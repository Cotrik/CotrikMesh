/*
 * DiagnalCollapseSimplifier.cpp
 *
 *  Created on: Jan 20, 2019
 *      Author: cotrik
 */

#include "DiagnalCollapseSimplifier.h"


struct collapsableDiagonal {
    size_t vid;
    size_t target_vid;
    size_t fid;
    double ranking;
};
bool opComp(const collapsableDiagonal& a, const collapsableDiagonal& b) {
    return a.ranking < b.ranking;
}

DiagnalCollapseSimplifier::DiagnalCollapseSimplifier(Mesh& mesh) : Simplifier(mesh) {
    // TODO Auto-generated constructor stub
    mu.SetMembers(mesh);

}

DiagnalCollapseSimplifier::~DiagnalCollapseSimplifier() {
    // TODO Auto-generated destructor stub
}

void DiagnalCollapseSimplifier::Run3(std::set<size_t>& canceledFids) {
    for (auto& f : mesh.F) {
        auto& v0 = mesh.V.at(f.Vids[0]);
        auto& v1 = mesh.V.at(f.Vids[1]);
        auto& v2 = mesh.V.at(f.Vids[2]);
        auto& v3 = mesh.V.at(f.Vids[3]);
        if (v0.type == CORNER || v1.type == CORNER || v2.type == CORNER || v3.type == CORNER) continue;

        if (v0.type == FEATURE && v1.type == FEATURE && v2.type == FEATURE && v0.label == v1.label && v1.label == v2.label) {
            canceledFids.insert(f.id);
            Collapse(v0.id, v2.id);
            return;
        } else if (v1.type == FEATURE && v2.type == FEATURE && v3.type == FEATURE && v1.label == v2.label && v2.label == v3.label) {
            canceledFids.insert(f.id);
            Collapse(v1.id, v3.id);
            return;
        } else if (v2.type == FEATURE && v3.type == FEATURE && v0.type == FEATURE && v2.label == v3.label && v3.label == v0.label) {
            canceledFids.insert(f.id);
            Collapse(v2.id, v0.id);
            return;
        } else if (v3.type == FEATURE && v0.type == FEATURE && v1.type == FEATURE && v3.label == v0.label && v0.label == v1.label) {
            canceledFids.insert(f.id);
            Collapse(v3.id, v1.id);
            return;
        }
    }
}

void DiagnalCollapseSimplifier::CollapseSinglets(std::set<size_t>& canceledFids) {
    for (auto& f: mesh.F) f.isVisited = false;
    for (auto& f : mesh.F) {
        if (f.isVisited) continue;
        f.isVisited = true;
        bool collapsed = false;
        auto& v0 = mesh.V.at(f.Vids[0]);
        auto& v1 = mesh.V.at(f.Vids[1]);
        auto& v2 = mesh.V.at(f.Vids[2]);
        auto& v3 = mesh.V.at(f.Vids[3]);
        if (v0.type == CORNER || v1.type == CORNER || v2.type == CORNER || v3.type == CORNER) continue;

        if (v0.type == FEATURE && v1.type == FEATURE && v2.type == FEATURE && v0.label == v1.label && v1.label == v2.label) {
            canceledFids.insert(f.id);
            Collapse(v0.id, v2.id);
            collapsed = true;
            // return;
        } else if (v1.type == FEATURE && v2.type == FEATURE && v3.type == FEATURE && v1.label == v2.label && v2.label == v3.label) {
            canceledFids.insert(f.id);
            collapsed = true;
            Collapse(v1.id, v3.id);
            // return;
        } else if (v2.type == FEATURE && v3.type == FEATURE && v0.type == FEATURE && v2.label == v3.label && v3.label == v0.label) {
            canceledFids.insert(f.id);
            collapsed = true;
            Collapse(v2.id, v0.id);
            // return;
        } else if (v3.type == FEATURE && v0.type == FEATURE && v1.type == FEATURE && v3.label == v0.label && v0.label == v1.label) {
            canceledFids.insert(f.id);
            Collapse(v3.id, v1.id);
            collapsed = true;
            // return;
        }

        if (collapsed) {
            for (auto nfid: f.N_Fids) {
                mesh.F.at(nfid).isVisited = true;
            }
        }
    }
}

void DiagnalCollapseSimplifier::Run(std::set<size_t>& canceledFids) {
    for (auto valence = Simplifier::maxValence - 1; valence > 4; --valence) {
        for (auto& f : mesh.F) {
            auto& v0 = mesh.V.at(f.Vids[0]);
            auto& v1 = mesh.V.at(f.Vids[1]);
            auto& v2 = mesh.V.at(f.Vids[2]);
            auto& v3 = mesh.V.at(f.Vids[3]);
            //if (v0.type == CORNER || v1.type == CORNER || v2.type == CORNER || v3.type == CORNER) continue;
            if (v1.N_Fids.size() >= valence && v3.N_Fids.size() >= valence && !v1.isBoundary && !v3.isBoundary) {
                if (v0.type == CORNER || v2.type == CORNER) continue;
                if (CanCollapseDiagnal(v0, v2)) {
                    Collapse(v0.id, v2.id);
                    canceledFids.insert(f.id);
                    return;
                } else if (CanCollapseDiagnal(v2, v0)) {
                    Collapse(v2.id, v0.id);
                    canceledFids.insert(f.id);
                    return;
                }
            }
            //else
                if (v0.N_Fids.size() >= valence && v2.N_Fids.size() >= valence && !v0.isBoundary && !v2.isBoundary) {
                if (v1.type == CORNER || v3.type == CORNER) continue;
                if (CanCollapseDiagnal(v1, v3)) {
                    Collapse(v1.id, v3.id);
                    canceledFids.insert(f.id);
                    return;
                } else if (CanCollapseDiagnal(v3, v1)) {
                    Collapse(v3.id, v1.id);
                    canceledFids.insert(f.id);
                    return;
                }
            }
        }
    }
}

void DiagnalCollapseSimplifier::RunCollective(std::set<size_t>& canceledFids) {
    std::vector<collapsableDiagonal> diagonals;
    /*for (auto& f: mesh.F) {
        auto& v0 = mesh.V.at(f.Vids[0]);
        auto& v1 = mesh.V.at(f.Vids[1]);
        auto& v2 = mesh.V.at(f.Vids[2]);
        auto& v3 = mesh.V.at(f.Vids[3]);
        if (v0.isBoundary || v1.isBoundary || v2.isBoundary || v3.isBoundary) continue;
        if (v0.type == FEATURE || v1.type == FEATURE || v2.type == FEATURE || v3.type == FEATURE) continue;
        
        collapsableDiagonal d1;
        d1.vid = v0.id;
        d1.target_vid = v2.id;
        d1.fid = f.id;
        d1.ranking = CalculateRanking(std::vector<size_t>{v0.id, v2.id, v1.id, v3.id});

        collapsableDiagonal d2;
        d2.vid = v1.id;
        d2.target_vid = v3.id;
        d2.fid = f.id;
        d2.ranking = CalculateRanking(std::vector<size_t>{v1.id, v3.id, v0.id, v2.id});

        if (d1.ranking < d2.ranking) {
            
            diagonals.push_back(d1);
        } else {
            diagonals.push_back(d2);
        }
    }
    int min_index = 0;
    double min_ranking = 1.0;
    for (int i = 0; i < diagonals.size(); i++) {
        auto& d = diagonals.at(i);
        if (d.ranking < min_ranking) {
            min_index = i;
            min_ranking = d.ranking;
        }
    }
    auto& d = diagonals.at(min_index);
    if (d.ranking > 0.55) return; 
    Collapse(d.vid, d.target_vid);
    canceledFids.insert(d.fid);
    return;*/
    for (auto& f: mesh.F) f.isVisited = false;
    for (auto valence = Simplifier::maxValence - 1; valence > 4; --valence) {
        for (auto& f: mesh.F) {
            if (f.isVisited) continue;
            // f.isVisited = true;
            auto& v0 = mesh.V.at(f.Vids[0]);
            auto& v1 = mesh.V.at(f.Vids[1]);
            auto& v2 = mesh.V.at(f.Vids[2]);
            auto& v3 = mesh.V.at(f.Vids[3]);
            //if (v0.type == CORNER || v1.type == CORNER || v2.type == CORNER || v3.type == CORNER) continue;
            if (v1.N_Fids.size() >= valence && v3.N_Fids.size() >= valence && !v1.isBoundary && !v3.isBoundary) {
                if (v0.type == CORNER || v2.type == CORNER) continue;
                if (CanCollapseDiagnal(v0, v2)) {
                    collapsableDiagonal d;
                    d.vid = v0.id;
                    d.target_vid = v2.id;
                    double min = mu.GetVertexEnergy(v0.id) + mu.GetVertexEnergy(v2.id);
                    double max = mu.GetVertexEnergy(v1.id) + mu.GetVertexEnergy(v3.id);
                    double normalizedArea = mu.GetMeshArea() / mu.GetFaceArea(f.id);
                    d.ranking = min / (max * normalizedArea);
                    diagonals.push_back(d);

                    // Collapse(v0.id, v2.id);
                    // CollapseDiagnal(v0, v2);
                    canceledFids.insert(f.id);
                    f.isVisited = true;
                    for (auto fid: f.N_Fids) {
                        mesh.F.at(fid).isVisited = true;
                    }
                    continue;
                    // return;
                } else if (CanCollapseDiagnal(v2, v0)) {
                    collapsableDiagonal d;
                    d.vid = v2.id;
                    d.target_vid = v0.id;
                    double min = mu.GetVertexEnergy(v0.id) + mu.GetVertexEnergy(v2.id);
                    double max = mu.GetVertexEnergy(v1.id) + mu.GetVertexEnergy(v3.id);
                    double normalizedArea = mu.GetMeshArea() / mu.GetFaceArea(f.id);
                    d.ranking = min / (max * normalizedArea);
                    diagonals.push_back(d);
                    
                    // Collapse(v2.id, v0.id);
                    // CollapseDiagnal(v2, v0);
                    canceledFids.insert(f.id);
                    f.isVisited = true;
                    for (auto fid: f.N_Fids) {
                        mesh.F.at(fid).isVisited = true;
                    }
                    continue;
                    // return;
                }
            }
            //else
                if (v0.N_Fids.size() >= valence && v2.N_Fids.size() >= valence && !v0.isBoundary && !v2.isBoundary) {
                if (v1.type == CORNER || v3.type == CORNER) continue;
                if (CanCollapseDiagnal(v1, v3)) {
                    collapsableDiagonal d;
                    d.vid = v1.id;
                    d.target_vid = v3.id;
                    double min = mu.GetVertexEnergy(v1.id) + mu.GetVertexEnergy(v3.id);
                    double max = mu.GetVertexEnergy(v0.id) + mu.GetVertexEnergy(v2.id);
                    double normalizedArea = mu.GetMeshArea() / mu.GetFaceArea(f.id);
                    d.ranking = min / (max * normalizedArea);
                    diagonals.push_back(d);
                    
                    // Collapse(v1.id, v3.id);
                    // CollapseDiagnal(v1, v3);
                    canceledFids.insert(f.id);
                    f.isVisited = true;
                    for (auto fid: f.N_Fids) {
                        mesh.F.at(fid).isVisited = true;
                    }
                    // return;
                } else if (CanCollapseDiagnal(v3, v1)) {
                    collapsableDiagonal d;
                    d.vid = v3.id;
                    d.target_vid = v1.id;
                    double min = mu.GetVertexEnergy(v1.id) + mu.GetVertexEnergy(v3.id);
                    double max = mu.GetVertexEnergy(v0.id) + mu.GetVertexEnergy(v2.id);
                    double normalizedArea = mu.GetMeshArea() / mu.GetFaceArea(f.id);
                    d.ranking = min / (max * normalizedArea);
                    diagonals.push_back(d);
                    
                    // Collapse(v3.id, v1.id);
                    // CollapseDiagnal(v3, v1);
                    canceledFids.insert(f.id);
                    f.isVisited = true;
                    for (auto fid: f.N_Fids) {
                        mesh.F.at(fid).isVisited = true;
                    }
                    // return;
                }
            }
        }
    }
    std::sort(diagonals.begin(), diagonals.end(), opComp);
    for (auto d: diagonals) {
        Collapse(d.vid, d.target_vid);
    }
}

void DiagnalCollapseSimplifier::CollapseDiagnal(const Vertex& v0, const Vertex& v2) {
    auto target_vid = v0.type == CORNER ? v0.id : v2.id;
    auto vid = v0.type == CORNER ? v2.id : v0.id;
    Collapse(vid, target_vid);
}

// bool DiagnalCollapseSimplifier::CanCollapseDiagnal(const Vertex& v, const Vertex& target_v) {
//    return !v.isBoundary && !target_v.isBoundary && v.N_Fids.size() == 3 && target_v.N_Fids.size() == 3;
// }
bool DiagnalCollapseSimplifier::CanCollapseDiagnal(const Vertex& v, const Vertex& target_v) {
    return !v.isBoundary && !target_v.isBoundary && v.N_Fids.size() == 3 && (target_v.N_Fids.size() >= Simplifier::minValence && target_v.N_Fids.size() <= Simplifier::maxValence);
}

void DiagnalCollapseSimplifier::GetDiagonalCollapseOps(std::multiset<SimplificationOperationStruct, bool(*)(SimplificationOperationStruct, SimplificationOperationStruct)>& SimplificationOps) {
    for (auto valence = Simplifier::maxValence - 1; valence > 4; --valence) {
        for (auto& f: mesh.F) {
            auto& v0 = mesh.V.at(f.Vids[0]);
            auto& v1 = mesh.V.at(f.Vids[1]);
            auto& v2 = mesh.V.at(f.Vids[2]);
            auto& v3 = mesh.V.at(f.Vids[3]);
            size_t source_id, target_id;
            bool OpExists = false;
            if (v1.N_Fids.size() >= valence && v3.N_Fids.size() >= valence && !v1.isBoundary && !v3.isBoundary) {
                if (v0.type == CORNER || v2.type == CORNER) continue;
                if (CanCollapseDiagnal(v0, v2)) {
                    source_id = v0.id;
                    target_id = v2.id;
                    OpExists = true;
                } else if (CanCollapseDiagnal(v2, v0)) {
                    source_id = v2.id;
                    target_id = v0.id;
                    OpExists = true;
                }
            } else if (v0.N_Fids.size() >= valence && v2.N_Fids.size() >= valence && !v0.isBoundary && !v2.isBoundary) {
                if (v1.type == CORNER || v3.type == CORNER) continue;
                if (CanCollapseDiagnal(v1, v3)) {
                    source_id = v1.id;
                    target_id = v3.id;
                    OpExists = true;
                } else if (CanCollapseDiagnal(v3, v1)) {
                    source_id = v3.id;
                    target_id = v1.id;
                    OpExists = true;
                }
            }
            if (OpExists) {
                SimplificationOperationStruct Op;
                Op.type = "Strict_Diagonal_Collapse";
                auto& source = mesh.V.at(source_id);
                auto& target = mesh.V.at(target_id);
                // Op.updateVertexIds.push_back(target_id);
                // Op.updatedVertexPos.push_back(0.5 * (source.xyz() + target.xyz()));
                // Op.profitability = mesh.GetQuadFaceArea(f.Vids) / mesh.totalArea;
                Op.profitability = glm::distance(source.xyz(), target.xyz());
                // Op.profitability = 1;
                for (auto fid: source.N_Fids) {
                    if (std::find(target.N_Fids.begin(), target.N_Fids.end(), fid) == target.N_Fids.end()) {
                        Face& n_f = mesh.F.at(fid);
                        Face newF;
                        for (int i = 0; i < n_f.Vids.size(); i++) {
                            n_f.Vids.at(i) == source.id ? newF.Vids.push_back(target.id) : newF.Vids.push_back(n_f.Vids.at(i));
                        }
                        Op.newFaces.push_back(newF);
                        Op.canceledFids.insert(n_f.id);
                    } else {
                        Op.canceledFids.insert(fid);
                    }
                }
                SimplificationOps.insert(Op);
            }
        }
    }
}

double DiagnalCollapseSimplifier::CalculateRanking(std::vector<size_t>& Vids) {
    double alpha_q = 0.05;
    double alpha_d = 0.45;
    double alpha_v = 0.5;

    glm::dvec4 newV = CalculateQEM(Vids.at(0), Vids.at(1));
    double E_q = newV[3];
    double E_d = glm::distance(mesh.V.at(Vids.at(0)).xyz(), mesh.V.at(Vids.at(1)).xyz());
    double E_v = 0.0;

    auto& v1 = mesh.V.at(Vids.at(0));
    auto& v2 = mesh.V.at(Vids.at(1));
    auto& v3 = mesh.V.at(Vids.at(2));
    auto& v4 = mesh.V.at(Vids.at(3));
    
    int target_val = v1.N_Vids.size() + v2.N_Vids.size() - 2;
    int val_a = v1.N_Vids.size() >= 4 ? v1.N_Vids.size() - 4 : 4 - v1.N_Vids.size();
    int val_c = v2.N_Vids.size() >= 4 ? v2.N_Vids.size() - 4 : 4 - v2.N_Vids.size();
    int val_b = v3.N_Vids.size() >= 4 ? v3.N_Vids.size() - 4 : 4 - v3.N_Vids.size();
    int val_d = v4.N_Vids.size() >= 4 ? v4.N_Vids.size() - 4 : 4 - v4.N_Vids.size();
    int val_b_ = v3.N_Vids.size() >= 5 ? v3.N_Vids.size() - 5 : 5 - v3.N_Vids.size();
    int val_d_ = v4.N_Vids.size() >= 5 ? v4.N_Vids.size() - 5 : 5 - v4.N_Vids.size();
    
    int val_term = 0;
    if (target_val > val_a) val_term += target_val - val_a;
    if (target_val > val_c) val_term += target_val - val_c;
    if (val_b_ > val_b_) val_term += val_b_ - val_b;
    if (val_d_ > val_d) val_term += val_d_ - val_d;

    E_v = (double) val_term;

    double r = (alpha_q * (1 - exp(-E_q))) + (alpha_d * (1 - exp(-E_d))) + (alpha_v * (1 - exp(-E_v)));
    // std::cout << "ranking: " << r << std::endl;
    return r;
}

glm::dvec4 DiagnalCollapseSimplifier::CalculateQEM(size_t v1_id, size_t v2_id) {
    auto& v1 = mesh.V.at(v1_id);
    auto& v2 = mesh.V.at(v2_id);

    // std::cout << "Vertex v1: " << v1.id << " " << v1.x << " " << v1.y << " " << v1.z << std::endl;
    // std::cout << "Vertex v2: " << v2.id << " " << v2.x << " " << v2.y << " " << v2.z << std::endl;

    glm::dmat4 Q1(0.0);
    // std::cout << "Calculating Q1" << std::endl;
    for (auto fid: v1.N_Fids) {
        auto& f = mesh.F.at(fid);
        // std::cout << "neighbor face: " << f.id << " " << f.Vids.at(0) << " " << f.Vids.at(1) << " " << f.Vids.at(2) << " " << f.Vids.at(3) << std::endl;
        int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v1.id));
        auto& v_a = mesh.V.at(f.Vids.at((idx+1)%f.Vids.size()));
        auto& v_b = mesh.V.at(f.Vids.at((idx+3)%f.Vids.size()));

        // std::cout << "Vertex v_a: " << v_a.id << " " << v_a.x << " " << v_a.y << " " << v_a.z << std::endl;
        // std::cout << "Vertex v_b: " << v_b.id << " " << v_b.x << " " << v_b.y << " " << v_b.z << std::endl;
        
        glm::dvec3 A = v_a.xyz() - v1.xyz();
        glm::dvec3 B = v_b.xyz() - v1.xyz();
        // std::cout << "A: " << A.x << " " << A.y << " " << A.z << std::endl;
        // std::cout << "B: " << B.x << " " << B.y << " " << B.z << std::endl;
        
        glm::dvec3 n = glm::normalize(glm::cross(A, B));
        // std::cout << "normal: " << n.x << " " << n.y << " " << n.z << std::endl;

        glm::dvec4 p(n, glm::dot(n,v1.xyz()));
        // std::cout << "plane: " << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << std::endl;
        
        glm::dmat4 Q;
        Q[0] = glm::dvec4(p[0]*p[0], p[0]*p[1], p[0]*p[2], p[0]*p[3]);
        Q[1] = glm::dvec4(p[1]*p[0], p[1]*p[1], p[1]*p[2], p[1]*p[3]);
        Q[2] = glm::dvec4(p[2]*p[0], p[2]*p[1], p[2]*p[2], p[2]*p[3]);
        Q[3] = glm::dvec4(p[3]*p[0], p[3]*p[1], p[3]*p[2], p[3]*p[3]);

        // std::cout << "Q: " << std::endl;
        // std::cout << Q[0][0] << " " << Q[1][0] << " " << Q[2][0] << " " << Q[3][0] << std::endl;
        // std::cout << Q[0][1] << " " << Q[1][1] << " " << Q[2][1] << " " << Q[3][1] << std::endl;
        // std::cout << Q[0][2] << " " << Q[1][2] << " " << Q[2][2] << " " << Q[3][2] << std::endl;
        // std::cout << Q[0][3] << " " << Q[1][3] << " " << Q[2][3] << " " << Q[3][3] << std::endl;

        Q1 += Q;

        // std::cout << "Q1: " << std::endl;
        // std::cout << Q1[0][0] << " " << Q1[1][0] << " " << Q1[2][0] << " " << Q1[3][0] << std::endl;
        // std::cout << Q1[0][1] << " " << Q1[1][1] << " " << Q1[2][1] << " " << Q1[3][1] << std::endl;
        // std::cout << Q1[0][2] << " " << Q1[1][2] << " " << Q1[2][2] << " " << Q1[3][2] << std::endl;
        // std::cout << Q1[0][3] << " " << Q1[1][3] << " " << Q1[2][3] << " " << Q1[3][3] << std::endl;
        // std::cout << "***********************************" << std::endl;
    }
    // std::cout << "Q1: " << std::endl;
    // std::cout << Q1[0][0] << " " << Q1[1][0] << " " << Q1[2][0] << " " << Q1[3][0] << std::endl;
    // std::cout << Q1[0][1] << " " << Q1[1][1] << " " << Q1[2][1] << " " << Q1[3][1] << std::endl;
    // std::cout << Q1[0][2] << " " << Q1[1][2] << " " << Q1[2][2] << " " << Q1[3][2] << std::endl;
    // std::cout << Q1[0][3] << " " << Q1[1][3] << " " << Q1[2][3] << " " << Q1[3][3] << std::endl;


    glm::dmat4 Q2(0.0);
    // std::cout << "Calculating Q2" << std::endl;
    for (auto fid: v2.N_Fids) {
        auto& f = mesh.F.at(fid);
        // std::cout << "neighbor face: " << f.id << " " << f.Vids.at(0) << " " << f.Vids.at(1) << " " << f.Vids.at(2) << " " << f.Vids.at(3) << std::endl;
        int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v2.id));
        auto& v_a = mesh.V.at(f.Vids.at((idx+1)%f.Vids.size()));
        auto& v_b = mesh.V.at(f.Vids.at((idx+3)%f.Vids.size()));

        // std::cout << "Vertex v_a: " << v_a.id << " " << v_a.x << " " << v_a.y << " " << v_a.z << std::endl;
        // std::cout << "Vertex v_b: " << v_b.id << " " << v_b.x << " " << v_b.y << " " << v_b.z << std::endl;
        
        glm::dvec3 A = v_a.xyz() - v2.xyz();
        glm::dvec3 B = v_b.xyz() - v2.xyz();
        // std::cout << "A: " << A.x << " " << A.y << " " << A.z << std::endl;
        // std::cout << "B: " << B.x << " " << B.y << " " << B.z << std::endl;
        
        glm::dvec3 n = glm::normalize(glm::cross(A, B));
        // std::cout << "normal: " << n.x << " " << n.y << " " << n.z << std::endl;
        glm::dvec4 p(n, glm::dot(n,v2.xyz()));
        // std::cout << "plane: " << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << std::endl;
        
        glm::dmat4 Q;
        Q[0] = glm::dvec4(p[0]*p[0], p[0]*p[1], p[0]*p[2], p[0]*p[3]);
        Q[1] = glm::dvec4(p[1]*p[0], p[1]*p[1], p[1]*p[2], p[1]*p[3]);
        Q[2] = glm::dvec4(p[2]*p[0], p[2]*p[1], p[2]*p[2], p[2]*p[3]);
        Q[3] = glm::dvec4(p[3]*p[0], p[3]*p[1], p[3]*p[2], p[3]*p[3]);

        // std::cout << "Q: " << std::endl;
        // std::cout << Q[0][0] << " " << Q[1][0] << " " << Q[2][0] << " " << Q[3][0] << std::endl;
        // std::cout << Q[0][1] << " " << Q[1][1] << " " << Q[2][1] << " " << Q[3][1] << std::endl;
        // std::cout << Q[0][2] << " " << Q[1][2] << " " << Q[2][2] << " " << Q[3][2] << std::endl;
        // std::cout << Q[0][3] << " " << Q[1][3] << " " << Q[2][3] << " " << Q[3][3] << std::endl;

        Q2 += Q;
        // std::cout << "Q2: " << std::endl;
        // std::cout << Q2[0][0] << " " << Q2[1][0] << " " << Q2[2][0] << " " << Q2[3][0] << std::endl;
        // std::cout << Q2[0][1] << " " << Q2[1][1] << " " << Q2[2][1] << " " << Q2[3][1] << std::endl;
        // std::cout << Q2[0][2] << " " << Q2[1][2] << " " << Q2[2][2] << " " << Q2[3][2] << std::endl;
        // std::cout << Q2[0][3] << " " << Q2[1][3] << " " << Q2[2][3] << " " << Q2[3][3] << std::endl;
        // std::cout << "***********************************" << std::endl;
    }

    // std::cout << "Q2: " << std::endl;
    // std::cout << Q2[0][0] << " " << Q2[1][0] << " " << Q2[2][0] << " " << Q2[3][0] << std::endl;
    // std::cout << Q2[0][1] << " " << Q2[1][1] << " " << Q2[2][1] << " " << Q2[3][1] << std::endl;
    // std::cout << Q2[0][2] << " " << Q2[1][2] << " " << Q2[2][2] << " " << Q2[3][2] << std::endl;
    // std::cout << Q2[0][3] << " " << Q2[1][3] << " " << Q2[2][3] << " " << Q2[3][3] << std::endl;

    glm::dmat4 Q_ = Q1+Q2;
    // std::cout << "Q_: " << std::endl;
    // std::cout << Q_[0][0] << " " << Q_[1][0] << " " << Q_[2][0] << " " << Q_[3][0] << std::endl;
    // std::cout << Q_[0][1] << " " << Q_[1][1] << " " << Q_[2][1] << " " << Q_[3][1] << std::endl;
    // std::cout << Q_[0][2] << " " << Q_[1][2] << " " << Q_[2][2] << " " << Q_[3][2] << std::endl;
    // std::cout << Q_[0][3] << " " << Q_[1][3] << " " << Q_[2][3] << " " << Q_[3][3] << std::endl;

    glm::dmat4 Q_v = Q1+Q2;
    Q_v[0][3] = 0.0;
    Q_v[1][3] = 0.0;
    Q_v[2][3] = 0.0;
    Q_v[3][3] = 1.0;
    // std::cout << "Q_v: " << std::endl;
    // std::cout << Q_v[0][0] << " " << Q_v[1][0] << " " << Q_v[2][0] << " " << Q_v[3][0] << std::endl;
    // std::cout << Q_v[0][1] << " " << Q_v[1][1] << " " << Q_v[2][1] << " " << Q_v[3][1] << std::endl;
    // std::cout << Q_v[0][2] << " " << Q_v[1][2] << " " << Q_v[2][2] << " " << Q_v[3][2] << std::endl;
    // std::cout << Q_v[0][3] << " " << Q_v[1][3] << " " << Q_v[2][3] << " " << Q_v[3][3] << std::endl;

    glm::dvec4 v_in(0.0, 0.0, 0.0, 1.0);

    // std::cout << "determinant: " << glm::determinant(Q_v) << std::endl; 
    double determinant = glm::determinant(Q_v);
    glm::dvec4 v_(0.5 * (v1.xyz() + v2.xyz()), 1.0);
    if (fabs(determinant) > 0.0) {
        glm::dmat4 inv_Q_v = glm::inverse(Q_v);
        // std::cout << "inv_Q_v: " << std::endl;
        // std::cout << inv_Q_v[0][0] << " " << inv_Q_v[1][0] << " " << inv_Q_v[2][0] << " " << inv_Q_v[3][0] << std::endl;
        // std::cout << inv_Q_v[0][1] << " " << inv_Q_v[1][1] << " " << inv_Q_v[2][1] << " " << inv_Q_v[3][1] << std::endl;
        // std::cout << inv_Q_v[0][2] << " " << inv_Q_v[1][2] << " " << inv_Q_v[2][2] << " " << inv_Q_v[3][2] << std::endl;
        // std::cout << inv_Q_v[0][3] << " " << inv_Q_v[1][3] << " " << inv_Q_v[2][3] << " " << inv_Q_v[3][3] << std::endl;

        v_ = inv_Q_v*v_in;
    }

    glm::dvec4 v_int = Q_*v_;
    // std::cout << "v_: " << v_[0] << " " << v_[1] << " " << v_[2] << " " << v_[3] << std::endl;
    // std::cout << "v_int: " << v_int[0] << " " << v_int[1] << " " << v_int[2] << " " << v_int[3] << std::endl;

    v_[3] = glm::dot(v_, v_int);
    return v_;
}


