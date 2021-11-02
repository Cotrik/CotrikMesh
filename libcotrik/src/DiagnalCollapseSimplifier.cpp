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
    double ranking;
};
bool opComp(const collapsableDiagonal& a, const collapsableDiagonal& b) {
    return a.ranking < b.ranking;
}

DiagnalCollapseSimplifier::DiagnalCollapseSimplifier(Mesh& mesh) : Simplifier(mesh) {
    // TODO Auto-generated constructor stub

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

void DiagnalCollapseSimplifier::GetDiagonalCollapseOps(std::multiset<SimplificationOperation, bool(*)(SimplificationOperation, SimplificationOperation)>& SimplificationOps) {
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
                SimplificationOperation Op;
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
