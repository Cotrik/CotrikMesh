#include <algorithm>
#include <map>
#include <time.h>
#include "ParallelFor.h"
#include "SemiGlobalSimplifier.h"

SemiGlobalSimplifier::SemiGlobalSimplifier() {}

SemiGlobalSimplifier::SemiGlobalSimplifier(Mesh& mesh_, MeshUtil& mu_, Smoother& smoother_) : mesh(mesh_), mu(mu_), smoother(smoother_) {
    // mu.SetMesh(mesh);
    // smoother.SetMesh(mesh);
}

SemiGlobalSimplifier::~SemiGlobalSimplifier() {}

void SemiGlobalSimplifier::CheckValidity() {
    if (mesh.V.size() == 0 || mesh.F.size() == 0 || mesh.C.size() == 0) {
        std::cout << "No mesh to use for Semi Global Simplifier." << std::endl;
        exit(0);
    }
}

void SemiGlobalSimplifier::SetMesh(Mesh& mesh_) {
    mesh = mesh_;
    mu.SetMesh(mesh);
    smoother.SetMesh(mesh);
}

void SemiGlobalSimplifier::SetIters(int iters_) {
    iters = iters_;
}

void SemiGlobalSimplifier::SetSimplificationOperations() {
    CheckValidity();
    
    // SetDiagonalCollapseOperations();
    SetDirectSeparatrixOperations(false);
}

void SemiGlobalSimplifier::FixBoundary() {
    CheckValidity();

    // int i = 0;
    // for (auto& f: mesh.F) {
    //     if (f.N_Fids.size() == 0) continue;
    //     bool skip = false;
    //     for (auto fvid: f.Vids) {
    //         auto& v = mesh.V.at(fvid);
    //         int featureCount = 0;
    //         for (auto vid: v.N_Vids) {
    //             if (mesh.V.at(vid).type == FEATURE) featureCount += 1;
    //         }
    //         if (featureCount > 2) {
    //             skip = true;
    //             break;
    //         }
    //     }
    //     if (skip) continue;
    //     for (int i = 0; i < f.Vids.size(); i++) {
    //         if ((mesh.V.at(f.Vids.at(i)).type == FEATURE || mesh.V.at(f.Vids.at(i)).isBoundary)
    //         && (mesh.V.at(f.Vids.at((i+1)%f.Vids.size())).type == FEATURE || mesh.V.at(f.Vids.at((i+1)%f.Vids.size())).isBoundary)
    //         && (mesh.V.at(f.Vids.at((i+2)%f.Vids.size())).type == FEATURE || mesh.V.at(f.Vids.at((i+2)%f.Vids.size())).isBoundary)) {
    //             std::cout << "fixing degenerate boundary" << std::endl;
    //             std::shared_ptr<SimplificationOperation> dc1 = std::make_shared<DiagonalCollapse>(mesh, mu, f.id, i, (i+2)%f.Vids.size());
    //             // dc1->SetRanking();
    //             dc1->PerformOperation();
    //             mesh.V.at(f.Vids.at((i+1)%f.Vids.size())).type = REGULAR;
    //             mesh.V.at(f.Vids.at((i+1)%f.Vids.size())).isBoundary = false;
    //             break;
    //         }
    //     }
    // }
    for (int i = 0; i < mesh.V.size(); i++) {
        // if (i >= iters) break;
        auto& v = mesh.V.at(i);
        int valence = v.N_Fids.size();
        bool failed = false;
        if (v.type == FEATURE && valence > 4) {
            int featureCount = 0;
            for (auto vid: v.N_Vids) {
                if (mesh.V.at(vid).type == FEATURE) featureCount += 1;
            }
            if (featureCount != 2) continue;
            bool performVertexSplit = true;
            for (int j = 0; j < v.N_Eids.size(); j++) {
                auto& e = mesh.E.at(v.N_Eids.at(j));
                if (mesh.V.at(e.Vids.at(0)).type == FEATURE && mesh.V.at(e.Vids.at(1)).type == FEATURE) continue;
                int count = 0;
                for (auto fid: e.N_Fids) {
                    auto& f = mesh.F.at(fid);
                    for (auto vid: f.Vids) {
                        if (mesh.V.at(vid).type == FEATURE) count += 1;
                    }
                }
                if (count == 4) performVertexSplit = false;
            }
            for (int j = 0; j < v.N_Fids.size(); j++) {
                int count = 0;
                auto& f = mesh.F.at(v.N_Fids.at(j));
                for (auto vid: f.Vids) {
                    if (mesh.V.at(vid).type == FEATURE) count += 1;
                }
                if (count > 2) performVertexSplit = false;
            }
            if (performVertexSplit && v.N_Vids.size() > 5) {
                std::shared_ptr<SimplificationOperation> s = std::make_shared<VertexSplit>(mesh, mu, v.id);
                s->PerformOperation();
                continue;
            }
            while (valence > 4) {
                bool breakLoop = false;
                for (int j = 0; j < v.N_Eids.size(); j++) {
                    auto& e = mesh.E.at(v.N_Eids.at(j));
                    if (mesh.V.at(e.Vids.at(0)).type == FEATURE && mesh.V.at(e.Vids.at(1)).type == FEATURE) continue;
                    int count = 0;
                    for (auto fid: e.N_Fids) {
                        auto& f = mesh.F.at(fid);
                        for (auto vid: f.Vids) {
                            if (mesh.V.at(vid).type == FEATURE) count += 1;
                        }
                    }
                    if (count == 4) continue;
                    bool clockwise = false;
                    for (auto fid: e.N_Fids) {
                        auto& f = mesh.F.at(fid);
                        int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
                        if (mesh.V.at(f.Vids.at((idx+1)%f.Vids.size())).type == FEATURE) {
                            clockwise = true;
                        }
                    }
                    std::shared_ptr<SimplificationOperation> s = std::make_shared<EdgeRotation>(mesh, mu, e.id, clockwise);
                    s->PerformOperation();
                }
                if (v.N_Fids.size() == valence) {
                    failed = true;
                    break;
                }
                valence = v.N_Fids.size();
            }
            // i += 1;
            if (failed) break;
        }
        if (v.isBoundary && v.N_Fids.size() != 2) {

        }
    }
}

void SemiGlobalSimplifier::SetDiagonalCollapseOperations() {
    CheckValidity();

    for (auto& f: mesh.F) {
        if (f.N_Fids.size() == 0 || f.Vids.empty()) continue;
        for (int i = 0; i < f.Vids.size(); i++) {
            if (mesh.V.at(f.Vids.at(i)).N_Fids.size() == 3 && mesh.V.at(f.Vids.at((i+2)%f.Vids.size())).N_Fids.size() == 3) {
                    std::shared_ptr<SimplificationOperation> dc1 = std::make_shared<DiagonalCollapse>(mesh, mu, f.id, i, (i+2)%f.Vids.size());
                    dc1->SetRanking();
                    Ops.push_back(dc1);
                break;
            }
        }
        
        // std::unique_ptr<SimplificationOperation> dc2 = std::make_unique<DiagonalCollapse>(mesh, mu, f.id, 1, 3);
        // dc2->SetRanking();
        // Ops.push_back(std::move(dc2));
    }
    std::cout << Ops.size() << std::endl;
    int i = 0;
    for (auto& op: Ops) {
        op->PerformOperation();
        i += 1;
        if (i >= iters) break;        
    }
}

void SemiGlobalSimplifier::SetBoundaryDirectSeparatrixOperations(bool looseCollapse) {
    CheckValidity();

    Op_Q.setMaxQueueOn();
    Op_Q.setSpecialComparisonOn();
    for (auto& v: mesh.V) {
        if (v.N_Vids.size() != 4 || v.type == FEATURE || v.isBoundary) continue;
        std::vector<size_t> c1, c2;
        for (auto vid: v.N_Vids) mesh.V.at(vid).N_Vids.size() == 3 ? c1.push_back(vid) : c2.push_back(vid);
        if (c1.size() < 2) continue;
        auto& s3_v1 = mesh.V.at(c1.at(0));
        auto& s3_v2 = mesh.V.at(c1.at(1));
        auto& sn_v1= mesh.V.at(c2.at(0));
        auto& sn_v2= mesh.V.at(c2.at(1));
        if (mu.GetDifference(s3_v1.N_Fids, s3_v2.N_Fids).size() != s3_v1.N_Fids.size()) continue;
        if (sn_v1.type == FEATURE || sn_v1.isBoundary || sn_v2.type == FEATURE || sn_v2.isBoundary) {
            bool skip = false;
            for (auto cvid: c2) {
                int count = 0;
                for (auto eid: v.N_Eids) {
                    auto& e = mesh.E.at(eid);
                    if ((e.Vids.at(0) == v.id && e.Vids.at(1) == cvid) || (e.Vids.at(1) == v.id && e.Vids.at(0) == cvid)) {
                        for (auto efid: e.N_Fids) {
                            auto& f = mesh.F.at(efid);
                            for (auto fvid: f.Vids) {
                                if (mesh.V.at(fvid).type == FEATURE || mesh.V.at(fvid).isBoundary) count += 1;
                            }
                        }
                        if (count == 4) skip = true;
                        break;
                    }
                }
            }
            if (skip) continue;

            std::shared_ptr<SimplificationOperation> ds = std::make_shared<DirectSeparatrixCollapse>(mesh, mu, v.id, c1, c2, looseCollapse);
            ds->SetRanking();
            if (ds->ranking < 0) continue;
            Op_Q.insert(ds->ranking, v.id, ds);
        }
    }
    int i = 0;
    std::cout << Op_Q.size() << " boundary direct separatrix operations" << std::endl;
    while (!Op_Q.empty()) {
        auto& op = Op_Q.pop();
        op->PerformOperation();
        // smoother.Smooth(op->smoothV);
        for (auto key: op->toUpdate) {
            auto nop = Op_Q.getByKey(key);
            if (!nop) continue;
            nop->SetRanking(op->GetLocation());
            Op_Q.update(nop->ranking, key);
        }
    }
    // std::cout << "Smoothing mesh" << std::endl;
    // std::vector<size_t> smoothv;
    // for (auto& v: mesh.V) smoothv.push_back(v.id);
    // smoother.Smooth(smoothv);
}


void SemiGlobalSimplifier::SetDirectSeparatrixOperations(bool looseCollapse) {
    CheckValidity();

    Op_Q.setMaxQueueOn();
    Op_Q.setSpecialComparisonOn();
    for (auto& v: mesh.V) {
        if (v.N_Vids.size() != 4 || v.type == FEATURE || v.isBoundary) continue;
        std::vector<size_t> c1, c2;
        for (auto vid: v.N_Vids) mesh.V.at(vid).N_Vids.size() == 3 ? c1.push_back(vid) : c2.push_back(vid);
        if (c1.size() < 2) continue;
        auto& s3_v1 = mesh.V.at(c1.at(0));
        auto& s3_v2 = mesh.V.at(c1.at(1));
        auto& sn_v1= mesh.V.at(c2.at(0));
        auto& sn_v2= mesh.V.at(c2.at(1));
        if (mu.GetDifference(s3_v1.N_Fids, s3_v2.N_Fids).size() != s3_v1.N_Fids.size()) continue;
        if (sn_v1.type == FEATURE || sn_v1.isBoundary || sn_v2.type == FEATURE || sn_v2.isBoundary) continue;

        std::shared_ptr<SimplificationOperation> ds = std::make_shared<DirectSeparatrixCollapse>(mesh, mu, v.id, c1, c2, looseCollapse);
        ds->SetRanking();
        if (ds->ranking < 0) continue;
        Op_Q.insert(ds->ranking, v.id, ds);
    }
    int i = 0;
    std::cout << Op_Q.size() << " direct separatrix operations" << std::endl;
    while (!Op_Q.empty()) {
        auto& op = Op_Q.pop();
        op->PerformOperation();
        // smoother.Smooth(op->smoothV);
        for (auto key: op->toUpdate) {
            auto nop = Op_Q.getByKey(key);
            if (!nop) continue;
            nop->SetRanking(op->GetLocation());
            Op_Q.update(nop->ranking, key);
        }
    }
    // std::cout << "Smoothing mesh" << std::endl;
    // std::vector<size_t> smoothv;
    // for (auto& v: mesh.V) smoothv.push_back(v.id);
    // smoother.Smooth(smoothv);

}

void SemiGlobalSimplifier::SetSeparatrixOperations() {
    CheckValidity();

    BaseComplexQuad bc(mesh);
    Op_Q.clear();
    Op_Q.setMinQueueOn();
    Op_Q.setSpecialComparisonOff();
    // Op_Q.setMaxQueueOn();
    for (auto& v: mesh.V) {
        if (!v.isSingularity || v.N_Fids.empty()) continue;
        std::vector<SingularityLink> links = TraceSingularityLinks(v, bc);
        for (int i = 0; i < links.size(); i++) {
            auto& linkV = links.at(i).linkVids;
            auto& linkE = links.at(i).linkEids;
            auto& v_front = mesh.V.at(linkV.front());
            auto& v_back = mesh.V.at(linkV.back());
            if (v_front.isBoundary || v_back.isBoundary) continue;
            if (v_front.N_Fids.size() != 3 || v_back.N_Fids.size() != 3) continue;
            if (mesh.V.at(GetDiagonalV(v_front.id, GetFaceId(v_front.id, linkV.at(1)))).N_Fids.size() <= 4 || mesh.V.at(GetDiagonalV(v_back.id, GetFaceId(v_back.id, linkV.at(linkV.size()-2)))).N_Fids.size() <= 4) continue;

            std::shared_ptr<SimplificationOperation> s = std::make_shared<SeparatrixCollapse>(mesh, mu, linkV, linkE);
            s->SetRanking();
            Op_Q.insert(s->ranking, s->GetCenterId(), s);
        }
    }
    std::cout << Op_Q.size() << " separatrix operations" << std::endl;
    int i = 0;
    while (!Op_Q.empty()) {
        auto& op = Op_Q.pop();
        op->PerformOperation();
        // std::cout << "Getting to update operations" << std::endl;
        for (auto vid: op->toUpdate) {
            if (!mesh.V.at(vid).isSingularity || mesh.V.at(vid).N_Fids.empty()) continue;
            // std::cout << mesh.V.at(vid).N_Vids.size() << " " << mesh.V.at(vid).N_Eids.size() << " " << mesh.V.at(vid).N_Fids.size() << std::endl;
            std::vector<SingularityLink> links = TraceSingularityLinks(mesh.V.at(vid), bc);
            // std::cout << "links " << links.size() << std::endl;
            for (int i = 0; i < links.size(); i++) {
                auto& linkV = links.at(i).linkVids;
                auto& linkE = links.at(i).linkEids;
                auto& v_front = mesh.V.at(linkV.front());
                auto& v_back = mesh.V.at(linkV.back());
                if (v_front.isBoundary || v_back.isBoundary) continue;
                if (v_front.N_Fids.size() != 3 || v_back.N_Fids.size() != 3) continue;
                if (mesh.V.at(GetDiagonalV(v_front.id, GetFaceId(v_front.id, linkV.at(1)))).N_Fids.size() <= 4 || mesh.V.at(GetDiagonalV(v_back.id, GetFaceId(v_back.id, linkV.at(linkV.size()-2)))).N_Fids.size() <= 4) continue;

                // std::cout << "making a separatrix operation" << std::endl;
                std::shared_ptr<SimplificationOperation> s = std::make_shared<SeparatrixCollapse>(mesh, mu, linkV, linkE);
                // std::cout << "made a separatrix operation" << std::endl;
                s->SetRanking();
                Op_Q.insert(s->ranking, s->GetCenterId(), s);
            }
            // std::cout << "got updated operation" << std::endl;
        }
        i += 1;
        // if (i > iters) break;
    }
    // std::cout << "Finished all separatrix collapse operations" << std::endl;
    // std::cout << "Smoothing mesh" << std::endl;
    // std::vector<size_t> smoothv;
    // for (auto& v: mesh.V) {
    //     if (v.N_Vids.empty() || v.N_Eids.empty() || v.N_Fids.empty()) continue;
    //     smoothv.push_back(v.id);
    // }
    // smoother.Smooth(smoothv);
}

void SemiGlobalSimplifier::SetBoundarySeparatrixOperations() {
    CheckValidity();

    BaseComplexQuad bc(mesh);
    Op_Q.clear();
    Op_Q.setMinQueueOn();
    Op_Q.setSpecialComparisonOff();
    // Op_Q.setMaxQueueOn();
    for (auto& v: mesh.V) {
        if (!v.isSingularity || v.N_Fids.empty()) continue;
        std::vector<SingularityLink> links = TraceSingularityLinks(v, bc);
        for (int i = 0; i < links.size(); i++) {
            auto& linkV = links.at(i).linkVids;
            auto& linkE = links.at(i).linkEids;
            if (linkV.size() != 3) continue;
            std::vector<size_t> verticesToCheck(linkV.begin(), linkV.end());
            for (auto linkVid: linkV) {
                mu.AddContents(verticesToCheck, mesh.V.at(linkVid).N_Vids);
                
            }
            bool skip = false;
            for (auto vid: verticesToCheck) {
                int featureCount = 0;
                for (auto nvid: mesh.V.at(vid).N_Vids) {
                    if (mesh.V.at(nvid).type == FEATURE || mesh.V.at(nvid).isBoundary) featureCount += 1;
                }
                if (featureCount > 2) {
                    skip = true;
                    break;
                }
            }
            if (skip) continue;
            auto& v_front = mesh.V.at(linkV.front());
            auto& v_back = mesh.V.at(linkV.back());
            if (v_front.isBoundary || v_front.type == FEATURE || v_back.isBoundary || v_back.type == FEATURE) continue;
            if (v_front.N_Fids.size() != 3 || v_back.N_Fids.size() != 3) continue;
            // if (mesh.V.at(GetDiagonalV(v_front.id, GetFaceId(v_front.id, linkV.at(1)))).N_Fids.size() <= 4 || mesh.V.at(GetDiagonalV(v_back.id, GetFaceId(v_back.id, linkV.at(linkV.size()-2)))).N_Fids.size() <= 4) continue;
            auto& midV = mesh.V.at(linkV.at(1));
            if (midV.type == FEATURE || midV.isBoundary) {
                std::shared_ptr<SimplificationOperation> s = std::make_shared<SeparatrixCollapse>(mesh, mu, linkV, linkE);
                s->SetRanking();
                Op_Q.insert(s->ranking, s->GetCenterId(), s);
            }
        }
    }
    std::cout << Op_Q.size() << " boundary separatrix operations" << std::endl;
    int i = 0;
    while (!Op_Q.empty()) {
        auto& op = Op_Q.pop();
        op->PerformOperation();
        i += 1;
        // if (i > iters) break;
    }
    // std::cout << "Finished all separatrix collapse operations" << std::endl;
    // std::cout << "Smoothing mesh" << std::endl;
    // std::vector<size_t> smoothv;
    // for (auto& v: mesh.V) {
    //     if (v.N_Vids.empty() || v.N_Eids.empty() || v.N_Fids.empty()) continue;
    //     smoothv.push_back(v.id);
    // }
    // smoother.Smooth(smoothv);
}

void SemiGlobalSimplifier::SetHalfSeparatrixOperations() {
    CheckValidity();

    BaseComplexQuad bc(mesh);
    Op_Q.clear();
    Op_Q.setMaxQueueOn();
    Op_Q.setSpecialComparisonOff();
    for (auto& v: mesh.V) {
        if (!v.isSingularity) continue;
        std::vector<SingularityLink> links = TraceSingularityLinks(v, bc);
        for (int i = 0; i < links.size(); i++) {
            auto& linkV = links.at(i).linkVids;
            auto& linkE = links.at(i).linkEids;
            auto& v_front = mesh.V.at(linkV.front());
            auto& v_back = mesh.V.at(linkV.back());
                
            if (!(v_front.isBoundary ^ v_back.isBoundary)) continue;
            if (v_front.isBoundary && !v_back.isBoundary) {
                std::reverse(linkV.begin(), linkV.end());
                std::reverse(linkE.begin(), linkE.end());
                auto& v_front = mesh.V.at(linkV.front());
                auto& v_back = mesh.V.at(linkV.back());
            }
            if (!(v_front.N_Fids.size() == 3 && v_back.N_Fids.size() == 2)) continue;
            std::shared_ptr<SimplificationOperation> s = std::make_shared<SeparatrixCollapse>(mesh, mu, linkV, linkE, true);
            s->SetRanking();
            Op_Q.insert(s->ranking, s->GetCenterId(), s);
        }
    }
    std::cout << Op_Q.size() << " half separatrix operations" << std::endl;
    int i = 0;
    while (!Op_Q.empty()) {
        auto& op = Op_Q.pop();
        op->PerformOperation();
        for (auto vid: op->toUpdate) {
            std::vector<SingularityLink> links = TraceSingularityLinks(mesh.V.at(vid), bc);
            for (int i = 0; i < links.size(); i++) {
                auto& linkV = links.at(i).linkVids;
                auto& linkE = links.at(i).linkEids;
                auto& v_front = mesh.V.at(linkV.front());
                auto& v_back = mesh.V.at(linkV.back());
                
                if (!(v_front.isBoundary ^ v_back.isBoundary)) continue;
                if (v_front.isBoundary && !v_back.isBoundary) {
                    std::reverse(linkV.begin(), linkV.end());
                    std::reverse(linkE.begin(), linkE.end());
                    auto& v_front = mesh.V.at(linkV.front());
                    auto& v_back = mesh.V.at(linkV.back());
                }
                if (!(v_front.N_Fids.size() == 3 && v_back.N_Fids.size() == 2)) continue;
                std::shared_ptr<SimplificationOperation> s = std::make_shared<SeparatrixCollapse>(mesh, mu, linkV, linkE, true);
                s->SetRanking();
                Op_Q.insert(s->ranking, s->GetCenterId(), s);
            }
        }
        i += 1;
        // if (i >= iters) break;
    }
    // std::cout << "Smoothing mesh" << std::endl;
    // std::vector<size_t> smoothv;
    // for (auto& v: mesh.V) {
    //     if (v.N_Vids.empty() || v.N_Eids.empty() || v.N_Fids.empty()) continue;
    //     smoothv.push_back(v.id);
    // }
    // smoother.Smooth(smoothv);
}

void SemiGlobalSimplifier::SetChordCollapseOperations() {
    CheckValidity();

    ChordExtractor ce(mesh);
    ce.Extract();
    // std::vector<size_t> chords = ce.SelectChords();
    // std::cout << chords.size() << std::endl;
    // ce.Write(std::vector<size_t>{chords.at(0)});
    // ce.Write(chords);
    // return;

    // int maxId = -1;
    // int rank = 0;
    for (auto& chord: ce.Chords) {
        std::shared_ptr<SimplificationOperation> s = std::make_shared<ChordCollapse>(mesh, mu, ce, chord.id);
        s->PerformOperation();
    }
    // for (auto& chord: ce.Chords) {
    //     std::cout << chord.Edges.size() << std::endl;
    //     if (maxId == -1 || (chord.Edges.size() > 500 && chord.Edges.size() < 1000)) {
    //         maxId = chord.id;
    //         rank = chord.Edges.size();
    //     }
    // }
    // auto& chord = ce.Chords.at(maxId);
    // std::cout << chord.Edges.size() << std::endl;
    // int colorValue = 0;
    // int ncolors = 15;
    // std::vector<size_t> c_indices;
    // std::vector<int> colors;

    // int prevId = ce.Edges.at(chord.Edges.at(0)).Vids.at(0);
    // for (auto eid: chord.Edges) {
    //     Edge& e = ce.Edges.at(eid);
    //     c_indices.push_back(prevId);
    //     prevId = e.Vids.at(0) == prevId ? e.Vids.at(1) : e.Vids.at(0);
    //     colors.push_back(colorValue%ncolors);
    //     colorValue += 1;
    // }
    // c_indices.push_back(prevId);
    // colors.push_back(colorValue%ncolors);
    
    // std::cout << "Writing output file" << std::endl;
    // std::ofstream ofs("chord.vtk");
    // ofs << "# vtk DataFile Version 3.0\n"
    //     << "singularityLinks" << ".vtk\n"
    //     << "ASCII\n\n"
    //     << "DATASET UNSTRUCTURED_GRID\n";
    // ofs << "POINTS " << mesh.V.size() << " double\n";
    // // for (auto& v: source.V) {
    // //     if (v.type == FEATURE) c_indices.push_back(v.id);
    // // }
    // // std::vector<size_t> c_indices = {12, 296};
    // // std::cout << c_indices.size() << std::endl;
    // for (size_t i = 0; i < mesh.V.size(); i++) {
    //     ofs << std::fixed << std::setprecision(7) <<  mesh.V.at(i).x << " " <<  mesh.V.at(i).y << " " <<  mesh.V.at(i).z << "\n";
    // }
    // ofs << "CELLS " << c_indices.size() << " " << 3 * c_indices.size() << std::endl;
    // for (size_t i = 0; i < c_indices.size(); i++) {
    //     auto& e = mesh.E.at(c_indices.at(i));
    //     ofs << "2 " << e.Vids.at(0) << " " << e.Vids.at(1) << std::endl;
    // }
    // ofs << "CELL_TYPES " << c_indices.size() << "\n";
    // for (size_t i = 0; i < c_indices.size(); i++) {
    //     ofs << "3" << std::endl;
    // }

    // ofs << "CELL_DATA " << c_indices.size() << "\n";
    // ofs << "SCALARS fixed int\n";
    // ofs << "LOOKUP_TABLE default\n";
    // for (auto c: colors) {
    //     ofs << c << "\n";
    // }

    // std::cout << "Smoothing mesh" << std::endl;
    // std::vector<size_t> smoothv;
    // for (auto& v: mesh.V) {
    //     if (v.N_Vids.empty() || v.N_Eids.empty() || v.N_Fids.empty()) continue;
    //     smoothv.push_back(v.id);
    // }
    // smoother.Smooth(smoothv);
}

void SemiGlobalSimplifier::SetEdgeRotationOperations() {
    CheckValidity();

    int i = 0;
    for (auto& e: mesh.E) {
        std::shared_ptr<SimplificationOperation> s = std::make_shared<EdgeRotation>(mesh, mu, e.id, true);
        s->PerformOperation();
        i += 1;
        // if (i > iters) break;
    }
    // std::cout << "Smoothing mesh" << std::endl;
    // std::vector<size_t> smoothv;
    // for (auto& v: mesh.V) {
    //     if (v.N_Vids.empty() || v.N_Eids.empty() || v.N_Fids.empty()) continue;
    //     smoothv.push_back(v.id);
    // }
    // smoother.Smooth(smoothv);
}

void SemiGlobalSimplifier::SetVertexRotationOperations() {
    CheckValidity();

    int i = 0;
    for (auto& v: mesh.V) {
        std::shared_ptr<SimplificationOperation> s = std::make_shared<VertexRotation>(mesh, mu, v.id);
        s->PerformOperation();
        // i += 1;
        // if (i > iters) break;
    }
    // std::cout << "Smoothing mesh" << std::endl;
    // std::vector<size_t> smoothv;
    // for (auto& v: mesh.V) {
    //     if (v.N_Vids.empty() || v.N_Eids.empty() || v.N_Fids.empty()) continue;
    //     smoothv.push_back(v.id);
    // }
    // smoother.Smooth(smoothv);
}

void SemiGlobalSimplifier::SetVertexSplitOperations() {
    CheckValidity();

    int i = 0;
    for (int j = 0; j < mesh.V.size(); j++) {
        auto& v = mesh.V.at(j);
        if (v.type != FEATURE) continue;
        if (v.N_Vids.size() < 6) continue;
        std::shared_ptr<SimplificationOperation> s = std::make_shared<VertexSplit>(mesh, mu, v.id);
        s->PerformOperation();
        // i += 1;
        // if (i > iters) break;
        // break;
    }
}

void SemiGlobalSimplifier::SetEdgeCollapseOperations() {
    CheckValidity();

    Edge& e = mesh.E.at(0);
    size_t source_id = e.Vids.at(0);
    size_t targte_id = e.Vids.at(1);
    std::shared_ptr<SimplificationOperation> vr = std::make_shared<VertexRotation>(mesh, mu, targte_id);
    std::shared_ptr<SimplificationOperation> dc = std::make_shared<DiagonalCollapse>(mesh, mu, -1, targte_id, source_id);
    std::shared_ptr<SimplificationOperation> s = std::make_shared<EdgeCollapse>(mesh, mu, vr, dc);
    s->PerformOperation();
}

void SemiGlobalSimplifier::SetQuadSplitOperations() {
    CheckValidity();

    size_t startId = -1;
    for (auto& v: mesh.V) {
        bool skip = false;
        if (v.type == FEATURE || v.isBoundary || v.N_Vids.size() != 5) skip = true;
        for (auto vid: v.N_Vids) if (mesh.V.at(vid).type == FEATURE || mesh.V.at(vid).isBoundary || mesh.V.at(vid).N_Vids.size() != 4) skip = true;
        if (skip) continue;
        int d = (v.N_Vids.size() / 2) + 1;
        bool clockwise = false;
        size_t mainV = v.N_Vids.at(0);
        std::vector<size_t> verticesToSplit;
        std::vector<size_t> verticesToChange;
        size_t startE;
        for (auto eid: v.N_Eids) {
            auto& e = mesh.E.at(eid);
            if (e.Vids.at(0) == mainV || e.Vids.at(1) == mainV) {
                startE = e.id;
            }
        }
        for (int j = 0; j < d; j++) {
            auto& e = mesh.E.at(startE);
            size_t ev = e.Vids.at(0) == v.id ? e.Vids.at(1) : e.Vids.at(0);
            for (auto fid: e.N_Fids) {
                auto& f = mesh.F.at(fid);
                int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
                if (!clockwise && f.Vids.at((idx+1)%f.Vids.size()) == ev) {
                    if (j == 0 || j == d-1) {
                        verticesToSplit.push_back(ev);
                    } else {
                        verticesToChange.push_back(ev);
                    }
                    startE = mu.GetDifference(mu.GetIntersection(f.Eids, v.N_Eids), std::vector<size_t>{e.id}).at(0);
                } else if (clockwise && f.Vids.at((idx+3)%f.Vids.size()) == ev) {
                    if (j == 0 || j == d-1) {
                        verticesToSplit.push_back(ev);
                    } else {
                        verticesToChange.push_back(ev);
                    }
                    startE = mu.GetDifference(mu.GetIntersection(f.Eids, v.N_Eids), std::vector<size_t>{e.id}).at(0);
                }
            }
        }
        startId = v.id;
        break;
    }
    // return;
    if (startId != -1) {
        auto& start = mesh.V.at(startId);
        size_t threeId = -1;
        size_t fiveId = -1;
        for (auto vid: start.N_Vids) {
            std::cout << mesh.V.at(vid).N_Vids.size() << std::endl; 
            if (mesh.V.at(vid).N_Vids.size() == 5) {
                fiveId = vid;
                break;
            }
        }
        if (fiveId == -1) return;
        auto& five = mesh.V.at(fiveId);
        for (auto vid: five.N_Vids) {
            std::cout << mesh.V.at(vid).N_Vids.size() << std::endl; 
            if (mesh.V.at(vid).N_Vids.size() == 3) {
                threeId = vid;
                break;
            }
        }
        if (threeId == -1) return;
        auto& three = mesh.V.at(threeId);
        std::cout << fiveId << " " << five.N_Vids.size() << " " << threeId << " " << three.N_Vids.size() << std::endl;
        std::unique_ptr<ThreeFivePair> tfp = std::make_unique<ThreeFivePair>(mesh, mu, threeId, fiveId);
        tfp->MoveLowerLeft();
        tfp->MoveLowerLeft();
        tfp->MoveLowerLeft();
    }
    return;

    int it = 0;
    for (int i = 0; i < mesh.V.size(); i++) {
        auto& v = mesh.V.at(i);
        if (v.type != FEATURE && v.N_Vids.size() > 4) {
            // std::cout << "Quad Split" << std::endl;
            // std::cout << "i: " << i << " " << mesh.V.size() << std::endl;
            int d = (v.N_Vids.size() / 2) + 1;
            bool clockwise = false;
            size_t mainV = v.N_Vids.at(0);
            std::vector<size_t> verticesToSplit;
            std::vector<size_t> verticesToChange;
            size_t startE;
            for (auto eid: v.N_Eids) {
                auto& e = mesh.E.at(eid);
                if (e.Vids.at(0) == mainV || e.Vids.at(1) == mainV) {
                    startE = e.id;
                }
            }
            for (int j = 0; j < d; j++) {
                auto& e = mesh.E.at(startE);
                size_t ev = e.Vids.at(0) == v.id ? e.Vids.at(1) : e.Vids.at(0);
                for (auto fid: e.N_Fids) {
                    auto& f = mesh.F.at(fid);
                    int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
                    if (!clockwise && f.Vids.at((idx+1)%f.Vids.size()) == ev) {
                        if (j == 0 || j == d-1) {
                            verticesToSplit.push_back(ev);
                        } else {
                            verticesToChange.push_back(ev);
                        }
                        startE = mu.GetDifference(mu.GetIntersection(f.Eids, v.N_Eids), std::vector<size_t>{e.id}).at(0);
                    } else if (clockwise && f.Vids.at((idx+3)%f.Vids.size()) == ev) {
                        if (j == 0 || j == d-1) {
                            verticesToSplit.push_back(ev);
                        } else {
                            verticesToChange.push_back(ev);
                        }
                        startE = mu.GetDifference(mu.GetIntersection(f.Eids, v.N_Eids), std::vector<size_t>{e.id}).at(0);
                    }
                }
            }
            // if (verticesToSplit.size() < 2) continue;
            std::shared_ptr<SimplificationOperation> qs = std::make_shared<QuadSplit>(mesh, mu, v.id, verticesToSplit, verticesToChange);
            qs->PerformOperation();
            it += 1;
            // if (it >= iters) break;
        }    
    }
}

void SemiGlobalSimplifier::ResolveSingularityPairs() {
    CheckValidity();

    for (int i = 0; i < mesh.V.size(); i++) {
        auto& v = mesh.V.at(i);
        if (v.type != FEATURE && v.N_Vids.size() == 3) {
            CheckAndResolveThreeFivePair(v.id);
        } else if (v.type != FEATURE && v.N_Vids.size() == 5) {
            CheckAndResolveFiveThreePair(v.id);
        }
        // if (i >= iters) break;
    }
}

void SemiGlobalSimplifier::CheckAndResolveThreeFivePair(size_t vid) {
    auto& v = mesh.V.at(vid);
    bool clockwise = false;
    bool isPair = false;
    for (auto fid: v.N_Fids) {
        auto& f = mesh.F.at(fid);
        int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), vid));
        if (mesh.V.at(f.Vids.at((idx+2)%f.Vids.size())).N_Vids.size() == 5 && mesh.V.at(f.Vids.at((idx+3)%f.Vids.size())).N_Vids.size() == 5 && mesh.V.at(f.Vids.at((idx+1)%f.Vids.size())).N_Vids.size() == 4) {
            isPair = true;
        } else if (mesh.V.at(f.Vids.at((idx+1)%f.Vids.size())).N_Vids.size() == 5 && mesh.V.at(f.Vids.at((idx+2)%f.Vids.size())).N_Vids.size() == 5 && mesh.V.at(f.Vids.at((idx+3)%f.Vids.size())).N_Vids.size() == 4) {
            isPair = true;
            clockwise = true;
        } else if (mesh.V.at(f.Vids.at((idx+1)%f.Vids.size())).N_Vids.size() == 5 && mesh.V.at(f.Vids.at((idx+2)%f.Vids.size())).N_Vids.size() == 4 && mesh.V.at(f.Vids.at((idx+3)%f.Vids.size())).N_Vids.size() == 5) {
            std::shared_ptr<SimplificationOperation> dc1 = std::make_shared<DiagonalCollapse>(mesh, mu, f.id, idx, (idx+2)%f.Vids.size());
            dc1->PerformOperation();
            break;
        }
        if (isPair) {
            for (auto eid: f.Eids) {
                auto& e = mesh.E.at(eid);
                if (mesh.V.at(e.Vids.at(0)).N_Vids.size() == 5 && mesh.V.at(e.Vids.at(1)).N_Vids.size() == 5) {
                    std::shared_ptr<SimplificationOperation> s = std::make_shared<EdgeRotation>(mesh, mu, e.id, clockwise);
                    s->PerformOperation();
                }
            }
            break;
        }
    }
}


void SemiGlobalSimplifier::CheckAndResolveFiveThreePair(size_t vid) {
    auto& v = mesh.V.at(vid);
    for (auto eid: v.N_Eids) {
        auto& e = mesh.E.at(eid);
        bool isPair = false;
        size_t neid = e.Vids.at(0) == vid ? e.Vids.at(1) : e.Vids.at(0);
        if (mesh.V.at(neid).N_Vids.size() != 4) continue;
        auto& f1 = mesh.F.at(e.N_Fids.at(0));
        auto& f2 = mesh.F.at(e.N_Fids.at(1));
        int idx1 = std::distance(f1.Vids.begin(), std::find(f1.Vids.begin(), f1.Vids.end(), vid));
        int idx2 = std::distance(f2.Vids.begin(), std::find(f2.Vids.begin(), f2.Vids.end(), vid));
        if ((mesh.V.at(f1.Vids.at((idx1+1)%f1.Vids.size())).N_Vids.size() == 3 && mesh.V.at(f2.Vids.at((idx2+2)%f2.Vids.size())).N_Vids.size() == 3)
        || (mesh.V.at(f1.Vids.at((idx1+2)%f1.Vids.size())).N_Vids.size() == 3 && mesh.V.at(f2.Vids.at((idx2+1)%f2.Vids.size())).N_Vids.size() == 3)) {
            std::shared_ptr<SimplificationOperation> s = std::make_shared<EdgeRotation>(mesh, mu, e.id, false);
            s->PerformOperation();
            isPair = true;
        } else if (mesh.V.at(f1.Vids.at((idx1+1)%f1.Vids.size())).N_Vids.size() == 3 && mesh.V.at(f2.Vids.at((idx2+3)%f2.Vids.size())).N_Vids.size() == 3) {
            std::vector<size_t> verticesToSplit = {f1.Vids.at((idx1+1)%f1.Vids.size()), f2.Vids.at((idx2+3)%f2.Vids.size())};
            std::vector<size_t> verticesToChange = {neid};
            std::shared_ptr<SimplificationOperation> qs = std::make_shared<QuadSplit>(mesh, mu, v.id, verticesToSplit, verticesToChange);
            qs->PerformOperation();
            isPair = true;
        } else if (mesh.V.at(f1.Vids.at((idx1+3)%f1.Vids.size())).N_Vids.size() == 3 && mesh.V.at(f2.Vids.at((idx2+1)%f2.Vids.size())).N_Vids.size() == 3) {
            std::vector<size_t> verticesToSplit = {f1.Vids.at((idx1+3)%f1.Vids.size()), f2.Vids.at((idx2+1)%f2.Vids.size())};
            std::vector<size_t> verticesToChange = {neid};
            std::shared_ptr<SimplificationOperation> qs = std::make_shared<QuadSplit>(mesh, mu, v.id, verticesToSplit, verticesToChange);
            qs->PerformOperation();
            isPair = true;
        } else if ((mesh.V.at(f1.Vids.at((idx1+2)%f1.Vids.size())).N_Vids.size() == 3 && mesh.V.at(f2.Vids.at((idx2+3)%f2.Vids.size())).N_Vids.size() == 3)
        || (mesh.V.at(f1.Vids.at((idx1+3)%f1.Vids.size())).N_Vids.size() == 3 && mesh.V.at(f2.Vids.at((idx2+2)%f2.Vids.size())).N_Vids.size() == 3)) {
            std::shared_ptr<SimplificationOperation> s = std::make_shared<EdgeRotation>(mesh, mu, e.id, true);
            s->PerformOperation();
            isPair = true;
        }
        if (isPair) break;
    }
}

void SemiGlobalSimplifier::GetSingularityPairs() {
    BaseComplexQuad bc(mesh);
    std::vector<SingularityLink> singularityLinks;
    // std::unordered_map<int, std::vector<int>> singularityMap;
    std::vector<std::vector<int>> singularityMap(mesh.V.size(), std::vector<int>{});
    clock_t start = clock();
    // auto cmp = [](SingularityLink left, SingularityLink right) {return left.a + left.b < right.a + right.b;};
    // std::vector<std::multiset<SingularityLink, decltype(cmp)>> singularityMap(mesh.V.size(), std::multiset<SingularityLink, decltype(cmp)>(cmp));
    std::mutex mtx;
    std::vector<size_t> Singularities;
    
    // std::priority_queue<SingularityLink, std::vector<SingularityLink>, decltype(cmp)> q(cmp);

    /*PARALLEL_FOR_BEGIN(mesh.V.size()) {
        auto& v = mesh.V.at(i);
        if (v.type == FEATURE || v.isBoundary || !v.isSingularity || v.N_Fids.empty()) continue;
        {
            std::lock_guard<std::mutex> lock(mtx);
            Singularities.push_back(v.id);
        }
    } PARALLEL_FOR_END();*/

    // int sid = 0;
    PARALLEL_FOR_BEGIN(mesh.V.size()) {
        // auto& v = mesh.V.at(Singularities.at(i));
        auto& v = mesh.V.at(i);
    // for (auto& v: mesh.V) {
        std::vector<SingularityLink> links = TraceSingularityLinks(v, bc);
        for (auto& l: links) {
            if (mesh.V.at(l.frontId).N_Vids.size() == 4 || mesh.V.at(l.backId).N_Vids.size() == 4) continue;
            {
                std::lock_guard<std::mutex> lock(mtx);
                // l.id = sid++;
                l.id = singularityLinks.size();
                l.a = l.linkEids.size();
                l.b = 0;
                singularityLinks.push_back(l);

                // singularityMap.at(v.id).insert(l);
                singularityMap.at(v.id).push_back(l.id);
            }
        }
    // }
    } PARALLEL_FOR_END();

    /*std::vector<SingularityLink> newLinks;
    PARALLEL_FOR_BEGIN(singularityLinks.size()) {
        auto& l = singularityLinks.at(i);
        if (doesCrossBoundary(l.linkVids, true)) continue;
        SingularityLink newL;
        for (int j = i+1; j < singularityLinks.size(); j++) {
            auto& l2 = singularityLinks.at(j);
            if (doesCrossBoundary(l2.linkVids, true)) continue;
            if (l.frontId == l2.frontId || l.frontId == l2.backId || l.backId == l2.frontId || l.backId == l2.backId) continue;
            std::vector<size_t> a = mu.GetIntersectionParallel(l.linkVids, l2.linkVids);
            if (a.size() != 1) continue;
            int vid = a.at(0);
            std::vector<size_t> b;
            for (auto eid: l.linkEids) {
                auto& e = mesh.E.at(eid);
                b.push_back(eid);
                if (e.Vids.at(1) == vid || e.Vids.at(0) == vid) break;
            }
            std::vector<size_t> c;
            int it = 0;
            for (auto eid: l2.linkEids) {
                auto& e = mesh.E.at(eid);
                c.push_back(eid);
                it += 1;
                if (e.Vids.at(1) == vid || e.Vids.at(0) == vid) break;
            }
            std::vector<size_t> c1(l2.linkEids.begin()+it, l2.linkEids.end());
            std::reverse(c.begin(), c.end());

            newL.a = b.size();
            bool isFirstHalfInvalid = doesCrossBoundary(c, false);
            bool isSecondHalfInvalid = doesCrossBoundary(c1, false);
            if (isFirstHalfInvalid && isSecondHalfInvalid) continue;
            if (c.size() <= c1.size()) {
                if (!isFirstHalfInvalid) {
                    newL.b = c.size();
                    b.insert(b.end(), c.begin(), c.end());
                } else if (!isSecondHalfInvalid) {
                    newL.b = c1.size();
                    b.insert(b.end(), c1.begin(), c1.end());
                }
            } else if (c1.size() < c.size()) {
                if (!isSecondHalfInvalid) {
                    newL.b = c1.size();
                    b.insert(b.end(), c1.begin(), c1.end());
                } else if (!isFirstHalfInvalid) {
                    newL.b = c.size();
                    b.insert(b.end(), c.begin(), c.end());
                }
            }

            std::vector<size_t> d;
            int prevId = l.frontId;
            for (auto el: b) {
                auto& e = mesh.E.at(el);
                d.push_back(prevId);
                prevId = e.Vids.at(0) == prevId ? e.Vids.at(1) : e.Vids.at(0);
            }
            d.push_back(prevId);
            newL.linkVids = d;
            newL.linkEids = b;
            newL.frontId = newL.linkVids.front();
            newL.backId = newL.linkVids.back();
            
            if (newL.linkVids.empty() || newL.linkEids.empty()) continue;
            SingularityLink newL2;
            newL2.a = newL.b;
            newL2.b = newL.a;
            std::vector<size_t> e = d;
            std::vector<size_t> f = b;
            std::reverse(e.begin(), e.end());
            std::reverse(f.begin(), f.end());
            newL2.linkVids = e;
            newL2.linkEids = f;
            newL2.frontId = newL2.linkVids.front();
            newL2.backId = newL2.linkVids.back();
            
            {
                std::lock_guard<std::mutex> lock(mtx);
                newL.id = sid++;
                newLinks.push_back(newL);
                newL2.id = sid++;
                newLinks.push_back(newL2);
                // singularityMap.at(v.id).insert(newL);
                // singularityMap.at(v.id).push_back(newL.id);
            }
        }
    } PARALLEL_FOR_END();*/

    /*PARALLEL_FOR_BEGIN(Singularities.size()) {
        auto& v = mesh.V.at(Singularities.at(i));
        for (auto vid: Singularities) {
            if (v.id == vid) continue;
            auto& v2 = mesh.V.at(vid);
            for (auto lid: singularityMap.at(v.id)) {
                auto& l = singularityLinks.at(lid);
                if (doesCrossBoundary(l.linkVids, true)) continue;
                int rank1 = -1; 
                int rank2 = -1;
                SingularityLink newL;
                for (auto lid2: singularityMap.at(v2.id)) {
                    auto& l2 = singularityLinks.at(lid2);
                    if (l.frontId == l2.frontId || l.frontId == l2.backId || l.backId == l2.frontId || l.backId == l2.backId) continue;
                    std::vector<size_t> a = mu.GetIntersectionParallel(l.linkVids, l2.linkVids);
                    if (a.size() != 1) continue;
                    int vid = a.at(0);
                    std::vector<size_t> b;
                    for (auto eid: l.linkEids) {
                        auto& e = mesh.E.at(eid);
                        b.push_back(eid);
                        if (e.Vids.at(1) == vid || e.Vids.at(0) == vid) break;
                    }
                    std::vector<size_t> c;
                    int it = 0;
                    for (auto eid: l2.linkEids) {
                        auto& e = mesh.E.at(eid);
                        c.push_back(eid);
                        it += 1;
                        if (e.Vids.at(1) == vid || e.Vids.at(0) == vid) break;
                    }
                    std::vector<size_t> c1(l2.linkEids.begin()+it, l2.linkEids.end());
                    std::reverse(c.begin(), c.end());

                    newL.a = b.size();
                    bool isFirstHalfInvalid = doesCrossBoundary(c, false);
                    bool isSecondHalfInvalid = doesCrossBoundary(c1, false);
                    if (isFirstHalfInvalid && isSecondHalfInvalid) continue;
                    if (c.size() <= c1.size()) {
                        if (!isFirstHalfInvalid) {
                            newL.b = c.size();
                            b.insert(b.end(), c.begin(), c.end());
                        } else if (!isSecondHalfInvalid) {
                            newL.b = c1.size();
                            b.insert(b.end(), c1.begin(), c1.end());
                        }
                    } else if (c1.size() < c.size()) {
                        if (!isSecondHalfInvalid) {
                            newL.b = c1.size();
                            b.insert(b.end(), c1.begin(), c1.end());
                        } else if (!isFirstHalfInvalid) {
                            newL.b = c.size();
                            b.insert(b.end(), c.begin(), c.end());
                        }
                    }

                    std::vector<size_t> d;
                    int prevId = l.frontId;
                    for (auto el: b) {
                        auto& e = mesh.E.at(el);
                        d.push_back(prevId);
                        prevId = e.Vids.at(0) == prevId ? e.Vids.at(1) : e.Vids.at(0);
                    }
                    d.push_back(prevId);
                    newL.linkVids = d;
                    newL.linkEids = b;
                    newL.frontId = newL.linkVids.front();
                    newL.backId = newL.linkVids.back();
                    if (newL.linkVids.empty() || newL.linkEids.empty()) continue;
                    {
                        std::lock_guard<std::mutex> lock(mtx);
                        newL.id = sid++;
                        newLinks.push_back(newL);
                        // singularityMap.at(v.id).insert(newL);
                        // singularityMap.at(v.id).push_back(newL.id);
                    }
                }
            }
        }
        
    } PARALLEL_FOR_END();

    PARALLEL_FOR_BEGIN(newLinks.size()) {
        auto& l = newLinks.at(i);
        {
            std::lock_guard<std::mutex> lock(mtx);
            l.id = singularityLinks.size();
            singularityLinks.push_back(l);
            // singularityMap.at(l.frontId).insert(l);
            singularityMap.at(l.frontId).push_back(l.id);
        }
    } PARALLEL_FOR_END();*/

    // auto cmpPair = [](SingularityGroup left, SingularityGroup right) {return left.rank > right.rank;};
    // std::priority_queue<SingularityGroup, std::vector<SingularityGroup>, decltype(cmpPair)> q(cmpPair);

    /*PARALLEL_FOR_BEGIN(Singularities.size()) {
        auto& v = mesh.V.at(Singularities.at(i));
        for (int j = 0; j < singularityMap.at(v.id).size(); j++) {
            std::multiset<SingularityLink, decltype(cmp)>::iterator it = singularityMap.at(v.id).begin();
            std::advance(it, j);
            auto& l = *it;
            size_t vid = l.frontId == v.id ? l.backId : l.frontId;
            auto& v2 = mesh.V.at(vid);
            if (v.N_Vids.size() > 5 || v2.N_Vids.size() > 5) continue;
            if (doesCrossBoundary(l.linkVids, true)) continue;
            bool skipSameSingularityType = v.N_Vids.size() == v2.N_Vids.size() ? true : false;    
            for (int k = j + 1; k < singularityMap.at(v.id).size(); k++) {
                it = singularityMap.at(v.id).begin();
                std::advance(it, k);
                auto& l2 = *it;
                size_t vid2 = l2.frontId == v.id ? l2.backId : l2.frontId;
                auto& v3 = mesh.V.at(vid2);
                if (v.N_Vids.size() > 5 || v3.N_Vids.size() > 5) continue;
                if (v.N_Vids.size() == v3.N_Vids.size() && skipSameSingularityType) continue;
                if (doesCrossBoundary(l2.linkVids, true)) continue;
                
                SingularityGroup sg;
                sg.l1 = l;
                sg.rank += l.a + l.b;
                sg.l2 = l2;
                sg.rank += l2.a + l2.b;
                {
                    std::lock_guard<std::mutex> lock(mtx);
                    q.push(sg);
                }
            }
        }
        for (auto& l: singularityMap.at(v.id)) {
            size_t vid = l.frontId == v.id ? l.backId : l.frontId;
            auto& v2 = mesh.V.at(vid);
            if (v2.N_Vids.size() < 3 || v2.N_Vids.size() > 5) continue;
            if (doesCrossBoundary(l.linkVids, true)) continue;
            if (v.N_Vids.size() == v2.N_Vids.size()) {
                if (skipSameSingularityType) {
                    continue;
                } else {
                    skipSameSingularityType = true;
                }
            }
            if (it == 0) {
                sg.l1 = l;
                sg.rank += l.a + l.b;
                it += 1;
            } else if (it == 1) {
                if ((sg.l1.frontId == l.frontId && sg.l1.backId == l.backId) || (sg.l1.backId == l.frontId && sg.l1.frontId == l.backId)) continue;
                sg.l2 = l;
                sg.rank += l.a + l.b;
                {
                    std::lock_guard<std::mutex> lock(mtx);
                    q.push(sg);
                    break;
                }
            }
        }
    } PARALLEL_FOR_END();*/

    int offset = 0;
    std::vector<SingularityLink> newLinks;
    PARALLEL_FOR_BEGIN(singularityMap.size()) {
    // for (int i = 0; i < singularityMap.size(); i++) {
        int n = singularityMap.at(i).size();
        for (int j = 0; j < n; j++) {
            auto& l = singularityLinks.at(singularityMap.at(i).at(j));
            if (doesCrossBoundary(l.linkVids, true)) continue;
            int rank1 = -1; 
            int rank2 = -1;
            SingularityLink newL;
            for (auto& l2: singularityLinks) {
                if (l.id == l2.id) continue;
                if (l.frontId == l2.frontId || l.frontId == l2.backId || l.backId == l2.frontId || l.backId == l2.backId) continue;
                std::vector<size_t> a = mu.GetIntersectionParallel(l.linkVids, l2.linkVids);
                if (a.size() != 1) continue;
                int vid = a.at(0);
                std::vector<size_t> b;
                for (auto eid: l.linkEids) {
                    auto& e = mesh.E.at(eid);
                    b.push_back(eid);
                    if (e.Vids.at(1) == vid || e.Vids.at(0) == vid) break;
                }
                std::vector<size_t> c;
                int it = 0;
                for (auto eid: l2.linkEids) {
                    auto& e = mesh.E.at(eid);
                    c.push_back(eid);
                    it += 1;
                    if (e.Vids.at(1) == vid || e.Vids.at(0) == vid) break;
                }
                std::vector<size_t> c1(l2.linkEids.begin()+it, l2.linkEids.end());
                std::reverse(c.begin(), c.end());

                newL.a = b.size();
                bool isFirstHalfInvalid = doesCrossBoundary(c, false);
                bool isSecondHalfInvalid = doesCrossBoundary(c1, false);
                if (isFirstHalfInvalid && isSecondHalfInvalid) continue;
                if (c.size() <= c1.size()) {
                    if (!isFirstHalfInvalid) {
                        newL.b = c.size();
                        b.insert(b.end(), c.begin(), c.end());
                    } else if (!isSecondHalfInvalid) {
                        newL.b = c1.size();
                        b.insert(b.end(), c1.begin(), c1.end());
                    }
                } else if (c1.size() < c.size()) {
                    if (!isSecondHalfInvalid) {
                        newL.b = c1.size();
                        b.insert(b.end(), c1.begin(), c1.end());
                    } else if (!isFirstHalfInvalid) {
                        newL.b = c.size();
                        b.insert(b.end(), c.begin(), c.end());
                    }
                }

                std::vector<size_t> d;
                int prevId = l.frontId;
                for (auto el: b) {
                    auto& e = mesh.E.at(el);
                    d.push_back(prevId);
                    prevId = e.Vids.at(0) == prevId ? e.Vids.at(1) : e.Vids.at(0);
                }
                d.push_back(prevId);
                newL.linkVids = d;
                newL.linkEids = b;
                newL.frontId = newL.linkVids.front();
                newL.backId = newL.linkVids.back();
                if (newL.linkVids.empty() || newL.linkEids.empty()) continue;
                {
                    std::lock_guard<std::mutex> lock(mtx);
                    newL.id = singularityLinks.size() + offset;
                    offset += 1;
                    newLinks.push_back(newL);
                    singularityMap.at(i).push_back(newL.id);
                }
            }
        }
    // }
    } PARALLEL_FOR_END();
    singularityLinks.insert(singularityLinks.end(), newLinks.begin(), newLinks.end());
    
    auto cmp = [](SingularityGroup left, SingularityGroup right) {return left.rank > right.rank;};
    std::priority_queue<SingularityGroup, std::vector<SingularityGroup>, decltype(cmp)> q(cmp);
    // std::multiset<SingularityGroup, decltype(cmp)> q(cmp);

    // std::vector<int> ThreeFivePairs;
    // PARALLEL_FOR_BEGIN(Singularities.size()) {
    //     auto& v = mesh.V.at(Singularities.at(i));
    // // for (auto& v: mesh.V) {
    //     bool skipSameSingularityType = false;
    //     int it = 0;
    //     for (auto lid: singularityMap.at(v.id)) {
    //         auto& l = singularityLinks.at(lid);
    //     }
    // // }
    // } PARALLEL_FOR_END();

    
    PARALLEL_FOR_BEGIN(singularityLinks.size()) {
        auto& l = singularityLinks.at(i);
    // for (auto& l: singularityLinks) {
        auto& v1 = mesh.V.at(l.frontId);
        auto& v2 = mesh.V.at(l.backId);
        if ((v1.N_Vids.size() == 3 && v2.N_Vids.size() == 5) || (v1.N_Vids.size() == 5 && v2.N_Vids.size() == 3)) {
            if (doesCrossBoundary(l.linkVids, true)) continue;
            std::vector<int> linksToCheck(singularityMap.at(l.frontId).begin(), singularityMap.at(l.frontId).end());
            linksToCheck.insert(linksToCheck.end(), singularityMap.at(l.backId).begin(), singularityMap.at(l.backId).end());
            int minLinkId = -1;
            int rank = -1;
            for (auto lid: linksToCheck) {
                auto& l2 = singularityLinks.at(lid);
                if ((l.frontId == l2.frontId && l.backId == l2.backId) || (l.frontId == l2.backId && l.backId == l2.frontId)) continue;
                if (mesh.V.at(l2.frontId).N_Vids.size() > 5 || mesh.V.at(l2.backId).N_Vids.size() > 5) continue;
                if (doesCrossBoundary(l2.linkVids, true)) continue;
                if (!mu.GetIntersection(l.linkEids, l2.linkEids).empty()) continue;
                // {
                //     std::lock_guard<std::mutex> lock(mtx);
                //     q.push(l);
                // }
                // SingularityGroup sg;
                // sg.l1 = l;
                // sg.l2 = l2;
                // sg.rank = l.a + l.b + l2.a + l2.b;
                // {
                //     std::lock_guard<std::mutex> lock(mtx);
                //     // q.insert(sg);
                //     q.push(sg);
                // }
                if (rank == -1 || rank > (l2.a + l2.b)) {
                    minLinkId = lid;
                    rank = l2.a + l2.b;
                }
                
            }
            if (minLinkId > -1) {
                auto& l2 = singularityLinks.at(minLinkId);
                SingularityGroup sg;
                sg.l1 = l;
                sg.l2 = l2;
                sg.rank = l.a + l.b + l2.a + l2.b;
                {
                    std::lock_guard<std::mutex> lock(mtx);
                    q.push(sg);
                }
            }
        }
    // }
    } PARALLEL_FOR_END();
    std::cout << "total singularity links: " << q.size() << std::endl;

    /*int qSize = q.size();
    // int qSize = ThreeFivePairs.size();
    std::cout << "Identified " << qSize << " groups" << std::endl;*/

    std::vector<bool> selected(mesh.V.size(), false);
    // std::vector<SingularityGroup> finalGroups;
    std::priority_queue<SingularityGroup, std::vector<SingularityGroup>, decltype(cmp)> finalGroups(cmp);
    while (!q.empty()) {
        auto& sg = q.top();
        q.pop();
        auto& l1 = sg.l1;
        auto& l2 = sg.l2;
        if (selected.at(l1.frontId) || selected.at(l1.backId) || selected.at(l2.frontId) || selected.at(l2.backId)) continue;
        selected.at(l1.frontId) = true;
        selected.at(l1.backId) = true;
        selected.at(l2.frontId) = true;
        selected.at(l2.backId) = true;
        finalGroups.push(sg);
    }
    /*while (!q.empty()) {
        auto& l = q.top();
        q.pop();
        int minLinkId = -1;
        int rank = -1;
        if (selected.at(l.frontId) || selected.at(l.backId)) {
            continue;
        }
        for (auto lid: singularityMap[l.frontId]) {
            auto& l2 = singularityLinks.at(lid);
            if ((l.frontId == l2.frontId && l.backId == l2.backId) || (l.frontId == l2.backId && l.backId == l2.frontId)) continue;
            if (selected.at(l2.frontId) || selected.at(l2.backId)) continue;
            if (mesh.V.at(l2.frontId).N_Vids.size() > 5 || mesh.V.at(l2.backId).N_Vids.size() > 5) continue;
            if (doesCrossBoundary(l2.linkVids, true)) continue;
            if (rank == -1 || rank > (l2.a + l2.b)) {
                minLinkId = lid;
                rank = l2.a + l2.b;
            }
        }
        for (auto lid: singularityMap[l.backId]) {
            auto& l2 = singularityLinks.at(lid);
            if ((l.frontId == l2.frontId && l.backId == l2.backId) || (l.frontId == l2.backId && l.backId == l2.frontId)) continue;
            if (selected.at(l2.frontId) || selected.at(l2.backId)) continue;
            if (mesh.V.at(l2.frontId).N_Vids.size() > 5 || mesh.V.at(l2.backId).N_Vids.size() > 5) continue;
            if (doesCrossBoundary(l2.linkVids, true)) continue;
            if (rank == -1 || rank > (l2.a + l2.b)) {
                minLinkId = lid;
                rank = l2.a + l2.b;
            }
        }
        if (minLinkId > -1) {
            auto& l2 = singularityLinks.at(minLinkId);
            selected.at(l.frontId) = true;
            selected.at(l.backId) = true;
            selected.at(l2.frontId) = true;
            selected.at(l2.backId) = true;
            linkPairs.push_back(std::make_pair(l.id, l2.id));
        }
    }*/

    // std::cout << finalGroups.size() << " link groups" << std::endl;

    int numSing = 0;
    for (auto& v: mesh.V) {
        if (v.isSingularity && v.N_Vids.size() > 0 && !selected.at(v.id)) {
            numSing += 1;
        }
    }
        std::cout << numSing << " non selected singularities" << std::endl;

    clock_t end = clock();
    double time = (double) (end-start) / CLOCKS_PER_SEC;
    std::cout << "execution took: " << time << " seconds" << std::endl;
    // std::cout << "identified singularities: " << singularityMap.size() << " identified links: " << singularityLinks.size() << std::endl;

    Smooth();

    int colorValue = 0;
    int ncolors = 15;
    std::vector<size_t> c_indices;
    std::vector<int> colors;
    // int it = 0;
    // for (int i = 0; i < singularityMap.size(); i++) {
    //     auto m = singularityMap.at(i);
    //     if (!m.empty()) it += 1;
    //     if (it == iters) {
    //         for (auto lid: m) {
    //             auto& l = singularityLinks.at(lid);
    //             std::vector<int> a(l.linkEids.size(), (colorValue%ncolors));
    //             c_indices.insert(c_indices.end(), l.linkEids.begin(), l.linkEids.end());
    //             colors.insert(colors.end(), a.begin(), a.end());
    //             colorValue += 1;
    //         }
    //         break;
    //     }
    // }
    // for (auto p: linkPairs) {
    //     auto& l1 = singularityLinks.at(std::get<0>(p));
    //     auto& l2 = singularityLinks.at(std::get<1>(p));
    //     std::vector<int> a(l1.linkEids.size(), (colorValue%ncolors));
    //     std::vector<int> b(l2.linkEids.size(), (colorValue%ncolors));
    //     c_indices.insert(c_indices.end(), l1.linkEids.begin(), l1.linkEids.end());
    //     c_indices.insert(c_indices.end(), l2.linkEids.begin(), l2.linkEids.end());
    //     colors.insert(colors.end(), a.begin(), a.end());
    //     colors.insert(colors.end(), b.begin(), b.end());
    //     colorValue += 1;
    // }
    // for (auto& sg: finalGroups) {
    //     std::cout << sg.rank << std::endl;
    //     auto& l1 = sg.l1;
    //     auto& l2 = sg.l2;
    //     std::vector<int> a(l1.linkEids.size(), (colorValue%ncolors));
    //     std::vector<int> b(l2.linkEids.size(), (colorValue%ncolors));
    //     c_indices.insert(c_indices.end(), l1.linkEids.begin(), l1.linkEids.end());
    //     c_indices.insert(c_indices.end(), l2.linkEids.begin(), l2.linkEids.end());
    //     colors.insert(colors.end(), a.begin(), a.end());
    //     colors.insert(colors.end(), b.begin(), b.end());
    //     colorValue += 1;
    // }
    while (!finalGroups.empty()) {
        auto& sg = finalGroups.top();
        finalGroups.pop();
        auto& l1 = sg.l1;
        auto& l2 = sg.l2;
        if (l1.linkVids.empty() || l2.linkVids.empty() || l1.linkEids.empty() || l2.linkEids.empty()) continue;
        size_t toMoveId = l1.frontId == l2.frontId || l1.frontId == l2.backId ? l1.frontId : l1.backId;
        size_t sourceId = l1.frontId == toMoveId ? l1.backId : l1.frontId;
        size_t secondaryId = l2.frontId == toMoveId ? l2.backId : l2.frontId;
        std::vector<size_t> mainPath = l1.linkVids;
        std::vector<size_t> secondaryPath = l2.linkVids;
        if (sourceId == mainPath.at(0)) {
            std::reverse(mainPath.begin(), mainPath.end());
        }
        if (secondaryId == secondaryPath.at(0)) {
            std::reverse(secondaryPath.begin(), secondaryPath.end());
        }

        if (mainPath.size() == 2) {
            std::cout << "rank: " << sg.rank << std::endl;
            std::cout << "l1 front: " << l1.frontId << " " << mesh.V.at(l1.frontId).N_Vids.size() << std::endl;
            std::cout << "l1 back: " << l1.backId << " " << mesh.V.at(l1.backId).N_Vids.size() << std::endl;
            std::cout << "l2 front: " << l2.frontId << " " << mesh.V.at(l2.frontId).N_Vids.size() << std::endl;
            std::cout << "l2 back: " << l2.backId << " " << mesh.V.at(l2.backId).N_Vids.size() << std::endl;
            auto& toMove = mesh.V.at(toMoveId);
            auto& source = mesh.V.at(sourceId);
            auto& secondary = mesh.V.at(secondaryId);
            std::cout << "toMove: " << toMove.id << std::endl;
            std::cout << "source: " << source.id << std::endl;
            std::cout << "secondary: " << secondary.id << std::endl;
            std::cout << "mainPath, ";
            for (auto id: mainPath) {
                std::cout << id << " ";
            }
            std::cout << std::endl;
            std::cout << "secondaryPath, ";
            for (auto id: secondaryPath) {
                std::cout << id << " ";
            }
            std::cout << std::endl;
            size_t threeId = toMove.N_Vids.size() == 3 ? toMoveId : sourceId;
            size_t fiveId = toMove.N_Vids.size() == 5 ? toMoveId : sourceId;                
            std::unique_ptr<ThreeFivePair> tfp = std::make_unique<ThreeFivePair>(mesh, mu, threeId, fiveId);
            for (int i = 1; i < secondaryPath.size(); i++) {
                size_t dest = secondaryPath.at(i);
                std::cout << "moving 3-5 pair to " << dest << std::endl;
                tfp->Move(dest);
            }
            std::vector<int> a(l1.linkEids.size(), (colorValue%ncolors));
            std::vector<int> b(l2.linkEids.size(), (colorValue%ncolors));
            c_indices.insert(c_indices.end(), l1.linkEids.begin(), l1.linkEids.end());
            c_indices.insert(c_indices.end(), l2.linkEids.begin(), l2.linkEids.end());
            colors.insert(colors.end(), a.begin(), a.end());
            colors.insert(colors.end(), b.begin(), b.end());
            colorValue += 1;    
        }
    
    }
    // while (!finalGroups.empty()) {
    //     auto& sg = finalGroups.top();
    //     finalGroups.pop();
    //     auto& l1 = sg.l1;
    //     auto& l2 = sg.l2;
    //     std::cout << "l1: " << mesh.V.at(l1.frontId).N_Vids.size() << " " << mesh.V.at(l1.backId).N_Vids.size() << std::endl;
    //     std::cout << "l2: " << mesh.V.at(l2.frontId).N_Vids.size() << " " << mesh.V.at(l2.backId).N_Vids.size() << std::endl;
    //     std::cout << "***********************" << std::endl;
    //     std::vector<int> a(l1.linkEids.size(), (colorValue%ncolors));
    //     colorValue += 7;
    //     std::vector<int> b(l2.linkEids.size(), (colorValue%ncolors));
    //     c_indices.insert(c_indices.end(), l1.linkEids.begin(), l1.linkEids.end());
    //     c_indices.insert(c_indices.end(), l2.linkEids.begin(), l2.linkEids.end());
    //     colors.insert(colors.end(), a.begin(), a.end());
    //     colors.insert(colors.end(), b.begin(), b.end());
    //     colorValue += 1;
    // }
    // for (auto vid: Singularities) {
    //     for (auto& l: singularityMap.at(vid)) {
    //         std::vector<int> a(l.linkEids.size(), (colorValue%ncolors));
    //         c_indices.insert(c_indices.end(), l.linkEids.begin(), l.linkEids.end());
    //         colors.insert(colors.end(), a.begin(), a.end());
    //         colorValue += 1;
    //         std::cout << l.a + l.b << std::endl;
    //     }
    //     break;
    // }

    // while (!q.empty()) {
    //     auto& l = q.top();
    //     q.pop();
    //     c_indices.insert(c_indices.end(), l.linkEids.begin(), l.linkEids.end());
    //     std::vector<int> a(l.linkEids.size(), (colorValue%ncolors));
    //     colors.insert(colors.end(), a.begin(), a.end());
    //     colorValue += 1;
    // }
    // while (!q.empty()) {
    //     auto& sg = q.top();
    //     std::cout << sg.rank << std::endl;
    //     q.pop();
    //     auto& l1 = singularityLinks.at(sg.l1_id);
    //     auto& l2 = singularityLinks.at(sg.l2_id);
    //     std::vector<int> a(l1.linkEids.size(), (colorValue%ncolors));
    //     std::vector<int> b(l2.linkEids.size(), (colorValue%ncolors));
    //     c_indices.insert(c_indices.end(), l1.linkEids.begin(), l1.linkEids.end());
    //     c_indices.insert(c_indices.end(), l2.linkEids.begin(), l2.linkEids.end());
    //     colors.insert(colors.end(), a.begin(), a.end());
    //     colors.insert(colors.end(), b.begin(), b.end());
    //     colorValue += 1;
    // }
    
    /*std::cout << "Writing output file" << std::endl;
    std::ofstream ofs("singularityLinks.vtk");
    ofs << "# vtk DataFile Version 3.0\n"
        << "singularityLinks" << ".vtk\n"
        << "ASCII\n\n"
        << "DATASET UNSTRUCTURED_GRID\n";
    ofs << "POINTS " << mesh.V.size() << " double\n";
    // for (auto& v: source.V) {
    //     if (v.type == FEATURE) c_indices.push_back(v.id);
    // }
    // std::vector<size_t> c_indices = {12, 296};
    // std::cout << c_indices.size() << std::endl;
    for (size_t i = 0; i < mesh.V.size(); i++) {
        ofs << std::fixed << std::setprecision(7) <<  mesh.V.at(i).x << " " <<  mesh.V.at(i).y << " " <<  mesh.V.at(i).z << "\n";
    }
    ofs << "CELLS " << c_indices.size() << " " << 3 * c_indices.size() << std::endl;
    for (size_t i = 0; i < c_indices.size(); i++) {
        auto& e = mesh.E.at(c_indices.at(i));
        ofs << "2 " << e.Vids.at(0) << " " << e.Vids.at(1) << std::endl;
    }
    ofs << "CELL_TYPES " << c_indices.size() << "\n";
    for (size_t i = 0; i < c_indices.size(); i++) {
        ofs << "3" << std::endl;
    }

    ofs << "CELL_DATA " << c_indices.size() << "\n";
    ofs << "SCALARS fixed int\n";
    ofs << "LOOKUP_TABLE default\n";
    for (auto c: colors) {
        ofs << c << "\n";
    }*/

}

std::vector<SingularityLink> SemiGlobalSimplifier::TraceSingularityLinks(Vertex& v, BaseComplexQuad& bc) {
    std::vector<bool> is_mesh_edge_visited(mesh.E.size(), false);
    std::vector<SingularityLink> links;
    for (auto edgeid : v.N_Eids) {
        const Edge& edge = mesh.E.at(edgeid);
        if (!is_mesh_edge_visited[edgeid]) {
            SingularityLink l;
            bc.TraceAlongEdge(v, edge, is_mesh_edge_visited, l.linkVids, l.linkEids);
            l.frontId = l.linkVids.front();
            l.backId = l.linkVids.back();
            links.push_back(l);
        }
    }
    return links;
}

void SemiGlobalSimplifier::PerformGlobalOperations() {
    CheckValidity();
}

size_t SemiGlobalSimplifier::GetFaceId(size_t vid, size_t exclude_vid) {
	auto& v = mesh.V.at(vid);
	size_t res;
	for (auto fid : v.N_Fids) {
		auto& f = mesh.F.at(fid);
		bool found_exclude_vid = false;
		for (auto fvid : f.Vids)
			if (fvid == exclude_vid) {
				found_exclude_vid = true;
				break;
			}
		if (!found_exclude_vid) {
			res = fid;
			break;
		}
	}
	return res;
}

bool SemiGlobalSimplifier::doesCrossBoundary(std::vector<size_t> in, bool isVertex) {
    for (auto id: in) {
        if (isVertex) {
            auto& v = mesh.V.at(id);
            if (v.type == FEATURE || v.isBoundary) {
                return true;
            }
        } else {
            auto& e = mesh.E.at(id);
            if (mesh.V.at(e.Vids.at(0)).type == FEATURE || mesh.V.at(e.Vids.at(0)).isBoundary || mesh.V.at(e.Vids.at(1)).type == FEATURE || mesh.V.at(e.Vids.at(1)).isBoundary) {
                return true;
            }
        }
    }
    return false;
}


size_t SemiGlobalSimplifier::GetDiagonalV(size_t vid, size_t fid) {
	auto& f = mesh.F.at(fid);
    int index = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), vid));
	return f.Vids.at((index+2)%f.Vids.size());
}

void SemiGlobalSimplifier::Smooth() {
    std::cout << "Smoothing mesh" << std::endl;
    // std::vector<size_t> smoothv;
    // for (auto& v: mesh.V) {
    //     if (v.N_Vids.empty() || v.N_Eids.empty() || v.N_Fids.empty()) continue;
    //     smoothv.push_back(v.id);
    // }
    // smoother.Smooth(smoothv);
    smoother.Smooth(std::vector<size_t>{});
}

// std::cout << "Writing output file" << std::endl;
// std::string outputf = "links.vtk";
// std::ofstream ofs(outputf.c_str());
// ofs << "# vtk DataFile Version 3.0\n"
//     << outputf.c_str() << "\n"
//     << "ASCII\n\n"
//     << "DATASET UNSTRUCTURED_GRID\n";
// ofs << "POINTS " << mesh.V.size() << " double\n";
// std::vector<size_t> c_indices;
// for (auto vid: target) {
//     c_indices.push_back(vid);
// }
// for (auto c: collapse) {
//     c_indices.push_back(c.at(0));
//     c_indices.push_back(c.at(1));
// }
// // std::vector<size_t> c_indices = {12, 296};
// // std::cout << c_indices.size() << std::endl;
// for (size_t i = 0; i < mesh.V.size(); i++) {
//     ofs << std::fixed << std::setprecision(7) <<  mesh.V.at(i).x << " " <<  mesh.V.at(i).y << " " <<  mesh.V.at(i).z << "\n";
// }
// ofs << "CELLS " << c_indices.size() << " " << 2 * c_indices.size() << std::endl;
// for (size_t i = 0; i < c_indices.size(); i++) {
//     ofs << "1 " << c_indices.at(i) << std::endl;
// }
// ofs << "CELL_TYPES " << c_indices.size() << "\n";
// for (size_t i = 0; i < c_indices.size(); i++) {
//     ofs << "1" << std::endl;
// }