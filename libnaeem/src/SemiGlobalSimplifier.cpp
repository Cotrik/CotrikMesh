#include <algorithm>
#include <map>
#include <time.h>
#include <queue>
#include <deque>
#include <utility>
#include "ParallelFor.h"
#include "SemiGlobalSimplifier.h"

SemiGlobalSimplifier::SemiGlobalSimplifier() {}

SemiGlobalSimplifier::SemiGlobalSimplifier(Mesh& mesh_, MeshUtil& mu_, Smoother& smoother_) : mesh(&mesh_), mu(&mu_), smoother(&smoother_) {}

SemiGlobalSimplifier::~SemiGlobalSimplifier() {}

void SemiGlobalSimplifier::CheckValidity() {
    if (mesh == NULL) {
        std::cout << "No mesh to use for Semi Global Simplifier." << std::endl;
        exit(0);
    }
    if (mu == NULL) {
        std::cout << "MeshUtil is not initialized for Semi Global Simplifier." << std::endl;
        exit(0);
    }
    if (smoother == NULL) {
        std::cout << "Smoother is not intitalized for Semi Global Simplifier." << std::endl;
        exit(0);
    }
    if (mesh->V.size() == 0 || mesh->F.size() == 0 || mesh->C.size() == 0) {
        std::cout << "No mesh to use for Semi Global Simplifier." << std::endl;
        exit(0);
    }
}

void SemiGlobalSimplifier::SetMembers(Mesh& mesh_, MeshUtil& mu_, Smoother& smoother_) {
    mesh = &mesh_;
    mu = &mu_;
    smoother = &smoother_;
}

void SemiGlobalSimplifier::SetIters(int iters_) {
    iters = iters_;
    smoother->SetIters(iters_);
}

void SemiGlobalSimplifier::SetSimplificationOperations() {
    CheckValidity();
    
    // SetDiagonalCollapseOperations();
    SetDirectSeparatrixOperations(false);
}

bool SemiGlobalSimplifier::FixBoundary() {
    CheckValidity();
    std::cout << "Fixing Boundary ..." << std::endl;
    bool res = false;
    // int i = 0;
    // for (auto& f: mesh->F) {
    //     if (f.N_Fids.size() == 0) continue;
    //     bool skip = false;
    //     for (auto fvid: f.Vids) {
    //         auto& v = mesh->V.at(fvid);
    //         int featureCount = 0;
    //         for (auto vid: v.N_Vids) {
    //             if (mesh->V.at(vid).type == FEATURE) featureCount += 1;
    //         }
    //         if (featureCount > 2) {
    //             skip = true;
    //             break;
    //         }
    //     }
    //     if (skip) continue;
    //     for (int i = 0; i < f.Vids.size(); i++) {
    //         if ((mesh->V.at(f.Vids.at(i)).type == FEATURE || mesh->V.at(f.Vids.at(i)).isBoundary)
    //         && (mesh->V.at(f.Vids.at((i+1)%f.Vids.size())).type == FEATURE || mesh->V.at(f.Vids.at((i+1)%f.Vids.size())).isBoundary)
    //         && (mesh->V.at(f.Vids.at((i+2)%f.Vids.size())).type == FEATURE || mesh->V.at(f.Vids.at((i+2)%f.Vids.size())).isBoundary)) {
    //             std::cout << "fixing degenerate boundary" << std::endl;
    //             std::shared_ptr<SimplificationOperation> dc1 = std::make_shared<DiagonalCollapse>(*mesh, *mu, *smoother, f.id, i, (i+2)%f.Vids.size());
    //             // dc1->SetRanking();
    //             dc1->PerformOperation();
    //             mesh->V.at(f.Vids.at((i+1)%f.Vids.size())).type = REGULAR;
    //             mesh->V.at(f.Vids.at((i+1)%f.Vids.size())).isBoundary = false;
    //             break;
    //         }
    //     }
    // }
    for (int i = 0; i < mesh->V.size(); i++) {
        // if (i >= iters) break;
        auto& v = mesh->V.at(i);
        int valence = v.N_Fids.size();
        bool failed = false;

        if (v.type == FEATURE && valence > v.idealValence) {
            int featureCount = 0;
            for (auto vid: v.N_Vids) {
                if (mesh->V.at(vid).type == FEATURE) featureCount += 1;
            }
            bool performVertexSplit = true;
            if (featureCount != 2) performVertexSplit = false;
            for (int j = 0; j < v.N_Eids.size(); j++) {
                auto& e = mesh->E.at(v.N_Eids.at(j));
                if (mesh->V.at(e.Vids.at(0)).type == FEATURE && mesh->V.at(e.Vids.at(1)).type == FEATURE) continue;
                int count = 0;
                for (auto fid: e.N_Fids) {
                    auto& f = mesh->F.at(fid);
                    for (auto vid: f.Vids) {
                        if (mesh->V.at(vid).type == FEATURE) count += 1;
                    }
                }
                if (count == 4) performVertexSplit = false;
            }
            for (int j = 0; j < v.N_Fids.size(); j++) {
                int count = 0;
                auto& f = mesh->F.at(v.N_Fids.at(j));
                for (auto vid: f.Vids) {
                    if (mesh->V.at(vid).type == FEATURE) count += 1;
                }
                if (count > 2) performVertexSplit = false;
            }
            if (performVertexSplit && v.N_Vids.size() > 5) {
                std::shared_ptr<SimplificationOperation> s = std::make_shared<VertexSplit>(*mesh, *mu, *smoother, v.id);
                s->PerformOperation();
                res = true;
                continue;
            }
            while (valence > v.idealValence) {
                bool breakLoop = false;
                for (int j = 0; j < v.N_Eids.size(); j++) {
                    auto& e = mesh->E.at(v.N_Eids.at(j));
                    if (mesh->V.at(e.Vids.at(0)).type == FEATURE && mesh->V.at(e.Vids.at(1)).type == FEATURE) continue;
                    int count = 0;
                    for (auto fid: e.N_Fids) {
                        auto& f = mesh->F.at(fid);
                        for (auto vid: f.Vids) {
                            if (mesh->V.at(vid).type == FEATURE) count += 1;
                        }
                    }
                    if (count >= 4 && v.idealValence >= 4) continue;
                    bool clockwise = false;
                    for (auto fid: e.N_Fids) {
                        auto& f = mesh->F.at(fid);
                        int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
                        if (mesh->V.at(f.Vids.at((idx+1)%f.Vids.size())).type == FEATURE) {
                            clockwise = true;
                        }
                    }
                    std::shared_ptr<SimplificationOperation> s = std::make_shared<EdgeRotation>(*mesh, *mu, *smoother, e.id, clockwise);
                    s->PerformOperation();
                    res = true;
                }
                if (v.N_Fids.size() == valence) {
                    failed = true;
                    break;
                }
                valence = v.N_Fids.size();
            }
            // i += 1;
            // if (failed) break;
        }
        if (v.isBoundary && valence > v.idealValence) {
            // while (valence > 2) {
                if (v.idealValence >= 2) {
                    int featureCount = 0;
                    for (auto vid: v.N_Vids) {
                        if (mesh->V.at(vid).isBoundary) featureCount += 1;
                    }
                    if (featureCount != 2) continue;
                    // std::cout << "Fixing boundary vertex" << std::endl;
                    std::shared_ptr<SimplificationOperation> s = std::make_shared<VertexSplit>(*mesh, *mu, *smoother, v.id);
                    s->PerformOperation();
                    if (mesh->V.at(i).N_Fids.size() != valence) res = true;
                    // res = true;
                    // std::cout << "After Fixing boundary vertex" << std::endl;
                    // valence = mesh->V.at(i).N_Fids.size();
                } else if (v.idealValence == 1) {
                    bool breakLoop = false;
                    for (int j = 0; j < v.N_Eids.size(); j++) {
                        auto& e = mesh->E.at(v.N_Eids.at(j));
                        if (mesh->V.at(e.Vids.at(0)).isBoundary && mesh->V.at(e.Vids.at(1)).isBoundary) continue;
                        bool clockwise = false;
                        for (auto fid: e.N_Fids) {
                            auto& f = mesh->F.at(fid);
                            int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
                            if (mesh->V.at(f.Vids.at((idx+1)%f.Vids.size())).isBoundary && mesh->V.at(f.Vids.at((idx+2)%f.Vids.size())).N_Fids.size() == 3) {
                                clockwise = true;
                            } else if (mesh->V.at(f.Vids.at((idx+3)%f.Vids.size())).isBoundary && mesh->V.at(f.Vids.at((idx+2)%f.Vids.size())).N_Fids.size() == 3) {
                                clockwise = false;
                            }
                        }
                        std::shared_ptr<SimplificationOperation> s = std::make_shared<EdgeRotation>(*mesh, *mu, *smoother, e.id, clockwise);
                        s->PerformOperation();
                        if (mesh->V.at(i).N_Fids.size() != valence) res = true;
                        // res = true;
                    }
                    // valence = mesh->V.at(i).N_Fids.size();
                }  
            // }
        } else if (v.isBoundary && valence == 1 && valence < v.idealValence) {
            double angle = mesh->GetAngle(v.id, v.N_Vids.at(0), v.N_Vids.at(1));
            if (fabs(180.0-angle) > 30.0) {
                continue;
            }
            auto& f = mesh->F.at(v.N_Fids.at(0));
            int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
            std::shared_ptr<SimplificationOperation> dc = std::make_shared<DiagonalCollapse>(*mesh, *mu, *smoother, f.id, (idx+1)%f.Vids.size(), (idx+3)%f.Vids.size());
            dc->PerformOperation();
        }
    }
    return res;
}

void SemiGlobalSimplifier::SetDiagonalCollapseOperations() {
    CheckValidity();

    for (auto& f: mesh->F) {
        if (f.N_Fids.size() == 0 || f.Vids.empty()) continue;
        for (int i = 0; i < f.Vids.size(); i++) {
            if (mesh->V.at(f.Vids.at(i)).N_Fids.size() == 3 && mesh->V.at(f.Vids.at((i+2)%f.Vids.size())).N_Fids.size() == 3) {
                    std::shared_ptr<SimplificationOperation> dc1 = std::make_shared<DiagonalCollapse>(*mesh, *mu, *smoother, f.id, i, (i+2)%f.Vids.size());
                    dc1->SetRanking();
                    Ops.push_back(dc1);
                break;
            }
        }
        
        // std::unique_ptr<SimplificationOperation> dc2 = std::make_unique<DiagonalCollapse>(*mesh, *mu, *smoother, f.id, 1, 3);
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
    for (auto& v: mesh->V) {
        if (v.N_Vids.size() != 4 || v.type == FEATURE || v.isBoundary) continue;
        std::vector<size_t> c1, c2;
        for (auto vid: v.N_Vids) mesh->V.at(vid).type != FEATURE && !mesh->V.at(vid).isBoundary && mesh->V.at(vid).N_Vids.size() == 3 ? c1.push_back(vid) : c2.push_back(vid);
        if (c1.size() != 2 && c2.size() != 2) continue;
        auto& s3_v1 = mesh->V.at(c1.at(0));
        auto& s3_v2 = mesh->V.at(c1.at(1));
        auto& sn_v1= mesh->V.at(c2.at(0));
        auto& sn_v2= mesh->V.at(c2.at(1));
        if (mu->GetDifference(s3_v1.N_Fids, s3_v2.N_Fids).size() != s3_v1.N_Fids.size()) continue;
        if ((sn_v1.type == FEATURE || sn_v1.isBoundary) && sn_v1.N_Fids.size() < sn_v1.idealValence) continue; 
        if ((sn_v2.type == FEATURE || sn_v2.isBoundary) && sn_v2.N_Fids.size() < sn_v2.idealValence) continue; 
        if (sn_v1.type == FEATURE || sn_v1.isBoundary || sn_v2.type == FEATURE || sn_v2.isBoundary) {
            bool skip = false;
            int featureCount = 0;
            for (auto nvid: v.N_Vids) {
                if (mesh->V.at(nvid).type == FEATURE || mesh->V.at(nvid).isBoundary) featureCount += 1;
            }
            for (auto cvid: c2) {
                for (auto eid: v.N_Eids) {
                    int count = 0;
                    auto& e = mesh->E.at(eid);
                    if ((e.Vids.at(0) == v.id && e.Vids.at(1) == cvid) || (e.Vids.at(1) == v.id && e.Vids.at(0) == cvid)) {
                        for (auto efid: e.N_Fids) {
                            auto& f = mesh->F.at(efid);
                            for (auto fvid: f.Vids) {
                                if (mesh->V.at(fvid).type == FEATURE || mesh->V.at(fvid).isBoundary) count += 1;
                            }
                        }
                        // if (count == 4 && featureCount == 2) skip = true;
                        if (count == 4 && mesh->V.at(cvid).idealValence == 4) skip = true;
                        break;
                    }
                }
            }
            if (skip) continue;
            std::shared_ptr<SimplificationOperation> ds = std::make_shared<DirectSeparatrixCollapse>(*mesh, *mu, *smoother, v.id, c1, c2, looseCollapse);
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
        // smoother->Smooth(op->smoothV);
        for (auto key: op->toUpdate) {
            auto nop = Op_Q.getByKey(key);
            if (!nop) continue;
            nop->SetRanking(op->GetLocation());
            Op_Q.update(nop->ranking, key);
        }
        // i += 1;
        // if (i >= iters) break;
    }
    // std::cout << "Smoothing mesh" << std::endl;
    // std::vector<size_t> smoothv;
    // for (auto& v: mesh->V) smoothv.push_back(v.id);
    // smoother->Smooth(smoothv);
}


void SemiGlobalSimplifier::SetDirectSeparatrixOperations(bool looseCollapse) {
    CheckValidity();

    Op_Q.setMaxQueueOn();
    Op_Q.setSpecialComparisonOn();
    for (auto& v: mesh->V) {
        if (v.N_Vids.size() != 4 || v.type == FEATURE || v.isBoundary) continue;
        std::vector<size_t> c1, c2;
        for (auto vid: v.N_Vids) mesh->V.at(vid).N_Vids.size() == 3 ? c1.push_back(vid) : c2.push_back(vid);
        if (c1.size() != 2 && c2.size() != 2) continue;
        auto& s3_v1 = mesh->V.at(c1.at(0));
        auto& s3_v2 = mesh->V.at(c1.at(1));
        auto& sn_v1= mesh->V.at(c2.at(0));
        auto& sn_v2= mesh->V.at(c2.at(1));
        if (mu->GetDifference(s3_v1.N_Fids, s3_v2.N_Fids).size() != s3_v1.N_Fids.size()) continue;
        if (sn_v1.type == FEATURE || sn_v1.isBoundary || sn_v2.type == FEATURE || sn_v2.isBoundary) continue;

        std::shared_ptr<SimplificationOperation> ds = std::make_shared<DirectSeparatrixCollapse>(*mesh, *mu, *smoother, v.id, c1, c2, looseCollapse);
        ds->SetRanking();
        if (ds->ranking < 0) continue;
        Op_Q.insert(ds->ranking, v.id, ds);
    }
    int i = 0;
    std::cout << Op_Q.size() << " direct separatrix operations" << std::endl;
    while (!Op_Q.empty()) {
        auto& op = Op_Q.pop();
        op->PerformOperation();
        // smoother->Smooth(op->smoothV);
        for (auto key: op->toUpdate) {
            auto nop = Op_Q.getByKey(key);
            if (!nop) continue;
            nop->SetRanking(op->GetLocation());
            Op_Q.update(nop->ranking, key);
        }
    }
    // std::cout << "Smoothing mesh" << std::endl;
    // std::vector<size_t> smoothv;
    // for (auto& v: mesh->V) smoothv.push_back(v.id);
    // smoother->Smooth(smoothv);

}

void SemiGlobalSimplifier::SetSeparatrixOperations() {
    CheckValidity();

    BaseComplexQuad bc(*mesh);
    Op_Q.clear();
    Op_Q.setMinQueueOn();
    Op_Q.setSpecialComparisonOff();
    // Op_Q.setMaxQueueOn();
    for (auto& v: mesh->V) {
        if (!v.isSingularity || v.N_Fids.empty()) continue;
        std::vector<SingularityLink> links = TraceSingularityLinks(v, bc);
        for (int i = 0; i < links.size(); i++) {
            auto& linkV = links.at(i).linkVids;
            auto& linkE = links.at(i).linkEids;
            auto& v_front = mesh->V.at(linkV.front());
            auto& v_back = mesh->V.at(linkV.back());
            if (v_front.isBoundary || v_back.isBoundary) continue;
            if (v_front.N_Fids.size() != 3 || v_back.N_Fids.size() != 3) continue;
            if (mesh->V.at(GetDiagonalV(v_front.id, GetFaceId(v_front.id, linkV.at(1)))).N_Fids.size() <= 4 || mesh->V.at(GetDiagonalV(v_back.id, GetFaceId(v_back.id, linkV.at(linkV.size()-2)))).N_Fids.size() <= 4) continue;

            std::shared_ptr<SimplificationOperation> s = std::make_shared<SeparatrixCollapse>(*mesh, *mu, *smoother, linkV, linkE);
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
            if (!mesh->V.at(vid).isSingularity || mesh->V.at(vid).N_Fids.empty()) continue;
            // std::cout << mesh->V.at(vid).N_Vids.size() << " " << mesh->V.at(vid).N_Eids.size() << " " << mesh->V.at(vid).N_Fids.size() << std::endl;
            std::vector<SingularityLink> links = TraceSingularityLinks(mesh->V.at(vid), bc);
            // std::cout << "links " << links.size() << std::endl;
            for (int i = 0; i < links.size(); i++) {
                auto& linkV = links.at(i).linkVids;
                auto& linkE = links.at(i).linkEids;
                auto& v_front = mesh->V.at(linkV.front());
                auto& v_back = mesh->V.at(linkV.back());
                if (v_front.isBoundary || v_back.isBoundary) continue;
                if (v_front.N_Fids.size() != 3 || v_back.N_Fids.size() != 3) continue;
                if (mesh->V.at(GetDiagonalV(v_front.id, GetFaceId(v_front.id, linkV.at(1)))).N_Fids.size() <= 4 || mesh->V.at(GetDiagonalV(v_back.id, GetFaceId(v_back.id, linkV.at(linkV.size()-2)))).N_Fids.size() <= 4) continue;

                // std::cout << "making a separatrix operation" << std::endl;
                std::shared_ptr<SimplificationOperation> s = std::make_shared<SeparatrixCollapse>(*mesh, *mu, *smoother, linkV, linkE);
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
    // for (auto& v: mesh->V) {
    //     if (v.N_Vids.empty() || v.N_Eids.empty() || v.N_Fids.empty()) continue;
    //     smoothv.push_back(v.id);
    // }
    // smoother->Smooth(smoothv);
}

void SemiGlobalSimplifier::SetBoundarySeparatrixOperations() {
    CheckValidity();

    BaseComplexQuad bc(*mesh);
    Op_Q.clear();
    Op_Q.setMinQueueOn();
    Op_Q.setSpecialComparisonOff();
    // Op_Q.setMaxQueueOn();
    for (auto& v: mesh->V) {
        if (!v.isSingularity || v.N_Fids.empty()) continue;
        std::vector<SingularityLink> links = TraceSingularityLinks(v, bc);
        for (int i = 0; i < links.size(); i++) {
            auto& linkV = links.at(i).linkVids;
            auto& linkE = links.at(i).linkEids;
            if (linkV.size() != 3) continue;
            std::vector<size_t> verticesToCheck(linkV.begin(), linkV.end());
            for (auto linkVid: linkV) {
                mu->AddContents(verticesToCheck, mesh->V.at(linkVid).N_Vids);
                
            }
            bool skip = false;
            for (auto vid: verticesToCheck) {
                int featureCount = 0;
                for (auto nvid: mesh->V.at(vid).N_Vids) {
                    if (mesh->V.at(nvid).type == FEATURE || mesh->V.at(nvid).isBoundary) featureCount += 1;
                }
                if (featureCount > 2) {
                    skip = true;
                    break;
                }
            }
            if (skip) continue;
            auto& v_front = mesh->V.at(linkV.front());
            auto& v_back = mesh->V.at(linkV.back());
            if (v_front.isBoundary || v_front.type == FEATURE || v_back.isBoundary || v_back.type == FEATURE) continue;
            if (v_front.N_Fids.size() != 3 || v_back.N_Fids.size() != 3) continue;
            // if (mesh->V.at(GetDiagonalV(v_front.id, GetFaceId(v_front.id, linkV.at(1)))).N_Fids.size() <= 4 || mesh->V.at(GetDiagonalV(v_back.id, GetFaceId(v_back.id, linkV.at(linkV.size()-2)))).N_Fids.size() <= 4) continue;
            auto& midV = mesh->V.at(linkV.at(1));
            if (midV.type == FEATURE || midV.isBoundary) {
                std::shared_ptr<SimplificationOperation> s = std::make_shared<SeparatrixCollapse>(*mesh, *mu, *smoother, linkV, linkE);
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
    // for (auto& v: mesh->V) {
    //     if (v.N_Vids.empty() || v.N_Eids.empty() || v.N_Fids.empty()) continue;
    //     smoothv.push_back(v.id);
    // }
    // smoother->Smooth(smoothv);
}

void SemiGlobalSimplifier::SetHalfSeparatrixOperations() {
    CheckValidity();

    BaseComplexQuad bc(*mesh);
    Op_Q.clear();
    Op_Q.setMaxQueueOn();
    Op_Q.setSpecialComparisonOff();
    for (auto& v: mesh->V) {
        if (!v.isSingularity) continue;
        std::vector<SingularityLink> links = TraceSingularityLinks(v, bc);
        for (int i = 0; i < links.size(); i++) {
            auto& linkV = links.at(i).linkVids;
            auto& linkE = links.at(i).linkEids;
            auto& v_front = mesh->V.at(linkV.front());
            auto& v_back = mesh->V.at(linkV.back());
                
            if (!(v_front.isBoundary ^ v_back.isBoundary)) continue;
            if (v_front.isBoundary && !v_back.isBoundary) {
                std::reverse(linkV.begin(), linkV.end());
                std::reverse(linkE.begin(), linkE.end());
                auto& v_front = mesh->V.at(linkV.front());
                auto& v_back = mesh->V.at(linkV.back());
            }
            if (!(v_front.N_Fids.size() == 3 && v_back.N_Fids.size() == 2)) continue;
            size_t fid = mu->GetDifference(v_front.N_Fids, std::vector<size_t>{mesh->E.at(linkE.at(0)).N_Fids}).at(0);
            auto& f = mesh->F.at(fid);
            int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v_front.id));
            int diagVid = f.Vids.at((idx+2)%f.Vids.size());
            auto& diagV = mesh->V.at(diagVid);
            if (diagV.N_Fids.size() <= 4) continue;
            std::shared_ptr<SimplificationOperation> s = std::make_shared<SeparatrixCollapse>(*mesh, *mu, *smoother, linkV, linkE, true);
            s->SetRanking();
            Op_Q.insert(s->ranking, s->GetCenterId(), s);
        }
    }
    std::cout << Op_Q.size() << " half separatrix operations" << std::endl;
    int i = 0;
    while (!Op_Q.empty()) {
        auto& op = Op_Q.pop();
        op->PerformOperation();
        // for (auto vid: op->toUpdate) {
        //     std::vector<SingularityLink> links = TraceSingularityLinks(mesh->V.at(vid), bc);
        //     for (int i = 0; i < links.size(); i++) {
        //         auto& linkV = links.at(i).linkVids;
        //         auto& linkE = links.at(i).linkEids;
        //         auto& v_front = mesh->V.at(linkV.front());
        //         auto& v_back = mesh->V.at(linkV.back());
                
        //         if (!(v_front.isBoundary ^ v_back.isBoundary)) continue;
        //         if (v_front.isBoundary && !v_back.isBoundary) {
        //             std::reverse(linkV.begin(), linkV.end());
        //             std::reverse(linkE.begin(), linkE.end());
        //             auto& v_front = mesh->V.at(linkV.front());
        //             auto& v_back = mesh->V.at(linkV.back());
        //         }
        //         if (!(v_front.N_Fids.size() == 3 && v_back.N_Fids.size() == 2)) continue;
        //         std::shared_ptr<SimplificationOperation> s = std::make_shared<SeparatrixCollapse>(*mesh, *mu, *smoother, linkV, linkE, true);
        //         s->SetRanking();
        //         Op_Q.insert(s->ranking, s->GetCenterId(), s);
        //     }
        // }
        i += 1;
        // if (i >= iters) break;
    }
    // std::cout << "Smoothing mesh" << std::endl;
    // std::vector<size_t> smoothv;
    // for (auto& v: mesh->V) {
    //     if (v.N_Vids.empty() || v.N_Eids.empty() || v.N_Fids.empty()) continue;
    //     smoothv.push_back(v.id);
    // }
    // smoother->Smooth(smoothv);
}

void SemiGlobalSimplifier::SetChordCollapseOperations() {
    CheckValidity();

    ChordExtractor ce(*mesh);
    ce.Extract();
    std::vector<size_t> chords = ce.SelectChords();
    // ce.Write(std::vector<size_t>{chords.at(0)});
    // ce.Write(chords);
    // return;

    // int maxId = -1;
    // int rank = 0;
    for (auto chordId: chords) {
        // std::shared_ptr<SimplificationOperation> s = std::make_shared<ChordCollapse>(*mesh, *mu, *smoother, ce, chord.id);
        std::shared_ptr<SimplificationOperation> s = std::make_shared<ChordCollapse>(*mesh, *mu, *smoother, ce, chordId);
        s->CalculateRanking();
        // s->PerformOperation();
        // break;
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
    // ofs << "POINTS " << mesh->V.size() << " double\n";
    // // for (auto& v: source.V) {
    // //     if (v.type == FEATURE) c_indices.push_back(v.id);
    // // }
    // // std::vector<size_t> c_indices = {12, 296};
    // // std::cout << c_indices.size() << std::endl;
    // for (size_t i = 0; i < mesh->V.size(); i++) {
    //     ofs << std::fixed << std::setprecision(7) <<  mesh->V.at(i).x << " " <<  mesh->V.at(i).y << " " <<  mesh->V.at(i).z << "\n";
    // }
    // ofs << "CELLS " << c_indices.size() << " " << 3 * c_indices.size() << std::endl;
    // for (size_t i = 0; i < c_indices.size(); i++) {
    //     auto& e = mesh->E.at(c_indices.at(i));
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
    // for (auto& v: mesh->V) {
    //     if (v.N_Vids.empty() || v.N_Eids.empty() || v.N_Fids.empty()) continue;
    //     smoothv.push_back(v.id);
    // }
    // smoother->Smooth(smoothv);
}

void SemiGlobalSimplifier::SetEdgeRotationOperations() {
    CheckValidity();

    int i = 0;
    for (auto& e: mesh->E) {
        std::shared_ptr<SimplificationOperation> s = std::make_shared<EdgeRotation>(*mesh, *mu, *smoother, e.id, true);
        s->PerformOperation();
        i += 1;
        // if (i > iters) break;
    }
    // std::cout << "Smoothing mesh" << std::endl;
    // std::vector<size_t> smoothv;
    // for (auto& v: mesh->V) {
    //     if (v.N_Vids.empty() || v.N_Eids.empty() || v.N_Fids.empty()) continue;
    //     smoothv.push_back(v.id);
    // }
    // smoother->Smooth(smoothv);
}

void SemiGlobalSimplifier::SetVertexRotationOperations() {
    CheckValidity();

    int i = 0;
    for (auto& v: mesh->V) {
        std::shared_ptr<SimplificationOperation> s = std::make_shared<VertexRotation>(*mesh, *mu, *smoother, v.id);
        s->PerformOperation();
        // i += 1;
        // if (i > iters) break;
    }
    // std::cout << "Smoothing mesh" << std::endl;
    // std::vector<size_t> smoothv;
    // for (auto& v: mesh->V) {
    //     if (v.N_Vids.empty() || v.N_Eids.empty() || v.N_Fids.empty()) continue;
    //     smoothv.push_back(v.id);
    // }
    // smoother->Smooth(smoothv);
}

void SemiGlobalSimplifier::SetVertexSplitOperations() {
    CheckValidity();

    int i = 0;
    for (int j = 0; j < mesh->V.size(); j++) {
        auto& v = mesh->V.at(j);
        if (v.type != FEATURE) continue;
        if (v.N_Vids.size() < 6) continue;
        std::shared_ptr<SimplificationOperation> s = std::make_shared<VertexSplit>(*mesh, *mu, *smoother, v.id);
        s->PerformOperation();
        // i += 1;
        // if (i > iters) break;
        // break;
    }
}

void SemiGlobalSimplifier::SetEdgeCollapseOperations() {
    CheckValidity();

    Edge& e = mesh->E.at(0);
    size_t source_id = e.Vids.at(0);
    size_t targte_id = e.Vids.at(1);
    std::shared_ptr<SimplificationOperation> vr = std::make_shared<VertexRotation>(*mesh, *mu, *smoother, targte_id);
    std::shared_ptr<SimplificationOperation> dc = std::make_shared<DiagonalCollapse>(*mesh, *mu, *smoother, -1, targte_id, source_id);
    std::shared_ptr<SimplificationOperation> s = std::make_shared<EdgeCollapse>(*mesh, *mu, *smoother, vr, dc);
    s->PerformOperation();
}

void SemiGlobalSimplifier::SetQuadSplitOperations() {
    CheckValidity();

    size_t startId = -1;
    for (auto& v: mesh->V) {
        bool skip = false;
        if (v.type == FEATURE || v.isBoundary || v.N_Vids.size() != 5) skip = true;
        for (auto vid: v.N_Vids) if (mesh->V.at(vid).type == FEATURE || mesh->V.at(vid).isBoundary || mesh->V.at(vid).N_Vids.size() != 4) skip = true;
        if (skip) continue;
        int d = (v.N_Vids.size() / 2) + 1;
        bool clockwise = false;
        size_t mainV = v.N_Vids.at(0);
        std::vector<size_t> verticesToSplit;
        std::vector<size_t> verticesToChange;
        size_t startE;
        for (auto eid: v.N_Eids) {
            auto& e = mesh->E.at(eid);
            if (e.Vids.at(0) == mainV || e.Vids.at(1) == mainV) {
                startE = e.id;
            }
        }
        for (int j = 0; j < d; j++) {
            auto& e = mesh->E.at(startE);
            size_t ev = e.Vids.at(0) == v.id ? e.Vids.at(1) : e.Vids.at(0);
            for (auto fid: e.N_Fids) {
                auto& f = mesh->F.at(fid);
                int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
                if (!clockwise && f.Vids.at((idx+1)%f.Vids.size()) == ev) {
                    if (j == 0 || j == d-1) {
                        verticesToSplit.push_back(ev);
                    } else {
                        verticesToChange.push_back(ev);
                    }
                    startE = mu->GetDifference(mu->GetIntersection(f.Eids, v.N_Eids), std::vector<size_t>{e.id}).at(0);
                } else if (clockwise && f.Vids.at((idx+3)%f.Vids.size()) == ev) {
                    if (j == 0 || j == d-1) {
                        verticesToSplit.push_back(ev);
                    } else {
                        verticesToChange.push_back(ev);
                    }
                    startE = mu->GetDifference(mu->GetIntersection(f.Eids, v.N_Eids), std::vector<size_t>{e.id}).at(0);
                }
            }
        }
        startId = v.id;
        break;
    }
    // return;
    if (startId != -1) {
        auto& start = mesh->V.at(startId);
        size_t threeId = -1;
        size_t fiveId = -1;
        for (auto vid: start.N_Vids) {
            std::cout << mesh->V.at(vid).N_Vids.size() << std::endl; 
            if (mesh->V.at(vid).N_Vids.size() == 5) {
                fiveId = vid;
                break;
            }
        }
        if (fiveId == -1) return;
        auto& five = mesh->V.at(fiveId);
        for (auto vid: five.N_Vids) {
            std::cout << mesh->V.at(vid).N_Vids.size() << std::endl; 
            if (mesh->V.at(vid).N_Vids.size() == 3) {
                threeId = vid;
                break;
            }
        }
        if (threeId == -1) return;
        auto& three = mesh->V.at(threeId);
        std::cout << fiveId << " " << five.N_Vids.size() << " " << threeId << " " << three.N_Vids.size() << std::endl;
        std::unique_ptr<ThreeFivePair> tfp = std::make_unique<ThreeFivePair>(*mesh, *mu, *smoother, threeId, fiveId);
        tfp->MoveLowerLeft();
        tfp->MoveLowerLeft();
        tfp->MoveLowerLeft();
    }
    return;

    int it = 0;
    for (int i = 0; i < mesh->V.size(); i++) {
        auto& v = mesh->V.at(i);
        if (v.type != FEATURE && v.N_Vids.size() > 4) {
            // std::cout << "Quad Split" << std::endl;
            // std::cout << "i: " << i << " " << mesh->V.size() << std::endl;
            int d = (v.N_Vids.size() / 2) + 1;
            bool clockwise = false;
            size_t mainV = v.N_Vids.at(0);
            std::vector<size_t> verticesToSplit;
            std::vector<size_t> verticesToChange;
            size_t startE;
            for (auto eid: v.N_Eids) {
                auto& e = mesh->E.at(eid);
                if (e.Vids.at(0) == mainV || e.Vids.at(1) == mainV) {
                    startE = e.id;
                }
            }
            for (int j = 0; j < d; j++) {
                auto& e = mesh->E.at(startE);
                size_t ev = e.Vids.at(0) == v.id ? e.Vids.at(1) : e.Vids.at(0);
                for (auto fid: e.N_Fids) {
                    auto& f = mesh->F.at(fid);
                    int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
                    if (!clockwise && f.Vids.at((idx+1)%f.Vids.size()) == ev) {
                        if (j == 0 || j == d-1) {
                            verticesToSplit.push_back(ev);
                        } else {
                            verticesToChange.push_back(ev);
                        }
                        startE = mu->GetDifference(mu->GetIntersection(f.Eids, v.N_Eids), std::vector<size_t>{e.id}).at(0);
                    } else if (clockwise && f.Vids.at((idx+3)%f.Vids.size()) == ev) {
                        if (j == 0 || j == d-1) {
                            verticesToSplit.push_back(ev);
                        } else {
                            verticesToChange.push_back(ev);
                        }
                        startE = mu->GetDifference(mu->GetIntersection(f.Eids, v.N_Eids), std::vector<size_t>{e.id}).at(0);
                    }
                }
            }
            // if (verticesToSplit.size() < 2) continue;
            std::shared_ptr<SimplificationOperation> qs = std::make_shared<QuadSplit>(*mesh, *mu, *smoother, v.id, verticesToSplit, verticesToChange);
            qs->PerformOperation();
            it += 1;
            // if (it >= iters) break;
        }    
    }
}

void SemiGlobalSimplifier::ResolveSingularityPairs() {
    CheckValidity();

    for (int i = 0; i < mesh->V.size(); i++) {
        auto& v = mesh->V.at(i);
        if (v.type != FEATURE && !v.isBoundary && v.N_Vids.size() == 3) {
            CheckAndResolveThreeFivePair(v.id);
        } else if (v.type != FEATURE && !v.isBoundary && v.N_Vids.size() == 5) {
            CheckAndResolveFiveThreePair(v.id);
        }
        // if (i >= iters) break;
    }
}

void SemiGlobalSimplifier::CheckAndResolveThreeFivePair(size_t vid) {
    auto& v = mesh->V.at(vid);
    bool clockwise = false;
    bool isPair = false;
    for (auto fid: v.N_Fids) {
        auto& f = mesh->F.at(fid);
        bool skip = false;
        for (auto fvid: f.Vids) {
            if (mesh->V.at(fvid).type == FEATURE || mesh->V.at(fvid).isBoundary) {
                skip = true;
                break;
            }
        }
        if (skip) continue;
        int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), vid));
        if (mesh->V.at(f.Vids.at((idx+2)%f.Vids.size())).N_Vids.size() == 5 && mesh->V.at(f.Vids.at((idx+3)%f.Vids.size())).N_Vids.size() == 5 && mesh->V.at(f.Vids.at((idx+1)%f.Vids.size())).N_Vids.size() == 4) {
            isPair = true;
        } else if (mesh->V.at(f.Vids.at((idx+1)%f.Vids.size())).N_Vids.size() == 5 && mesh->V.at(f.Vids.at((idx+2)%f.Vids.size())).N_Vids.size() == 5 && mesh->V.at(f.Vids.at((idx+3)%f.Vids.size())).N_Vids.size() == 4) {
            isPair = true;
            clockwise = true;
        } else if (mesh->V.at(f.Vids.at((idx+1)%f.Vids.size())).N_Vids.size() == 5 && (mesh->V.at(f.Vids.at((idx+2)%f.Vids.size())).N_Vids.size() == 4 || mesh->V.at(f.Vids.at((idx+2)%f.Vids.size())).N_Vids.size() == 3) && mesh->V.at(f.Vids.at((idx+3)%f.Vids.size())).N_Vids.size() == 5) {
            std::shared_ptr<SimplificationOperation> dc1 = std::make_shared<DiagonalCollapse>(*mesh, *mu, *smoother, f.id, idx, (idx+2)%f.Vids.size());
            dc1->PerformOperation();
            break;
        }
        if (isPair) {
            for (auto eid: f.Eids) {
                auto& e = mesh->E.at(eid);
                if (mesh->V.at(e.Vids.at(0)).N_Vids.size() == 5 && mesh->V.at(e.Vids.at(1)).N_Vids.size() == 5) {
                    std::shared_ptr<SimplificationOperation> s = std::make_shared<EdgeRotation>(*mesh, *mu, *smoother, e.id, clockwise);
                    s->PerformOperation();
                }
            }
            break;
        }
    }
}


void SemiGlobalSimplifier::CheckAndResolveFiveThreePair(size_t vid) {
    auto& v = mesh->V.at(vid);
    for (auto eid: v.N_Eids) {
        auto& e = mesh->E.at(eid);
        bool isPair = false;
        size_t neid = e.Vids.at(0) == vid ? e.Vids.at(1) : e.Vids.at(0);
        if (mesh->V.at(neid).N_Vids.size() != 4) continue;
        auto& f1 = mesh->F.at(e.N_Fids.at(0));
        auto& f2 = mesh->F.at(e.N_Fids.at(1));
        std::vector<size_t> fvids = f1.Vids;
        fvids.insert(fvids.end(), f2.Vids.begin(), f2.Vids.end());
        bool skip = false;
        for (auto fvid: fvids) {
            if (mesh->V.at(fvid).type == FEATURE || mesh->V.at(fvid).isBoundary) {
                skip = true;
                break;
            }
        }
        if (skip) continue;
        int idx1 = std::distance(f1.Vids.begin(), std::find(f1.Vids.begin(), f1.Vids.end(), vid));
        int idx2 = std::distance(f2.Vids.begin(), std::find(f2.Vids.begin(), f2.Vids.end(), vid));
        if ((mesh->V.at(f1.Vids.at((idx1+1)%f1.Vids.size())).N_Vids.size() == 3 && mesh->V.at(f2.Vids.at((idx2+2)%f2.Vids.size())).N_Vids.size() == 3)
        || (mesh->V.at(f1.Vids.at((idx1+2)%f1.Vids.size())).N_Vids.size() == 3 && mesh->V.at(f2.Vids.at((idx2+1)%f2.Vids.size())).N_Vids.size() == 3)) {
            std::shared_ptr<SimplificationOperation> s = std::make_shared<EdgeRotation>(*mesh, *mu, *smoother, e.id, false);
            s->PerformOperation();
            isPair = true;
        } else if (mesh->V.at(f1.Vids.at((idx1+1)%f1.Vids.size())).N_Vids.size() == 3 && mesh->V.at(f2.Vids.at((idx2+3)%f2.Vids.size())).N_Vids.size() == 3) {
            std::vector<size_t> verticesToSplit = {f1.Vids.at((idx1+1)%f1.Vids.size()), f2.Vids.at((idx2+3)%f2.Vids.size())};
            std::vector<size_t> verticesToChange = {neid};
            std::shared_ptr<SimplificationOperation> qs = std::make_shared<QuadSplit>(*mesh, *mu, *smoother, v.id, verticesToSplit, verticesToChange);
            qs->PerformOperation();
            isPair = true;
        } else if (mesh->V.at(f1.Vids.at((idx1+3)%f1.Vids.size())).N_Vids.size() == 3 && mesh->V.at(f2.Vids.at((idx2+1)%f2.Vids.size())).N_Vids.size() == 3) {
            std::vector<size_t> verticesToSplit = {f1.Vids.at((idx1+3)%f1.Vids.size()), f2.Vids.at((idx2+1)%f2.Vids.size())};
            std::vector<size_t> verticesToChange = {neid};
            std::shared_ptr<SimplificationOperation> qs = std::make_shared<QuadSplit>(*mesh, *mu, *smoother, v.id, verticesToSplit, verticesToChange);
            qs->PerformOperation();
            isPair = true;
        } else if ((mesh->V.at(f1.Vids.at((idx1+2)%f1.Vids.size())).N_Vids.size() == 3 && mesh->V.at(f2.Vids.at((idx2+3)%f2.Vids.size())).N_Vids.size() == 3)
        || (mesh->V.at(f1.Vids.at((idx1+3)%f1.Vids.size())).N_Vids.size() == 3 && mesh->V.at(f2.Vids.at((idx2+2)%f2.Vids.size())).N_Vids.size() == 3)) {
            std::shared_ptr<SimplificationOperation> s = std::make_shared<EdgeRotation>(*mesh, *mu, *smoother, e.id, true);
            s->PerformOperation();
            isPair = true;
        }
        if (isPair) break;
    }
}

bool SemiGlobalSimplifier::ResolveHighValences() {
    CheckValidity();
    std::cout << "Resolving high valences ..." << std::endl;
    bool res = false;
    for (int vid = 0; vid < mesh->V.size(); vid++) {
        auto& v = mesh->V.at(vid);
        if (v.N_Vids.size() == 6 && v.type != FEATURE && !v.isBoundary) {
            size_t startE = v.N_Eids.at(0);
            size_t secE;
            bool foundCandidate = false;
            std::vector<size_t> avoidableE;
            for (int i = 0; i < v.N_Eids.size(); i++) {
                auto& e = mesh->E.at(startE);
                size_t ev = e.Vids.at(0) == v.id ? e.Vids.at(1) : e.Vids.at(0);
                for (auto fid: e.N_Fids) {
                    auto& f = mesh->F.at(fid);
                    int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
                    if (f.Vids.at((idx+1)%f.Vids.size()) == ev) {
                        if (mesh->V.at(f.Vids.at((idx+2)%f.Vids.size())).type != FEATURE && !mesh->V.at(f.Vids.at((idx+2)%f.Vids.size())).isBoundary && mesh->V.at(f.Vids.at((idx+2)%f.Vids.size())).N_Vids.size() == 3) {
                            avoidableE.push_back(startE);
                            secE = startE;
                            foundCandidate = true;
                        }    
                        startE = mu->GetDifference(mu->GetIntersection(f.Eids, v.N_Eids), std::vector<size_t>{e.id}).at(0);
                    }
                }
                if (foundCandidate) {
                    avoidableE.push_back(startE);
                    break;
                }
            }
            bool canSplit = false;
            if (foundCandidate) {
                startE = secE;
                for (int i = 0; i < 4; i++) {
                    auto& e = mesh->E.at(startE);
                    size_t ev = e.Vids.at(0) == v.id ? e.Vids.at(1) : e.Vids.at(0);
                    for (auto fid: e.N_Fids) {
                        auto& f = mesh->F.at(fid);
                        int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
                        if (f.Vids.at((idx+1)%f.Vids.size()) == ev) {
                            if (i == 3 && mesh->V.at(f.Vids.at((idx+2)%f.Vids.size())).type != FEATURE && !mesh->V.at(f.Vids.at((idx+2)%f.Vids.size())).isBoundary && mesh->V.at(f.Vids.at((idx+2)%f.Vids.size())).N_Vids.size() == 3) {
                                canSplit = true;
                                avoidableE.push_back(startE);
                            }    
                            startE = mu->GetDifference(mu->GetIntersection(f.Eids, v.N_Eids), std::vector<size_t>{e.id}).at(0);
                        }
                    }
                    if (canSplit) {
                        avoidableE.push_back(startE);
                    }
                }
            }
            if (canSplit) {
                std::vector<size_t> edgesToSplit = mu->GetDifference(v.N_Eids, avoidableE);
                std::shared_ptr<SimplificationOperation> s = std::make_shared<VertexSplit>(*mesh, *mu, *smoother, v.id, edgesToSplit);
                s->PerformOperation();
                res = true;
                continue;
            }
        }
        if (v.N_Vids.size() > 5 && v.type != FEATURE && !v.isBoundary) {
            int d = (v.N_Vids.size() / 2) + 1;
            size_t mainV;
            bool foundThree = false;
            for (auto nvid: v.N_Vids) {
                if (mesh->V.at(nvid).N_Vids.size() == 3 && mesh->V.at(nvid).type != FEATURE && !mesh->V.at(nvid).isBoundary) {
                    mainV = nvid;
                    foundThree = true;
                    break;
                }
            }
            if (!foundThree) {
                for (auto nvid: v.N_Vids) {
                    if (mesh->V.at(nvid).type != FEATURE && !mesh->V.at(nvid).isBoundary) {
                        mainV = nvid;
                        break;
                    }
                }
            }
            std::vector<size_t> verticesToSplit;
            std::vector<size_t> verticesToChange;
            size_t startE;
            for (auto eid: v.N_Eids) {
                auto& e = mesh->E.at(eid);
                if (e.Vids.at(0) == mainV || e.Vids.at(1) == mainV) {
                    startE = e.id;
                }
            }
            for (int j = 0; j < d; j++) {
                auto& e = mesh->E.at(startE);
                size_t ev = e.Vids.at(0) == v.id ? e.Vids.at(1) : e.Vids.at(0);
                for (auto fid: e.N_Fids) {
                    auto& f = mesh->F.at(fid);
                    int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
                    if (f.Vids.at((idx+1)%f.Vids.size()) == ev) {
                        if (j == 0 || j == d-1) {
                            verticesToSplit.push_back(ev);
                        } else {
                            verticesToChange.push_back(ev);
                        }
                        startE = mu->GetDifference(mu->GetIntersection(f.Eids, v.N_Eids), std::vector<size_t>{e.id}).at(0);
                    }
                }
            }
            std::shared_ptr<SimplificationOperation> qs = std::make_shared<QuadSplit>(*mesh, *mu, *smoother, v.id, verticesToSplit, verticesToChange);
            qs->PerformOperation();
            res = true;
        }
    }
    return res;
}

SingularityLink SemiGlobalSimplifier::PrototypeGetLink(size_t vid, BaseComplexQuad& bc, size_t vertexToSkip, std::vector<size_t> edgesToCheck, bool checkValence, bool boundary) {
    auto& v = mesh->V.at(vid);
    std::priority_queue<SingularityLink, std::vector<SingularityLink>, LinkComparator> q;
    std::vector<SingularityLink> links = GetLinks(v.id, bc, checkValence);
    int valenceToCheck = mesh->V.at(vid).N_Vids.size() == 3 ? 5 : 3;
    for (auto& l: links) {
        // if (PrototypeCheckBoundarySingularity(l.frontId)) continue;
        // if (PrototypeCheckBoundarySingularity(l.frontId) || PrototypeCheckBoundarySingularity(l.backId)) continue;
        if (!edgesToCheck.empty() && (l.backId == vertexToSkip || !IsExclusive(l.frontId, l.linkEids, edgesToCheck))) continue;
        if (!checkValence) {
            if (PrototypeCheckBoundarySingularity(l.backId)) continue;
            if (l.backId == vertexToSkip) continue;
        //     // if (PrototypeCheckBoundarySingularity(l.frontId) || PrototypeCheckBoundarySingularity(l.backId)) l.rank += 3.0;
        //     if (PrototypeCheckBoundarySingularity(l.frontId) || PrototypeCheckBoundarySingularity(l.backId)) continue;
            // if (mesh->V.at(l.backId).N_Vids.size() != valenceToCheck) l.rank += 100;
        // } else if (PrototypeCheckBoundarySingularity(l.frontId) || PrototypeCheckBoundarySingularity(l.backId)) {
        } else if (PrototypeCheckBoundarySingularity(l.backId)) {
            l.rank = 0.0;
        }
        // if (boundary && (PrototypeCheckBoundarySingularity(l.frontId) || PrototypeCheckBoundarySingularity(l.backId))) {
        //     q.push(l);
        //     continue;
        // }
        q.push(l);
    }
    SingularityLink l;
    if (!q.empty()) l = q.top();
    return l;
}

void SemiGlobalSimplifier::PrototypeE() {
    CheckValidity();

    BaseComplexQuad bc(*mesh);
    std::priority_queue<SingularityGroup, std::vector<SingularityGroup>, GroupComparator> q;
    for (auto& v: mesh->V) {
        if (v.type == FEATURE || v.isBoundary || PrototypeCheckBoundarySingularity(v.id) || (v.N_Vids.size() != 3 && v.N_Vids.size() != 5)) continue;
        SingularityLink l1 = PrototypeGetLink(v.id, bc, 0, std::vector<size_t>{}, true, false);
        if (l1.linkVids.empty()) continue;
        SingularityLink l2 = PrototypeGetLink(v.id, bc, l1.backId, l1.linkEids, true, false);
        if (l2.linkVids.empty()) continue;
        SingularityGroup s;
        s.l1 = l1;
        s.l2 = l2;
        s.rank = s.l1.rank + s.l2.rank + (s.l1.a + s.l1.b) + (s.l2.a + s.l2.b);
        q.push(s);
    }
    std::cout << q.size() << std::endl;
    int i = 0;
    
    while (!q.empty()) {
        SingularityGroup s = q.top();
        SingularityLink l1 = s.l1;
        SingularityLink l2 = s.l2;
        if (i >= iters && ((l1.a+l1.b > 2 && l1.b > 0) || (l2.a+l2.b > 2 && l2.b > 0))) {
            std::cout << "l1: " << l1.a+l1.b << " l1 rots: " << l1.rot << std::endl;
            std::cout << "l2: " << l2.a+l2.b << " l2 rots: " << l2.rot << std::endl;
            std::cout << "link rotation from l2 to l1: " << PrototypeGetRotations(l1.frontId, l2.linkEids.front(), l1.linkEids.front()) << std::endl;
            std::cout << "link rotation from l1 to l2: " << PrototypeGetRotations(l1.frontId, l1.linkEids.front(), l2.linkEids.front()) << std::endl;
            PrototypeSaveMesh(l1, l2, "test");
            break;
        }
        // if ((s.l2.a+s.l2.b) < (s.l1.a+s.l1.b)) {
        //     l1 = s.l2;
        //     l2 = s.l1;
        // }
        // PrototypeResolveGroup(l1, l2);
        q.pop();
        i += 1;
        // if (i >= iters) {
        //     break;
        // }
    }
    // while (FixValences());
}

int SemiGlobalSimplifier::PrototypeGetRotations(size_t vid, size_t start, size_t end) {
    auto& v = mesh->V.at(vid);
    int nRots = 0;
    for (int i = 0; i < v.N_Eids.size(); i++) {
        auto& e = mesh->E.at(start);
        size_t evid = e.Vids.at(0) == v.id ? e.Vids.at(1) : e.Vids.at(0);
        for (auto fid: e.N_Fids) {
            auto& f = mesh->F.at(fid);
            int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
            if (f.Vids.at((idx+1)%f.Vids.size()) == evid) {
                start = mu->GetDifference(mu->GetIntersection(v.N_Eids, f.Eids), std::vector<size_t>{e.id}).at(0);
                break;
            }
        }
        nRots += 1;
        if (start == end) break;
    }
    return nRots;
}

void SemiGlobalSimplifier::PrototypeSaveMesh(SingularityLink& l1, SingularityLink& l2, std::string in) {
    int colorValue = 0;
    int ncolors = 15;
    std::vector<size_t> c_indices;
    std::vector<int> colors;
    
    c_indices.insert(c_indices.end(), l1.linkEids.begin(), l1.linkEids.end());
    c_indices.insert(c_indices.end(), l2.linkEids.begin(), l2.linkEids.end());
    std::vector<int> a(l1.linkEids.size(), (colorValue%ncolors));
    colors.insert(colors.end(), a.begin(), a.end());
    colorValue += 1;
    a.clear();
    a.resize(l2.linkEids.size(), (colorValue%ncolors));
    colors.insert(colors.end(), a.begin(), a.end());    
    
    std::ofstream ofs(in+"_singularityLinks.vtk");
    ofs << "# vtk DataFile Version 3.0\n"
        << "singularityLinks" << ".vtk\n"
        << "ASCII\n\n"
        << "DATASET UNSTRUCTURED_GRID\n";
    ofs << "POINTS " << mesh->V.size() << " double\n";
    for (size_t i = 0; i < mesh->V.size(); i++) {
        ofs << std::fixed << std::setprecision(7) <<  mesh->V.at(i).x << " " <<  mesh->V.at(i).y << " " <<  mesh->V.at(i).z << "\n";
    }
    ofs << "CELLS " << c_indices.size() << " " << 3 * c_indices.size() << std::endl;
    for (size_t i = 0; i < c_indices.size(); i++) {
        auto& e = mesh->E.at(c_indices.at(i));
        if (e.Vids.empty()) continue;
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
    }

    ofs.close();
    ofs.clear();
    ofs.open(in+"_out.vtk");

    c_indices.clear();
    for (auto& f: mesh->F) {
        if (f.Vids.empty() || f.N_Fids.empty()) continue;
        c_indices.push_back(f.id);
    }
    ofs << "# vtk DataFile Version 3.0\n"
        << "Mesh" << ".vtk\n"
        << "ASCII\n\n"
        << "DATASET UNSTRUCTURED_GRID\n";
    ofs << "POINTS " << mesh->V.size() << " double\n";
    for (size_t i = 0; i < mesh->V.size(); i++) {
        ofs << std::fixed << std::setprecision(7) <<  mesh->V.at(i).x << " " <<  mesh->V.at(i).y << " " <<  mesh->V.at(i).z << "\n";
    }
    ofs << "CELLS " << c_indices.size() << " " << 5 * c_indices.size() << std::endl;
    for (size_t i = 0; i < c_indices.size(); i++) {
        auto& f = mesh->F.at(c_indices.at(i));
        ofs << "4 " << f.Vids.at(0) << " " << f.Vids.at(1) << " " << f.Vids.at(2) << " " << f.Vids.at(3) << " " << std::endl;
    }
    ofs << "CELL_TYPES " << c_indices.size() << "\n";
    for (size_t i = 0; i < c_indices.size(); i++) {
        ofs << "9" << std::endl;
    }
}

void SemiGlobalSimplifier::PrototypeResolveGroup(SingularityLink& l1, SingularityLink& l2) {
    // std::cout << "l1: " << mesh->V.at(l1.frontId).N_Vids.size() << " " << mesh->V.at(l1.backId).N_Vids.size() << std::endl;
    // std::cout << "l2: " << mesh->V.at(l2.frontId).N_Vids.size() << " " << mesh->V.at(l2.backId).N_Vids.size() << std::endl;
    // std::cout << "Before l1: " << l1.a + l1.b << " l2: " << l2.a + l2.b << std::endl;
    while (true) {
        // std::cout << "l1.a " << l1.a << " l1.b " << l1.b << std::endl;
        // std::cout << "l1 front: " << mesh->V.at(l1.frontId).N_Vids.size() << " l1 back: " << mesh->V.at(l1.backId).N_Vids.size() << std::endl;
        // std::cout << "l2.a " << l2.a << " l2.b " << l2.b << std::endl;
        // std::cout << "l2 front: " << mesh->V.at(l2.frontId).N_Vids.size() << " l2 back: " << mesh->V.at(l2.backId).N_Vids.size() << std::endl;
        if (mesh->V.at(l1.frontId).N_Vids.size() != 3 && mesh->V.at(l1.frontId).N_Vids.size() != 5) break;
        if (mesh->V.at(l1.backId).N_Vids.size() != 3 && mesh->V.at(l1.backId).N_Vids.size() != 5) break;
        if (mesh->V.at(l2.backId).N_Vids.size() != 3 && mesh->V.at(l2.backId).N_Vids.size() != 5) break;
        if (PrototypeCancelThreeFivePair(l1, l2) > 0) break;
        if (PullSingularity(l1, l2) == -1) break;
        if (PrototypeCancelThreeFivePair(l1, l2) > 0) break;
        if (PullSingularity(l2, l1) == -1) break;
        if (PrototypeCancelThreeFivePair(l1, l2) > 0) break;
    }
    PrototypeCancelThreeFivePair(l1, l2);
    // std::cout << "After l1: " << l1.a + l1.b << " l2: " << l2.a + l2.b << std::endl;
}

int SemiGlobalSimplifier::PrototypeCancelThreeFivePair(SingularityLink& l1, SingularityLink& l2) {
    if (l1.a + l1.b == 1) {
        if (!ValidatePath(l1.linkVids) || !ValidatePath(l2.linkVids)) return 0;
        std::vector<size_t> threeFiveIds;
        if (mesh->V.at(l1.linkVids.at(0)).N_Vids.size() == 3 && mesh->V.at(l1.linkVids.at(1)).N_Vids.size() == 5) {
            threeFiveIds.push_back(l1.linkVids.at(0));
            threeFiveIds.push_back(l1.linkVids.at(1));
        } else if (mesh->V.at(l1.linkVids.at(0)).N_Vids.size() == 5 && mesh->V.at(l1.linkVids.at(1)).N_Vids.size() == 3) {
            threeFiveIds.push_back(l1.linkVids.at(1));
            threeFiveIds.push_back(l1.linkVids.at(0));
        }
        if (threeFiveIds.size() < 2) return 0;
        MovePair(threeFiveIds, l2.linkVids);
        return 1;
    }
    return 0;
}

void SemiGlobalSimplifier::PrototypeC() {
    CheckValidity();

    std::queue<size_t> q;
    BaseComplexQuad bc(*mesh);
    int i = 0;
    bool breakLoop = false;
    int valenceToCancel = 3;
    while (!breakLoop) {
        if (i++ >= iters) breakLoop = true;
        for (auto& v: mesh->V) {
            if (v.type != FEATURE && !v.isBoundary && v.N_Vids.size() == valenceToCancel) {
                if (PrototypeCancelSingularity(v.id, bc) == -1) {
                    valenceToCancel = valenceToCancel == 3 ? 5 : 3;
                    continue;
                }
                // PrototypeCancelSingularity(q.front(), bc);
                // q.push(v.id);
                break;
            }
        }
        valenceToCancel = valenceToCancel == 3 ? 5 : 3;
        // while (!q.empty()) {
        //     int vid = PrototypeCancelSingularity(q.front(), bc);
        //     q.pop();
        //     if (vid != -1 && mesh->V.at(vid).type != FEATURE && !mesh->V.at(vid).isBoundary && (mesh->V.at(vid).N_Vids.size() == 3 || mesh->V.at(vid).N_Vids.size() == 5)) {
        //         std::cout << "vid: " << vid << " " << mesh->V.at(vid).N_Vids.size() << std::endl;
        //         q.push(vid);
        //     }
        //     // break;
        // }
        while(FixValences());
    }
}

int SemiGlobalSimplifier::PrototypeCancelSingularity(size_t vid, BaseComplexQuad& bc) {
    auto& v = mesh->V.at(vid);
    SingularityLink l = PrototypeGetLink(v.id, bc);
    if (l.linkVids.empty()) return -1;
    return PrototypeCancelSingularityPair(l, bc);
}

void SemiGlobalSimplifier::PrototypeD() {
    CheckValidity();
    BaseComplexQuad bc(*mesh);
    std::priority_queue<SingularityLink, std::vector<SingularityLink>, LinkComparator> q;
    for (auto& v: mesh->V) {
        if (v.type == FEATURE || v.isBoundary || (v.N_Vids.size() != 3 && v.N_Vids.size() != 5)) continue;
        SingularityLink l = PrototypeGetLink(v.id, bc);
        if (l.linkVids.empty()) continue;
        // if (PrototypeCheckBoundarySingularity(l.frontId) && mesh->V.at(l.frontId).N_Vids.size() == 3) continue;
        // if (PrototypeCheckBoundarySingularity(l.frontId) && PrototypeCheckBoundarySingularity(l.backId)) continue;
        // if (PrototypeCheckBoundarySingularity(l.frontId) && v.N_Vids.size() == 3) continue;
        // if (!PrototypeCheckBoundarySingularity(l.backId)) continue;
        q.push(l);
    }

    int colorValue = 0;
    int ncolors = 15;
    std::vector<size_t> c_indices;
    std::vector<int> colors;
    int i = 1;
    std::cout << "q size: " << q.size() << std::endl;
    while (!q.empty()) {
        SingularityLink l = q.top();
        // SingularityLink l2 = PrototypeGetLink(l.frontId, bc, l.backId, l.linkEids, false, false);
        // if (i == iters) {
            // std::cout << "l rank: " << l.rank << " a+b: " << l.a+l.b << " " << mesh->V.at(l.frontId).N_Vids.size() << " " << mesh->V.at(l.backId).N_Vids.size() << std::endl;
            PrototypeCancelSingularityPair(l, bc);
            // std::cout << "l2 rank: " << l2.rank << " a+b: " << l2.a+l2.b << std::endl;
            // c_indices.insert(c_indices.end(), l.linkEids.begin(), l.linkEids.end());
            // c_indices.insert(c_indices.end(), l2.linkEids.begin(), l2.linkEids.end());
            // std::vector<int> a(l.linkEids.size(), (colorValue%ncolors));
            // colors.insert(colors.end(), a.begin(), a.end());
            // colorValue += 1;
            // a.clear();
            // a.resize(l2.linkEids.size(), (colorValue%ncolors));
            // colors.insert(colors.end(), a.begin(), a.end());
            // break;
        // }
        // i += 1;
        q.pop();
    }
    while (FixValences());
    /*std::cout << "Writing output file" << std::endl;
    std::ofstream ofs("singularityLinks.vtk");
    ofs << "# vtk DataFile Version 3.0\n"
        << "singularityLinks" << ".vtk\n"
        << "ASCII\n\n"
        << "DATASET UNSTRUCTURED_GRID\n";
    ofs << "POINTS " << mesh->V.size() << " double\n";
    // for (auto& v: source.V) {
    //     if (v.type == FEATURE) c_indices.push_back(v.id);
    // }
    // std::vector<size_t> c_indices = {12, 296};
    // std::cout << c_indices.size() << std::endl;
    for (size_t i = 0; i < mesh->V.size(); i++) {
        ofs << std::fixed << std::setprecision(7) <<  mesh->V.at(i).x << " " <<  mesh->V.at(i).y << " " <<  mesh->V.at(i).z << "\n";
    }
    ofs << "CELLS " << c_indices.size() << " " << 3 * c_indices.size() << std::endl;
    for (size_t i = 0; i < c_indices.size(); i++) {
        auto& e = mesh->E.at(c_indices.at(i));
        if (e.Vids.empty()) continue;
        ofs << "2 " << e.Vids.at(0) << " " << e.Vids.at(1) << std::endl;

        // ofs << "1 " << c_indices.at(i) << std::endl;
    }
    ofs << "CELL_TYPES " << c_indices.size() << "\n";
    for (size_t i = 0; i < c_indices.size(); i++) {
        ofs << "3" << std::endl;
        // ofs << "1" << std::endl;
    }

    ofs << "CELL_DATA " << c_indices.size() << "\n";
    ofs << "SCALARS fixed int\n";
    ofs << "LOOKUP_TABLE default\n";
    for (auto c: colors) {
        ofs << c << "\n";
    }*/
}

void SemiGlobalSimplifier::PrototypeB() {
    CheckValidity();

    BaseComplexQuad bc(*mesh);
    std::priority_queue<SingularityLink, std::vector<SingularityLink>, LinkComparator> q;
    for (auto& v: mesh->V) {
        if (v.type == FEATURE || v.isBoundary || (v.N_Vids.size() != 3 && v.N_Vids.size() != 5)) continue;
        SingularityLink l = PrototypeGetLink(v.id, bc);
        if (!l.linkVids.empty()) q.push(l);
    }
    
    std::cout << q.size() << std::endl;
    int i = 0;
    while (!q.empty()) {
        // if (i >= iters) break;
        SingularityLink l = q.top();
        PrototypeCancelSingularityPair(l, bc);
        q.pop();
        i += 1;
    }
    while(FixValences());
}

int SemiGlobalSimplifier::PrototypeCancelSingularityPair(SingularityLink& l, BaseComplexQuad& bc) {
    auto& front = mesh->V.at(l.frontId);
    auto& back = mesh->V.at(l.backId);

    int nextVid = -1;
    if ((front.N_Vids.size() != 3 && front.N_Vids.size() != 5) || front.type == FEATURE || front.isBoundary) return nextVid;
    if ((back.N_Vids.size() != 3 && back.N_Vids.size() != 5) || back.type == FEATURE || back.isBoundary) return nextVid;

    
    // std::cout << "rank: " << l.rank << " " << l.a + l.b << std::endl;
    // std::cout << "valences: " << mesh->V.at(l.frontId).N_Vids.size() << " " << mesh->V.at(l.backId).N_Vids.size() << std::endl;
    // c_indices.insert(c_indices.end(), l.linkEids.begin(), l.linkEids.end());
    // std::vector<int> a(l.linkEids.size(), (colorValue%ncolors));
    // colors.insert(colors.end(), a.begin(), a.end());
    // colorValue += 1;


    std::vector<size_t> mainPath = l.linkVids;
    if (!ValidatePath(mainPath)) return -1;
    // std::cout << "main path: ";
    // for (auto id: mainPath) std::cout << id << " ";
    // std::cout << std::endl;
    for (int i = 0; i < mainPath.size()-1; i++) {
        size_t toMove = mainPath.at(i);
        if (mesh->V.at(toMove).N_Vids.empty()) return -1;
        // std::cout << "main path step: " << i << " to move: " << toMove << std::endl;
        SingularityLink l2 = PrototypeGetLink(toMove, bc, l.backId, std::vector<size_t>(l.linkEids.begin()+i, l.linkEids.end()), false, false);
        if (l2.linkVids.empty()) return -1;
        // std::cout << "rank: " << l2.rank << " " << l2.a + l2.b << std::endl;
        // std::cout << "valences: " << mesh->V.at(l2.frontId).N_Vids.size() << " " << mesh->V.at(l2.backId).N_Vids.size() << std::endl;
        // int colorValue = 0;
        // int ncolors = 15;
        // std::vector<size_t> c_indices;
        // std::vector<int> colors;
        // c_indices.insert(c_indices.end(), l2.linkEids.begin(), l2.linkEids.end());
        // std::vector<int> a(l2.linkEids.size(), (colorValue%ncolors));
        // colors.insert(colors.end(), a.begin(), a.end());
        // colorValue += 1;
        std::vector<size_t> secondaryPath = l2.linkVids;
        // std::cout << "secondary path: ";
        // for (auto id: secondaryPath) std::cout << id << " ";
        // std::cout << std::endl;        
        size_t sourceDir = mainPath.at(i+1);
        size_t secondaryDir = secondaryPath.at(1);
        std::vector<size_t> threeFiveIds = GetThreeFivePairIds(toMove, sourceDir, secondaryDir);
        if (threeFiveIds.size() < 2) return -1;
        MovePair(threeFiveIds, secondaryPath);
        // break;
        // if (MovePair(threeFiveIds, secondaryPath, true) == -1) return;
    }
    return nextVid;
}

bool SemiGlobalSimplifier::PrototypeCheckBoundarySingularity(size_t vid) {
    auto& v = mesh->V.at(vid);
    for (auto nvid: v.N_Vids) {
        if (mesh->V.at(nvid).type == FEATURE || mesh->V.at(nvid).isBoundary) {
            return true;
        }
    }
    return false;
}

void SemiGlobalSimplifier::PrototypeBoundary(bool checkValence) {
    CheckValidity();

    BaseComplexQuad bc(*mesh);
    std::queue<size_t> Singularities;
    std::vector<bool> isAvailable(mesh->V.size(), false);
    for (auto& v: mesh->V) {
        if (v.type != FEATURE && !v.isBoundary && (v.N_Vids.size() == 3 || v.N_Vids.size() == 5)) {
            Singularities.push(v.id);
            isAvailable.at(v.id) = true;
        }
    }
    while (!Singularities.empty()) {
        size_t vid = Singularities.front();
        Singularities.pop();
        auto& v = mesh->V.at(vid);
        if (v.type == FEATURE || v.isBoundary || (v.N_Vids.size() != 3 && v.N_Vids.size() != 5)) continue;
        PrototypeA(vid, bc);
    }
    while(FixValences());
}

void SemiGlobalSimplifier::PrototypeA(size_t vid, BaseComplexQuad& bc, bool checkValence) {
    auto& v = mesh->V.at(vid);
    auto cmp = [](SingularityLink left, SingularityLink right) {return left.a + left.b > right.a + right.b;};
    auto cmp_pair = [](std::pair<SingularityLink, SingularityLink> left, std::pair<SingularityLink, SingularityLink> right) {
        return (left.first.a + left.first.b + left.second.a + left.second.b > right.first.a + right.first.b + right.second.a + right.second.b);
    };
    std::priority_queue<std::pair<SingularityLink, SingularityLink>, std::vector<std::pair<SingularityLink, SingularityLink>>, decltype(cmp_pair)> pair_q(cmp_pair);
    std::vector<SingularityLink> links = GetLinks(v.id, bc, checkValence);
    std::priority_queue<SingularityLink, std::vector<SingularityLink>, decltype(cmp)> q(cmp);
    for (auto& l: links) {
        q.push(l);
    }
    std::vector<SingularityLink> linkPair;
    bool foundBoundarySingularity = false;
    while (!q.empty()) {
        auto& l = q.top();
        auto& lv = mesh->V.at(l.backId);
        for (auto nvid: lv.N_Vids) {
            if (mesh->V.at(nvid).type == FEATURE || mesh->V.at(nvid).isBoundary) {
                foundBoundarySingularity = true;
                break;
            }
        }
        if (linkPair.empty()) {
            linkPair.push_back(l);
        } else if (foundBoundarySingularity) {
            auto& l2 = linkPair.at(0);
            bool isDisjoint = true;
            for (int i = 1; i < l2.linkVids.size(); i++) {
                if (std::find(l.linkVids.begin()+1, l.linkVids.end(), l2.linkVids.at(i)) != l.linkVids.end()) {
                    isDisjoint = false;
                    break;
                }
            }
            if (isDisjoint) {
                linkPair.push_back(l);
                break;
            }
        }
        q.pop();
    }
    if (linkPair.size() < 2) return;
    auto& l1 = linkPair.at(0);
    auto& l2 = linkPair.at(1);
    GenerateSingularityPair(l1, l2);
    if (l1.a + l1.b > 1) return;
    std::vector<size_t> threeFiveIds;
    if (mesh->V.at(l1.frontId).N_Vids.size() == 3) {
        threeFiveIds.push_back(l1.frontId);
        threeFiveIds.push_back(l1.backId);
    } else {
        threeFiveIds.push_back(l1.backId);
        threeFiveIds.push_back(l1.frontId);
    }
    std::vector<size_t> secondaryPath = l2.linkVids;
    MovePair(threeFiveIds, secondaryPath);
}

void SemiGlobalSimplifier::AlignAndResolveSingularities(bool checkValence) {
    CheckValidity();
    
    BaseComplexQuad bc(*mesh);
    std::queue<size_t> Singularities;
    for (auto& v: mesh->V) {
        if (v.type != FEATURE && !v.isBoundary && v.N_Vids.size() != 4) Singularities.push(v.id);
    }

    std::vector<bool> isAvailable(mesh->V.size(), true);
    while (!Singularities.empty()) {
        size_t vid = Singularities.front();
        Singularities.pop();
        auto& v = mesh->V.at(vid);
        if (v.type == FEATURE || v.isBoundary || !isAvailable.at(v.id) || (v.N_Vids.size() != 3 && v.N_Vids.size() != 5)) continue;
        AlignSingularities(vid, Singularities, isAvailable, bc, checkValence);
    }
    while (ResolveIsolatedSingularities(bc)) {
        std::cout << "resolving ..." << std::endl;
        FixValences();
        // break;
    }
}

bool SemiGlobalSimplifier::ResolveIsolatedSingularities(BaseComplexQuad& bc) {
    bool res = false;
    for (auto& v: mesh->V) {
        if (v.type == FEATURE || v.isBoundary || (v.N_Vids.size() != 3 && v.N_Vids.size() != 5)) continue;
        if (IsPair(v.id)) continue;
        std::deque<size_t> dq;
        dq.push_front(v.id);
        while (!dq.empty()) {
            size_t vid = dq.front();
            dq.pop_front();
            auto& s = mesh->V.at(vid);
            if (s.type == FEATURE || s.isBoundary || (s.N_Vids.size() != 3 && s.N_Vids.size() != 5)) continue;
            if (IsPair(s.id)) continue;
            std::vector<SingularityLink> links = GetLinks(s.id, bc, false);
            int minId = -1;
            int minValue = -1;
            for (int i = 0; i < links.size(); i++) {
                auto& l = links.at(i);
                if (!IsPair(l.backId)) continue;
                if (minId == -1 || l.a + l.b < minValue) {
                    minValue = l.a + l.b;
                    minId = i;
                }
            }
            if (minId > -1) {
                auto& l = links.at(minId);
                std::vector<size_t> threeFiveIds = GetPair(l.backId);
                std::cout << "singularity valence: " << s.N_Vids.size() << std::endl;
                std::cout << "l back id: " << mesh->V.at(l.backId).N_Vids.size() << std::endl;
                std::cout << "threeFiveIds: " << mesh->V.at(threeFiveIds.at(0)).N_Vids.size() << " " << mesh->V.at(threeFiveIds.at(1)).N_Vids.size() << std::endl;
                std::vector<size_t> path = l.linkVids;
                std::cout << "path: " << path.size() << std::endl;
                std::reverse(path.begin(), path.end());
                int newId = MovePair(threeFiveIds, path, true);
                // std::cout << "sid: " << s.id << std::endl;
                std::cout << "l backId: " << l.backId << std::endl;
                std::cout << "newId: " << newId << std::endl;
                std::cout << "threeFiveIds: " << threeFiveIds.at(0) << " " << threeFiveIds.at(1) << std::endl;
                if (newId != -1 && newId != threeFiveIds.at(0) && newId != threeFiveIds.at(1)) {
                    res = true;
                    dq.push_front(newId);
                }
            }
        }
        if (res) break;
    }
    return res;
}

std::vector<size_t> SemiGlobalSimplifier::GetPair(size_t vid) {
    std::vector<size_t> res = {vid};
    auto& v = mesh->V.at(vid);
    int valenceToCheck = v.N_Vids.size() == 3 ? 5 : 3;
    for (auto nvid: v.N_Vids) {
        if (mesh->V.at(nvid).N_Vids.size() == valenceToCheck) {
            res.push_back(nvid);
            break;
        }
    }
    if (mesh->V.at(res.at(0)).N_Vids.size() != 3) {
        std::reverse(res.begin(), res.end());
    }
    return res;
}

bool SemiGlobalSimplifier::IsPair(size_t vid) {
    auto& v = mesh->V.at(vid);
    if (v.N_Vids.size() != 3 && v.N_Vids.size() != 5) return false;
    int valenceToCheck = v.N_Vids.size() == 3 ? 5 : 3;
    for (auto nvid: v.N_Vids) {
        auto& nv = mesh->V.at(nvid);
        if (nv.type != FEATURE && !nv.isBoundary && nv.N_Vids.size() == valenceToCheck) return true;
    }
    return false;
}

void SemiGlobalSimplifier::AlignSingularities(size_t vid, std::queue<size_t>& Singularities, std::vector<bool>& isAvailable, BaseComplexQuad& bc, bool checkValence) {
    CheckValidity();

    auto& v = mesh->V.at(vid);
    std::vector<SingularityLink> links = GetLinks(v.id, bc, checkValence);
    auto cmp = [](SingularityLink left, SingularityLink right) {return left.a + left.b > right.a + right.b;};
    std::priority_queue<SingularityLink, std::vector<SingularityLink>, decltype(cmp)> q(cmp);
    for (auto& l: links) {
        if (!isAvailable.at(l.backId)) continue;
        q.push(l);
    }
    std::vector<SingularityLink> linkPair;
    while (!q.empty()) {
        auto& l = q.top();
        if (linkPair.empty()) {
            linkPair.push_back(l);
        } else {
            auto& l2 = linkPair.at(0);
            bool isDisjoint = true;
            for (int i = 1; i < l2.linkVids.size(); i++) {
                if (std::find(l.linkVids.begin()+1, l.linkVids.end(), l2.linkVids.at(i)) != l.linkVids.end()) {
                    isDisjoint = false;
                    break;
                }
            }
            if (isDisjoint) {
                linkPair.push_back(l);
                break;
            }
        }
        q.pop();
    }
    if (linkPair.size() < 2) return;
    auto& l1 = linkPair.at(0);
    auto& l2 = linkPair.at(1);
        
    GenerateSingularityPair(l1, l2);
    isAvailable.resize(mesh->V.size(), true);
    // Singularities.push(l1.frontId);
    // Singularities.push(l1.backId);
    if (l1.a + l1.b == 1) {
        isAvailable.at(l1.frontId) = false;
        isAvailable.at(l1.backId) = false;
        Singularities.push(l2.backId);
        // std::vector<size_t> threeFiveIds;
        // if (mesh->V.at(l1.frontId).N_Vids.size() == 3) {
        //     threeFiveIds.push_back(l1.frontId);
        //     threeFiveIds.push_back(l1.backId);
        // } else {
        //     threeFiveIds.push_back(l1.backId);
        //     threeFiveIds.push_back(l1.frontId);
        // }
        // std::vector<size_t> secondaryPath = l2.linkVids;
        // MovePair(threeFiveIds, secondaryPath);
    }
}

void SemiGlobalSimplifier::AlignSingularities() {
    CheckValidity();

    BaseComplexQuad bc(*mesh);
    std::vector<bool> vertexVisited(mesh->V.size(), false);
    std::vector<SingularityLink> links;
    std::vector<SingularityLink> DirectPairs;
    size_t sid = 0;
    int it = 0;
    int nSingularities = 0;
    for (auto& v: mesh->V) {
        if (vertexVisited.at(v.id)) continue;
        if (v.isBoundary || v.type == FEATURE) continue;
        if (v.N_Vids.size() == 4) continue;
        links = GetLinks(v.id, bc);
        if (links.size() < 2) continue;
        auto cmp = [](SingularityLink left, SingularityLink right) {return left.a + left.b > right.a + right.b;};
        std::priority_queue<SingularityLink, std::vector<SingularityLink>, decltype(cmp)> q(cmp);
        for (auto& l: links) {
            q.push(l);
        }
        auto& l = q.top();
        if (l.a + l.b == 1 && !vertexVisited.at(l.frontId) && !vertexVisited.at(l.backId)) {
            vertexVisited.at(l.frontId) = true;
            vertexVisited.at(l.backId) = true;
            DirectPairs.push_back(l);
        }
    }
    for (auto& v: mesh->V) {
        // std::cout << "v: " << v.id << std::endl;
        if (vertexVisited.at(v.id)) continue;
        if (v.isBoundary || v.type == FEATURE) continue;
        if (v.N_Vids.size() == 4) continue;
        // if (v.N_Vids.size() != 3) continue;
        // if (it++ == iters) break;
        // std::cout << "iter: " << it++ << std::endl; 
        links = GetLinks(v.id, bc);
        if (links.size() < 2) continue;
        // std::cout << "number links: " << links.size() << std::endl;
        auto cmp = [](SingularityLink left, SingularityLink right) {return left.a + left.b > right.a + right.b;};
        std::priority_queue<SingularityLink, std::vector<SingularityLink>, decltype(cmp)> q(cmp);
        // std::cout << "before link operations" << std::endl;
        for (auto& l: links) {
            // std::cout << "l.a: " << l.a << " l.b: " << l.b << std::endl;
            // std::cout << "l frontId: " << l.frontId << " l backId: " << l.backId << std::endl;
            // std::cout << "l vids: " << l.linkVids.size() << std::endl;
            // std::cout << vertexVisited.size() << std::endl; 
            if (vertexVisited.at(l.backId)) continue;
            q.push(l);
        }
        std::vector<SingularityLink> linkPair;
        // std::cout << "befor queue operations" << std::endl;
        while (!q.empty()) {
            auto& l = q.top();
            // std::cout << "l vids: " << l.linkVids.size() << std::endl;
            if (linkPair.empty()) {
                linkPair.push_back(l);
            } else {
                auto& l2 = linkPair.at(0);
                bool isDisjoint = true;
                for (int i = 1; i < l2.linkVids.size(); i++) {
                    if (std::find(l.linkVids.begin()+1, l.linkVids.end(), l2.linkVids.at(i)) != l.linkVids.end()) {
                        isDisjoint = false;
                        break;
                    }
                }
                if (isDisjoint) {
                    linkPair.push_back(l);
                    break;
                }
            }
            q.pop();
        }
        // std::cout << "link Pair size: " << linkPair.size() << std::endl;
        if (linkPair.size() < 2) continue;
        auto& l1 = linkPair.at(0);
        auto& l2 = linkPair.at(1);
        // std::cout << "l1.a " << l1.a << " l1.b: " << l1.b << std::endl;
        // if (l1.a + l1.b == 1) {
        //     vertexVisited.at(l1.frontId) = true;
        //     vertexVisited.at(l1.backId) = true;
        // } else {
            GenerateSingularityPair(l1, l2);
        // }
        vertexVisited.resize(mesh->V.size(), false);
        // std::cout << "l1.a " << l1.a << " l1.b: " << l1.b << std::endl;
        // std::cout << "l1 frontId: " << l1.frontId << " l1 backId: " << l1.backId << std::endl;
        // std::cout << vertexVisited.size() << std::endl; 
        if (l1.a + l1.b == 1 && !vertexVisited.at(l1.backId) && !vertexVisited.at(l1.backId)) {
            vertexVisited.at(l1.frontId) = true;
            vertexVisited.at(l1.backId) = true;
            DirectPairs.push_back(l1);
        }
        // links = linkPair;
        // links.clear();
        // links.push_back(l1);
        // links.push_back(l2);
        // sid = v.id;
        // break;
    }
    std::cout << "Got " << DirectPairs.size() << " pairs" << std::endl;

    /*std::cout << "Singularity Links for vertex: " << sid << std::endl;
    
    int colorValue = 0;
    int ncolors = 15;
    std::vector<size_t> c_indices;
    std::vector<int> colors;
    for (auto& l: links) {
        std::cout << "a: " << l.a << " b: " << l.b << std::endl;
        c_indices.insert(c_indices.end(), l.linkVids.begin(), l.linkVids.end());
        std::vector<int> a(l.linkVids.size(), (colorValue%ncolors));
        colors.insert(colors.end(), a.begin(), a.end());
        colorValue += 1;
    }
    
    std::cout << "Writing output file" << std::endl;
    std::ofstream ofs("singularityLinks.vtk");
    ofs << "# vtk DataFile Version 3.0\n"
        << "singularityLinks" << ".vtk\n"
        << "ASCII\n\n"
        << "DATASET UNSTRUCTURED_GRID\n";
    ofs << "POINTS " << mesh->V.size() << " double\n";
    // for (auto& v: source.V) {
    //     if (v.type == FEATURE) c_indices.push_back(v.id);
    // }
    // std::vector<size_t> c_indices = {12, 296};
    // std::cout << c_indices.size() << std::endl;
    for (size_t i = 0; i < mesh->V.size(); i++) {
        ofs << std::fixed << std::setprecision(7) <<  mesh->V.at(i).x << " " <<  mesh->V.at(i).y << " " <<  mesh->V.at(i).z << "\n";
    }
    ofs << "CELLS " << c_indices.size() << " " << 2 * c_indices.size() << std::endl;
    for (size_t i = 0; i < c_indices.size(); i++) {
        // auto& e = mesh->E.at(c_indices.at(i));
        // if (e.Vids.empty()) continue;
        // ofs << "2 " << e.Vids.at(0) << " " << e.Vids.at(1) << std::endl;

        ofs << "1 " << c_indices.at(i) << std::endl;
    }
    ofs << "CELL_TYPES " << c_indices.size() << "\n";
    for (size_t i = 0; i < c_indices.size(); i++) {
        // ofs << "3" << std::endl;
        ofs << "1" << std::endl;
    }

    ofs << "CELL_DATA " << c_indices.size() << "\n";
    ofs << "SCALARS fixed int\n";
    ofs << "LOOKUP_TABLE default\n";
    for (auto c: colors) {
        ofs << c << "\n";
    }*/
}

void SemiGlobalSimplifier::GenerateSingularityPair(SingularityLink& l1, SingularityLink& l2) {
    // std::cout << "Inside GenerateSingularityPair" << std::endl;
    while (true) {
        // std::cout << "l1.a " << l1.a << " l1.b " << l1.b << std::endl;
        // std::cout << "l1 front: " << mesh->V.at(l1.frontId).N_Vids.size() << " l1 back: " << mesh->V.at(l1.backId).N_Vids.size() << std::endl;
        // std::cout << "l2.a " << l2.a << " l2.b " << l2.b << std::endl;
        // std::cout << "l2 front: " << mesh->V.at(l2.frontId).N_Vids.size() << " l2 back: " << mesh->V.at(l2.backId).N_Vids.size() << std::endl;
        if (l1.a + l1.b == 1) break;
        if (mesh->V.at(l1.frontId).N_Vids.size() != 3 && mesh->V.at(l1.frontId).N_Vids.size() != 5) break;
        if (mesh->V.at(l1.backId).N_Vids.size() != 3 && mesh->V.at(l1.backId).N_Vids.size() != 5) break;
        if (mesh->V.at(l2.backId).N_Vids.size() != 3 && mesh->V.at(l2.backId).N_Vids.size() != 5) break;
        if (PullSingularity(l1, l2) == -1) break;
        if (PullSingularity(l2, l1) == -1) break;
    }
    // std::cout << "**********************************" << std::endl;
}

int SemiGlobalSimplifier::PullSingularity(SingularityLink& l1, SingularityLink& l2) {
    // PrototypeSaveMesh(l1, l2, "Before");
    std::cout << "Inside PullSingularity" << std::endl;
    std::vector<size_t> mainPath = l1.linkVids;
    std::vector<size_t> secondaryPath = l2.linkVids;
    if (!ValidatePath(mainPath) || !ValidatePath(secondaryPath)) return -1;
    std::cout << "mainPath: ";
    for (auto vid: mainPath) {
        std::cout << vid << " ";
    }
    std::cout  << std::endl;
    
    std::cout << "secondaryPath: ";
    for (auto vid: secondaryPath) {
        std::cout << vid << " ";
    }
    std::cout  << std::endl;

    std::vector<int> Rots = GetTraverseInfo(l1, l2);
    if (Rots.at(0) == 0) return -1;
    std::cout << "Got traverse info" << std::endl;
    std::cout << "Rots: ";
    for (auto r: Rots) {
        std::cout << r << " ";
    }
    std::cout  << std::endl;
    // return -1;
    auto& toMove = mesh->V.at(mainPath.at(0));
    size_t sourceDir = mainPath.at(1);
    size_t secondaryDir = secondaryPath.at(1);

    std::cout << "toMove: " << toMove.id << " sourceDir: " << sourceDir << " secondaryDir: " << secondaryDir << std::endl;

    std::vector<size_t> threeFiveIds = GetThreeFivePairIds(toMove.id, sourceDir, secondaryDir);
    if (threeFiveIds.size() < 2) return -1;
    std::cout << "threeID: " << threeFiveIds.at(0) << " " << mesh->V.at(threeFiveIds.at(0)).N_Vids.size() << " fiveID: " << threeFiveIds.at(1) << " " << mesh->V.at(threeFiveIds.at(1)).N_Vids.size() << std::endl;
    auto& sourceDirV = mesh->V.at(sourceDir);
    // size_t startId;
    // if (sourceDirV.N_Vids.size() == 3) {
    //     startId = threeFiveIds.at(1);
    // } else if (sourceDirV.N_Vids.size() == 5) {
    //     size_t startE = 0;
    //     for (auto eid: sourceDirV.N_Eids) {
    //         auto& e = mesh->E.at(eid);
    //         if (e.Vids.at(0) == mainPath.at(2) || e.Vids.at(1) == mainPath.at(2)) {
    //             startE = eid;
    //             break;
    //         }
    //     }
    //     for (int i = 0; i < nRot; i++) {
    //         auto& e = mesh->E.at(startE);
    //         for (auto fid: e.N_Fids) {
    //             auto& f = mesh->F.at(fid);
    //             size_t neid = mu->GetDifference(mu->GetIntersection(sourceDirV.N_Eids, f.Eids), std::vector<size_t>{startE}).at(0);
    //             auto& ne = mesh->E.at(neid);
    //             size_t evid = e.Vids.at(0) == sourceDirV.id ? e.Vids.at(1) : e.Vids.at(0);
    //             int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), sourceDirV.id));
    //             if (f.Vids.at((idx+3)%f.Vids.size()) == evid) {
    //                 startE = neid;
    //                 startId = evid;
    //             }
    //         }
    //     }
    // }
    // if (std::find(sourceDirV.N_Vids.begin(), sourceDirV.N_Vids.end(), threeFiveIds.at(0)) != sourceDirV.N_Vids.end()) {
    //     startId = threeFiveIds.at(0);
    // } else if (std::find(sourceDirV.N_Vids.begin(), sourceDirV.N_Vids.end(), threeFiveIds.at(1)) != sourceDirV.N_Vids.end()) {
    //     startId = threeFiveIds.at(1);
    // }
    
    // std::shared_ptr<ThreeFivePair> tfp = std::make_shared<ThreeFivePair>(*mesh, *mu, *smoother, threeFiveIds.at(0), threeFiveIds.at(1));
    // bool skipCheck = false;
    // for (int i = 1; i < secondaryPath.size(); i++) {
    //     if (mesh->V.at(tfp->GetPairIds().at(1)).N_Vids.size() == 6) skipCheck = true;
    //     tfp->Move(secondaryPath.at(i), skipCheck);
    //     // return -1;
    // }
    // if (MovePair(threeFiveIds, secondaryPath) == -1) return -1;
    MovePair(threeFiveIds, secondaryPath);

    std::cout << "After moving singularity" << std::endl;
    l1.linkVids = std::vector<size_t>(mainPath.begin()+1, mainPath.end());
    l1.frontId = l1.linkVids.front();
    l1.backId = l1.linkVids.back();

    std::cout << "l1 vids: ";
    for (auto vid: l1.linkVids) {
        std::cout << vid << " ";
    }
    std::cout  << std::endl;

    l1.a -= 1;
    if (l1.a == 0) {
        l1.a = l1.b;
        l1.b = 0;
        // Rots.at(0) = Rots.at(0) - 1;
    }

    // std::vector<size_t> newPath = TraversePath(l1.linkVids.at(1), l1.linkVids.at(0), Rots);
    l2.linkVids = TraversePath(l1.linkVids.at(1), l1.linkVids.at(0), Rots);
    l2.frontId = l2.linkVids.front();
    l2.backId = l2.linkVids.back();

    std::cout << "l2 vids: ";
    for (auto vid: l2.linkVids) {
        std::cout << vid << " ";
    }
    std::cout  << std::endl;
    // return -1;

    /*std::vector<size_t> newPath1 = {sourceDir, startId};

    std::cout << "startId: " << startId << " " << mesh->V.at(startId).N_Vids.size() << std::endl;
    int vid1 = newPath1.at(0);
    int vid2 = newPath1.at(1);
    for (int i = 1; i < l2.a; i++) {
        std::cout << "vid1: " << vid1 << " " << mesh->V.at(vid1).N_Vids.size() << " vid2: " << vid2 << " " << mesh->V.at(vid2).N_Vids.size() << std::endl;
        auto& v1 = mesh->V.at(vid1);
        auto& v2 = mesh->V.at(vid2);
        if (v2.N_Vids.size() != 4) continue;
        for (auto nvid: v2.N_Vids) {
            if (nvid == vid1) continue;
            if (mu->GetIntersection(mesh->V.at(nvid).N_Fids, v1.N_Fids).empty()) {
                newPath1.push_back(nvid);
                vid1 = vid2;
                vid2 = nvid;
            }
        }
    }
    std::cout << "newPath1: ";
    for (auto vid: newPath1) {
        std::cout << vid << " ";
    }
    std::cout  << std::endl;
    auto& v1 = mesh->V.at(vid1);
    auto& v2 = mesh->V.at(vid2);
    std::vector<int> vids;
    std::vector<size_t> newPath2;
    std::vector<size_t> newPath3;
    for (auto nvid: v2.N_Vids) {
        if (nvid == vid1) continue;
        if (!mu->GetIntersection(mesh->V.at(nvid).N_Fids, v1.N_Fids).empty()) {
            vids.push_back(nvid);
        }
    }
    std::cout << "vids: " << vids.size() << std::endl;
    vid1 = vid2;
    int vid3 = vids.at(0);
    int vid4 = vids.at(1);
    for (int i = 0; i < l2.b; i++) {
        auto& v1 = mesh->V.at(vid1);
        auto& v2 = mesh->V.at(vid2);
        auto& v3 = mesh->V.at(vid3);
        auto& v4 = mesh->V.at(vid4);
        if (v3.N_Vids.size() == 4) {
            for (auto nvid: v3.N_Vids) {
                if (nvid == vid1) continue;
                if (mu->GetIntersection(mesh->V.at(nvid).N_Fids, v1.N_Fids).empty()) {
                    newPath2.push_back(nvid);
                    vid1 = vid3;
                    vid3 = nvid;
                }
            }
        }
        if (v4.N_Vids.size() == 4) {
            for (auto nvid: v4.N_Vids) {
                if (nvid == vid2) continue;
                if (mu->GetIntersection(mesh->V.at(nvid).N_Fids, v2.N_Fids).empty()) {
                    newPath3.push_back(nvid);
                    vid2 = vid4;
                    vid4 = nvid;
                }
            }
        }
    }

    int vertexToCheck = mesh->V.at(l2.frontId).N_Vids.size() == 3 ? 5 : 3;
    if (mesh->V.at(newPath2.back()).N_Vids.size() == vertexToCheck && newPath2.size() == l2.b) {
        newPath1.insert(newPath1.end(), newPath2.begin(), newPath2.end());
    } else if (mesh->V.at(newPath3.back()).N_Vids.size() == vertexToCheck && newPath3.size() == l2.b) {
        newPath1.insert(newPath1.end(), newPath3.begin(), newPath3.end());
    }
    if (l2.linkVids.size() == newPath1.size()) {
        l2.linkVids = newPath1;
        l2.frontId = l2.linkVids.front();
        l2.backId = l2.linkVids.back();
    }*/
    // PrototypeSaveMesh(l1, l2, "After");
    if (l1.a + l1.b == 1) return -1;
    return 0;
}

int SemiGlobalSimplifier::MovePair(std::vector<size_t> threeFiveIds, std::vector<size_t>& secondaryPath, bool checkValence) {
    if (mesh->V.at(threeFiveIds.at(0)).N_Vids.size() != 3 && mesh->V.at(threeFiveIds.at(1)).N_Vids.size() != 5) return -1;
    int valenceToCheck = mesh->V.at(secondaryPath.back()).N_Vids.size();
    std::shared_ptr<ThreeFivePair> tfp = std::make_shared<ThreeFivePair>(*mesh, *mu, *smoother, threeFiveIds.at(0), threeFiveIds.at(1));
    bool skipCheck = false;
    for (int i = 1; i < secondaryPath.size(); i++) {
        bool skipCheck = i < secondaryPath.size() -1 && mesh->V.at(tfp->GetPairIds().at(1)).N_Vids.size() == 6 ? true : false;
        // if (i < secondaryPath.size() -1 && mesh->V.at(tfp->GetPairIds().at(1)).N_Vids.size() == 6) {
        //     skipCheck = true;
        // } else {
        //     skipCheck = false;
        // }
        tfp->Move(secondaryPath.at(i), skipCheck);
        // skipCheck = false;
    }
    if (checkValence) {
        std::vector<size_t> ids = tfp->GetPairIds();
        if (mesh->V.at(ids.at(0)).N_Vids.size() == valenceToCheck) {
            return ids.at(0);
        } else if (mesh->V.at(ids.at(1)).N_Vids.size() == valenceToCheck) {
            return ids.at(1);
        }    
    }

    return -1;
}

std::vector<int> SemiGlobalSimplifier::GetTraverseInfo(SingularityLink& l1, SingularityLink& l2) {
    size_t start = l1.linkVids.at(1);
    std::vector<size_t> path = l2.linkVids;
    std::vector<int> rots;

    size_t prev = start;
    size_t current = 0;
    size_t next = 0;

    for (int i = 0; i < path.size() - 1; i++) {
        current = path.at(i);
        next = path.at(i+1);
        // std::cout << i << " prev: " << prev << " current: " << current << " next: " << next << std::endl;
        auto& v = mesh->V.at(current);
        // for (auto id: v.N_Eids) {
        //     auto& e = mesh->E.at(id);
        //     std::cout << e.id << ": " << e.Vids.at(0) << " " << e.Vids.at(1) << std::endl;
        // }
        int eid = 0;
        for (int j = 0; j < v.N_Eids.size(); j++) {
            auto& e = mesh->E.at(v.N_Eids.at(j));
            if (e.Vids.at(0) == prev || e.Vids.at(1) == prev) {
                eid = e.id;
                break;
            }
        }
        int nRot = 0;
        for (int j = 0; j < v.N_Eids.size(); j++) {
            auto& e = mesh->E.at(eid);
            // std::cout << "edge: " << e.id << " " << e.Vids.at(0) << " " << e.Vids.at(1) << std::endl;
            if (e.Vids.at(0) == next || e.Vids.at(1) == next) break;
            for (auto fid: e.N_Fids) {
                auto& f = mesh->F.at(fid);
                // std::cout << "face: " << f.Vids.at(0) << " " << f.Vids.at(1) << " " << f.Vids.at(2) << " " << f.Vids.at(3) << std::endl;
                // std::cout << "face edges: " << f.Eids.at(0) << " " << f.Eids.at(1) << " " << f.Eids.at(2) << " " << f.Eids.at(3) << std::endl;
                int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), current));
                // std::cout << "idx: " << idx << " " << (idx+3)%f.Vids.size() << std::endl;
                if (f.Vids.at((idx+1)%f.Vids.size()) == prev) {
                    eid = mu->GetDifference(mu->GetIntersection(v.N_Eids, f.Eids), std::vector<size_t>{e.id}).at(0);
                    prev = f.Vids.at((idx+3)%f.Vids.size());
                    break;
                }
            }
            nRot += 1;
        }
        rots.push_back(nRot);
        prev = current;
    }
    if (l1.a == 1 && l1.b > 0) {
        current = l1.linkVids.at(1);
        prev = l1.linkVids.at(2);
        next = l1.linkVids.at(0);

        auto& v = mesh->V.at(current);
        int eid = 0;
        for (int j = 0; j < v.N_Eids.size(); j++) {
            auto& e = mesh->E.at(v.N_Eids.at(j));
            if (e.Vids.at(0) == prev || e.Vids.at(1) == prev) {
                eid = e.id;
                break;
            }
        }
        int nRot = 0;
        for (int j = 0; j < v.N_Eids.size(); j++) {
            auto& e = mesh->E.at(eid);
            if (e.Vids.at(0) == next || e.Vids.at(1) == next) break;
            for (auto fid: e.N_Fids) {
                auto& f = mesh->F.at(fid);
                int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), current));
                if (f.Vids.at((idx+1)%f.Vids.size()) == prev) {
                    eid = mu->GetDifference(mu->GetIntersection(v.N_Eids, f.Eids), std::vector<size_t>{e.id}).at(0);
                    prev = f.Vids.at((idx+3)%f.Vids.size());
                    break;
                }
            }
            nRot += 1;
        }
        rots.at(0) = (rots.at(0) + nRot) - 2;
    }
    return rots;
}

std::vector<size_t> SemiGlobalSimplifier::TraversePath(size_t prev, size_t current, std::vector<int> Rots) {
    std::vector<size_t> res;
    for (auto nRot: Rots) {
        // std::cout << nRot << std::endl;
        size_t eid = 0;
        res.push_back(current);
        auto& v = mesh->V.at(current);
        for (auto id: v.N_Eids) {
            auto& e = mesh->E.at(id);
            if (e.Vids.at(0) == prev || e.Vids.at(1) == prev) {
                eid = id;
                break;
            }
        }
        for (int i = 0; i < nRot; i++) {
            auto& e = mesh->E.at(eid);
            for (auto fid: e.N_Fids) {
                auto& f = mesh->F.at(fid);
                int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), current));
                if (f.Vids.at((idx+1)%f.Vids.size()) == prev) {
                    eid = mu->GetDifference(mu->GetIntersection(v.N_Eids, f.Eids), std::vector<size_t>{e.id}).at(0);
                    prev = f.Vids.at((idx+3)%f.Vids.size());
                    break;
                }
            }
        }
        size_t temp = current;
        current = prev;
        prev = temp;
        // auto& e = mesh->E.at(eid);
        // prev = current;
        // current = e.Vids.at(0) == current ? e.Vids.at(1) : e.Vids.at(0);
    }
    res.push_back(current);
    return res;
}

void SemiGlobalSimplifier::ResolveSingularities() {
    CheckValidity();

    // Smooth(*mesh, std::vector<size_t>{});
    clock_t start = clock();
    // Smooth();
    BaseComplexQuad bc(*mesh);

    std::vector<size_t> Singularities;
    int nSingularities = 0;
    int nThreeSingularities = 0;
    int nFiveSingularities = 0;
    PARALLEL_FOR_BEGIN(mesh->V.size()) {
        auto& v = mesh->V.at(i);
    // for (auto& v: mesh->V) {
        if (v.N_Vids.size() > 0 && v.N_Fids.size() > 0 && v.N_Vids.size() != 4 && !v.isBoundary) nSingularities += 1;
        if (v.N_Vids.size() > 0 && v.N_Fids.size() > 0 && v.N_Vids.size() != 3 && v.isBoundary) nSingularities += 1;
        if ((v.N_Vids.size() == 3 || v.N_Vids.size() == 5) && v.type != FEATURE && !v.isBoundary) {
            Singularities.push_back(v.id);
            nThreeSingularities += (int) (v.N_Vids.size() == 3);
            nFiveSingularities += (int) (v.N_Vids.size() == 5);
        }
    // }
    } PARALLEL_FOR_END();

    std::cout << "SINGULARITIES: " << nSingularities << " DETECTED SINGULARITIES: " << Singularities.size() << std::endl;
    std::cout << "THREE SINGULARITIES: " << nThreeSingularities << " FIVE SINGULARITIES: " << nFiveSingularities << std::endl;
    return;
    GetSingularityGroups(Singularities, bc);

    std::cout << "END OF SIMPLIFICATION LOOP" << std::endl;
    // int it = 0;
    // for (auto sid: Singularities) {
    //     int mSize = mesh->V.size();
    //     auto& s = mesh->V.at(sid);
    //     std::cout << "Resolving Singularity: " << s.id << std::endl;
    //     // if ((s.N_Vids.size() != 3 && s.N_Vids.size() != 5) || s.type == FEATURE || s.isBoundary) continue;
    //     int id = ResolveSingularity(sid, bc);
    //     while (id != -1) {
    //         id = ResolveSingularity(id, bc);
    //         // std::cout << "NEW ID: " << id << std::endl;
    //     }
    //     // if (mesh->V.size() != mSize) {
    //     //     Smooth();
    //     // }
    //     it += 1;
    //     if (it == iters) {
    //         break;
    //     }
    // }
    // Smooth();
    clock_t end = clock();
    double time = (double) (end-start) / CLOCKS_PER_SEC;
    std::cout << "execution took: " << time << " seconds" << std::endl;
}

std::vector<SingularityLink> SemiGlobalSimplifier::GetLinks(size_t sid, BaseComplexQuad& bc, bool checkValence) {
    auto& s = mesh->V.at(sid);
    std::vector<SingularityLink> links = TraceSingularityLinks(s, bc);
    std::vector<SingularityLink> tempLinks;
    tempLinks.insert(tempLinks.end(), links.begin(), links.end());
    links = SelectLinks(links, s.N_Vids.size(), checkValence);
    for (auto& l: links) {
        l.a = l.linkEids.size();
        l.b = 0;
    }
    for (auto l: tempLinks) {
        std::vector<SingularityLink> newLinks = GetCrossLinks(l, bc, checkValence);
        links.insert(links.end(), newLinks.begin(), newLinks.end());    
    }
    return links;
}

std::vector<SingularityLink> SemiGlobalSimplifier::SelectLinks(std::vector<SingularityLink> links, int valence, bool checkValence) {
    std::vector<SingularityLink> res;
    for (auto l: links) {
        int valenceToCheck = valence == 3 ? 5 : 3;
        auto& back = mesh->V.at(l.backId);
        // if (l.frontId == l.backId || back.N_Vids.size() != valenceToCheck || back.type == FEATURE || back.isBoundary || !ValidateLink(l)) continue;
        if (l.frontId == l.backId || back.type == FEATURE || back.isBoundary) continue;
        if (doesCrossBoundary(l.linkVids, true)) continue;
        // if (doesCrossBoundary(l.linkVids, true) || !ValidateLink(l)) continue;
        // if (!ValidateLink(l)) continue;
        if (checkValence && back.N_Vids.size() != valenceToCheck) continue;
        res.push_back(l);
    }
    return res;
}

void SemiGlobalSimplifier::GetSingularityGroups(std::vector<size_t> Singularities, BaseComplexQuad& bc) {
    auto cmp = [](SingularityLink left, SingularityLink right) {return left.a + left.b > right.a + right.b;};
    // auto cmp = [](SingularityLink left, SingularityLink right) {return 1;};
    std::priority_queue<SingularityLink, std::vector<SingularityLink>, decltype(cmp)> q1(cmp);
    std::mutex mtx;
    
    std::vector<SingularityLink> tempLinks;
    std::vector<bool> SingularitiesToAvoid(mesh->V.size(), false);
    PARALLEL_FOR_BEGIN(Singularities.size()) {
        auto& s = mesh->V.at(Singularities.at(i));
    // for (auto sid: Singularities) {
    //     auto& s = mesh->V.at(sid);
        // int valenceToCheck = s.N_Vids.size() == 3 ? 5 : 3;
        if ((s.N_Vids.size() == 3 && s.N_Vids.size() != 5) || s.type == FEATURE || s.isBoundary) continue;
        std::vector<SingularityLink> links = GetLinks(s.id, bc);
        // std::vector<SingularityLink> links = TraceSingularityLinks(s, bc);
        // {
        //     std::lock_guard<std::mutex> lock(mtx);    
        //     tempLinks.insert(tempLinks.end(), links.begin(), links.end());
        // }
        // links = SelectLinks(links);
        for (auto& l: links) {
            {
                std::lock_guard<std::mutex> lock(mtx);    
                q1.push(l);
            }
        }
    // }
    } PARALLEL_FOR_END();
    // for (auto l: tempLinks) {
    //     std::vector<SingularityLink> newLinks = GetCrossLinks(l, bc);
    //     for (auto newL: newLinks) q.push(newL);    
    // }
    
    std::priority_queue<SingularityLink, std::vector<SingularityLink>, decltype(cmp)> q(cmp);
    std::vector<bool> selected(mesh->V.size(), false);
    while (!q1.empty()) {
        auto l = q1.top();
        if (!selected.at(l.frontId) && !selected.at(l.backId)) {
            q.push(l);
            selected.at(l.frontId) = true;
            selected.at(l.backId) = true;
        }
        q1.pop();
    }

    // std::cout << q.size() << std::endl;

    /*int it = 0;
    while (!q.empty()) {
        // std::cout << q.size() << std::endl; 
        auto& l = q.top();
        // std::cout << l.a << " " << l.b << std::endl;
        auto& front = mesh->V.at(l.frontId);
        auto& back = mesh->V.at(l.backId);
        std::vector<size_t> mainPath = l.linkVids;
        if (!ValidatePath(mainPath)) {
            q.pop();
            continue;
        }
        int id = mainPath.at(0);
        // std::cout << "NEW MOVING ITERATION" << std::endl;
        while (id != -1) {
            // std::cout << "to move: " << id << std::endl;
            int offset = std::distance(mainPath.begin(), std::find(mainPath.begin(), mainPath.end(), id));
            // std::cout << "Offset: " << offset << " " << " id:" << id << std::endl;
            std::vector<size_t> secPath = GetSecondaryPath(offset, mainPath, bc);
            // std::cout << "mainPath: " << mainPath.size() << std::endl;
            // for (auto id: mainPath) {
            //     std::cout << id << " ";
            // }
            // std::cout << std::endl;
            // std::cout << "secPath: " << secPath.size() << std::endl;
            // for (auto id: secPath) {
            //     std::cout << id << " ";
            // }
            // std::cout << std::endl;
            if (!ValidatePath(mainPath) || !ValidatePath(secPath) || mainPath.empty() || secPath.empty()) break;
            if (mainPath.size() > 2 && secPath.size() == 2) {   
            // if (mainPath.size() > secPath.size()) {
                std::swap(mainPath, secPath);
            }
            id = MoveSingularity(offset, mainPath, secPath);
        }
        // std::cout << "****************************" << std::endl;

        q.pop();
        // it += 1;
        // if (it == iters) {
        //     break;
        // }
    }*/
    
    int colorValue = 0;
    int ncolors = 15;
    std::vector<size_t> c_indices;
    std::vector<int> colors;
    while (!q.empty()) {
        auto& l = q.top();
        q.pop();
        c_indices.insert(c_indices.end(), l.linkEids.begin(), l.linkEids.end());
        std::vector<int> a(l.linkEids.size(), (colorValue%ncolors));
        colors.insert(colors.end(), a.begin(), a.end());
        colorValue += 1;
    }
    
    std::cout << "Writing output file" << std::endl;
    std::ofstream ofs("singularityLinks.vtk");
    ofs << "# vtk DataFile Version 3.0\n"
        << "singularityLinks" << ".vtk\n"
        << "ASCII\n\n"
        << "DATASET UNSTRUCTURED_GRID\n";
    ofs << "POINTS " << mesh->V.size() << " double\n";
    // for (auto& v: source.V) {
    //     if (v.type == FEATURE) c_indices.push_back(v.id);
    // }
    // std::vector<size_t> c_indices = {12, 296};
    // std::cout << c_indices.size() << std::endl;
    for (size_t i = 0; i < mesh->V.size(); i++) {
        ofs << std::fixed << std::setprecision(7) <<  mesh->V.at(i).x << " " <<  mesh->V.at(i).y << " " <<  mesh->V.at(i).z << "\n";
    }
    ofs << "CELLS " << c_indices.size() << " " << 3 * c_indices.size() << std::endl;
    for (size_t i = 0; i < c_indices.size(); i++) {
        auto& e = mesh->E.at(c_indices.at(i));
        if (e.Vids.empty()) continue;
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
    }
}

std::vector<size_t> SemiGlobalSimplifier::GetSecondaryPath(int offset, std::vector<size_t>& mainPath, BaseComplexQuad& bc) {
    mainPath = std::vector<size_t>(mainPath.begin()+offset, mainPath.end());
    size_t frontId = mainPath.at(0);
    size_t backId = mainPath.at(mainPath.size()-1);
    std::vector<SingularityLink> links_a = GetLinks(frontId, bc);
    std::vector<SingularityLink> links_b = GetLinks(backId, bc);

    int rank = -1;
    std::vector<size_t> res;
    for (auto l: links_a) {
        if (l.backId == frontId || l.backId == backId) continue;
        if (!mu->GetIntersection(std::vector<size_t>(mainPath.begin()+1,mainPath.end()-1), std::vector<size_t>(l.linkVids.begin()+1, l.linkVids.end()-1)).empty()) continue;
        if (rank == -1 || rank > l.a + l.b) {
            rank = l.a + l.b;
            res = l.linkVids;
        } 
    }
    for (auto l: links_b) {
        if (l.backId == frontId || l.backId == backId) continue;
        if (!mu->GetIntersection(std::vector<size_t>(mainPath.begin()+1,mainPath.end()-1), std::vector<size_t>(l.linkVids.begin()+1, l.linkVids.end()-1)).empty()) continue;
        if (rank == -1 || rank > l.a + l.b) {
            rank = l.a + l.b;
            res = l.linkVids;
        } 
    }
    if (rank > -1 && res.at(0) == backId) {
        std::reverse(mainPath.begin(), mainPath.end());
    }
    return res;
}

std::vector<SingularityLink> SemiGlobalSimplifier::GetCrossLinks(SingularityLink& l, BaseComplexQuad& bc, bool checkValence) {
    std::vector<SingularityLink> newLinks;
    std::vector<size_t> verticesToCheck;
    for (int i = 1; i < l.linkVids.size() - 1; i++) {
        auto& v = mesh->V.at(l.linkVids.at(i));
        if (v.isBoundary || v.type == FEATURE) break;
        verticesToCheck.push_back(l.linkVids.at(i));
    }
    std::mutex mtx;
    PARALLEL_FOR_BEGIN(verticesToCheck.size()) {
        auto& v = mesh->V.at(verticesToCheck.at(i));
    // for (auto vid: verticesToCheck) {
    //     auto& v = mesh->V.at(vid);
        std::vector<size_t> a;
        for (auto vid: l.linkVids) {
            a.push_back(vid);
            if (vid == v.id) break;
        }
        std::vector<size_t> b;
        for (auto eid: l.linkEids) {
            auto& e = mesh->E.at(eid);
            b.push_back(eid);
            if (e.Vids.at(1) == v.id || e.Vids.at(0) == v.id) break;
        }
        std::vector<SingularityLink> tempLinks = TraceSingularityLinks(v, bc);
        tempLinks = SelectLinks(tempLinks, mesh->V.at(l.frontId).N_Vids.size(), checkValence);
        for (auto tempL: tempLinks) {
            if (tempL.backId == l.frontId || tempL.backId == l.backId) continue;
            SingularityLink newL;
            newL.linkVids = a;
            newL.linkVids.insert(newL.linkVids.end(), tempL.linkVids.begin()+1, tempL.linkVids.end());
            newL.linkEids = b;
            newL.linkEids.insert(newL.linkEids.end(), tempL.linkEids.begin(), tempL.linkEids.end());
            newL.frontId = newL.linkVids.front();
            newL.backId = newL.linkVids.back();
            newL.a = b.size();
            newL.b = tempL.linkEids.size();
            newL.rot = PrototypeGetRotations(tempL.linkVids.front(), tempL.linkEids.front(), b.back());
            {
                std::lock_guard<std::mutex> lock(mtx);
                newLinks.push_back(newL);
            }
        }
    // }
    } PARALLEL_FOR_END(); 
    return newLinks;
}

bool SemiGlobalSimplifier::ValidateLink(SingularityLink& l) {
    // int start = mesh->V.at(l.linkVids.at(0)).N_Vids.size() == 4 ? 0 : 1;
    for (int i = 1; i < l.linkVids.size() - 1; i++) {
        auto& v = mesh->V.at(l.linkVids.at(i));
        for (auto nvid: v.N_Vids) {
            if (nvid == l.frontId || nvid == l.backId) continue;
            auto& nv = mesh->V.at(nvid);
            if (nv.N_Vids.size() != 4 || nv.type == FEATURE || nv.isBoundary) return false;
            // if (nv.N_Vids.size() != 4 && nv.type != FEATURE && !nv.isBoundary) return false;
            // if (nv.N_Vids.size() != nv.idealValence && (nv.type == FEATURE || nv.isBoundary)) return false;
        }
    }
    return true;
}

bool SemiGlobalSimplifier::ValidatePath(std::vector<size_t> p) {
    for (auto id: p) if (mesh->V.at(id).N_Vids.empty()) return false;
    for (int i = 1; i < p.size(); i++) {
        auto& v1 = mesh->V.at(p.at(i-1));
        auto& v2 = mesh->V.at(p.at(i));
        if (v1.N_Vids.empty() || v2.N_Vids.empty()) return false;
        if (std::find(v1.N_Vids.begin(), v1.N_Vids.end(), v2.id) == v1.N_Vids.end()) return false;
        if (std::find(v2.N_Vids.begin(), v2.N_Vids.end(), v1.id) == v2.N_Vids.end()) return false;
    }
    return true;
}

int SemiGlobalSimplifier::ResolveSingularity(size_t sid, BaseComplexQuad& bc) {
    auto& s = mesh->V.at(sid);
    if ((s.N_Vids.size() != 3 && s.N_Vids.size() != 5) || s.type == FEATURE || s.isBoundary) return -1;
    std::mutex mtx;
    auto cmp = [](SingularityLink left, SingularityLink right) {return left.a + left.b > right.a + right.b;};
    std::priority_queue<SingularityLink, std::vector<SingularityLink>, decltype(cmp)> q(cmp);
    
    int valenceToCheck = s.N_Vids.size() == 3 ? 5 : 3;
    std::vector<SingularityLink> links = TraceSingularityLinks(s, bc);
    for (auto& l: links) {
        // if (mesh->V.at(l.backId).N_Vids.size() != valenceToCheck) continue;
        
        std::vector<size_t> verticesToCheck;
        bool skipLink = false;
        for (int i = 1; i < l.linkVids.size() - 1; i++) {
            auto& v = mesh->V.at(l.linkVids.at(i));
            for (auto nvid: v.N_Vids) {
                if (mesh->V.at(nvid).isBoundary || mesh->V.at(nvid).type == FEATURE) {
                    skipLink = true;
                    break;
                }
            }
            if (v.isBoundary || v.type == FEATURE) break;
            verticesToCheck.push_back(l.linkVids.at(i));
        }
        if (!skipLink && (mesh->V.at(l.backId).N_Vids.size() == 3 || mesh->V.at(l.backId).N_Vids.size() == 5)) {
            l.a = l.linkEids.size();
            l.b = 0;
            q.push(l);
        }
        
        PARALLEL_FOR_BEGIN(verticesToCheck.size()) {
            auto& v = mesh->V.at(verticesToCheck.at(i));
            std::vector<size_t> a;
            for (auto vid: l.linkVids) {
                a.push_back(vid);
                if (vid == v.id) break;
            }
            std::vector<size_t> b;
            for (auto eid: l.linkEids) {
                auto& e = mesh->E.at(eid);
                b.push_back(eid);
                if (e.Vids.at(1) == v.id || e.Vids.at(0) == v.id) break;
            }
            std::vector<SingularityLink> tempLinks = TraceSingularityLinks(v, bc);
            for (auto& tempL: tempLinks) {
                // if (mesh->V.at(tempL.backId).N_Vids.size() != valenceToCheck) continue;
                if (mesh->V.at(tempL.backId).N_Vids.size() != 3 && mesh->V.at(tempL.backId).N_Vids.size() != 5) continue;
                // if (mesh->V.at(tempL.backId).type == FEATURE || mesh->V.at(tempL.backId).isBoundary) continue;
                if (tempL.backId == l.frontId || tempL.backId == l.backId || doesCrossBoundary(tempL.linkVids, true)) continue;
                SingularityLink newL;
                newL.linkVids = a;
                newL.linkVids.insert(newL.linkVids.end(), tempL.linkVids.begin()+1, tempL.linkVids.end());
                newL.linkEids = b;
                newL.linkEids.insert(newL.linkEids.end(), tempL.linkEids.begin(), tempL.linkEids.end());
                newL.frontId = newL.linkVids.front();
                newL.backId = newL.linkVids.back();
                newL.a = b.size();
                newL.b = tempL.linkEids.size();
                {
                    std::lock_guard<std::mutex> lock(mtx);
                    q.push(newL);
                }
            }
        } PARALLEL_FOR_END();
    }
    if (q.size() > 1) {
        SingularityGroup sg;
        std::priority_queue<SingularityLink, std::vector<SingularityLink>, decltype(cmp)> q2 = q;
        bool foundFirstLink = false;
        while (!q2.empty()) {
            auto l = q2.top();
            if (l.frontId != l.backId && mesh->V.at(l.backId).N_Vids.size() == valenceToCheck) {
                sg.l1 = l;
                foundFirstLink = true;
                break;
            }
            q2.pop();
        }
        if (foundFirstLink) {
            bool foundSecondLink = false;
            while (!q.empty()) {
                auto l = q.top();
                // if (mesh->V.at(l.frontId).N_Vids.size() != 3 && !IsExclusive(l, sg.l1)) {
                //     q.pop();
                //     continue;
                // }
                // if (l.backId != sg.l1.backId && IsExclusive(l, sg.l1)) {
                if (l.frontId != l.backId && l.backId != sg.l1.backId) {
                    sg.l2 = l;
                    foundSecondLink = true;
                    break;
                }
                q.pop();
            }
            if (foundSecondLink) {
                /*std::cout << "Singularity Links for vertex: " << sid << std::endl;
                std::cout << "Link 1: " << sg.l1.frontId << " " << mesh->V.at(sg.l1.frontId).N_Vids.size() << " " << sg.l1.backId << " " << mesh->V.at(sg.l1.backId).N_Vids.size() << std::endl; 
                std::cout << "Link 2: " << sg.l2.frontId << " " << mesh->V.at(sg.l2.frontId).N_Vids.size() << " " << sg.l2.backId << " " << mesh->V.at(sg.l2.backId).N_Vids.size() << std::endl; 
                std::cout << "*******************************************************" << std::endl;
                std::vector<size_t> c_indices;
                c_indices.insert(c_indices.end(), sg.l1.linkEids.begin(), sg.l1.linkEids.end());
                c_indices.insert(c_indices.end(), sg.l2.linkEids.begin(), sg.l2.linkEids.end());  
                std::cout << "Writing links output file" << std::endl;
                std::string fname = "SingularityLinks_"+std::to_string(sid)+".vtk";
                std::ofstream ofs(fname);
                ofs << "# vtk DataFile Version 3.0\n"
                    << fname << "\n"
                    << "ASCII\n\n"
                    << "DATASET UNSTRUCTURED_GRID\n";
                ofs << "POINTS " << mesh->V.size() << " double\n";
                // for (auto& v: source.V) {
                //     if (v.type == FEATURE) c_indices.push_back(v.id);
                // }
                // std::vector<size_t> c_indices = {12, 296};
                // std::cout << c_indices.size() << std::endl;
                for (size_t i = 0; i < mesh->V.size(); i++) {
                    ofs << std::fixed << std::setprecision(7) <<  mesh->V.at(i).x << " " <<  mesh->V.at(i).y << " " <<  mesh->V.at(i).z << "\n";
                }
                ofs << "CELLS " << c_indices.size() << " " << 3 * c_indices.size() << std::endl;
                for (size_t i = 0; i < c_indices.size(); i++) {
                    auto& e = mesh->E.at(c_indices.at(i));
                    if (e.Vids.empty()) continue;
                    ofs << "2 " << e.Vids.at(0) << " " << e.Vids.at(1) << std::endl;
                }
                ofs << "CELL_TYPES " << c_indices.size() << "\n";
                for (size_t i = 0; i < c_indices.size(); i++) {
                    ofs << "3" << std::endl;
                }
                c_indices.clear();
                for (size_t i = 0; i < mesh->E.size(); i++) {
                    auto& e = mesh->E.at(i);
                    if (e.Vids.empty()) continue;
                    c_indices.push_back(e.id);
                }
                std::cout << "Writing mesh output file" << std::endl;
                std::string fname2 = "Mesh_"+std::to_string(sid)+".vtk";
                std::ofstream ofs2(fname2);
                ofs2 << "# vtk DataFile Version 3.0\n"
                    << fname2 << "\n"
                    << "ASCII\n\n"
                    << "DATASET UNSTRUCTURED_GRID\n";
                ofs2 << "POINTS " << mesh->V.size() << " double\n";
                // for (auto& v: source.V) {
                //     if (v.type == FEATURE) c_indices.push_back(v.id);
                // }
                // std::vector<size_t> c_indices = {12, 296};
                // std::cout << c_indices.size() << std::endl;
                for (size_t i = 0; i < mesh->V.size(); i++) {
                    ofs2 << std::fixed << std::setprecision(7) <<  mesh->V.at(i).x << " " <<  mesh->V.at(i).y << " " <<  mesh->V.at(i).z << "\n";
                }
                ofs2 << "CELLS " << c_indices.size() << " " << 3 * c_indices.size() << std::endl;
                for (size_t i = 0; i < c_indices.size(); i++) {
                    auto& e = mesh->E.at(c_indices.at(i));
                    ofs2 << "2 " << e.Vids.at(0) << " " << e.Vids.at(1) << std::endl;
                }
               ofs2 << "CELL_TYPES " << c_indices.size() << "\n";
                for (size_t i = 0; i < c_indices.size(); i++) {
                    ofs2 << "3" << std::endl;
                }*/
                return MoveSingularity(sg);
                // return -1;
            }
        }
    }
    /*if (q.size() > 1) {
        SingularityGroup sg;
        sg.l1 = q.top();
        q.pop();
        bool foundLink = false;
        while (!q.empty()) {
            auto& l = q.top();
            if (!mu->GetIntersection(std::vector<size_t>(l.linkVids.begin()+1, l.linkVids.end()-1), std::vector<size_t>(sg.l1.linkVids.begin()+1, sg.l1.linkVids.end()-1)).empty()) {
                q.pop();
                continue;
            }
            sg.l2 = l;
            foundLink = true;
            q.pop();
            break;
        }
        if (foundLink) {
            return MoveSingularity(sg);
        }
    }*/
    return -1;
}

int SemiGlobalSimplifier::MoveSingularity(int offset, std::vector<size_t>& mainPath, std::vector<size_t>& secondaryPath) {
    auto& toMove = mesh->V.at(mainPath.at(0));
    auto& source = mesh->V.at(mainPath.at(mainPath.size()-1));
    auto& secondary = mesh->V.at(secondaryPath.at(secondaryPath.size()-1));

    if ((toMove.N_Vids.size() != 3 && toMove.N_Vids.size() != 5) || toMove.type == FEATURE || toMove.isBoundary) return -1;
    if ((source.N_Vids.size() != 3 && source.N_Vids.size() != 5) || source.type == FEATURE || source.isBoundary) return -1;
    if ((secondary.N_Vids.size() != 3 && secondary.N_Vids.size() != 5) || secondary.type == FEATURE || secondary.isBoundary) return -1;

    // std::cout << "Offset: " << offset << std::endl;
    size_t sourceDir = mainPath.at(1);
    size_t secondaryDir = secondaryPath.at(1);

    // std::cout << "toMove: " << toMove.id << " " << toMove.N_Vids.size() << std::endl;
    // std::cout <<  " source: " << source.id << " " << source.N_Vids.size() << std::endl;
    // std::cout << " secondary: " << secondary.id << " " << secondary.N_Vids.size() << std::endl;
    // std::cout << "sourceDir: " << sourceDir << " " << mesh->V.at(sourceDir).N_Vids.size() << std::endl;
    // std::cout << " secondaryDir: " << secondaryDir << " " << mesh->V.at(secondaryDir).N_Vids.size() << std::endl;

    int res = -1;

    size_t threeId, fiveId;
    if (source.id == sourceDir) {
        threeId = toMove.N_Vids.size() == 3 ? toMove.id : source.id;
        fiveId = toMove.N_Vids.size() == 5 ? toMove.id : source.id;
    } else {
        std::vector<size_t> threeFiveIds = GetThreeFivePairIds(toMove.id, sourceDir, secondaryDir);
        if (threeFiveIds.size() < 2) return res;
        threeId = threeFiveIds.at(0);
        fiveId = threeFiveIds.at(1);
        res = sourceDir;
    }
    // std::cout << "threeID: " << threeId << ": " << mesh->V.at(threeId).N_Vids.size() << " fiveID: " << fiveId << ": " << mesh->V.at(fiveId).N_Vids.size() << std::endl;
    std::shared_ptr<ThreeFivePair> tfp = std::make_shared<ThreeFivePair>(*mesh, *mu, *smoother, threeId, fiveId);
    bool skipCheck = false;
    for (int i = 1; i < secondaryPath.size(); i++) {
        // bool skipCheck = tfp->GetPairIds().at(1) == secondaryPath.at(secondaryPath.size()-1) && mesh->V.at(tfp->GetPairIds().at(1)).N_Vids.size() == 6 ? true : false;
        if (mesh->V.at(tfp->GetPairIds().at(1)).N_Vids.size() == 6) skipCheck = true;
        tfp->Move(secondaryPath.at(i), skipCheck);
    }
    return res;
}

int SemiGlobalSimplifier::MoveSingularity(SingularityGroup& sg) {
    auto& l1 = sg.l1;
    auto& l2 = sg.l2;

    auto& toMove = mesh->V.at(l1.frontId);
    auto& source = mesh->V.at(l1.backId);
    auto& secondary = mesh->V.at(l2.backId);
    size_t sourceDir = l1.linkVids.at(1);
    size_t secondaryDir = l2.linkVids.at(1);
    std::vector<size_t> secondaryPath(l2.linkVids.begin()+1, l2.linkVids.end());

    int res = -1;

    size_t threeId, fiveId;
    if (source.id == sourceDir) {
        threeId = toMove.N_Vids.size() == 3 ? toMove.id : source.id;
        fiveId = toMove.N_Vids.size() == 5 ? toMove.id : source.id;
    } else {
        std::vector<size_t> threeFiveIds = GetThreeFivePairIds(toMove.id, sourceDir, secondaryDir);
        if (threeFiveIds.size() < 2) return res;
        threeId = threeFiveIds.at(0);
        fiveId = threeFiveIds.at(1);
        res = sourceDir;
    }

    std::shared_ptr<ThreeFivePair> tfp = std::make_shared<ThreeFivePair>(*mesh, *mu, *smoother, threeId, fiveId);
    for (int i = 0; i < secondaryPath.size(); i++) {
        tfp->Move(secondaryPath.at(i));
    }
    return res;
}

bool SemiGlobalSimplifier::IsExclusive(SingularityLink& l1, SingularityLink& l2) {
    for (int i = 0; i < l1.linkEids.size(); i++) {
        for (int j = 0; j < l2.linkEids.size(); j++) {
            auto& e1 = mesh->E.at(l1.linkEids.at(i));
            auto& e2 = mesh->E.at(l2.linkEids.at(j));
            if (!mu->GetIntersection(e1.N_Fids, e2.N_Fids).empty()) return false;
        }
    }
    return true;
}

bool SemiGlobalSimplifier::IsExclusive(size_t vid, std::vector<size_t> a, std::vector<size_t> b) {
    auto& e1 = mesh->E.at(a.at(0));
    auto& e2 = mesh->E.at(b.at(0));
    if (e1.id == e2.id) return false;
    std::vector<size_t> commonF = mu->GetIntersection(e1.N_Fids, e2.N_Fids);
    if (!commonF.empty() && a.size() > 1 && b.size() > 1) {
        auto& f = mesh->F.at(commonF.at(0));
        if (mu->GetIntersection(f.Eids, std::vector<size_t>{a.at(0), a.at(1), b.at(0), b.at(1)}).size() != 2) return false;
    }
    return true;
}

void SemiGlobalSimplifier::GetSingularityPairs() {
    BaseComplexQuad bc(*mesh);
    std::vector<SingularityLink> singularityLinks;
    std::vector<std::vector<size_t>> singularityMap(mesh->V.size());
    SetSingularityLinks(singularityLinks, singularityMap, bc);
    std::cout << "Singularity Links: " << singularityLinks.size() << std::endl;
    std::vector<SingularityGroup> groups;
    SelectSingularityGroups(groups, singularityLinks, singularityMap);
    ResolveSingularityGroups(groups, bc);
    // for (auto& sg: groups) {
    //     std::cout << "rank: " << sg.rank << " " << sg.l1.id << " " << sg.l2.id << " ";
    //     std::cout << mesh->V.at(sg.l1.backId).N_Vids.size() << " ";
    //     std::cout << mesh->V.at(sg.l1.frontId).N_Vids.size() << " ";
    //     std::cout << mesh->V.at(sg.l2.backId).N_Vids.size() << " " << std::endl;;
    // }
    
    std::cout << "Singularity groups: " << groups.size() << std::endl;
    Smooth();
    return;
    clock_t start = clock();
    // auto cmp = [](SingularityLink left, SingularityLink right) {return left.a + left.b < right.a + right.b;};
    // std::vector<std::multiset<SingularityLink, decltype(cmp)>> singularityMap(mesh->V.size(), std::multiset<SingularityLink, decltype(cmp)>(cmp));
    std::mutex mtx;
    std::vector<size_t> Singularities;
    
    // std::priority_queue<SingularityLink, std::vector<SingularityLink>, decltype(cmp)> q(cmp);

    /*PARALLEL_FOR_BEGIN(mesh->V.size()) {
        auto& v = mesh->V.at(i);
        if (v.type == FEATURE || v.isBoundary || !v.isSingularity || v.N_Fids.empty()) continue;
        {
            std::lock_guard<std::mutex> lock(mtx);
            Singularities.push_back(v.id);
        }
    } PARALLEL_FOR_END();*/

    // int sid = 0;
    /*PARALLEL_FOR_BEGIN(mesh->V.size()) {
        // auto& v = mesh->V.at(Singularities.at(i));
        auto& v = mesh->V.at(i);
    // for (auto& v: mesh->V) {
        std::vector<SingularityLink> links = TraceSingularityLinks(v, bc);
        for (auto& l: links) {
            if (mesh->V.at(l.frontId).N_Vids.size() == 4 || mesh->V.at(l.backId).N_Vids.size() == 4) continue;
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
    } PARALLEL_FOR_END();*/
    /*std::vector<SingularityLink> newLinks;
    PARALLEL_FOR_BEGIN(singularityLinks.size()) {
        auto& l = singularityLinks.at(i);
        if (doesCrossBoundary(l.linkVids, true)) continue;
        SingularityLink newL;
        for (int j = i+1; j < singularityLinks.size(); j++) {
            auto& l2 = singularityLinks.at(j);
            if (doesCrossBoundary(l2.linkVids, true)) continue;
            if (l.frontId == l2.frontId || l.frontId == l2.backId || l.backId == l2.frontId || l.backId == l2.backId) continue;
            std::vector<size_t> a = mu->GetIntersectionParallel(l.linkVids, l2.linkVids);
            if (a.size() != 1) continue;
            int vid = a.at(0);
            std::vector<size_t> b;
            for (auto eid: l.linkEids) {
                auto& e = mesh->E.at(eid);
                b.push_back(eid);
                if (e.Vids.at(1) == vid || e.Vids.at(0) == vid) break;
            }
            std::vector<size_t> c;
            int it = 0;
            for (auto eid: l2.linkEids) {
                auto& e = mesh->E.at(eid);
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
                auto& e = mesh->E.at(el);
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
        auto& v = mesh->V.at(Singularities.at(i));
        for (auto vid: Singularities) {
            if (v.id == vid) continue;
            auto& v2 = mesh->V.at(vid);
            for (auto lid: singularityMap.at(v.id)) {
                auto& l = singularityLinks.at(lid);
                if (doesCrossBoundary(l.linkVids, true)) continue;
                int rank1 = -1; 
                int rank2 = -1;
                SingularityLink newL;
                for (auto lid2: singularityMap.at(v2.id)) {
                    auto& l2 = singularityLinks.at(lid2);
                    if (l.frontId == l2.frontId || l.frontId == l2.backId || l.backId == l2.frontId || l.backId == l2.backId) continue;
                    std::vector<size_t> a = mu->GetIntersectionParallel(l.linkVids, l2.linkVids);
                    if (a.size() != 1) continue;
                    int vid = a.at(0);
                    std::vector<size_t> b;
                    for (auto eid: l.linkEids) {
                        auto& e = mesh->E.at(eid);
                        b.push_back(eid);
                        if (e.Vids.at(1) == vid || e.Vids.at(0) == vid) break;
                    }
                    std::vector<size_t> c;
                    int it = 0;
                    for (auto eid: l2.linkEids) {
                        auto& e = mesh->E.at(eid);
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
                        auto& e = mesh->E.at(el);
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
        auto& v = mesh->V.at(Singularities.at(i));
        for (int j = 0; j < singularityMap.at(v.id).size(); j++) {
            std::multiset<SingularityLink, decltype(cmp)>::iterator it = singularityMap.at(v.id).begin();
            std::advance(it, j);
            auto& l = *it;
            size_t vid = l.frontId == v.id ? l.backId : l.frontId;
            auto& v2 = mesh->V.at(vid);
            if (v.N_Vids.size() > 5 || v2.N_Vids.size() > 5) continue;
            if (doesCrossBoundary(l.linkVids, true)) continue;
            bool skipSameSingularityType = v.N_Vids.size() == v2.N_Vids.size() ? true : false;    
            for (int k = j + 1; k < singularityMap.at(v.id).size(); k++) {
                it = singularityMap.at(v.id).begin();
                std::advance(it, k);
                auto& l2 = *it;
                size_t vid2 = l2.frontId == v.id ? l2.backId : l2.frontId;
                auto& v3 = mesh->V.at(vid2);
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
            auto& v2 = mesh->V.at(vid);
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

    /*int offset = 0;
    std::vector<SingularityLink> newLinks;
    PARALLEL_FOR_BEGIN(singularityMap.size()) {
    // for (int s = 0; s < singularityMap.size(); s++) {
        int n = singularityMap.at(i).size();
        for (int j = 0; j < n; j++) {
            auto& l = singularityLinks.at(singularityMap.at(i).at(j));
            if (doesCrossBoundary(l.linkVids, true)) continue;
            int rank1 = -1; 
            int rank2 = -1;
            for (auto& l2: singularityLinks) {
                if (l.id == l2.id) continue;
                if (l.frontId == l2.frontId || l.frontId == l2.backId || l.backId == l2.frontId || l.backId == l2.backId) continue;
                std::vector<size_t> a = mu->GetIntersectionParallel(l.linkVids, l2.linkVids);
                if (a.size() != 1) continue;
                int vid = a.at(0);
                std::vector<size_t> b;
                for (auto eid: l.linkEids) {
                    auto& e = mesh->E.at(eid);
                    b.push_back(eid);
                    if (e.Vids.at(1) == vid || e.Vids.at(0) == vid) break;
                }
                std::vector<size_t> c;
                int it = 0;
                for (auto eid: l2.linkEids) {
                    auto& e = mesh->E.at(eid);
                    c.push_back(eid);
                    it += 1;
                    if (e.Vids.at(1) == vid || e.Vids.at(0) == vid) break;
                }
                std::vector<size_t> c1(l2.linkEids.begin()+it, l2.linkEids.end());
                std::reverse(c.begin(), c.end());

                SingularityLink newL;
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
                    auto& e = mesh->E.at(el);
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
    //     auto& v = mesh->V.at(Singularities.at(i));
    // // for (auto& v: mesh->V) {
    //     bool skipSameSingularityType = false;
    //     int it = 0;
    //     for (auto lid: singularityMap.at(v.id)) {
    //         auto& l = singularityLinks.at(lid);
    //     }
    // // }
    // } PARALLEL_FOR_END();

    
    std::priority_queue<SingularityGroup, std::vector<SingularityGroup>, decltype(cmp)> finalGroups(cmp);
    PARALLEL_FOR_BEGIN(singularityLinks.size()) {
        auto& l = singularityLinks.at(i);
    // for (auto& l: singularityLinks) {
        auto& v1 = mesh->V.at(l.frontId);
        auto& v2 = mesh->V.at(l.backId);
        if ((v1.N_Vids.size() == 3 && v2.N_Vids.size() == 5) || (v1.N_Vids.size() == 5 && v2.N_Vids.size() == 3)) {
            if (doesCrossBoundary(l.linkVids, true)) continue;
            std::vector<int> linksToCheck(singularityMap.at(l.frontId).begin(), singularityMap.at(l.frontId).end());
            linksToCheck.insert(linksToCheck.end(), singularityMap.at(l.backId).begin(), singularityMap.at(l.backId).end());
            int minLinkId = -1;
            int rank = -1;
            for (auto lid: linksToCheck) {
                auto& l2 = singularityLinks.at(lid);
                if ((l.frontId == l2.frontId && l.backId == l2.backId) || (l.frontId == l2.backId && l.backId == l2.frontId)) continue;
                if (mesh->V.at(l2.frontId).N_Vids.size() > 5 || mesh->V.at(l2.backId).N_Vids.size() > 5) continue;
                if (doesCrossBoundary(l2.linkVids, true)) continue;
                if (!mu->GetIntersection(l.linkEids, l2.linkEids).empty()) continue;
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
                    // q.push(sg);
                    finalGroups.push(sg);
                }
            }
        }
    // }
    } PARALLEL_FOR_END();
    std::cout << "total singularity links: " << q.size() << std::endl;

    /*int qSize = q.size();
    // int qSize = ThreeFivePairs.size();
    std::cout << "Identified " << qSize << " groups" << std::endl;*/

    // std::vector<bool> selected(mesh->V.size(), false);
    // std::vector<SingularityGroup> finalGroups;
    // std::priority_queue<SingularityGroup, std::vector<SingularityGroup>, decltype(cmp)> finalGroups(cmp);
    // while (!q.empty()) {
    //     auto sg = q.top();
    //     auto l1 = sg.l1;
    //     auto l2 = sg.l2;
    //     if (selected.at(l1.frontId) || selected.at(l1.backId) || selected.at(l2.frontId) || selected.at(l2.backId)) {
    //         q.pop();
    //         continue;
    //     }
    //     selected.at(l1.frontId) = true;
    //     selected.at(l1.backId) = true;
    //     selected.at(l2.frontId) = true;
    //     selected.at(l2.backId) = true;
    //     finalGroups.push(sg);
    //     q.pop();
    // }*/
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
            if (mesh->V.at(l2.frontId).N_Vids.size() > 5 || mesh->V.at(l2.backId).N_Vids.size() > 5) continue;
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
            if (mesh->V.at(l2.frontId).N_Vids.size() > 5 || mesh->V.at(l2.backId).N_Vids.size() > 5) continue;
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

    // int numSing = 0;
    // for (auto& v: mesh->V) {
    //     if (v.isSingularity && v.N_Vids.size() > 0 && !selected.at(v.id)) {
    //         numSing += 1;
    //     }
    // }
    //     std::cout << numSing << " non selected singularities" << std::endl;

    /*clock_t end = clock();
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
    std::vector<bool> selected(mesh->V.size(), false);
    while (!finalGroups.empty()) {
        auto sg = finalGroups.top();
        auto l1 = sg.l1;
        auto l2 = sg.l2;
        if (selected.at(l1.frontId) || selected.at(l1.backId) || selected.at(l2.frontId) || selected.at(l2.backId)) {
            finalGroups.pop();
            continue;
        }
        if (l1.linkVids.empty() || l2.linkVids.empty() || l1.linkEids.empty() || l2.linkEids.empty()) {
            finalGroups.pop();
            // std::cout << "Queue size: " << finalGroups.size() << std::endl;
            continue;
        }
        std::vector<int> a(l1.linkEids.size(), (colorValue%ncolors));
        std::vector<int> b(l2.linkEids.size(), (colorValue%ncolors));
        c_indices.insert(c_indices.end(), l1.linkEids.begin(), l1.linkEids.end());
        c_indices.insert(c_indices.end(), l2.linkEids.begin(), l2.linkEids.end());
        colors.insert(colors.end(), a.begin(), a.end());
        colors.insert(colors.end(), b.begin(), b.end());
        colorValue += 1;
        size_t toMoveId = l1.frontId == l2.frontId || l1.frontId == l2.backId ? l1.frontId : l1.backId;
        size_t sourceId = l1.frontId == toMoveId ? l1.backId : l1.frontId;
        size_t secondaryId = l2.frontId == toMoveId ? l2.backId : l2.frontId;
        if (mesh->V.at(toMoveId).N_Vids.size() == 4 || mesh->V.at(sourceId).N_Vids.size() ==4 || mesh->V.at(secondaryId).N_Vids.size() == 4) {
            finalGroups.pop();
            continue;
        }
        if (mesh->V.at(toMoveId).N_Vids.empty() || mesh->V.at(sourceId).N_Vids.empty() || mesh->V.at(secondaryId).N_Vids.empty()) {
            finalGroups.pop();
            continue;
        }
        std::vector<size_t> mainPath = l1.linkVids;
        std::vector<size_t> secondaryPath = l2.linkVids;
        if (sourceId == mainPath.at(0)) {
            std::reverse(mainPath.begin(), mainPath.end());
        }
        if (secondaryId == secondaryPath.at(0)) {
            std::reverse(secondaryPath.begin(), secondaryPath.end());
        }
        
        std::cout << "rank: " << sg.rank << std::endl;
        for (int i = 1; i < mainPath.size(); i++) {
            std::cout << "secondary path: ";
            for (auto vid: secondaryPath) {
                std::cout << vid << " ";
            }
            std::cout << std::endl;
            if (mesh->V.at(toMoveId).N_Vids.empty() || mesh->V.at(sourceId).N_Vids.empty() || mesh->V.at(secondaryId).N_Vids.empty()) break;
            size_t newSecondary = MoveSingularities(toMoveId, sourceId, secondaryId, mainPath.at(i), secondaryPath.at(1), std::vector<size_t>(secondaryPath.begin()+1, secondaryPath.end()));
            if (i == mainPath.size()-1) continue;
            secondaryId = newSecondary;
            if (!secondaryId) break;
            std::cout << "new secondary: " << mesh->V.at(secondaryId).N_Vids.size() << std::endl;
            if (mesh->V.at(secondaryId).N_Vids.size() == 4) break;
            std::cout << "toMove: " << mesh->V.at(toMoveId).N_Vids.size() << " source: " << mesh->V.at(sourceId).N_Vids.size() << " secondary: " << mesh->V.at(secondaryId).N_Vids.size() << std::endl; 
            
            SetSecondaryPath(secondaryId, toMoveId, sourceId, secondaryPath, bc);
        }
        if (sg.rank > 6) {
            finalGroups.pop();
            break;
        }
        /*if (mainPath.size() == 2) {
            // std::cout << "l1 front: " << l1.frontId << " " << mesh->V.at(l1.frontId).N_Vids.size() << std::endl;
            // std::cout << "l1 back: " << l1.backId << " " << mesh->V.at(l1.backId).N_Vids.size() << std::endl;
            // std::cout << "l2 front: " << l2.frontId << " " << mesh->V.at(l2.frontId).N_Vids.size() << std::endl;
            // std::cout << "l2 back: " << l2.backId << " " << mesh->V.at(l2.backId).N_Vids.size() << std::endl;
            auto& toMove = mesh->V.at(toMoveId);
            auto& source = mesh->V.at(sourceId);
            auto& secondary = mesh->V.at(secondaryId);
            // std::cout << "toMove: " << toMove.id << std::endl;
            // std::cout << "source: " << source.id << std::endl;
            // std::cout << "secondary: " << secondary.id << std::endl;
            // std::cout << "mainPath, ";
            // for (auto id: mainPath) {
            //     std::cout << id << " ";
            // }
            // std::cout << std::endl;
            // std::cout << "secondaryPath, ";
            // for (auto id: secondaryPath) {
            //     std::cout << id << " ";
            // }
            // std::cout << std::endl;
            size_t threeId = toMove.N_Vids.size() == 3 ? toMoveId : sourceId;
            size_t fiveId = toMove.N_Vids.size() == 5 ? toMoveId : sourceId;                
            std::shared_ptr<ThreeFivePair> tfp = std::make_shared<ThreeFivePair>(*mesh, *mu, *smoother, threeId, fiveId);
            for (int i = 1; i < secondaryPath.size(); i++) {
                size_t dest = secondaryPath.at(i);
                // std::cout << "moving 3-5 pair to " << dest << std::endl;
                tfp->Move(dest);
            }
            // std::cout << "Done with moving pair" << std::endl;
            // std::vector<int> a(l1.linkEids.size(), (colorValue%ncolors));
            // std::vector<int> b(l2.linkEids.size(), (colorValue%ncolors));
            // c_indices.insert(c_indices.end(), l1.linkEids.begin(), l1.linkEids.end());
            // c_indices.insert(c_indices.end(), l2.linkEids.begin(), l2.linkEids.end());
            // colors.insert(colors.end(), a.begin(), a.end());
            // colors.insert(colors.end(), b.begin(), b.end());
            // colorValue += 1;    
        } else {
            std::cout << "l1 front: " << l1.frontId << " " << mesh->V.at(l1.frontId).N_Vids.size() << std::endl;
            std::cout << "l1 back: " << l1.backId << " " << mesh->V.at(l1.backId).N_Vids.size() << std::endl;
            std::cout << "l2 front: " << l2.frontId << " " << mesh->V.at(l2.frontId).N_Vids.size() << std::endl;
            std::cout << "l2 back: " << l2.backId << " " << mesh->V.at(l2.backId).N_Vids.size() << std::endl;
            auto& toMove = mesh->V.at(toMoveId);
            auto& source = mesh->V.at(sourceId);
            auto& secondary = mesh->V.at(secondaryId);
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
            std::vector<size_t> threeFiveIds = GetThreeFivePairIds(toMove.id, mainPath.at(1), secondaryPath.at(1));
        }
        selected.at(l1.frontId) = true;
        selected.at(l1.backId) = true;
        selected.at(l2.frontId) = true;
        selected.at(l2.backId) = true;
        finalGroups.pop();
        // break;
        // std::cout << "Queue size: " << finalGroups.size() << std::endl;
    }
    std::cout << "End of GetSingularityPairs function" << std::endl;
    std::cout << c_indices.size() << std::endl;
    // while (!finalGroups.empty()) {
    //     auto& sg = finalGroups.top();
    //     finalGroups.pop();
    //     auto& l1 = sg.l1;
    //     auto& l2 = sg.l2;
    //     std::cout << "l1: " << mesh->V.at(l1.frontId).N_Vids.size() << " " << mesh->V.at(l1.backId).N_Vids.size() << std::endl;
    //     std::cout << "l2: " << mesh->V.at(l2.frontId).N_Vids.size() << " " << mesh->V.at(l2.backId).N_Vids.size() << std::endl;
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
    
    std::cout << "Writing output file" << std::endl;
    std::ofstream ofs("singularityLinks.vtk");
    ofs << "# vtk DataFile Version 3.0\n"
        << "singularityLinks" << ".vtk\n"
        << "ASCII\n\n"
        << "DATASET UNSTRUCTURED_GRID\n";
    ofs << "POINTS " << mesh->V.size() << " double\n";
    // for (auto& v: source.V) {
    //     if (v.type == FEATURE) c_indices.push_back(v.id);
    // }
    // std::vector<size_t> c_indices = {12, 296};
    // std::cout << c_indices.size() << std::endl;
    for (size_t i = 0; i < mesh->V.size(); i++) {
        ofs << std::fixed << std::setprecision(7) <<  mesh->V.at(i).x << " " <<  mesh->V.at(i).y << " " <<  mesh->V.at(i).z << "\n";
    }
    ofs << "CELLS " << c_indices.size() << " " << 3 * c_indices.size() << std::endl;
    for (size_t i = 0; i < c_indices.size(); i++) {
        auto& e = mesh->E.at(c_indices.at(i));
        if (e.Vids.empty()) continue;
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

void SemiGlobalSimplifier::SetSingularityLinks(std::vector<SingularityLink>& SingularityLinks, std::vector<std::vector<size_t>>& SingularityMap, BaseComplexQuad& bc) {
    std::mutex mtx;
    PARALLEL_FOR_BEGIN(mesh->V.size()) {
        auto& v = mesh->V.at(i);
    // for (auto& v: mesh->V) {
        if (!v.isSingularity || v.type == FEATURE || v.isBoundary) continue;
        std::vector<SingularityLink> links = TraceSingularityLinks(v, bc);
        for (auto& l: links) {
            if (mesh->V.at(l.frontId).N_Vids.size() == 4 || mesh->V.at(l.backId).N_Vids.size() == 4) continue;
            {
                std::lock_guard<std::mutex> lock(mtx);
                l.id = SingularityLinks.size();
                l.a = l.linkEids.size();
                l.b = 0;
                SingularityLinks.push_back(l);
                SingularityMap.at(l.frontId).push_back(l.id);
            }
        }
    // }
    } PARALLEL_FOR_END();

    int n = SingularityLinks.size();
    int offset = n;
    std::vector<SingularityLink> newLinks;
    PARALLEL_FOR_BEGIN(n) {
        auto l = SingularityLinks.at(i);
        for (int j = i+1; j < n; j++) {
            auto l2 = SingularityLinks.at(j);
            if (l.frontId == l2.frontId || l.frontId == l2.backId || l.backId == l2.frontId || l.backId == l2.backId) continue;
            std::vector<size_t> a = mu->GetIntersectionParallel(l.linkVids, l2.linkVids);
            if (a.size() != 1) continue;
            int vid = a.at(0);
            std::vector<size_t> b;
            for (auto eid: l.linkEids) {
                auto& e = mesh->E.at(eid);
                b.push_back(eid);
                if (e.Vids.at(1) == vid || e.Vids.at(0) == vid) break;
            }
            if (doesCrossBoundary(b, false)) continue;
            std::vector<size_t> c;
            int it = 0;
            for (auto eid: l2.linkEids) {
                auto& e = mesh->E.at(eid);
                c.push_back(eid);
                it += 1;
                if (e.Vids.at(1) == vid || e.Vids.at(0) == vid) break;
            }
            std::vector<size_t> c1(l2.linkEids.begin()+it, l2.linkEids.end());
            std::reverse(c.begin(), c.end());
            if (!doesCrossBoundary(c, false)) {
                SingularityLink newL;
                newL.a = b.size();
                newL.b = c.size();
                newL.linkEids = b;
                newL.linkEids.insert(newL.linkEids.end(), c.begin(), c.end());
                std::vector<size_t> d;
                int prevId = l.frontId;
                for (auto el: newL.linkEids) {
                    auto& e = mesh->E.at(el);
                    newL.linkVids.push_back(prevId);
                    prevId = e.Vids.at(0) == prevId ? e.Vids.at(1) : e.Vids.at(0);
                }
                newL.linkVids.push_back(prevId);
                newL.frontId = newL.linkVids.front();
                newL.backId = newL.linkVids.back();
                SingularityLink newL2;
                newL2.a = newL.b;
                newL2.b = newL.a;
                newL2.linkEids = newL.linkEids;
                newL2.linkVids = newL.linkVids;
                std::reverse(newL2.linkEids.begin(), newL2.linkEids.end());
                std::reverse(newL2.linkVids.begin(), newL2.linkVids.end());
                newL2.frontId = newL2.linkVids.front();
                newL2.backId = newL2.linkVids.back();
                {
                    std::lock_guard<std::mutex> lock(mtx);
                    newL.id = offset;
                    offset += 1;
                    newLinks.push_back(newL);
                    SingularityMap.at(newL.frontId).push_back(newL.id);
                    
                    newL2.id = offset;
                    offset += 1;
                    newLinks.push_back(newL2);
                    SingularityMap.at(newL2.frontId).push_back(newL2.id);
                }
            }
            if (!doesCrossBoundary(c1, false)) {
                SingularityLink newL;
                newL.a = b.size();
                newL.b = c1.size();
                newL.linkEids = b;
                newL.linkEids.insert(newL.linkEids.end(), c1.begin(), c1.end());
                std::vector<size_t> d;
                int prevId = l.frontId;
                for (auto el: newL.linkEids) {
                    auto& e = mesh->E.at(el);
                    newL.linkVids.push_back(prevId);
                    prevId = e.Vids.at(0) == prevId ? e.Vids.at(1) : e.Vids.at(0);
                }
                newL.linkVids.push_back(prevId);
                newL.frontId = newL.linkVids.front();
                newL.backId = newL.linkVids.back();
                SingularityLink newL2;
                newL2.a = newL.b;
                newL2.b = newL.a;
                newL2.linkEids = newL.linkEids;
                newL2.linkVids = newL.linkVids;
                std::reverse(newL2.linkEids.begin(), newL2.linkEids.end());
                std::reverse(newL2.linkVids.begin(), newL2.linkVids.end());
                newL2.frontId = newL2.linkVids.front();
                newL2.backId = newL2.linkVids.back();
                {
                    std::lock_guard<std::mutex> lock(mtx);
                    newL.id = offset;
                    offset += 1;
                    newLinks.push_back(newL);
                    SingularityMap.at(newL.frontId).push_back(newL.id);
                    
                    newL2.id = offset;
                    offset += 1;
                    newLinks.push_back(newL2);
                    SingularityMap.at(newL2.frontId).push_back(newL2.id);
                }
            }
        }
    } PARALLEL_FOR_END();
    SingularityLinks.insert(SingularityLinks.end(), newLinks.begin(), newLinks.end());
}

void SemiGlobalSimplifier::SelectSingularityGroups(std::vector<SingularityGroup>& Groups, std::vector<SingularityLink>& SingularityLinks, std::vector<std::vector<size_t>>& SingularityMap) {
    std::mutex mtx;
    auto cmp = [](SingularityGroup left, SingularityGroup right) {return left.rank > right.rank;};
    std::priority_queue<SingularityGroup, std::vector<SingularityGroup>, decltype(cmp)> q(cmp);
    
    PARALLEL_FOR_BEGIN(mesh->V.size()) { 
        auto& v = mesh->V.at(i);
        if (!v.isSingularity || v.type == FEATURE || v.isBoundary) continue;
        if (v.N_Vids.size() != 3 && v.N_Vids.size() != 5) continue;
        auto& lids = SingularityMap.at(v.id);
        int valenceToCheck = v.N_Vids.size() == 3 ? 5 : 3;
        int l1_id = -1;
        int l2_id = -1;
        int rank = -1;
        for (int j = 0; j < lids.size(); j++) {
            auto& l = SingularityLinks.at(lids.at(j));
            if (doesCrossBoundary(l.linkVids, true) || mesh->V.at(l.backId).N_Vids.size() != valenceToCheck) continue;
            if (l1_id == -1 || rank > l.a + l.b) {
                l1_id = l.id;
                rank = l.a + l.b;
            }
        }
        if (l1_id == -1) continue;
        bool skipValenceCheck = rank == 1 ? true : false;
        rank = -1;
        for (int j = 0; j < lids.size(); j++) {
            auto& l = SingularityLinks.at(lids.at(j));
            auto& l1 = SingularityLinks.at(l1_id);
            if ((l.frontId == l1.frontId && l.backId == l1.backId) || (l.backId == l1.frontId && l.frontId == l1.backId)) continue;
            if (doesCrossBoundary(l.linkVids, true)) continue;
            if (!skipValenceCheck && mesh->V.at(l.backId).N_Vids.size() != valenceToCheck) continue;
            if (!mu->GetIntersection(std::vector<size_t>(l.linkVids.begin()+1, l.linkVids.end()-1), std::vector<size_t>(l1.linkVids.begin()+1, l1.linkVids.end()-1)).empty()) continue;
            if (l2_id == -1 || rank > l.a + l.b) {
                l2_id = l.id;
                rank = l.a + l.b;
            }
        }
        if (l2_id == -1) continue;
        auto& l1 = SingularityLinks.at(l1_id);
        auto& l2 = SingularityLinks.at(l2_id);
        SingularityGroup sg;
        sg.l1 = (l1.a + l1.b) <= (l2.a + l2.b) ? l1 : l2;
        sg.l2 = sg.l1.id == l1.id ? l2 : l1;
        sg.rank = (sg.l1.a + sg.l1.b + sg.l2.a + sg.l2.b);
        sg.rank += (double)(sg.l1.a + sg.l1.b) / (double)(sg.l1.a + sg.l1.b + sg.l2.a + sg.l2.b);
        {
            std::lock_guard<std::mutex> lock(mtx);
            q.push(sg);
        }
    } PARALLEL_FOR_END();

    std::vector<bool> selected(mesh->V.size(), false);
    while (!q.empty()) {
        auto sg = q.top();
        auto l1 = sg.l1;
        auto l2 = sg.l2;
        if (selected.at(l1.frontId) || selected.at(l1.backId) || selected.at(l2.frontId) || selected.at(l2.backId)) {
            q.pop();
            continue;
        }
        selected.at(l1.frontId) = true;
        selected.at(l1.backId) = true;
        selected.at(l2.frontId) = true;
        selected.at(l2.backId) = true;
        Groups.push_back(sg);
        q.pop();
    }
}

void SemiGlobalSimplifier::ResolveSingularityGroups(std::vector<SingularityGroup>& Groups, BaseComplexQuad& bc) {
    int colorValue = 0;
    int ncolors = 15;
    std::vector<size_t> c_indices;
    std::vector<int> colors;
    
    for (auto& sg: Groups) {
        auto l1 = sg.l1;
        auto l2 = sg.l2;
        // if (l1.a + l1.b > 7) continue;
        if (l1.linkVids.empty() || l2.linkVids.empty() || l1.linkEids.empty() || l2.linkEids.empty()) continue;
        std::cout << l1.frontId << " " << l1.backId << " " << l2.frontId << " " << l2.backId << std::endl;
        size_t toMoveId = l1.frontId;
        size_t sourceId = l1.backId;
        size_t secondaryId = l2.backId;
        if (mesh->V.at(toMoveId).N_Vids.size() == 4 || mesh->V.at(sourceId).N_Vids.size() ==4 || mesh->V.at(secondaryId).N_Vids.size() == 4) continue;
        if (mesh->V.at(toMoveId).N_Vids.empty() || mesh->V.at(sourceId).N_Vids.empty() || mesh->V.at(secondaryId).N_Vids.empty()) continue;
        std::vector<size_t> mainPath = l1.linkVids;
        std::vector<size_t> secondaryPath = l2.linkVids;
        if (sourceId == mainPath.at(0)) {
            std::reverse(mainPath.begin(), mainPath.end());
        }
        if (secondaryId == secondaryPath.at(0)) {
            std::reverse(secondaryPath.begin(), secondaryPath.end());
        }
        
        std::cout << "rank: " << sg.rank << std::endl;
        for (int i = 1; i < mainPath.size(); i++) {
            std::cout << "secondary path: ";
            for (auto el: secondaryPath) {
                std::cout << el << " ";
            }
            std::cout << std::endl;
            if (secondaryPath.empty()) break;
            // size_t newSecondary = MoveSingularities(toMoveId, sourceId, secondaryId, mainPath.at(i), secondaryPath.at(1), std::vector<size_t>(secondaryPath.begin()+1, secondaryPath.end()));
            MoveSingularities(toMoveId, sourceId, secondaryId, mainPath.at(i), secondaryPath.at(1), std::vector<size_t>(secondaryPath.begin()+1, secondaryPath.end()));
            if (i == mainPath.size()-1) continue;
            // secondaryId = newSecondary;
            std::cout << "new Secondary " << secondaryId << " " << mesh->V.at(secondaryId).N_Vids.size() << std::endl;
            if (!secondaryId) break;
            if (mesh->V.at(secondaryId).N_Vids.size() == 4) break;            
            SetSecondaryPath(secondaryId, toMoveId, sourceId, mainPath, secondaryPath, bc);
        }
        std::cout << "at the end of singularity move" << std::endl;
        // std::vector<int> a(l1.linkEids.size(), (colorValue%ncolors));
        // std::vector<int> b(l2.linkEids.size(), (colorValue%ncolors));
        // c_indices.insert(c_indices.end(), l1.linkEids.begin(), l1.linkEids.end());
        // c_indices.insert(c_indices.end(), l2.linkEids.begin(), l2.linkEids.end());
        // colors.insert(colors.end(), a.begin(), a.end());
        // colors.insert(colors.end(), b.begin(), b.end());
        // colorValue += 1;
        // if (l1.a + l1.b > 7 || l2.a + l2.b > 7) break;
    }
    // Smooth();
    /*std::cout << "Writing output file" << std::endl;
    std::ofstream ofs("singularityLinks.vtk");
    ofs << "# vtk DataFile Version 3.0\n"
        << "singularityLinks" << ".vtk\n"
        << "ASCII\n\n"
        << "DATASET UNSTRUCTURED_GRID\n";
    ofs << "POINTS " << mesh->V.size() << " double\n";
    // for (auto& v: source.V) {
    //     if (v.type == FEATURE) c_indices.push_back(v.id);
    // }
    // std::vector<size_t> c_indices = {12, 296};
    // std::cout << c_indices.size() << std::endl;
    for (size_t i = 0; i < mesh->V.size(); i++) {
        ofs << std::fixed << std::setprecision(7) <<  mesh->V.at(i).x << " " <<  mesh->V.at(i).y << " " <<  mesh->V.at(i).z << "\n";
    }
    ofs << "CELLS " << c_indices.size() << " " << 3 * c_indices.size() << std::endl;
    for (size_t i = 0; i < c_indices.size(); i++) {
        auto& e = mesh->E.at(c_indices.at(i));
        if (e.Vids.empty()) continue;
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

void SemiGlobalSimplifier::SetSecondaryPath(size_t& secondaryId, size_t& toMoveId, size_t& sourceId, std::vector<size_t>& mainPath, std::vector<size_t>& secondaryPath, BaseComplexQuad& bc) {
    std::vector<SingularityLink> secondaryLinks = TraceSingularityLinks(mesh->V.at(secondaryId), bc);
    // std::vector<SingularityLink> moveLinks = TraceSingularityLinks(mesh->V.at(toMoveId), bc);
    std::vector<SingularityLink> cmpLinks = TraceSingularityLinks(mesh->V.at(toMoveId), bc);
    // std::vector<SingularityLink> cmpLinks = TraceSingularityLinks(mesh->V.at(sourceId), bc);
    // cmpLinks.insert(cmpLinks.end(), moveLinks.begin(), moveLinks.end());
    auto cmp = [](SingularityLink left, SingularityLink right) {return left.a + left.b < right.a + right.b;};
    std::priority_queue<SingularityLink, std::vector<SingularityLink>, decltype(cmp)> q(cmp);
    
    for (auto l: secondaryLinks) {
        if (l.backId == toMoveId) {
            q.push(l);
            continue;
        }
        for (auto l2: cmpLinks) {
            if (l.frontId == l2.frontId || l.frontId == l2.backId || l.backId == l2.frontId || l.backId == l2.backId) continue;
            std::vector<size_t> a = mu->GetIntersectionParallel(l.linkVids, l2.linkVids);
            if (a.size() != 1) continue;
            int vid = a.at(0);
            std::vector<size_t> b;
            for (auto eid: l.linkEids) {
                auto& e = mesh->E.at(eid);
                b.push_back(eid);
                if (e.Vids.at(1) == vid || e.Vids.at(0) == vid) break;
            }
            if (doesCrossBoundary(b, false)) continue;
            std::vector<size_t> c;
            int it = 0;
            for (auto eid: l2.linkEids) {
                auto& e = mesh->E.at(eid);
                c.push_back(eid);
                it += 1;
                if (e.Vids.at(1) == vid || e.Vids.at(0) == vid) break;
            }
            std::reverse(c.begin(), c.end());
            if (!doesCrossBoundary(c, false)) {
                SingularityLink newL;
                newL.a = b.size();
                newL.b = c.size();
                newL.linkEids = b;
                newL.linkEids.insert(newL.linkEids.end(), c.begin(), c.end());
                std::vector<size_t> d;
                int prevId = l.frontId;
                for (auto el: newL.linkEids) {
                    auto& e = mesh->E.at(el);
                    newL.linkVids.push_back(prevId);
                    prevId = e.Vids.at(0) == prevId ? e.Vids.at(1) : e.Vids.at(0);
                }
                newL.linkVids.push_back(prevId);
                newL.frontId = newL.linkVids.front();
                newL.backId = newL.linkVids.back();
                if (!mu->GetIntersection(std::vector<size_t>(mainPath.begin()+1, mainPath.end()-1), std::vector<size_t>(newL.linkVids.begin()+1, newL.linkVids.end()-1)).empty()) continue;
                q.push(newL);
            }
        }   
    }
    secondaryPath.clear();
    if (!q.empty()) {
        auto& l = q.top();
        secondaryPath = l.linkVids;
        std::reverse(secondaryPath.begin(), secondaryPath.end());
    }
}

void SemiGlobalSimplifier::MoveSingularities(size_t& toMoveId, size_t& sourceId, size_t& secondaryId, size_t& sourceDir, size_t& secondaryDir, std::vector<size_t>& secondaryPath) {
    auto& toMove = mesh->V.at(toMoveId);
    std::cout << "toMove: " << toMoveId << " " << toMove.N_Vids.size() << " source: " << sourceId << " " << mesh->V.at(sourceId).N_Vids.size() << " secondary: " << secondaryId << " " << mesh->V.at(secondaryId).N_Vids.size() << std::endl; 
    std::cout << "toMove: " << toMoveId << " " << toMove.N_Vids.size() << " source dir: " << sourceDir << " " << mesh->V.at(sourceDir).N_Vids.size() << " secondary dir: " << secondaryDir << " " << mesh->V.at(secondaryDir).N_Vids.size() << std::endl; 
    size_t threeId, fiveId;
    if (sourceId == sourceDir) {
        threeId = toMove.N_Vids.size() == 3 ? toMoveId : sourceId;
        fiveId = toMove.N_Vids.size() == 5 ? toMoveId : sourceId;
    } else {
        std::vector<size_t> threeFiveIds = GetThreeFivePairIds(toMove.id, sourceDir, secondaryDir);
        threeId = threeFiveIds.at(0);
        fiveId = threeFiveIds.at(1);
    }
    toMoveId = sourceDir;
    std::cout << "Three Five Ids: " << mesh->V.at(threeId).N_Vids.size() << " " << mesh->V.at(fiveId).N_Vids.size() << std::endl;
    if (threeId == secondaryId || fiveId == secondaryId) {
        if (mesh->V.at(threeId).N_Vids.size() == 4 && mesh->V.at(fiveId).N_Vids.size() == 5) {
            // return fiveId;
            secondaryId = fiveId;
            return;
        } else if (mesh->V.at(threeId).N_Vids.size() == 3 && mesh->V.at(fiveId).N_Vids.size() == 4) {
            // return threeId;
            secondaryId = threeId;
            return;
        } else if (mesh->V.at(threeId).N_Vids.size() == 3 && mesh->V.at(fiveId).N_Vids.size() == 6) {
            std::shared_ptr<ThreeFivePair> tfp = std::make_shared<ThreeFivePair>(*mesh, *mu, *smoother, threeId, fiveId);
            tfp->SplitSixSingularity();
            // return tfp->GetPairIds().at(1);
            secondaryId = tfp->GetPairIds().at(1);
            return;
        }
    }
    std::shared_ptr<ThreeFivePair> tfp = std::make_shared<ThreeFivePair>(*mesh, *mu, *smoother, threeId, fiveId);
    std::vector<size_t> pids;
    for (int i = 0; i < secondaryPath.size(); i++) {
        std::cout << "Moving to: " << secondaryPath.at(i) << " with valence: " << mesh->V.at(secondaryPath.at(i)).N_Vids.size() << std::endl;
        tfp->Move(secondaryPath.at(i));
        threeId = tfp->GetPairIds().at(0);
        fiveId = tfp->GetPairIds().at(1);
        if (threeId == secondaryId || fiveId == secondaryId || mesh->V.at(secondaryId).N_Vids.empty()) {
            if (mesh->V.at(threeId).N_Vids.size() == 4 && mesh->V.at(fiveId).N_Vids.size() == 5) {
                // return fiveId;
                secondaryId = fiveId;
                return;
            } else if (mesh->V.at(threeId).N_Vids.size() == 3 && mesh->V.at(fiveId).N_Vids.size() == 4) {
                // return threeId;
                secondaryId = threeId;
                return;
            } else if (mesh->V.at(threeId).N_Vids.size() == 3 && mesh->V.at(fiveId).N_Vids.size() == 6) {
                std::shared_ptr<ThreeFivePair> tfp = std::make_shared<ThreeFivePair>(*mesh, *mu, *smoother, threeId, fiveId);
                tfp->SplitSixSingularity();
                // return tfp->GetPairIds().at(1);
                secondaryId = tfp->GetPairIds().at(1);
                return;
            }
        }
        std::cout << "Three Five Ids: " << mesh->V.at(tfp->GetPairIds().at(0)).N_Vids.size() << " " << mesh->V.at(tfp->GetPairIds().at(1)).N_Vids.size() << std::endl;
    }
    size_t newSecondary = tfp->GetResolvedSingularity();
    if (newSecondary != -1) {
        secondaryId = newSecondary;
    }
    // std::cout << "newSecondary: " << newSecondary << std::endl;
    // toMoveId = sourceDir;
    // return newSecondary;
}

std::vector<size_t> SemiGlobalSimplifier::GetThreeFivePairIds(size_t vid, size_t mainId, size_t secondaryId) {
    std::vector<size_t> res;
    auto& v = mesh->V.at(vid);
    auto& main = mesh->V.at(mainId);
    auto& secondary = mesh->V.at(secondaryId);
    // std::cout << "Inside GetThreeFivePairIds, v: " << v.id << " main: " << main.id << " secondary: " << secondary.id << std::endl; 
    if (v.N_Vids.size() == 3) {
        for (auto fid: v.N_Fids) {
            auto& f = mesh->F.at(fid);
            // for (auto fvid: f.Vids) {
            //     if (mesh->V.at(fvid).isBoundary || mesh->V.at(fvid).type == FEATURE) return res;
            // }
            if (std::find(f.Vids.begin(), f.Vids.end(), main.id) != f.Vids.end()
            && std::find(f.Vids.begin(), f.Vids.end(), secondary.id) != f.Vids.end()) {
                int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
                if (secondary.N_Vids.size() == 3) {
                    size_t id = mu->GetDifference(secondary.N_Vids, std::vector<size_t>{vid, f.Vids.at((idx+2)%f.Vids.size())}).at(0);
                    res.insert(res.begin(), id);
                } else {
                    res.insert(res.begin(), secondary.id);
                }
                std::shared_ptr<SimplificationOperation> dc = std::make_shared<DiagonalCollapse>(*mesh, *mu, *smoother, f.id, (idx+2)%f.Vids.size(), idx);
                dc->PerformOperation();
                res.push_back(f.Vids.at((idx+2)%f.Vids.size()));
                // for (auto nvid: mesh->V.at(res.at(0)).N_Vids) {
                //     if (nvid == main.id) continue;
                //     if (mesh->V.at(nvid).N_Vids.size() == 3) res.insert(res.begin(), nvid);
                // }
                break;
            }
        }
    } else if (v.N_Vids.size() == 5) {
        std::vector<size_t> verticesToSplit;
        std::vector<size_t> verticesToChange;
        size_t startE;
        // std::cout << "v nvids: ";
        // for (auto id: v.N_Vids) {
        //     std::cout << id << " ";
        // }
        // std::cout << std::endl;
        // std::cout << "v neids: " << std::endl;
        for (auto eid: v.N_Eids) {
            auto& e = mesh->E.at(eid);
            // std::cout << e.Vids.at(0) << " " << e.Vids.at(1) << std::endl;
            if (e.Vids.at(0) == main.id || e.Vids.at(1) == main.id) {
                startE = e.id;
            }
        }
        
        int d = (v.N_Vids.size() / 2) + 1;
        // std::cout << "d: " << d << std::endl;
        bool addSplit = false;
        for (int j = 0; j < d; j++) {
            auto& e = mesh->E.at(startE);
            size_t ev = e.Vids.at(0) == v.id ? e.Vids.at(1) : e.Vids.at(0);
            for (auto fid: e.N_Fids) {
                auto& f = mesh->F.at(fid);
                int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
                // std::cout << "idx: " << idx << std::endl;
                if (f.Vids.at((idx+1)%f.Vids.size()) == ev) {
                    if (ev == secondary.id && j <= d-1) {
                        addSplit = true;
                    }
                    if (j == 0 || (j == d-1 && addSplit)) {
                        verticesToSplit.push_back(ev);
                    } else {
                        verticesToChange.push_back(ev);
                    }
                    startE = mu->GetDifference(mu->GetIntersection(f.Eids, v.N_Eids), std::vector<size_t>{e.id}).at(0);
                }
            }
        }
        auto& e = mesh->E.at(startE);
        size_t ev = e.Vids.at(0) == v.id ? e.Vids.at(1) : e.Vids.at(0);
        if (verticesToSplit.size() < 2) verticesToSplit.push_back(ev);
        // std::cout << "Before quad split" << std::endl;
        // std::cout << verticesToSplit.size() << " " << verticesToChange.size() << std::endl;
        if (mesh->V.at(verticesToSplit.at(0)).isBoundary || mesh->V.at(verticesToSplit.at(0)).type == FEATURE
        || mesh->V.at(verticesToSplit.at(1)).isBoundary || mesh->V.at(verticesToSplit.at(1)).type == FEATURE) return res;
        
        // std::cout << "Mesh vertices: " << mesh->V.size() << std::endl;
        
        std::shared_ptr<SimplificationOperation> qs = std::make_shared<QuadSplit>(*mesh, *mu, *smoother, v.id, verticesToSplit, verticesToChange);
        qs->PerformOperation();

        auto& updatedV = mesh->V.at(vid);
        // std::cout << "Mesh vertices: " << mesh->V.size() << std::endl;

        // std::cout << "After quad split" << std::endl;
        // std::cout << "v: " << updatedV.N_Vids.size() << std::endl;
        // std::cout << "Last v: " << mesh->V.at(mesh->V.size()-1).N_Vids.size() << std::endl;
        if (updatedV.N_Vids.size() == 3) {
            res.insert(res.begin(), updatedV.id);
        } else if (mesh->V.at(mesh->V.size()-1).N_Vids.size() == 3) {
            res.insert(res.begin(), mesh->V.at(mesh->V.size()-1).id);
        }
        // std::cout << "verticesToSplit at 1: " << mesh->V.at(verticesToSplit.at(1)).N_Vids.size() << std::endl;
        res.push_back(verticesToSplit.at(1));
        // for (auto nvid: mesh->V.at(res.at(0)).N_Vids) {
        //     if (nvid == main.id) continue;
        //     if (mesh->V.at(nvid).N_Vids.size() == 5) {
        //         res.push_back(nvid);
        //         break;
        //     }
        // }
    }

    return res;
}

std::vector<SingularityLink> SemiGlobalSimplifier::TraceSingularityLinks(Vertex& v, BaseComplexQuad& bc) {
    std::vector<bool> is_mesh_edge_visited(mesh->E.size(), false);
    std::vector<SingularityLink> links;
    for (auto edgeid : v.N_Eids) {
        const Edge& edge = mesh->E.at(edgeid);
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
	auto& v = mesh->V.at(vid);
	size_t res;
	for (auto fid : v.N_Fids) {
		auto& f = mesh->F.at(fid);
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
            auto& v = mesh->V.at(id);
            if (v.type == FEATURE || v.isBoundary) {
                return true;
            }
        } else {
            auto& e = mesh->E.at(id);
            if (mesh->V.at(e.Vids.at(0)).type == FEATURE || mesh->V.at(e.Vids.at(0)).isBoundary || mesh->V.at(e.Vids.at(1)).type == FEATURE || mesh->V.at(e.Vids.at(1)).isBoundary) {
                return true;
            }
        }
    }
    return false;
}


size_t SemiGlobalSimplifier::GetDiagonalV(size_t vid, size_t fid) {
	auto& f = mesh->F.at(fid);
    int index = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), vid));
	return f.Vids.at((index+2)%f.Vids.size());
}

void SemiGlobalSimplifier::Smooth() {
    std::cout << "Smoothing mesh" << std::endl;
    // std::vector<size_t> smoothv;
    // for (auto& v: mesh->V) {
    //     if (v.N_Vids.empty() || v.N_Eids.empty() || v.N_Fids.empty()) continue;
    //     smoothv.push_back(v.id);
    // }
    // smoother->Smooth(smoothv);
    // smoother->SetMesh(*mesh);
    smoother->Smooth(std::vector<size_t>{});
}

bool SemiGlobalSimplifier::FixValences() {
    bool res = false;
    for (auto& v: mesh->V) {
        VertexSplit s(*mesh, *mu, *smoother, v.id);
        s.FixDoublet(v.id);
    }

    ResolveSingularityPairs();
    for (auto& v: mesh->V) {
        if ((v.type == FEATURE || v.isBoundary) && !v.N_Fids.empty()) mesh->SetIdealValence(v.id);
    }
    res = FixBoundary();
    res = ResolveHighValences();
    return res;
}

// std::cout << "Writing output file" << std::endl;
// std::string outputf = "links.vtk";
// std::ofstream ofs(outputf.c_str());
// ofs << "# vtk DataFile Version 3.0\n"
//     << outputf.c_str() << "\n"
//     << "ASCII\n\n"
//     << "DATASET UNSTRUCTURED_GRID\n";
// ofs << "POINTS " << mesh->V.size() << " double\n";
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
// for (size_t i = 0; i < mesh->V.size(); i++) {
//     ofs << std::fixed << std::setprecision(7) <<  mesh->V.at(i).x << " " <<  mesh->V.at(i).y << " " <<  mesh->V.at(i).z << "\n";
// }
// ofs << "CELLS " << c_indices.size() << " " << 2 * c_indices.size() << std::endl;
// for (size_t i = 0; i < c_indices.size(); i++) {
//     ofs << "1 " << c_indices.at(i) << std::endl;
// }
// ofs << "CELL_TYPES " << c_indices.size() << "\n";
// for (size_t i = 0; i < c_indices.size(); i++) {
//     ofs << "1" << std::endl;
// }