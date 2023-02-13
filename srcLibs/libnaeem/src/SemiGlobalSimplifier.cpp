#include <algorithm>
#include <map>
#include <time.h>
#include <queue>
#include <deque>
#include <utility>
#include <memory>
#include <cmath>
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
            // for (int j = 0; j < v.N_Fids.size(); j++) {
            //     int count = 0;
            //     auto& f = mesh->F.at(v.N_Fids.at(j));
            //     for (auto vid: f.Vids) {
            //         if (mesh->V.at(vid).type == FEATURE) count += 1;
            //     }
            //     if (count > 2) performVertexSplit = false;
            // }
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
        if (f.N_Fids.empty() || f.Vids.empty()) continue;
        for (int i = 0; i < f.Vids.size(); i++) {
            // if (mesh->V.at(f.Vids.at(i)).N_Fids.size() == 3 && mesh->V.at(f.Vids.at((i+2)%f.Vids.size())).N_Fids.size() == 3) {
            if (mesh->V.at(f.Vids.at((i+1)%f.Vids.size())).N_Fids.size() >= 4 || mesh->V.at(f.Vids.at((i+3)%f.Vids.size())).N_Fids.size() >= 4) {
                    std::shared_ptr<SimplificationOperation> dc1 = std::make_shared<DiagonalCollapse>(*mesh, *mu, *smoother, f.id, i, (i+2)%f.Vids.size());
                    // dc1->SetRanking();
                    Ops.push_back(dc1);
                // break;
            }
        }
        
        // std::unique_ptr<SimplificationOperation> dc2 = std::make_unique<DiagonalCollapse>(*mesh, *mu, *smoother, f.id, 1, 3);
        // dc2->SetRanking();
        // Ops.push_back(std::move(dc2));
    }
    std::cout << Ops.size() << std::endl;
    int i = 0;
    for (auto op: Ops) {
        op->PerformOperation();
        i += 1;
        // if (i >= iters) break;        
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
        auto op = Op_Q.pop();
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
        auto op = Op_Q.pop();
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
            if (v_front.isBoundary || v_back.isBoundary || v_front.type == FEATURE || v_back.type == FEATURE) continue;
            if (v_front.id == v_back.id) continue;
            if (v_front.N_Fids.size() != 3 || v_back.N_Fids.size() != 3) continue;
            if (mesh->V.at(GetDiagonalV(v_front.id, GetFaceID(v_front.id, linkV.at(1)))).N_Fids.size() <= 4 || mesh->V.at(GetDiagonalV(v_back.id, GetFaceID(v_back.id, linkV.at(linkV.size()-2)))).N_Fids.size() <= 4) continue;

            std::shared_ptr<SimplificationOperation> s = std::make_shared<SeparatrixCollapse>(*mesh, *mu, *smoother, linkV, linkE);
            s->SetRanking();
            Op_Q.insert(s->ranking, s->GetCenterId(), s);
        }
    }
    std::cout << Op_Q.size() << " separatrix operations" << std::endl;
    int i = 0;
    while (!Op_Q.empty()) {
        auto op = Op_Q.pop();
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
                if (v_front.isBoundary || v_back.isBoundary || v_front.type == FEATURE || v_back.type == FEATURE) continue;
                if (v_front.id == v_back.id) continue;
                if (v_front.N_Fids.size() != 3 || v_back.N_Fids.size() != 3) continue;
                if (mesh->V.at(GetDiagonalV(v_front.id, GetFaceID(v_front.id, linkV.at(1)))).N_Fids.size() <= 4 || mesh->V.at(GetDiagonalV(v_back.id, GetFaceID(v_back.id, linkV.at(linkV.size()-2)))).N_Fids.size() <= 4) continue;

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
            // if (mesh->V.at(GetDiagonalV(v_front.id, GetFaceID(v_front.id, linkV.at(1)))).N_Fids.size() <= 4 || mesh->V.at(GetDiagonalV(v_back.id, GetFaceID(v_back.id, linkV.at(linkV.size()-2)))).N_Fids.size() <= 4) continue;
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
        auto op = Op_Q.pop();
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
        auto op = Op_Q.pop();
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
        s->PerformOperation();
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
        // tfp->MoveLowerLeft();
        // tfp->MoveLowerLeft();
        // tfp->MoveLowerLeft();
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
        std::vector<size_t> edgesToCheck;
        if (v.N_Vids.size() == 6 && v.type != FEATURE && !v.isBoundary) {
            // std::cout << "vid: " << v.id << " " << v.N_Vids.size() << " " << v.N_Eids.size() << " " << v.N_Fids.size() << std::endl;
            bool split = false;
            for (auto fid: v.N_Fids) {
                auto& f = mesh->F.at(fid);
                int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
                if (mesh->V.at(f.Vids.at((idx+2)%f.Vids.size())).type != FEATURE && !mesh->V.at(f.Vids.at((idx+2)%f.Vids.size())).isBoundary && mesh->V.at(f.Vids.at((idx+2)%f.Vids.size())).N_Vids.size() == 3) {
                    split= true;
                    std::vector<size_t> edgesToSplit;
                    for (auto eid: mu->GetIntersection(f.Eids, v.N_Eids)) {
                        auto& e = mesh->E.at(eid);
                        mu->AddContents(edgesToSplit, mu->GetDifference(mu->GetIntersection(mesh->F.at(mu->GetDifference(e.N_Fids, std::vector<size_t>{fid})[0]).Eids, v.N_Eids), std::vector<size_t>{eid}));
                    }
                    // std::cout << "Fixing six singularity with 3 singularity" << std::endl;
                    std::shared_ptr<SimplificationOperation> s = std::make_shared<VertexSplit>(*mesh, *mu, *smoother, v.id, edgesToSplit);
                    s->PerformOperation();
                    res = true;
                    break;
                }
            }
            if (split) continue;
            /*size_t startE = v.N_Eids.at(0);
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
            }*/
        }
        if (v.N_Vids.size() > 5 && v.type != FEATURE && !v.isBoundary) {
            // std::cout << "vid: " << v.id << " " << v.N_Vids.size() << " " << v.N_Eids.size() << " " << v.N_Fids.size() << std::endl;
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
            // std::cout << "Fixing six singularity with quad split" << std::endl;
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
        if (l.backId == vertexToSkip) continue;
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
            l1.volt = l1.rots > 1 ? 1 : 0;
            l2.volt = l2.rots > 1 ? 1 : 0;
            s.l1l2Rots = PrototypeGetRotations(l1.frontId, l1.linkEids.front(), l2.linkEids.front());
            s.l2l1Rots = PrototypeGetRotations(l1.frontId, l2.linkEids.front(), l1.linkEids.front());
            if (mesh->V.at(l1.frontId).N_Vids.size() == 3) {
                s.l1l2Volt = s.l1l2Rots > 1 ? 1 : 0;
                s.l2l1Volt = s.l2l1Rots > 1 ? 1 : 0;
            } else if (mesh->V.at(l1.frontId).N_Vids.size() == 5) {
                s.l1l2Volt = s.l1l2Rots <= 2 ? 0 : 1;
                s.l2l1Volt = s.l2l1Rots <= 2 ? 0 : 1;
                // s.l1l2Volt = (s.l1l2Rots == 2 && s.l2l1Rots == 3) || (s.l1l2Rots > 2) ? 1 : 0;
                // s.l2l1Volt = (s.l1l2Rots == 3 && s.l2l1Rots == 2) || (s.l2l1Rots > 1) ? 1 : 0;
            }
            std::string link1Volt = l1.volt == 1 ? "High" : "Low";
            std::string link2Volt = l2.volt == 1 ? "High" : "Low";
            std::string l1l2Volt = s.l1l2Volt == 1 ? "High" : "Low";
            std::string l2l1Volt = s.l2l1Volt == 1 ? "High" : "Low";
            
            std::cout << "l1: " << l1.a+l1.b << " l1 volt: " << link1Volt << std::endl;
            std::cout << "l2: " << l2.a+l2.b << " l2 volt: " << link2Volt << std::endl;
            std::cout << "link volt from l1 to l2: " << l1l2Volt << std::endl;
            std::cout << "link volt from l2 to l1: " << l2l1Volt << std::endl;
            std::string valid = PrototypeIsLinkValid(s) == true ? "Yes" : "No";
            std::cout << "Is link Valid? " << valid << std::endl;
            if (PrototypeIsLinkValid(s)) std::cout << "Element change: " << PrototypeGetElementPrediction(s) << std::endl;

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

void SemiGlobalSimplifier::PrototypeF(int idxOffset) {
    std::cout << "Prototype F" << std::endl;
    BaseComplexQuad bc(*mesh);
    int n = mesh->V.size();
    for (int vid = 0; vid < n; vid++) {
        auto& v = mesh->V.at(vid);
        if (v.type == FEATURE || v.isBoundary || (v.N_Vids.size() != 3 && v.N_Vids.size() != 5)) continue;
        int valenceToCheck = v.N_Vids.size() == 5 ? 3 : 5;
        bool breakLoop = false;
        for (auto fid: v.N_Fids) {
            auto& f = mesh->F.at(fid);
            int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
            auto& fv = mesh->V.at(f.Vids.at((idx+idxOffset)%f.Vids.size()));
            if (fv.type != FEATURE && !fv.isBoundary && fv.N_Vids.size() == valenceToCheck) {
                std::vector<size_t> threeFiveIds = valenceToCheck == 5 ? std::vector<size_t>{v.id, fv.id} : std::vector<size_t>{fv.id, v.id};
                // std::cout << mesh->V.at(threeFiveIds.at(0)).id << " " << mesh->V.at(threeFiveIds.at(1)).id << std::endl;
                // std::cout << mesh->V.at(threeFiveIds.at(0)).N_Vids.size() << " " << mesh->V.at(threeFiveIds.at(1)).N_Vids.size() << std::endl;
                SingularityLink l1 = PrototypeGetLink(v.id, bc, fv.id, std::vector<size_t>{}, false, false);
                SingularityLink l2 = PrototypeGetLink(fv.id, bc, v.id, std::vector<size_t>{}, false, false);
                if (l1.linkVids.empty() || l2.linkVids.empty()) continue;
                std::vector<size_t> toMove = l1.a+l1.b <= l2.a+l2.b ? l1.linkVids : l2.linkVids;
                // size_t step = toMove.at(1);
                // if (mu->GetDifference(mesh->V.at(threeFiveIds.at(0)).N_Vids, mesh->V.at(threeFiveIds.at(1)).N_Vids).at(0) != step) {
                //     continue;
                // }
                // std::cout << "Got a down direction" << std::endl;
                // std::cout << "toMove size: " << toMove.size() << std::endl;
                if (idxOffset == 2) {
                    std::shared_ptr<DiagonalThreeFivePair> tfp = std::make_shared<DiagonalThreeFivePair>(*mesh, *mu, *smoother, threeFiveIds.at(0), threeFiveIds.at(1));
                    for (int i = 1; i < toMove.size(); i++) {
                        tfp->Move(toMove.at(i), delta);
                        if (!CheckMeshValidity()) return;
                    }
                } else {
                    std::shared_ptr<ThreeFivePair> tfp = std::make_shared<ThreeFivePair>(*mesh, *mu, *smoother, threeFiveIds.at(0), threeFiveIds.at(1));
                    for (int i = 1; i < toMove.size(); i++) {
                        tfp->Move(toMove.at(i), delta);
                    }
                }
                break;
            }
        }
        // if (breakLoop) break;
    }
    while (FixValences());
}

void SemiGlobalSimplifier::PrototypeG(int vid, BaseComplexQuad& bc) {
    auto& v = mesh->V.at(vid);
    std::vector<SingularityLink> links;
    TraceSingularityLinks(v.id, links);
    PrototypeSaveMesh(links, "test");
    std::cout << "vertex: " << v.id << " " << v.N_Vids.size() << " " << links.size() << std::endl;
}

bool SemiGlobalSimplifier::PrototypeH(int idxOffset) {
    bool res = false;
    int it = 0;
    std::vector<std::shared_ptr<DiagonalThreeFivePair>> DiagonalPairs;
    std::vector<std::shared_ptr<ThreeFivePair>> DirectPairs;
    std::priority_queue<SingularityLink, std::vector<SingularityLink>, LinkComparator> q;

    // std::cout << "delta: " << delta << std::endl;
    // for (auto& v: mesh->V) {
    //     if (v.isBoundary || v.type == FEATURE || v.N_Vids.empty() || (v.N_Vids.size() != 3 && v.N_Vids.size() != 5)) continue;
    //     SingularityLink l = PrototypeGetLink(v.id);
    //     PrototypeMoveSingularity(l);
    //     it += 1;
    //     if (it >= iters) break;
    // }
    // return true;
    for (auto& v: mesh->V) {
        if (v.isBoundary || v.type == FEATURE || (v.N_Vids.size() != 3 && v.N_Vids.size() != 5)) continue;
        q.push(PrototypeGetLink(v.id));
        /*for (auto fid: v.N_Fids) {
            int valenceToCheck = v.N_Vids.size() == 5 ? 3 : 5;
            for (int n = 1; n <= 3; n++) {
                auto& fv = mesh->V.at(GetFaceV(v.id, fid, n));
                if (mu->IsSharpFeature(fv.id) || fv.N_Vids.size() != valenceToCheck) continue;
                std::vector<size_t> threeFiveIds = valenceToCheck == 5 ? std::vector<size_t>{v.id, fv.id} : std::vector<size_t>{fv.id, v.id};
                if (n%2==0) {
                    std::shared_ptr<DiagonalThreeFivePair> tfp = std::make_shared<DiagonalThreeFivePair>(*mesh, *mu, *smoother, threeFiveIds.at(0), threeFiveIds.at(1));
                    DiagonalPairs.push_back(tfp);
                } else {
                    std::shared_ptr<ThreeFivePair> tfp = std::make_shared<ThreeFivePair>(*mesh, *mu, *smoother, threeFiveIds.at(0), threeFiveIds.at(1));
                    DirectPairs.push_back(tfp);
                }
            }
            if (fv.type != FEATURE && !fv.isBoundary && fv.N_Vids.size() == valenceToCheck) {
            // if (fv.N_Vids.size() == valenceToCheck) {
                it += 1;
                // if (it < iters) continue;
                std::vector<size_t> threeFiveIds = valenceToCheck == 5 ? std::vector<size_t>{v.id, fv.id} : std::vector<size_t>{fv.id, v.id};
                if (idxOffset == 2) {
                    res = true;
                    std::shared_ptr<DiagonalThreeFivePair> tfp = std::make_shared<DiagonalThreeFivePair>(*mesh, *mu, *smoother, threeFiveIds.at(0), threeFiveIds.at(1));
                    std::vector<SingularityLink> links;
                    TraceSingularityLinks(threeFiveIds.at(0), links);
                    TraceSingularityLinks(threeFiveIds.at(1), links);
                    SelectDiagonalPairLink(threeFiveIds, links);
                    // if (it == iters) {
                    //     PrototypeSaveMesh(links, "test");
                    // }
                    std::vector<size_t> path(links.at(0).linkVids.begin()+1, links.at(0).linkVids.end());
                    for (int i = 0; i < path.size(); i++) {
                        tfp->Move(path.at(i));
                        // if (it == iters && i == 0) break;
                        // if (it == it2 && i == iters) break;
                    }
                    // for (auto& l: links) std::cout << l.rank << " ";
                    // std::cout << std::endl;
                } else {
                    res = true;
                    std::shared_ptr<ThreeFivePair> tfp = std::make_shared<ThreeFivePair>(*mesh, *mu, *smoother, threeFiveIds.at(0), threeFiveIds.at(1));
                    std::vector<SingularityLink> links;
                    TraceSingularityLinks(threeFiveIds.at(0), links);
                    TraceSingularityLinks(threeFiveIds.at(1), links);
                    // std::cout << "direct pair links: " << links.size() << std::endl;
                    SelectDirectPairLink(threeFiveIds, links);
                    // if (it == iters) {
                    //     PrototypeSaveMesh(links, "test");
                    // }
                    std::vector<size_t> path(links.at(0).linkVids.begin()+1, links.at(0).linkVids.end());
                    // std::cout << "linkVids size: " << links.at(0).linkVids.size() << " path size: " << path.size() << std::endl;
                    // std::cout << "three: " << threeFiveIds.at(0) << " five: " << threeFiveIds.at(1) << std::endl;
                    // for (auto id: links.at(0).linkVids) std::cout << id << " ";
                    // std::cout << std::endl;
                    // for (auto id: path) std::cout << id << " ";
                    // std::cout << std::endl;
                    for (int i = 0; i < path.size(); i++) {
                        tfp->Move(path.at(i));
                        // std::cout << "move iteration: " << i << std::endl;
                    }
                    // for (auto& l: links) std::cout << l.rank << " ";
                    // std::cout << std::endl;
                    // PrototypeSaveMesh(links, "test");
                }
                break;
            }
        }*/
    }
    std::cout << "Singularity Links: " << q.size() << std::endl;
    // std::vector<SingularityLink> links;
    while (!q.empty()) {
        // links.push_back(q.top());
        auto l = q.top();
        q.pop();
        it += 1;
        if (PrototypeIsLinkValid((SingularityLink&) l)) {
            std::cout << "working on link of rank: " << l.rank << std::endl;
            PrototypeMoveSingularity((SingularityLink&) l);
            while (FixValences());
            if (it == iters) {
                break;
            }
            // std::cout << "After singularity movement" << std::endl;
        }
    }
    // PrototypeSaveMesh(links, "test");

    // std::cout << "Diagonal Pairs: " << DiagonalPairs.size() << std::endl;
    // std::cout << "Direct Pairs: " << DirectPairs.size() << std::endl;
    return res;
}

bool SemiGlobalSimplifier::PrototypeI() {
    bool res = false;
    res = PrototypeH();
    res = PrototypeH(1);
    // res = ResolveHighValences();
    // while (FixValences());
    return res;
}

void SemiGlobalSimplifier::PrototypeJ() {
    int it = 0;
    for (auto& v: mesh->V) {
        if (v.type != FEATURE && !v.isBoundary && (v.N_Vids.size() == 3 || v.N_Vids.size() == 5)) {
            std::vector<SingularityLink> links;
            TraceSingularityLinks(v.id, links, true);
            it += 1;
            if (it == iters) {
                SelectLinks(links);
                PrototypeSaveMesh(links, "test");
                PrototypeMoveSingularity(links.at(0));
                std::cout << "rank: " << links.at(0).rank << " a: " << links.at(0).a << " b: " << links.at(0).b << std::endl;
                std::cout << "linkVids: " << links.at(0).linkVids.size() << std::endl;
                break;
            }
        }
    }
}

void SemiGlobalSimplifier::PrototypeMoveSingularity(SingularityLink& l) {
    int i = 0;
    int j = l.linkVids.size()-1;
    bool toggleDirection = true;
    int iter = 0;
    // PrototypeSaveMesh(std::vector<SingularityLink> {l}, "test");
    /*while (i < j) {
        std::cout << "i: " << i << " j: " << j << std::endl; 
        // int n = j - i;
        std::vector<size_t> tfp, vToAvoid;
        bool diagPair;
        SingularityLink pairLink;
        if (i+1 == j) {
            tfp = {l.linkVids.at(i), l.linkVids.at(i+1)};
            vToAvoid = tfp;
            if (mesh->V.at(tfp.at(0)).N_Vids.size() != 3) std::reverse(tfp.begin(), tfp.end());
            diagPair = Contains(mesh->V.at(tfp.at(0)).N_Vids, tfp.at(1)) ? false : true;
            pairLink = PrototypePairLink(tfp, vToAvoid, diagPair);
            // PrototypeSaveMesh(std::vector<SingularityLink> {pairLink}, "test2");
            PrototypeMovePair(tfp, pairLink, diagPair); 
            i += 1;
            std::cout << "delta inside move singularity loop: " << delta << std::endl;
            continue;
        }
        int toMove = i;
        int dest = i+1;
        vToAvoid = {l.linkVids.at(i+1), l.linkVids.at(j)};
        // if (!toggleDirection) {
        //     toMove = j;
        //     dest = j-1;
        //     vToAvoid = {l.linkVids.at(j-1), l.linkVids.at(i)};
        // }
        diagPair = PrototypePairIds(tfp, l.linkVids, toMove, dest);
        mu->AddContents(vToAvoid, tfp);
        pairLink = PrototypePairLink(tfp, vToAvoid, diagPair);
        PrototypeMovePair(tfp, pairLink, diagPair);
        delta += mesh->V.at(dest).N_Vids.size() == 3 ? -2 : 2;
        // if (toggleDirection) {
            i += 1;
        // } else {
            // j -= 1;
        // }
        // toggleDirection = !toggleDirection;
        // iter += 1;
        // if (iter%2 == 0 && j-i == n) break;
        std::cout << "delta inside move singularity loop: " << delta << std::endl;
    }
    return;*/
    while (i < j) {
        std::cout << "i: " << i << " j: " << j << std::endl; 
        std::vector<size_t> tfp_a, tfp_b, verticesToAvoid;
        SingularityLink l_a, l_b;
        bool diagPair_a, diagPair_b;
        PrototypeSaveMesh(std::vector<SingularityLink>{l}, "test");
        if (i+1 == j) {
            // std::cout << "Getting pair link" << std::endl;
            tfp_a.push_back(l.linkVids.at(i));
            tfp_a.push_back(l.linkVids.at(i+1));
            verticesToAvoid = tfp_a;
            if (mesh->V.at(tfp_a.at(0)).N_Vids.size() != 3) std::reverse(tfp_a.begin(), tfp_a.end());
            diagPair_a = Contains(mesh->V.at(tfp_a.at(0)).N_Vids, tfp_a.at(1)) ? false : true;
            l_a = PrototypePairLink(tfp_a, verticesToAvoid, diagPair_a);
            PrototypeSaveMesh(std::vector<SingularityLink>{l_a}, "test2");
            // std::cout << "Got pair ids a link" << std::endl;
            PrototypeMovePair(tfp_a, l_a, diagPair_a);
            // std::cout << "Moving pair a" << std::endl;
            // links.push_back(l_a);
            i += 1;
            continue;
        }
        diagPair_a = PrototypePairIds(tfp_a, l.linkVids, i, i+1);
        std::cout << "Got pair ids a, NVids at i+1: " << l.linkVids.at(i+1) << " " << mesh->V.at(l.linkVids.at(i+1)).N_Vids.size() << std::endl;
        diagPair_b = PrototypePairIds(tfp_b, l.linkVids, j, j-1);
        std::cout << "Got pair ids b, NVids at j-1: " <<  l.linkVids.at(j-1) << " " << mesh->V.at(l.linkVids.at(j-1)).N_Vids.size() << std::endl;
        mu->AddContents(verticesToAvoid, std::vector<size_t>{l.linkVids.at(i+1), l.linkVids.at(j-1)});
        mu->AddContents(verticesToAvoid, tfp_a);
        mu->AddContents(verticesToAvoid, tfp_b);
        // std::cout << "verticesToAvoid: ";
        // for (auto aid: verticesToAvoid) std::cout << aid << " ";
        // std::cout << std::endl;
        if (!tfp_a.empty()) {
            l_a = PrototypePairLink(tfp_a, verticesToAvoid, diagPair_a);
            std::cout << "Got pair ids a link" << std::endl;
        }
        if (!tfp_b.empty()) {
            l_b = PrototypePairLink(tfp_b, verticesToAvoid, diagPair_b);
            std::cout << "Got pair ids b link" << std::endl;
        }
        // PrototypeSaveMesh(std::vector<SingularityLink>{l_a, l_b}, "test2");
        if (i+2 == j) {
            // std::cout << "Moving pair a, l_a: " << l_a.linkVids.front() << " " << l_a.linkVids.back() << std::endl;
            PrototypeMovePair(tfp_a, l_a, diagPair_a);
            std::cout << "Moved pair a" << std::endl;
            // std::cout << "Moving pair b, l_b: " << l_b.linkVids.front() << " " << l_b.linkVids.back() << std::endl;
            PrototypeMovePair(tfp_b, l_b, diagPair_b);
            std::cout << "Moved pair b" << std::endl;
            i += 1;
            j -= 1;
        } else if (!l_a.linkVids.empty() && !l_b.linkVids.empty() && fabs(l_a.rank) <= fabs(l_b.rank)) {
            // links.push_back(l_a);
            // std::cout << "Resolving pair ids b" << std::endl;
            // std::cout << "Three Five ids to resolve: " << tfp_b.at(0) << " " << mesh->V.at(tfp_b.at(0)).N_Vids.size() << " " << tfp_b.at(1) << " " << mesh->V.at(tfp_b.at(1)).N_Vids.size() << std::endl;
            // std::cout << "id at j = " << j << " " << l.linkVids.at(j) << " id at j-1 = " << j-1 << " " << l.linkVids.at(j-1) << std::endl;
            PrototypeResolvePairIds(tfp_b, l.linkVids, j, j-1);
            // std::cout << "id at j = " << j << " " << l.linkVids.at(j) << " " << mesh->V.at(l.linkVids.at(j)).N_Vids.size() << std::endl;
            // std::cout << "Moving pair a" << std::endl;
            // std::cout << "Three Five pair to move: " << tfp_a.at(0) << " " << mesh->V.at(tfp_a.at(0)).N_Vids.size() << " " << tfp_a.at(1) << " " << mesh->V.at(tfp_a.at(1)).N_Vids.size() << std::endl;
            // std::cout << "link front: " << l_a.linkVids.front() << " " << l_a.linkVids.back() << " " << mesh->V.at(l_a.linkVids.back()).N_Vids.size() << std::endl;
            PrototypeMovePair(tfp_a, l_a, diagPair_a); 
            i += 1;
        } else if (!l_a.linkVids.empty() && !l_b.linkVids.empty() && fabs(l_a.rank) > fabs(l_b.rank)) {
            // links.push_back(l_b);
            // std::cout << "Resolving pair ids a" << std::endl;
            // std::cout << "Three Five ids to resolve: " << tfp_a.at(0) << " " << mesh->V.at(tfp_a.at(0)).N_Vids.size() << " " << tfp_a.at(1) << " " << mesh->V.at(tfp_a.at(1)).N_Vids.size() << std::endl;
            // std::cout << "id at i = " << i << " " << l.linkVids.at(i) << " id at i+1 = " << i+1 << " " << l.linkVids.at(i+1) << std::endl;
            PrototypeResolvePairIds(tfp_a, l.linkVids, i, i+1);
            // std::cout << "id at i: " << i << " " << l.linkVids.at(i) << " " << mesh->V.at(l.linkVids.at(i)).N_Vids.size() << std::endl;
            // std::cout << "Moving pair b" << std::endl;
            // std::cout << "Three Five pair to move: " << tfp_b.at(0) << " " << mesh->V.at(tfp_b.at(0)).N_Vids.size() << " " << tfp_b.at(1) << " " << mesh->V.at(tfp_b.at(1)).N_Vids.size() << std::endl;
            // std::cout << "link front: " << l_b.linkVids.front() << " " << l_b.linkVids.back() << " " << mesh->V.at(l_a.linkVids.back()).N_Vids.size() << std::endl;
            PrototypeMovePair(tfp_b, l_b, diagPair_b);
            j -= 1;
        } else if (!l_a.linkVids.empty() && l_b.linkVids.empty()) {
            // std::cout << "Moving pair a" << std::endl;
            PrototypeMovePair(tfp_a, l_a, diagPair_a);
            i += 1;
        } else if (l_a.linkVids.empty() && !l_b.linkVids.empty()) {
            // std::cout << "Moving pair b" << std::endl;
            PrototypeMovePair(tfp_b, l_b, diagPair_b);
            j -= 1;
        } else if (l_a.linkVids.empty() && l_b.linkVids.empty()) {
            break;
        }
    }
    return;
    size_t frontSingularity = l.linkVids.front();
    size_t backSingularity = l.linkVids.back();
    // std::cout << "link size: " << l.linkVids.size() << std::endl;
    for (int i = 1; i < l.linkVids.size(); i++) {
        auto& toMove = mesh->V.at(l.linkVids.at(i-1));
        auto& dest = mesh->V.at(l.linkVids.at(i));
        // std::cout << "toMove: " << toMove.id << " " << toMove.N_Vids.size() << " dest: " << dest.id << " " << dest.N_Vids.size() << std::endl;
        bool moved = false;
        for (auto eid: toMove.N_Eids) {
            auto& e = mesh->E.at(eid);
            // std::cout << "edge vids: " << e.Vids.at(0) << " " << e.Vids.at(1) << std::endl;
            if (Contains(e.Vids, dest.id)) {
                // std::cout << "dest is in edge of singularity" << std::endl;
                if (toMove.N_Vids.size() == 3) {
                    auto& f = mesh->F.at(mu->GetDifference(toMove.N_Fids, e.N_Fids).at(0));
                    int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), toMove.id));
                    size_t threeId = f.Vids.at((idx+2)%f.Vids.size());
                    size_t fiveId = f.Vids.at((idx+1)%f.Vids.size());
                    std::shared_ptr<SimplificationOperation> dc = std::make_shared<DiagonalCollapse>(*mesh, *mu, *smoother, f.id, (idx+1)%f.Vids.size(), (idx+3)%f.Vids.size());
                    dc->PerformOperation();
                    
                    std::shared_ptr<ThreeFivePair> tfp = std::make_shared<ThreeFivePair>(*mesh, *mu, *smoother, threeId, fiveId);
                    std::vector<SingularityLink> links;
                    TraceSingularityLinks(threeId, links);
                    TraceSingularityLinks(fiveId, links);
                    SelectDirectPairLink(std::vector<size_t>{threeId, fiveId}, links, std::vector<size_t>{threeId, fiveId, dest.id, l.linkVids.back()});
                    std::vector<size_t> path(links.at(0).linkVids.begin()+1, links.at(0).linkVids.end());
                    for (int pid = 0; pid < path.size(); pid++) {
                        tfp->Move(path.at(pid), delta);
                    }
                    moved = true;
                    break;
                } else if (toMove.N_Vids.size() == 5) {
                    // std::cout << "dest is present in five N_Vids " << std::endl;
                    std::vector<size_t> tempEdges;
                    for (auto fid: e.N_Fids) {
                        mu->AddContents(tempEdges, mu->GetDifference(mu->GetIntersection(mesh->F.at(fid).Eids, toMove.N_Eids), std::vector<size_t>{eid}));
                    }
                    std::vector<size_t> faces = e.N_Fids;
                    for (auto tmp: tempEdges) mu->AddContents(faces, mesh->E.at(tmp).N_Fids);
                    auto& f = mesh->F.at(mu->GetDifference(toMove.N_Fids, faces).at(0));
                    size_t threeId = toMove.id;
                    size_t fiveId = GetFaceV(toMove.id, f.id, 2);
                    // std::cout << "tempEdges: " << tempEdges.size() << std::endl;
                    std::shared_ptr<SimplificationOperation> vs = std::make_shared<VertexSplit>(*mesh, *mu, *smoother, toMove.id, tempEdges);
                    vs->PerformOperation();
                    // std::cout << "performed vertex split" << std::endl;

                    std::shared_ptr<ThreeFivePair> tfp = std::make_shared<ThreeFivePair>(*mesh, *mu, *smoother, threeId, fiveId);
                    // std::cout << "made direct pair" << std::endl;
                    std::vector<SingularityLink> links;
                    TraceSingularityLinks(threeId, links);
                    PrototypeSaveMesh(links, "test_a");
                    // std::cout << "traced three links" << std::endl;
                    TraceSingularityLinks(fiveId, links);
                    PrototypeSaveMesh(links, "test_b");
                    // std::cout << "traced five links" << std::endl;
                    SelectDirectPairLink(std::vector<size_t>{threeId, fiveId}, links, std::vector<size_t>{threeId, fiveId, dest.id, l.linkVids.back()});
                    PrototypeSaveMesh(links, "test");
                    // std::cout << "selected links" << std::endl;
                    std::vector<size_t> path(links.at(0).linkVids.begin()+1, links.at(0).linkVids.end());
                    for (int pid = 0; pid < path.size(); pid++) {
                        tfp->Move(path.at(pid), delta);
                    }
                    moved = true;
                    break;
                }
            }
        }
        if (moved) continue;
        for (auto fid: toMove.N_Fids) {
            auto& f = mesh->F.at(fid);
            if (GetFaceV(toMove.id, f.id, 2) == dest.id) {
                if (toMove.N_Vids.size() == 3) {
                    int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), toMove.id));
                    size_t threeId = mu->GetDifference(toMove.N_Vids, mu->GetIntersection(toMove.N_Vids, f.Vids)).at(0);
                    size_t fiveId = f.Vids.at((idx+1)%f.Vids.size());
                    std::shared_ptr<SimplificationOperation> dc = std::make_shared<DiagonalCollapse>(*mesh, *mu, *smoother, f.id, (idx+1)%f.Vids.size(), (idx+3)%f.Vids.size());
                    dc->PerformOperation();
                    
                    std::shared_ptr<DiagonalThreeFivePair> tfp = std::make_shared<DiagonalThreeFivePair>(*mesh, *mu, *smoother, threeId, fiveId);
                    std::vector<SingularityLink> links;
                    TraceSingularityLinks(threeId, links);
                    TraceSingularityLinks(fiveId, links);
                    SelectDiagonalPairLink(std::vector<size_t>{threeId, fiveId}, links, std::vector<size_t>{threeId, fiveId, dest.id, l.linkVids.back()});
                    std::vector<size_t> path(links.at(0).linkVids.begin()+1, links.at(0).linkVids.end());
                    for (int pid = 0; pid < path.size(); pid++) {
                        tfp->Move(path.at(pid), delta);
                    }
                    break;
                } else if (toMove.N_Vids.size() == 5) {
                    std::vector<size_t> diffEdges = mu->GetIntersection(toMove.N_Eids, f.Eids);
                    for (auto diffE: diffEdges) {
                        mu->AddContents(diffEdges, mu->GetIntersection(toMove.N_Eids, mesh->F.at(mesh->E.at(diffE).N_Fids.at(0)).Eids));
                        mu->AddContents(diffEdges, mu->GetIntersection(toMove.N_Eids, mesh->F.at(mesh->E.at(diffE).N_Fids.at(1)).Eids));
                    }
                    auto& e = mesh->E.at(mu->GetDifference(toMove.N_Eids, diffEdges).at(0));
                    std::vector<size_t> tempEdges;
                    for (auto efid: e.N_Fids) {
                        mu->AddContents(tempEdges, mu->GetDifference(mu->GetIntersection(mesh->F.at(efid).Eids, toMove.N_Eids), std::vector<size_t>{e.id}));
                    }
                    std::vector<size_t> faces = e.N_Fids;
                    for (auto tmp: tempEdges) mu->AddContents(faces, mesh->E.at(tmp).N_Fids);
                    size_t threeId = toMove.id;
                    size_t fiveId = e.Vids.at(0) == toMove.id ? e.Vids.at(1) : e.Vids.at(0);
                    std::shared_ptr<SimplificationOperation> vs = std::make_shared<VertexSplit>(*mesh, *mu, *smoother, toMove.id, tempEdges);
                    vs->PerformOperation();

                    std::shared_ptr<DiagonalThreeFivePair> tfp = std::make_shared<DiagonalThreeFivePair>(*mesh, *mu, *smoother, threeId, fiveId);
                    std::vector<SingularityLink> links;
                    TraceSingularityLinks(threeId, links);
                    TraceSingularityLinks(fiveId, links);
                    SelectDiagonalPairLink(std::vector<size_t>{threeId, fiveId}, links, std::vector<size_t>{threeId, fiveId, dest.id, l.linkVids.back()});
                    std::vector<size_t> path(links.at(0).linkVids.begin()+1, links.at(0).linkVids.end());
                    for (int pid = 0; pid < path.size(); pid++) {
                        tfp->Move(path.at(pid), delta);
                    }
                    break;
                }
            }
        }
    }
}

SingularityLink SemiGlobalSimplifier::PrototypeGetLink(size_t vid) {
    std::vector<SingularityLink> links;
    TraceSingularityLinks(vid, links, true);
    SelectLinks(links);
    return links.at(0);   
}

bool SemiGlobalSimplifier::PrototypePairIds(std::vector<size_t>& threeFiveIds, std::vector<size_t>& vids, size_t toMoveIdx, size_t destIdx) {
    std::cout << "Getting pair ids" << std::endl; 
    auto& toMove = mesh->V.at(vids.at(toMoveIdx));
    auto& dest = mesh->V.at(vids.at(destIdx));
    if (toMove.N_Vids.size() != 3 && toMove.N_Vids.size() != 5) return false;
    std::cout << "toMove: " << toMove.id << " " << toMove.N_Vids.size() << std::endl;
    std::cout << "dest: " << dest.id << " " << dest.N_Vids.size() << std::endl;
     bool gotDirectPair = false;
    for (auto eid: toMove.N_Eids) {
        auto& e = mesh->E.at(eid);
        std::cout << "Looking for direct pair" << std::endl;
        if (Contains(e.Vids, dest.id)) {
            if (toMove.N_Vids.size() == 3) {
                auto& f = mesh->F.at(mu->GetDifference(toMove.N_Fids, e.N_Fids).at(0));
                int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), toMove.id));
                threeFiveIds.push_back(f.Vids.at((idx+2)%f.Vids.size()));
                threeFiveIds.push_back(f.Vids.at((idx+1)%f.Vids.size()));
                std::cout << "Identified three five pair: " << threeFiveIds.at(0) << " " << threeFiveIds.at(1) << std::endl;
                if (mesh->V.at(threeFiveIds.at(0)).N_Vids.size() < 4 || mu->IsSharpFeature(threeFiveIds.at(0))) {
                    threeFiveIds.clear();
                    return false;
                }
                vids.at(toMoveIdx) = threeFiveIds.at(0);
                std::shared_ptr<SimplificationOperation> dc = std::make_shared<DiagonalCollapse>(*mesh, *mu, *smoother, f.id, (idx+1)%f.Vids.size(), (idx+3)%f.Vids.size());
                dc->PerformOperation();
                gotDirectPair = true;
                break;
            } else if (toMove.N_Vids.size() == 5) {
                std::vector<size_t> tempEdges;
                for (auto fid: e.N_Fids) {
                    mu->AddContents(tempEdges, mu->GetDifference(mu->GetIntersection(mesh->F.at(fid).Eids, toMove.N_Eids), std::vector<size_t>{eid}));
                }
                std::vector<size_t> faces = e.N_Fids;
                for (auto tmp: tempEdges) mu->AddContents(faces, mesh->E.at(tmp).N_Fids);
                auto& f = mesh->F.at(mu->GetDifference(toMove.N_Fids, faces).at(0));
                threeFiveIds.push_back(toMove.id);
                threeFiveIds.push_back(GetFaceV(toMove.id, f.id, 2));
                std::cout << "Identified three five pair: " << threeFiveIds.at(0) << " " << threeFiveIds.at(1) << std::endl;
                if (mesh->V.at(threeFiveIds.at(1)).N_Vids.size() >= 5 || mu->IsSharpFeature(threeFiveIds.at(1))) {
                    threeFiveIds.clear();
                    return false;
                }
                vids.at(toMoveIdx) = threeFiveIds.at(1);
                std::shared_ptr<SimplificationOperation> vs = std::make_shared<VertexSplit>(*mesh, *mu, *smoother, toMove.id, tempEdges);
                vs->PerformOperation();
                gotDirectPair = true;
                break;
            }
        }
    }
    if (gotDirectPair) {
        if (mesh->V.at(threeFiveIds.at(0)).N_Vids.size() == 4 || mesh->V.at(threeFiveIds.at(1)).N_Vids.size() == 4) {
            threeFiveIds.clear();
        }
        return false;
    }
    for (auto fid: toMove.N_Fids) {
        auto& f = mesh->F.at(fid);
        std::cout << "Looking for diagonal pair" << std::endl;
        std::cout << "f Vids: ";
        for (auto fvid: f.Vids) std::cout << fvid << " ";
        std::cout << std::endl;
        std::cout << "f Eids: ";
        for (auto feid: f.Eids) std::cout << feid << " ";
        std::cout << std::endl;
        if (GetDiagonalV(toMove.id, f.id) == dest.id) {
            if (toMove.N_Vids.size() == 3) {
                int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), toMove.id));
                threeFiveIds.push_back(mu->GetDifference(toMove.N_Vids, mu->GetIntersection(toMove.N_Vids, f.Vids)).at(0));
                threeFiveIds.push_back(f.Vids.at((idx+1)%f.Vids.size()));
                std::cout << "Identified three five pair: " << threeFiveIds.at(0) << " " << threeFiveIds.at(1) << std::endl;
                if (mesh->V.at(threeFiveIds.at(0)).N_Vids.size() < 4 || mu->IsSharpFeature(threeFiveIds.at(0))) {
                    threeFiveIds.clear();
                    return false;
                }
                std::shared_ptr<SimplificationOperation> dc = std::make_shared<DiagonalCollapse>(*mesh, *mu, *smoother, f.id, (idx+1)%f.Vids.size(), (idx+3)%f.Vids.size());
                dc->PerformOperation();
                std::cout << "toMove: " << toMove.id << " " << toMove.N_Vids.size() << " threeFiveIds at 1: " << threeFiveIds.at(1) << " " << mesh->V.at(threeFiveIds.at(1)).N_Vids.size() << std::endl;
                if (mesh->V.at(threeFiveIds.at(1)).N_Vids.empty()) threeFiveIds.at(1) = toMove.id;
                vids.at(toMoveIdx) = threeFiveIds.at(0);
                break;
            } else if (toMove.N_Vids.size() == 5) {
                std::cout << "five Edges: " << std::endl;
                for (auto teid: toMove.N_Eids) {
                    std::cout << teid << " " << mesh->E.at(teid).Vids.at(0) << " " << mesh->E.at(teid).Vids.at(1) << std::endl;
                }
                std::vector<size_t> diffEdges = mu->GetIntersection(toMove.N_Eids, f.Eids);
                std::cout << "diffEdges: ";
                for (auto teid: diffEdges) std::cout << teid << " ";
                std::cout << std::endl;
                for (auto diffEid: diffEdges) {
                    auto& diffE = mesh->E.at(diffEid);
                    auto diffFid = diffE.N_Fids.at(0) == f.id ? diffE.N_Fids.at(1) : diffE.N_Fids.at(0);
                    auto& diffF = mesh->F.at(diffFid);
                    std::cout << "diffF Eids: ";
                    for (auto feid: diffF.Eids) std::cout << feid << " ";
                    std::cout << std::endl;
                    mu->AddContents(diffEdges, mu->GetIntersection(toMove.N_Eids, diffF.Eids));
                    // mu->AddContents(diffEdges, mu->GetIntersection(toMove.N_Eids, diffE.N_Fids.at(0).Eids));
                    // mu->AddContents(diffEdges, mu->GetIntersection(toMove.N_Eids, mesh->E.at(diffE).N_Fids.at(1).Eids));
                }
                std::cout << "diffEdges: ";
                for (auto teid: diffEdges) std::cout << teid << " ";
                std::cout << std::endl;
                auto& e = mesh->E.at(mu->GetDifference(toMove.N_Eids, diffEdges).at(0));
                std::vector<size_t> tempEdges;
                for (auto efid: e.N_Fids) {
                    mu->AddContents(tempEdges, mu->GetDifference(mu->GetIntersection(mesh->F.at(efid).Eids, toMove.N_Eids), std::vector<size_t>{e.id}));
                }
                std::vector<size_t> faces = e.N_Fids;
                for (auto tmp: tempEdges) mu->AddContents(faces, mesh->E.at(tmp).N_Fids);
                std::cout << "tempEdges: ";
                for (auto teid: tempEdges) std::cout << teid << " ";
                std::cout << std::endl;
                std::cout << "Edge: " << e.id << " " << e.Vids.at(0) << " " << e.Vids.at(1) << std::endl;
                threeFiveIds.push_back(toMove.id);
                threeFiveIds.push_back((e.Vids.at(0) == toMove.id ? e.Vids.at(1) : e.Vids.at(0)));
                std::cout << "Identified three five pair: " << threeFiveIds.at(0) << " " << threeFiveIds.at(1) << std::endl;
                if (mesh->V.at(threeFiveIds.at(1)).N_Vids.size() >= 5 || mu->IsSharpFeature(threeFiveIds.at(1))) {
                    threeFiveIds.clear();
                    return false;
                }
                vids.at(toMoveIdx) = threeFiveIds.at(1);
                std::shared_ptr<SimplificationOperation> vs = std::make_shared<VertexSplit>(*mesh, *mu, *smoother, toMove.id, tempEdges);
                vs->PerformOperation();
                break;
            }
        }
    }
    if (!threeFiveIds.empty() && (mesh->V.at(threeFiveIds.at(0)).N_Vids.size() == 4 || mesh->V.at(threeFiveIds.at(1)).N_Vids.size() == 4)) {
        threeFiveIds.clear();
    }
    return true;  
}

void SemiGlobalSimplifier::PrototypeResolvePairIds(std::vector<size_t>& threeFiveIds, std::vector<size_t>& vids, size_t toMoveIdx, size_t destIdx) {
    if (threeFiveIds.empty()) return;
    if (mesh->V.at(threeFiveIds.at(0)).N_Vids.empty() || mesh->V.at(threeFiveIds.at(1)).N_Vids.empty()) return;
    auto& three = mesh->V.at(threeFiveIds.at(0));
    auto& five = mesh->V.at(threeFiveIds.at(1));
    auto& dest = mesh->V.at(vids.at(destIdx));
    std::cout << "Resolving pair ids: " << three.id << " " << three.N_Vids.size() << " " << five.id << " " << five.N_Vids.size() << std::endl;
    std::cout << "dest: " << dest.id << " " << dest.N_Vids.size() << std::endl;
    if (dest.N_Vids.size() == 5) {
        for (auto fid: three.N_Fids) {
            std::cout << "Checking three NFids" << std::endl;
            auto& f = mesh->F.at(fid);
            int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), three.id));
            size_t diagId = GetDiagonalV(three.id, f.id);
            if (diagId == dest.id || diagId == five.id) {
                std::cout << "found diag id: " << diagId << std::endl;
                std::cout << "setting toMoveIdx: " << toMoveIdx << std::endl;
                vids.at(toMoveIdx) = f.Vids.at((idx+1)%f.Vids.size());
                std::shared_ptr<SimplificationOperation> dc = std::make_shared<DiagonalCollapse>(*mesh, *mu, *smoother, f.id, (idx+1)%f.Vids.size(), (idx+3)%f.Vids.size());
                dc->PerformOperation();
                break;
            }
        }
    } else if (dest.N_Vids.size() == 3) {
        std::vector<size_t> edges = mu->GetIntersection(five.N_Eids, three.N_Eids);
        mu->AddContents(edges, mu->GetIntersection(five.N_Eids, dest.N_Eids));
        auto& e = mesh->E.at(edges.at(0));
        std::vector<size_t> tempEdges = mu->GetDifference(mu->GetIntersection(mesh->F.at(e.N_Fids.at(0)).Eids, five.N_Eids), edges);
        mu->AddContents(tempEdges, mu->GetDifference(mu->GetIntersection(mesh->F.at(e.N_Fids.at(1)).Eids, five.N_Eids), edges));
        std::cout << "setting toMoveIdx: " << toMoveIdx << std::endl;
        vids.at(toMoveIdx) = five.id;
        std::shared_ptr<SimplificationOperation> vs = std::make_shared<VertexSplit>(*mesh, *mu, *smoother, five.id, tempEdges);
        vs->PerformOperation();
    }
}

void SemiGlobalSimplifier::PrototypeMovePair(std::vector<size_t>& threeFiveIds, SingularityLink& l, bool isDiagonal) {
    if (threeFiveIds.empty() || l.linkVids.empty()) return;
    if (mesh->V.at(threeFiveIds.at(0)).N_Vids.empty() || mesh->V.at(threeFiveIds.at(1)).N_Vids.empty()) return;
    // if (!PrototypeIsLinkValid(l)) return;
    // PrototypeSaveMesh(std::vector<SingularityLink>{l}, "test_move_pair");
    std::cout << "threeFiveIds: " << threeFiveIds.size() << " three: " << threeFiveIds.at(0) << " " << mesh->V.at(threeFiveIds.at(0)).N_Vids.size() << " five: " << threeFiveIds.at(1) << " " << mesh->V.at(threeFiveIds.at(1)).N_Vids.size() << std::endl;
    std::cout << "link is valid, front: " << l.linkVids.front() << " back: " << l.linkVids.back() << std::endl;
    if (isDiagonal) {
        std::cout << "moving diagonal pair" << std::endl;
        std::shared_ptr<DiagonalThreeFivePair> tfp = std::make_shared<DiagonalThreeFivePair>(*mesh, *mu, *smoother, threeFiveIds.at(0), threeFiveIds.at(1));
        for (int i = 1; i < l.linkVids.size(); i++) {
            std::cout << i << " ";
            tfp->Move(l.linkVids.at(i), delta);
            SingularityLink t_l;
            t_l.linkVids = std::vector<size_t>(l.linkVids.begin()+i+1,l.linkVids.end());
            PrototypeSaveMesh(std::vector<SingularityLink>{t_l}, "test2"+std::to_string(i+1));
            // PrototypeSaveMesh(std::vector<SingularityLink>{l}, "test_move_pair");
        }
        std::cout << std::endl;
    } else {
        std::cout << "moving direct pair" << std::endl;
        std::shared_ptr<ThreeFivePair> tfp = std::make_shared<ThreeFivePair>(*mesh, *mu, *smoother, threeFiveIds.at(0), threeFiveIds.at(1));
        for (int i = 1; i < l.linkVids.size(); i++) {
            std::cout << i << " ";
            tfp->Move(l.linkVids.at(i), delta);
            SingularityLink t_l;
            t_l.linkVids = std::vector<size_t>(l.linkVids.begin()+i+1,l.linkVids.end());
            PrototypeSaveMesh(std::vector<SingularityLink>{t_l}, "test2"+std::to_string(i+1));
            // PrototypeSaveMesh(std::vector<SingularityLink>{l}, "test_move_pair");
        }
        std::cout << std::endl;
    }
}

bool SemiGlobalSimplifier::PrototypeIsLinkValid(SingularityLink& l, std::vector<size_t> verticesToAvoid) {
    if (l.linkVids.empty()) return false;
    if ((mesh->V.at(l.linkVids.front()).N_Vids.size() != 3 && mesh->V.at(l.linkVids.front()).N_Vids.size() != 5) ||
        (mesh->V.at(l.linkVids.back()).N_Vids.size() != 3 && mesh->V.at(l.linkVids.back()).N_Vids.size() != 5)) {
            return false;
    }
    if (l.linkVids.size() < 3) return true;
    for (int i = 1; i < l.linkVids.size()-1; i++) {
        auto& v = mesh->V.at(l.linkVids.at(i));
        if (v.N_Vids.size() != 4) return false;
        // if (!verticesToAvoid.empty() && Contains(v.N_Vids, verticesToAvoid)) return false;
        bool foundPrev, foundNext = false;
        for (auto fid: v.N_Fids) {
            if (Contains(mesh->F.at(fid).Vids, l.linkVids.at(i-1))) {
                foundPrev = true;
            }
            if (Contains(mesh->F.at(fid).Vids, l.linkVids.at(i+1))) {
                foundNext = true;
            }
        }
        if (!foundPrev || !foundNext) return false;
    }
    return true;
}

SingularityLink SemiGlobalSimplifier::PrototypePairLink(std::vector<size_t> threeFiveIds, std::vector<size_t>& verticesToAvoid, bool isDiagonal) {
    std::vector<SingularityLink> links;
    SingularityLink l;
    if (threeFiveIds.empty()) return l;
    if (mesh->V.at(threeFiveIds.at(0)).N_Vids.empty() || mesh->V.at(threeFiveIds.at(1)).N_Vids.empty()) return l;
    // std::cout << "three: " << threeFiveIds.at(0) << " " << mesh->V.at(threeFiveIds.at(0)).N_Vids.size() << " five: " << threeFiveIds.at(1) << " " << mesh->V.at(threeFiveIds.at(1)).N_Vids.size() << std::endl; 
    TraceSingularityLinks(threeFiveIds.at(0), links);
    TraceSingularityLinks(threeFiveIds.at(1), links);
    PrototypeSaveMesh(links, "test_pair_links");
    if (isDiagonal) {
        // std::cout << "Selecting Diagonal Pair link" << std::endl;
        SelectDiagonalPairLink(threeFiveIds, links, verticesToAvoid);
    } else {
        // std::cout << "Selecting Direct Pair link" << std::endl;
        SelectDirectPairLink(threeFiveIds, links, verticesToAvoid);
    }
    if (!links.empty()) {
        l = links.at(0);
    }
    if (!l.linkVids.empty()) mu->AddContents(verticesToAvoid, std::vector<size_t>{l.linkVids.back()});
    return l;
}

void SemiGlobalSimplifier::SelectLinks(std::vector<SingularityLink>& links) {
    SingularityLink link;
    std::priority_queue<SingularityLink, std::vector<SingularityLink>, LinkComparator> q;
    for (auto& l: links) {
        int valenceToCheck = mesh->V.at(l.linkVids.front()).N_Vids.size() == 3 ? 5 : 3;
        if ((mesh->V.at(l.linkVids.back()).N_Vids.size() != 3 && mesh->V.at(l.linkVids.back()).N_Vids.size() != 5) || 
        mesh->V.at(l.linkVids.back()).N_Vids.size() != valenceToCheck || doesCrossBoundary(l.linkVids, true)) continue;
        // if (!PrototypeIsLinkValid(l)) continue;
        l.rank = mesh->V.at(l.linkVids.front()).N_Vids.size() == 3 ? -1 * l.rank : 1 * l.rank;
        q.push(l);
    }
    links.clear();
    if (!q.empty()) {
        link = q.top();
    }
    links.push_back(link);
    // q.pop();
    // while (!q.empty()) {
    //     auto& l = q.top();
    //     if (l.linkVids.back() != links.at(0).linkVids.back()) {
    //         std::cout << links.at(0).linkVids.size() << " " << l.linkVids.size() << std::endl;
    //         links.push_back(l);
    //         break;
    //     }
    //     q.pop();
    // }
}

void SemiGlobalSimplifier::SelectDirectPairLink(std::vector<size_t> threeFiveIds, std::vector<SingularityLink>& links, std::vector<size_t> verticesToAvoid) {
    auto& three = mesh->V.at(threeFiveIds.at(0));
    auto& five = mesh->V.at(threeFiveIds.at(1));
    size_t fid = mu->GetIntersection(three.N_Fids, five.N_Fids).at(0);
    auto& pairEdge = mesh->E.at(mu->GetIntersection(three.N_Eids, five.N_Eids).at(0));
    std::vector<size_t> threeEdges = mu->GetDifference(three.N_Eids, std::vector<size_t>{pairEdge.id});
    std::vector<size_t> fiveEdges = mu->GetDifference(five.N_Eids, std::vector<size_t>{pairEdge.id});

    std::priority_queue<SingularityLink, std::vector<SingularityLink>, LinkComparator> q;
    for (auto& l: links) {
        // int valenceToCheck = mesh->V.at(l.linkVids.front()).N_Vids.size() == 3 ? 5 : 3;
        // if (Contains(verticesToAvoid, l.linkVids.back()) || mesh->V.at(l.linkVids.back()).N_Vids.size() != valenceToCheck) continue;
        if (Contains(verticesToAvoid, l.linkVids.back()) || mesh->V.at(l.linkVids.back()).isBoundary) continue;
        // if (!PrototypeIsLinkValid(l, verticesToAvoid)) continue;
        if ((mesh->V.at(l.linkVids.back()).N_Vids.size() != 3 && mesh->V.at(l.linkVids.back()).N_Vids.size() != 5)
        || doesCrossBoundary(l.linkVids, true)) continue;
        if (l.b > 0) GetDiagonalPath(l);
        if (l.linkVids.empty()) continue;
        if (l.rots > 0) {
            if (Contains(threeEdges, l.linkEids.at(0))) {
                if ((l.rots == 1 && GetEdgeRots(pairEdge.id, l.linkEids.at(0)) > 1) || (l.rots > 1 && GetEdgeRots(pairEdge.id, l.linkEids.at(0)) == 1)) {
                    if (Contains(three.N_Vids, l.linkVids.at(1))) {
                        l.rank = -1*(l.a+(2*l.b));
                    } else {
                        l.rank = -1*((2*l.a)+l.b);
                    }
                } else if ((l.rots == 1 && GetEdgeRots(pairEdge.id, l.linkEids.at(0)) == 1) || (l.rots > 1 && GetEdgeRots(pairEdge.id, l.linkEids.at(0)) > 1)) {
                    if (Contains(three.N_Vids, l.linkVids.at(1))) {
                        l.rank = -1*l.a;
                    } else {
                        l.rank = l.b-1;
                    }
                }
            } else if (Contains(fiveEdges, l.linkEids.at(0))) {
                if ((l.rots == 1 && GetEdgeRots(pairEdge.id, l.linkEids.at(0)) == 1) || (l.rots > 1 && GetEdgeRots(pairEdge.id, l.linkEids.at(0)) == 4)) {
                    if (Contains(five.N_Vids, l.linkVids.at(1))) {
                        l.rank = -1*((l.a-1)+(2*l.b));
                    } else {
                        l.rank = -1*(1+(2*(l.a-1))+l.b);
                    }
                } else if ((l.rots > 1 && GetEdgeRots(pairEdge.id, l.linkEids.at(0)) == 1) || (l.rots == 1 && GetEdgeRots(pairEdge.id, l.linkEids.at(0)) == 4)) {
                    if (Contains(five.N_Vids, l.linkVids.at(1))) {
                        l.rank = -1*(l.a-1);
                    } else {
                        l.rank = l.b-1;
                    }
                } else if ((l.rots == 1 && GetEdgeRots(pairEdge.id, l.linkEids.at(0)) == 2) || (l.rots > 1 && GetEdgeRots(pairEdge.id, l.linkEids.at(0)) == 3)) {
                    if (Contains(five.N_Vids, l.linkVids.at(1))) {
                        l.rank = l.a;
                    } else {
                        l.rank = -1*(l.b-1);
                    }
                } else if ((l.rots > 1 && GetEdgeRots(pairEdge.id, l.linkEids.at(0)) == 2) || (l.rots == 1 && GetEdgeRots(pairEdge.id, l.linkEids.at(0)) == 3)) {
                    if (Contains(five.N_Vids, l.linkVids.at(1))) {
                        l.rank = l.a + (2*l.b);
                    } else {
                        l.rank = (2*l.a) + l.b;
                    }
                }
            }
        } else {
            if (Contains(three.N_Vids, l.linkVids.at(1))) {
                l.rank = -1*l.a;    
            } else if (Contains(five.N_Vids, l.linkVids.at(1))) {
                for (auto eid: fiveEdges) {
                    if ((GetEdgeRots(pairEdge.id, l.linkEids.at(0)) == 1 || GetEdgeRots(pairEdge.id, l.linkEids.at(0)) == 4) && (mesh->E.at(eid).Vids.at(0) == l.linkVids.at(1) || mesh->E.at(eid).Vids.at(1) == l.linkVids.at(1))) {
                        l.rank = -1*(l.a-1);
                    } else if ((GetEdgeRots(pairEdge.id, l.linkEids.at(0)) == 2 || GetEdgeRots(pairEdge.id, l.linkEids.at(0)) == 3) && (mesh->E.at(eid).Vids.at(0) == l.linkVids.at(1) || mesh->E.at(eid).Vids.at(1) == l.linkVids.at(1))) {
                        l.rank = l.a;
                    }
                }
            }
        }
        // l.rank = mesh->F.size() + l.rank;
        q.push(l);
    }
    links.clear();
    if (!q.empty()) links.push_back(q.top());
}

void SemiGlobalSimplifier::SelectDiagonalPairLink(std::vector<size_t> threeFiveIds, std::vector<SingularityLink>& links, std::vector<size_t> verticesToAvoid) {
    auto& three = mesh->V.at(threeFiveIds.at(0));
    auto& five = mesh->V.at(threeFiveIds.at(1));
    size_t fid = mu->GetIntersection(three.N_Fids, five.N_Fids).at(0);
    std::vector<size_t> rotEdges = mu->GetIntersection(mesh->F.at(fid).Eids, three.N_Eids);
    size_t threeEdge = mu->GetDifference(three.N_Eids, rotEdges).at(0);
    std::vector<size_t> nRotEdges = mu->GetIntersection(mesh->F.at(fid).Eids, five.N_Eids);
    std::vector<size_t> fiveEdges;
    for (auto eid: nRotEdges) {
        size_t id = mesh->E.at(eid).N_Fids.at(0) == fid ? mesh->E.at(eid).N_Fids.at(1) : mesh->E.at(eid).N_Fids.at(0);
        mu->AddContents(rotEdges, mu->GetDifference(mu->GetIntersection(five.N_Eids, mesh->F.at(id).Eids), nRotEdges));
        mu->AddContents(fiveEdges, mu->GetIntersection(five.N_Eids, mesh->F.at(id).Eids));
    }
    size_t fiveEdge = mu->GetDifference(five.N_Eids, fiveEdges).at(0);
    mu->AddContents(nRotEdges, mu->GetDifference(five.N_Eids, rotEdges));
    mu->AddContents(nRotEdges, mu->GetDifference(three.N_Eids, rotEdges));
    std::priority_queue<SingularityLink, std::vector<SingularityLink>, LinkComparator> q;
    for (auto& l: links) {
        // int valenceToCheck = mesh->V.at(l.linkVids.front()).N_Vids.size() == 3 ? 5 : 3;
        // if (Contains(verticesToAvoid, l.linkVids.back()) || mesh->V.at(l.linkVids.back()).N_Vids.size() != valenceToCheck || mesh->V.at(l.linkVids.back()).isBoundary) continue;
        if (Contains(verticesToAvoid, l.linkVids.back()) || mesh->V.at(l.linkVids.back()).isBoundary) continue;
        // if (!PrototypeIsLinkValid(l, verticesToAvoid)) continue;
        if (l.b > 0) l.rots = GetEdgeRots(l.linkEids.at(l.a-1), l.linkEids.at(l.a));
        int start = 0;
        int end = l.a;
        if (Contains(nRotEdges, l.linkEids.at(0))) {
            if (l.b > 0 ) {
                start = l.a;
                end = l.linkVids.size();
            } else {
                start = 0;
                end = 0;
            }
            if (Contains(five.N_Eids, l.linkEids.at(0))) {
                l.rank = l.linkEids.at(0) == fiveEdge ? 2*l.a : -1*(2*(l.a-1));
            } else {
                l.rank = -1*(2*l.a);
            }
        }
        if ((mesh->V.at(l.linkVids.back()).N_Vids.size() != 3 && mesh->V.at(l.linkVids.back()).N_Vids.size() != 5)
        || doesCrossBoundary(std::vector<size_t>(l.linkVids.begin()+start, l.linkVids.begin()+end), true)) continue;
        if (Contains(rotEdges, l.linkEids.at(0))) {
            l.rank = 0;
            if (Contains(three.N_Eids, l.linkEids.at(0))) {
                if ((l.rots > 1 && GetEdgeRots(threeEdge, l.linkEids.at(0)) == 1) || (l.rots == 1 && GetEdgeRots(threeEdge, l.linkEids.at(0)) > 1)) {
                    l.rank = -1*(2*l.b);
                } else if ((l.rots == 1 && GetEdgeRots(threeEdge, l.linkEids.at(0)) == 1) || (l.rots > 1 && GetEdgeRots(threeEdge, l.linkEids.at(0)) > 1)) {
                    l.rank = 2*l.b;
                }
            } else if (Contains(five.N_Eids, l.linkEids.at(0))) {
                if ((l.rots > 1 && GetEdgeRots(fiveEdge, l.linkEids.at(0)) > 1) || (l.rots == 1 && GetEdgeRots(fiveEdge, l.linkEids.at(0)) == 1)) {
                    l.rank = -1*(2*l.b);
                } else if ((l.rots == 1 && GetEdgeRots(fiveEdge, l.linkEids.at(0)) > 1) || (l.rots > 1 && GetEdgeRots(fiveEdge, l.linkEids.at(0)) == 1)) {
                    l.rank = 2*l.b;
                }
            }
        }
        // l.rank = mesh->F.size() + l.rank;
        q.push(l);
    }
    links.clear();
    if (!q.empty()) links.push_back(q.top());
    std::vector<SingularityLink> selectedLinks;
    while (!q.empty()) {
        selectedLinks.push_back(q.top());
        q.pop();
    }
    PrototypeSaveMesh(selectedLinks, "test_diagonal_links");
}

bool SemiGlobalSimplifier::PrototypeIsLinkValid(SingularityGroup& s) {
    bool res = true;

    auto& l1 = s.l1;
    auto& l2 = s.l2;
    if (mesh->V.at(l1.frontId).N_Vids.size() == 3) {
        if (l1.b > 0 && l2.b > 0) {
            if (s.l1l2Volt == l1.volt && s.l2l1Volt == l2.volt) res = false;
            if (s.l1l2Volt != l1.volt && s.l2l1Volt != l2.volt) res = false;
            if (s.l1l2Volt == l1.volt && s.l2l1Volt != l2.volt && l1.a <= l2.a) res = false;
            if (s.l1l2Volt != l1.volt && s.l2l1Volt == l2.volt && l2.a <= l1.a) res = false;
        }
        else if (l1.b > 0 && l2.b == 0 && s.l1l2Volt == l1.volt && l1.a <= l2.a) res = false;
        else if (l1.b == 0 && l2.b > 0 && s.l2l1Volt == l2.volt && l2.a <= l1.a) res = false; 
    } else if (mesh->V.at(l1.frontId).N_Vids.size() == 5) {
        if (l1.b > 0 && l2.b > 0) {
            if (s.l1l2Volt == l1.volt && s.l2l1Volt == l2.volt) res = false;
            if (s.l1l2Volt == l1.volt && s.l2l1Volt != l2.volt && (l1.rots == 1 || l1.rots == 4) && l1.a <= l2.a) res = false;
        }
        else if (l1.b > 0 && l2.b == 0 && s.l1l2Volt == l1.volt && (l1.rots == 1 || l1.rots == 4) && l1.a <= l2.a) res = false;
        else if (l1.b == 0 && l2.b > 0 && s.l2l1Volt == l2.volt && (l1.rots == 1 || l1.rots == 4) && l2.a <= l1.a) res = false;
    }
    return res;
}

int SemiGlobalSimplifier::PrototypeGetElementPrediction(SingularityGroup& s) {
    int n = 0;
    
    auto& l1 = s.l1;
    auto& l2 = s.l2;
    if (mesh->V.at(l1.frontId).N_Vids.size() == 3) {
        n += l1.a + l1.b + l2.a + l2.b + (0.5 * (l1.a * (l1.a + 1))) + (0.5 * (l2.a * (l2.a + 1)));
        n *= -1;
        if (l1.volt == s.l1l2Volt && l1.b > 0) n += ((l1.a * l1.b) - (0.5 * (l1.b * (l1.b + 1))));
        if (l1.volt != s.l1l2Volt && l1.b > 0) n -= ((l1.a * l1.b) + (0.5 * (l1.b * (l1.b + 1))));
        if (l2.volt == s.l2l1Volt && l2.b > 0) n += ((l2.a * l2.b) - (0.5 * (l1.b * (l1.b + 1))));
        if (l2.volt != s.l2l1Volt && l2.b > 0) n -= ((l2.a * l2.b) + (0.5 * (l1.b * (l1.b + 1))));
    } else if (mesh->V.at(l1.frontId).N_Vids.size() == 5) {
        n += l1.a + l1.b + l2.a + l2.b;
        if (s.l1l2Rots == 1 || s.l1l2Rots == 4) {
            n -= (0.5 * (l1.a * (l1.a + 1)));
        } else {
            n += (0.5 * (l1.a * (l1.a + 1)));
        }
        if (l1.volt == s.l1l2Volt && l1.b > 0) n += ((l1.a * l1.b) - (0.5 * (l1.b * (l1.b + 1))));
        if (l1.volt != s.l1l2Volt && l1.b > 0) n += ((l1.a * l1.b) + (0.5 * (l1.b * (l1.b + 1))));
        if (l2.volt == s.l2l1Volt && l2.b > 0) n += ((l2.a * l2.b) - (0.5 * (l1.b * (l1.b + 1))));
        if (l2.volt != s.l2l1Volt && l2.b > 0) n += ((l2.a * l2.b) + (0.5 * (l1.b * (l1.b + 1))));
    }

    return n;
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

void SemiGlobalSimplifier::PrototypeSaveMesh(std::vector<SingularityLink>& links, std::string in) {
    int colorValue = 0;
    int ncolors = 15;
    std::vector<size_t> c_indices;
    std::vector<int> colors;
    
    std::vector<Edge> edges;
    for (auto& l: links) {
        if (l.linkVids.empty()) continue;
        for (int i = 1; i < l.linkVids.size(); i++) {
            Edge e;
            e.Vids = {l.linkVids.at(i-1), l.linkVids.at(i)};
            edges.push_back(e);
        }
        // c_indices.insert(c_indices.end(), l.linkEids.begin(), l.linkEids.end());
        // c_indices.insert(c_indices.end(), l.linkVids.begin(), l.linkVids.end());
        std::vector<int> a(l.linkVids.size()-1, (colorValue%ncolors));
        // std::vector<int> a(l.linkVids.size(), (colorValue%ncolors));
        colors.insert(colors.end(), a.begin(), a.end());
        colorValue += 1;
        a.clear();        
    }
    
    std::ofstream ofs(in+"_singularityLinks.vtk");
    ofs << "# vtk DataFile Version 3.0\n"
        << "singularityLinks" << ".vtk\n"
        << "ASCII\n\n"
        << "DATASET UNSTRUCTURED_GRID\n";
    ofs << "POINTS " << mesh->V.size() << " double\n";
    for (size_t i = 0; i < mesh->V.size(); i++) {
        ofs << std::fixed << std::setprecision(7) <<  mesh->V.at(i).x << " " <<  mesh->V.at(i).y << " " <<  mesh->V.at(i).z << "\n";
    }
    // ofs << "CELLS " << c_indices.size() << " " << 3 * c_indices.size() << std::endl;
    // for (size_t i = 0; i < c_indices.size(); i++) {
    //     auto& e = mesh->E.at(c_indices.at(i));
    //     if (e.Vids.empty()) continue;
    //     ofs << "2 " << e.Vids.at(0) << " " << e.Vids.at(1) << std::endl;
    // }
    ofs << "CELLS " << edges.size() << " " << 3 * edges.size() << std::endl;
    for (size_t i = 0; i < edges.size(); i++) {
        ofs << "2 " << edges.at(i).Vids.at(0) << " " << edges.at(i).Vids.at(1) << std::endl;
    }
    // ofs << "CELL_TYPES " << c_indices.size() << "\n";
    // for (size_t i = 0; i < c_indices.size(); i++) {
    //     ofs << "3" << std::endl;
    // }
    ofs << "CELL_TYPES " << edges.size() << "\n";
    for (size_t i = 0; i < edges.size(); i++) {
        ofs << "3" << std::endl;
    }

    // ofs << "CELL_DATA " << c_indices.size() << "\n";
    ofs << "CELL_DATA " << edges.size() << "\n";
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
        // std::cout << "Moving Pair" << std::endl;
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
    // std::cout << "three: " << mesh->V.at(threeFiveIds.at(0)).N_Vids.size() << " " << mesh->V.at(threeFiveIds.at(1)).N_Vids.size() << std::endl;
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
        tfp->Move(secondaryPath.at(i), delta, skipCheck);
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
    PARALLEL_FOR_BEGIN(0, mesh->V.size()) {
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
    // return;
    // GetSingularityGroups(Singularities, bc);

    std::cout << "END OF SIMPLIFICATION LOOP" << std::endl;
    // int it = 0;
    for (auto sid: Singularities) {
        int mSize = mesh->V.size();
        auto& s = mesh->V.at(sid);
        std::cout << "Resolving Singularity: " << s.id << std::endl;
        // if ((s.N_Vids.size() != 3 && s.N_Vids.size() != 5) || s.type == FEATURE || s.isBoundary) continue;
        int id = ResolveSingularity(sid, bc);
        while (id != -1) {
            id = ResolveSingularity(id, bc);
            // std::cout << "NEW ID: " << id << std::endl;
        }
        // if (mesh->V.size() != mSize) {
        //     Smooth();
        // }
        // it += 1;
        // if (it == iters) {
        //     break;
        // }
    }
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
    // std::recursive_mutex mtx;
    
    std::vector<SingularityLink> tempLinks;
    std::vector<bool> SingularitiesToAvoid(mesh->V.size(), false);
    PARALLEL_FOR_BEGIN(0, Singularities.size()) {
        auto& s = mesh->V.at(Singularities.at(i));
    // for (auto sid: Singularities) {
    //     auto& s = mesh->V.at(sid);
        // int valenceToCheck = s.N_Vids.size() == 3 ? 5 : 3;
        if ((s.N_Vids.size() == 3 && s.N_Vids.size() != 5) || s.type == FEATURE || s.isBoundary) continue;
        std::vector<SingularityLink> links = GetLinks(s.id, bc);
        // std::vector<SingularityLink> links = TraceSingularityLinks(s, bc);
        // {
        //     std::lock_guard<std::recursive_mutex> lock(mtx);    
        //     tempLinks.insert(tempLinks.end(), links.begin(), links.end());
        // }
        // links = SelectLinks(links);
        for (auto& l: links) {
            {
                std::lock_guard<std::recursive_mutex> lock(mtx);    
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

    int it = 0;
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
    }
    
    /*int colorValue = 0;
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
    }*/
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
    // std::recursive_mutex mtx;
    PARALLEL_FOR_BEGIN(0, verticesToCheck.size()) {
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
            newL.rots = PrototypeGetRotations(tempL.linkVids.front(), tempL.linkEids.front(), b.back());
            {
                std::lock_guard<std::recursive_mutex> lock(mtx);
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

bool SemiGlobalSimplifier::ValidatePath(std::vector<size_t> p, bool checkValences) {
    if (checkValences) {
        for (int i = 1; i < p.size()-1; i++) {
            auto& v = mesh->V.at(p.at(i));
            if (v.N_Vids.empty()) return false;
            for (auto nvid: v.N_Vids) {
                if (nvid == p.front() || nvid == p.back()) continue;
                auto& nv = mesh->V.at(nvid);
                // if (nv.N_Vids.size() != 4 || nv.type == FEATURE || nv.isBoundary) return false;
                if (nv.N_Vids.size() != 4) return false;
            }
        }
    }
    
    // for (auto id: p) if (mesh->V.at(id).N_Vids.empty()) return false;
    for (int i = 1; i < p.size(); i++) {
        auto& v1 = mesh->V.at(p.at(i-1));
        auto& v2 = mesh->V.at(p.at(i));
        // if (v1.N_Vids.empty() || v2.N_Vids.empty()) return false;
        if (std::find(v1.N_Vids.begin(), v1.N_Vids.end(), v2.id) == v1.N_Vids.end()) return false;
        if (std::find(v2.N_Vids.begin(), v2.N_Vids.end(), v1.id) == v2.N_Vids.end()) return false;
    }
    return true;
}

int SemiGlobalSimplifier::ResolveSingularity(size_t sid, BaseComplexQuad& bc) {
    auto& s = mesh->V.at(sid);
    if ((s.N_Vids.size() != 3 && s.N_Vids.size() != 5) || s.type == FEATURE || s.isBoundary) return -1;
    // std::recursive_mutex mtx;
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
        
        PARALLEL_FOR_BEGIN(0, verticesToCheck.size()) {
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
                    std::lock_guard<std::recursive_mutex> lock(mtx);
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
        tfp->Move(secondaryPath.at(i), delta, skipCheck);
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
        tfp->Move(secondaryPath.at(i), delta);
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
    // std::recursive_mutex mtx;
    std::vector<size_t> Singularities;
    
    // std::priority_queue<SingularityLink, std::vector<SingularityLink>, decltype(cmp)> q(cmp);

    /*PARALLEL_FOR_BEGIN(mesh->V.size()) {
        auto& v = mesh->V.at(i);
        if (v.type == FEATURE || v.isBoundary || !v.isSingularity || v.N_Fids.empty()) continue;
        {
            std::lock_guard<std::recursive_mutex> lock(mtx);
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
                std::lock_guard<std::recursive_mutex> lock(mtx);
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
                std::lock_guard<std::recursive_mutex> lock(mtx);
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
                        std::lock_guard<std::recursive_mutex> lock(mtx);
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
            std::lock_guard<std::recursive_mutex> lock(mtx);
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
                    std::lock_guard<std::recursive_mutex> lock(mtx);
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
                    std::lock_guard<std::recursive_mutex> lock(mtx);
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
                    std::lock_guard<std::recursive_mutex> lock(mtx);
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
                //     std::lock_guard<std::recursive_mutex> lock(mtx);
                //     q.push(l);
                // }
                // SingularityGroup sg;
                // sg.l1 = l;
                // sg.l2 = l2;
                // sg.rank = l.a + l.b + l2.a + l2.b;
                // {
                //     std::lock_guard<std::recursive_mutex> lock(mtx);
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
                    std::lock_guard<std::recursive_mutex> lock(mtx);
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
    // std::recursive_mutex mtx;
    PARALLEL_FOR_BEGIN(0, mesh->V.size()) {
        auto& v = mesh->V.at(i);
    // for (auto& v: mesh->V) {
        if (!v.isSingularity || v.type == FEATURE || v.isBoundary) continue;
        std::vector<SingularityLink> links = TraceSingularityLinks(v, bc);
        for (auto& l: links) {
            if (mesh->V.at(l.frontId).N_Vids.size() == 4 || mesh->V.at(l.backId).N_Vids.size() == 4) continue;
            {
                std::lock_guard<std::recursive_mutex> lock(mtx);
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
    PARALLEL_FOR_BEGIN(0, n) {
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
                    std::lock_guard<std::recursive_mutex> lock(mtx);
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
                    std::lock_guard<std::recursive_mutex> lock(mtx);
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
    // std::recursive_mutex mtx;
    auto cmp = [](SingularityGroup left, SingularityGroup right) {return left.rank > right.rank;};
    std::priority_queue<SingularityGroup, std::vector<SingularityGroup>, decltype(cmp)> q(cmp);
    
    PARALLEL_FOR_BEGIN(0, mesh->V.size()) { 
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
            std::lock_guard<std::recursive_mutex> lock(mtx);
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

void SemiGlobalSimplifier::MoveSingularities(size_t& toMoveId, size_t& sourceId, size_t& secondaryId, size_t& sourceDir, size_t& secondaryDir, std::vector<size_t> secondaryPath) {
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
        tfp->Move(secondaryPath.at(i), delta);
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
                // std::cout << "Collapsing 3-singularity for 3-5 pair generation" << std::endl;
                res.push_back(f.Vids.at((idx+2)%f.Vids.size()));
                std::shared_ptr<SimplificationOperation> dc = std::make_shared<DiagonalCollapse>(*mesh, *mu, *smoother, f.id, (idx+2)%f.Vids.size(), idx);
                dc->PerformOperation();
                // std::cout << "After collapsing" << std::endl;
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
            if (l.linkVids.front() == l.linkVids.back() || mesh->V.at(l.linkVids.back()).type == FEATURE || mesh->V.at(l.linkVids.back()).isBoundary || mesh->V.at(l.linkVids.back()).N_Vids.size() > 6) continue;
            if (doesCrossBoundary(l.linkVids, true)) continue;
        
            l.frontId = l.linkVids.front();
            l.backId = l.linkVids.back();
            links.push_back(l);
        }
    }
    return links;
}

void SemiGlobalSimplifier::TraceSingularityLinks(size_t vid, std::vector<SingularityLink>& links, bool traceDiagonals) {
    auto& v = mesh->V.at(vid);
    std::vector<SingularityLink> newLinks;
    for (auto eid: v.N_Eids) {
        TraceLink(v, mesh->E.at(eid), newLinks);
    }
    std::vector<SingularityLink> tmpLinks;
    for (auto& l: newLinks) {
        PARALLEL_FOR_BEGIN(1, l.linkVids.size()-1) {
        // for (int i = 1; i < l.linkVids.size()-1; i++) {
            auto& nv = mesh->V.at(l.linkVids.at(i));
            for (auto neid: nv.N_Eids) {
                auto& ne = mesh->E.at(neid);
                if (ne.Vids.at(0) == l.linkVids.at(i-1) || ne.Vids.at(1) == l.linkVids.at(i-1)
                || ne.Vids.at(0) == l.linkVids.at(i+1) || ne.Vids.at(1) == l.linkVids.at(i+1)) continue;
                if ((mesh->V.at(ne.Vids.at(0)).type == FEATURE || mesh->V.at(ne.Vids.at(0)).isBoundary)
                && (mesh->V.at(ne.Vids.at(1)).type == FEATURE || mesh->V.at(ne.Vids.at(1)).isBoundary)) continue;
                TraceLink(nv, ne, tmpLinks, std::vector<size_t>(l.linkVids.begin(), l.linkVids.begin()+i), std::vector<size_t>(l.linkEids.begin(), l.linkEids.begin()+i));
            }
        // }
        } PARALLEL_FOR_END();
    }
    links.insert(links.end(), newLinks.begin(), newLinks.end());
    links.insert(links.end(), tmpLinks.begin(), tmpLinks.end());
    if (!traceDiagonals) return;
    std::vector<SingularityLink> newLinksD;
    for (auto fid: v.N_Fids) {
        TraceLink(v, mesh->F.at(fid), newLinksD);
    }
    std::vector<SingularityLink> tmpLinksD;
    for (auto& l: newLinksD) {
        PARALLEL_FOR_BEGIN(1, l.linkVids.size()-1) {
        // for (int i = 1; i < l.linkVids.size()-1; i++) {
            auto& nv = mesh->V.at(l.linkVids.at(i));
            for (auto nfid: nv.N_Fids) {
                auto& nf = mesh->F.at(nfid);
                size_t nextV = GetFaceV(nv.id, nf.id, 2);
                if (nextV == l.linkVids.at(i-1) || nextV == l.linkVids.at(i+1)) continue;
                if ((mesh->V.at(nextV).type == FEATURE || mesh->V.at(nextV).isBoundary)) continue;
                TraceLink(nv, nf, tmpLinksD, std::vector<size_t>(l.linkVids.begin(), l.linkVids.begin()+i));
            }
        // }
        } PARALLEL_FOR_END();
    }
    links.insert(links.end(), newLinksD.begin(), newLinksD.end());
    links.insert(links.end(), tmpLinksD.begin(), tmpLinksD.end());
}

void SemiGlobalSimplifier::TraceLink(const Vertex& v, const Edge& edge, std::vector<SingularityLink>& links, std::vector<size_t> vids, std::vector<size_t> eids) {
    SingularityLink l;
    TraceAlongEdge(v, edge, l.linkVids, l.linkEids);
    if (l.linkVids.empty()) return;
    if (mesh->V.at(l.linkVids.back()).N_Vids.size() > 6) return;
    if (vids.empty()) {
        l.a = l.linkEids.size();
        l.rank = l.linkEids.size();
    } else {
        l.a = eids.size();
        l.b = l.linkEids.size();
        l.linkVids.insert(l.linkVids.begin(), vids.begin(), vids.end());
        l.linkEids.insert(l.linkEids.begin(), eids.begin(), eids.end());
    }
    l.rank = l.a + l.b;
    l.frontId = l.linkVids.front();
    l.backId = l.linkVids.back();
    l.delta = &delta;
    {
        std::lock_guard<std::recursive_mutex> lock(mtx);
        links.push_back(l);
    }
}

void SemiGlobalSimplifier::TraceLink(const Vertex& v, const Face& face, std::vector<SingularityLink>& links, std::vector<size_t> vids) {
    SingularityLink l;
    TraceAlongDiagonal(v, face, l.linkVids);
    if (l.linkVids.empty()) return;
    if (mesh->V.at(l.linkVids.back()).N_Vids.size() > 6) return;
    if (vids.empty()) {
        l.a = l.linkVids.size() - 1;
    } else {
        l.a = vids.size();
        l.b = l.linkVids.size()-1;
        l.linkVids.insert(l.linkVids.begin(), vids.begin(), vids.end());
    }
    l.rank = l.a + l.b;
    l.frontId = l.linkVids.front();
    l.backId = l.linkVids.back();
    l.diagonal = true;
    l.delta = &delta;
    {
        std::lock_guard<std::recursive_mutex> lock(mtx);
        links.push_back(l);
    }
}

void SemiGlobalSimplifier::TraceAlongEdge(const Vertex& start_vertex, const Edge& start_edge, std::vector<size_t>& vids_link, std::vector<size_t>& eids_link) {
    Vertex* currentVertex = (Vertex*) &start_vertex;
    Edge* currentEdge = (Edge*) &start_edge;
    vids_link.push_back(currentVertex->id);
    while (currentEdge != NULL) {
        eids_link.push_back(currentEdge->id);
        currentVertex = (Vertex*)&mesh->V.at(currentEdge->Vids[0] == currentVertex->id ? currentEdge->Vids[1] : currentEdge->Vids[0]);
        vids_link.push_back(currentVertex->id);
        if (currentVertex->N_Vids.size() != 4 || currentVertex->id == start_vertex.id) break;
        std::vector<size_t> edgesToAvoid;
        mu->AddContents(edgesToAvoid, mu->GetIntersection(currentVertex->N_Eids, mesh->F.at(currentEdge->N_Fids.at(0)).Eids));
        mu->AddContents(edgesToAvoid, mu->GetIntersection(currentVertex->N_Eids, mesh->F.at(currentEdge->N_Fids.at(1)).Eids));
        currentEdge = (Edge*)&mesh->E.at(mu->GetDifference(currentVertex->N_Eids, edgesToAvoid).at(0));
    }
}

void SemiGlobalSimplifier::TraceAlongEdge(const Vertex& start_vertex, const Edge& start_edge, std::vector<size_t>& vids_link, bool checkBoundary) {
    Vertex* currentVertex = (Vertex*) &start_vertex;
    Edge* currentEdge = (Edge*) &start_edge;
    vids_link.push_back(currentVertex->id);
    while (currentEdge != NULL) {
        currentVertex = (Vertex*)&mesh->V.at(currentEdge->Vids[0] == currentVertex->id ? currentEdge->Vids[1] : currentEdge->Vids[0]);
        vids_link.push_back(currentVertex->id);
        if (currentVertex->N_Vids.size() != 4 || (checkBoundary && (currentVertex->isBoundary || currentVertex->type == FEATURE)) || currentVertex->id == start_vertex.id) break;
        std::vector<size_t> edgesToAvoid;
        mu->AddContents(edgesToAvoid, mu->GetIntersection(currentVertex->N_Eids, mesh->F.at(currentEdge->N_Fids.at(0)).Eids));
        mu->AddContents(edgesToAvoid, mu->GetIntersection(currentVertex->N_Eids, mesh->F.at(currentEdge->N_Fids.at(1)).Eids));
        currentEdge = (Edge*)&mesh->E.at(mu->GetDifference(currentVertex->N_Eids, edgesToAvoid).at(0));
    }
}

void SemiGlobalSimplifier::TraceAlongEdge(const Vertex& start_vertex, const Edge& start_edge, std::vector<size_t>& vids_link, int b1, int b2, int rotOffset, bool checkBoundary) {
    Vertex* currentVertex = (Vertex*) &start_vertex;
    Edge* currentEdge = (Edge*) &start_edge;
    Edge* prevEdge = (Edge*) &start_edge;
    vids_link.push_back(currentVertex->id);
    while (b1-- > 0) {
        currentVertex = (Vertex*)&mesh->V.at(currentEdge->Vids[0] == currentVertex->id ? currentEdge->Vids[1] : currentEdge->Vids[0]);
        vids_link.push_back(currentVertex->id);
        if (currentVertex->N_Vids.size() != 4 || (checkBoundary && (currentVertex->isBoundary || currentVertex->type == FEATURE)) || currentVertex->id == start_vertex.id) break;
        std::vector<size_t> edgesToAvoid;
        mu->AddContents(edgesToAvoid, mu->GetIntersection(currentVertex->N_Eids, mesh->F.at(currentEdge->N_Fids.at(0)).Eids));
        mu->AddContents(edgesToAvoid, mu->GetIntersection(currentVertex->N_Eids, mesh->F.at(currentEdge->N_Fids.at(1)).Eids));
        prevEdge = currentEdge;
        currentEdge = (Edge*)&mesh->E.at(mu->GetDifference(currentVertex->N_Eids, edgesToAvoid).at(0));
    }
    currentEdge = (Edge*)&mesh->E.at(GetCCedgeAt(currentVertex->id, currentEdge->id, rotOffset));
    while (b2-- > 0) {
        currentVertex = (Vertex*)&mesh->V.at(currentEdge->Vids[0] == currentVertex->id ? currentEdge->Vids[1] : currentEdge->Vids[0]);
        vids_link.push_back(currentVertex->id);
        if (currentVertex->N_Vids.size() != 4 || currentVertex->id == start_vertex.id) break;
        std::vector<size_t> edgesToAvoid;
        mu->AddContents(edgesToAvoid, mu->GetIntersection(currentVertex->N_Eids, mesh->F.at(currentEdge->N_Fids.at(0)).Eids));
        mu->AddContents(edgesToAvoid, mu->GetIntersection(currentVertex->N_Eids, mesh->F.at(currentEdge->N_Fids.at(1)).Eids));
        currentEdge = (Edge*)&mesh->E.at(mu->GetDifference(currentVertex->N_Eids, edgesToAvoid).at(0));
    }
}

void SemiGlobalSimplifier::TraceAlongDiagonal(const Vertex& start_vertex, const Face& start_face, std::vector<size_t>& vids_link) {
    Vertex* currentVertex = (Vertex*)&start_vertex;
    Face* currentFace = (Face*)&start_face;
    vids_link.push_back(currentVertex->id);
    while (currentVertex != NULL) {
        currentVertex = (Vertex*)&mesh->V.at(GetFaceV(currentVertex->id, currentFace->id, 2));
        vids_link.push_back(currentVertex->id);
        if (currentVertex->N_Vids.size() != 4 || currentVertex->id == start_vertex.id) break;
        std::vector<size_t> diffFaces;
        for (auto eid: currentFace->Eids) mu->AddContents(diffFaces, mesh->E.at(eid).N_Fids);
        diffFaces = mu->GetDifference(currentVertex->N_Fids, diffFaces);
        if (diffFaces.empty()) return;
        currentFace = (Face*)&mesh->F.at(diffFaces.at(0)); 
    }
}

void SemiGlobalSimplifier::GetDiagonalPath(SingularityLink& l) {
    l.rots = GetEdgeRots(l.linkEids.at(l.a-1), l.linkEids.at(l.a));
    size_t startIdx = 0;
    size_t endIdx = 2*l.a;
    int steps = l.a;
    if (l.a > l.b) {
        startIdx = l.a - l.b;
        endIdx = l.linkVids.size()-1;
        steps = l.b;
    }
    auto& v = mesh->V.at(l.linkVids.at(startIdx));
    std::vector<size_t> temp;
    for (auto fid: v.N_Fids) {
        size_t currentVid = l.linkVids.at(startIdx);
        size_t currentFid = fid;
        for (int i = 0; i < steps; i++) {
            auto& f = mesh->F.at(currentFid);
            currentVid = GetFaceV(currentVid, f.id, 2);
            temp.push_back(currentVid);
            std::vector<size_t> diffFaces;
            for (auto eid: f.Eids) mu->AddContents(diffFaces, mesh->E.at(eid).N_Fids);
            if (mesh->V.at(currentVid).N_Vids.size() != 4) break;
            currentFid = mu->GetDifference(mesh->V.at(currentVid).N_Fids, diffFaces).at(0);
        }
        if (currentVid == l.linkVids.at(endIdx)) {   
            std::vector<size_t> newVids;
            if (l.a > l.b) {
                newVids.insert(newVids.end(), l.linkVids.begin(), l.linkVids.begin()+startIdx+1);
                newVids.insert(newVids.end(), temp.begin(), temp.end());
                l.a -= steps;
                l.b = steps;
            } else {
                newVids = {l.linkVids.at(startIdx)};
                newVids.insert(newVids.end(), temp.begin(), temp.end());
                newVids.insert(newVids.end(), l.linkVids.begin()+endIdx, l.linkVids.end());
                l.a = steps;
                l.b -= steps;    
            }
            l.linkVids = newVids;
            l.rank = l.a + l.b;
            break;
        }
        temp.clear();
    }
    if (temp.empty()) {
        l.linkVids.clear();
    }
}

int SemiGlobalSimplifier::GetEdgeRots(size_t eid1, size_t eid2) {
    auto& e1 = mesh->E.at(eid1);
    auto& e2 = mesh->E.at(eid2);
    auto& v = mesh->V.at(mu->GetIntersection(e1.Vids, e2.Vids).at(0));

    int nRots = 0;
    size_t eid = e1.id;
    for (int i = 0; i < v.N_Eids.size(); i++) {
        if (eid == e2.id) break;
        auto& e = mesh->E.at(eid);
        size_t prev = e.Vids.at(0) == v.id ? e.Vids.at(1) : e.Vids.at(0);
        for (auto fid: e.N_Fids) {
            auto& f = mesh->F.at(fid);
            if (GetFaceV(v.id, f.id, 1) == prev) {
                eid = mu->GetDifference(mu->GetIntersection(v.N_Eids, f.Eids), std::vector<size_t>{e.id}).at(0);
                break;
            }
        }
        nRots += 1;
    }
    return nRots;
}

int SemiGlobalSimplifier::GetFaceRots(size_t fid1, size_t fid2, size_t vid) {
    auto& f1 = mesh->E.at(fid1);
    auto& f2 = mesh->E.at(fid2);
    auto& v = mesh->V.at(vid);

    int nRots = 0;
    size_t fid = f1.id;
    for (int i = 0; i < v.N_Fids.size(); i++) {
        if (fid == f2.id) break;
        auto& f = mesh->F.at(fid);
        size_t prev = GetFaceV(v.id, f.id, 3);
        for (auto vfid: v.N_Fids) {
            auto& vf = mesh->F.at(vfid);
            if (GetFaceV(v.id, vf.id, 1) == prev) {
                fid = vfid;
                break;
            }
        }
        nRots += 1;
    }
    return nRots;
}

void SemiGlobalSimplifier::PerformGlobalOperations() {
    CheckValidity();
}

size_t SemiGlobalSimplifier::GetFaceID(size_t vid, size_t exclude_vid) {
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
    if (in.empty()) return false;
    if (mu->IsSharpFeature(in.back())) return true;
    for (int id = 0; id < in.size(); id++) {
        // if (isVertex) {
            auto& v = mesh->V.at(in.at(id));
            if (v.type == FEATURE || v.isBoundary) {
                return true;
            }
        // } else {
        //     auto& e = mesh->E.at(in.at(id));
        //     if (mesh->V.at(e.Vids.at(0)).type == FEATURE || mesh->V.at(e.Vids.at(0)).isBoundary || mesh->V.at(e.Vids.at(1)).type == FEATURE || mesh->V.at(e.Vids.at(1)).isBoundary) {
        //         return true;
        //     }
        // }
    }
    return false;
}

size_t SemiGlobalSimplifier::GetDiagonalV(size_t vid, size_t fid) {
	auto& f = mesh->F.at(fid);
    int index = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), vid));
	return f.Vids.at((index+2)%f.Vids.size());
}

size_t SemiGlobalSimplifier::GetFaceV(size_t vid, size_t fid, int offset) {
	auto& f = mesh->F.at(fid);
    int index = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), vid));
	return f.Vids.at((index+offset)%f.Vids.size());
}

bool SemiGlobalSimplifier::Contains(std::vector<size_t> v, size_t val) {
    return (std::find(v.begin(), v.end(), val) != v.end());
}

bool SemiGlobalSimplifier::Contains(std::vector<size_t> v, std::vector<size_t> v2) {
    return !mu->GetIntersection(v, v2).empty();
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

bool SemiGlobalSimplifier::RemoveDoublets() {
    bool res = false;
    for (auto& v: mesh->V) {
        if (!v.isBoundary && v.type != FEATURE && v.N_Vids.size() == 2) res = true;
        VertexSplit s(*mesh, *mu, *smoother, v.id);
        s.FixDoublet(v.id);
    }
    return res;
}

bool SemiGlobalSimplifier::FixValences() {
    bool res = false;

    res = RemoveDoublets();
    // ResolveSingularityPairs();
    for (auto& v: mesh->V) {
        if ((v.type == FEATURE || v.isBoundary) && !v.N_Fids.empty()) mesh->SetIdealValence(v.id);
    }
    res = FixBoundary();
    res = ResolveHighValences();
    return res;
}

void SemiGlobalSimplifier::PrototypeK() {
    std::priority_queue<Separatrix, std::vector<Separatrix>, SeparatrixComparator> q;
    std::vector<std::vector<size_t>> vec;
    std::vector<std::vector<size_t>> vec2;
    int it = 0;
    for (auto& v: mesh->V) {
        if ((v.N_Vids.size() != 3 && v.N_Vids.size() != 5)|| v.isBoundary || v.type == FEATURE) continue;
        Separatrix s;
        for (auto eid: v.N_Eids) {
            PrototypeTraceSeparatrix(v.id, eid, s);
        }
        if (!s.empty) {
            q.push(s);
            // it += 1;
            // if (it == iters) {
                // vec.push_back(s.vidsA);
                // vec.push_back(s.vidsB);
                // std::cout << "separatrix boundary: " << s.vidsA.size() << std::endl;
            // }
        }
        // if (it == iters){
            // break;
        // }
    }
    std::cout << "Extracted " << q.size() << " separatrices" << std::endl;
    while (!q.empty()) {
        auto& s = (Separatrix&) q.top();
        if (s.b1 == 1 && s.b2 == 0) {
        //     it += 1;
        //     // bool includeIters = it == 7 ? true : false;
            vec.push_back(s.vidsA);
            bool includeIters = false;
            PrototypeTransportSeparatrix(s, includeIters, vec2);
        }
        if (s.b1 == 1 && s.b2 == 1) {
        //     it += 1;
        //     // bool includeIters = it == 7 ? true : false;
            vec.push_back(s.vidsA);
            bool includeIters = false;
            PrototypeTransportSeparatrix(s, includeIters, vec2);
        }
        if (s.b1 > 1 && s.b2 == 0) {
            // it += 1;
            // bool includeIters = it == 7 ? true : false;
            vec.push_back(s.vidsA);
            bool includeIters = false;
            PrototypeTransportSeparatrix(s, includeIters, vec2);
        }
        q.pop();
        // std::cout << "it: " << it << " iters: " << iters << std::endl;
        if (it == iters) {
            // break;
        }
    }
    // PrototypeSaveSeparatrices(vec, "primary");
    // PrototypeSaveSeparatrices(vec2, "secondary");
}

void SemiGlobalSimplifier::PrototypeTraceSeparatrix(size_t vid, size_t eid, Separatrix& s, bool checkValence) {
    auto& v = mesh->V.at(vid);
    int valenceToCheck = v.N_Vids.size() == 3 ? 5 : 3;
    std::vector<size_t> vids;
    TraceAlongEdge(mesh->V.at(vid), mesh->E.at(eid), vids);
    if (vids.empty()) return;
    if (checkValence && mesh->V.at(vids.back()).N_Vids.size() != valenceToCheck) return;
    if (!checkValence && (mesh->V.at(vids.back()).N_Vids.size() != 3 && mesh->V.at(vids.back()).N_Vids.size() != 5)) return;
    if (!ValidatePath(vids)) return;
    s.vidsA = vids;
    s.frontId = s.vidsA.front();
    s.backId = s.vidsA.back();
    s.b1 = vids.size()-1;
    s.b2 = 0;
    s.empty = false;
    s.minSize = s.b1+s.b2;
    // PrototypeSetDirMap(s);

    int n = vids.size()-1;
    for (int i = 1; i < n; i++) {
        int b1 = i;
        int b2 = 0;
        auto& nv = mesh->V.at(vids.at(i));
        for (auto eid: nv.N_Eids) {
            auto& e = mesh->E.at(eid);
            if (Contains(e.Vids, vids.at(i+1)) || Contains(e.Vids, vids.at(i-1))) continue;
            std::vector<size_t> secVids;
            TraceAlongEdge(v, e, secVids);
            if (secVids.empty()) continue;
            if (checkValence && mesh->V.at(secVids.back()).N_Vids.size() != valenceToCheck) return;
            if (!checkValence && (mesh->V.at(secVids.back()).N_Vids.size() != 3 && mesh->V.at(secVids.back()).N_Vids.size() != 5)) return;
            if (!ValidatePath(secVids)) return;
    
            b2 = secVids.size()-1;
            size_t rotEdge = GetEdgeId(secVids.at(0), secVids.at(1));
            s.rots = GetEdgeRots(rotEdge, eid);
            secVids.insert(secVids.begin(), vids.begin(), vids.begin()+i);
            
            s.vidsA = secVids;
            s.frontId = secVids.front();
            s.backId = secVids.back();
            s.b1 = b1;
            s.b2 = b2;
            s.empty = false;
            s.minSize = s.b1+s.b2;
            // PrototypeSetDirMap(s);
        }
    }
}

void SemiGlobalSimplifier::PrototypeTransportSeparatrix(Separatrix& s, bool includeIters, std::vector<std::vector<size_t>>& vec) {
    std::cout << "vidsA: " << s.vidsA.size() << " vidsB: " << s.vidsB.size() << std::endl; 
    std::cout << "three: " << s.threeId << " " << mesh->V.at(s.threeId).N_Vids.size() << " five: " << s.fiveId << " " << mesh->V.at(s.fiveId).N_Vids.size() << std::endl;
    std::cout << "b1: " << s.b1 << " b2: " << s.b2 << std::endl;
    if (mesh->V.at(s.threeId).N_Vids.size() != 3 || mesh->V.at(s.fiveId).N_Vids.size() != 5) return;
    if (!ValidatePath(s.vidsA)) return;
    std::cout << "Paths validated" << std::endl;
    if (s.b1 == 1 && s.b2 == 0) {
        PrototypeTransportDirectPair(s.threeId, s.fiveId, includeIters);
    } else if (s.b1 == 1 && s.b2 == 1) {
        PrototypeTransportDiagonalPair(s.threeId, s.fiveId, includeIters);
    } 
    // else {
        // PrototypeTransportLink(s);
    // } 
    else if (s.b1 > 1 && s.b2 == 0) {
        PrototypeTransportDirectLink(s, vec);
    }
}

void SemiGlobalSimplifier::PrototypeTransportDiagonalPair(size_t threeId, size_t fiveId, bool includeIters, std::vector<size_t> verticesToAvoid) {
    Separatrix s;
    auto& three = mesh->V.at(threeId);
    auto& five = mesh->V.at(fiveId);

    mu->AddContents(verticesToAvoid, std::vector<size_t>{threeId, fiveId});
    auto& f = mesh->F.at(mu->GetIntersection(three.N_Fids, five.N_Fids).at(0));
    // std::cout << "f Vids: ";
    // for (auto id: f.Vids) std::cout << id << " ";
    // std::cout << std::endl;
    // std::cout << "f Eids: ";
    // for (auto id: f.Eids) std::cout << id << " ";
    // std::cout << std::endl;
    // std::cout << "five N_Eids: ";
    // for (auto id: five.N_Eids) std::cout << id << " ";
    // std::cout << std::endl;
    std::vector<size_t> nRotEdges = mu->GetIntersection(f.Eids, five.N_Eids);
    std::vector<size_t> rotEdges;
    for (auto eid: nRotEdges) {
        auto& e = mesh->E.at(eid);
        size_t efid = e.N_Fids.at(0) == f.id ? e.N_Fids.at(1) : e.N_Fids.at(0);
        mu->AddContents(rotEdges, mu->GetDifference(mu->GetIntersection(five.N_Eids, mesh->F.at(efid).Eids), std::vector<size_t>{eid}));
    }
    mu->AddContents(nRotEdges, mu->GetDifference(five.N_Eids, mu->GetUnion(nRotEdges, rotEdges)));
    mu->AddContents(nRotEdges, mu->GetDifference(three.N_Eids, mu->GetIntersection(f.Eids, three.N_Eids)));
    mu->AddContents(rotEdges, mu->GetIntersection(f.Eids, three.N_Eids));
    // std::cout << "nRotEdges: " << nRotEdges.size() << " rotEdges: " << rotEdges.size() << std::endl;
    for (auto eid: rotEdges) {
        auto& e = mesh->E.at(eid);
        size_t vid = e.Vids.at(0) == threeId || e.Vids.at(0) == fiveId ? e.Vids.at(0) : e.Vids.at(1);
        size_t vertexToAvoid = vid == threeId ? fiveId : threeId;
        // std::cout << "Tracing rotEdge" << std::endl;
        PrototypeTraceDiagSeparatrix(vid, eid, s, true, verticesToAvoid);
    }
    for (auto eid: nRotEdges) {
        auto& e = mesh->E.at(eid);
        size_t vid = e.Vids.at(0) == threeId || e.Vids.at(0) == fiveId ? e.Vids.at(0) : e.Vids.at(1);
        size_t vertexToAvoid = vid == threeId ? fiveId : threeId;
        // std::cout << "Tracing nRotEdge" << std::endl;
        PrototypeTraceDiagSeparatrix(vid, eid, s, false, verticesToAvoid);
    }
    if (!s.empty) {
        // std::cout << "separatrix: " << s.vidsA.size() << " front: " << mesh->V.at(s.vidsA.front()).N_Vids.size() << " back: " << mesh->V.at(s.vidsA.back()).N_Vids.size() << std::endl;
        // std::cout << "separatrix b1: " << s.b1 << " b2: " << s.b2 << std::endl;
        // PrototypeSaveSeparatrices(std::vector<std::vector<size_t>>{s.vidsA});
        std::shared_ptr<DiagonalThreeFivePair> tfp = std::make_shared<DiagonalThreeFivePair>(*mesh, *mu, *smoother, threeId, fiveId);
        // int n = iters < s.vidsA.size() ? iters : s.vidsA.size();
        int n = s.vidsA.size();
        if (includeIters) n = iters < s.vidsA.size() ? iters : s.vidsA.size();
        for (int i = 1; i < n; i++) {
            tfp->Move(s.vidsA.at(i), delta);
        }
    }
}

void SemiGlobalSimplifier::PrototypeTraceDiagSeparatrix(size_t vid, size_t eid, Separatrix& s, bool checkBoundary, std::vector<size_t> verticesToAvoid) {
    auto& v = mesh->V.at(vid);
    std::vector<size_t> vids;
    TraceAlongEdge(mesh->V.at(vid), mesh->E.at(eid), vids, checkBoundary);
    if (!vids.empty()) {
        int b1 = checkBoundary ? 0 : vids.size()-1;
        int b2 = 0;
        if (!PrototypeCheckBoundarySingularity(vids.back()) && ValidatePath(vids) && !Contains(verticesToAvoid, vids.back()) && (mesh->V.at(vids.back()).N_Vids.size() == 5 || mesh->V.at(vids.back()).N_Vids.size() == 3) && (s.minSize == -1 || s.minSize > b1+b2)) {
            // std::cout << "checkBoundary: " << checkBoundary << std::endl;
            s.frontId = vids.front();
            s.backId = vids.back();
            s.threeId = s.frontId;
            s.fiveId = s.backId;
            s.b1 = b1;
            s.b2 = b2;
            s.vidsA = vids;
            s.vidsB.clear();
            s.empty = false;
            s.minSize = s.b1+s.b2;
        }
        int n = vids.size()-1;
        // std::cout << "n: " << n << std::endl;
        for (int i = 1; i < n; i++) {
            b1 = checkBoundary ? 0 : i;
            auto& v = mesh->V.at(vids.at(i));
            for (auto eid: v.N_Eids) {
                auto& e = mesh->E.at(eid);
                if (Contains(e.Vids, vids.at(i+1)) || Contains(e.Vids, vids.at(i-1))) continue;
                // std::cout << "GetEdgeId" << std::endl;
                // size_t edgeId = GetEdgeId(vids.at(i), vids.at(i-1));
                // std::cout << "GetCCedgeAt" << std::endl;
                // edgeId = GetCCedgeAt(vids.at(i), edgeId, 1);
                std::vector<size_t> secVids;
                // std::cout << "Tracing secondary link" << std::endl;
                TraceAlongEdge(v, e, secVids, !checkBoundary);
                b2 = checkBoundary ? secVids.size()-1 : 0;
                secVids.insert(secVids.begin(), vids.begin(), vids.begin()+i);
                if (!PrototypeCheckBoundarySingularity(secVids.back()) && ValidatePath(secVids) && !Contains(verticesToAvoid, secVids.back()) && (mesh->V.at(secVids.back()).N_Vids.size() == 5 || mesh->V.at(secVids.back()).N_Vids.size() == 3) && (s.minSize == -1 || s.minSize > b1+b2)) {
                    // std::cout << "checkBoundary: " << checkBoundary << std::endl;
                //     size_t edgeId = GetEdgeId(secVids.at(secVids.size()-1), secVids.at(secVids.size()-2));    
                //     edgeId = GetCCedgeAt(secVids.back(), edgeId, mesh->V.at(secVids.back()).N_Vids.size()-1);
                //     std::vector<size_t> secVidsB;
                    // std::cout << "b1: " << b1 << " b2: " << b2 << std::endl;
                    // TraceAlongEdge(mesh->V.at(secVids.back()), mesh->E.at(edgeId), secVidsB, b1, b2);
                    // std::cout << "secVids: " << secVids.size() << " " << " front: " << secVids.front() << std::endl;
                    // std::cout << "secVidsB: " << secVidsB.size() << " " << "back: " << secVidsB.back() << std::endl;
                    // if (secVids.front() == secVidsB.back()) {
                        s.frontId = secVids.front();
                        s.backId = secVids.back();
                        s.threeId = s.frontId;
                        s.fiveId = s.backId;
                        s.b1 = b1;
                        s.b2 = b2;
                        s.vidsA = secVids;
                        // s.vidsB = secVidsB;
                        s.empty = false;
                        s.minSize = s.b1+s.b2;
                    // }
                }
            }
        }
    }
}

void SemiGlobalSimplifier::PrototypeTransportDirectPair(size_t threeId, size_t fiveId, bool includeIters, std::vector<size_t> verticesToAvoid) {
    Separatrix s;
    auto& three = mesh->V.at(threeId);
    auto& five = mesh->V.at(fiveId);
    mu->AddContents(verticesToAvoid, std::vector<size_t>{threeId, fiveId});

    size_t fiveEdge = mu->GetIntersection(three.N_Eids, five.N_Eids).at(0);
    std::vector<size_t> threeEdges = mu->GetDifference(three.N_Eids, std::vector<size_t>{fiveEdge});
    std::vector<size_t> rotEdges;
    for (auto fid: mesh->E.at(fiveEdge).N_Fids) {
        mu->AddContents(rotEdges, mu->GetDifference(mu->GetIntersection(five.N_Eids, mesh->F.at(fid).Eids), std::vector<size_t>{fiveEdge}));
    }
    std::vector<size_t> nRotEdges = mu->GetDifference(five.N_Eids, mu->GetUnion(rotEdges, std::vector<size_t>{fiveEdge}));
    int rotOffset = 1;
    for (auto eid: threeEdges) {
        if (GetEdgeRots(eid, fiveEdge) > 1) {
            rotOffset = 3;
        } else {
            rotOffset = 1;
        }
        auto& e = mesh->E.at(eid);
        size_t vid = e.Vids.at(0) == threeId || e.Vids.at(0) == fiveId ? e.Vids.at(0) : e.Vids.at(1);
        size_t vertexToAvoid = vid == threeId ? fiveId : threeId;
        PrototypeTraceDirectSeparatrix(three.id, eid, s, verticesToAvoid, rotOffset, false);
    }
    for (auto eid: rotEdges) {
        if (GetEdgeRots(eid, fiveEdge) > 1) {
            rotOffset = 3;
        } else {
            rotOffset = 1;
        }
        auto& e = mesh->E.at(eid);
        size_t vid = e.Vids.at(0) == threeId || e.Vids.at(0) == fiveId ? e.Vids.at(0) : e.Vids.at(1);
        size_t vertexToAvoid = vid == threeId ? fiveId : threeId;
        PrototypeTraceDirectSeparatrix(five.id, eid, s, verticesToAvoid, rotOffset, true);
    }
    for (auto eid: nRotEdges) {
        if (GetEdgeRots(eid, fiveEdge) > 2) {
            rotOffset = 3;
        } else {
            rotOffset = 1;
        }
        auto& e = mesh->E.at(eid);
        size_t vid = e.Vids.at(0) == threeId || e.Vids.at(0) == fiveId ? e.Vids.at(0) : e.Vids.at(1);
        size_t vertexToAvoid = vid == threeId ? fiveId : threeId;
        PrototypeTraceDirectSeparatrix(five.id, eid, s, verticesToAvoid, rotOffset, false);
    }

    // std::cout << "f Vids: ";
    // for (auto id: f.Vids) std::cout << id << " ";
    // std::cout << std::endl;
    // std::cout << "f Eids: ";
    // for (auto id: f.Eids) std::cout << id << " ";
    // std::cout << std::endl;
    // std::cout << "five N_Eids: ";
    // for (auto id: five.N_Eids) std::cout << id << " ";
    // std::cout << std::endl;
    
    if (!s.empty) {
        // std::cout << "separatrix: " << s.vidsA.size() << " front: " << mesh->V.at(s.vidsA.front()).N_Vids.size() << " back: " << mesh->V.at(s.vidsA.back()).N_Vids.size() << std::endl;
        // std::cout << "separatrix b1: " << s.b1 << " b2: " << s.b2 << std::endl;
        // PrototypeSaveSeparatrices(std::vector<std::vector<size_t>>{s.vidsA});
        std::shared_ptr<ThreeFivePair> tfp = std::make_shared<ThreeFivePair>(*mesh, *mu, *smoother, threeId, fiveId);
        // int n = iters < s.vidsA.size() ? iters : s.vidsA.size();
        int n = s.vidsA.size();
        // if (includeIters) n = iters < s.vidsA.size() ? iters : s.vidsA.size();
        for (int i = 1; i < n; i++) {
            tfp->Move(s.vidsA.at(i), delta);
        }
    }
}

void SemiGlobalSimplifier::PrototypeTraceDirectSeparatrix(size_t vid, size_t eid, Separatrix& s, std::vector<size_t> verticesToAvoid, int rotOffset, bool isRotEdge) {
    auto& v = mesh->V.at(vid);
    std::vector<size_t> vids;
    TraceAlongEdge(mesh->V.at(vid), mesh->E.at(eid), vids);
    if (!vids.empty()) {
        int b1 = isRotEdge && vids.size() > 2 ? vids.size()-2 : vids.size()-1;
        int b2 = 0;
        if (!PrototypeCheckBoundarySingularity(vids.back()) && ValidatePath(vids) && !Contains(verticesToAvoid, vids.back()) && (mesh->V.at(vids.back()).N_Vids.size() == 5 || mesh->V.at(vids.back()).N_Vids.size() == 3) && (s.minSize == -1 || s.minSize > b1+b2)) {
            // std::cout << "checkBoundary: " << checkBoundary << std::endl;
            s.frontId = vids.front();
            s.backId = vids.back();
            s.threeId = s.frontId;
            s.fiveId = s.backId;
            s.b1 = b1;
            s.b2 = b2;
            s.vidsA = vids;
            s.vidsB.clear();
            s.empty = false;
            s.minSize = s.b1+s.b2;
        }
        int n = vids.size()-1;
        // std::cout << "n: " << n << std::endl;
        for (int i = 1; i < n; i++) {
            b1 = i;
            auto& v = mesh->V.at(vids.at(i));
            size_t eid = GetEdgeId(vids.at(i), vids.at(i-1));
            eid = GetCCedgeAt(v.id, eid, rotOffset);

            // for (auto eid: v.N_Eids) {
                auto& e = mesh->E.at(eid);
                // if (Contains(e.Vids, vids.at(i+1)) || Contains(e.Vids, vids.at(i-1))) continue;
                // std::cout << "GetEdgeId" << std::endl;
                // size_t edgeId = GetEdgeId(vids.at(i), vids.at(i-1));
                // std::cout << "GetCCedgeAt" << std::endl;
                // edgeId = GetCCedgeAt(vids.at(i), edgeId, 1);
                std::vector<size_t> secVids;
                // std::cout << "Tracing secondary link" << std::endl;
                TraceAlongEdge(v, e, secVids);
                b2 = secVids.size()-1;
                if (isRotEdge) {
                    if (b1 > b2) {
                        b1 = b1-b2-1;
                        b2 = 0;
                    } else {
                        b1 = b2-b1;
                        b2 = 0;
                    }
                }
                secVids.insert(secVids.begin(), vids.begin(), vids.begin()+i);
                if (!PrototypeCheckBoundarySingularity(secVids.back()) && ValidatePath(secVids) && !Contains(verticesToAvoid, secVids.back()) && (mesh->V.at(secVids.back()).N_Vids.size() == 5 || mesh->V.at(secVids.back()).N_Vids.size() == 3) && (s.minSize == -1 || s.minSize > b1+b2)) {
                    // std::cout << "checkBoundary: " << checkBoundary << std::endl;
                //     size_t edgeId = GetEdgeId(secVids.at(secVids.size()-1), secVids.at(secVids.size()-2));    
                //     edgeId = GetCCedgeAt(secVids.back(), edgeId, mesh->V.at(secVids.back()).N_Vids.size()-1);
                //     std::vector<size_t> secVidsB;
                    // std::cout << "b1: " << b1 << " b2: " << b2 << std::endl;
                    // TraceAlongEdge(mesh->V.at(secVids.back()), mesh->E.at(edgeId), secVidsB, b1, b2);
                    // std::cout << "secVids: " << secVids.size() << " " << " front: " << secVids.front() << std::endl;
                    // std::cout << "secVidsB: " << secVidsB.size() << " " << "back: " << secVidsB.back() << std::endl;
                    // if (secVids.front() == secVidsB.back()) {
                        s.frontId = secVids.front();
                        s.backId = secVids.back();
                        s.threeId = s.frontId;
                        s.fiveId = s.backId;
                        s.b1 = b1;
                        s.b2 = b2;
                        s.vidsA = secVids;
                        // s.vidsB = secVidsB;
                        s.empty = false;
                        s.minSize = s.b1+s.b2;
                    // }
                }
            // }
        }
    }
}

void SemiGlobalSimplifier::PrototypeTransportLink(Separatrix& primary) {
    for (int i = 1; i < primary.vidsA.size(); i++) {
        size_t singularityId = primary.vidsA.at(i-1);
        size_t moveDir = primary.vidsA.at(i);
        std::vector<size_t> verticesToAvoid = {singularityId, primary.vidsA.back()};
        PrototypeMoveLinkPair(singularityId, moveDir, verticesToAvoid);
    }

}

void SemiGlobalSimplifier::PrototypeMoveLinkPair(size_t singularityId, size_t moveDir, std::vector<size_t> verticesToAvoid) {
    auto& move = mesh->V.at(moveDir);
    auto& singularity = mesh->V.at(singularityId);

    std::vector<size_t> threeFiveIds;
    bool gotDirectPair = false;
    for (auto eid: singularity.N_Eids) {
        auto& e = mesh->E.at(eid);
        if (Contains(e.Vids, move.id)) {
            if (singularity.N_Vids.size() == 3) {
                auto& f = mesh->F.at(mu->GetDifference(singularity.N_Fids, e.N_Fids).at(0));
                int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), singularity.id));
                threeFiveIds.push_back(f.Vids.at((idx+2)%f.Vids.size()));
                threeFiveIds.push_back(f.Vids.at((idx+1)%f.Vids.size()));
                if (mesh->V.at(threeFiveIds.at(0)).N_Vids.size() < 4 || mu->IsSharpFeature(threeFiveIds.at(0))) {
                    threeFiveIds.clear();
                    continue;
                }
                std::shared_ptr<SimplificationOperation> dc = std::make_shared<DiagonalCollapse>(*mesh, *mu, *smoother, f.id, (idx+1)%f.Vids.size(), (idx+3)%f.Vids.size());
                dc->PerformOperation();
                gotDirectPair = true;
                break;
            } else if (singularity.N_Vids.size() == 5) {
                std::vector<size_t> tempEdges;
                for (auto fid: e.N_Fids) {
                    mu->AddContents(tempEdges, mu->GetDifference(mu->GetIntersection(mesh->F.at(fid).Eids, singularity.N_Eids), std::vector<size_t>{eid}));
                }
                std::vector<size_t> faces = e.N_Fids;
                for (auto tmp: tempEdges) mu->AddContents(faces, mesh->E.at(tmp).N_Fids);
                auto& f = mesh->F.at(mu->GetDifference(singularity.N_Fids, faces).at(0));
                threeFiveIds.push_back(singularity.id);
                threeFiveIds.push_back(GetFaceV(singularity.id, f.id, 2));
                if (mesh->V.at(threeFiveIds.at(1)).N_Vids.size() >= 5 || mu->IsSharpFeature(threeFiveIds.at(1))) {
                    threeFiveIds.clear();
                    continue;
                }
                std::shared_ptr<SimplificationOperation> vs = std::make_shared<VertexSplit>(*mesh, *mu, *smoother, singularity.id, tempEdges);
                vs->PerformOperation();
                gotDirectPair = true;
                break;
            }
        }
    }
    if (!gotDirectPair) {
        for (auto fid: singularity.N_Fids) {
            auto& f = mesh->F.at(fid);
            if (GetDiagonalV(singularity.id, f.id) == move.id) {
                if (singularity.N_Vids.size() == 3) {
                    int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), singularity.id));
                    threeFiveIds.push_back(mu->GetDifference(singularity.N_Vids, mu->GetIntersection(singularity.N_Vids, f.Vids)).at(0));
                    threeFiveIds.push_back(f.Vids.at((idx+1)%f.Vids.size()));
                    if (mesh->V.at(threeFiveIds.at(0)).N_Vids.size() < 4 || mu->IsSharpFeature(threeFiveIds.at(0))) {
                        threeFiveIds.clear();
                        continue;
                    }
                    std::shared_ptr<SimplificationOperation> dc = std::make_shared<DiagonalCollapse>(*mesh, *mu, *smoother, f.id, (idx+1)%f.Vids.size(), (idx+3)%f.Vids.size());
                    dc->PerformOperation();
                    break;
                } else if (singularity.N_Vids.size() == 5) {
                    std::vector<size_t> diffEdges = mu->GetIntersection(singularity.N_Eids, f.Eids);
                    for (auto diffEid: diffEdges) {
                        auto& diffE = mesh->E.at(diffEid);
                        auto diffFid = diffE.N_Fids.at(0) == f.id ? diffE.N_Fids.at(1) : diffE.N_Fids.at(0);
                        auto& diffF = mesh->F.at(diffFid);
                        mu->AddContents(diffEdges, mu->GetIntersection(singularity.N_Eids, diffF.Eids));
                    }
                    auto& e = mesh->E.at(mu->GetDifference(singularity.N_Eids, diffEdges).at(0));
                    std::vector<size_t> tempEdges;
                    for (auto efid: e.N_Fids) {
                        mu->AddContents(tempEdges, mu->GetDifference(mu->GetIntersection(mesh->F.at(efid).Eids, singularity.N_Eids), std::vector<size_t>{e.id}));
                    }
                    std::vector<size_t> faces = e.N_Fids;
                    threeFiveIds.push_back(singularity.id);
                    threeFiveIds.push_back((e.Vids.at(0) == singularity.id ? e.Vids.at(1) : e.Vids.at(0)));
                    if (mesh->V.at(threeFiveIds.at(1)).N_Vids.size() >= 5 || mu->IsSharpFeature(threeFiveIds.at(1))) {
                        threeFiveIds.clear();
                        continue;
                    }
                    std::shared_ptr<SimplificationOperation> vs = std::make_shared<VertexSplit>(*mesh, *mu, *smoother, singularity.id, tempEdges);
                    vs->PerformOperation();
                    break;
                }
            }
        }
    }
    if (threeFiveIds.empty()) return;
    if (gotDirectPair) {
        PrototypeTransportDirectPair(threeFiveIds.at(0), threeFiveIds.at(1), false, verticesToAvoid);
    } else {
        PrototypeTransportDiagonalPair(threeFiveIds.at(0), threeFiveIds.at(1), false, verticesToAvoid);
    }
}

void SemiGlobalSimplifier::PrototypeTransportDirectLink(Separatrix& primary, std::vector<std::vector<size_t>>& vec) {
    Separatrix s;
    auto& three = mesh->V.at(primary.threeId);
    auto& five = mesh->V.at(primary.fiveId);

    size_t threeEdge = GetEdgeId(primary.vidsA.at(0), primary.vidsA.at(1));
    size_t fiveEdge = GetEdgeId(primary.vidsA.at(primary.vidsA.size()-1), primary.vidsA.at(primary.vidsA.size()-2));
    std::vector<size_t> threeEdges = mu->GetDifference(three.N_Eids, std::vector<size_t>{threeEdge});
    std::vector<size_t> rotEdges;
    for (auto fid: mesh->E.at(fiveEdge).N_Fids) {
        mu->AddContents(rotEdges, mu->GetDifference(mu->GetIntersection(five.N_Eids, mesh->F.at(fid).Eids), std::vector<size_t>{fiveEdge}));
    }
    std::vector<size_t> nRotEdges = mu->GetDifference(five.N_Eids, mu->GetUnion(rotEdges, std::vector<size_t>{fiveEdge}));
    int rotOffset = 1;
    for (auto eid: threeEdges) {
        if (GetEdgeRots(eid, threeEdge) > 1) {
            rotOffset = 3;
        } else {
            rotOffset = 1;
        }
        auto& e = mesh->E.at(eid);
        size_t vid = e.Vids.at(0) == s.threeId || e.Vids.at(0) == s.fiveId ? e.Vids.at(0) : e.Vids.at(1);
        size_t vertexToAvoid = vid == s.threeId ? s.fiveId : s.threeId;
        PrototypeTraceDirectLinkSeparatrix(three.id, eid, s, vertexToAvoid, rotOffset, false);
    }
    for (auto eid: rotEdges) {
        if (GetEdgeRots(eid, fiveEdge) > 1) {
            rotOffset = 3;
        } else {
            rotOffset = 1;
        }
        auto& e = mesh->E.at(eid);
        size_t vid = e.Vids.at(0) == s.threeId || e.Vids.at(0) == s.fiveId ? e.Vids.at(0) : e.Vids.at(1);
        size_t vertexToAvoid = vid == s.threeId ? s.fiveId : s.threeId;
        PrototypeTraceDirectLinkSeparatrix(five.id, eid, s, vertexToAvoid, rotOffset, true);
    }
    for (auto eid: nRotEdges) {
        if (GetEdgeRots(eid, fiveEdge) > 2) {
            rotOffset = 3;
        } else {
            rotOffset = 1;
        }
        auto& e = mesh->E.at(eid);
        size_t vid = e.Vids.at(0) == s.threeId || e.Vids.at(0) == s.fiveId ? e.Vids.at(0) : e.Vids.at(1);
        size_t vertexToAvoid = vid == s.threeId ? s.fiveId : s.threeId;
        PrototypeTraceDirectLinkSeparatrix(five.id, eid, s, vertexToAvoid, rotOffset, false);
    }
    if (!s.empty) {
        std::cout << "separatrix: " << s.vidsA.size() << " front: " << mesh->V.at(s.vidsA.front()).N_Vids.size() << " back: " << mesh->V.at(s.vidsA.back()).N_Vids.size() << std::endl;
        std::cout << "separatrix b1: " << s.b1 << " b2: " << s.b2 << std::endl;
        // PrototypeSaveSeparatrices(std::vector<std::vector<size_t>>{s.vidsA});
        // std::shared_ptr<ThreeFivePair> tfp = std::make_shared<ThreeFivePair>(*mesh, *mu, *smoother, threeId, fiveId);
        // int n = iters < s.vidsA.size() ? iters : s.vidsA.size();
        vec.push_back(s.vidsA);
        // int n = s.vidsA.size();
        // for (int i = 1; i < s.vidsA.size(); i++) {
        //     size_t moveDirection = s.vidsA.at(i);
        //     size_t sourceDirection = primary.vidsA.at(1);
        //     PrototypeMoveLinkPair(primary, s, s.vidsA.at(i-1), moveDirection, sourceDirection);
        // }
        // if (includeIters) n = iters < s.vidsA.size() ? iters : s.vidsA.size();
        // for (int i = 1; i < n; i++) {
        //     tfp->Move(s.vidsA.at(i), delta);
        // }
    }
}

void SemiGlobalSimplifier::PrototypeTraceDirectLinkSeparatrix(size_t vid, size_t eid, Separatrix& s, size_t vertexToAvoid, int rotOffset, bool isRotEdge) {
    auto& v = mesh->V.at(vid);
    std::vector<size_t> vids;
    TraceAlongEdge(mesh->V.at(vid), mesh->E.at(eid), vids);
    int valenceToCheck = mesh->V.at(vid).N_Vids.size() == 3 ? 5 : 3;
    if (!vids.empty()) {
        int b1 = vids.size()-1;
        int b2 = 0;
        if (!PrototypeCheckBoundarySingularity(vids.back()) && ValidatePath(vids) && vids.back() != vertexToAvoid && (mesh->V.at(vids.back()).N_Vids.size() == valenceToCheck) && (s.minSize == -1 || s.minSize > b1+b2)) {
            // std::cout << "checkBoundary: " << checkBoundary << std::endl;
            s.frontId = vids.front();
            s.backId = vids.back();
            s.threeId = s.frontId;
            s.fiveId = s.backId;
            s.b1 = b1;
            s.b2 = b2;
            s.vidsA = vids;
            s.vidsB.clear();
            s.empty = false;
            s.minSize = s.b1+s.b2;
        }
        int n = vids.size()-1;
        // std::cout << "n: " << n << std::endl;
        for (int i = 1; i < n; i++) {
            b1 = i;
            auto& v = mesh->V.at(vids.at(i));
            size_t eid = GetEdgeId(vids.at(i), vids.at(i-1));
            eid = GetCCedgeAt(v.id, eid, rotOffset);

            // for (auto eid: v.N_Eids) {
                auto& e = mesh->E.at(eid);
                // if (Contains(e.Vids, vids.at(i+1)) || Contains(e.Vids, vids.at(i-1))) continue;
                // std::cout << "GetEdgeId" << std::endl;
                // size_t edgeId = GetEdgeId(vids.at(i), vids.at(i-1));
                // std::cout << "GetCCedgeAt" << std::endl;
                // edgeId = GetCCedgeAt(vids.at(i), edgeId, 1);
                std::vector<size_t> secVids;
                // std::cout << "Tracing secondary link" << std::endl;
                TraceAlongEdge(v, e, secVids);
                b2 = secVids.size()-1;
                if (isRotEdge) {
                    if (b1 > b2) {
                        b1 = b1-b2;
                        b2 = 0;
                    } else {
                        b1 = b2-b1;
                        b2 = 0;
                    }
                }
                secVids.insert(secVids.begin(), vids.begin(), vids.begin()+i);
                if (!PrototypeCheckBoundarySingularity(secVids.back()) && ValidatePath(secVids) && secVids.back() != vertexToAvoid && (mesh->V.at(secVids.back()).N_Vids.size() == valenceToCheck) && (s.minSize == -1 || s.minSize > b1+b2)) {
                    // std::cout << "checkBoundary: " << checkBoundary << std::endl;
                //     size_t edgeId = GetEdgeId(secVids.at(secVids.size()-1), secVids.at(secVids.size()-2));    
                //     edgeId = GetCCedgeAt(secVids.back(), edgeId, mesh->V.at(secVids.back()).N_Vids.size()-1);
                //     std::vector<size_t> secVidsB;
                    // std::cout << "b1: " << b1 << " b2: " << b2 << std::endl;
                    // TraceAlongEdge(mesh->V.at(secVids.back()), mesh->E.at(edgeId), secVidsB, b1, b2);
                    // std::cout << "secVids: " << secVids.size() << " " << " front: " << secVids.front() << std::endl;
                    // std::cout << "secVidsB: " << secVidsB.size() << " " << "back: " << secVidsB.back() << std::endl;
                    // if (secVids.front() == secVidsB.back()) {
                        s.frontId = secVids.front();
                        s.backId = secVids.back();
                        s.threeId = s.frontId;
                        s.fiveId = s.backId;
                        s.b1 = b1;
                        s.b2 = b2;
                        s.vidsA = secVids;
                        // s.vidsB = secVidsB;
                        s.empty = false;
                        s.minSize = s.b1+s.b2;
                    // }
                }
            // }
        }
    }
}

void SemiGlobalSimplifier::PrototypeMoveLinkPair(Separatrix& p, Separatrix& s, size_t singularityId, size_t moveDir, size_t sourceDir) {
    auto& move = mesh->V.at(moveDir);
    auto& source = mesh->V.at(sourceDir);
    auto& singularity = mesh->V.at(singularityId);
    
    bool reversePath = singularity.N_Vids.size() == 5 ? true : false;
    if (reversePath) {
        std::reverse(p.vidsA.begin(), p.vidsA.end());
    }

    std::vector<size_t> threeFiveIds;
    bool gotDirectPair = false;
    for (auto eid: singularity.N_Eids) {
        auto& e = mesh->E.at(eid);
        if (Contains(e.Vids, move.id)) {
            if (singularity.N_Vids.size() == 3) {
                auto& f = mesh->F.at(mu->GetDifference(singularity.N_Fids, e.N_Fids).at(0));
                int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), singularity.id));
                threeFiveIds.push_back(f.Vids.at((idx+2)%f.Vids.size()));
                threeFiveIds.push_back(f.Vids.at((idx+1)%f.Vids.size()));
                // if (mesh->V.at(threeFiveIds.at(0)).N_Vids.size() < 4 || mu->IsSharpFeature(threeFiveIds.at(0))) {
                //     threeFiveIds.clear();
                //     return false;
                // }
                std::shared_ptr<SimplificationOperation> dc = std::make_shared<DiagonalCollapse>(*mesh, *mu, *smoother, f.id, (idx+1)%f.Vids.size(), (idx+3)%f.Vids.size());
                dc->PerformOperation();
                gotDirectPair = true;
                break;
            } else if (singularity.N_Vids.size() == 5) {
                std::vector<size_t> tempEdges;
                for (auto fid: e.N_Fids) {
                    mu->AddContents(tempEdges, mu->GetDifference(mu->GetIntersection(mesh->F.at(fid).Eids, singularity.N_Eids), std::vector<size_t>{eid}));
                }
                std::vector<size_t> faces = e.N_Fids;
                for (auto tmp: tempEdges) mu->AddContents(faces, mesh->E.at(tmp).N_Fids);
                auto& f = mesh->F.at(mu->GetDifference(singularity.N_Fids, faces).at(0));
                threeFiveIds.push_back(singularity.id);
                threeFiveIds.push_back(GetFaceV(singularity.id, f.id, 2));
                // if (mesh->V.at(threeFiveIds.at(1)).N_Vids.size() >= 5 || mu->IsSharpFeature(threeFiveIds.at(1))) {
                //     threeFiveIds.clear();
                //     return false;
                // }
                std::shared_ptr<SimplificationOperation> vs = std::make_shared<VertexSplit>(*mesh, *mu, *smoother, singularity.id, tempEdges);
                vs->PerformOperation();
                gotDirectPair = true;
                break;
            }
        }
    }
    if (!gotDirectPair) {
        for (auto fid: singularity.N_Fids) {
            auto& f = mesh->F.at(fid);
            if (GetDiagonalV(singularity.id, f.id) == move.id) {
                if (singularity.N_Vids.size() == 3) {
                    int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), singularity.id));
                    threeFiveIds.push_back(mu->GetDifference(singularity.N_Vids, mu->GetIntersection(singularity.N_Vids, f.Vids)).at(0));
                    threeFiveIds.push_back(f.Vids.at((idx+1)%f.Vids.size()));
                    // if (mesh->V.at(threeFiveIds.at(0)).N_Vids.size() < 4 || mu->IsSharpFeature(threeFiveIds.at(0))) {
                    //     threeFiveIds.clear();
                    //     return false;
                    // }
                    std::shared_ptr<SimplificationOperation> dc = std::make_shared<DiagonalCollapse>(*mesh, *mu, *smoother, f.id, (idx+1)%f.Vids.size(), (idx+3)%f.Vids.size());
                    dc->PerformOperation();
                    break;
                } else if (singularity.N_Vids.size() == 5) {
                    std::vector<size_t> diffEdges = mu->GetIntersection(singularity.N_Eids, f.Eids);
                    for (auto diffEid: diffEdges) {
                        auto& diffE = mesh->E.at(diffEid);
                        auto diffFid = diffE.N_Fids.at(0) == f.id ? diffE.N_Fids.at(1) : diffE.N_Fids.at(0);
                        auto& diffF = mesh->F.at(diffFid);
                        mu->AddContents(diffEdges, mu->GetIntersection(singularity.N_Eids, diffF.Eids));
                    }
                    auto& e = mesh->E.at(mu->GetDifference(singularity.N_Eids, diffEdges).at(0));
                    std::vector<size_t> tempEdges;
                    for (auto efid: e.N_Fids) {
                        mu->AddContents(tempEdges, mu->GetDifference(mu->GetIntersection(mesh->F.at(efid).Eids, singularity.N_Eids), std::vector<size_t>{e.id}));
                    }
                    std::vector<size_t> faces = e.N_Fids;
                    threeFiveIds.push_back(singularity.id);
                    threeFiveIds.push_back((e.Vids.at(0) == singularity.id ? e.Vids.at(1) : e.Vids.at(0)));
                    // if (mesh->V.at(threeFiveIds.at(1)).N_Vids.size() >= 5 || mu->IsSharpFeature(threeFiveIds.at(1))) {
                    //     threeFiveIds.clear();
                    //     return false;
                    // }
                    std::shared_ptr<SimplificationOperation> vs = std::make_shared<VertexSplit>(*mesh, *mu, *smoother, singularity.id, tempEdges);
                    vs->PerformOperation();
                    break;
                }
            }
        }
    }
    if (gotDirectPair) {
        std::shared_ptr<ThreeFivePair> tfp = std::make_shared<ThreeFivePair>(*mesh, *mu, *smoother, threeFiveIds.at(0), threeFiveIds.at(1));
        int n = p.vidsA.size();
        for (int i = 1; i < n; i++) {
            tfp->Move(p.vidsA.at(i), delta);
        }
    } else {
        std::shared_ptr<DiagonalThreeFivePair> tfp = std::make_shared<DiagonalThreeFivePair>(*mesh, *mu, *smoother, threeFiveIds.at(0), threeFiveIds.at(1));
        int n = p.vidsA.size();
        // if (includeIters) n = iters < s.vidsA.size() ? iters : s.vidsA.size();
        for (int i = 1; i < n; i++) {
            tfp->Move(p.vidsA.at(i), delta);
        }
    }
    for (auto eid: move.N_Eids) {
        std::vector<size_t> vids;
        if (reversePath) {
            TraceAlongEdge(move, mesh->E.at(eid), vids, p.b2, p.b1, p.rots);
            if (!vids.empty() && vids.back() == p.threeId) {
                p.threeId = vids.back();
                p.fiveId = vids.front();
                p.frontId = vids.front();
                p.backId = vids.back();
                p.vidsA = vids;
                std::reverse(p.vidsA.begin(), p.vidsA.end());
                break;
            }
        } else {
            TraceAlongEdge(move, mesh->E.at(eid), vids, p.b1, p.b2, p.rots);
            if (!vids.empty() && vids.back() == p.fiveId) {
                p.threeId = vids.front();
                p.fiveId = vids.back();
                p.frontId = vids.front();
                p.backId = vids.back();
                p.vidsA = vids;
                break;
            }
        }
    }
}

void SemiGlobalSimplifier::PrototypeSaveSeparatrices(std::vector<std::vector<size_t>> vec, std::string in) {
    
    int colorValue = 0;
    int ncolors = 15;
    std::vector<int> colors;
    
    std::vector<Edge> edges;
    for (auto vids: vec) {
        if (vids.empty()) continue;
        for (int i = 1; i < vids.size(); i++) {
            Edge e;
            e.Vids = {vids.at(i-1), vids.at(i)};
            edges.push_back(e);
        }
        std::vector<int> a(vids.size()-1, (colorValue%ncolors));
        colors.insert(colors.end(), a.begin(), a.end());
        colorValue += 1;
        a.clear();
    }
    
    std::ofstream ofs(in+"_singularityLinks.vtk");
    ofs << "# vtk DataFile Version 3.0\n"
        << "singularityLinks" << ".vtk\n"
        << "ASCII\n\n"
        << "DATASET UNSTRUCTURED_GRID\n";
    ofs << "POINTS " << mesh->V.size() << " double\n";
    for (size_t i = 0; i < mesh->V.size(); i++) {
        ofs << std::fixed << std::setprecision(7) <<  mesh->V.at(i).x << " " <<  mesh->V.at(i).y << " " <<  mesh->V.at(i).z << "\n";
    }
    ofs << "CELLS " << edges.size() << " " << 3 * edges.size() << std::endl;
    for (size_t i = 0; i < edges.size(); i++) {
        ofs << "2 " << edges.at(i).Vids.at(0) << " " << edges.at(i).Vids.at(1) << std::endl;
    }
    ofs << "CELL_TYPES " << edges.size() << "\n";
    for (size_t i = 0; i < edges.size(); i++) {
        ofs << "3" << std::endl;
    }

    ofs << "CELL_DATA " << edges.size() << "\n";
    ofs << "SCALARS fixed int\n";
    ofs << "LOOKUP_TABLE default\n";
    for (auto c: colors) {
        ofs << c << "\n";
    }

    ofs.close();
    ofs.clear();
    ofs.open(in+"_out.vtk");

    std::vector<size_t> c_indices;
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

size_t SemiGlobalSimplifier::GetEdgeId(size_t vid, size_t vid2) {
    size_t res;
    auto& v = mesh->V.at(vid);
    for (auto eid: v.N_Eids) {
        auto& e = mesh->E.at(eid);
        if (e.Vids.at(0) == vid2 || e.Vids.at(1) == vid2) return eid;
    }
    return res;
}

size_t SemiGlobalSimplifier::GetFaceId(size_t vid, size_t vid2) {
    size_t res;
    auto& v = mesh->V.at(vid);
    for (auto fid: v.N_Fids) {
        auto& f = mesh->F.at(fid);
        if (Contains(f.Vids, vid2)) return fid;
    }
    return res;
}

size_t SemiGlobalSimplifier::GetCCedgeAt(size_t vid, size_t eid, int counter) {
    auto& v = mesh->V.at(vid);
    for (int i = 0; i < counter; i++) {
        auto& e = mesh->E.at(eid);
        size_t prev = e.Vids.at(0) == v.id ? e.Vids.at(1) : e.Vids.at(0);
        for (auto fid: e.N_Fids) {
            auto& f = mesh->F.at(fid);
            if (GetFaceV(v.id, f.id, 1) == prev) {
                eid = mu->GetDifference(mu->GetIntersection(v.N_Eids, f.Eids), std::vector<size_t>{e.id}).at(0);
                break;
            }
        }
    }
    return eid;
}

size_t SemiGlobalSimplifier::GettCCFaceAt(size_t vid, size_t eid, int counter) {
    size_t res;
    auto& v = mesh->V.at(vid);
    for (int i = 0; i < counter; i++) {
        auto& e = mesh->E.at(eid);
        size_t prev = e.Vids.at(0) == v.id ? e.Vids.at(1) : e.Vids.at(0);
        for (auto fid: e.N_Fids) {
            auto& f = mesh->F.at(fid);
            if (GetFaceV(v.id, f.id, 1) == prev) {
                res = fid;
                eid = GetEdgeId(v.id, GetFaceV(v.id, f.id, 3));
                break;
            }
        }
    }
    return res;
}

bool SemiGlobalSimplifier::CheckMeshValidity() {
    bool res = true;
    for (auto& v: mesh->V) {
        if ((v.N_Vids.size() + v.N_Eids.size() + v.N_Fids.size()) / 3 != v.N_Vids.size()) {
            std::cout << v.id << " : " << v.N_Vids.size() << " " << v.N_Eids.size() << " " << v.N_Fids.size() << std::endl;
            return false;
        }
        for (auto nvid: v.N_Vids) {
            auto& nv = mesh->V.at(nvid);
            if (std::find(nv.N_Vids.begin(), nv.N_Vids.end(), v.id) == nv.N_Vids.end()) {
                std::cout << "v neighbors: " << v.id << " ";
                for (auto id: v.N_Vids) std::cout << id << " ";
                std::cout << std::endl;
                std::cout << "nv neighbors: " << nv.id << " ";
                for (auto id: nv.N_Vids) std::cout << id << " ";
                std::cout << std::endl;
                return false;
            }
        }
    }
    for (auto& f: mesh->F) {
        if (f.Vids.empty() || f.N_Fids.empty()) continue;
        for (auto eid: f.Eids) {
            auto& e = mesh->E.at(eid);
            if (e.Vids.size() != 2) {
                res = false;
                std::cout << f.id << " face vertices: " << f.Vids.size() << std::endl;
                std::cout << f.id << " face edges: " << f.Eids.size() << std::endl;
                std::cout << e.id << " vertices: " << e.Vids.size() << std::endl;
            }
        }
    }
    for (auto& e: mesh->E) {
        if (e.Vids.empty()) continue;
        for (auto fid: e.N_Fids) {
            auto& f = mesh->F.at(fid);
            // if (f.Vids.empty()) {
            //     std::cout << "edge face is empty" << std::endl;
            // }
            for (auto feid: f.Eids) {
                auto& fe = mesh->E.at(feid);
                if (fe.Vids.size() != 2) {
                    res = false;
                    std::cout << f.id << " face vertices: " << f.Vids.size() << std::endl;
                    std::cout << f.id << " face edges: " << f.Eids.size() << std::endl;
                    std::cout << fe.id << " vertices: " << fe.Vids.size() << std::endl;
                }
            }
        }
    }
    for (auto& v: mesh->V) {
        for (auto eid: v.N_Eids) {
            auto& e = mesh->E.at(eid);
            // if (e.N_Fids.size() != 2) std::cout << "edge faces are not 2" << std::endl;
            for (auto fid: e.N_Fids) {
                auto& f = mesh->F.at(fid);
                for (auto feid: f.Eids) {
                    auto& fe = mesh->E.at(feid);
                    if (fe.Vids.size() != 2) {
                        res = false;
                        std::cout << f.id << " face vertices: " << f.Vids.size() << std::endl;
                        std::cout << f.id << " face edges: " << f.Eids.size() << std::endl;
                        std::cout << fe.id << " vertices: " << fe.Vids.size() << std::endl;
                    }
                }
            }    
        }
    }
    return res;
}

void SemiGlobalSimplifier::SaveEdgeMesh() {
    std::string outputf = "edge_mesh.obj";
    std::ofstream ofs(outputf.c_str());
    std::vector<size_t> c_indices;
    for (auto& e: mesh->E) {
        if (e.Vids.empty()) continue;
        c_indices.push_back(e.id);
    }
    for (size_t i = 0; i < mesh->V.size(); i++) {
        ofs << std::fixed << "v " <<  mesh->V.at(i).x << " " <<  mesh->V.at(i).y << " " <<  mesh->V.at(i).z << "\n";
    }
    for (size_t i = 0; i < c_indices.size(); i++) {
        auto& e = mesh->E.at(c_indices.at(i));
        ofs << "l " << e.Vids.at(0)+1 << " " << e.Vids.at(1)+1 << std::endl;
    }
}

void SemiGlobalSimplifier::PrototypeExecute() {
    std::vector<bool> isSingularity = PrototypeSetSingularities();
    std::vector<Separatrix> separatrices;
    // PARALLEL_FOR_BEGIN(0, mesh->V.size()) {
    for (int i = 0; i < mesh->V.size(); i++) {
        if (!isSingularity.at(i) || mesh->V.at(i).N_Vids.size() == 3) continue;
        auto& v = mesh->V.at(i);
        {
            std::lock_guard<std::recursive_mutex> lock(mtx);
            PrototypeExtractSingularityLinks(v.id, separatrices);
            break;
        }
    }
    auto& s = separatrices.at(1);
    std::cout << "b1: " << s.b1 << " b2: " << s.b2 << " rots: " << s.rots << std::endl;
    size_t eid = GetEdgeId(s.vidsA.at(0), s.vidsA.at(1));
    size_t vid;
    std::shared_ptr<SingularityPair> sp;
    for (int i = 0; i < iters; i++) {
        size_t fid = GettCCFaceAt(s.vidsA.at(0), eid, 1);
        auto& e = mesh->E.at(eid);
        auto& f = mesh->F.at(fid);
        vid = GetFaceV(s.vidsA.at(0), f.id, 2);
        eid = GetCCedgeAt(s.vidsA.at(0), eid, 1);
    }
    
    std::vector<size_t> tfIds;
    std::cout << "vid: " << vid << std::endl;
    int secVid = -1;
    if (vid == s.vidsA.at(2)) {
        for (auto fid: mesh->V.at(s.vidsA.at(0)).N_Fids) {
            auto& f = mesh->F.at(fid);
            if (Contains(f.Vids, vid)) {
                secVid = mu->GetDifference(f.Vids, std::vector<size_t>(s.vidsA.begin(), s.vidsA.begin()+3)).at(0);
                break;
            }
        }
    }

    bool isDiagonal = PrototypeGetPairIds(s.vidsA.at(0), vid, tfIds, secVid);
    std::cout << "Diagonal Pair? " << isDiagonal << " toMove: " << vid << std::endl;
    if (isDiagonal) {
        sp = std::make_shared<DiagonalThreeFivePair>(*mesh, *mu, *smoother, tfIds.at(0), tfIds.at(1));
    } else {
        sp = std::make_shared<ThreeFivePair>(*mesh, *mu, *smoother, tfIds.at(0), tfIds.at(1));
    }
    int startIdx = 0;
    std::vector<size_t> tfVids = mu->GetUnion(mesh->V.at(tfIds.at(0)).N_Vids, mesh->V.at(tfIds.at(1)).N_Vids);
    tfVids = mu->GetDifference(tfVids, std::vector<size_t>{vid});
    
    for (int i = 0; i < s.vidsA.size(); i++) {
        for (auto nvid: tfVids) {
            if (Contains(mesh->V.at(nvid).N_Vids, s.vidsA.at(i)) && i > startIdx) {
                s.vidsA.at(i-1) = nvid;
                startIdx = i-1;
            }
        }
    }
    std::cout << "separatrix size: " << s.vidsA.size() << " startIdx: " << startIdx << std::endl;
    for (int i = startIdx; i < s.vidsA.size()-1; i++) {
        sp->Move(s.vidsA.at(i), delta);
    }
    sp->Move(s.vidsA.back(), delta, false);

    // PrototypeSaveSeparatrices(std::vector<std::vector<size_t>>{s.vidsA}, "test");
    // } PARALLEL_FOR_END();
}

void SemiGlobalSimplifier::PrototypeExtractSingularityLinks(size_t vid, std::vector<Separatrix>& separatrices) {
    auto& v = mesh->V.at(vid);
    std::cout << "Extracting Links for " << v.N_Vids.size() << "-singularity: " << v.id << std::endl;
    PrototypeTraceLinks(v.id, separatrices);
    std::cout << "Extracted " << separatrices.size() << " links" << std::endl;
}

void SemiGlobalSimplifier::PrototypeTraceLinks(size_t vid, std::vector<Separatrix>& seps) {
    auto& v = mesh->V.at(vid);
    for (auto eid: v.N_Eids) {
        Separatrix s;
        TraceAlongEdge(v, mesh->E.at(eid), s.vidsA);
        if (s.vidsA.empty()) continue;
        for (int i = 1; i < s.vidsA.size()-1; i++) {
            auto& sv = mesh->V.at(s.vidsA.at(i));
            for (auto sveid: sv.N_Eids) {
                auto& e = mesh->E.at(sveid);
                if (Contains(e.Vids, s.vidsA.at(i-1)) || Contains(e.Vids, s.vidsA.at(i+1))) continue;
                Separatrix es;
                TraceAlongEdge(sv, e, es.vidsA);
                if (es.vidsA.empty()) continue;
                if (mesh->V.at(es.vidsA.back()).isBoundary && mesh->V.at(es.vidsA.back()).N_Vids.size() == 3) continue;
                if (mesh->V.at(es.vidsA.back()).N_Vids.size() != 3 && mesh->V.at(es.vidsA.back()).N_Vids.size() != 5) continue;
                es.b1 = i;
                es.b2 = es.vidsA.size()-1;
                size_t rotEdge = GetEdgeId(s.vidsA.at(i), s.vidsA.at(i-1));
                es.rots = GetEdgeRots(rotEdge, e.id);
                es.vidsA.insert(es.vidsA.begin(), s.vidsA.begin(), s.vidsA.begin()+i);
                es.frontId = es.vidsA.front();
                es.backId = es.vidsA.back();
                PrototypeSetDirMap(es);
                seps.push_back(es);   
            }    
        }
        if (mesh->V.at(s.vidsA.back()).isBoundary && mesh->V.at(s.vidsA.back()).N_Vids.size() == 3) continue;
        if (mesh->V.at(s.vidsA.back()).N_Vids.size() != 3 && mesh->V.at(s.vidsA.back()).N_Vids.size() != 5) continue;
        s.frontId = s.vidsA.front();
        s.backId = s.vidsA.back();
        s.b1 = s.vidsA.size()-1;
        PrototypeSetDirMap(s);
        seps.push_back(s);        
    }
    /*for (auto fid: v.N_Fids) {
        Separatrix s;
        TraceAlongDiagonal(v, mesh->F.at(fid), s.vidsA);
        if (s.vidsA.empty()) continue;
        for (int i = 1; i < s.vidsA.size()-1; i++) {
            auto& sv = mesh->V.at(s.vidsA.at(i));
            for (auto svfid: sv.N_Fids) {
                auto& f = mesh->F.at(svfid);
                if (Contains(f.Vids, s.vidsA.at(i-1)) || Contains(f.Vids, s.vidsA.at(i+1))) continue;
                Separatrix es;
                TraceAlongDiagonal(sv, f, es.vidsA);
                if (es.vidsA.empty()) continue;
                if (mesh->V.at(es.vidsA.back()).isBoundary && mesh->V.at(es.vidsA.back()).N_Vids.size() == 3) continue;
                if (mesh->V.at(es.vidsA.back()).N_Vids.size() != 3 && mesh->V.at(es.vidsA.back()).N_Vids.size() != 5) continue;
                es.b1 = i;
                es.b2 = es.vidsA.size()-1;
                size_t rotFace = GetFaceId(s.vidsA.at(i), s.vidsA.at(i-1));
                es.rots = GetFaceRots(rotFace, f.id, sv.id);
                es.vidsA.insert(es.vidsA.begin(), s.vidsA.begin(), s.vidsA.begin()+i);
                es.frontId = es.vidsA.front();
                es.backId = es.vidsA.back();
                PrototypeSetDirMap(es);
                seps.push_back(es);
            }
        }
        if (mesh->V.at(s.vidsA.back()).isBoundary && mesh->V.at(s.vidsA.back()).N_Vids.size() == 3) continue;
        if (mesh->V.at(s.vidsA.back()).N_Vids.size() != 3 && mesh->V.at(s.vidsA.back()).N_Vids.size() != 5) continue;
        s.frontId = s.vidsA.front();
        s.backId = s.vidsA.back();
        s.b1 = s.vidsA.size()-1;
        PrototypeSetDirMap(s);
        seps.push_back(s);
    }*/
    // std::vector<std::vector<size_t>> sepv;
    // for (auto s: seps) {
    //     sepv.push_back(s.vidsA);
    // }
    // auto& s = seps.at(iters);
    // std::cout << "b1: " << s.b1 << " b2: " << s.b2 << " rots: " << s.rots << std::endl;
    // PrototypeSaveSeparatrices(std::vector<std::vector<size_t>>{s.vidsA}, "test");
}

bool SemiGlobalSimplifier::PrototypeGetPairIds(size_t singularityId, size_t moveDir, std::vector<size_t>& threeFiveIds, int secVid) {
    bool res = false;
    auto& move = mesh->V.at(moveDir);
    auto& singularity = mesh->V.at(singularityId);

    for (auto eid: singularity.N_Eids) {
        auto& e = mesh->E.at(eid);
        if (Contains(e.Vids, move.id)) {
            if (singularity.N_Vids.size() == 3) {
                auto& f = mesh->F.at(mu->GetDifference(singularity.N_Fids, e.N_Fids).at(0));
                int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), singularity.id));
                threeFiveIds.push_back(f.Vids.at((idx+2)%f.Vids.size()));
                threeFiveIds.push_back(f.Vids.at((idx+1)%f.Vids.size()));
                if (mesh->V.at(threeFiveIds.at(0)).N_Vids.size() < 4 || mu->IsSharpFeature(threeFiveIds.at(0))) {
                    threeFiveIds.clear();
                    continue;
                }
                std::shared_ptr<SimplificationOperation> dc = std::make_shared<DiagonalCollapse>(*mesh, *mu, *smoother, f.id, (idx+1)%f.Vids.size(), (idx+3)%f.Vids.size());
                dc->PerformOperation();
                break;
            } else if (singularity.N_Vids.size() == 5) {
                std::vector<size_t> tempEdges;
                for (auto fid: e.N_Fids) {
                    mu->AddContents(tempEdges, mu->GetDifference(mu->GetIntersection(mesh->F.at(fid).Eids, singularity.N_Eids), std::vector<size_t>{eid}));
                }
                std::vector<size_t> faces = e.N_Fids;
                for (auto tmp: tempEdges) mu->AddContents(faces, mesh->E.at(tmp).N_Fids);
                auto& f = mesh->F.at(mu->GetDifference(singularity.N_Fids, faces).at(0));
                threeFiveIds.push_back(singularity.id);
                threeFiveIds.push_back(GetFaceV(singularity.id, f.id, 2));
                if (mesh->V.at(threeFiveIds.at(1)).N_Vids.size() >= 5 || mu->IsSharpFeature(threeFiveIds.at(1))) {
                    threeFiveIds.clear();
                    continue;
                }
                std::shared_ptr<SimplificationOperation> vs = std::make_shared<VertexSplit>(*mesh, *mu, *smoother, singularity.id, tempEdges);
                vs->PerformOperation();
                break;
            }
        }
    }
    if (threeFiveIds.empty()) {
        for (auto fid: singularity.N_Fids) {
            auto& f = mesh->F.at(fid);
            if (GetDiagonalV(singularity.id, f.id) == move.id) {
                if (singularity.N_Vids.size() == 3) {
                    if (secVid != -1) { 
                        bool clockwise = GetFaceV((size_t) secVid, f.id, 1) == singularity.id ? false : true;
                        size_t eid = GetEdgeId(move.id, (size_t) secVid);
                        size_t secFid = mesh->E.at(eid).N_Fids.at(0) == f.id ? mesh->E.at(eid).N_Fids.at(1) : mesh->E.at(eid).N_Fids.at(0);
                        threeFiveIds.push_back(secVid);
                        threeFiveIds.push_back(GetFaceV((size_t) secVid, secFid, 2));
                        std::shared_ptr<SimplificationOperation> s = std::make_shared<EdgeRotation>(*mesh, *mu, *smoother, eid, clockwise);
                        s->PerformOperation();
                        res = true;
                        break;
                    }
                    int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), singularity.id));
                    threeFiveIds.push_back(mu->GetDifference(singularity.N_Vids, mu->GetIntersection(singularity.N_Vids, f.Vids)).at(0));
                    threeFiveIds.push_back(f.Vids.at((idx+1)%f.Vids.size()));
                    if (mesh->V.at(threeFiveIds.at(0)).N_Vids.size() < 4 || mu->IsSharpFeature(threeFiveIds.at(0))) {
                        threeFiveIds.clear();
                        continue;
                    }
                    std::shared_ptr<SimplificationOperation> dc = std::make_shared<DiagonalCollapse>(*mesh, *mu, *smoother, f.id, (idx+1)%f.Vids.size(), (idx+3)%f.Vids.size());
                    dc->PerformOperation();
                    res = true;
                    break;
                } else if (singularity.N_Vids.size() == 5) {
                    if (secVid != -1) { 
                        bool clockwise = GetFaceV((size_t) secVid, f.id, 1) == singularity.id ? true : false;
                        size_t eid = GetEdgeId(singularity.id, (size_t) secVid);
                        size_t secFid = mesh->E.at(eid).N_Fids.at(0) == f.id ? mesh->E.at(eid).N_Fids.at(1) : mesh->E.at(eid).N_Fids.at(0);
                        threeFiveIds.push_back(secVid);
                        threeFiveIds.push_back(GetFaceV((size_t) secVid, secFid, 2));
                        std::shared_ptr<SimplificationOperation> s = std::make_shared<EdgeRotation>(*mesh, *mu, *smoother, eid, clockwise);
                        s->PerformOperation();
                        res = true;
                        break;
                    }
                    std::vector<size_t> diffEdges = mu->GetIntersection(singularity.N_Eids, f.Eids);
                    for (auto diffEid: diffEdges) {
                        auto& diffE = mesh->E.at(diffEid);
                        auto diffFid = diffE.N_Fids.at(0) == f.id ? diffE.N_Fids.at(1) : diffE.N_Fids.at(0);
                        auto& diffF = mesh->F.at(diffFid);
                        mu->AddContents(diffEdges, mu->GetIntersection(singularity.N_Eids, diffF.Eids));
                    }
                    auto& e = mesh->E.at(mu->GetDifference(singularity.N_Eids, diffEdges).at(0));
                    std::vector<size_t> tempEdges;
                    for (auto efid: e.N_Fids) {
                        mu->AddContents(tempEdges, mu->GetDifference(mu->GetIntersection(mesh->F.at(efid).Eids, singularity.N_Eids), std::vector<size_t>{e.id}));
                    }
                    std::vector<size_t> faces = e.N_Fids;
                    threeFiveIds.push_back(singularity.id);
                    threeFiveIds.push_back((e.Vids.at(0) == singularity.id ? e.Vids.at(1) : e.Vids.at(0)));
                    if (mesh->V.at(threeFiveIds.at(1)).N_Vids.size() >= 5 || mu->IsSharpFeature(threeFiveIds.at(1))) {
                        threeFiveIds.clear();
                        continue;
                    }
                    std::shared_ptr<SimplificationOperation> vs = std::make_shared<VertexSplit>(*mesh, *mu, *smoother, singularity.id, tempEdges);
                    vs->PerformOperation();
                    res = true;
                    break;
                }
            }
        }
    }
    return res;
}

void SemiGlobalSimplifier::PrototypeSetDirMap(Separatrix& s) {
    auto& v1 = mesh->V.at(s.frontId);
    auto& v2 = mesh->V.at(s.backId);

    if (v1.N_Vids.size() == 3 && v2.N_Vids.size() == 3) {
        PrototypeSet33DirMap(s);
    } else if (v1.N_Vids.size() == 5 && v2.N_Vids.size() == 5) {
        PrototypeSet55DirMap(s);
    } else {
        PrototypeSet35DirMap(s);
    }
}

void SemiGlobalSimplifier::PrototypeSet33DirMap(Separatrix& s) {
    std::cout << "Setting dirMap for 3-3 separatrix" << std::endl;
    auto& v1 = mesh->V.at(s.frontId);
    auto& v2 = mesh->V.at(s.backId);

    size_t startEid1 = GetEdgeId(v1.id, s.vidsA.at(1));
    size_t startEid2 = GetEdgeId(v2.id, s.vidsA.at(s.vidsA.size()-2));

    std::vector<int> counts = {0,-1*s.b1,-2*s.b1,-1*s.b1,0};

    if (s.rots == 1) {
        counts[0] += (-2*s.b2);
        counts[1] += (-1*s.b2);
        counts[3] += s.b2;
        counts[4] += s.b1 > 1 ? 2*s.b2 : 0;
        startEid2 = GetCCedgeAt(v2.id, startEid2, 1);
    } else if (s.rots > 1) {
        counts[0] += s.b1 > 1 ? 2*s.b2 : 0;
        counts[1] += s.b2;
        counts[3] += (-1*s.b2);
        counts[4] += (-2*s.b2);
        startEid2 = GetCCedgeAt(v2.id, startEid2, 2);
    }
    size_t eid1 = startEid1;
    size_t eid2 = startEid2;
    
    std::vector<size_t> a,b;
    for (int i = 0; i < 3; i++) {
        size_t fid1 = GettCCFaceAt(v1.id, eid1, 1);
        size_t fid2 = GettCCFaceAt(v2.id, eid2, 1);
        auto& e1 = mesh->E.at(eid1);
        auto& e2 = mesh->E.at(eid2);
        auto& f1 = mesh->F.at(fid1);
        auto& f2 = mesh->F.at(fid2);

        a.push_back(GetFaceV(v1.id, f1.id, 1));
        a.push_back(GetFaceV(v1.id, f1.id, 2));
        b.push_back(GetFaceV(v2.id, f2.id, 1));
        b.push_back(GetFaceV(v2.id, f2.id, 2));
        eid1 = GetCCedgeAt(v1.id, eid1, 1);
        eid2 = GetCCedgeAt(v2.id, eid2, 1);
    }
    
    std::vector<size_t> dir1(a.begin()+1,a.end());
    std::vector<size_t> dir2(b.begin()+1,b.end());
    for (int i = 0; i < dir1.size(); i++) {
        if (!mesh->V.at(dir2.at(i)).isBoundary && mesh->V.at(dir2.at(i)).type != FEATURE) {
            s.dirMap[std::to_string(dir1.at(i))] = counts.at(i);
            std::cout << "dir index: " << std::to_string(dir1.at(i)) << " count: " << s.dirMap[std::to_string(dir1.at(i))] << std::endl;
        } else {
            s.dirMap[std::to_string(dir1.at(i))] = -1;
        }
    }
}

void SemiGlobalSimplifier::PrototypeSet55DirMap(Separatrix& s) {
    std::cout << "Setting dirMap for 5-5 separatrix" << std::endl;
    auto& v1 = mesh->V.at(s.frontId);
    auto& v2 = mesh->V.at(s.backId);

    size_t startEid1 = GetEdgeId(v1.id, s.vidsA.at(1));
    size_t startEid2 = GetEdgeId(v2.id, s.vidsA.at(s.vidsA.size()-2));

    std::vector<int> counts = {-2*s.b1,-s.b1,0,s.b1,2*s.b1,s.b1,0,-s.b1,-2*s.b1};
    if (s.rots == 1) {
        counts[1] += s.b2;
        counts[2] += (2*s.b2);
        counts[3] += s.b2;
        counts[5] += (-1*s.b2);
        counts[6] += (-2*s.b2);
        counts[7] += (-1*s.b2);
        counts[8] = s.b1 == 1 ? (-2*s.b2) : counts[8];
        startEid2 = GetCCedgeAt(v2.id, startEid2, 1);
    } else if (s.rots > 1) {
        counts[0] = s.b1 == 1 ? (-2*s.b2) : counts[0];
        counts[1] += (-1*s.b2);
        counts[2] += (-2*s.b2);
        counts[3] += (-1*s.b2);
        counts[5] += s.b2;
        counts[6] += (2*s.b2);
        counts[7] += s.b2;
        startEid2 = GetCCedgeAt(v2.id, startEid2, 4);
    }
    size_t eid1 = startEid1;
    size_t eid2 = startEid2;

    std::vector<size_t> a,b;
    for (int i = 0; i < 5; i++) {
        size_t fid1 = GettCCFaceAt(v1.id, eid1, 1);
        size_t fid2 = GettCCFaceAt(v2.id, eid2, 1);
        auto& e1 = mesh->E.at(eid1);
        auto& e2 = mesh->E.at(eid2);
        auto& f1 = mesh->F.at(fid1);
        auto& f2 = mesh->F.at(fid2);

        a.push_back(GetFaceV(v1.id, f1.id, 2));
        a.push_back(GetFaceV(v1.id, f1.id, 3));
        b.push_back(GetFaceV(v2.id, f2.id, 2));
        b.push_back(GetFaceV(v2.id, f2.id, 3));
        
        eid1 = GetCCedgeAt(v1.id, eid1, 1);
        eid2 = GetCCedgeAt(v2.id, eid2, 1);
    }
    
    std::vector<size_t> dir1(a.begin(),a.end()-1);
    std::vector<size_t> dir2(b.begin(),b.end()-1);
    for (int i = 0; i < dir1.size(); i++) {
        if (!mesh->V.at(dir2.at(i)).isBoundary && mesh->V.at(dir2.at(i)).type != FEATURE) {
            s.dirMap[std::to_string(dir1.at(i))] = counts.at(i);
            std::cout << "dir index: " << std::to_string(dir1.at(i)) << " count: " << s.dirMap[std::to_string(dir1.at(i))] << std::endl;
        } else {
            s.dirMap[std::to_string(dir1.at(i))] = -1;
        }
    }
}

void SemiGlobalSimplifier::PrototypeSet35DirMap(Separatrix& s) {
    std::cout << "Setting dirMap for 3-5 separatrix" << std::endl;
    auto& v1 = mesh->V.at(s.frontId);
    auto& v2 = mesh->V.at(s.backId);

    size_t startEid1 = GetEdgeId(v1.id, s.vidsA.at(1));
    size_t startEid2 = GetEdgeId(v2.id, s.vidsA.at(s.vidsA.size()-2));

    if (v1.N_Vids.size() == 3) {
        std::vector<int> counts = {0,-1*s.b1,-2*s.b1,-1*s.b1,0};
        if (s.rots == 1) {
            counts[0] += (-2*s.b2);
            counts[1] += (-1*s.b2);
            counts[3] += s.b2;
            counts[4] += s.b1 > 1 ? 2*s.b2 : 0;
            startEid2 = GetCCedgeAt(v2.id, startEid2, 1);
        } else if (s.rots > 1) {
            counts[0] += s.b1 > 1 ? 2*s.b2 : 0;
            counts[1] += s.b2;
            counts[3] += (-1*s.b2);
            counts[4] += (-2*s.b2);
            startEid2 = GetCCedgeAt(v2.id, startEid2, 4);
        }
        std::vector<size_t> a,b;
        size_t eid1 = startEid1;
        for (int i = 0; i < 3; i++) {
            size_t fid1 = GettCCFaceAt(v1.id, eid1, 1);
            auto& e1 = mesh->E.at(eid1);
            auto& f1 = mesh->F.at(fid1);
            a.push_back(GetFaceV(v1.id, f1.id, 1));
            a.push_back(GetFaceV(v1.id, f1.id, 2));
            eid1 = GetCCedgeAt(v1.id, eid1, 1);
        }
        size_t eid2 = startEid2;
        for (int i = 0; i < 5; i++) {
            size_t fid2 = GettCCFaceAt(v2.id, eid2, 1);
            auto& e2 = mesh->E.at(eid2);
            auto& f2 = mesh->F.at(fid2);
            b.push_back(GetFaceV(v2.id, f2.id, 1));
            b.push_back(GetFaceV(v2.id, f2.id, 2));
            eid2 = GetCCedgeAt(v2.id, eid2, 1);
        }
        std::vector<size_t> dir1(a.begin()+1,a.end());
        std::vector<size_t> dir2 = {b[6],b[7],b[8],b[1],b[2]};
        
        for (int i = 0; i < dir1.size(); i++) {
            if (!mesh->V.at(dir2.at(i)).isBoundary && mesh->V.at(dir2.at(i)).type != FEATURE) {
                s.dirMap[std::to_string(dir1.at(i))] = counts.at(i);
            std::cout << "dir index: " << std::to_string(dir1.at(i)) << " count: " << s.dirMap[std::to_string(dir1.at(i))] << std::endl;
            } else {
                s.dirMap[std::to_string(dir1.at(i))] = -1;
            }
        }
    } else {
        std::vector<int> counts = {-2*s.b1,-1*s.b1,0,s.b1,2*s.b1,s.b1,0,-1*s.b1,-2*s.b1};
        if (s.rots == 1) {
            counts[1] += s.b2;
            counts[2] += (2*s.b2);
            counts[3] += s.b2;
            counts[5] += (-1*s.b2);
            counts[6] += (-2*s.b2);
            counts[7] += (-1*s.b2);
            counts[8] = s.b1 == 1 ? (-2*s.b2) : counts[8];
            startEid2 = GetCCedgeAt(v2.id, startEid2, 1);
        } else if (s.rots > 1) {
            counts[0] = s.b1 == 1 ? (-2*s.b2) : counts[0];
            counts[1] += (-1*s.b2);
            counts[2] += (-2*s.b2);
            counts[3] += (-1*s.b2);
            counts[5] += s.b2;
            counts[6] += (2*s.b2);
            counts[7] += s.b2;
            startEid2 = GetCCedgeAt(v2.id, startEid2, 2);
        }
        std::vector<size_t> a,b;
        size_t eid1 = startEid1;
        for (int i = 0; i < 5; i++) {
            size_t fid1 = GettCCFaceAt(v1.id, eid1, 1);
            auto& e1 = mesh->E.at(eid1);
            auto& f1 = mesh->F.at(fid1);
            a.push_back(GetFaceV(v1.id, f1.id, 2));
            a.push_back(GetFaceV(v1.id, f1.id, 3));
            eid1 = GetCCedgeAt(v1.id, eid1, 1);
        }
        size_t eid2 = startEid2;
        for (int i = 0; i < 3; i++) {
            size_t fid2 = GettCCFaceAt(v2.id, eid2, 1);
            auto& e2 = mesh->E.at(eid2);
            auto& f2 = mesh->F.at(fid2);
            b.push_back(GetFaceV(v2.id, f2.id, 1));
            b.push_back(GetFaceV(v2.id, f2.id, 2));
            eid2 = GetCCedgeAt(v2.id, eid2, 1);
        }
        std::vector<size_t> dir1(a.begin(),a.end()-1);
        std::vector<size_t> dir2 = {b[2],b[3],b[4],b[4],b[4],b[0],b[0],b[1],b[2]};

        for (int i = 0; i < dir1.size(); i++) {
            if (!mesh->V.at(dir2.at(i)).isBoundary && mesh->V.at(dir2.at(i)).type != FEATURE) {
                s.dirMap[std::to_string(dir1.at(i))] = counts.at(i);
            std::cout << "dir index: " << std::to_string(dir1.at(i)) << " count: " << s.dirMap[std::to_string(dir1.at(i))] << std::endl;
            } else {
                s.dirMap[std::to_string(dir1.at(i))] = -1;
            }
        }
    }
}

std::vector<bool> SemiGlobalSimplifier::PrototypeSetSingularities() {
    std::vector<bool> res(mesh->V.size(), false);
    PARALLEL_FOR_BEGIN(0, mesh->V.size()) {
        auto& v = mesh->V.at(i);
        if ((v.N_Vids.size() == 3 || v.N_Vids.size() == 5) && (!v.isBoundary && v.type != FEATURE)) {
            res.at(i) = true;
        }
    } PARALLEL_FOR_END();
    return res;
}