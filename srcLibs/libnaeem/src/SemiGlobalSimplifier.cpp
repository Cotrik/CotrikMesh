#include <algorithm>
#include <map>
#include <time.h>
#include <queue>
#include <deque>
#include <utility>
#include <memory>
#include <cmath>
#include <random>
#include <chrono>
#include "ThreadPool.h"
#include "verdict.h"
#include "ParallelFor.h"
#include "SemiGlobalSimplifier.h"

SemiGlobalSimplifier::SemiGlobalSimplifier() {}

SemiGlobalSimplifier::SemiGlobalSimplifier(Mesh& mesh_, MeshUtil& mu_, Smoother& smoother_, KDTree& kdtree_) {
    mesh = &mesh_;
    mu = &mu_;
    smoother = &smoother_;
    kd = &kdtree_;
    sp = std::make_unique<SurfaceProjector>(mesh_);
    SetFaceMetrics();
    // renderer.SetMesh(mesh_);
}

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
    SetFaceMetrics();
}

void SemiGlobalSimplifier::SetIters(int iters_) {
    iters = iters_;
    smoother->SetIters(iters_);
}

void SemiGlobalSimplifier::SetFaceMetrics() {
    CheckValidity();
    Smooth();
    PARALLEL_FOR_BEGIN(0, mesh->F.size()) {
        auto& f = mesh->F.at(i);
        auto qV_arr = [&] (Face& f) {
            double coords[4][3];
            for (int i = 0; i < 4; i++) {
                auto& v = mesh->V.at(f.Vids.at(i));
                coords[i][0] = v.x; coords[i][1] = v.y; coords[i][2] = v.z;
            }
            return coords;
        }(f);
        f.area = v_quad_area(4, qV_arr);
        // f.shape = v_quad_shape(4, qV_arr);
        f.shape = v_quad_aspect_ratio(4, qV_arr);
    } PARALLEL_FOR_END();
    PARALLEL_FOR_BEGIN(0, mesh->F.size()) {
        auto& f = mesh->F.at(i);
        auto neighborhood = [&] () {
            std::unordered_set<size_t> neighborhood;
            auto setNeighborhood = [&] (Face& f) {
                for (auto vid: f.Vids) {
                    auto& v = mesh->V.at(vid);
                    neighborhood.insert(v.N_Fids.begin(), v.N_Fids.end());
                }
            };
            setNeighborhood(f);
            for (auto vid: f.Vids) {
                auto& v = mesh->V.at(vid);
                for (auto qid: v.N_Fids) {
                    setNeighborhood(mesh->F.at(qid));
                }
            }
            return neighborhood;
        }();
        double sum_area = 0.0;
        double sum_shape = 0.0;
        for (auto neighbor: neighborhood) {
            sum_area += mesh->F.at(neighbor).area;
            sum_shape += mesh->F.at(neighbor).shape;
        }
        f.avg_area = sum_area / neighborhood.size();
        f.threshold_size = std::pow(std::min(f.area/f.avg_area, f.avg_area/f.area), 2);
        f.threshold_shape = sum_shape / neighborhood.size();
        // f.threshold_shape = f.shape;
    } PARALLEL_FOR_END();
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

bool SemiGlobalSimplifier::SetBoundaryDirectSeparatrixOperations(bool looseCollapse) {
    CheckValidity();

    Op_Q.setMaxQueueOn();
    Op_Q.setSpecialComparisonOn();
    for (auto& v: mesh->V) {
        // std::cout << "on " << v.id << std::endl;
        if (v.N_Fids.empty() || mesh->virtualValence(v) != 4 || v.type == FEATURE || v.isBoundary) continue;
        std::vector<size_t> c1, c2;
        for (auto vid: v.N_Vids) mesh->V.at(vid).type != FEATURE && !mesh->V.at(vid).isBoundary && mesh->valence(mesh->V.at(vid)) == 3 ? c1.push_back(vid) : c2.push_back(vid);
        if (c1.size() != 2 && c2.size() != 2) continue;
        auto& s3_v1 = mesh->V.at(c1.at(0));
        auto& s3_v2 = mesh->V.at(c1.at(1));
        auto& sn_v1= mesh->V.at(c2.at(0));
        auto& sn_v2= mesh->V.at(c2.at(1));
        if (mu->GetDifference(s3_v1.N_Fids, s3_v2.N_Fids).size() != s3_v1.N_Fids.size()) continue;
        // std::cout << "here" << std::endl;
        if ((sn_v1.type == FEATURE || sn_v1.isBoundary) && mesh->virtualValence(sn_v1, v.N_Fids) < 5) continue; 
        if ((sn_v2.type == FEATURE || sn_v2.isBoundary) && mesh->virtualValence(sn_v2, v.N_Fids) < 5) continue; 
        if (sn_v1.type == FEATURE || sn_v1.isBoundary || sn_v2.type == FEATURE || sn_v2.isBoundary) {
            // bool skip = false;
            // if (mesh->virtualValence())
            // int featureCount = 0;
            // for (auto nvid: v.N_Vids) {
            //     if (mesh->V.at(nvid).type == FEATURE || mesh->V.at(nvid).isBoundary) featureCount += 1;
            // }
            // for (auto cvid: c2) {
            //     for (auto eid: v.N_Eids) {
            //         int count = 0;
            //         auto& e = mesh->E.at(eid);
            //         if ((e.Vids.at(0) == v.id && e.Vids.at(1) == cvid) || (e.Vids.at(1) == v.id && e.Vids.at(0) == cvid)) {
            //             for (auto efid: e.N_Fids) {
            //                 auto& f = mesh->F.at(efid);
            //                 for (auto fvid: f.Vids) {
            //                     if (mesh->V.at(fvid).type == FEATURE || mesh->V.at(fvid).isBoundary) count += 1;
            //                 }
            //             }
            //             // if (count == 4 && featureCount == 2) skip = true;
            //             if (count == 4 && mesh->V.at(cvid).idealValence == 4) skip = true;
            //             break;
            //         }
            //     }
            // }
            // if (skip) continue;
            std::shared_ptr<SimplificationOperation> ds = std::make_shared<DirectSeparatrixCollapse>(*mesh, *mu, *smoother, v.id, c1, c2, looseCollapse);
            ds->SetRanking();
            // std::cout << "ranking: " << ds->ranking << std::endl;
            if (ds->ranking < 0) continue;
            Op_Q.insert(ds->ranking, v.id, ds);
        }
    }
    int i = 0;
    std::cout << Op_Q.size() << " boundary direct separatrix operations" << std::endl;
    bool res = Op_Q.size() > 0;
    while (!Op_Q.empty()) {
        // std::cout << "Performing Operation" << std::endl;
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
    std::cout << "Smoothing mesh" << std::endl;
    // std::vector<size_t> smoothv;
    // for (auto& v: mesh->V) smoothv.push_back(v.id);
    // smoother->Smooth(smoothv);
    return res;
}


bool SemiGlobalSimplifier::SetDirectSeparatrixOperations(bool looseCollapse) {
    CheckValidity();

    Op_Q.setMaxQueueOn();
    Op_Q.setSpecialComparisonOn();
    for (auto& v: mesh->V) {
        if (mesh->virtualValence(v) != 4 || v.type == FEATURE || v.isBoundary) continue;
        std::vector<size_t> c1, c2;
        for (auto vid: v.N_Vids) mesh->valence(mesh->V.at(vid)) == 3 ? c1.push_back(vid) : c2.push_back(vid);
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
    bool res = Op_Q.size() > 0;
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
    return res;
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
    std::vector<bool> visited(mesh->V.size(), false);
    for (auto& v: mesh->V) {
        bool skip = true;
        // if (visited.at(v.id)) continue;
        // for (auto nfid: v.N_Fids) {
        //     auto&f = mesh->F.at(nfid);
        //     for (auto fvid: f.Vids) {
        //         visited.at(fvid) = true;
        //     }
        // }
        if (v.N_Vids.size() != 4) continue;
        int count = 0;
        for (auto nvid: v.N_Vids) {
            if (mesh->V.at(nvid).N_Vids.size() == 4) count++; 
        }
        if (count != 4) continue;
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

void SemiGlobalSimplifier::PrototypeSaveMesh(const std::vector<SingularityLink>& links, std::string in) {
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

void SemiGlobalSimplifier::TraceAlongEdge(const Vertex& start_vertex, const Edge& start_edge, bool& crossBoundary, int& boundary_distance, std::vector<size_t>& vids_link, bool checkBoundary) {
    Vertex* currentVertex = (Vertex*) &start_vertex;
    Edge* currentEdge = (Edge*) &start_edge;
    vids_link.push_back(currentVertex->id);
    int dist = 0;
    while (currentEdge != NULL) {
        currentVertex = (Vertex*)&mesh->V.at(currentEdge->Vids[0] == currentVertex->id ? currentEdge->Vids[1] : currentEdge->Vids[0]);
        vids_link.push_back(currentVertex->id);
        if ((currentVertex->isBoundary || currentVertex->type == FEATURE) && boundary_distance == 0) {
            crossBoundary = true;
            boundary_distance = dist;
        }
        if (currentVertex->N_Vids.size() != 4 || (checkBoundary && crossBoundary) || currentVertex->id == start_vertex.id) break;
        std::vector<size_t> edgesToAvoid;
        mu->AddContents(edgesToAvoid, mu->GetIntersection(currentVertex->N_Eids, mesh->F.at(currentEdge->N_Fids.at(0)).Eids));
        mu->AddContents(edgesToAvoid, mu->GetIntersection(currentVertex->N_Eids, mesh->F.at(currentEdge->N_Fids.at(1)).Eids));
        currentEdge = (Edge*)&mesh->E.at(mu->GetDifference(currentVertex->N_Eids, edgesToAvoid).at(0));
        dist++;
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

void SemiGlobalSimplifier::TraceAlongDiagonal(const Vertex& start_vertex, const Face& start_face, int& boundary_distance, std::vector<size_t>& vids_link) {
    Vertex* currentVertex = (Vertex*)&start_vertex;
    Face* currentFace = (Face*)&start_face;
    vids_link.push_back(currentVertex->id);
    int dist = 0;
    while (currentVertex != NULL) {
        currentVertex = (Vertex*)&mesh->V.at(GetFaceV(currentVertex->id, currentFace->id, 2));
        vids_link.push_back(currentVertex->id);
        if ((currentVertex->isBoundary || currentVertex->type == FEATURE) && boundary_distance == 0) boundary_distance = dist;
        if (currentVertex->N_Vids.size() != 4 || currentVertex->id == start_vertex.id) break;
        std::vector<size_t> diffFaces;
        for (auto eid: currentFace->Eids) mu->AddContents(diffFaces, mesh->E.at(eid).N_Fids);
        diffFaces = mu->GetDifference(currentVertex->N_Fids, diffFaces);
        if (diffFaces.empty()) return;
        currentFace = (Face*)&mesh->F.at(diffFaces.at(0)); 
        dist++;
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

void SemiGlobalSimplifier::Smooth(vMesh* m) {
    // std::cout << "Smoothing mesh" << std::endl;
    // std::vector<size_t> smoothv;
    // for (auto& v: mesh->V) {
    //     if (v.N_Vids.empty() || v.N_Eids.empty() || v.N_Fids.empty()) continue;
    //     smoothv.push_back(v.id);
    // }
    // smoother->Smooth(smoothv);
    // smoother->SetMesh(*mesh);
    // smoother->Smooth(std::vector<size_t>{});
    int iters = 10;
    std::vector<size_t> V;
    bool useVM = m != nullptr;
    std::vector<glm::dvec3> new_xyz;
    if (useVM) {
        for (auto it_map = m->vmap.begin(); it_map != m->vmap.end(); it_map++) {
            V.push_back(it_map->second.id);
            new_xyz.push_back(it_map->second.xyz());
        }
    } else {
        V.resize(mesh->V.size());
        new_xyz.resize(mesh->V.size());
        PARALLEL_FOR_BEGIN(0, mesh->V.size()) {
            V.at(i) = i;
            new_xyz.at(i) = mesh->V.at(i).xyz();
        } PARALLEL_FOR_END();
    }
    if (m == nullptr) m = new vMesh(mesh);
    auto normalize = [&] (glm::dvec3 v) {
        if (glm::length(v) == 0.) return glm::dvec3(0.0, 0.0, 0.0);
        return glm::normalize(v);
    };
    auto getVertexNormal = [&] (Vertex& v) {
        glm::dvec3 normal(0.0, 0.0, 0.0);
        auto n = [&] (size_t vid_1, size_t vid_2, size_t vid_3) {
            glm::dvec3 a = m->getVertex(vid_2, useVM).xyz() - m->getVertex(vid_1, useVM).xyz();
            glm::dvec3 b = m->getVertex(vid_3, useVM).xyz() - m->getVertex(vid_1, useVM).xyz();
            return normalize(glm::cross(a, b));
        };
        for (auto fid: v.N_Fids) {
            auto& f = m->getFace(fid, useVM);
            normal += n(f.Vids.at(0), f.Vids.at(1), f.Vids.at(2));
            normal += n(f.Vids.at(0), f.Vids.at(2), f.Vids.at(3));
            normal += n(f.Vids.at(1), f.Vids.at(2), f.Vids.at(3));
            normal += n(f.Vids.at(1), f.Vids.at(3), f.Vids.at(0));
        }
        normal = normalize(normal);
        return normal;
    };
    auto projectToNormal = [&] (Vertex& v, glm::dvec3 p, glm::dvec3 n) {
        return p - (glm::dot(p-v.xyz(), n) * n);
    };
    auto SetPosition = [&] (Vertex& v, int idx) {
        vInfo info_v(mesh, v.id, m, useVM);
        auto nvids = info_v.vids();
        glm::dvec3 move(0.0, 0.0, 0.0);
        auto normal = getVertexNormal(v);
        double count = 0.0;
        for (int i = 0; i < nvids.size(); i++) {
            auto& v_n = m->getVertex(nvids.at(i), useVM);
            auto& v_n_prev = m->getVertex(nvids.at((i-1+nvids.size())%nvids.size()), useVM);
            auto& v_n_next = m->getVertex(nvids.at((i+1)%nvids.size()), useVM);

            auto A = v_n.xyz();
            auto B = v_n_prev.xyz();
            auto C = v_n_next.xyz();

            glm::dvec3 AB = B - A;
            glm::dvec3 BC = C - B;
            glm::dvec3 CA = A - C;

            double a = glm::length(BC);
            double b = glm::length(CA);
            double c = glm::length(AB);
            
            glm::dvec3 incenter = glm::dvec3((a*A.x)+(b*B.x)+(c*C.x), (a*A.y)+(b*B.y)+(c*C.y), (a*A.z)+(b*B.z)+(c*C.z)) / (a+b+c);
            auto dir = normalize(incenter - A);
            dir = (A + (dir * glm::length(v.xyz() - A))) - v.xyz();
            auto target = v.xyz() + dir;
            move += projectToNormal(v, target, normal);
            // move += target;
            count += 1.0;
        }
        if (count > 0.0) move /= count;
        new_xyz.at(idx) = move;
        // v.xyz(move);
    };
    auto SetPositionBoundary = [&] (Vertex& v, int idx) {
        // std::cout << "inside SetPositionBoundary" << std::endl;
        vInfo info_v(mesh, v.id, m, useVM);
        // std::cout << "Got vInfo " << std::endl;
        auto nvids = info_v.vids();
        std::vector<size_t> boundaryVertices;
        for (auto nvid: nvids) {
            auto& nv = m->getVertex(nvid, useVM);
            if (nv.isBoundary || nv.type == FEATURE) boundaryVertices.push_back(nvid);
        }
        // std::cout << "boundaryVertices size: " << boundaryVertices.size() << std::endl;
        if (boundaryVertices.empty() || boundaryVertices.size() > 2) return;
        if (boundaryVertices.size() == 1) {
            SetPosition(v, idx);
            return;
        }
        auto& b1 = m->getVertex(boundaryVertices.at(0), useVM);
        auto& b2 = m->getVertex(boundaryVertices.at(1), useVM);
        glm::dvec3 move(0.0, 0.0, 0.0);
        auto normal = getVertexNormal(v);
        int k = 0;
        for (auto nvid: nvids) {
            if (nvid == b1.id || nvid == b2.id) continue;
            auto& nv = m->getVertex(nvid, useVM);
            glm::dvec3 A(v.xyz() - nv.xyz());
            glm::dvec3 B(b1.xyz() - nv.xyz());
            glm::dvec3 C(b2.xyz() - nv.xyz());
            double theta_1 = atan2(glm::length(glm::cross(A, B)), glm::dot(A, B));
            double theta_2 = atan2(glm::length(glm::cross(A, C)), glm::dot(A, C));
            double theta = (theta_2 - theta_1) / 2.0;
            glm::dvec3 r(v.xyz());
            double l = 0.0;
            if (theta > 0) {
                r = b2.xyz() - v.xyz();
                l += fabs(theta/theta_2) * glm::length(r);
            } else if (theta < 0) {
                r = b1.xyz() - v.xyz();
                l += fabs(theta/theta_1) * glm::length(r);
            }
            move += projectToNormal(v, v.xyz()+l*normalize(r), normal);
            // move += v.xyz()+l*normalize(r);
            k += 1;
        }
        if (k > 0) move /= k;
        new_xyz.at(idx) = move;
        // v.xyz(move);
    };
    auto SetOptimizedPositions = [&] (int idx) {
        // std::cout << "getting vertex at: " << idx << std::endl;
        auto& v = m->getVertex(V.at(idx), useVM);
        // std::cout << "got vertex at: " << idx << std::endl;
        if (v.N_Fids.empty()) return;
        if (v.isBoundary || v.type == FEATURE) {
            // std::cout << "got vertex ";
            SetPositionBoundary(v, idx);
        } else {
            SetPosition(v, idx);
        }
    };
    auto SetProjectedPositions = [&] (int idx) {
        auto& v = m->getVertex(V.at(idx), useVM);
        if (v.N_Fids.empty()) return;
        if (v.isBoundary || v.type == FEATURE) {
            v.xyz(v.xyz()+sp->projectToBoundary(v.xyz()));
        } else {
            auto normal = getVertexNormal(v);
            v.xyz(v.xyz()+sp->projectToSurface(v.xyz(), normal));
        }
    };
    while (iters--) {
        // std::cout << "**********smoothing iter " << iters+1 << "*******************" << std::endl;
        // std::cout << "Getting new positions" << std::endl;
        // std::cout << "vertices: " << mesh->V.size() << " V: " << V.size() << std::endl;
        if (useVM) {
            for (int i = 0; i < V.size(); i++) {
                SetOptimizedPositions(i);
            }
        } else {
            // for (int i = 0; i < mesh->V.size(); i++) {
            PARALLEL_FOR_BEGIN(0, mesh->V.size()) {
                SetOptimizedPositions(i);
            } PARALLEL_FOR_END();
            // }
        }
        // std::cout << "Setting new positions" << std::endl;
        if (useVM) {
            for (int i = 0; i < V.size(); i++) {                
                m->getVertex(V.at(i), useVM).xyz(new_xyz.at(i));
            }
        } else {
            PARALLEL_FOR_BEGIN(0, mesh->V.size()) {
                m->getVertex(V.at(i), useVM).xyz(new_xyz.at(i));
            } PARALLEL_FOR_END();
        }
        // std::cout << "Projecting to surface" << std::endl;
        if (useVM) {
            for (int i = 0; i < V.size(); i++) {
                SetProjectedPositions(i);
            }
        } else {
            PARALLEL_FOR_BEGIN(0, mesh->V.size()) {
                SetProjectedPositions(i);
            } PARALLEL_FOR_END();
        }
    }
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
        if (f.Vids.empty()) continue;
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
    // PARALLEL_FOR_BEGIN(0, mesh->V.size()) {
    // std::map<size_t, double> singularityDirMap; 
    std::priority_queue<Separatrix, std::vector<Separatrix>, SeparatrixComparator> q;
    std::vector<std::vector<size_t>> sepv;
    for (int i = 0; i < mesh->V.size(); i++) {
        // if (!isSingularity.at(i) || mesh->V.at(i).N_Vids.size() == 5) continue;
        if (!isSingularity.at(i)) continue;
        auto& v = mesh->V.at(i);
        {
            std::lock_guard<std::recursive_mutex> lock(mtx);
            std::vector<Separatrix> separatrices;
            PrototypeExtractSingularityLinks(v.id, separatrices);
            std::pair<size_t,double> dirMap = PrototypeGetDir(v.id, separatrices);
            if (dirMap.second > 0) {
                Separatrix* s = PrototypeGetSeparatrix(v.id, separatrices, dirMap);
                if (s && (*s).vidsA.size() > 0) {
                    // PrototypeMove((Separatrix&)s);
                    // q.push(*s);
                    caretaker.saveState(mesh);
                    sepv.push_back((*s).vidsA);
                    sepv.push_back(std::vector<size_t>{(*s).vidsA.front(), (*s).moveDir});
                    PrototypeSaveSeparatrices(sepv, "test");
                    PrototypeMove(*s);
                    PrototypeSaveSeparatrices(sepv, "test2");
                    caretaker.restoreState(mesh);    
                }
            }
            // std::cout << "dir: " << dirMap.first << " value: " << dirMap.second << std::endl;
            
            // for (auto fid: v.N_Fids) {
            //     std::cout << singularityDirMap[GetFaceV(v.id, fid, 1)] << std::endl;
            //     std::cout << singularityDirMap[GetFaceV(v.id, fid, 2)] << std::endl;
            // }
            break;
        }
    }
    // std::cout << "Queue size: " << q.size() << std::endl;
    
    /*while(!q.empty()) {
        auto& s = q.top();
        std::cout << "Separatrix score: " << s.score << " toMove: " << s.moveDir << " rots: " << s.rots << " b1: " << s.b1 << " b2: " << s.b2 << " dirMap: " << s.dirMap.at(s.moveDir) << std::endl;
        // for (auto fid: mesh->V.at(s.vidsA.front()).N_Fids) {
        //     std::cout << GetFaceV(mesh->V.at(s.vidsA.front()).id, fid, 1) << " singularityDirMap: " << singularityDirMap[GetFaceV(mesh->V.at(s.vidsA.front()).id, fid, 1)] << " s.dirMap: " << s.dirMap.at(GetFaceV(mesh->V.at(s.vidsA.front()).id, fid, 1)) << std::endl;
        //     std::cout << GetFaceV(mesh->V.at(s.vidsA.front()).id, fid, 2) << " singularityDirMap: " << singularityDirMap[GetFaceV(mesh->V.at(s.vidsA.front()).id, fid, 2)] << " s.dirMap: " << s.dirMap.at(GetFaceV(mesh->V.at(s.vidsA.front()).id, fid, 2)) << std::endl;
        // }
        PrototypeMove((Separatrix&)s);
        q.pop();
        // break;
    }*/
    // while(FixValences());
    
    /*auto& s = separatrices.at(1);
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
    sp->Move(s.vidsA.back(), delta, false);*/

    // PrototypeSaveSeparatrices(std::vector<std::vector<size_t>>{s.vidsA}, "test");
    // } PARALLEL_FOR_END();
}

void SemiGlobalSimplifier::PrototypeExtractSingularityLinks(size_t vid, std::vector<Separatrix>& separatrices) {
    auto& v = mesh->V.at(vid);
    // std::cout << "Extracting Links for " << v.N_Vids.size() << "-singularity: " << v.id << std::endl;
    PrototypeTraceLinks(v.id, separatrices);
    // std::cout << "Extracted " << separatrices.size() << " links" << std::endl;
}

void SemiGlobalSimplifier::PrototypeTraceLinks(size_t vid, std::vector<Separatrix>& seps) {
    auto& v = mesh->V.at(vid);
    for (auto eid: v.N_Eids) {
        Separatrix s;
        TraceAlongEdge(v, mesh->E.at(eid), s.crossb1, s.b1, s.vidsA);
        if (s.vidsA.empty()) continue;
        bool crossb1 = false;
        for (int i = 1; i < s.vidsA.size()-1; i++) {
            auto& sv = mesh->V.at(s.vidsA.at(i));
            crossb1 = (sv.isBoundary || sv.type == FEATURE);
            for (auto sveid: sv.N_Eids) {
                auto& e = mesh->E.at(sveid);
                if (Contains(e.Vids, s.vidsA.at(i-1)) || Contains(e.Vids, s.vidsA.at(i+1))) continue;
                Separatrix es;
                TraceAlongEdge(sv, e, es.crossb2, es.b1, es.vidsA);
                if (es.vidsA.empty()) continue;
                if (mu->IsSharpFeature(es.vidsA.back())) continue;
                // if (mesh->V.at(es.vidsA.back()).isBoundary && mesh->V.at(es.vidsA.back()).N_Vids.size() != 3) continue;
                if (!mesh->V.at(es.vidsA.back()).isBoundary && mesh->V.at(es.vidsA.back()).N_Vids.size() != 3 && mesh->V.at(es.vidsA.back()).N_Vids.size() != 5) continue;
                es.b1 = i;
                es.b2 = es.vidsA.size()-1;
                size_t rotEdge = GetEdgeId(s.vidsA.at(i), s.vidsA.at(i-1));
                es.rots = GetEdgeRots(rotEdge, e.id);
                es.vidsA.insert(es.vidsA.begin(), s.vidsA.begin(), s.vidsA.begin()+i);
                es.frontId = es.vidsA.front();
                es.backId = es.vidsA.back();
                es.crossb1 = crossb1;
                PrototypeSetDirMap(es);
                seps.push_back(es); 
                // break;
            }    
        }
        // break;
        if (mu->IsSharpFeature(s.vidsA.back())) continue;
        // if (mesh->V.at(s.vidsA.back()).isBoundary && mesh->V.at(s.vidsA.back()).N_Vids.size() != 3) continue;
        if (!mesh->V.at(s.vidsA.back()).isBoundary && mesh->V.at(s.vidsA.back()).N_Vids.size() != 3 && mesh->V.at(s.vidsA.back()).N_Vids.size() != 5) continue;
        s.frontId = s.vidsA.front();
        s.backId = s.vidsA.back();
        s.b1 = s.vidsA.size()-1;
        PrototypeSetDirMap(s);
        seps.push_back(s);        
    }
    for (auto fid: v.N_Fids) {
        Separatrix s;
        TraceAlongDiagonal(v, mesh->F.at(fid), s.b1, s.vidsA);
        if (s.vidsA.empty()) continue;
        for (int i = 1; i < s.vidsA.size()-1; i++) {
            auto& sv = mesh->V.at(s.vidsA.at(i));
            for (auto svfid: sv.N_Fids) {
                auto& f = mesh->F.at(svfid);
                if (Contains(f.Vids, s.vidsA.at(i-1)) || Contains(f.Vids, s.vidsA.at(i+1))) continue;
                Separatrix es;
                TraceAlongDiagonal(sv, f, es.b1, es.vidsA);
                if (es.vidsA.empty()) continue;
                if (mu->IsSharpFeature(es.vidsA.back())) continue;
                // if (mesh->V.at(es.vidsA.back()).isBoundary && mesh->V.at(es.vidsA.back()).N_Vids.size() != 3) continue;
                if (!mesh->V.at(es.vidsA.back()).isBoundary && mesh->V.at(es.vidsA.back()).N_Vids.size() != 3 && mesh->V.at(es.vidsA.back()).N_Vids.size() != 5) continue;
                es.b1 = i;
                es.b2 = es.vidsA.size()-1;
                size_t rotFace = GetFaceId(s.vidsA.at(i), s.vidsA.at(i-1));
                es.rots = GetFaceRots(rotFace, f.id, sv.id);
                es.vidsA.insert(es.vidsA.begin(), s.vidsA.begin(), s.vidsA.begin()+i);
                es.frontId = es.vidsA.front();
                es.backId = es.vidsA.back();
                es.diagonal = true;
                // PrototypeSetDirMap(es);
                seps.push_back(es);
            }
        }
        if (mu->IsSharpFeature(s.vidsA.back())) continue;
        // if (mesh->V.at(s.vidsA.back()).isBoundary && mesh->V.at(s.vidsA.back()).N_Vids.size() != 3) continue;
        if (!mesh->V.at(s.vidsA.back()).isBoundary && mesh->V.at(s.vidsA.back()).N_Vids.size() != 3 && mesh->V.at(s.vidsA.back()).N_Vids.size() != 5) continue;
        s.frontId = s.vidsA.front();
        s.backId = s.vidsA.back();
        s.b1 = s.vidsA.size()-1;
        s.diagonal = true;
        // PrototypeSetDirMap(s);
        seps.push_back(s);
    }
    // std::vector<std::vector<size_t>> sepv;
    // for (auto s: seps) {
    //     sepv.push_back(s.vidsA);
    // }
    // auto& s = seps.at(iters);
    // std::cout << "b1: " << s.b1 << " b2: " << s.b2 << " rots: " << s.rots << std::endl;
    // PrototypeSaveSeparatrices(std::vector<std::vector<size_t>>{s.vidsA}, "test");
    // PrototypeSaveSeparatrices(sepv, "test");
}

bool SemiGlobalSimplifier::PrototypeGetPairIds(size_t singularityId, size_t moveDir, std::vector<size_t>& threeFiveIds, int secVid) {
    bool res = false;
    auto& move = mesh->V.at(moveDir);
    auto& singularity = mesh->V.at(singularityId);
    // std::cout << "Inside PrototypeGetPairIds" << std::endl;
    for (auto eid: singularity.N_Eids) {
        auto& e = mesh->E.at(eid);
        if (Contains(e.Vids, move.id)) {
            if (singularity.N_Vids.size() == 3) {
                // std::cout << "3-singularity and dest is adjacent" << std::endl;
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
                // std::cout << "5-singularity and dest is adjacent" << std::endl;
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
                    // if (secVid != -1) { 
                    //     bool clockwise = GetFaceV((size_t) secVid, f.id, 1) == singularity.id ? false : true;
                    //     size_t eid = GetEdgeId(move.id, (size_t) secVid);
                    //     size_t secFid = mesh->E.at(eid).N_Fids.at(0) == f.id ? mesh->E.at(eid).N_Fids.at(1) : mesh->E.at(eid).N_Fids.at(0);
                    //     threeFiveIds.push_back(secVid);
                    //     threeFiveIds.push_back(GetFaceV((size_t) secVid, secFid, 2));
                    //     std::shared_ptr<SimplificationOperation> s = std::make_shared<EdgeRotation>(*mesh, *mu, *smoother, eid, clockwise);
                    //     s->PerformOperation();
                    //     res = true;
                    //     break;
                    // }
                    // std::cout << "3-singularity and dest is diagonal" << std::endl;
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
                    // if (secVid != -1) { 
                    //     bool clockwise = GetFaceV((size_t) secVid, f.id, 1) == singularity.id ? true : false;
                    //     size_t eid = GetEdgeId(singularity.id, (size_t) secVid);
                    //     size_t secFid = mesh->E.at(eid).N_Fids.at(0) == f.id ? mesh->E.at(eid).N_Fids.at(1) : mesh->E.at(eid).N_Fids.at(0);
                    //     threeFiveIds.push_back(secVid);
                    //     threeFiveIds.push_back(GetFaceV((size_t) secVid, secFid, 2));
                    //     std::shared_ptr<SimplificationOperation> s = std::make_shared<EdgeRotation>(*mesh, *mu, *smoother, eid, clockwise);
                    //     s->PerformOperation();
                    //     res = true;
                    //     break;
                    // }
                    // std::cout << "5-singularity and dest is diagonal" << std::endl;
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

    // if (v1.isBoundary || v1.type == FEATURE || v2.isBoundary || v2.type == FEATURE) return;
    // if (v1.N_Vids.size() != 3 && v1.N_Vids.size() != 5 && v2.N_Vids.size() != 3 && v2.N_Vids.size() != 5) return;
    for (auto fid: v1.N_Fids) {
        s.dirMap[GetFaceV(v1.id, fid, 1)] = 0;
        s.dirMap[GetFaceV(v1.id, fid, 2)] = 0;
    }

    if (v1.N_Vids.size() == 3 && v2.N_Vids.size() == 3) {
        PrototypeSet33DirMap(s);
    } else if (v1.N_Vids.size() == 5 && v2.N_Vids.size() == 5) {
        PrototypeSet55DirMap(s);
    } else {
        PrototypeSet35DirMap(s);
    }
}

void SemiGlobalSimplifier::PrototypeSet33DirMap(Separatrix& s) {
    // std::cout << "Setting dirMap for 3-3 separatrix" << std::endl;
    auto& v1 = mesh->V.at(s.frontId);
    auto& v2 = mesh->V.at(s.backId);

    size_t startEid1 = GetEdgeId(v1.id, s.vidsA.at(1));
    size_t startEid2 = GetEdgeId(v2.id, s.vidsA.at(s.vidsA.size()-2));

    double b1 = (double)s.b1;
    double b2 = (double)s.b2;
    std::vector<double> counts = {0,-1*b1,-2*b1,-1*b1,0};

    if (s.rots == 1) {
        counts[0] += (-2*b2);
        counts[1] += (-1*b2);
        counts[3] += b2;
        counts[4] += b1 > 1 ? 2*b2 : 0;
        startEid2 = GetCCedgeAt(v2.id, startEid2, 1);
    } else if (s.rots > 1) {
        counts[0] += b1 > 1 ? 2*b2 : 0;
        counts[1] += b2;
        counts[3] += (-1*b2);
        counts[4] += (-2*b2);
        startEid2 = GetCCedgeAt(v2.id, startEid2, 2);
    }
    double min = counts.at(0);
    for (auto el: counts) {
        if (el < min) min = el;
    }
    if (s.crossb1) {
        counts[0] = min-1; counts[1] = min-1; counts[3] = min-1; counts[4] = min-1;
    }
    if (s.crossb2) {
        counts[1] = min-1; counts[3] = min-1;
        counts[0] = counts[0] > 0 ? counts[0] : min-1; counts[2] = min-1; counts[4] = counts[4] > 0 ? counts[4] : min-1;
    }
    size_t eid1 = startEid1;
    size_t eid2 = startEid2;
    
    std::vector<size_t> a,b;
    if (v2.isBoundary) {
        for (int i = 0; i < 3; i++) {
            size_t fid1 = GettCCFaceAt(v1.id, eid1, 1);
            auto& e1 = mesh->E.at(eid1);
            auto& f1 = mesh->F.at(fid1);

            a.push_back(GetFaceV(v1.id, f1.id, 1));
            a.push_back(GetFaceV(v1.id, f1.id, 2));
            eid1 = GetCCedgeAt(v1.id, eid1, 1);
        }
        std::vector<size_t> dir1(a.begin()+1,a.end());
        for (int i = 0; i < dir1.size(); i++) {
            if (counts.at(i) > min-1) {
                s.dirMap[dir1.at(i)] = 1 / (fabs(mu->delta + counts.at(i)) + 1);
            }
        }
        return;
    }
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
        if (!mesh->V.at(dir2.at(i)).isBoundary && mesh->V.at(dir2.at(i)).type != FEATURE && counts.at(i) > min-1) {
            s.dirMap[dir1.at(i)] = 1 / (fabs(mu->delta + counts.at(i)) + 1);
            // std::cout << "dir index: " << std::to_string(dir1.at(i)) << " count: " << s.dirMap[std::to_string(dir1.at(i))] << std::endl;
        }
        // else {
        //     s.dirMap[std::to_string(dir1.at(i))] = 0;
        // }
    }
}

void SemiGlobalSimplifier::PrototypeSet55DirMap(Separatrix& s) {
    // std::cout << "Setting dirMap for 5-5 separatrix" << std::endl;
    auto& v1 = mesh->V.at(s.frontId);
    auto& v2 = mesh->V.at(s.backId);

    size_t startEid1 = GetEdgeId(v1.id, s.vidsA.at(1));
    size_t startEid2 = GetEdgeId(v2.id, s.vidsA.at(s.vidsA.size()-2));

    double b1 = (double)s.b1;
    double b2 = (double)s.b2;
    std::vector<double> counts = {-2*b1,-b1,0,b1,2*b1,b1,0,-b1,-2*b1};
    if (s.rots == 1) {
        counts[1] += b2;
        counts[2] += (2*b2);
        counts[3] += b2;
        counts[5] += (-1*b2);
        counts[6] += (-2*b2);
        counts[7] += (-1*b2);
        counts[8] = (int)b1 == 1 ? (-2*b2) : counts[8];
        startEid2 = GetCCedgeAt(v2.id, startEid2, 1);
    } else if (s.rots > 1) {
        counts[0] = (int)b1 == 1 ? (-2*b2) : counts[0];
        counts[1] += (-1*b2);
        counts[2] += (-2*b2);
        counts[3] += (-1*b2);
        counts[5] += b2;
        counts[6] += (2*b2);
        counts[7] += b2;
        startEid2 = GetCCedgeAt(v2.id, startEid2, 4);
    }
    double min = counts.at(0);
    for (auto el: counts) {
        if (el < min) min = el;
    }
    // std::cout << min << std::endl;
    if (s.crossb1) {
        counts[1] = min-1; counts[2] = min-1; counts[3] = min-1; counts[5] = min-1; counts[7] = min-1; counts[6] = min-1;
    }
    if (s.crossb2) {
        counts[1] = min-1; counts[3] = min-1; counts[5] = min-1; counts[7] = min-1;
        counts[0] = (s.rots > 1 && (int)b1 == 1) ? counts[0] : min-1; counts[4] = min-1; counts[8] = (s.rots == 1 && (int)b1 == 1) ? counts[8] : min-1; 
    }
    size_t eid1 = startEid1;
    size_t eid2 = startEid2;

    std::vector<size_t> a,b;
    if (v2.isBoundary) {
        for (int i = 0; i < 5; i++) {
            size_t fid1 = GettCCFaceAt(v1.id, eid1, 1);
            auto& e1 = mesh->E.at(eid1);
            auto& f1 = mesh->F.at(fid1);

            a.push_back(GetFaceV(v1.id, f1.id, 2));
            a.push_back(GetFaceV(v1.id, f1.id, 3));
            
            eid1 = GetCCedgeAt(v1.id, eid1, 1);
        }
        std::vector<size_t> dir1(a.begin(),a.end()-1);
        for (int i = 0; i < dir1.size(); i++) {
            if (counts.at(i) > min-1) {
                s.dirMap[dir1.at(i)] = 1 / (fabs(mu->delta + counts.at(i)) + 1);
            }
        }
        return;
    }
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
        if (!mesh->V.at(dir2.at(i)).isBoundary && mesh->V.at(dir2.at(i)).type != FEATURE && counts.at(i) > min-1) {
            s.dirMap[dir1.at(i)] = 1 / (fabs(mu->delta + counts.at(i)) + 1);
            // std::cout << "dir index: " << std::to_string(dir1.at(i)) << " count: " << s.dirMap[std::to_string(dir1.at(i))] << std::endl;
        } 
        // else {
        //     s.dirMap[std::to_string(dir1.at(i))] = 0;
        // }
    }
}

void SemiGlobalSimplifier::PrototypeSet35DirMap(Separatrix& s) {
    // std::cout << "Setting dirMap for 3-5 separatrix" << std::endl;
    auto& v1 = mesh->V.at(s.frontId);
    auto& v2 = mesh->V.at(s.backId);

    size_t startEid1 = GetEdgeId(v1.id, s.vidsA.at(1));
    size_t startEid2 = GetEdgeId(v2.id, s.vidsA.at(s.vidsA.size()-2));

    double b1 = (double)s.b1;
    double b2 = (double)s.b2;
    if (v1.N_Vids.size() == 3) {
        std::vector<double> counts = {0,-1*b1,-2*b1,-1*b1,0};
        if (s.rots == 1) {
            counts[0] += (-2*b2);
            counts[1] += (-1*b2);
            counts[3] += b2;
            counts[4] += b1 > 1 ? 2*b2 : 0;
            startEid2 = GetCCedgeAt(v2.id, startEid2, 1);
        } else if (s.rots > 1) {
            counts[0] += b1 > 1 ? 2*b2 : 0;
            counts[1] += b2;
            counts[3] += (-1*b2);
            counts[4] += (-2*b2);
            startEid2 = GetCCedgeAt(v2.id, startEid2, 4);
        }
        double min = counts.at(0);
        for (auto el: counts) {
            if (el < min) min = el;
        }
        if (s.crossb1) {
            counts[0] = min-1; counts[1] = min-1; counts[3] = min-1; counts[4] = min-1;
        }
        if (s.crossb2) {
            counts[1] = min-1; counts[3] = min-1;
            counts[0] = counts[0] > 0 ? counts[0] : min-1; counts[2] = min-1; counts[4] = counts[4] > 0 ? counts[4] : min-1;
        }
        std::vector<size_t> a,b;
        size_t eid1 = startEid1;
        if (v2.isBoundary) {
            for (int i = 0; i < 3; i++) {
                size_t fid1 = GettCCFaceAt(v1.id, eid1, 1);
                auto& e1 = mesh->E.at(eid1);
                auto& f1 = mesh->F.at(fid1);

                a.push_back(GetFaceV(v1.id, f1.id, 1));
                a.push_back(GetFaceV(v1.id, f1.id, 2));
                eid1 = GetCCedgeAt(v1.id, eid1, 1);
            }
            std::vector<size_t> dir1(a.begin()+1,a.end());
            for (int i = 0; i < dir1.size(); i++) {
                if (counts.at(i) > min-1) {
                    s.dirMap[dir1.at(i)] = 1 / (fabs(mu->delta + counts.at(i)) + 1);
                }
            }
            return;
        }
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
            if (!mesh->V.at(dir2.at(i)).isBoundary && mesh->V.at(dir2.at(i)).type != FEATURE && counts.at(i) > min-1) {
                s.dirMap[dir1.at(i)] = 1 / (fabs(mu->delta + counts.at(i)) + 1);
                // std::cout << "dir index: " << dir1.at(i) << " count: " << s.dirMap[dir1.at(i)] << std::endl;
            }
            // else {
            //     s.dirMap[std::to_string(dir1.at(i))] = 0;
            // }
        }
    } else {
        std::vector<double> counts = {-2*b1,-1*b1,0,b1,2*b1,b1,0,-1*b1,-2*b1};
        if (s.rots == 1) {
            counts[1] += b2;
            counts[2] += (2*b2);
            counts[3] += b2;
            counts[5] += (-1*b2);
            counts[6] += (-2*b2);
            counts[7] += (-1*b2);
            counts[8] = b1 == 1 ? (-2*b2) : counts[8];
            startEid2 = GetCCedgeAt(v2.id, startEid2, 1);
        } else if (s.rots > 1) {
            counts[0] = b1 == 1 ? (-2*b2) : counts[0];
            counts[1] += (-1*b2);
            counts[2] += (-2*b2);
            counts[3] += (-1*b2);
            counts[5] += b2;
            counts[6] += (2*b2);
            counts[7] += b2;
            startEid2 = GetCCedgeAt(v2.id, startEid2, 2);
        }
        // std::cout << "s.crossb1: " << s.crossb1 << " s.crossb2: " << s.crossb2 << std::endl;
        double min = counts.at(0);
        for (auto el: counts) {
            if (el < min) min = el;
        }
        // std::cout << min << std::endl;
        if (s.crossb1) {
            counts[1] = min-1; counts[2] = min-1; counts[3] = min-1; counts[5] = min-1; counts[7] = min-1; counts[6] = min-1;
        }
        if (s.crossb2) {
            counts[1] = min-1; counts[3] = min-1; counts[5] = min-1; counts[7] = min-1;
            counts[0] = (s.rots > 1 && (int)b1 == 1) ? counts[0] : min-1; counts[4] = min-1; counts[8] = (s.rots == 1 && (int)b1 == 1) ? counts[8] : min-1; 
        }
        // for (auto el: counts) {
        //     std::cout << el << " ";
        // }
        // std::cout << std::endl;
        std::vector<size_t> a,b;
        size_t eid1 = startEid1;
        if (v2.isBoundary) {
            for (int i = 0; i < 5; i++) {
                size_t fid1 = GettCCFaceAt(v1.id, eid1, 1);
                auto& e1 = mesh->E.at(eid1);
                auto& f1 = mesh->F.at(fid1);

                a.push_back(GetFaceV(v1.id, f1.id, 2));
                a.push_back(GetFaceV(v1.id, f1.id, 3));
                
                eid1 = GetCCedgeAt(v1.id, eid1, 1);
            }
            std::vector<size_t> dir1(a.begin(),a.end()-1);
            for (int i = 0; i < dir1.size(); i++) {
                if (counts.at(i) > min-1) {
                    s.dirMap[dir1.at(i)] = 1 / (fabs(mu->delta + counts.at(i)) + 1);
                    // std::cout << "dir index: " << dir1.at(i) << " count: " << s.dirMap[dir1.at(i)] << std::endl;
                }
            }
            return;
        }
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
            if (!mesh->V.at(dir2.at(i)).isBoundary && mesh->V.at(dir2.at(i)).type != FEATURE && counts.at(i) > min-1) {
                s.dirMap[dir1.at(i)] = 1 / (fabs(mu->delta + counts.at(i)) + 1);
                // std::cout << "dir index: " << dir1.at(i) << " count: " << s.dirMap[dir1.at(i)] << std::endl;
            } 
            // else {
            //     s.dirMap[std::to_string(dir1.at(i))] = 0;
            // }
        }
    }
}

std::pair<size_t,double> SemiGlobalSimplifier::PrototypeGetDir(size_t vid, std::vector<Separatrix>& separatrices) {
    std::pair<size_t,double> res = std::make_pair<size_t,double>(0,0);
    auto& v = mesh->V.at(vid);
    int oppValence = v.N_Vids.size() == 3 ? 5 : 3;
    std::map<size_t, std::vector<double>> dirValueMap;
    std::vector<size_t> nVids;
    for (auto fid: v.N_Fids) {
        nVids.push_back(GetFaceV(vid, fid, 1));
        nVids.push_back(GetFaceV(vid, fid, 2));
    }
    for (auto vid: nVids) {
        dirValueMap[vid] = {0.0,0.0,0.0,0.0};
        // singularityDirMap[vid] = 0.0;
    }
    for (auto s: separatrices) {
        auto& backV = mesh->V.at(s.vidsA.back());
        int boundaryValence = backV.isBoundary && ((backV.N_Vids.size() == 4 && oppValence == 5) || (backV.N_Vids.size() == 2 && !mu->IsSharpFeature(backV.id) && oppValence == 3)) ? backV.N_Vids.size() : -1;
        // if (boundaryValence > -1) {
        //     std::cout << "backV boundary: " << backV.isBoundary << " N_Vids: " << backV.N_Vids.size() << std::endl;
        // }
        size_t dirV = s.vidsA.at(1);
        if ((backV.isBoundary && boundaryValence == -1) || (!backV.isBoundary && backV.N_Vids.size() != oppValence)) {
            double n = ++dirValueMap[dirV].at(0);
            double alpha = dirValueMap[dirV].at(1);
            double a = 1-(1/((double)s.b1+(double)s.b2));
            dirValueMap[dirV].at(1) = (((n-1)*alpha)+a)/n;
        } else if ((backV.isBoundary && boundaryValence > -1) || (!backV.isBoundary && backV.N_Vids.size() == oppValence)) {
            double n = ++dirValueMap[dirV].at(2);
            double beta = dirValueMap[dirV].at(3);
            double b = 1/((double)s.b1+(double)s.b2);
            dirValueMap[dirV].at(3) = (((n-1)*beta)+b)/n;
        }
        
    }
    double score = -1;
    for (auto vid: nVids) {
        double value = (int) dirValueMap[vid].at(3) == 1 ? dirValueMap[vid].at(3) : dirValueMap[vid].at(1)*dirValueMap[vid].at(3);
        // std::cout << "dirValueMap[vid].at(1): " << dirValueMap[vid].at(1) << " dirValueMap[vid].at(3): " << dirValueMap[vid].at(3) << std::endl;
        // std::cout << "value: " << value << std::endl;
        if (value > score) {
            score = value;
            res.first = vid;
            res.second = score;
        } 
        // singularityDirMap[vid] = ;
        // std::cout << "singularityDirMap at vid " << vid << " : " << singularityDirMap[vid] << std::endl;
    }
    return res;
}

Separatrix* SemiGlobalSimplifier::PrototypeGetSeparatrix(size_t vid, std::vector<Separatrix>& separatrices, std::pair<size_t,double>& dirMap) {
    Separatrix* res;
    // std::vector<size_t> nVids;
    // for (auto fid: mesh->V.at(vid).N_Fids) {
    //     nVids.push_back(GetFaceV(vid, fid, 1));
    //     nVids.push_back(GetFaceV(vid, fid, 2));
    // }
    double maxScore = 0.0;
    for (auto& s: separatrices) {
        if (s.diagonal) continue;
        // std::cout << "dirMap size: " << s.dirMap.size() << std::endl;
        // for (auto nvid: nVids) {
        //     std::cout << "nvid: " << nvid << " dirMap: " << s.dirMap[nvid] << std::endl;
        // }
        // std::cout << "singularityDirMap size: " << singularityDirMap.size() << std::endl;
        // for (auto nvid: nVids) {
        //     std::cout << "nvid: " << nvid << " singularityDirMap: " << singularityDirMap[nvid] << std::endl;
        // }
        // for (auto nvid: nVids) {
            // double score = singularityDirMap[nvid] * s.dirMap[nvid];
            double score = s.dirMap[dirMap.first];
            // double score = (1/((double)s.b1+(double)s.b2)) * dirMap.second * s.dirMap[dirMap.first];
            // std::cout << "nvid: " << nvid << " singularityDirMap[nvid]: " << singularityDirMap[nvid] << " s.dirMap[nvid]: " << s.dirMap[nvid] << " score: " << score << " s.score: " << s.score << std::endl;
            
            if (score > maxScore) {
                // bool sameSeparatrixEnd = false;
                // for (auto& s2: separatrices) {
                //     if (Contains(s2.vidsA, dirMap.first) && s2.vidsA.back() == s.vidsA.back() && mesh->V.at(s.vidsA.front()).N_Vids.size() != mesh->V.at(s.vidsA.back()).N_Vids.size()) {
                //         sameSeparatrixEnd = true;
                //         break;
                //     }
                // }
                // if (sameSeparatrixEnd) continue;
                s.moveDir = dirMap.first;
                s.score = dirMap.second;
                // s.score = score;
                maxScore = s.score;
                res = &s;
            }
        // }
    }
    return res;
}

void SemiGlobalSimplifier::PrototypeMove(Separatrix& s) {
    auto& front = mesh->V.at(s.vidsA.front());
    auto& back = mesh->V.at(s.vidsA.back());
    if (front.N_Vids.empty() || (front.N_Vids.size() != 3 && front.N_Vids.size() != 5)) return;
    size_t vid = s.moveDir;
    // bool foundDir = false;
    // std::cout << "moveDir: " << vid << std::endl;
    // for (auto fid: front.N_Fids) {
        // std::cout << "f Vids: ";
        // for (auto fvid: mesh->F.at(fid).Vids) std::cout << fvid << " ";
        // std::cout << std::endl;
        // if (Contains(mesh->F.at(fid).Vids, vid)) foundDir = true;
    // }
    // if (!foundDir) return;
    // std::cout << "front: " << front.id << " N_Vids: " << front.N_Vids.size() << std::endl;
    // std::cout << "back: " << back.id << " N_Vids: " << back.N_Vids.size() << std::endl;
    if (!ValidatePath(s.vidsA, false)) return;
    std::shared_ptr<SingularityPair> sp;
    
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

    std::vector<size_t> tfIds;
    bool isDiagonal = PrototypeGetPairIds(s.vidsA.at(0), vid, tfIds, secVid);
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
    for (int i = startIdx; i < s.vidsA.size()-1; i++) {
        std::cout << "Moving to " << s.vidsA.at(i) << " Vids: " << mesh->V.at(s.vidsA.at(i)).N_Vids.size() << std::endl;
        PrototypeSaveSeparatrices(std::vector<std::vector<size_t>>{}, "test2");
        sp->Move(s.vidsA.at(i), delta);
    }
    sp->Move(s.vidsA.back(), delta, false);
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

void SemiGlobalSimplifier::func1() {
    // mu->SetV_Scores();
    std::vector<Path> paths = PrototypeGetPaths();
    // std::cout << "Got " << paths.size() << " paths" << std::endl;
    // PrototypeSaveSeparatrices(std::vector<std::vector<size_t>>{paths.at(iters).vids}, "test");
    // std::vector<std::vector<size_t>> paths_;
    // for (auto path: paths) paths_.push_back(path.vids);
    // PrototypeSaveSeparatrices(paths_, "test");
}

// a function named GetPaths() that retrieves all paths between singularities

std::vector<Path> SemiGlobalSimplifier::PrototypeGetPaths() {
    // std::vector<bool> singularities = PrototypeSetSingularities();
    std::vector<Path> paths;
    std::vector<vQEM> vQEMs;
    /*for (size_t i = 0; i < singularities.size(); i++) {
        if (singularities.at(i)) {
            std::cout << "Getting paths for singularity at: " << i << std::endl;
            PrototypeGetPathsAt(i, paths);
            // auto score = mu->GetVertexScore(i);
            // std::cout << "Same score: " << score->same_singularity_score << " Diff score: " << score->diff_singularity_score << " Boundary score: " << score->boundary_score << std::endl;
            // std::cout << "Total score: " << (score->same_singularity_score+score->diff_singularity_score+score->boundary_score) / 3 << std::endl;
            break;
        }
    }*/
    struct path_v {
        size_t curr, prev;
        int dist = 0, dist_b = 0;
        std::vector<size_t> vids;
        bool continuous = false;
    };

    // lambda function to get vertex neighboring vertices in order

    auto getVertexNeighbors = [&] (size_t vid, size_t start, bool getDiagonals = false) -> std::vector<size_t> {
        auto& v = this->mesh->V.at(vid);
        std::vector<size_t> res;
        size_t eid = this->GetEdgeId(v.id, start);
        auto getFace_V = [] (Face& f, int idx, int offset) -> size_t {
            return f.Vids.at((idx+offset)%f.Vids.size());
        };
        for (int i = 0; i < v.N_Eids.size(); i++) {
            auto& e = mesh->E.at(eid);
            size_t prev = e.Vids.at(0) == v.id ? e.Vids.at(1) : e.Vids.at(0);
            res.push_back(prev);
            for (auto fid: e.N_Fids) {
                auto& f = this->mesh->F.at(fid);
                int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
                if (getFace_V(f, idx, 1) == prev) {
                    if (getDiagonals) res.push_back(getFace_V(f, idx, 2));
                    eid = this->mu->GetDifference(this->mu->GetIntersection(v.N_Eids, f.Eids), std::vector<size_t>{e.id}).at(0);
                    break;
                }
            }
        }
        return res;
    };

    // lambda function to determine type of terminal vertex of a path

    auto terminalVertexType = [&] (size_t initial_, size_t terminal_, bool& same, bool& diff, bool& boundary) {
        auto& initial = this->mesh->V.at(initial_);
        auto& terminal = this->mesh->V.at(terminal_);
        int initital_valence = initial.N_Vids.size();
        int terminal_valence = terminal.N_Vids.size();
        if (initital_valence == 3) {
            if (terminal.isBoundary) {
                if (terminal.N_Vids.size() == 4 || (terminal.N_Vids.size() == 3 && terminal.idealValence == 1)) diff = true;
                if (terminal.N_Vids.size() == 2 && terminal.idealValence != 1) same = true;
            } else {
                if (terminal.N_Vids.size() == 5) diff = true;
                if (terminal.N_Vids.size() == 3) same = true;
            }
        } else if (initital_valence == 5) {
            if (terminal.isBoundary) {
                if (terminal.N_Vids.size() == 2 && terminal.idealValence != 1) diff = true;
                if (terminal.N_Vids.size() == 4 || (terminal.N_Vids.size() == 3 && terminal.idealValence == 1)) same = true;
            } else {
                if (terminal.N_Vids.size() == 3) diff = true;
                if (terminal.N_Vids.size() == 5) same = true;
            }
        }
        if (terminal.isBoundary || terminal.type == FEATURE) boundary = true;
    };

    // lambda function to determine singularity
    auto isSingularity = [&] (Vertex& v) {
        bool res = false;
        if (v.N_Vids.size() == 3 || v.N_Vids.size() == 5) {
            if (!v.isBoundary && v.type != FEATURE) res = true;
        }
        // if (v.isBoundary && (v.N_Vids.size() == 4 || (!this->mu->IsSharpFeature(v.id) && v.N_Vids.size() == 2))) res = true;
        if (v.isBoundary && v.N_Vids.size() != v.idealValence+1) res = true;
        return res;
    };

    // lambda function to determine non essential irregular vertex

    auto isIrregular = [&] (Vertex& v) {
        bool res = false;
        if (!v.isBoundary) {
            if (v.type != FEATURE && (v.N_Vids.size() > 5 || v.N_Vids.size() < 3)) res = true;
            if (this->mu->IsSharpFeature(v.id)) res = true;
        }
        if (v.isBoundary) {
            // if (v.N_Vids.size() > 4 || this->mu->IsSharpFeature(v.id)) res = true;
            if ((v.idealValence == 2 && (v.N_Vids.size() > 4 || v.N_Vids.size() < 2))) res = true;
        }
        return res;
    };

    // lambda function to set vertex QEMs
    auto setvQEMs = [&] (std::vector<vQEM>& vQEMs) {
        vQEMs.resize(this->mesh->V.size());
        PARALLEL_FOR_BEGIN(0, this->mesh->V.size()) {
            double qem = this->mu->CalculateQEM(i);
            {
                std::lock_guard<std::mutex> lock(mx);
                vQEMs.at(i).qem = qem;
            }
        } PARALLEL_FOR_END();
        
        PARALLEL_FOR_BEGIN(0, this->mesh->V.size()) {
            auto& v = this->mesh->V.at(i);
            std::unordered_set<size_t> vids;
            vids.insert(v.N_Vids.begin(), v.N_Vids.end());
            for (auto vid: v.N_Vids) {
                auto& nv = this->mesh->V.at(vid);
                vids.insert(nv.N_Vids.begin(), nv.N_Vids.end());
            }
            double sum_qem = 0, mean_qem = 0.0, std_qem = 0.0;
            for (auto vid: vids) {
                sum_qem += vQEMs.at(vid).qem;
            }
            mean_qem = sum_qem / vids.size();
            for (auto vid: vids) {
                std_qem += std::pow(vQEMs.at(vid).qem - mean_qem, 2);
            }
            std_qem = std::sqrt(std_qem / vids.size());
            {
                std::lock_guard<std::mutex> lock(mx);
                vQEMs.at(i).lower_bound = mean_qem - std_qem;
                vQEMs.at(i).upper_bound = mean_qem + std_qem;
                double val = vQEMs.at(i).upper_bound - vQEMs.at(i).qem;
                vQEMs.at(i).qem = val < 0 ? exp(val) - 1.0 : 1.0 - exp(-val);
            }
        } PARALLEL_FOR_END();
    };

    // lambda function to get vertex QEMs
    auto getVQEMs = [&] (Vertex& v) {
        std::unordered_set<size_t> vids;
        vids.insert(v.N_Vids.begin(), v.N_Vids.end());
        for (auto vid: v.N_Vids) {
            auto& nv = this->mesh->V.at(vid);
            vids.insert(nv.N_Vids.begin(), nv.N_Vids.end());
        }
        std::vector<double> qems;
        for (auto vid: vids) {
            qems.push_back(this->mu->CalculateQEM(vid));
        }
        double sum_qem = 0, mean_qem = 0.0, std_qem = 0.0;
        for (auto qem: qems) {
            sum_qem += qem;
        }
        mean_qem = sum_qem / qems.size();
        for (auto qem: qems) {
            std_qem += std::pow(qem - mean_qem, 2);
        }
        std_qem = std::sqrt(std_qem / qems.size());
        vQEM vqem;
        // vqem.qem = this->mu->CalculateQEM(v.id);
        vqem.lower_bound = mean_qem - std_qem;
        vqem.upper_bound = mean_qem + std_qem;
        double val = vqem.upper_bound - this->mu->CalculateQEM(v.id);
        vqem.qem = val < 0 ? exp(val) - 1.0 : 1.0 - exp(-val);
        return vqem;
    };

    // lambda function to get singularity score
    auto getSingularityScore = [&] (size_t vid) {
        auto& v = this->mesh->V.at(vid);
        double score = 0;
        std::vector<size_t> neighbors = getVertexNeighbors(v.id, v.N_Vids.at(0));
        // double boundary_score = 0;
        // int boundary_min_dist = std::numeric_limits<int>::max();
        // std::cout << "Getting singularity score for: " << v.id << " " << v.N_Vids.size() << std::endl;
        if (v.isBoundary || v.type == FEATURE) score -= exp(2.0);
        std::priority_queue<double> same_score_dist, diff_score_dist, boundary_score_dist;
        std::vector<bool> visited(this->mesh->V.size(), false);
        for (auto nid: neighbors) {
            path_v pv; pv.curr = nid; pv.prev = v.id; pv.dist = 1; pv.continuous = true;
            std::queue<path_v> q;
            q.push(pv);
            int same_count = 0, diff_count = 0;
            while (!q.empty()) {
                auto pv_= q.front();
                q.pop();
                auto& v_ = this->mesh->V.at(pv_.curr);
                if (visited.at(v_.id) || isIrregular(v_)) {
                    visited.at(v_.id) = true;
                    continue;
                } 
                if (v_.N_Vids.size() != 4 || v_.isBoundary || v_.type == FEATURE) {
                    // std::cout << "encountered vertex: " << v_.N_Vids.size() << " boundary: " << v_.isBoundary << std::endl;
                    visited.at(v_.id) = true;
                    bool same = false; bool diff = false; bool boundary = false;
                    terminalVertexType(v.id, v_.id, same, diff, boundary);
                    double diag_dist = (pv_.dist_b > 0) ? std::max(pv_.dist, pv_.dist_b) : -1;
                    // if (pv_.dist_b > 0) {
                    //     diag_dist = (pv_.dist + pv_.dist_b)%2 == 0 ? std::max(pv_.dist, pv_.dist_b) : std::min(pv_.dist, pv_.dist_b) + fabs(pv_.dist - pv_.dist_b);
                    // }
                    double dir_dist = pv_.dist + pv_.dist_b;
                    if (same) {
                        // score_values.push(-exp(1-pv_.dist));
                        // score -= exp(1-pv_.dist);
                        // same_score_dist.push(((dir_dist-4)-fabs(dir_dist-4)));
                        same_score_dist.push(((4-dir_dist)+fabs(4-dir_dist)));
                    }
                    if (diff) {
                        // score_values.push(exp(1-pv_.dist));
                        // score += exp(1-pv_.dist);
                        // diff_score_dist.push((1.0/(double)pv_.dist));
                        // std::cout << "encountered opposit valence" << std::endl;
                        // std::cout << "diag_dist: " << diag_dist << " dir_dist: " << dir_dist << std::endl;
                        diag_dist > 0 && diag_dist < dir_dist ? diff_score_dist.push((1.0/diag_dist)) : diff_score_dist.push((1.0/dir_dist));
                    }
                    if (boundary) {
                        // boundary_score_dist.push(((dir_dist-4)-fabs(dir_dist-4)));
                        boundary_score_dist.push(((4-dir_dist)+fabs(4-dir_dist)));
                        // if (pv_.dist < boundary_min_dist) {
                        //     boundary_min_dist = pv_.dist;
                        //     boundary_score = exp(1-pv_.dist);
                        // }
                        // break;
                    }
                    continue;
                }
                std::vector<size_t> neighbors_ = getVertexNeighbors(pv_.curr, pv_.prev);
                // std::cout << "v: " << v_.id << " neighbors(" << neighbors_.size() << "): " << neighbors_.at(1) << " " << neighbors_.at(3) << " " << neighbors_.at(2) << std::endl;
                path_v pv1, pv2, pv3;
                pv1.curr = neighbors_.at(1); pv2.curr = neighbors_.at(2); pv3.curr = neighbors_.at(3);
                pv1.prev = pv_.curr; pv2.prev = pv_.curr; pv3.prev = pv_.curr;
                // pv1.dist = pv_.dist+1; pv2.dist = pv_.dist+1; pv3.dist = pv_.dist+1;
                if (pv_.continuous) {
                    visited.at(v_.id) = true;
                    pv1.dist = pv_.dist; pv1.dist_b = 1;
                    pv2.dist = pv_.dist+1;
                    pv3.dist = pv_.dist; pv3.dist_b = 1;
                    pv2.continuous = true;
                    q.push(pv1); q.push(pv3); q.push(pv2);
                } else {
                    pv2.dist = pv_.dist; pv2.dist_b = pv_.dist_b+1;
                    q.push(pv2);
                }
                
                // visited.at(v_.id) = true;
            }
        }
        int it = 0;
        // std::cout << "same_score_dist: " << same_score_dist.size() << " diff_score_dist: " << diff_score_dist.size() << " boundary_score_dist: " << boundary_score_dist.size() << std::endl;
        while (it++ < 3) {
            if (!same_score_dist.empty()) {
                score -= (1 - exp(-same_score_dist.top())); // 1 - e^((x-4)-|x-4|)
                // std::cout << "same score dist: " << same_score_dist.top() << std::endl;
        // std::cout << "same score decrease: -" << (1 - exp(-same_score_dist.top())) << std::endl;
                same_score_dist.pop();
            }
            if (!diff_score_dist.empty()) {
                // std::cout << "diff_score_dist: " << diff_score_dist.top() << std::endl;
                score += exp(diff_score_dist.top()); // e^(1/x)
        // std::cout << "diff score increase: +" << exp(diff_score_dist.top()) << std::endl;
                diff_score_dist.pop();
            }
        }
        if (!boundary_score_dist.empty()) {
            // std::cout << "boundary score dist: " << boundary_score_dist.top() << std::endl;
            score -= (1 - exp(-boundary_score_dist.top())); // 1 - e^((x-4)-|x-4|)
        // std::cout << "boundary score decrease : -" << (1 - exp(-boundary_score_dist.top())) << std::endl;
        }
        // std::cout << "score: " << score << std::endl;
        return score;
    };

    // lambda function to get valence score
    auto getValenceScore = [&] (Vertex& v) {
        if (v.N_Vids.empty()) return 0.0;
        int ideal_valence = 4;
        if (v.isBoundary || v.type == FEATURE) {
            this->mesh->SetIdealValence(v.id);
            ideal_valence = v.idealValence + 1;
        }
        // std::cout << "N_vids: " << v.N_Vids.size() << std::endl;
        // std::cout << "ideal_valence: " << ideal_valence << std::endl;
        // std::cout << "isBoundary: " << v.isBoundary << std::endl;
        // std::cout << "difference: " << (int)v.N_Vids.size() - (int)ideal_valence << std::endl;
        // return 1.0 - exp(fabs((int)v.N_Vids.size() - (int)ideal_valence));
        return fabs((int)v.N_Vids.size() - (int)ideal_valence);
    };

    // lambda function to determine valid path
    auto isPathValid = [&] (Path& p) {
        // std::cout << "path vids: " << p.vids.size() << std::endl;
        // this->PrototypeSaveSeparatrices(std::vector<std::vector<size_t>>{p.vids}, "test");
        // std::lock_guard<std::mutex> lock(this->mx);
        struct pScore {
            // double singularity_score = 0.0;
            // double valence_score = 0.0;
            double qem_score = 0.0;
        };
        bool res = false;
        // std::unordered_set<size_t> vids(p.vids.begin(), p.vids.end());
        double singularity_score_before = 0.0;
        double valence_score_before = 0.0;
        // double qem_score_before = 0.0;
        // int n = 0;
        // this->PrototypeSaveSeparatrices(std::vector<std::vector<size_t>>{p.vids}, "test_before");
        // std::cout << "Calculating scores before..." << std::endl;
        // std::cout << "QEMs size: " << vQEMs.size() << std::endl;
        std::unordered_set<size_t> pVids;
        for (auto vid: p.vids) {
            // pVids.insert(vid);
            auto& pv = this->mesh->V.at(vid);
            for (auto fid: pv.N_Fids) {
                auto& pf = this->mesh->F.at(fid);
                pVids.insert(pf.Vids.begin(), pf.Vids.end());
            }
            // pVids.insert(pv.N_Vids.begin(), pv.N_Vids.end());
            // for (auto nvid: pv.N_Vids) {
                // auto& nv = this->mesh->V.at(nvid);
                // pVids.insert(nv.N_Vids.begin(), nv.N_Vids.end());
            // }
        }
        std::map<size_t, pScore> vid_scores;
        int n_singularities_before = 0;
        for (auto vid: pVids) {
            auto& v = this->mesh->V.at(vid);
            // std::cout << "vertex: " << v.id << "(" << v.N_Vids.size() << ") ";
            pScore ps;
            if (isSingularity(v)) {
                singularity_score_before += getSingularityScore(vid);
                ++n_singularities_before;
            }
            // ps.valence_score = getValenceScore(v);
            ps.qem_score = vQEMs.at(vid).qem;
            vid_scores.insert(std::make_pair(vid, ps));
            // if (vid < vQEMs.size() && !this->mesh->V.at(vid).N_Vids.empty()) {
            //     // ++n;
            //     // qem_score_curr = (((n-1)*qem_score_curr) + vQEMs.at(vid)) / n;
            //     double qem = vQEMs.at(vid).upper_bound - vQEMs.at(vid).qem;
            //     std::cout << "qem upper bound: " << vQEMs.at(vid).upper_bound << " v qem: " << vQEMs.at(vid).qem << " ";
            //     double score = qem < 0 ? exp(qem) - 1.0 : 1.0 - exp(-qem);
            //     qem_score_before += score;
            // }
            valence_score_before += getValenceScore(v);
            // std::cout << valence_score_before << " " << singularity_score_before << " " << qem_score_before << std::endl;
        }
        n_singularities_before > 0 ? singularity_score_before /= n_singularities_before : 1;
        // double score_before = singularity_score_before + qem_score_before + valence_score_before;
        // std::cout << "valence_score_before: " << valence_score_before << " singularity_score_before: " << singularity_score_before << std::endl;
        double best_score = 0.0;
        int best_dir = -1;
        std::vector<size_t> dirs = getVertexNeighbors(p.vids.at(0), p.vids.at(1), true);
        this->mu->UpdateContents(dirs, std::vector<size_t>{p.vids.at(1)});
        int it = 0;
        for (auto dir: dirs) {
            {
                // std::lock_guard<std::mutex> lock(this->mx);
                
                // std::cout << "DIRECTION ITERATION: " << it++ << std::endl;
                this->caretaker.saveState(this->mesh);
                double singularity_score_after = std::numeric_limits<double>::min();
                double valence_score_after = 0.0;
                double qem_score_after = 0.0;
                p.move_dir = dir;
                // try {
                    this->PrototypeMoveAlongPath(p);

                // } catch (const std::exception& e) {
                    // singularity_score_after = -exp(5.0);
                    // this->caretaker.restoreState(this->mesh);
                    // continue;
                // }
                // this->PrototypeSaveSeparatrices(std::vector<std::vector<size_t>>{p.vids}, "test_after_"+std::to_string(it));
                // std::cout << "After moving along path: " << std::endl;
                // n = 0;
                // std::cout << "\n****************\nCalculating scores after..." << std::endl;
                int prev_vSize = this->caretaker.vSize();
                int curr_vSize = this->mesh->V.size();
                int n_singularities_after = 0;
                if (prev_vSize < curr_vSize) {
                    int diff = curr_vSize - prev_vSize;
                    for (int i = 0; i < diff; i++) {
                        auto& v = this->mesh->V.at((size_t)(prev_vSize + i));
                        if (isSingularity(v)) {
                            singularity_score_after += getSingularityScore(v.id);
                            // std::cout << "singularity score after: " << singularity_score_after << std::endl;
                            ++n_singularities_after;
                        }
                        double curr_valence_score = getValenceScore(v);
                        valence_score_after += curr_valence_score;
                    }
                }
                for (auto vid: pVids) {
                    auto& v = this->mesh->V.at(vid);
                    if (v.N_Vids.empty()) continue;
                    // std::cout << v.N_Vids.size() << " ";
                    if (isSingularity(v)) {
                        singularity_score_after += getSingularityScore(vid);
                        // std::cout << "singularity score after: " << singularity_score_after << std::endl;
                        ++n_singularities_after;
                    }
                    double curr_valence_score = getValenceScore(v);
                    // valence_score_before += vid_scores.at(vid).valence_score;
                    valence_score_after += curr_valence_score;
                    // std::cout << "valence: " << v.N_Vids.size() <<  " curr valence score: " << valence_score_after << " valence score before: " << vid_scores.at(vid).valence_score << std::endl;
                    double val = vQEMs.at(vid).upper_bound - this->mu->CalculateQEM(vid);
                    qem_score_after = val < 0 && !(qem_score_after < 0) ? -1.0 : 1.0;
                    // double qem = val < 0 ? exp(val) - 1.0 : 1.0 - exp(-val);
                    // qem_score_after += (qem - vid_scores.at(vid).qem_score);
                    // qem_score_after += qem;
                    // std::cout << "-----------------" << std::endl;
                }
                n_singularities_after > 0 ? singularity_score_after /= n_singularities_after : 1;
                // std::cout << "n_singularities_after: " << n_singularities_after << " " << singularity_score_after << std::endl;
                // std::cout << std::endl;
                // std::cout << "\n****************\n" << std::endl;
                // for (auto vid: p.vids) {
                //     auto& v = this->mesh->V.at(vid);
                //     std::cout << "vertex: " << v.id << "(" << v.N_Vids.size() << ") ";
                //     if (isSingularity(v)) singularity_score_after += getSingularityScore(vid);
                //     if (vid < vQEMs.size() && !this->mesh->V.at(vid).N_Vids.empty()) {
                //         // ++n;
                //         double qem = vQEMs.at(vid).upper_bound - this->mu->CalculateQEM(vid);
                //         std::cout << "qem upper bound: " << vQEMs.at(vid).upper_bound << " v qem: " << this->mu->CalculateQEM(vid) << " ";
                //         double score = qem > 0 ? 1.0 - exp(-qem) : exp(qem) - 1.0;
                //         // qem_score_curr = (((n-1)*qem_score_curr) + score) / n;
                //         qem_score_after += score;
                //     }
                //     if (!v.N_Vids.empty()) valence_score_after += getValenceScore(v);
                //     std::cout << valence_score_after << " " << singularity_score_after << " " << qem_score_after << std::endl;
                // }
                // this->caretaker.restoreState(this->mesh);
                // double score_after = (singularity_score_after - singularity_score_before) + qem_score_after + valence_score_after;
                // std::cout << "valence_score_after: " << valence_score_after << " valence_score_before: " << valence_score_before << std::endl;
                // std::cout << "singularity_score_after: " << singularity_score_after << " singularity_score_before: " << singularity_score_before << std::endl;
                double sscore = singularity_score_after == std::numeric_limits<double>::min() ? 1 : singularity_score_after - singularity_score_before;
                // double sscore = singularity_score_after - singularity_score_before;
                double vscore = 2.0 - exp(valence_score_after - valence_score_before); // 2-e^(v_after-v_before)
                double score_after = sscore < 0 && vscore < 0 ? -(sscore*vscore) : (sscore*vscore);
                score_after = score_after < 0 && qem_score_after < 0 ? -(qem_score_after*score_after) : (qem_score_after*score_after);
                // std::cout << "sscore: " << sscore << "   vscore: " << vscore << " score_after: " << score_after << std::endl;
                // double score_after = () - fabs(valence_score_after - valence_score_before);
                // std::cout << "singularity score: before: " << singularity_score_before << " singularity score: after: " << singularity_score_after << " valence_score_diff: " << -fabs(valence_score_after - valence_score_before) << std::endl;
                // std::cout << "score before: " << score_before << " score after: " << score_after << std::endl;
                if (score_after > best_score) {
                    // std::cout << "DIRECTION ITERATION: " << it++ << std::endl;
                    // this->PrototypeSaveSeparatrices(std::vector<std::vector<size_t>>{p.vids}, "test_after_"+std::to_string(it));
                    // std::cout << "valence_score_after: " << valence_score_after << " valence_score_before: " << valence_score_before << std::endl;
                    // std::cout << "singularity_score_after: " << singularity_score_after << " singularity_score_before: " << singularity_score_before << std::endl;
                    // std::cout << "sscore: " << sscore << " vscore: " << vscore << " score_after: " << score_after << std::endl;
                    // std::cout << "*******************************" << std::endl;
                    
                    best_score = score_after;
                    // std::cout << "BEST SCORE: " << best_score << std::endl;
                    best_dir = dir;
                }
                this->caretaker.restoreState(this->mesh);
                // if (it == iters) break;
                // std::cout << "*******************************" << std::endl;
            }
            // std::cout << "Restored mesh" << std::endl;
        }
        if (best_score > 0 && best_dir != -1) {
            res = true;
            p.move_dir = best_dir;
            p.score = best_score;
        }
        return res;
    };

    auto getPath = [&] (path_v& pv, std::priority_queue<Path, std::vector<Path>, p_comp>& paths) {
        std::queue<path_v> q;
        q.push(pv);
        while (!q.empty()) {
            auto pv_= q.front();
            // std::cout << "pv_ curr: " << pv_.curr << " pv_ prev: " << pv_.prev << std::endl;
            q.pop();
            auto& v_ = this->mesh->V.at(pv_.curr);
            // std::cout << "vertex: " << v_.id << "(" << v_.N_Vids.size() << ") ";
            if (isIrregular(v_)) continue;
            if ((!v_.isBoundary && v_.type != FEATURE && (v_.N_Vids.size() == 3 || v_.N_Vids.size() == 5)) ||
                (v_.isBoundary && (v_.N_Vids.size() == 3 || (!this->mu->IsSharpFeature(v_.id) && (v_.N_Vids.size() == 4 || v_.N_Vids.size() == 2)))) || 
                v_.id == pv_.vids.front()) {
                Path p;
                p.vids = pv_.vids;
                if (isPathValid(p)) {
                    // {
                        // std::lock_guard<std::mutex> lock(this->mx);
                        paths.push(p);
                    // }
                }
                continue;
            }
            std::vector<size_t> neighbors_ = getVertexNeighbors(pv_.curr, pv_.prev);
            // std::cout << "neighbors: " << neighbors_.size() << std::endl;
            if (neighbors_.size() != 4) continue;
            path_v pv1, pv2, pv3;
            pv1.prev = pv_.curr; pv2.prev = pv_.curr; pv3.prev = pv_.curr;
            pv1.curr = neighbors_.at(1); pv2.curr = neighbors_.at(2); pv3.curr = neighbors_.at(3);
            pv1.vids = pv_.vids; pv2.vids = pv_.vids; pv3.vids = pv_.vids;
            pv1.vids.push_back(pv1.curr); pv2.vids.push_back(pv2.curr); pv3.vids.push_back(pv3.curr);
            if (pv_.continuous) {
                pv2.continuous = true;
                if (!this->mesh->V.at(pv1.curr).isBoundary && this->mesh->V.at(pv1.curr).type != FEATURE) q.push(pv1);
                if (!this->mesh->V.at(pv3.curr).isBoundary && this->mesh->V.at(pv3.curr).type != FEATURE) q.push(pv3);
            }
            q.push(pv2);
        }
        // std::cout << std::endl;
    };

    auto getPaths = [&] (path_v& pv, std::vector<Path>& paths) {
        std::queue<path_v> q;
        q.push(pv);
        while (!q.empty()) {
            auto pv_= q.front();
            // std::cout << "pv_ curr: " << pv_.curr << " pv_ prev: " << pv_.prev << std::endl;
            q.pop();
            auto& v_ = this->mesh->V.at(pv_.curr);
            // std::cout << "vertex: " << v_.id << "(" << v_.N_Vids.size() << ") ";
            if (isIrregular(v_)) continue;
            if (v_.N_Vids.size() != 4 || v_.isBoundary) {
            // if ((!v_.isBoundary && v_.type != FEATURE && (v_.N_Vids.size() == 3 || v_.N_Vids.size() == 5)) ||
            //     (v_.isBoundary && (v_.N_Vids.size() == 3 || (!this->mu->IsSharpFeature(v_.id) && (v_.N_Vids.size() == 4 || v_.N_Vids.size() == 2)))) || 
            //     v_.id == pv_.vids.front()) {
                Path p;
                p.vids = pv_.vids;
                // if (isPathValid(p)) {
                    // {
                        // std::lock_guard<std::mutex> lock(this->mx);
                        paths.push_back(p);
                    // }
                // }
                continue;
            }
            std::vector<size_t> neighbors_ = getVertexNeighbors(pv_.curr, pv_.prev);
            // std::cout << "neighbors: " << neighbors_.size() << std::endl;
            if (neighbors_.size() != 4) continue;
            path_v pv1, pv2, pv3;
            pv1.prev = pv_.curr; pv2.prev = pv_.curr; pv3.prev = pv_.curr;
            pv1.curr = neighbors_.at(1); pv2.curr = neighbors_.at(2); pv3.curr = neighbors_.at(3);
            pv1.vids = pv_.vids; pv2.vids = pv_.vids; pv3.vids = pv_.vids;
            pv1.vids.push_back(pv1.curr); pv2.vids.push_back(pv2.curr); pv3.vids.push_back(pv3.curr);
            if (pv_.continuous) {
                pv2.continuous = true;
                if (!this->mesh->V.at(pv1.curr).isBoundary && this->mesh->V.at(pv1.curr).type != FEATURE) q.push(pv1);
                if (!this->mesh->V.at(pv3.curr).isBoundary && this->mesh->V.at(pv3.curr).type != FEATURE) q.push(pv3);
            }
            q.push(pv2);
        }
        // std::cout << std::endl;
    };

    // lambda function to get exclusive paths
    auto isPathExclusive = [&] (Path& p, std::vector<bool>& visited) {
        bool res = true;
        std::unordered_set<size_t> vids;
        for (auto vid: p.vids) {
            auto& pv = this->mesh->V.at(vid);
            for (auto fid: pv.N_Fids) {
                auto& pf = this->mesh->F.at(fid);
                vids.insert(pf.Vids.begin(), pf.Vids.end());
            }
        }
        for (auto vid: vids) {
            if (visited.at(vid)) {
                res = false;
                break;
            }
        }
        if (res) {
            for (auto vid: vids) visited.at(vid) = true;
        }
        return res;
    };

    std::vector<path_v> pvs;
    // PARALLEL_FOR_BEGIN(0, mesh->V.size()) {
    for (int i = 0; i < mesh->V.size(); i++) {
        auto& v = mesh->V.at(i);
        if ((v.N_Vids.size() == 3 || v.N_Vids.size() == 5) && (!v.isBoundary && v.type != FEATURE)) {
        // if ((v.N_Vids.size() == 3) && (!v.isBoundary && v.type != FEATURE)) {
            std::vector<size_t> neighbors = getVertexNeighbors(i, v.N_Vids.at(0));
            for (auto nid: neighbors) {
                path_v pv; pv.curr = nid; pv.prev = v.id; pv.vids = {v.id, nid}; pv.continuous = true;
                {
                    std::lock_guard<std::mutex> lock(mx);
                    pvs.push_back(pv);
                }
            }
        }
    }
    // } PARALLEL_FOR_END();
    
    std::cout << "path initiators: " << pvs.size() << std::endl;

    // PARALLEL_FOR_BEGIN(0, pvs.size()) {
    // setvQEMs(vQEMs);
    // std::priority_queue<Path, std::vector<Path>, p_comp> fPaths;
    for (int i = 0; i < pvs.size(); i++) {
        // std::cout << "GETTING Path: " << i << std::endl;
        getPaths(pvs.at(i), paths);
        // getPath(pvs.at(i), fPaths);
        // getPaths(pvs.at(i), paths);
    }
    // } PARALLEL_FOR_END();

    std::cout << "paths: " << paths.size() << std::endl;
    
    // std::vector<std::vector<size_t>> pv;
    // for (auto p: paths) {
    //     pv.push_back(p.vids);
    // }
    // pv.push_back(paths.at(1).vids);
    // PrototypeSaveSeparatrices(pv, "paths");

    // std::vector<Path> fPaths;
    std::priority_queue<Path, std::vector<Path>, p_comp> fPaths;
    setvQEMs(vQEMs);
    // PARALLEL_FOR_BEGIN(0, paths.size()) {
    for (int i = 0; i < paths.size(); i++) {
    // for (int i = iters; i < iters+1; i++) {
    // for (int i = 0; i < 1; i++) {
        std::cout << "CHECKING PATH: " << i << std::endl;
        
        PrototypeSaveSeparatrices(std::vector<std::vector<size_t>>{paths.at(i).vids}, "saving_path_"+std::to_string(i));
        bool res = false;
        // {
            // std::lock_guard<std::mutex> lock(mx);
            res = isPathValid(paths.at(i));
        // }
        if (res) {
            // std::lock_guard<std::mutex> lock(mx);
            // fPaths.push_back(paths.at(i));
            fPaths.push(paths.at(i));
        }
        // break;
    }
    // } PARALLEL_FOR_END();

    // PrototypeMoveAlongPath(fPaths.at(0));
    std::cout << "Final paths: " << fPaths.size() << std::endl;
    int it = 0;
    std::priority_queue<Path, std::vector<Path>, p_comp> finalPaths;
    std::vector<bool> visited(mesh->V.size(), false);
    while (!fPaths.empty()) {
        auto& p = (Path&) fPaths.top();
        std::cout << "score: " << p.score << " length: " << p.vids.size()-1 << std::endl;
        // caretaker.saveState(mesh);
        // PrototypeSaveSeparatrices(std::vector<std::vector<size_t>>{p.vids}, "saving_path_"+std::to_string(++it));
        // std::cout << "it: " << it << std::endl;
        // PrototypeMoveAlongPath(p);
        // PrototypeSaveSeparatrices(std::vector<std::vector<size_t>>{}, "executed_path_"+std::to_string(it));
        // caretaker.restoreState(mesh); 
        if (isPathExclusive(p, visited)) {
            PrototypeSaveSeparatrices(std::vector<std::vector<size_t>>{p.vids}, "saving_path_"+std::to_string(++it));
            finalPaths.push(p);
        }
        fPaths.pop();
        // break;
    }
    std::cout << "Final exclusive paths: " << finalPaths.size() << std::endl;
    it = 0;
    while (!finalPaths.empty()) {
        auto& p = (Path&) finalPaths.top();
        std::cout << "score: " << p.score << std::endl;
        // caretaker.saveState(mesh);
        // PrototypeSaveSeparatrices(std::vector<std::vector<size_t>>{p.vids}, "saving_path_"+std::to_string(++it));
        // std::cout << "it: " << it << std::endl;
        
        PrototypeMoveAlongPath(p);
        RenderMesh();
        PrototypeSaveSeparatrices(std::vector<std::vector<size_t>>{}, "executed_path_"+std::to_string(++it));
        // caretaker.restoreState(mesh); 
        finalPaths.pop();
        // break;
    }
    return paths;
}

// PrototypeGetPathsAt(vid, paths, singularities) that performs Depth First Search to find multiple paths between vid and singularities
// it also takes into account the diagonal vertices in faces during search

void SemiGlobalSimplifier::PrototypeGetPathsAt(size_t vid, std::vector<Path>& paths) {
    
    // lambda function to set path parameters
    auto setPathParams = [&] (Path& p, int b1, int b2, int rots, bool intersect_b1, bool is_diagonal, std::vector<size_t>& sec_vids) {
        if (p.vids.empty()) return false;
        if (this->mu->IsSharpFeature(p.vids.back())) return false;
        p.b1 = b1;
        p.b2 = b2;
        p.rots = rots;
        p.intersect_b1 = intersect_b1;
        p.is_diagonal = is_diagonal;
        if (b2 > 0) p.vids.insert(p.vids.begin(), sec_vids.begin(), sec_vids.begin()+b1);
        return true;
    };

    // lambda function to get vertex neighboring vertices in order

    auto getVertexNeighbors = [&] (size_t vid, size_t start) -> std::vector<size_t> {
        auto& v = this->mesh->V.at(vid);
        std::vector<size_t> res;
        size_t eid = this->GetEdgeId(v.id, start);
        auto getFace_V = [] (Face& f, int idx, int offset) -> size_t {
            return f.Vids.at((idx+offset)%f.Vids.size());
        };
        for (int i = 0; i < v.N_Eids.size(); i++) {
            auto& e = mesh->E.at(eid);
            size_t prev = e.Vids.at(0) == v.id ? e.Vids.at(1) : e.Vids.at(0);
            res.push_back(prev);
            for (auto fid: e.N_Fids) {
                auto& f = this->mesh->F.at(fid);
                int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
                if (getFace_V(f, idx, 1) == prev) {
                    eid = this->mu->GetDifference(this->mu->GetIntersection(v.N_Eids, f.Eids), std::vector<size_t>{e.id}).at(0);
                    break;
                }
            }
        }
        return res;
    };

    auto getIndirectNeighbors = [&] (size_t vid, size_t start) -> std::vector<size_t> {
        auto& v = this->mesh->V.at(vid);
        std::vector<size_t> res;
        size_t fid = this->GetFaceId(v.id, start);
        auto getFace_V = [] (Face& f, int idx, int offset) -> size_t {
            return f.Vids.at((idx+offset)%f.Vids.size());
        };
        for (int i = 0; i < v.N_Fids.size(); i++) {
            auto& f = mesh->F.at(fid);
            int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
            size_t prev = getFace_V(f, idx, 3);
            res.push_back(getFace_V(f, idx, 2));
            for (auto nfid: v.N_Fids) {
                auto& nf = this->mesh->F.at(nfid);
                idx = std::distance(nf.Vids.begin(), std::find(nf.Vids.begin(), nf.Vids.end(), v.id));
                if (getFace_V(f, idx, 1) == prev) {
                    fid = nfid;
                    break;
                }
            }
        }
        return res;
    };

    // lambda function to determine type of terminal vertex of a path

    auto terminalVertexType = [&] (size_t initial_, size_t terminal_, bool& same, bool& diff, bool& boundary) {
        auto& initial = this->mesh->V.at(initial_);
        auto& terminal = this->mesh->V.at(terminal_);
        int initital_valence = initial.N_Vids.size();
        int terminal_valence = terminal.N_Vids.size();
        if (initital_valence == 3) {
            if (terminal.isBoundary) {
                if (terminal.N_Vids.size() == 4) diff = true;
                if (terminal.N_Vids.size() == 2 && !this->mu->IsSharpFeature(terminal.id)) same = true;
            } else {
                if (terminal.N_Vids.size() == 5) diff = true;
                if (terminal.N_Vids.size() == 3) same = true;
            }
        } else if (initital_valence == 5) {
            if (terminal.isBoundary) {
                if (terminal.N_Vids.size() == 2 && !this->mu->IsSharpFeature(terminal.id)) diff = true;
                if (terminal.N_Vids.size() == 4) same = true;
            } else {
                if (terminal.N_Vids.size() == 3) diff = true;
                if (terminal.N_Vids.size() == 5) same = true;
            }
        }
        if (terminal.isBoundary || terminal.type == FEATURE) boundary = true;
    };

    // lambda function to determine singularity
    auto isSingularity = [&] (size_t vid) {
        bool res = false;
        auto& v = this->mesh->V.at(vid);
        if (v.N_Vids.size() == 3 || v.N_Vids.size() == 5) {
            if (!v.isBoundary && v.type != FEATURE) res = true;
        }
        if (v.isBoundary && (v.N_Vids.size() == 4 || (!mu->IsSharpFeature(vid) && v.N_Vids.size() == 2))) res = true;
        return res;
    };

    // lambda function to set singularity score

    auto getPaths = [&] (size_t vid, bool incluePaths,  std::vector<Path>& paths) {
        // std::cout << "Inside get Paths" << std::endl;
        struct path_v {
            size_t id;
            std::vector<size_t> vids;
            bool continuous = false;
            int rots = 0;
            int b1 = 0;
            int b2 = 0;
            bool intersect_b1 = false;
            bool intersect_b2 = false;
            int boundary_distance = 0;
        };
        auto& v = this->mesh->V.at(vid);
        double score = 0.0;
        if ((v.N_Vids.size() != 3 && v.N_Vids.size() != 5) || v.isBoundary) return score;
        auto v_score = this->mu->GetVertexScore(v.id);
        int same_count = 0; int diff_count = 0; int same_singularity_min_dist = 0;
        int diff_singularity_max_dist = std::numeric_limits<int>::max(); int boundary_min_dist = 0;
        double same_singularity_score = 0.0; double diff_singularity_score = 0.0; double boundary_score = 0.0;
        std::vector<size_t> neighbors = getVertexNeighbors(v.id, v.N_Vids.at(0));
        // std::cout << "neighbors: " << neighbors.size() << std::endl;
        int i = 0;
        for (auto nid: neighbors) {
            // std::cout << "Start: " << v.id << "(" << v.N_Vids.size() << ") -> ";
            // std::cout << "nid: " << nid << std::endl;
            path_v pv;
            pv.id = nid;
            pv.vids = {v.id, nid};
            pv.continuous = true;
            pv.b1 = 1;
            std::queue<path_v> q;
            q.push(pv);
            while (!q.empty()) {
                // std::cout << "Queue size: " << q.size() << std::endl;
                auto pv_= q.front();
                q.pop();
                auto& v_ = this->mesh->V.at(pv_.id);
                if (v_.isBoundary || v_.type == FEATURE) {
                    if (pv_.continuous && !pv_.intersect_b1) {
                        pv_.intersect_b1 = true; pv_.boundary_distance = pv_.b1;
                    } else if (!pv.continuous && !pv.intersect_b1 && !pv.intersect_b2) {
                        pv_.intersect_b2 = true; pv_.boundary_distance = pv_.b1 + pv_.b2;
                    }
                    pv_.continuous ? pv_.intersect_b1 = true : pv_.intersect_b2 = true;

                }
                // std::cout << v_.id << "(" << v_.N_Vids.size() << ") -> ";
                if (v_.N_Vids.size() != 4 || v_.id == vid) {
                    // std::cout << "Non regular vertex" << v_.id << std::endl;
                    if (incluePaths) {
                        Path p;
                        p.vids = pv_.vids; p.b1 = pv_.b1; p.b2 = pv_.b2; p.rots = pv_.rots;
                        p.intersect_b1 = pv_.intersect_b1; p.intersect_b2 = pv_.intersect_b2; p.boundary_distance = pv_.boundary_distance;
                        paths.push_back(p);
                    } else {
                        int dist = pv_.b1 + pv_.b2;
                        bool same, diff, boundary;
                        // std::cout << "Getting terminal vertex type" << std::endl;
                        terminalVertexType(v.id, v_.id, same, diff, boundary);
                        // std::cout << "Got terminal vertex type" << std::endl;
                        if (same && dist <= same_singularity_min_dist && ++same_count <= 3) {
                            same_singularity_score = (((same_count-1)*same_singularity_score)+(1-(1/(double)dist)))/same_count;
                        }
                        if (diff && dist <= diff_singularity_max_dist && ++diff_count <= 3) {
                            diff_singularity_max_dist = dist + (dist/diff_count);
                            diff_singularity_score = (((diff_count-1)*diff_singularity_score)+(1/(double)dist))/diff_count;
                        }
                        if (pv_.boundary_distance > 0 && pv_.boundary_distance <= boundary_min_dist) {
                            boundary_min_dist = pv_.boundary_distance;
                            boundary_score = 1-(1/(double)pv_.boundary_distance);
                        }
                    }
                    continue;
                }
                // std::cout << pv_.vids.size() << " ";
                // if (++i == this->iters) break;
                // std::cout << "Getting neighbors of v_" << v_.id << "(" << v_.N_Vids.size() << ")" << std::endl;
                std::vector<size_t> neighbors_ = getVertexNeighbors(v_.id, pv_.vids.at(pv_.vids.size()-2));
                // std::cout << "Got neighbors of v_: " << neighbors_.size() << std::endl;
                // std::cout << neighbors_.size() << std::endl;
                path_v pv1, pv2, pv3;
                pv1.id = neighbors_.at(1); pv2.id = neighbors_.at(2); pv3.id = neighbors_.at(3);
                pv1.vids = pv_.vids; pv2.vids = pv_.vids; pv3.vids = pv_.vids;
                pv1.vids.push_back(pv1.id); pv2.vids.push_back(pv2.id); pv3.vids.push_back(pv3.id);
                // std::cout << "pv1 vids: " << pv1.vids.size() << " pv2 vids: " << pv2.vids.size() << " pv3 vids: " << pv3.vids.size() << " ";
                if (pv_.continuous) {
                    pv2.continuous = true;
                    pv1.b1 = pv_.b1; pv2.b1 = pv_.b1+1; pv3.b1 = pv_.b1;
                    pv1.b2 = 1; pv3.b2 = 1; pv1.rots = 1; pv3.rots = 3;
                    q.push(pv1); q.push(pv3); q.push(pv2);
                } else {
                    pv2.b2 += 1;
                    q.push(pv2);
                }
            }
            // std::cout << std::endl;
        }
        score = (0.25 * same_singularity_score) + (0.5 * diff_singularity_score) + (0.25 * boundary_score);
        return score;
    };

    // lambda function to set path singularity score
    auto getBestPaths = [&] (std::vector<Path>& paths) {
        std::priority_queue<Path, std::vector<Path>, p_comp> pq;
        for (auto& p: paths) {
            for (auto vid: p.vids) {
                p.vset.insert(vid);
                p.vset.insert(this->mesh->V.at(vid).N_Vids.begin(), this->mesh->V.at(vid).N_Vids.end());
                for (auto nvid: this->mesh->V.at(vid).N_Vids) p.vset.insert(this->mesh->V.at(nvid).N_Vids.begin(), this->mesh->V.at(nvid).N_Vids.end());
            }
            double path_qem_score_before = 0;
            double path_singularity_score_before = 0;
            int n1 = 0;
            int n2 = 0;
            for (auto vid: p.vset) {
                ++n1;
                // std::cout << "Getting vertex score for: " << vid << "(" << this->mesh->V.at(vid).N_Vids.size() << ")" << std::endl;
                auto v_score = this->mu->GetVertexScore(vid, true);
                // std::cout << "Got vertex score for: " << vid << "(" << this->mesh->V.at(vid).N_Vids.size() << ")" << std::endl;
                // this->mu->CalculateQEM(vid);
                path_qem_score_before = (((n1-1)*path_qem_score_before)+v_score->qem_score)/n1;
                if (isSingularity(vid)) {
                    ++n2;
                    double cur_sing_score = getPaths(vid, false, paths);
                    // auto v_score = this->mu->GetVertexScore(vid);
                    path_singularity_score_before = (((n2-1)*path_singularity_score_before)+cur_sing_score)/n2;
                }
            }
            double score_before = (path_qem_score_before + path_singularity_score_before)/2;
            int best_dir = -1;
            double best_score = 0;
            std::vector<size_t> dirs = getVertexNeighbors(p.vids.at(0), p.vids.at(1));
            for (int i = 1; i < dirs.size(); i++) {
                // std::cout << "Testing dir: " << dirs.at(i) << std::endl;
                size_t dir = dirs.at(i);
                this->caretaker.saveState(this->mesh);
                p.move_dir = dir;
                // std::cout << "Moving transient pair" << std::endl;
                // this->PrototypeSaveSeparatrices(std::vector<std::vector<size_t>>{p.vids}, "test3");
                this->PrototypeMoveAlongPath(p);
                // std::cout << "After moving transient pair" << std::endl;
                double path_qem_score_after = 0;
                double path_singularity_score_after = 0;
                n1 = 0; n2 = 0;
                for (auto vid: p.vset) {
                    ++n1;
                    // std::cout << "Getting vertex score for: " << vid << std::endl;
                    auto v_score = this->mu->GetVertexScore(vid, true);
                    // std::cout << "Got vertex score for: " << vid << std::endl;
                    // this->mu->CalculateQEM(vid);
                    path_qem_score_after = (((n1-1)*path_qem_score_after)+v_score->qem_score)/n1;
                    if (isSingularity(vid)) {
                        ++n2;
                        // std::cout << "Getting singularity score for: " << vid << "(" << this->mesh->V.at(vid).N_Vids.size() << ")" << std::endl;
                        
                        double cur_sing_score = getPaths(vid, false, paths);
                        // std::cout << "Got singularity score for: " << vid << "(" << this->mesh->V.at(vid).N_Vids.size() << ")" << std::endl;
                        // auto v_score = this->mu->GetVertexScore(vid);
                        path_singularity_score_after = (((n2-1)*path_singularity_score_after)+cur_sing_score)/n2;
                    }
                }
                double score_after = (path_qem_score_after + path_singularity_score_after)/2;
                double score = score_before - score_after;
                std::cout << "score: " << score << std::endl;
                if (score > best_score) {
                    best_score = score;
                    best_dir = dir;
                }
                // std::cout << "Before restoring state" << std::endl;
                this->caretaker.restoreState(this->mesh);
                // std::cout << "After restoring state" << std::endl;
            }
            if (best_dir >= 0 && best_score > 0) {
                p.move_dir = (size_t) best_dir;
                p.score = best_score;
                pq.push(p);
            }
        }
        return pq;
    };

    // getPaths(vid, true, paths);
    // std::priority_queue<Path, std::vector<Path>, p_comp> pq = getBestPaths(paths);
    // std::cout << pq.size() << std::endl;
    // std::cout << "move_dir: " << pq.top().move_dir << std::endl;
    // std::cout << "first connection: " << pq.top().vids.at(1) << std::endl;
    // std::vector<size_t> vids = {pq.top().move_dir};
    // for (auto vid: pq.top().vids) {
        // vids.push_back(vid);
    // }
    // for (auto vid: vids) std::cout << vid << " ";
    // std::cout << std::endl;
    // PrototypeSaveSeparatrices(std::vector<std::vector<size_t>>{vids}, "test");

    // separate path extraction and score calculation

    auto getSingularityScore = [&] (size_t vid) {
        struct v_trace {
            size_t curr;
            size_t prev;
            int dist = 0;
            bool continuous = false;
        };
        auto& v = this->mesh->V.at(vid);
        double score = 0;
        if ((v.N_Vids.size() != 3 && v.N_Vids.size() != 5) || v.isBoundary) return score;
        std::vector<size_t> neighbors = getVertexNeighbors(v.id, v.N_Vids.at(0));
        double boundary_score = 0;
        int boundary_min_dist = std::numeric_limits<int>::max();
        std::priority_queue<int> same_score_dist;
        std::priority_queue<int> diff_score_dist;
        std::vector<bool> visited(this->mesh->V.size(), false);
        
        for (auto nid: neighbors) {
            v_trace pv;
            pv.curr = nid;
            pv.prev = v.id;
            pv.dist = 1;
            pv.continuous = true;
            std::queue<v_trace> q;
            q.push(pv);
            while (!q.empty()) {
                auto pv_= q.front();
                q.pop();
                auto& v_ = this->mesh->V.at(pv_.curr);
                if (visited.at(v_.id)) continue;
                if (v_.N_Vids.size() != 4 || v_.isBoundary || v_.type == FEATURE) {
                    bool same = false; bool diff = false; bool boundary = false;
                    terminalVertexType(v.id, v_.id, same, diff, boundary);
                    if (same) {
                        // same_score_values.push(-exp(1-pv_.dist));
                        same_score_dist.push(pv_.dist);
                    }
                    if (diff) {
                        diff_score_dist.push(pv_.dist);
                    }
                    if (boundary) {
                        boundary_min_dist = pv_.dist;
                        boundary_score = exp(1-pv_.dist);
                    }
                    visited.at(v_.id) = true;
                    continue;
                }
                std::vector<size_t> neighbors_ = getVertexNeighbors(pv_.curr, pv_.prev);
                v_trace pv1, pv2, pv3;
                pv1.curr = neighbors_.at(1); pv2.curr = neighbors_.at(2); pv3.curr = neighbors_.at(3);
                pv1.prev = pv_.curr; pv2.prev = pv_.curr; pv3.prev = pv_.curr;
                pv1.dist = pv_.dist+1; pv2.dist = pv_.dist+1; pv3.dist = pv_.dist+1;
                if (pv_.continuous) {
                    pv2.continuous = true;
                    q.push(pv1); q.push(pv3); q.push(pv2);
                } else {
                    q.push(pv2);
                }
                visited.at(v_.id) = true;
            }
        }
        // int it = 0;
        // while (it++ < 3 && !score_values.empty()) {
        //     score += score_values.top();
        //     score_values.pop();
        // }
        score -= boundary_score;
        return score;
    };
    // parallelize paths extraction
    // parallelize path selection
    // add valence badness


    /*auto& v = mesh->V.at(vid);
    for (auto eid: v.N_Eids) {
        Path p;
        TraceAlongEdge(v, mesh->E.at(eid), p.intersect_b1, p.boundary_distance, p.vids);
        if (p.vids.empty()) continue;
        bool intersect_b1 = false;
        for (int i = 1; i < p.vids.size()-1; i++) {
            auto& sv = mesh->V.at(p.vids.at(i));
            intersect_b1 = (sv.isBoundary || sv.type == FEATURE);
            for (auto sveid: sv.N_Eids) {
                auto& e = mesh->E.at(sveid);
                if (Contains(e.Vids, p.vids.at(i-1)) || Contains(e.Vids, p.vids.at(i+1))) continue;
                Path sec_p;
                TraceAlongEdge(sv, e, sec_p.intersect_b2, sec_p.boundary_distance, sec_p.vids);
                if (p.boundary_distance > 0) sec_p.boundary_distance = p.boundary_distance;
                if (sec_p.vids.empty()) continue;
                if (!setPathParams(sec_p, i, sec_p.vids.size()-1, GetEdgeRots(GetEdgeId(p.vids.at(i), p.vids.at(i-1)), e.id), intersect_b1, false, p.vids)) continue;
                pq.push(sec_p);
                // paths.push_back(sec_p);
            }    
        }
        if (!setPathParams(p, p.vids.size()-1, 0, 0, p.intersect_b1, false, std::vector<size_t>{})) continue;
        pq.push(p);
        // paths.push_back(p);
    }
    for (auto fid: v.N_Fids) {
        Path p;
        TraceAlongDiagonal(v, mesh->F.at(fid), p.boundary_distance, p.vids);
        if (p.vids.empty()) continue;
        for (int i = 1; i < p.vids.size()-1; i++) {
            auto& sv = mesh->V.at(p.vids.at(i));
            for (auto svfid: sv.N_Fids) {
                auto& f = mesh->F.at(svfid);
                if (Contains(f.Vids, p.vids.at(i-1)) || Contains(f.Vids, p.vids.at(i+1))) continue;
                Path sec_p;
                TraceAlongDiagonal(sv, f, sec_p.boundary_distance, sec_p.vids);
                if (p.boundary_distance > 0) sec_p.boundary_distance = p.boundary_distance;
                if (!setPathParams(sec_p, i, sec_p.vids.size()-1, GetFaceRots(GetFaceId(p.vids.at(i), p.vids.at(i-1)), f.id, sv.id), false, true, p.vids)) continue;
                pq.push(sec_p);
                // paths.push_back(sec_p);
            }
        }
        if (!setPathParams(p, p.vids.size()-1, 0, 0, false, true, std::vector<size_t>{})) continue;
        pq.push(p);
        // paths.push_back(p);
    }*/
}

// PrototypeMoveAlongPath(Path& p) moves starting singularity in the dir and traces transient 3-5 pair along the path

void SemiGlobalSimplifier::PrototypeMoveAlongPath(Path& p) {
    // std::cout << "Inside PrototypeMoveAlongPath" << std::endl;
    // std::lock_guard<std::mutex> lock(mx);
    std::vector<size_t> vids = p.vids;
    auto& starting_singularity = mesh->V.at(vids.front());
    auto& terminal_singularity = mesh->V.at(vids.back());
    // std::cout << "start: " << starting_singularity.id << "(" << starting_singularity.N_Vids.size() << ") terminal: " << terminal_singularity.id << "(" << terminal_singularity.N_Vids.size() << ")" << std::endl;
    // std::cout << "vids: " << vids.size() << " ";
    // for (auto vid: vids) std::cout << vid << "(" << mesh->V.at(vid).N_Vids.size() << ")->";
    // std::cout << "move dir: " << p.move_dir << std::endl;
    // PrototypeSaveSeparatrices(std::vector<std::vector<size_t>>{vids}, "test");
    size_t vid = p.move_dir;

    std::shared_ptr<SingularityPair> sp;
    
    int secVid = -1;
    // if (vid == vids.at(2)) {
    //     for (auto fid: mesh->V.at(vids.at(0)).N_Fids) {
    //         auto& f = mesh->F.at(fid);
    //         if (Contains(f.Vids, vid)) {
    //             secVid = mu->GetDifference(f.Vids, std::vector<size_t>(vids.begin(), vids.begin()+3)).at(0);
    //             break;
    //         }
    //     }
    // }

    std::vector<size_t> tfIds;
    bool isDiagonal = PrototypeGetPairIds(vids.at(0), vid, tfIds, secVid);
    // std::cout << "isDiagonal: " << isDiagonal << std::endl;
    // std::cout << "tFIds: " << tfIds.size() << std::endl;
    if (tfIds.size() != 2) return;
    if (isDiagonal) {
        sp = std::make_shared<DiagonalThreeFivePair>(*mesh, *mu, *smoother, tfIds.at(0), tfIds.at(1));
        // std::cout << "DIAGONAL PAIR " << std::endl;
    } else {
        sp = std::make_shared<ThreeFivePair>(*mesh, *mu, *smoother, tfIds.at(0), tfIds.at(1));
        // std::cout << "DIRECT PAIR " << std::endl;
    }
    int startIdx = 0;
    std::vector<size_t> tfVids = mu->GetUnion(mesh->V.at(tfIds.at(0)).N_Vids, mesh->V.at(tfIds.at(1)).N_Vids);
    tfVids = mu->GetDifference(tfVids, std::vector<size_t>{vid});
    
    for (int i = 0; i < vids.size(); i++) {
        for (auto nvid: tfVids) {
            if (Contains(mesh->V.at(nvid).N_Vids, vids.at(i)) && i > startIdx) {
                vids.at(i-1) = nvid;
                startIdx = i-1;
            }
        }
    }
    RenderMesh();

    // PrototypeSaveSeparatrices(std::vector<std::vector<size_t>>{}, "test_begin");
    for (int i = startIdx; i < vids.size()-1; i++) {
        // std::cout << "Movement iteration: " << i << std::endl;
    // for (int i = startIdx; i < iters; i++) {
        // std::cout << "Moving to " << vids.at(i) << " Vids: " << mesh->V.at(vids.at(i)).N_Vids.size() << std::endl;
        // PrototypeSaveSeparatrices(std::vector<std::vector<size_t>>{}, "test_"+std::to_string(i)+"_before");
        sp->Move(vids.at(i), delta);
        RenderMesh();
        // PrototypeSaveSeparatrices(std::vector<std::vector<size_t>>{}, "test_"+std::to_string(i)+"_after");
    }
    sp->Move(vids.back(), delta, false);
    RenderMesh();
    // PrototypeSaveSeparatrices(std::vector<std::vector<size_t>>{}, "test_end");
}

bool SemiGlobalSimplifier::TestFlips() {
    // auto testFlipQuads = [&]() {
    //     auto qV_arr = [&] (Face& f) {
    //         double coords[4][3];
    //         for (int i = 0; i < 4; i++) {
    //             auto& v = mesh->V.at(f.Vids.at(i));
    //             coords[i][0] = v.x; coords[i][1] = v.y; coords[i][2] = v.z;
    //         }
    //         return coords;
    //     };
    //     for (auto& f: mesh->F) {
    //         std::cout << "aspect_ratio_score for: " << f.id << " " << v_quad_aspect_ratio(4, qV_arr(f)) << std::endl; 
    //         std::cout << "shape_score for: " << f.id << " " << v_quad_shape(4, qV_arr(f)) << std::endl; 
    //         std::cout << "---------------------" << std::endl;
    //     }
    //     return 0;
    // }();
    // return;
    // for (int mainIter = 0; mainIter < iters; mainIter++) {
    // for (int mainIter = 0; mainIter < 1; mainIter++) {
        // std::cout << "MAIN ITER: " << mainIter << std::endl;
    // std::vector<std::vector<size_t>> paths;
    // for (auto& v: mesh->V) {
    //     vMesh m(mesh);
    //     PerformOperation(Operation("Rotate", v.id, std::vector<size_t>{}, false), &m);
    //     m.Update();
    // }
    // return 0;
    int mainIter = 0;
    int singularity_idx = 175;
    int path_idx = 991;
    int iteration_idx = 2;
    int mainIter_idx = 0;
    int singularity_it = 0;
    int path_it = 0;
    int iteration_it = 0;
    auto isSingularity = [&](Vertex& v) {
        if (v.isBoundary || v.type == FEATURE) return false;
        // if (v.N_Fids.size() == 3 || v.N_Fids.size() == 5) return true;
        // return false;
        return (mesh->virtualValence(v) == 3 || mesh->virtualValence(v) == 5);
    };
    int it = 0;
    const int MAIN = 1, MAIN_BRANCH = 2, BRANCH = 3, DEFLECT = 4, DEFLECT_SKIP = 5, SKIP = 6;
    struct vPath {
        int id;
        int flag = 1;
        vPath(int id_, int flag_) : id(id_), flag(flag_) {}
    };
    struct executablePath {
        vMesh m;
        std::vector<size_t> path;
        double score = -1.0;
        double valence_score = 0.0;
        double singularity_score = 0.0;
        double element_score = 0.0;
        bool operator<(const executablePath& p) const {return score < p.score;}
        int singularity_it = 0; int iteration_it = 0; int path_it = 0; int mainIter = 0;
        executablePath() {}
        executablePath& operator=(const executablePath& other) {
            if (this == &other) return *this;
            m = other.m;
            path = other.path;
            score = other.score;
            valence_score = other.valence_score;
            singularity_score = other.singularity_score;
            element_score = other.element_score;
            singularity_it = other.singularity_it;
            iteration_it = other.iteration_it;
            path_it = other.path_it;
            mainIter = other.mainIter;
            return *this;
        }
    };
    struct pathQueueItem {
        pathQueueItem() {}
        pathQueueItem(std::vector<size_t> path_, int pth_id_) {
            path = path_;
            pth_id = pth_id_;
        }
        std::vector<size_t> path;
        int pth_id;
    };
    std::priority_queue<executablePath> paths_;
        
    auto corner = [&] (size_t vid, vMesh* m = nullptr, bool useVM = true) {
        if (m == nullptr) m = new vMesh(mesh);
        auto& v = m->getVertex(vid, useVM);
        if (!v.isBoundary && v.type != FEATURE) return false;
        int ideal_valence = m->getIdealValence(v.id, useVM);
        if (v.isBoundary && ideal_valence != 2) return true;
        if (v.type == FEATURE && ideal_valence != 4) return true;
        vInfo info_v(mesh, v.id, m, useVM);
        auto nvids = info_v.vids();
        std::vector<size_t> features;
        for (auto nvid: nvids) {
            auto& nv = m->getVertex(nvid, useVM);
            if (nv.isBoundary || nv.type == FEATURE) features.push_back(nvid);
        }
        if (features.size() == 2) {
            double angle = m->getAngle(v.id, features.at(0), features.at(1), useVM);
            if (angle < 165 || angle > 195) return true;
        } else if (features.size() > 2) {
            features.push_back(features.at(0));
            std::vector<glm::dvec3> normals;
            for (int i = 0; i < features.size()-1; i++) {
                auto& v1 = m->getVertex(features.at(i), useVM);
                auto& v2 = m->getVertex(features.at(i+1), useVM);
                normals.push_back(glm::normalize(glm::cross(v1.xyz()-v.xyz(), v2.xyz()-v.xyz())));
            }
            for (int i = 0; i < normals.size()-1; i++) {
                auto A = normals.at(i);
                auto B = normals.at(i+1);
                double angle = atan2(glm::length(glm::cross(A, B)), glm::dot(A, B));
                if (angle < 0) angle += 2*PI;
                angle = angle * 180.0 / PI;
                if (angle < 165 || angle > 195) return true; 
            }
        }
        return false;
    };
    auto renewMesh = [&] () {
        std::vector<Vertex> V;
        std::vector<Edge> E;
        std::vector<Face> F;
        std::unordered_map<size_t, size_t> vmap;
        size_t vidx = 0;
        for (auto& v: mesh->V) {
            if (v.N_Fids.empty()) continue;
            Vertex v_(v);
            v_.id = V.size();
            v_.N_Fids.clear();
            V.push_back(v_);
            vmap[v.id] = v_.id;
        }
        for (auto& f: mesh->F) {
            if (f.Vids.empty()) continue;
            for (int i = 0; i < f.Vids.size(); i++) {
                f.Vids.at(i) = vmap[f.Vids.at(i)];
            }
            Face f_(f);
            f_.id = F.size();
            F.push_back(f_);
        }
        for (auto& f: F) {
            for (int i = 0; i < f.Vids.size(); i++) {
                auto& v = V.at(f.Vids.at(i));
                mu->AddContents(v.N_Fids, std::vector<size_t>{f.id});
            }
        }
        mesh->V.swap(V);
        mesh->E.swap(E);
        mesh->F.swap(F);
    };
    auto smooth = [&] (vMesh* m = nullptr, int smoothIters = 5) {
        int it = 0;
        bool useVM = m != nullptr;
        std::vector<size_t> V;
        // std::unordered_map<int, glm::dvec3> new_xyz;
        std::vector<glm::dvec3> new_xyz;
        if (useVM) {
            for (auto it_map = m->vmap.begin(); it_map != m->vmap.end(); it_map++) {
                V.push_back(it_map->second.id);
                new_xyz.push_back(it_map->second.xyz());
            }
        } else {
            V.resize(mesh->V.size());
            new_xyz.resize(mesh->V.size());
            PARALLEL_FOR_BEGIN(0, mesh->V.size()) {
                V.at(i) = i;
                new_xyz.at(i) = mesh->V.at(i).xyz();
            } PARALLEL_FOR_END();
        }
        if (m == nullptr) m = new vMesh(mesh);
        bool log = mainIter == mainIter_idx && path_it == path_idx && singularity_it == singularity_idx && iteration_it == iteration_idx;
        log = false;
        auto getSmoothPoint = [&] (Vertex& v, Vertex& nv, Vertex& prev, Vertex& next) {
            auto A = v.xyz() - nv.xyz();
            auto B = prev.xyz() - nv.xyz();
            auto C = next.xyz() - nv.xyz();

            double theta_1 = atan2(glm::length(glm::cross(A, B)), glm::dot(A, B));
            double theta_2 = atan2(glm::length(glm::cross(A, C)), glm::dot(A, C));
            double theta = (theta_2 - theta_1) / 2.0;

            glm::dvec3 r(v.xyz());
            double l = 0.0;
            if (theta > 0) {
                r = next.xyz() - v.xyz();
                l = fabs(theta/theta_2) * glm::length(r);
            } else if (theta < 0) {
                r = prev.xyz() - v.xyz();
                l = fabs(theta/theta_1) * glm::length(r);
            }
            return v.xyz()+(l*glm::normalize(r));
        };
        auto normalize = [&] (glm::dvec3 v) {
            if (glm::length(v) == 0.0) return glm::dvec3(0.0, 0.0, 0.0);
            return glm::normalize(v);
        };
        auto vertexNormal = [&] (Vertex& v) {
            glm::dvec3 normal(0.0, 0.0, 0.0);
            auto n = [&] (size_t vid_1, size_t vid_2, size_t vid_3) {
                glm::dvec3 a = m->getVertex(vid_2, useVM).xyz() - m->getVertex(vid_1, useVM).xyz();
                glm::dvec3 b = m->getVertex(vid_3, useVM).xyz() - m->getVertex(vid_1, useVM).xyz();
                return normalize(glm::cross(a, b));
            };
            for (auto fid: v.N_Fids) {
                auto& f = m->getFace(fid, useVM);
                normal += n(f.Vids.at(0), f.Vids.at(1), f.Vids.at(2));
                normal += n(f.Vids.at(0), f.Vids.at(2), f.Vids.at(3));
                normal += n(f.Vids.at(1), f.Vids.at(2), f.Vids.at(3));
                normal += n(f.Vids.at(1), f.Vids.at(3), f.Vids.at(0));
            }
            return normalize(normal);
        };
        auto projectToNormal = [&] (glm::dvec3 x, glm::dvec3 p, glm::dvec3 n) {
            // std::cout << "dot p,n: " << glm::dot(p-x, n) << std::endl;
            // auto toSub = glm::dot(p-x, n) * n;
            // std::cout << "to sub from p: " << toSub.x << " " << toSub.y << " " << toSub.z << std::endl;
            return p - (glm::dot(p-x, n) * n);
        };
        auto smoothInterior = [&] (Vertex& v, int id) {
            vInfo info_v(mesh, v.id, m, useVM);
            auto nvids = info_v.vids();
            double polyArea = 0.0;
            glm::dvec3 centroid(0.0, 0.0, 0.0);
            for (int nvid = 0; nvid < nvids.size(); nvid++) {
                auto& v2 = m->getVertex(nvids.at(nvid), useVM); // A
                auto& v3 = m->getVertex(nvids.at((nvid+nvids.size()-1)%nvids.size()), useVM); // B
                auto& v4 = m->getVertex(nvids.at((nvid+1)%nvids.size()), useVM); // C
                
                auto A = v2.xyz();
                auto B = v3.xyz();
                auto C = v4.xyz();

                glm::dvec3 AB = B - A;
                glm::dvec3 BC = C - B;
                glm::dvec3 CA = A - C;
                glm::dvec3 AC = C - A;

                double a = glm::length(BC);
                double b = glm::length(CA);
                double c = glm::length(AB);

                auto diag_1 = v.xyz() - A;
                auto diag_2 = B - C;
                
                // glm::dvec3 incenter = ((a * v2.xyz()) + (b * v3.xyz()) + (c * v4.xyz())) / (a + b + c);
                glm::dvec3 incenter = glm::dvec3((a*A.x)+(b*B.x)+(c*C.x), (a*A.y)+(b*B.y)+(c*C.y), (a*A.z)+(b*B.z)+(c*C.z));
                incenter /= (a+b+c);
                // double area = 0.5 * glm::length(glm::cross(AB, AC));
                double area = 0.5 * glm::length(glm::cross(diag_1, diag_2));
                // centroid += (area * incenter); 
                centroid += (A + normalize(incenter - A) * glm::length(v.xyz() - A)) * area;
                polyArea += area;
            }
            if (polyArea == 0.0) return;
            // std::cout << "v " << v.id << " pos calculation" << std::endl;
            // std::cout << "current pos: " << v.x << " " << v.y << " " << v.z << std::endl;
            // auto new_pos = (centroid/polyArea);
            // std::cout << "new pos: " << new_pos.x << " " << new_pos.y << " " << new_pos.z << std::endl;
            // auto normal = vertexNormal(v);
            // std::cout << "normal: " << normal.x << " " << normal.y << " " << normal.z << std::endl;
            // auto projected_point = projectToNormal(v.xyz(), new_pos, normal);
            // std::cout << "projected point: " << projected_point.x << " " << projected_point.y << " " << projected_point.z << std::endl;
            // centroid /= polyArea;
            v.xyz(v.xyz() + (1.0 * ((centroid/polyArea) - v.xyz())));
            // auto normal = vertexNormal(v);
            // v.xyz(projectToNormal(v.xyz(), centroid/polyArea, vertexNormal(v)));
            // new_xyz[id] = centroid/polyArea;
        };
        auto smoothBoundary = [&] (Vertex& v, int id) {
            if (corner(v.id, m, useVM)) return;
            vInfo info_v(mesh, v.id, m, useVM);
            auto nvids = info_v.vids();
            std::vector<size_t> boundaryVertices;
            for (auto nvid: nvids) {
                auto& nv = m->getVertex(nvid, useVM);
                if (nv.isBoundary || nv.type == FEATURE) boundaryVertices.push_back(nvid);
            }
            // if (v.N_Fids.size() == 6) std::cout << "boundary vertices: " << boundaryVertices.size() << std::endl;
            if (boundaryVertices.size() > 2) return;
            if (boundaryVertices.size() == 1) {
                // std::cout << "boundary vertex has one boundary neighbor: " << v.id << std::endl;
                smoothInterior(v, id);
                return;
            }
            auto& b1 = m->getVertex(boundaryVertices.at(0), useVM);
            auto& b2 = m->getVertex(boundaryVertices.at(1), useVM);
            glm::dvec3 centroid(0.0, 0.0, 0.0);
            int k = 0;
            auto r1 = b1.xyz() - v.xyz();
            auto r2 = b2.xyz() - v.xyz();
            double l1 = 0.0;
            double l2 = 0.0;
            for (auto nvid: nvids) {
                if (nvid == b1.id || nvid == b2.id) continue;
                auto& nv = m->getVertex(nvid, useVM);
                glm::dvec3 A(v.xyz() - nv.xyz());
                glm::dvec3 B(b1.xyz() - nv.xyz());
                glm::dvec3 C(b2.xyz() - nv.xyz());
                double theta_1 = atan2(glm::length(glm::cross(A, B)), glm::dot(A, B));
                double theta_2 = atan2(glm::length(glm::cross(A, C)), glm::dot(A, C));
                double theta = (theta_2 - theta_1) / 2.0;
                glm::dvec3 r(v.xyz());
                double l = 0.0;
                if (theta > 0) {
                    r = b2.xyz() - v.xyz();
                    l += fabs(theta/theta_2) * glm::length(r);
                } else if (theta < 0) {
                    r = b1.xyz() - v.xyz();
                    l += fabs(theta/theta_1) * glm::length(r);
                }
                v.xyz(projectToNormal(v.xyz(), v.xyz()+l*normalize(r), vertexNormal(v)));
                // v.xyz(v.xyz()+l*glm::normalize(r));
                // centroid += (v.xyz()+l*glm::normalize(r));
                // k++;
            }
            // if (l1 > l2) {
            //     v.xyz(v.xyz()+((l1-l2)*glm::normalize(r1)));
            // } else {
            //     v.xyz(v.xyz()+((l2-l1)*glm::normalize(r2)));
            // }
            // if (k > 0) v.xyz(centroid/(double) k);
        };
        
        auto calculatePos = [&] (int id, Vertex& v) {
            // if (useVM) std::cout << "calculatePos point " << id << std::endl;
            if (v.N_Fids.empty()) return;
            // vInfo info_v(mesh, v.id, m, useVM);
            // auto nvids = info_v.vids();
            if (v.isBoundary || v.type == FEATURE) {
                // smoothBoundary(v, id);
                // if (corner(v.id, m, useVM)) return;
                // if (useVM) std::cout << "Boundary Vertex" << std::endl;
                // std::vector<size_t> boundaryVertices;
                // for (auto nvid: nvids) {
                //     auto& nv = m->getVertex(nvid, useVM);
                //     if (nv.isBoundary || nv.type == FEATURE) boundaryVertices.push_back(nvid);
                // }
                // if (boundaryVertices.size() != 2) return;
                // auto& b1 = m->getVertex(boundaryVertices.at(0), useVM);
                // auto& b2 = m->getVertex(boundaryVertices.at(1), useVM);
                // std::vector<glm::dvec3> points;
                // for (auto nvid: nvids) {
                //     if (nvid == boundaryVertices.at(0) || nvid == boundaryVertices.at(1)) continue;
                //     points.push_back(getSmoothPoint(v, m->getVertex(nvid, useVM), b1, b2));
                // }
                // glm::dvec3 p = points.at(0);
                // for (int pid = 1; pid < points.size(); pid++) {
                //     if (glm::length(p-v.xyz()) > glm::length(points.at(pid)-v.xyz())) p = points.at(pid);
                // }
                // kd->SearchAndProject(p);
                // new_xyz.at(id) = p;
            } else {
                smoothInterior(v, id);
                // if (useVM) std::cout << "Normal Vertex" << std::endl;
                // glm::dvec3 p(0.0);
                // int np = 0;
                // for (int nvid = 0; nvid < nvids.size(); nvid++) {
                //     p += getSmoothPoint(v, m->getVertex(nvids.at(nvid), useVM), 
                //         m->getVertex(nvids.at((nvid+nvids.size()-1)%nvids.size()), useVM), m->getVertex(nvids.at((nvid+1)%nvids.size()), useVM));
                //     np++;
                // }
                // if (np > 0) p /= np;
                // if (useVM) std::cout << "calculating projection" << std::endl;
                // kd->SearchAndProject(p);
                // if (useVM) std::cout << "calculated projection" << std::endl;
                // new_xyz.at(id) = p;
            }
        };
        auto setPos = [&] (int id, Vertex& v) {
            if (v.type == FEATURE || v.isBoundary) return;
            auto p = v.xyz();
            glm::dvec3 n(0.0, 0.0, 0.0);
            for (auto fid: v.N_Fids) {
                auto& f = m->getFace(fid, useVM);
                int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
                const auto& AB = m->getVertex(f.Vids.at((idx+1)%f.Vids.size()), useVM).xyz() - v.xyz();
                const auto& AC = m->getVertex(f.Vids.at((idx+2)%f.Vids.size()), useVM).xyz() - v.xyz();
                const auto& AD = m->getVertex(f.Vids.at((idx+3)%f.Vids.size()), useVM).xyz() - v.xyz();
                auto T_area_a = glm::cross(AB, AC);
                auto T_area_b = glm::cross(AC, AD);
                n += (0.5 * glm::length(T_area_a) * T_area_a);
                n += (0.5 * glm::length(T_area_b) * T_area_b);
            }
            n = glm::normalize(n);
            
            kd->SearchAndProject(p, n);
            // v.xyz(new_xyz.at(id));
            v.xyz(p);
        };
        while (it < smoothIters) {
            // std::cout << "iter: " << it << std::endl;
            // std::cout << "V: " << V.size() << std::endl;
            // if (useVM) {
            //     for (int i = 0; i < V.size(); i++) {
            //         setPos(i, m->getVertex(V.at(i), useVM));
            //     }
            // } else {
            //     PARALLEL_FOR_BEGIN(0, V.size()) {
            //         setPos(i, m->getVertex(V.at(i), useVM));
            //     } PARALLEL_FOR_END();
            // }

            if (useVM) {
                for (int i = 0; i < V.size(); i++) {
                    // std::cout << "Calculating point " << i << std::endl;
                    calculatePos(i, m->getVertex(V.at(i), useVM));
                    // std::cout << "Calculated point " << i << std::endl;
                }
            } else {
                // PARALLEL_FOR_BEGIN(0, V.size()) {
                //     calculatePos(i, m->getVertex(V.at(i), useVM));
                // } PARALLEL_FOR_END();
                int it_ = 0;
                for (int i = 0; i < V.size(); i++) {
                    // std::cout << "Calculating point " << i << std::endl;
                    // auto& v = m->getVertex(V.at(i), useVM);
                    // if (v.isBoundary || v.type == FEATURE) continue;
                    calculatePos(i, m->getVertex(V.at(i), useVM));
                    // if (++it_ > 9) break;
                    // std::cout << "Calculated point " << i << std::endl;
                }
                // PARALLEL_FOR_BEGIN(0, V.size()) {
                //     m->getVertex(V.at(i), useVM).xyz(new_xyz.at(i));
                // } PARALLEL_FOR_END();
            }
            if (useVM) {
                // for (int i = 0; i < V.size(); i++) {
                //     setPos(i, m->getVertex(V.at(i), useVM));
                // }
            } else {
                PARALLEL_FOR_BEGIN(0, V.size()) {
                    setPos(i, m->getVertex(V.at(i), useVM));
                } PARALLEL_FOR_END();
            }
            // for (auto vid: V) {
            /*for (int i = 0; i < V.size(); i++) {
                auto vid = V.at(i);
                // auto& v = m->getVertex(vid);
                // auto& v = it_map->second;
                auto& v = m->getVertex(vid, useVM);
                if (v.N_Fids.empty()) continue;
                // if (path_it == path_idx && iteration_it == iteration_idx) std::cout  << "smoothing " << vid << std::endl;
                if (v.isBoundary || v.type == FEATURE) {
                    continue;
                    // if (path_it == path_idx && iteration_it == iteration_idx) continue;
                    // if (iteration_it != iteration_idx || path_it != path_idx) continue;
                    // if (path_it == path_idx && iteration_it == iteration_idx) std::cout << "smoothing boundary" << std::endl;
                    // bool log = path_it == path_idx && iteration_it == iteration_idx;
                    // int ideal_valence = m->getIdealValence(v.id, useVM, log);1
                    // int ideal_valence = m->getIdealValence(v.id, useVM);
                    // if (path_it == path_idx && iteration_it == iteration_idx) std::cout << "ideal valence: " << ideal_valence << std::endl;
                    // if ((v.isBoundary && ideal_valence != 2) || (v.type == FEATURE && ideal_valence != 4)) continue;
                    if (corner(v.id, m, useVM)) continue;                
                    vInfo info_v(mesh, v.id, m, useVM);
                    // std::cout << "boundary: " << v.id << (v.isBoundary ? " yes" : " no") << std::endl;
                    auto nvids = info_v.vids();
                    std::vector<size_t> boundaryVertices;
                    for (auto nvid: nvids) {
                        auto& nv = m->getVertex(nvid, useVM);
                        if (nv.isBoundary || nv.type == FEATURE) boundaryVertices.push_back(nvid);
                    }
                    if (boundaryVertices.size() != 2) continue;
                    // std::cout << "boundaryVertices: " << boundaryVertices[0] << " " << boundaryVertices[1] << std::endl;
                    // std::cout << "nvids: "; for (auto nvid: nvids) std::cout << nvid << " "; std::cout << std::endl;
                    auto& b1 = m->getVertex(boundaryVertices.at(0), useVM);
                    auto& b2 = m->getVertex(boundaryVertices.at(1), useVM);
                    // if ([&] () {
                    //     double angle = m->getAngle(v.id, b1.id, b2.id, useVM);
                    //     if (angle < 165 || angle > 195) return true;
                    //     return false;
                    // }()) continue;
                    for (auto nvid: nvids) {
                        // std::cout << "nvid: " << nvid << std::endl;
                        // std::cout << (nvid == boundaryVertices.at(0)) << " " << (nvid == boundaryVertices.at(1)) << std::endl;
                        if (nvid == boundaryVertices.at(0) || nvid == boundaryVertices.at(1)) continue;
                        auto& nv = m->getVertex(nvid, useVM);
                        // std::cout << "nv: " << nv.id << std::endl;
                        auto A = v.xyz() - nv.xyz();
                        auto B = b1.xyz() - nv.xyz();
                        auto C = b2.xyz() - nv.xyz();

                        double theta_1 = atan2(glm::length(glm::cross(A, B)), glm::dot(A, B));
                        double theta_2 = atan2(glm::length(glm::cross(A, C)), glm::dot(A, C));

                        double theta = (theta_2 - theta_1) / 2.0;

                        glm::dvec3 r(v.xyz());
                        double l = 0.0;
                        if (theta > 0) {
                            r = b2.xyz() - v.xyz();
                            l = fabs(theta/theta_2) * glm::length(r);
                        } else if (theta < 0) {
                            r = b1.xyz() - v.xyz();
                            l = fabs(theta/theta_1) * glm::length(r);
                        }
                        // std::cout << "l: " << l << std::endl;
                        // v.xyz(v.xyz()+(l*glm::normalize(r)));
                        new_xyz.at(i) = v.xyz()+(l*glm::normalize(r));
                        break;
                    }
                     
                } else {
                    // continue;
                    vInfo info_v(mesh, v.id, m, useVM);
                    auto nvids = info_v.vids();
                    double polyArea = 0.0;
                    glm::dvec3 centroid(0.0, 0.0, 0.0);
                    if (log) {
                        std::cout << "v: " << v.x << " " << v.y << " " << v.z << std::endl;
                        std::cout << "nfids: " << v.N_Fids.size() << " " << info_v.fids().size() << std::endl;
                        std::cout << "nvids: " << nvids.size() << std::endl;
                    }
                    for (int i = 0; i < nvids.size(); i++) {
                        auto v2 = m->getVertex(nvids.at(i), useVM).xyz();
                        auto v3 = m->getVertex(nvids.at((i+nvids.size()-1)%nvids.size()), useVM).xyz();
                        auto v4 = m->getVertex(nvids.at((i+1)%nvids.size()), useVM).xyz();
                        
                        glm::dvec3 AB = v3-v2; glm::dvec3 BC = v4-v3; glm::dvec3 CA = v2-v4; glm::dvec3 AC = v4-v2;

                        double a = glm::length(BC); double b = glm::length(CA); double c = glm::length(AB);
                        glm::dvec3 incenter = ((a * v2) + (b * v3) + (c * v4)) / (a + b + c);
                        
                        double area = 0.5 * glm::length(glm::cross(AB, AC));
                        if (log) {
                            std::cout << "v2: " << v2.x << " " << v2.y << " " << v2.z << std::endl;
                            std::cout << "v3: " << v3.x << " " << v3.y << " " << v3.z << std::endl;
                            std::cout << "v4: " << v4.x << " " << v4.y << " " << v4.z << std::endl;
                            std::cout << "AB: " << AB.x << " " << AB.y << " " << AB.z << std::endl;
                            std::cout << "BC: " << BC.x << " " << BC.y << " " << BC.z << std::endl;
                            std::cout << "CA: " << CA.x << " " << CA.y << " " << CA.z << std::endl;
                            std::cout << "AC: " << AC.x << " " << AC.y << " " << AC.z << std::endl;
                            std::cout << "a: " << a << " b: " << b << " c: " << c << std::endl;
                            std::cout << "incenter: " << incenter.x << " " << incenter.y << " " << incenter.z << std::endl;
                            std::cout << "area: " << area << std::endl;
                        }
                        centroid += (area * incenter); 
                        polyArea += area;
                    }
                    if (polyArea == 0.0) continue;
                    if (log) {
                        std::cout << "centroid: " << centroid.x << " " << centroid.y << " " << centroid.z << std::endl;
                        std::cout << "polyArea: " << polyArea << std::endl;
                        std::cout << "***************************" << std::endl;
                    }
                    centroid /= polyArea;
                    double min_distance = std::numeric_limits<double>::max();
                    glm::dvec3 intersection(0.0, 0.0, 0.0);
                    for (auto fid: info_v.fids()) {
                        auto& f = m->getFace(fid);
                        int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
                        Plane p(v, m->getVertex(f.Vids.at((idx+1)%f.Vids.size()), useVM), m->getVertex(f.Vids.at((idx+3)%f.Vids.size()), useVM));
                        glm::dvec3 intersect(0.0, 0.0, 0.0);
                        double distance = p.DistanseFromPoint(centroid, intersect);
                        if (distance < min_distance) {
                            min_distance = distance;
                            intersection = intersect;
                        }
                    }
                    // v.xyz(intersection);
                    new_xyz.at(i) = intersection;
                }
            }
            for (int i = 0; i < V.size(); i++) {
                auto& v = m->getVertex(V.at(i), useVM);
                v.xyz(new_xyz.at(i));
            }
            for (int i = 0; i < V.size(); i++) {
                auto vid = V.at(i);
                auto& v = m->getVertex(vid, useVM);
                if (v.N_Fids.empty()) continue;
                if (v.isBoundary || v.type == FEATURE) {
                    // if (path_it == path_idx && iteration_it == iteration_idx) continue;
                    // if (iteration_it != iteration_idx || path_it != path_idx) continue;
                    // if (path_it == path_idx && iteration_it == iteration_idx) std::cout << "smoothing boundary" << std::endl;
                    // bool log = path_it == path_idx && iteration_it == iteration_idx;
                    // int ideal_valence = m->getIdealValence(v.id, useVM, log);1
                    // int ideal_valence = m->getIdealValence(v.id, useVM);
                    // if (path_it == path_idx && iteration_it == iteration_idx) std::cout << "ideal valence: " << ideal_valence << std::endl;
                    // if ((v.isBoundary && ideal_valence != 2) || (v.type == FEATURE && ideal_valence != 4)) continue;
                    if (corner(v.id, m, useVM)) continue;                
                    vInfo info_v(mesh, v.id, m, useVM);
                    if (log) std::cout << "boundary: " << v.id << (v.isBoundary ? " yes" : " no") << std::endl;
                    auto nvids = info_v.vids();
                    std::vector<size_t> boundaryVertices;
                    for (auto nvid: nvids) {
                        auto& nv = m->getVertex(nvid, useVM);
                        if (nv.isBoundary || nv.type == FEATURE) boundaryVertices.push_back(nvid);
                    }
                    if (boundaryVertices.size() != 2) continue;
                    // std::cout << "boundaryVertices: " << boundaryVertices[0] << " " << boundaryVertices[1] << std::endl;
                    // std::cout << "nvids: "; for (auto nvid: nvids) std::cout << nvid << " "; std::cout << std::endl;
                    auto& b1 = m->getVertex(boundaryVertices.at(0), useVM);
                    auto& b2 = m->getVertex(boundaryVertices.at(1), useVM);
                    // if ([&] () {
                    //     double angle = m->getAngle(v.id, b1.id, b2.id, useVM);
                    //     if (angle < 165 || angle > 195) return true;
                    //     return false;
                    // }()) continue;
                    for (auto nvid: nvids) {
                        // std::cout << "nvid: " << nvid << std::endl;
                        // std::cout << (nvid == boundaryVertices.at(0)) << " " << (nvid == boundaryVertices.at(1)) << std::endl;
                        if (nvid == boundaryVertices.at(0) || nvid == boundaryVertices.at(1)) continue;
                        auto& nv = m->getVertex(nvid, useVM);
                        // std::cout << "nv: " << nv.id << std::endl;
                        auto A = v.xyz() - nv.xyz();
                        auto B = b1.xyz() - nv.xyz();
                        auto C = b2.xyz() - nv.xyz();

                        double theta_1 = atan2(glm::length(glm::cross(A, B)), glm::dot(A, B));
                        double theta_2 = atan2(glm::length(glm::cross(A, C)), glm::dot(A, C));

                        double theta = (theta_2 - theta_1) / 2.0;

                        glm::dvec3 r(v.xyz());
                        double l = 0.0;
                        if (theta > 0) {
                            r = b2.xyz() - v.xyz();
                            l = fabs(theta/theta_2) * glm::length(r);
                        } else if (theta < 0) {
                            r = b1.xyz() - v.xyz();
                            l = fabs(theta/theta_1) * glm::length(r);
                        }
                        // std::cout << "l: " << l << std::endl;
                        // v.xyz(v.xyz()+(l*glm::normalize(r)));
                        new_xyz.at(i) = v.xyz()+(l*glm::normalize(r));
                        break;
                    }
                     
                }
            }
            for (int i = 0; i < V.size(); i++) {
                auto& v = m->getVertex(V.at(i), useVM);
                v.xyz(new_xyz.at(i));
            }*/
            it++;
        }
        if (!useVM) {
            // PARALLEL_FOR_BEGIN(0, V.size()) {
            //     setPos(i, m->getVertex(V.at(i), useVM));
            // } PARALLEL_FOR_END();
        }
    };
    auto FixDoublets = [&] () {
        // std::cout << "Fixing doublets" << std::endl;
        vMesh m(mesh);
        bool res = false;
        std::vector<size_t> b_irr;
        for (auto& v: mesh->V) {
            // std::cout << "checking vertex " << v.id << "(" << v.N_Fids.size() << ")" << std::endl;
            // if (v.N_Fids.size() == 6) {
            //     vInfo info_v(mesh, v.id, &m);
            //     PrototypeSaveSeparatrices(std::vector<std::vector<size_t>> {std::vector<size_t>{v.id, info_v.vids().at(0)}}, "test");
            // }
            if (v.N_Fids.empty() || v.isBoundary || v.type == FEATURE) continue;
            if (m.valence(v) == 2) {
                b_irr.push_back(v.id);
            }
            // std::cout << "after valence calculation" << std::endl;
            // if (v.N_Fids.size() == 2) b_irr.push_back(v.id);
        }
        // vMesh m(mesh);
        // std::cout << "# doublets " << b_irr.size() << std::endl;
        for (auto id: b_irr) {
            m.Clear();
            auto& v = m.getVertex(id);
            if (v.N_Fids.empty()) continue;
            auto& f = m.getFace(v.N_Fids.at(0));
            int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
            // std::cout << "Before doublet collapse" << std::endl;
            if (!PerformOperation(Operation("Collapse", v.id, std::vector<size_t>{f.Vids.at((idx+2)%f.Vids.size()), v.id}, false), &m)) continue;
            // std::cout << "After doublet collapse" << std::endl;
            m.Update();
            res = true;
        }
        return res;
    };
    auto fixFeatures = [&] () {
        std::cout << "Fixing features" << std::endl;
        while (FixDoublets());
        bool res = false;
        std::unordered_set<size_t> b_irr;
        for (auto& v: mesh->V) {
            if (v.N_Fids.empty()) continue;
            if (v.isBoundary || v.type == FEATURE) {
                if (mesh->valence(v) == 2) {
                    b_irr.insert(v.id);
                    continue;
                }
                auto plane = mesh->planes(v);
                for (auto p: plane) {
                    if (mesh->virtualValence(v, p) != 4) {
                        b_irr.insert(v.id);
                    }
                }
            }
        }
        vMesh m(mesh);
        int it = 0;
        for (auto id: b_irr) {
            // std::cout << "id: " << id << std::endl;
            auto& v = m.getVertex(id);
            if (v.N_Fids.empty()) continue;
            if (m.valence(v) == 2) {
                // std::cout << "Fixing Boundary Doublet" << std::endl;
                // m.Clear();            
                // vInfo info_v(mesh, v.id, &m);
                // PerformOperation(Operation("Rotate", v.id, std::vector<size_t>{}, false), &m);
                // info_v = vInfo(mesh, v.id, &m);
                // for (auto fid: info_v.fids()) {
                //     auto& f = m.getFace(fid);
                //     int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
                //     PerformOperation(Operation("Collapse", v.id, std::vector<size_t>{v.id, f.Vids.at((idx+2)%f.Vids.size())}, false), &m);
                // }
                // m.Update();
                // res = true;
                // continue;
            }
            auto plane = m.planes(v);
            // std::cout << "planes: " << plane.size() << std::endl;
            for (auto p: plane) {
                int val = m.virtualValence(v, true, p);
                if (val == 2) {
                    // std::cout << "Fixing Boundary Doublet" << std::endl;
                    // m.Clear();            
                    // vInfo info_v(mesh, v.id, &m);
                    // PerformOperation(Operation("Rotate", v.id, std::vector<size_t>{}, false), &m);
                    // info_v = vInfo(mesh, v.id, &m);
                    // for (auto fid: info_v.fids()) {
                    //     auto& f = m.getFace(fid);
                    //     int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
                    //     PerformOperation(Operation("Collapse", v.id, std::vector<size_t>{f.Vids.at((idx+2)%f.Vids.size()), v.id}, false), &m);
                    // }
                    // m.Update();
                    // res = true;
                    // break;
                } else if (val == 3) {
                    // continue;
                    /*std::cout << "Fixing val 3 " << v.id << std::endl;
                    m.Clear();            
                    auto& f = m.getFace(p[0]);
                    int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
                    // std::cout << "Let me pass" << std::endl;
                    if (!PerformOperation(Operation("Collapse", f.Vids.at((idx+1)%f.Vids.size()), std::vector<size_t>{f.Vids.at((idx+1)%f.Vids.size()), f.Vids.at((idx+3)%f.Vids.size())}, false), &m)) break;
                    // std::cout << "After Performing operation" << std::endl;
                    // std::cout << "v nfids: " << v.N_Fids.size() << std::endl;
                    // vInfo info_v(mesh, v.id, &m);
                    // std::cout << "got info v" << std::endl;
                    // if (info_v.fids().size() != 2) break;
                    // auto& toCollapse = m.getFace(info_v.fids().at(0));
                    // int toCollapse_idx = std::distance(toCollapse.Vids.begin(), std::find(toCollapse.Vids.begin(), toCollapse.Vids.end(), v.id));
                    // std::vector<size_t> path = {f.Vids.at((idx+1)%f.Vids.size()), f.Vids.at((idx+3)%f.Vids.size())};
                    // std::cout << "here" << std::endl;
                    // PerformOperation(Operation("Collapse", v.id, std::vector<size_t>{v.id, toCollapse.Vids.at((toCollapse_idx+2)%toCollapse.Vids.size())}, false), &m);
                    // PrototypeSaveSeparatrices(std::vector<std::vector<size_t>>{path}, "test");
                    // std::cout << "m fmap: " << m.fmap.size() << " m vmap; " << m.vmap.size() << std::endl;
                    // std::cout << "setting vertex feature: " << v.id << "(" << v.N_Fids.size() << ")" << std::endl;
                    m.getVertex(v.id).type = REGULAR;
                    m.getVertex(v.id).isBoundary = false;
                    // std::cout << "regular: " << (v.type == REGULAR) <<  " feature: " << (v.type == FEATURE) << " isBoundary: " << v.isBoundary << std::endl;
                    m.Update();
                    // std::cout << "v new fids: " << v.N_Fids.size() << std::endl;
                    // std::cout << "regular: " << (v.type == REGULAR) <<  " feature: " << (v.type == FEATURE) << " isBoundary: " << v.isBoundary << std::endl;
                    // PrototypeSaveSeparatrices(std::vector<std::vector<size_t>>{}, "test2");
                    res = true;
                    break;
                    // std::cout << "Fixed val 3" << std::endl;
                    // if (++it == iters) return res;*/
                } else if (val > 4) {
                    // continue;
                    // std::cout << "Fixing val > 4: " << v.id << std::endl;
                    m.Clear();
                    vInfo info_v(mesh, v.id, &m);
                    for (auto nvid: info_v.vids()) {
                        auto& nv = m.getVertex(nvid);
                        if (nv.isBoundary || nv.type == FEATURE || nv.N_Fids.size() == 2) continue;
                        // std::cout << "nv " << nv.id << "(" << nv.N_Fids.size() << ")" << std::endl;
                        vInfo info_nv(mesh, nv.id, &m);
                        std::vector<size_t> edge = {v.id, nv.id};
                        // std::cout << "edge: " << edge.at(0) << " " << edge.at(1) << std::endl;
                        if (plane.size() > 2 && PerformOperation(Operation("Flip", v.id, edge, false), &m)) {
                            // std::cout << "Fixing corner" << std::endl;
                            // PerformOperation(Operation("Flip", v.id, edge, false), &m);
                            m.Update();
                            m.Clear();
                            res = true;
                            // if (++it == iters) return res;
                            break;
                        }
                        auto faces = mu->GetIntersection(info_v.fids(), info_nv.fids());
                        // std::cout << "faces size: " << faces.size() << std::endl;
                        if (faces.size() != 2) continue;
                        auto& f1 = m.getFace(faces.at(0));
                        auto& f2 = m.getFace(faces.at(1));
                        // std::cout << "f1: "; for (auto fvid: f1.Vids) std::cout << fvid << " "; std::cout << std::endl;
                        // std::cout << "f2: "; for (auto fvid: f2.Vids) std::cout << fvid << " "; std::cout << std::endl;
                        std::vector<size_t> edgeccw;
                        std::vector<size_t> edgecw;
                        int idx_v_f1 = std::distance(f1.Vids.begin(), std::find(f1.Vids.begin(), f1.Vids.end(), v.id));
                        int idx_v_f2 = std::distance(f2.Vids.begin(), std::find(f2.Vids.begin(), f2.Vids.end(), v.id));
                        if (f1.Vids.at((idx_v_f1+1)%4) == nv.id) {
                            edgeccw = {f1.Vids.at((idx_v_f1+2)%4), f2.Vids.at((idx_v_f2+1)%4)};
                            edgecw = {f1.Vids.at((idx_v_f1+3)%4), f2.Vids.at((idx_v_f2+2)%4)};
                        } else {
                            edgeccw = {f1.Vids.at((idx_v_f1+1)%4), f2.Vids.at((idx_v_f2+2)%4)};
                            edgecw = {f1.Vids.at((idx_v_f1+2)%4), f2.Vids.at((idx_v_f2+3)%4)};
                        }
                        // std::cout << "edgeccw: " << edgeccw.at(0) << " " << edgeccw.at(1) << std::endl;
                        // std::cout << "edgecw: " << edgecw.at(0) << " " << edgecw.at(1) << std::endl;
                        auto checkEdge = [&] (std::vector<size_t> edge) {
                            auto& v1 = m.getVertex(edge.at(0));
                            auto& v2 = m.getVertex(edge.at(1));
                            // std::cout << "v1 boundary: " << v1.isBoundary << " " << v1.type << std::endl;
                            // std::cout << "v2 boundary: " << v2.isBoundary << " " << v2.type << std::endl;
                            if (v1.isBoundary || v2.isBoundary || v1.type == FEATURE || v2.type == FEATURE) return false;
                            return true;
                        };
                        if (checkEdge(edgeccw) && PerformOperation(Operation("Flip", v.id, edge, false), &m)) {
                            // PerformOperation(Operation("Flip", v.id, edge, false), &m);
                            m.Update();
                            res = true;
                            // std::cout << "Performed Operation" << std::endl;
                            break;
                        } else if (checkEdge(edgecw) && PerformOperation(Operation("Flip", v.id, edge, true), &m)) {
                            // PerformOperation(Operation("Flip", v.id, edge, true), &m);
                            m.Update();
                            // std::cout << "Performed operation" << std::endl;
                            res = true;
                            break;
                        }
                    }
                    if (res) break;
                }
                // std::cout << "Done with the plane" << std::endl;
            }
            // if (++it == iters) {
            //     PrototypeSaveSeparatrices(std::vector<std::vector<size_t>> {}, "test");
            //     exit(0);
            // }
        }
        while (FixDoublets());
        return res;
    };
    // renewMesh();
    // if (mainIter == 3) {
    //     PrototypeSaveSeparatrices(std::vector<std::vector<size_t>>{}, "test");
    // }
    std::vector<executablePath> paths;
    std::vector<pathQueueItem> pathQueue;

    // std::cout << "Inside TestFlips" << std::endl;
    // smoother->vtkSmoother(10);
    // PrototypeSaveSeparatrices(std::vector<std::vector<size_t>>{}, "test");
    // while (fixFeatures());
    // std::cout << "Fixed features" << std::endl;
    Smooth(nullptr);
    // return 0;
    // fixFeatures();

    // return 0;
    {
        std::mutex mtx;
        ThreadPool pool;
        std::cout << "Made pool" << std::endl;
        auto push = [&] (std::unordered_map<size_t, int>& parent, std::queue<vPath*>& q, size_t vid, int flag, int parentId) {
            if (parent.find(vid) != parent.end()) return;
            parent[vid] = parentId;
            q.push(new vPath(vid, flag));
            // auto path = {(size_t) vid, (size_t) parentId};
            // paths.push_back(path);
        };
        auto isTarget = [&] (Vertex& v, std::vector<size_t> qids = {}) {
            return (mesh->virtualValence(v, qids) == 3 || mesh->virtualValence(v, qids) == 5);
            // if (v.isBoundary || v.type == FEATURE) {
            //     vMesh m(mesh);
            //     int ideal_valence = m.getIdealValence(v.id, false);
            //     if (v.N_Fids.size() != ideal_valence) return true;
            // }
            // if ((!v.isBoundary && v.type != FEATURE) && v.N_Fids.size() != 4) return true;
            // return false;
        };
        auto countSingularities = [&] (const vInfo& info_v) {
            std::unordered_set<size_t> res;
            for (auto fid: info_v.fids()) {
                auto& f = mesh->F.at(fid);
                for (auto fvid: f.Vids) {
                    if (!mesh->V.at(fvid).isBoundary && mesh->V.at(fvid).N_Fids.size() != 4) {
                        res.insert(fvid);
                    }
                }
            }
            return res;
        };
        auto validV = [&] (Vertex& v, vPath* cur, size_t vid, int flag, std::unordered_map<size_t, int>& parent, std::queue<vPath*>& q) {
            if (parent.find(vid) != parent.end()) return false;
            // std::cout << "validating vertex: " << vid << "(" << mesh->V.at(vid).N_Fids.size() << ")" << std::endl;
            vInfo info_v(mesh, vid);
            // std::cout << "got info v" << std::endl;
            if ((mesh->V.at(cur->id).isBoundary || mesh->V.at(cur->id).type == FEATURE) && (mesh->V.at(vid).isBoundary || mesh->V.at(vid).type == FEATURE)) return false;
            // std::cout << "after boundary check" << std::endl;
            auto singularities = countSingularities(info_v);
            bool res = (singularities.size() == 0);
            // std::cout << "counted singularities and res is " << res << std::endl;
            if (singularities.size() == 1) {
                if (*begin(singularities) == v.id) return true;
                // vInfo info_s(mesh, *begin(singularities));
                auto nvids = info_v.vids(cur->id);
                // std::cout << "nvids.at(2): " << nvids.at(2) << " v.id: " << v.id << std::endl;
                // std::cout << "nvid n: " << mesh->V.at(nvids.at(2)).N_Fids.size() << std::endl;
                if (parent.find(nvids.at(2)) == parent.end() && !mesh->V.at(nvids.at(2)).isBoundary && mesh->V.at(nvids.at(2)).N_Fids.size() != 4) {
                    push(parent, q, vid, SKIP, cur->id);
                    push(parent, q, nvids.at(2), SKIP, vid);
                }
            }
            // std::cout << "res: " << res << " " << singularities.size() << std::endl;
            if (res && (flag == DEFLECT)) {
                // std::cout << "res: " << res << " " << flag << std::endl;
                auto nvids = info_v.vids(cur->id);
                for (int i = 1; i < nvids.size(); i++) {
                    if (parent.find(nvids.at(i)) != parent.end() || countSingularities(vInfo(mesh, nvids.at(i))).size() == 0) continue;
                    // auto singularities = countSingularities(vInfo(mesh, nvids.at(i)));
                    // std::cout << "singularities: " << singularities.size() << std::endl;
                    // if (singularities.size() > 0) {
                        // auto nvid = info_v.vids(nvids.at(i)).at(2);
                        // if (parent.find(nvid) != parent.end() || countSingularities(vInfo(mesh, nvid)).size() > 0) continue;
                        // std::cout << "Adding deflect vertex: " << vid << " and its branch " << nvid << std::endl;
                        push(parent, q, vid, DEFLECT, cur->id);
                        // push(parent, q, nvid, BRANCH, vid);
                        // if (!DEFLECT_SKIP) push(parent, q, nvid, BRANCH, vid);
                        res = false;
                    // }
                }
            }
            return res;
        };
        auto pathVs = [&] (Vertex& v, vPath* vp, std::unordered_map<size_t, int>& parent, std::queue<vPath*>& q) {
            // std::cout << "Getting path's next vertices for " << vp->id << "(" << mesh->V.at(vp->id).N_Fids.size() << ") " 
            // << " parent: " << parent[vp->id] << "(" << mesh->V.at(parent[vp->id]).N_Fids.size() << ")" << std::endl;
            vInfo info_v(mesh, vp->id);
            auto nvids = info_v.vids(parent[vp->id]);
            // for (int i = 1; i < nvids.size(); i++) {
            //     if (parent.find(nvids.at(i)) != parent.end()) continue;
            //     if ((v.isBoundary || v.type == FEATURE) && (mesh->V.at(nvids.at(i)).isBoundary || mesh->V.at(nvids.at(i)).type == FEATURE)) continue;
            //     push(parent, q, nvids.at(i), MAIN, vp->id);
            // }
            // return;
            if (vp->flag == MAIN) {
                // std::cout << "MAIN vertex neighbors" << std::endl;
                if (validV(v, vp, nvids.at(2), vp->flag, parent, q)) push(parent, q, nvids.at(2), MAIN, vp->id);
                int flag = (parent[vp->id] == v.id) ? MAIN_BRANCH : BRANCH;
                // std::cout << "parent: " << parent[vp->id] << " v.id: " << v.id << " flag " << flag <<  std::endl;
                if (validV(v, vp, nvids.at(1), vp->flag, parent, q)) push(parent, q, nvids.at(1), flag, vp->id);
                // std::vector<size_t> path = {(size_t)parent[vp->id], (size_t)vp->id , nvids.at(1)};
                // PrototypeSaveSeparatrices(std::vector<std::vector<size_t>>{path}, "test");
                // std::cout << "nvids size: " << nvids.size() << std::endl;
                if (nvids.size() > 3 && validV(v, vp, nvids.at(3), vp->flag, parent, q)) push(parent, q, nvids.at(3), flag, vp->id);
            } else if (vp->flag == BRANCH || vp->flag == MAIN_BRANCH) {
                // std::cout << "BRANCH vertex neighbors" << std::endl;
                if (validV(v, vp, nvids.at(2), vp->flag, parent, q)) push(parent, q, nvids.at(2), vp->flag, vp->id);
                if (vp->flag == MAIN_BRANCH) return;
                // std::cout << "checking for DEFLECT and flag is: " << vp->flag << std::endl;
                validV(v, vp, nvids.at(1), DEFLECT, parent, q);
                nvids.size() > 3 ? validV(v, vp, nvids.at(3), DEFLECT, parent, q) : 1;
            } else if (vp->flag == DEFLECT) {
                // std::cout << "DEFLECT vertex neighbors" << std::endl;
                validV(v, vp, nvids.at(2), DEFLECT, parent, q);
                if (validV(v, vp, nvids.at(1), vp->flag, parent, q)) push(parent, q, nvids.at(1), BRANCH, vp->id);
                if (nvids.size() > 3 && validV(v, vp, nvids.at(3), vp->flag, parent, q)) push(parent, q, nvids.at(3), BRANCH, vp->id);
            }
        };
        auto checkPath_ = [&] (std::vector<size_t>& path, int pth_id) {
            // return;
            iteration_it = 0;
            auto& s = mesh->V.at(path.at(0));
            vInfo info_s(mesh, s.id);
            auto nvids = info_s.vids(path.at(1));
            auto nfids = info_s.fids();
            // executablePath p;
            double max_score = -1.0;
            // if (path_it == path_idx) {
            //     PrototypeSaveSeparatrices(std::vector<std::vector<size_t>> {path}, "test");
            // }
            vMesh vmesh(mesh);
            auto m = &vmesh;
            for (int idx_ = 0; idx_ < info_s.nvids()*2; idx_++) {
                // std::cout << "at the vert start" << std::endl;
                // std::cout << "iteration_it: " << iteration_it << " idx_: " << idx_ <<
                // " path_it: " << path_it << " path_idx: " << path_idx << " singularity_it: " << singularity_it <<
                // " singularity_idx: " << singularity_idx << " mainIter: " << mainIter << " info_s.nvids size: " << info_s.nvids() << std::endl;
                iteration_it++;
                // if (path_it != path_idx || iteration_it != iteration_idx) continue;
                // if (path_it != path_idx) continue;
                bool log = path_it == path_idx && iteration_it == iteration_idx && mainIter == mainIter_idx && singularity_it == singularity_idx;
                // bool log = iteration_it == 7;
                // log = false;
                // log = true;
                // std::cout << "log: " << log << std::endl;
                if (log) std::cout << "idx_: " << idx_ << std::endl;
                // vMesh* m = new vMesh(mesh);
                m->Clear();
                // p.m = m;
                // p.path = path;
                // p.score = 1;
                int maxIdx = 0;
                if (log) std::cout << "path_it: " << path_it << " path_idx: " << path_idx << " idx: " << idx_ << std::endl;
                if (log) {
                    // m->Update();
                    // PrototypeSaveSeparatrices(std::vector<std::vector<size_t>> {path}, "test");
                    // exit(0);
                }
                if (log) std::cout << "GETTING THREE FIVE PAIR" << std::endl;
                auto tfp = [&] (bool three) {
                    tfPair tfp(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), false);
                    // auto m = vm.at(idx_);
                    int pos = idx_;
                    if (info_s.nvids() != 3 && info_s.nvids() != 5) return tfp;
                    if (three) {
                        tfp.fId = s.id;
                        if (pos < 3) {
                            auto& f = m->getFace(nfids.at((pos+1)%3));
                            int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), s.id));
                            tfp.tId = f.Vids.at((idx+2)%4);
                        } else {
                            pos = pos%3;
                            tfp.diag = true;
                            tfp.tId = nvids.at(pos);
                        }
                        auto dest = nvids.at(pos);
                        PerformOperation(Operation("Rotate", s.id, std::vector<size_t>{}, false), m);
                        std::vector<size_t> faces;
                        vInfo info_t(mesh, s.id, m);
                        for (auto fid: info_t.fids()) {
                            auto& f = m->getFace(fid);
                            if (std::find(f.Vids.begin(), f.Vids.end(), dest) != f.Vids.end()) continue;
                            int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), s.id));
                            PerformOperation(Operation("Collapse", s.id, std::vector<size_t>{s.id, f.Vids.at((idx+2)%4)}, false), m);
                        }
                    } else {
                        tfp.tId = s.id;
                        if (pos < 5) {
                            auto& f = m->getFace(nfids.at((pos+2)%5));
                            int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), s.id));
                            tfp.fId = f.Vids.at((idx+2)%4);
                        } else {
                            pos = pos%5;
                            tfp.diag = true;
                            tfp.fId = nvids.at(pos);
                        }
                        auto edge1 = {s.id, nvids.at((pos+1)%nvids.size())};
                        auto edge2 = {s.id, nvids.at((pos+4)%nvids.size())};
                        PerformOperation(Operation("Split", s.id, edge1, false), m);
                        PerformOperation(Operation("Split", s.id, edge2, false), m);
                        PerformOperation(Operation("Rotate", s.id, std::vector<size_t>{}, false), m);
                    }
                    auto setMaxIdx = [&] (size_t vid) {
                        vInfo info_t(mesh, vid, m);
                        for (auto nvid: info_t.vids()) {
                            auto it = std::find(path.begin(), path.end(), nvid);
                            if (it == path.end()) continue;
                            int idx = std::distance(path.begin(), it);
                            maxIdx = std::max(maxIdx, idx);
                        }
                    };
                    setMaxIdx(tfp.tId); setMaxIdx(tfp.fId); setMaxIdx(nvids.at(pos));
                    return tfp;
                }(info_s.nvids() == 3);
                if (tfp.tId == std::numeric_limits<size_t>::max() || tfp.fId == std::numeric_limits<size_t>::max()) continue;
                if (log) std::cout << "UPDATING PATH" << std::endl;
                auto nPath = [&] (tfPair& tfp) {
                    std::vector<size_t> tPath(path.begin()+maxIdx, path.end());
                    std::unordered_map<size_t, int> parent;
                    std::queue<size_t> q;
                    // std::queue<size_t> q2;
                    auto skip = [&] (Vertex& v) {
                        if (v.id == tfp.tId || v.id == tfp.fId || v.N_Fids.size() == 3 || v.N_Fids.size() == 5) return true;
                        return false;
                    };
                    parent[tfp.tId] = -1; parent[tfp.fId] = -1;
                    for (auto vid: mu->GetDifference(vInfo(mesh, tfp.tId, m).vids(), std::vector<size_t>{tfp.fId})) {
                        parent[vid] = tfp.tId;
                        q.push(vid);
                        // q2.push(vid);
                    }
                    for (auto vid: mu->GetDifference(vInfo(mesh, tfp.fId, m).vids(), std::vector<size_t>{tfp.tId})) {
                        parent[vid] = tfp.fId;
                        q.push(vid);
                        // q2.push(vid);
                    }
                    // if (log){
                    //     std::cout << "q2: ";
                    //     while (!q2.empty()) {
                    //         std::cout << q2.front() << " ";
                    //         q2.pop();
                    //     }
                    //     std::cout << std::endl;
                    // }
                    while (!q.empty()) {
                        auto vid = q.front();
                        q.pop();

                        auto& v = m->getVertex(vid);
                        if (skip(v)) continue;
                        auto it = std::find(tPath.begin(), tPath.end(), vid);
                        if (it != tPath.end()) {
                            std::vector<size_t> newPath;
                            auto path_id = vid;
                            while (true) {
                                if (parent[path_id] == -1) break;
                                // CHANGE THIS
                                // if (path_id > mesh->V.size()) {
                                    // newPath.push_back(path_id-1);
                                // } else {
                                    newPath.push_back(path_id);
                                // }
                                path_id = parent[path_id];
                            }
                            std::reverse(newPath.begin(), newPath.end());
                            newPath.insert(newPath.end(), it+1, tPath.end());
                            return newPath;
                        }
                        vInfo info_v(mesh, vid, m);
                        auto nvids = info_v.vids(parent[vid]);
                        for (auto nvid: nvids) {
                            if (parent.find(nvid) != parent.end()) continue;
                            parent[nvid] = vid;
                            q.push(nvid);
                        }
                    }
                    return std::vector<size_t>{};
                }(tfp);
                if (log) {
                    // std::cout << "tfp.tId: " << tfp.tId << " tfp.fId: " << tfp.fId << std::endl;
                    // auto& t = m->getVertex(tfp.tId);
                    // auto& f = m->getVertex(tfp.fId);
                    // std::cout << "three fids: " << t.N_Fids.size() << " five fids: " << f.N_Fids.size() << std::endl;
                    // m->Update();
                    // std::cout << "path: ";
                    // for (auto id: path) {
                    //     std::cout << id << " ";
                    // }
                    // std::cout << std::endl;
                    // std::cout << "nPath: ";
                    // for (int pid = 0; pid < nPath.size(); pid++) {
                    //     std::cout << nPath.at(pid) << " ";
                    //     if (nPath.at(pid) > m->maxvid) nPath.at(pid) -= 1;
                    // }
                    // std::cout << std::endl;
                    // PrototypeSaveSeparatrices(std::vector<std::vector<size_t>> {nPath}, "test");
                }
                if (log) {
                    // m->Update();
                    // PrototypeSaveSeparatrices(std::vector<std::vector<size_t>> {path}, "test");
                    // exit(0);
                }
                auto movePair = [&] (tfPair& tfp, size_t dest) {
                    // std::cout << "inside movePair" << std::endl;
                    auto& three = m->getVertex(tfp.tId);
                    auto& five = m->getVertex(tfp.fId);
                    // if (three.N_Fids.size() != 3 || three.isBoundary || three.type == FEATURE || 
                    //     five.N_Fids.size() != 5 || five.isBoundary || five.type == FEATURE) return;
                        if (m->valence(three) != 3 || three.isBoundary || 
                        m->valence(five) != 5 || five.isBoundary) return 0;
                    // std::cout << "three: " << three.id << "(" << three.N_Fids.size() << ") " << " five: " << five.id << "(" << five.N_Fids.size() << ") dest: " << dest << std::endl;
                    if (log) {
                        std::cout << "three: " << three.id << "(" << three.N_Fids.size() << ") " << " five: " << five.id << "(" << five.N_Fids.size() << ") dest: " << dest << std::endl;
                        std::cout << "three fids: ";
                        for (auto fid: three.N_Fids) std::cout << fid << " ";
                        std::cout << std::endl;
                        for (auto fid: three.N_Fids) {
                            std::cout << "fid: " << fid << " ";
                            auto& f = m->getFace(fid);
                            std::cout << "vids: ";
                            for (auto vid: f.Vids) std::cout << vid << " ";
                            std::cout << std::endl;
                        }
                        std::cout << "five fids: ";
                        for (auto fid: five.N_Fids) std::cout << fid << " ";
                        std::cout << std::endl;
                    }
                    vInfo info_three(mesh, tfp.tId, m);
                    vInfo info_five(mesh, tfp.fId, m);
                    if (log) {
                        std::cout << "info_three fids: ";
                        for (auto fid: info_three.fids()) std::cout << fid << " ";
                        std::cout << std::endl;
                        std::cout << "info_five fids: ";
                        for (auto fid: info_five.fids()) std::cout << fid << " ";
                        std::cout << std::endl;
                    }
                    std::vector<size_t> vids = info_three.vids(five.id);
                    std::vector<size_t> fids = info_three.fids();
                    if (log) {
                        std::cout << "info_three fids: ";
                        for (auto fid: info_three.fids()) std::cout << fid << " ";
                        std::cout << std::endl;
                        std::cout << "info_five fids: ";
                        for (auto fid: info_five.fids()) std::cout << fid << " ";
                        std::cout << std::endl;
                    }
                    auto collapseEdge = [&] (size_t vid, size_t vid2) {
                        std::vector<size_t> edge = {vid, vid2};
                        for (int i = 0; i < path.size(); i++) {
                            if (path.at(i) == vid2) return std::vector<size_t>{vid2, vid};
                        }
                        // auto& t_v = m->getVertex(edge[1]);
                        // if (t_v.isBoundary && m->getIdealValence(t_v.id) != 2) return std::vector<size_t>{edge[1], edge[0]};
                        // if (t_v.type == FEATURE && m->getIdealValence(t_v.id) != 4) return std::vector<size_t>{edge[1], edge[0]};
                        return edge;
                    };

                    auto addVirtualElements = [&] (size_t vid) {
                        // std::cout << "adding virtual elements" << std::endl;
                        auto& v = m->getVertex(vid);
                        vInfo info_v(mesh, vid, m);
                        auto nvids = info_v.vids();
                        auto& v1 = m->AddVertex(v.xyz()); auto& v2 = m->AddVertex(v.xyz()); auto& v3 = m->AddVertex(v.xyz());
                        v1.isBoundary = true; v2.isBoundary = true; v3.isBoundary = true;
                        // std::cout << "added vertices: " << v1.id << " " << v2.id << " " << v3.id << std::endl;
                        auto& f1 = m->AddFace({vid, v2.id, v1.id, nvids.at(0)});
                        auto& f2 = m->AddFace({vid, nvids.at(info_v.nvids()-1), v3.id, v2.id});
                        // std::cout << "faces: " << f1.id << " " << f2.id << std::endl;
                        auto addNFace = [&] (Face& f) {
                            // std::cout << "face: " << f.id << " vids: ";
                            for (auto fvid: f.Vids) {
                                // std::cout << fvid << " ";
                                m->setVertex(fvid);
                                auto& fv = m->getVertex(fvid);
                                mu->AddContents(fv.N_Fids, std::vector<size_t>{f.id});
                            }
                            // std::cout << std::endl;
                        };
                        addNFace(f1); addNFace(f2);
                        return std::vector<size_t>{v1.id, v2.id, v3.id};
                    };

                    auto removeVirtualElements = [&] (std::vector<size_t> vids) {
                        // std::cout << "removing virtual elements" << std::endl;
                        auto& v = m->getVertex(vids.at(1));
                        // std::cout << "virtual v: " << v.id << std::endl;
                        // std::cout << "N_Fids: " << v.N_Fids.size() << std::endl;
                        // std::cout << "v.N_Fids: " << v.N_Fids.at(0) << " " << v.N_Fids.at(1) << std::endl; 
                        vInfo info_v(mesh, vids.at(1), m);
                        auto nfids = info_v.fids();
                        // std::cout << "nfids: ";
                        // for (auto fid: nfids) std::cout << fid << " ";
                        // std::cout << std::endl;
                        for (auto fid: nfids) {
                            // std::cout << "removing face: " << fid << std::endl;
                            auto& f = m->getFace(fid);
                            for (auto fvid: f.Vids) {
                                auto& fv = m->getVertex(fvid);
                                mu->UpdateContents(fv.N_Fids, std::vector<size_t>{f.id});
                            }
                            m->fmap.erase(fid);
                            // std::cout << "removed " << fid << std::endl;
                        }
                        m->vmap.erase(vids.at(0)); m->vmap.erase(vids.at(1)); m->vmap.erase(vids.at(2));
                        // std::cout << "removed: ";
                        // for (auto vid: vids) std::cout << vid << " ";
                        // std::cout << std::endl;
                    };
                    
                    if (tfp.diag) {
                        if (log) std::cout << "Moving diagonal pair" << std::endl;
                        if (log) {
                            std::cout << "info_three fids: ";
                            for (auto fid: info_three.fids()) std::cout << fid << " ";
                            std::cout << std::endl;
                            std::cout << "info_five fids: ";
                            for (auto fid: info_five.fids()) std::cout << fid << " ";
                            std::cout << std::endl;
                        }
                        size_t face_id = mu->GetIntersection(info_three.fids(), info_five.fids()).at(0);
                        if (log) std::cout << "face_id: " << face_id << std::endl;
                        auto& f = m->getFace(face_id);
                        // std::cout << "face id: " << face_id << std::endl;
                        std::vector<size_t> fids = info_five.fids(face_id);
                        std::vector<size_t> vids = info_five.vids();
                        // std::cout << "fids: " << fids.size() << " vids: " << vids.size() << std::endl;
                        // std::cout << "vids: ";
                        // for (auto vid: vids) std::cout << vid << " ";
                        // std::cout << std::endl;
                        // if (vids.at(3) != dest && (m->getVertex(dest).isBoundary || m->getVertex(dest).type == FEATURE)) return;
                        if (log) std::cout << "vids: " << vids.size() << " info_five nvids: " << info_five.nvids() << std::endl; 
                        if (vids.at(0) == dest || vids.at(1) == dest || vids.at(2) == dest || vids.at(4) == dest ||
                        (info_five.nvids() == 6 && vids.at(5) == dest)) {
                            if (m->getVertex(dest).isBoundary || m->getVertex(dest).type == FEATURE) return 0;
                            if (log) std::cout << "performing flip" << std::endl;
                            int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), tfp.fId));
                            bool clockwise = (vids.at(1) == dest || vids.at(2) == dest);
                            std::vector<size_t> edge = {tfp.fId, dest};
                            if (vids.at(2) == dest) edge.at(1) = vids.at(1);
                            if (vids.at(4) == dest || (info_five.nvids() == 6 && vids.at(5) == dest)) edge.at(1) = vids.at(0);
                            // std::cout << "moving pair " << (clockwise ? "clockwise" : "counter-clockwise") << std::endl;
                            PerformOperation(Operation("Flip", tfp.fId, edge, clockwise), m);
                            tfp.tId = clockwise ? vids.at(1) : vids.at(0); tfp.fId = clockwise ? vids.at(2) : info_five.nvids() == 6 ? vids.at(5) : vids.at(4);
                            [this, &m, &tfp, &collapseEdge, &log] () {
                                if (!m->getVertex(tfp.tId).isBoundary && m->getVertex(tfp.tId).N_Fids.size() == 2) {
                                    auto& v = m->getVertex(tfp.tId);
                                    if (log) std::cout << "Performing rotate on doublet: " << v.id << "(" << v.N_Fids.size() << ")" << std::endl;
                                    PerformOperation(Operation("Rotate", tfp.tId, std::vector<size_t>{}, false), m);
                                    vInfo info_t(mesh, tfp.tId, m);
                                    for (auto fid: info_t.fids()) {
                                        auto& f = m->getFace(fid);
                                        int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), tfp.tId));
                                        auto edge = collapseEdge(tfp.tId, f.Vids.at((idx+2)%4)); tfp.tId = edge.at(0);
                                        if (log) {
                                            std::cout << "Performing collapse on three: " << edge.at(0) << "(" << m->getVertex(edge.at(0)).N_Fids.size() << ")"
                                            << " " << edge.at(1) << "(" << m->getVertex(edge.at(1)).N_Fids.size() << ")" << std::endl;
                                        }
                                        PerformOperation(Operation("Collapse", tfp.tId, edge, false), m);
                                    }
                                    if (log) {
                                        std::cout << "vmap size: " << m->vmap.size() << " fmap size: " << m->fmap.size() << std::endl;
                                        std::cout << "vMesh vids: " << std::endl;
                                        for (auto el = m->vmap.begin(); el != m->vmap.end(); el++) {
                                            std::cout << el->first << " " << el->second.id << std::endl;
                                        }
                                        std::cout << std::endl;
                                        std::cout << "vMesh faces: " << std::endl;
                                        for (auto el = m->fmap.begin(); el != m->fmap.end(); el++) {
                                            std::cout << "face id: " << el->second.id << " face vids: ";
                                            for (auto fvid: el->second.Vids) {
                                                std::cout << fvid << " ";
                                            }
                                            std::cout << std::endl;
                                        }
                                        // exit(0);
                                    }
                                } else if (m->getVertex(tfp.fId).N_Fids.size() == 6 && !(m->getVertex(tfp.fId).isBoundary || m->getVertex(tfp.fId).type == FEATURE)) {
                                    auto& v = m->getVertex(tfp.fId);
                                    vInfo info_t(mesh, tfp.fId, m);
                                    auto fid = mu->GetIntersection(info_t.fids(), vInfo(mesh, tfp.tId, m).fids()).at(0);
                                    auto fids = info_t.fids(fid);
                                    auto vids = info_t.vids();
                                    std::vector<size_t> edge1 = {tfp.fId, vids.at(2)}; std::vector<size_t> edge2 = {tfp.fId, vids.at(5)};
                                    PerformOperation(Operation("Split", tfp.fId, edge1, false), m);
                                    PerformOperation(Operation("Split", tfp.fId, edge2, false), m);
                                    PerformOperation(Operation("Rotate", tfp.fId, std::vector<size_t>{}, false), m);
                                }
                            }();
                            return 1;
                        }
                        if (vids.at(3) == dest) {
                            if (log) std::cout << "moving up" << std::endl;
                            if ((five.isBoundary || five.type == FEATURE) && (m->getVertex(dest).isBoundary || m->getVertex(dest).type == FEATURE)) return 0;
                            std::vector<size_t> nvs = {vids.at(2), vids.at(4)};
                            for (auto id: nvs) {
                                PerformOperation(Operation("Split", tfp.fId, std::vector<size_t>{tfp.fId, id}, false), m);
                            }
                            PerformOperation(Operation("Rotate", tfp.fId, std::vector<size_t>{}, false), m);
                            tfp.tId = tfp.fId; tfp.fId = dest;
                            [this, &m, &tfp, &addVirtualElements, &removeVirtualElements] () {
                                auto& v = m->getVertex(tfp.fId);
                                // std::cout << "v nfids: " << v.N_Fids.size() << " isBoundary: " << v.isBoundary << std::endl;
                                if (v.N_Fids.size() == 6 && !(v.isBoundary || v.type == FEATURE)) {
                                    vInfo info_t(mesh, tfp.fId, m);
                                    auto fid = mu->GetIntersection(info_t.fids(), vInfo(mesh, tfp.tId, m).fids()).at(0);
                                    auto fids = info_t.fids(fid);
                                    auto vids = info_t.vids();
                                    std::vector<size_t> edge1 = {tfp.fId, vids.at(2)}; std::vector<size_t> edge2 = {tfp.fId, vids.at(5)};
                                    PerformOperation(Operation("Split", tfp.fId, edge1, false), m);
                                    PerformOperation(Operation("Split", tfp.fId, edge2, false), m);
                                    PerformOperation(Operation("Rotate", tfp.fId, std::vector<size_t>{}, false), m);
                                }
                                if (v.N_Fids.size() == 3 && v.isBoundary) {
                                    vInfo info_t(mesh, tfp.fId, m);
                                    auto t_nvids = info_t.vids();
                                    std::vector<size_t> edge_vids = {t_nvids.at(0), t_nvids.at(info_t.nvids()-1)};
                                    auto virtual_vids = addVirtualElements(tfp.fId);
                                    for (auto id: edge_vids) {
                                        PerformOperation(Operation("Split", tfp.fId, std::vector<size_t>{tfp.fId, id}, false), m);
                                    }
                                    PerformOperation(Operation("Rotate", tfp.fId, std::vector<size_t>{}, false), m);
                                    removeVirtualElements(virtual_vids);
                                }
                            }();
                            return 1;
                        }
                        fids = info_three.fids(face_id);
                        vids = info_three.vids();
                        if (vids.at(2) == dest) {
                            if (log) std::cout << "moving down" << std::endl;
                            if (log) std::cout << "three boundary: " << three.isBoundary << " feature: " << (three.type == FEATURE) << std::endl;
                            if (log) std::cout << "dest boundary: " << m->getVertex(dest).isBoundary << " feature: " << (m->getVertex(dest).type == FEATURE) << std::endl;
                            if ((three.isBoundary || three.type == FEATURE) && (m->getVertex(dest).isBoundary || m->getVertex(dest).type == FEATURE)) return 0;
                            if (log) std::cout << "moving to " << dest << std::endl;
                            if (log) {
                                auto& tfpt = m->getVertex(tfp.tId);
                                std::cout << "three: " << tfpt.id << "(" << tfpt.N_Fids.size() << ")" << std::endl;
                                for (auto id: tfpt.N_Fids) {
                                    auto& tfptf = m->getFace(id);
                                    std::cout << id << ": ";
                                    for (auto fvid: tfptf.Vids) {
                                        std::cout << fvid << " ";
                                    }
                                    std::cout << std::endl;
                                }
                            }
                            PerformOperation(Operation("Rotate", tfp.tId, std::vector<size_t>{}, false), m);
                            if (log) std::cout << "Rotated" << std::endl;
                            std::vector<size_t> faces;
                            vInfo info_t(mesh, tfp.tId, m);
                            if (log) {
                                auto& tfpt = m->getVertex(tfp.tId);
                                std::cout << "three: " << tfpt.id << "(" << tfpt.N_Fids.size() << ")" << std::endl;
                                for (auto id: tfpt.N_Fids) {
                                    auto& tfptf = m->getFace(id);
                                    std::cout << id << ": ";
                                    for (auto fvid: tfptf.Vids) {
                                        std::cout << fvid << " ";
                                    }
                                    std::cout << std::endl;
                                }
                            }
                            for (auto fid: info_t.fids()) {
                                if (log) std::cout << "fid: " << fid << std::endl;
                                auto& f = m->getFace(fid);
                                if (std::find(f.Vids.begin(), f.Vids.end(), dest) != f.Vids.end()) continue;
                                int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), tfp.tId));
                                auto edge = collapseEdge(tfp.tId, f.Vids.at((idx+2)%4)); tfp.tId = edge.at(0);
                                glm::dvec3 coords = {0.0, 0.0, 0.0};
                                for (auto e: edge) {
                                    if (corner(e, m)) coords = m->getVertex(e).xyz();
                                }
                                PerformOperation(Operation("Collapse", tfp.tId, edge, false, coords), m);
                                if (log) std::cout << "Collapsed" << std::endl;
                            }
                            tfp.fId = tfp.tId; tfp.tId = dest;
                            [this, &m, &tfp, &collapseEdge, &addVirtualElements, &removeVirtualElements, &corner, &log, &path] () {
                                if (log) std::cout << "checking stuff" << std::endl;
                                auto& v = m->getVertex(tfp.tId);
                                if (!v.isBoundary && v.N_Fids.size() == 2) {
                                    if (log) std::cout << "le bhai trouble" << std::endl;
                                    PerformOperation(Operation("Rotate", tfp.tId, std::vector<size_t>{}, false), m);
                                    if (log) std::cout << "Rotated for doublet" << std::endl;
                                    vInfo info_t(mesh, tfp.tId, m);
                                    for (auto fid: info_t.fids()) {
                                        auto& f = m->getFace(fid);
                                        int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), tfp.tId));
                                        auto edge = collapseEdge(tfp.tId, f.Vids.at((idx+2)%4)); tfp.tId = edge.at(0);
                                        glm::dvec3 coords = {0.0, 0.0, 0.0};
                                        for (auto e: edge) {
                                            if (corner(e, m)) coords = m->getVertex(e).xyz();
                                        }
                                        PerformOperation(Operation("Collapse", tfp.tId, edge, false, coords), m);
                                        if (log) std::cout << "Collapsed for doublet" << std::endl;
                                    }
                                }
                                if (v.N_Fids.size() == 1 && v.isBoundary) {
                                    if (log) std::cout << "I see" << std::endl;
                                    // glm::dvec3 coords = {0.0, 0.0, 0.0};
                                    // if (m->getIdealValence(v.id) != 2) {
                                        // std::cout << "ITS A CORNER" << std::endl;
                                        // coords = v.xyz();
                                        // std::cout << "coords: " << coords.x << " " << coords.y << " " << coords.z << std::endl;
                                        // std::cout << "glm::length(): " << glm::length(coords) << std::endl;
                                    // }
                                    auto& f = m->getFace(v.N_Fids.at(0));
                                    glm::dvec3 coords = {0.0, 0.0, 0.0};
                                    for (auto fvid: f.Vids) {
                                        if (corner(fvid, m)) {
                                            coords = m->getVertex(fvid).xyz();
                                            // std::cout << "Got Corner" << std::endl;
                                            // std::cout << "coords: " << coords.x << " " << coords.y << " " << coords.z << std::endl;
                                            break;
                                        }
                                    }
                                    auto virtual_vids = addVirtualElements(tfp.tId);
                                    PerformOperation(Operation("Rotate", tfp.tId, std::vector<size_t>{}, false), m);
                                    
                                    vInfo virtual_info(mesh, tfp.tId, m);
                                    auto virtual_fids = m->getVertex(tfp.tId).N_Fids;
                                    // int itskip = 0;
                                    for (auto fid: virtual_fids) {
                                        auto& f = m->getFace(fid);
                                        if (std::find(f.Vids.begin(), f.Vids.end(), virtual_vids.at(1)) != f.Vids.end()) continue;
                                        int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), tfp.tId));
                                        auto edge = collapseEdge(tfp.tId, f.Vids.at((idx+2)%4)); tfp.tId = edge.at(0);
                                        PerformOperation(Operation("Collapse", tfp.tId, edge, false, coords), m);
                                        // if (log && itskip == 1) {
                                        //     m->Update();
                                        //     PrototypeSaveSeparatrices(std::vector<std::vector<size_t>> {path}, "test");
                                        // }
                                        // itskip++;
                                    }
                                    removeVirtualElements(virtual_vids);
                                }
                            }();
                            return 1;
                        }
                    } else {
                        if (log) std::cout << "moving direct pair" << std::endl;
                        // std::cout << "vids: ";
                        // for (auto vid: vids) std::cout << vid << " ";
                        // std::cout << std::endl;
                        // std::cout << "five vids: ";
                        // for (auto vid: info_five.vids()) std::cout << vid << " ";
                        // std::cout << std::endl;
                        if (m->getVertex(dest).isBoundary || m->getVertex(dest).type == FEATURE) return 0;

                        if (vids.at(1) == dest || vids.at(2) == dest) {
                            if (log) std::cout << "Performing Collapse" << std::endl;
                            size_t fidx = (vids.at(1) == dest) ? fids.at(0) : fids.at(2);
                            auto& f = m->getFace(fidx);
                            int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), tfp.tId));
                            auto edge = collapseEdge(tfp.tId, f.Vids.at((idx+2)%4)); tfp.tId = edge.at(0);
                            if (m->getVertex(edge[1]).isBoundary || m->getVertex(edge[1]).type == FEATURE) return 0;
                            tfp.fId = tfp.tId; tfp.tId = dest;
                            PerformOperation(Operation("Collapse", tfp.tId, edge, false), m);
                            [this, &m, &tfp, &collapseEdge] () {
                                auto& v = m->getVertex(tfp.tId);
                                if (!v.isBoundary && v.N_Fids.size() == 2) {
                                    auto& f = m->getFace(v.N_Fids.at(0));
                                    int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), tfp.tId));
                                    auto edge = collapseEdge(tfp.tId, f.Vids.at((idx+2)%4)); tfp.tId = edge.at(0);
                                    PerformOperation(Operation("Collapse", tfp.tId, edge, false), m);
                                }
                            }();
                            return 1;
                        }
                        vids = info_five.vids(three.id);
                        fids = info_five.fids();
                        // std::cout << "vids: ";
                        // for (auto vid: vids) std::cout << vid << " ";
                        // std::cout << std::endl;
                        // std::cout << "fids: " << std::endl;
                        // for (auto fid: fids) {
                        //     std::cout << "fid: " << fid << " ";
                        //     auto& f = m->getFace(fid);
                        //     for (auto fvid: f.Vids) std::cout << fvid << " "; std::cout << std::endl;
                        // }
                        // std::cout << std::endl;
                        if (vids.at(1) == dest || vids.at(4) == dest) {
                            if (log) std::cout << "Performing Flip" << std::endl;
                            bool clockwise = (vids.at(4) == dest);
                            // std::cout << "clockwise: " << clockwise << std::endl;
                            int fidx = (vids.at(1) == dest) ? fids.at(1) : fids.at(3);
                            // std::cout << "fidx: " << fidx << std::endl;
                            auto& f = m->getFace(fidx);
                            // std::cout << "f vids: "; for (auto fvid: f.Vids) std::cout << fvid << " "; std::cout << std::endl;
                            int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), tfp.fId));
                            // std::cout << "new five: " << m->getVertex(tfp.fId).N_Fids.size() << std::endl;
                            size_t temp = f.Vids.at((idx+2)%4);
                            PerformOperation(Operation("Flip", tfp.fId, std::vector<size_t>{tfp.fId, dest}, clockwise), m);
                            tfp.tId = dest; tfp.fId = temp;
                            // std::cout << "new five: " << tfp.fId << " " << m->getVertex(tfp.fId).N_Fids.size() << std::endl;
                            [this, &m, &tfp, &collapseEdge] () {
                                auto& v = m->getVertex(tfp.tId);
                                if (!v.isBoundary && v.N_Fids.size() == 2) {
                                    auto& f = m->getFace(v.N_Fids.at(0));
                                    int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), tfp.tId));
                                    auto edge = collapseEdge(tfp.tId, f.Vids.at((idx+2)%4)); tfp.tId = edge.at(0);
                                    PerformOperation(Operation("Collapse", tfp.tId, edge, false), m);
                                }
                            }();
                            return 1;
                        } else if (vids.at(2) == dest || vids.at(3) == dest) {
                            if (log) std::cout << "Performing Split" << std::endl;
                            std::vector<size_t> edge = {tfp.fId, vids.at(2) == dest ? vids.at(1) : vids.at(4)};
                            PerformOperation(Operation("Split", tfp.fId, edge, false), m);
                            tfp.tId = m->max_vid; tfp.fId = dest;
                            [this, &m, &tfp] () {
                                auto& v = m->getVertex(tfp.fId);
                                if (v.N_Fids.size() == 6) {
                                    vInfo info_six(mesh, tfp.fId, m);
                                    auto vids = info_six.vids(tfp.tId);
                                    auto fids = info_six.fids();
                                    std::vector<size_t> edge1 = {tfp.fId, vids.at(1)};
                                    std::vector<size_t> edge2 = {tfp.fId, vids.at(2)};
                                    PerformOperation(Operation("Split", tfp.fId, edge1, false), m);
                                    PerformOperation(Operation("Flip", tfp.fId, edge2, true), m);
                                }
                            }();
                            return 1;
                        }
                    }
                };
                if (log) std::cout << "MOVING PAIR" << std::endl;
                // std::cout << "nPath: " << nPath.size() << std::endl;
                bool ok = true;
                for (int i = 0; i < nPath.size(); i++) {
                    if (log) std::cout << "i " << i << std::endl;
                    if (i+1 < nPath.size()-1 && (nPath.at(i+1) == tfp.tId || nPath.at(i+1) == tfp.fId)) continue;
                    auto vid = nPath.at(i);
                    ok = movePair(tfp, vid);
                    if (!ok) break;
                    if (log && i == iters) {
                        // MeshCaretaker c;
                        // c.saveState(mesh);
                        m->Update();
                        PrototypeSaveSeparatrices(std::vector<std::vector<size_t>> {path}, "test");
                        exit(0);
                        // c.restoreState(mesh);
                    }
                }
                if (log) std::cout << "AFTER MOVING PAIR" << std::endl;
                if (!ok) continue;
                if (log) {
                    // m->Update();
                    // PrototypeSaveSeparatrices(std::vector<std::vector<size_t>> {path}, "test");
                    // exit(0);
                }
                if (log) {
                    // m->Update();
                    // PrototypeSaveSeparatrices(std::vector<std::vector<size_t>> {path}, "test");
                    // exit(0);
                }
                if (log) std::cout << "GETTING VALENCE SCORE" << std::endl;
                auto valenceScore = [&] () {
                    std::unordered_set<size_t> qids;
                    for (auto it = m->fmap.begin(); it != m->fmap.end(); it++) qids.insert(it->second.id);
                    for (auto vid: path) {
                        auto& v = mesh->V.at(vid);
                        qids.insert(v.N_Fids.begin(), v.N_Fids.end());
                    }
                    
                    auto ideal_valence = [&] (Vertex& v, bool useVM = true) {
                        if (v.isBoundary || v.type == FEATURE) {
                            // bool log = false;
                            // if (log) log = true;
                            int ideal_valence = m->getIdealValence(v.id, useVM, log);
                            if (log) std::cout << "boundary vertex: " << v.N_Fids.size() <<  " ideal_valence: " << ideal_valence << std::endl;
                            return ideal_valence;
                        }
                        return 4;
                    };

                    auto compareValences = [&] (Vertex& v, bool useVM = true) {
                        if (v.N_Fids.empty()) return 0.0;
                        double score = 0.0;
                        // std::cout << "comparing valences" << std::endl;
                        if (v.isBoundary || v.type == FEATURE) {
                            // std::cout << "v is boundary or feature" << std::endl;
                            auto plane = m->planes(v, useVM);
                            for (auto p: plane) {
                                int val = m->valence(v, useVM, p);
                                int idealVal = m->idealValence(v, useVM, p);
                                int valence_diff = std::abs(val - idealVal);
                                if (valence_diff > 0) {
                                    score -= exp(10*valence_diff);
                                }
                            }
                        } else {
                            // std::cout << "non boundary v " << v.id << "(" << v.N_Fids.size() << ")" << std::endl;
                            // std::cout << "useVM_: " << useVM << std::endl;
                            int val = m->valence(v, useVM);
                            // std::cout << "val: " << val << std::endl;
                            // std::cout << "useVM_: " << useVM << std::endl;
                            int idealVal = m->idealValence(v, useVM);
                            // std::cout << "idealVal " << idealVal << std::endl;
                            int valence_diff = std::abs(val - idealVal);
                            if (valence_diff > 0) {
                                score -= exp(10*valence_diff);
                            }
                            // std::cout << "score " << score << std::endl;    
                        }
                        return score;
                        /*if (v.isBoundary || v.type == FEATURE) {
                            vInfo info_v(mesh, v.id, m, useVM);
                            double score = 0.0;
                            auto fids = info_v.fids();
                            for (int i = 0; i < fids.size(); i++) {
                                auto& f = m->getFace(fids.at(i), useVM);
                                int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
                                auto& fv = m->getVertex(f.Vids.at((idx+1)%4), useVM);
                                if (fv.isBoundary || fv.type == FEATURE) {
                                    fids = info_v.fids(fids.at(i));
                                    break;
                                }
                            }
                            double total_angle = 0.0;
                            int n = 0;
                            int total_valence = 0;
                            if (log) std::cout << "boundary score for: " << v.id << "(" << fids.size() << ")" << std::endl;
                            for (auto fid: fids) {
                                n++;
                                auto& f = m->getFace(fid, useVM);
                                int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
                                auto& v1 = m->getVertex(f.Vids.at((idx+1)%4), useVM);
                                auto& v2 = m->getVertex(f.Vids.at((idx+2)%4), useVM);
                                auto& v3 = m->getVertex(f.Vids.at((idx+3)%4), useVM);
                                total_angle += m->getAngle(v.id, v1.id, v2.id, useVM);
                                total_angle += m->getAngle(v.id, v2.id, v3.id, useVM);
                                if (v3.isBoundary || v3.type == FEATURE) {
                                    int idealValence = (floor(total_angle / 90.0) + floor(std::fmod(total_angle, 90.0) / 45));
                                    int valence_diff = std::abs((int) n - idealValence);
                                    total_valence += idealValence;
                                    score -= valence_diff > 0 ? exp(10*valence_diff) : 0.0;
                                    if (log) std::cout << "n " << n << " idealValence: " << idealValence << " valence_diff " << valence_diff << " score " << score << std::endl;  
                                    n = 0;
                                    total_angle = 0.0;
                                } 
                                if (useVM && !v1.isBoundary && v1.type != FEATURE && (v2.isBoundary || v2.type == FEATURE) && !v3.isBoundary && v3.type != FEATURE) {
                                    if (log) std::cout << "v2 is boundary or feature" << std::endl;
                                    score -= exp(50);
                                }
                            }
                            int total_valence_diff = std::abs((int) v.N_Fids.size() - total_valence);
                            score -= total_valence_diff > 0 ? exp(20*total_valence_diff) : 0.0;
                            return score;
                        } else {
                            // if (log) std::cout << "score for: " << v.id << "(" << v.N_Fids.size() << ")" << std::endl;
                            int valence_diff = std::abs((int) v.N_Fids.size() - 4);
                            // if (log) std::cout << "valence_diff " << valence_diff << " score " << -exp(10*valence_diff) << std::endl;  
                            return valence_diff > 0 ? -exp(10*valence_diff) : 0.0;
                        }*/
                    };

                    double score_before = [&] () {
                        if (log) std::cout << "CALCULATING SCORE BEFORE: " << std::endl;
                        double score = 0.0;
                        std::unordered_set<size_t> vids;
                        for (auto qid: qids) {
                            if (qid > m->maxfid) continue;
                            auto& q = mesh->F.at(qid);
                            for (auto vid: q.Vids) vids.insert(vid);
                        }
                        for (auto vid: vids) {
                            auto& v = mesh->V.at(vid);
                            // int valence_diff = std::abs((int) v.N_Fids.size() - ideal_valence(v, false));
                            // if (log) std::cout << "v: " << v.id << " v N_Fids: " << v.N_Fids.size() << " ideal_valence: " << ideal_valence(v, false) << " valence_diff: " << valence_diff << std::endl;
                            // std::cout << valence_diff << std::endl;
                            // if (valence_diff > 0) score -= exp(10*valence_diff);
                            score += compareValences(v, false);
                        }
                        if (log) std::cout << score << std::endl;
                        return score;
                    }();
                    double score_after = [&] () {
                        if (log) std::cout << "CALCULATING SCORE AFTER " << std::endl;
                        double score = 0.0;
                        std::unordered_set<size_t> vids;
                        for (auto qid: qids) {
                            auto& q = m->getFace(qid);
                            for (auto vid: q.Vids) vids.insert(vid);
                        }
                        for (auto vid: vids) {
                            auto& v = m->getVertex(vid);
                            // int valence_diff = std::abs((int) v.N_Fids.size() - ideal_valence(v));
                            // if (log) std::cout << "v: " << v.id << " v N_Fids: " << v.N_Fids.size() << " ideal_valence: " << ideal_valence(v) << " valence_diff: " << valence_diff << std::endl;
                            // std::cout << valence_diff << std::endl;
                            // if (valence_diff > 0) score -= exp(10*valence_diff);
                            score += compareValences(v);
                        }
                        if (log) std::cout << score << std::endl;
                        return score;
                    }();
                    
                    if (log) std::cout << "valence_score_before: " << score_before << " valence_score_after: " << score_after << std::endl;
                    return score_after - score_before;
                }();
                if (valenceScore < 0.0) continue;
                if (log) std::cout << "GETTING SINGULARITY SCORE" << std::endl;
                auto singularityScore = [&] () {
                    std::unordered_set<size_t> vids;
                    for (auto it = m->vmap.begin(); it != m->vmap.end(); it++) vids.insert(it->second.id); 
                    auto trace = [&] (size_t vid, bool useVM = false) {
                        struct node {
                            size_t id;
                            int starter;
                            size_t parent;
                            bool continuous = false;
                            double a = 0; double b = 0;
                            node(size_t id_, size_t parent_, bool continuous_ = false, double a_ = 0, double b_ = 0, int starter_ = -1) {
                                id = id_;
                                parent = parent_;
                                continuous = continuous_;
                                starter = starter_;
                                a = a_;
                                b = b_;
                            }
                        };
                        auto& v = m->getVertex(vid, useVM);
                        int valence_v = m->virtualValence(v, useVM);
                        bool skip = valence_v != 3 && valence_v != 5;
                        if (log) std::cout << "Inside trace for: " << vid << "(" << v.N_Fids.size() << ")" << std::endl;
                        double score = 0.0;
                        // if (v.isBoundary || v.type == FEATURE) score -= exp(30);
                        std::queue<node> q;
                        std::vector<bool> visited(m->max_vid+1);
                        vInfo info_v(mesh, v.id, m, useVM);
                        for (auto nvid: info_v.vids()) {
                            auto& nv = m->getVertex(nvid, useVM);
                            int starter = -1;
                            if (v.isBoundary || v.type == FEATURE) {
                                if (nv.isBoundary || nv.type == FEATURE) continue;
                                int b_valence = m->virtualValence(v, useVM, nv.N_Fids);
                                if (b_valence == 3 || b_valence == 5) skip = false;
                                if (b_valence != 3 && b_valence != 5) continue;
                                starter = nvid;
                            }
                            q.push(node(nvid, v.id, false, 1, starter));
                            visited.at(nvid) = true;
                        }
                        if (skip) return 0.0;
                        bool boundary_found = false;
                        bool opposite_found = false;
                        double min_boundary_dist = 4;
                        double min_separatrix_dist = 4;
                        double min_same_dist = 4;
                        while (!q.empty()) {
                            // std::cout << "q size: " << q.size() << std::endl;
                            std::vector<node> nodes;
                            while (!q.empty()) {
                                nodes.push_back(q.front());
                                q.pop();
                            }
                            std::vector<node> q_nodes;
                            // std::cout << "nodes: " << nodes.size() << std::endl;
                            for (auto n: nodes) {
                                // std::cout << "n.a: " << n.a << " n.b: " << n.b << std::endl;

                                // int dist = (n.a + n.b) - std::min(n.a, n.b);
                                double dist = std::max(n.a, n.b);
                                auto& qv = m->getVertex(n.id, useVM);
                                // std::cout << "qv: " << qv.id << "(" << qv.N_Fids.size() << ")" << std::endl; 
                                // if ((qv.isBoundary || qv.type == FEATURE) && dist < min_boundary_dist) {
                                if ((qv.isBoundary || qv.type == FEATURE) && !boundary_found) {
                                    if (log) std::cout << "v.id: " << v.id << " qv.id: " << qv.id << " setting min dist to: " << dist << " n.a: " << n.a << " n.b: " << n.b << std::endl;
                                    min_boundary_dist = dist;
                                    // score += -(1.0/dist+1);
                                    // score += -(1.0 - exp((min_boundary_dist-4.0) - fabs(min_boundary_dist-4.0)));
                                    double max = 5.0;
                                    // double factor = 2.0 * max;
                                    double min = min_boundary_dist < max ? min_boundary_dist : 1.0;
                                    score += -exp(max/min);
                                    // score += -(1.0 - exp((min_boundary_dist-4.0) - fabs(min_boundary_dist-4.0)));
                                    if (log) std::cout << "Found boundary and score now is: " << score << std::endl;
                                    boundary_found = true;
                                }
                                int valence = m->virtualValence(qv, useVM, m->getVertex(n.parent, useVM).N_Fids);
                                if (valence != 4 && !qv.isBoundary && qv.type != FEATURE) {
                                    bool singularity_3_or_5 = valence == 3 || valence >= 5;
                                    if (singularity_3_or_5) {
                                        if (v.isBoundary || v.type == FEATURE) {
                                            if (n.starter != -1) {
                                                valence_v = m->virtualValence(v, useVM, m->getVertex(n.starter, useVM).N_Fids);
                                            }
                                        }
                                        opposite_found = !((valence_v == 3 && valence == 3) || (valence_v == 5 && valence >= 5));
                                        if (log) {
                                            std::cout << "v.id: " << v.id << " qv.id: " << qv.id <<  " singularity dist: " << dist << " n.a: " << n.a << " n.b: " << n.b << std::endl;
                                            std::cout << "qv N_fids: " << qv.N_Fids.size() << " v N_fids: " << v.N_Fids.size() << std::endl;
                                        }
                                        // score += opposite_found ? -(1.0 - exp((dist-4.0) - fabs(dist-4.0))) : exp(1.0/dist);
                                        // score += opposite_found ? exp(1.0/dist) : -(1.0 - exp((dist-4.0) - fabs(dist-4.0)));
                                        double max = 5.0;
                                        double min = dist < max ? dist : 1.0;
                                        // score += opposite_found ? exp(5.0/dist) : -exp(max/min);
                                        score += opposite_found ? exp(5.0/dist) : -exp(5.0/dist);
                                        if (log) std::cout << "Found " << (opposite_found ? "opposite" : "same") << " singularity and score now is: " << score << std::endl;
                                        // score += same ? -exp(1.0/dist) : exp(1.0/dist);
                                        // double min = std::min(n.a, n.b);
                                        // if (n.b > 0 && min < min_separatrix_dist) {
                                        //     min_separatrix_dist = min;
                                        // }
                                        // double min = fabs((std::min(n.a, n.b) - 4.0) - fabs(std::min(n.a, n.b) - 4.0));
                                        // double min = (std::min(n.a, n.b) - 4.0) - fabs(std::min(n.a, n.b) - 4.0);
                                        // score += n.b == 0 ? 0 : (1.0 - exp(min/(min+0.1)));
                                        // score += n.b == 0 ? 0 : -(1.0 - exp(min));
                                    }
                                    if (opposite_found) {
                                        return score;
                                    } else {
                                        continue;
                                    }
                                }
                                vInfo info_qv(mesh, qv.id, m, useVM);
                                if (qv.isBoundary || qv.type == FEATURE || info_qv.nvids() != 4) continue;
                                auto qv_nvids = info_qv.vids(n.parent);
                                if (!n.continuous) {
                                    // std::cout << "n is not contiuous" << std::endl;
                                    // std::cout << "n.a: " << n.a << " n.a+1: " << n.a+1 << std::endl;
                                    q_nodes.push_back(node(qv_nvids.at(1), n.id, true, n.a, n.b+1, n.starter));
                                    q_nodes.push_back(node(qv_nvids.at(2), n.id, false, n.a+1, 0, n.starter));
                                    q_nodes.push_back(node(qv_nvids.at(3), n.id, true, n.a, n.b+1, n.starter));
                                } else {
                                    // std::cout << "n is contiuous" << std::endl;
                                    q_nodes.push_back(node(qv_nvids.at(2), n.id, true, n.a, n.b+1, n.starter));
                                }
                            }
                            for (auto n: q_nodes) {
                                auto& nv = m->getVertex(n.id, useVM);
                                // std::cout << "checking visit status" << std::endl;
                                // std::cout << "nv: " << nv.id << "(" << nv.N_Fids.size() << ")" << std::endl;
                                // std::cout << "m max_vid: " << m->max_vid << std::endl;
                                if (visited.at(nv.id)) continue;
                                // std::cout << "checked visit status" << std::endl;
                                // if (nv.isBoundary || nv.type == FEATURE) continue;
                                visited.at(nv.id) = true;
                                q.push(n);
                            }
                        }
                        // if (min_boundary_dist < min_separatrix_dist) {
                            // score += -(1.0 - exp((min_boundary_dist-4.0) - fabs(min_boundary_dist-4.0)));
                        // } else {
                            // score += -(1.0 - exp((min_separatrix_dist-4.0) - fabs(min_separatrix_dist-4.0)));
                        // }
                        return score;
                    };

                    auto ideal_valence = [&] (Vertex& v, bool useVM = true) {
                        if (v.isBoundary || v.type == FEATURE) {
                            // bool log = false;
                            // if (log) log = true;
                            int ideal_valence = m->getIdealValence(v.id, useVM, log);
                            // if (log) std::cout << "boundary vertex: " << v.N_Fids.size() <<  " ideal_valence: " << ideal_valence << std::endl;
                            return ideal_valence;
                        }
                        return 4;
                    };

                    auto score_before = [&] () {
                        if (log) std::cout << "CALCULATING SINGULARITY SCORE BEFORE" << std::endl;
                        double score = 0.0;
                        for (auto vid: vids) {
                            if (vid > m->maxvid) continue;
                            auto& v = m->getVertex(vid, false);
                            if (v.N_Fids.empty()) continue;
                            // bool skip = true;
                            // int valence = m->virtualValence(v, false);
                            if (v.isBoundary || v.type == FEATURE) {
                                auto plane = m->planes(v, false);
                                for (auto p: plane) {
                                    if (m->virtualValence(v, false, p) != 4) {
                                        score -= exp(30);
                                        // skip = false;
                                    }
                                }
                            }
                            // if (valence != 3 && valence != 5) skip = continue;
                            // if (v.N_Fids.size() == ideal_valence(v, false)) continue; 
                            // vInfo info_v(mesh, vid, m, false);
                            // if (info_v.nfids() != 3 && info_v.nfids() != 5) continue;
                            score += trace(vid);
                        }
                        if (log) std::cout << "SINGULARITY SCORE BEFORE: " << score << std::endl;
                        return score;
                    }();
                    auto score_after = [&] () {
                        if (log) std::cout << "CALCULATING SINGULARITY SCORE AFTER" << std::endl;
                        double score = 0.0;
                        for (auto vid: vids) {
                            auto& v = m->getVertex(vid);
                            if (v.N_Fids.empty()) continue;
                            // if (v.N_Fids.size() == ideal_valence(v)) continue; 
                            // bool skip = true;
                            if (v.isBoundary || v.type == FEATURE) {
                                auto plane = m->planes(v);
                                for (auto p: plane) {
                                    if (m->virtualValence(v, true, p) != 4) {
                                        score -= exp(30);
                                        // skip = false;
                                    }
                                }
                            }
                            // else if (m->virtualValence(v) != 4) {
                            //     skip = false;
                            // }
                            // if (skip) continue;
                            // vInfo info_v(mesh, vid, m);
                            // if (info_v.nfids() != 3 && info_v.nfids() != 5) continue;
                            score += trace(vid, true);
                        }
                        if (log) std::cout << "SINGULARITY SCORE AFTER: " << score << std::endl;
                        return score;
                    }();

                    // std::cout << "singularity score before: " << score_before << " singularity score after: " << score_after << std::endl;
                    return score_after - score_before;
                }();
                if (!(valenceScore > 0.0) && singularityScore < 0.0) continue;
                if (log) std::cout << "SMOOTHING" << std::endl;
                Smooth(m);
                if (log) std::cout << "GETTING ELEMENT SCORE" << std::endl;
                auto elementScore = [&] () {
                    // struct sizeNshape {
                    //     double size = 0.0;
                    //     double shape = 0.0;
                    //     double area = 0.0;
                    //     double avg_area = 0.0;
                    // };
                    double score = 0.0;
                    // std::unordered_map<size_t, sizeNshape> sNs;
                    std::unordered_set<size_t> qids;
                    double elementChange = 1.0;
                    for (auto it = m->fmap.begin(); it != m->fmap.end(); it++) {
                        if (it->first > m->maxfid || it->second.Vids.empty()) {
                            elementChange += 1;
                            // continue;
                        }
                        if (it->second.Vids.empty()) continue;
                        qids.insert(it->second.id);
                    }
                    auto qV_arr = [&] (size_t fid, bool useVM = false) {
                        double coords[4][3];
                        auto& q = m->getFace(fid, useVM);
                        for (int i = 0; i < 4; i++) {
                            auto& v = m->getVertex(q.Vids.at(i), useVM);
                            coords[i][0] = v.x; coords[i][1] = v.y; coords[i][2] = v.z;
                        }
                        return coords;
                    };
                    double shape_score = 0.0;
                    double size_score = 0.0;
                    for (auto qid: qids) {
                        auto& q = m->getFace(qid);
                        double threshold_shape = q.threshold_shape;
                        double threshold_size = q.threshold_size;
                        double shape = 0.0;
                        double size = 0.0;
                        if (!q.Vids.empty()) {
                            auto coords = qV_arr(qid, true);
                            double area = v_quad_area(4, coords);
                            size = std::pow(std::min(area/q.avg_area, q.avg_area/area), 2);
                            // shape = v_quad_shape(4, coords);
                            if (v_quad_scaled_jacobian(4, coords) < 0) return -exp(40);
                            shape = v_quad_aspect_ratio(4, coords); 
                        }
                        // shape_score += std::max(1.0, 0.75 * threshold_shape) - shape;
                        shape_score += (1.75 * threshold_shape) - shape;
                        size_score += size - (0.5 * threshold_size);
                    
                        if (log) {
                            // std::cout << "face metrics: " << std::endl;
                            // std::cout << "area: " << area << " avg area: " << q.avg_area << std::endl;
                            // std::cout << "threshold shape: " << threshold_shape << " shape: " << shape << std::endl;
                            // std::cout << "(1.5 * threshold_shape) - shape: " << (1.5 * threshold_shape) - shape << std::endl;
                            // std::cout << "threshold size: " << threshold_size << " size: " << size << std::endl;
                            // std::cout << "size - (0.1 * threshold_size): " << size - (0.1 * threshold_size) << std::endl;
                            // std::cout << "shape score: " << shape_score << " size score: " << size_score << std::endl;                    
                        }
                    }
                    score = (shape_score + size_score) / elementChange;
                    if (log) {
                        // std::cout << "element change: " << elementChange << std::endl;
                    }
                    return score;
                }();

                double path_score = -2.0;
                if (singularityScore > 0.0 && elementScore > 0.0) path_score = (singularityScore + elementScore);
                if (path_score < 0.0 && valenceScore > 0.0) {
                // if (path_score < 0.0) {
                    path_score = valenceScore + singularityScore + elementScore; 
                } else {
                    path_score += valenceScore;
                }
                // if (mainIter == 4) {
                // std::cout << "path_it: " << path_it << " singularity_it: " << singularity_it << " mainIter: " << mainIter << std::endl;
                // std::cout << "iteration idx: " << idx_ << std::endl;
                // std::cout << "valence score: " << valenceScore << std::endl;
                // std::cout << "element score: " << elementScore << std::endl;
                // std::cout << "singularity score: " << singularityScore << std::endl;
                // std::cout << "path_score: " << path_score << std::endl;
                // std::cout << "*****************************" << std::endl;
                // }
                    
                if (log) {
                    // std::cout << "Mesh->V.size(): " << mesh->V.size() << std::endl;
                    // std::cout << "iteration idx: " << idx_ << std::endl;
                    // std::cout << "valence score: " << valenceScore << std::endl;
                    // std::cout << "element score: " << elementScore << std::endl;
                    // std::cout << "singularity score: " << singularityScore << std::endl;
                    // std::cout << "path_score: " << path_score << std::endl;
                    // std::cout << "*****************************" << std::endl;
                    // m->Update();
                    // PrototypeSaveSeparatrices(std::vector<std::vector<size_t>>{path}, "test"); 
                    // exit(0);
                }
                {
                    std::lock_guard<std::mutex> lock(mtx);
                    if (path_score > paths.at(pth_id).score) {
                        // std::cout << "iteration idx: " << idx_ << std::endl;
                        // std::cout << "valence score: " << valenceScore << std::endl;
                        // std::cout << "element score: " << elementScore << std::endl;
                        // std::cout << "singularity score: " << singularityScore << std::endl;
                        // std::cout << "path_score: " << path_score << std::endl;
                        // std::cout << "*****************************" << std::endl;    
                        if (log) std::cout << "path_score is greater than max_score" << std::endl;
                        // max_score = path_score;
                        paths.at(pth_id).m = vmesh;
                        paths.at(pth_id).score = path_score;
                        paths.at(pth_id).valence_score = valenceScore;
                        paths.at(pth_id).element_score = elementScore;
                        paths.at(pth_id).singularity_score = singularityScore;
                        paths.at(pth_id).singularity_it = singularity_it;
                        paths.at(pth_id).path_it = path_it;
                        paths.at(pth_id).iteration_it = iteration_it;
                        paths.at(pth_id).mainIter = mainIter;
                        paths.at(pth_id).path = path;
                        // for (auto pid: path) {
                        //     p.path.push_back(pid);
                        // }
                        // p.path.swap(path);
                        if (log) std::cout << "path: " << paths.at(pth_id).path.size() << " " << path.size() << std::endl;
                        // std::cout << "path: " << p.path.size() << " " << path.size() << std::endl;
                    }

                }
                if (log) std::cout << "Ending iteration" << std::endl;
                // if (log && p.m) {
                //     std::cout << "vmap size: " << p.m->vmap.size() << " fmap size: " << p.m->fmap.size() << std::endl;
                //     std::cout << "vMesh vids: " << std::endl;
                //     for (auto el = p.m->vmap.begin(); el != p.m->vmap.end(); el++) {
                //         std::cout << el->first << " " << el->second.id << std::endl;
                //     }
                //     std::cout << std::endl;
                //     std::cout << "vMesh faces: " << std::endl;
                //     for (auto el = p.m->fmap.begin(); el != p.m->fmap.end(); el++) {
                //         std::cout << "face id: " << el->second.id << " face vids: ";
                //         for (auto fvid: el->second.Vids) {
                //             std::cout << fvid << " ";
                //         }
                //         std::cout << std::endl;
                //     }
                //     // exit(0);
                // }
                // if (log) 
                // std::cout << "At the very end " << "iteration_it: " << iteration_it << " idx_: " << idx_ <<
                // " path_it: " << path_it << " path_idx: " << path_idx << " singularity_it: " << singularity_it <<
                // " singularity_idx: " << singularity_idx << " mainIter: " << mainIter << " info_s.nvids size: " << info_s.nvids() << std::endl;
            }
            // if (path_it == path_idx && iteration_it == iteration_idx) std::cout << "did all iterations" << std::endl;
            // if (p.m && p.score >= 0.0) {
            // if (p.m) {
                // std::lock_guard<std::mutex> lock(mtx);
                // std::cout << "Adding to paths" << std::endl;
                // p.path = path;
                // paths_.push(p);
                // std::cout << "paths size: " << paths_.size() << std::endl;
            // }
        };
        auto handlePathQueue = [&] (pathQueueItem item, bool terminate = false) {
            auto executeLoop = [&] () {
                std::cout << "Loop execution start " << pathQueue.size() << std::endl;
                PARALLEL_FOR_BEGIN(0, pathQueue.size()) {
                    checkPath_(pathQueue.at(i).path, pathQueue.at(i).pth_id);
                } PARALLEL_FOR_END();
                std::cout << "Loop execution end" << std::endl;
            };
            if (terminate && pathQueue.size() > 0) {
                executeLoop();
                return;
            }
            pathQueue.push_back(item);
            if (pathQueue.size() >= 10000) {
                executeLoop();
                pathQueue.clear();
            }
        };
        auto getPaths = [&] (Vertex& v, int pth_id, const std::function<bool(Vertex&, std::vector<size_t>)>& isTarget) {
            std::cout << "Getting paths for " << v.id << "(" << v.N_Fids.size() << ")" << std::endl;
            singularity_it++;
            // std::cout << "singularity idx: " << singularity_it << std::endl;
            // if (singularity_it != singularity_idx) return;
            // if (singularity_it != 4) return;
            std::vector<std::vector<size_t>> pathsToSave;
            struct Node {
                size_t id;
                bool continuous = false;
                Node(size_t id_, bool continuous_ = false) {id = id_; continuous = continuous_;}
                Node& operator=(const Node& n) {
                    id = n.id;
                    continuous = n.continuous;
                    return *this;
                }
            };
            vInfo info_v(mesh, v.id);
            auto nvids = info_v.vids();
            for (auto nvid: nvids) {
                // std::cout << "starting new direction" << std::endl;
                int boundaryCount = 0;
                int singularityCount = 0;
                int nPaths = 0;
                std::unordered_map<size_t, int> parent;
                std::queue<Node> q;
                parent[v.id] = -1;
                q.push(Node(nvid));
                parent[nvid] = v.id;
                while (!q.empty()) {
                    auto curV = q.front();
                    // std::cout << "curV: " << curV.id << ", ";
                    q.pop();
                    std::vector<size_t> qids = {};
                    if (mesh->V.at(curV.id).isBoundary || mesh->V.at(curV.id).type == FEATURE) {
                        // continue;
                        qids = mesh->V.at(parent[curV.id]).N_Fids;
                        if (mesh->V.at(parent[curV.id]).isBoundary || mesh->V.at(parent[curV.id]).type == FEATURE) {
                            qids = mesh->V.at(parent[parent[curV.id]]).N_Fids;
                        }
                    }
                    if (isTarget(mesh->V.at(curV.id), qids)) {
                        if (boundaryCount > 10) continue;
                        if (mesh->V.at(curV.id).isBoundary) {
                            boundaryCount++;
                        } else {
                            singularityCount++;
                        }
                        std::vector<size_t> path;
                        auto path_id = curV.id;
                        // std::cout << "got path, boundary count: " << boundaryCount <<  ", singularityCount: " << singularityCount << ", nPaths: " << nPaths << std::endl;
                        while (path_id != -1) {
                            path.push_back(path_id);
                            path_id = parent[path_id];
                        }
                        std::reverse(path.begin(), path.end());
                        path_it++;
                        // if (nPaths == iters) {
                            // pathsToSave.push_back(path);
                            // std::cout << "Got path: " << path_it << std::endl;
                            // checkPath_(path, pth_id);
                        // }
                        pool.submit(checkPath_, path, pth_id);
                        if (++nPaths >= 10 && singularityCount > 1) break;
                        continue;
                    }
                    vInfo info_curV(mesh, curV.id);
                    auto curV_nvids = info_curV.vids(parent[curV.id]);
                    // std::vector<Node> nodes;
                    // for (int i = 1; i < curV_nvids.size(); i++) {
                    //     nodes.push_back(Node(curV_nvids.at(i), i == 2));
                    // }
                    // std::cout << "curV_nvids: " << curV_nvids.size() << std::endl;
                    // for (auto vid: curV_nvids) {
                    //     std::cout << vid << " ";
                    // }
                    // std::cout << std::endl;
                    if (!curV.continuous) {
                        for (int i = 1; i < curV_nvids.size(); i++) {
                            if (parent.find(curV_nvids.at(i)) == parent.end()) {
                                parent[curV_nvids.at(i)] = curV.id;
                                q.push(Node(curV_nvids.at(i), i != 2));
                            }    
                        }
                        // if (parent.find(curV_nvids.at(1)) == parent.end()) {
                        //     parent[curV_nvids.at(1)] = curV.id; 
                        //     q.push(Node(curV_nvids.at(1), true));
                        // }
                        // std::cout << "q front: " << q.front().id << std::endl;
                        // if (parent.find(curV_nvids.at(2)) == parent.end()) {
                        //     parent[curV_nvids.at(2)] = curV.id; 
                        //     q.push(Node(curV_nvids.at(2)));
                        // }
                        // std::cout << "q front: " << q.front().id << std::endl;
                        // if (parent.find(curV_nvids.at(3)) == parent.end()) {
                        //     parent[curV_nvids.at(3)] = curV.id;
                        //     q.push(Node(curV_nvids.at(3), true));
                        // std::cout << "q front: " << q.front().id << std::endl;
                    } else {
                        if (parent.find(curV_nvids.at(2)) == parent.end()) {
                            parent[curV_nvids.at(2)] = curV.id;
                            if (!mesh->V.at(curV_nvids.at(2)).isBoundary && mesh->V.at(curV_nvids.at(2)).type != FEATURE) {
                                q.push(Node(curV_nvids.at(2), true));
                            }
                        }
                        // std::cout << "q front: " << q.front().id << std::endl;
                    }
                }
                break;
            }
            /*int it2 = 0;
            vInfo info_v(mesh, v.id);
            auto nvids = info_v.vids();
            std::unordered_set<size_t> oneRing;
            for (auto fid: info_v.fids()) {
                auto& f = mesh->F.at(fid);
                for (auto vid: f.Vids) {
                    if (vid == v.id) continue;
                    oneRing.insert(vid);
                }
            }
            // std::unordered_map<size_t, int> parent;
            // std::queue<vPath*> q;
            // parent[v.id] = -1;
            // for (auto nvid: nvids) {
            //     push(parent, q, nvid, MAIN, v.id);
            // }
            // ThreadPool pool();
            // executablePath p;
            for (auto nvid: nvids) {
            // for (int n = 4; n < 5; n++) {
                // auto nvid = nvids.at(n);
                int nPaths = 0;
                std::unordered_map<size_t, int> parent;
                std::queue<vPath*> q;
                parent[v.id] = -1;
                push(parent, q, nvid, MAIN, v.id);
                // std::vector<size_t> path;
                // int it2 = 0;
                int boundaryCount = 0;
                int singularityCount = 0;
                while (!q.empty()) {
                    // if (it2++ >= iters) break;
                    auto curV = q.front();
                    // path = {(size_t) curV->id, (size_t) parent[curV->id]};
                    // std::cout << "curV: " << curV->id << "(" << mesh->V.at(curV->id).N_Fids.size() << ") " << curV->flag << std::endl;
                    q.pop();
                    // std::cout << "q size: " << q.size() << std::endl;
                    if (isTarget(mesh->V.at(curV->id))) {
                        if (boundaryCount > 10) continue;
                        if (mesh->V.at(curV->id).isBoundary) {
                            boundaryCount++;
                        } else {
                            singularityCount++;
                        }
                        std::vector<size_t> path;
                        auto path_id = curV->id;
                        int it = 0;
                        while (path_id != -1) {
                            path.push_back(path_id);
                            path_id = parent[path_id];
                        }
                        std::reverse(path.begin(), path.end());
                        // int pathCrossings = 0;
                        // for (auto pvid: path) {
                        //     if (std::find(oneRing.begin(), oneRing.end(), pvid) != oneRing.end()) pathCrossings++;
                        // }
                        // if (pathCrossings > 1) continue;
                        // CheckPath(path);
                        std::cout << "Got path: " << path_it << std::endl;
                        // if (mainIter == 5) {
                            pathsToSave.push_back(path);
                            checkPath_(path, pth_id);
                            // PrototypeSaveSeparatrices(pathsToSave, "test");
                            // return;
                        // } else {
                            // while (!pool.available());
                            // {std::cout << "waiting for queue to be available" << std::endl;}
                            // handlePathQueue(pathQueueItem(path, pth_id));
                            // pool.submit(checkPath_, path, pth_id);
                        // }
                        // std::cout << "After checking path" << std::endl;
                        path_it++;
                        // std::cout << "nPaths: " << nPaths << " singularityCount: " << singularityCount << std::endl;
                        if (++nPaths >= 10 && singularityCount > 1) break;
                        // return;
                        // pool.submit(CheckPath, path);
                        // paths.push_back(path);
                        continue;
                    }
                    if (curV->flag == SKIP) {delete curV; continue;}
                    pathVs(v, curV, parent, q);
                    delete curV;
                }
                // paths.push_back(path);
            }*/
            // std::cout << "pathsToSave: " << pathsToSave.size() << std::endl;
            // PrototypeSaveSeparatrices(pathsToSave, "test");
            // if (p.score >= 0.0) {
            //     std::lock_guard<std::mutex> lock(mtx);
            //     paths_.push(p);
            // }
        };

        std::vector<size_t> singularities;
        for (auto& v: mesh->V) {
            if (isSingularity(v)) {
                singularities.push_back(v.id);
            }
        }
        
        // size_t seed = std::hash<size_t>{}();
        // size_t seed = std::hash<size_t>{}(std::chrono::high_resolution_clock::now().time_since_epoch().count() * singularities.size());
        // std::default_random_engine e{unsigned(seed)};
        // std::shuffle(singularities.begin(), singularities.end(), e);
        // int num_singularities = ceil(0.5*singularities.size());
        int num_singularities = ceil(1*singularities.size());
        paths.resize(num_singularities);
        // int num_singularities = singularities.size();
        // std::cout << "num singularities: " << num_singularities << std::endl;
        int breaker = 0;
        for (int i = 0; i < num_singularities; i++) {
            auto& v = mesh->V.at(singularities.at(i));
            getPaths(v, i, isTarget);
            // if (++breaker > iters) break;
        }
        // handlePathQueue(pathQueueItem(), true);
        /*for (auto& v: mesh->V) {
            if (v.isBoundary || v.type == FEATURE) continue;
            if (isSingularity(v)) {
                // if (it++ == 0) continue;
                // if (singularity_it++ != singularity_idx) continue;
                getPaths(v, paths, isTarget);
                // std::random_device rd;
                // std::mt19937 g(rd());
                // std::uniform_int_distribution<size_t> dist(0, paths.size()-1);
                // int idx = dist(g);
                // idx = 32;
                // auto& path = paths.at(idx);
                // PrototypeSaveSeparatrices(std::vector<std::vector<size_t>>{path}, "test");
                // PrototypeSaveSeparatrices(paths, "test");
                // return;
                // std::cout << "Chose path at idx: " << idx << std::endl;
                // std::cout << "path size: " << path.size() << std::endl;
                // for (auto path: paths) {
                // for (int i = iters; i < iters+1; i++) {
                // std::cout << "paths: " << paths.size() << std::endl;
                // for (int i = 0; i < paths.size(); i++) {
                //     std::cout << "path idx: " << i << std::endl;
                //     auto path = paths.at(i);
                //     std::cout << "path: ";
                //     for (auto id: path) std::cout << id << " ";
                //     std::cout << std::endl;
                //    std::cout << "Checking path" << std::endl;
                    // checkPath_(path);
                    // break;
                // }
                // PrototypeSaveSeparatrices(paths, "test");
                // std::cout << "Found " << paths.size() << " paths for vertex: " << v.id << std::endl;
                // break;
            }
        }*/
    }
    std::vector<std::vector<size_t>> fpaths;
    std::vector<executablePath> final_paths;
    for (auto p: paths) {
        if (p.score >= 0.0) paths_.push(p);
    }
    std::cout << "paths: " << paths_.size() << std::endl;
    std::vector<bool> visited(mesh->V.size(), false);
    // if (paths_.empty()) break;
    while (!paths_.empty()) {
        auto p = paths_.top();
        // std::cout << "Got path with score: " << p.score << std::endl;
        // fpaths.push_back(p.path);
        // PrototypeSaveSeparatrices(fpaths, "test");
        // p.m->Update();
        paths_.pop();
        // std::cout << "Popped path out" << std::endl;
        std::set<size_t> path_vids = [&] () {
            std::set<size_t> vids;
            for (auto vid: p.path) {
                auto& v = mesh->V.at(vid);
                for (auto fid: v.N_Fids) {
                    auto& f = mesh->F.at(fid);
                    vids.insert(f.Vids.begin(), f.Vids.end());
                }
            }
            // for (auto el = p.m->vmap.begin(); el != p.m->vmap.end(); el++) {
            //     if (el->second.id > p.m->maxvid) continue;
            //     vids.insert(el->second.id);
            // }
            for (auto it = p.m.vmap.begin(); it != p.m.vmap.end(); it++) {
                if (it->first > p.m.maxvid) continue;
                // std::cout << "inserting " << it->first << std::endl;
                vids.insert(it->first);
            }
            for (auto it = p.m.fmap.begin(); it != p.m.fmap.end(); it++) {
                if (it->first <= p.m.maxfid && it->second.Vids.empty()) {
                    // std::cout << "face Vids are empty: " << it->first << " " << it->second.id << " mesh F size: " << mesh->F.size() << std::endl;
                    auto& f = mesh->F.at(it->first);
                    vids.insert(f.Vids.begin(), f.Vids.end());
                } else {
                    for (auto vid: it->second.Vids) {
                        if (vid > p.m.maxvid) continue;
                        vids.insert(vid);
                    }
                }
            }
            return vids;
        }();
        // std::cout << "got path vids" << std::endl;
        bool already_taken = false;
        for (auto vid: path_vids) {
            // std::cout << "checking " << vid << std::endl;
            if (visited.at(vid)) {
                already_taken = true;
                break;
            }
        }
        // std::cout << "checked path vids" << std::endl;
            // for (auto it = p.m->vmap.begin(); it != p.m->vmap.end(); it++) {
            //     if (it->second.id > p.m->maxvid) continue;
            //     if (visited.at(it->second.id)) {
                    // std::cout << "checking vertex: " << it->second.id << " visited: " << (visited.at(it->second.id) ? " yes " : " no ") << std::endl;
                    // already_taken = true;
                    // break;
                // }
                // std::cout << "checking vertex: " << it->second.id << std::endl;
                // vInfo info_v(mesh, it->second.id, p.m);
                // for (auto id: info_v.vids()) {
                    // if (id > p.m->max_vid) continue;
                    // if (visited.at(id)) return;
                // }
                // for (auto fid: it->second.N_Fids) {
                //     auto& f = p.m->getFace(fid);
                //     for (auto vid: f.Vids) {
                //         if (vid > p.m->maxvid) continue;
                //         if (visited.at(vid)) return;
                //     }
                // }
            // }
        if (already_taken) continue;
        // std::cout << "path: " << p.path.size() << std::endl;
        // for (auto vid: p.path) {
        //     auto& v = mesh->V.at(vid);
        //     std::cout << "v: " << v.id << "(" << v.N_Fids.size() << ")" << std::endl;
        // }
        // std::cout << "path not already crossed" << std::endl;
        for (auto vid: path_vids) {
            // std::cout << "setting " << vid << std::endl;
            visited.at(vid) = true;
        }   
        // std::cout << "-----------------------" << std::endl;
            // for (auto it = p.m->vmap.begin(); it != p.m->vmap.end(); it++) {
            //     if (it->second.id > p.m->maxvid) continue;
                // std::cout << "setting vertex: " << it->second.id << std::endl;
                // visited.at(it->second.id) = true;
                // vInfo info_v(mesh, it->second.id, p.m);
                // for (auto id: info_v.vids()) {
                    // if (id > p.m->max_vid) continue;
                    // visited.at(id) = true;
                // }
                // std::cout << "visited vertex: " << it->second.id << std::endl;
                // for (auto fid: it->second.N_Fids) {
                //     auto& f = p.m->getFace(fid);
                //     for (auto vid: f.Vids) {
                //         if (vid > p.m->maxvid) continue;
                //         visited.at(vid) = true;
                //     }
                // }
            // }
        // std::cout << "inserting path with score: " << p.score << " in final paths" << std::endl;
        final_paths.push_back(p);
            // fpaths.push_back(p.path);
        // }();
        // break;
    }
    std::cout << "final paths: " << final_paths.size() << std::endl;
    // for (auto p: fpaths) {

    // }
    int i = 0;
    int diff_vids = 0;
    int diff_fids = 0;
    int new_maxvid = 0;
    int new_maxfid = 0;
    // for (auto& p: final_paths) {
    //     i++;
    for (int i = 0; i < ceil(1*final_paths.size()); i++) {    
        // if (mainIter == 0 && i == iters) break;
        // PrototypeSaveSeparatrices(std::vector<std::vector<size_t>>{p.path}, "test");
        auto& p = final_paths.at(i);
        std::cout << "i: " << i <<  " mainIter_idx: " << p.mainIter << " singularity_idx: " << p.singularity_it << " path_idx: " << p.path_it << " iteration_idx: " << p.iteration_it << std::endl;
        std::cout << "path score: " << p.score << " valence score: " << p.valence_score << " element score: " << p.element_score << " singularity score: " << p.singularity_score << std::endl;
        
        // std::cout << "vMesh vids: ";
        int least_diff = 4;
        int current_size = mesh->V.size();
        int min_id = mesh->V.size();
        for (auto el = p.m.vmap.begin(); el != p.m.vmap.end(); el++) {
            if (el->first > p.m.maxvid && el->first < min_id) {
                min_id = el->first;
            }
        }
        for (auto el = p.m.fmap.begin(); el != p.m.fmap.end(); el++) {
            for (int j = 0; j < el->second.Vids.size(); j++) {
                if (el->second.Vids.at(j) > p.m.maxvid && el->second.Vids.at(j) < min_id) {
                    min_id = el->second.Vids.at(j);
                }
            }
        }
        if (min_id < current_size) {
            least_diff += (current_size - min_id);
            std::map<size_t, Vertex> new_vmap;
            for (auto el = p.m.vmap.begin(); el != p.m.vmap.end(); el++) {
                if (el->first > p.m.maxvid) {
                    el->second.id += least_diff;
                }
                // new_vmap[el->second.id] = el->second;
                new_vmap[el->second.id] = Vertex(el->second.xyz());
                new_vmap[el->second.id].id = el->second.id; new_vmap[el->second.id].isBoundary = el->second.isBoundary; new_vmap[el->second.id].type = el->second.type;
                new_vmap[el->second.id].N_Fids = el->second.N_Fids;
            }
            for (auto el = p.m.fmap.begin(); el != p.m.fmap.end(); el++) {
                for (int j = 0; j < el->second.Vids.size(); j++) {
                    if (el->second.Vids.at(j) > p.m.maxvid) {
                        el->second.Vids.at(j) += least_diff;
                    }
                }
            }
            p.m.vmap = new_vmap;
        }
        // for (auto el = p.m->vmap.begin(); el != p.m->vmap.end(); el++) {
            // std::cout << el->second.id << " ";
        //     if (el->first > p.m->maxvid) {
        //         diff_vids++;
        //     }
        // }
        // std::cout << std::endl;
        // std::cout << "vMesh faces: " << std::endl;
        // for (auto el = p.m->fmap.begin(); el != p.m->fmap.end(); el++) {
        //     if (el->first > p.m->maxfid) {
        //         diff_fids++;
        //     }
        //     std::cout << "face id: " << el->second.id << " face vids: ";
        //     for (auto fvid: el->second.Vids) {
        //         std::cout << fvid << " ";
        //     }
        //     std::cout << std::endl;
        // }
        // if (mainIter == -1) {
        //     std::cout << "vmap size: " << p.m->vmap.size() << " fmap size: " << p.m->fmap.size() << std::endl;
        //     std::cout << "vMesh vids: " << std::endl;
        //     for (auto el = p.m->vmap.begin(); el != p.m->vmap.end(); el++) {
        //         std::cout << el->first << " " << el->second.id << std::endl;
        //     }
        //     std::cout << std::endl;
        //     std::cout << "vMesh faces: " << std::endl;
        //     for (auto el = p.m->fmap.begin(); el != p.m->fmap.end(); el++) {
        //         std::cout << "face id: " << el->second.id << " face vids: ";
        //         for (auto fvid: el->second.Vids) {
        //             std::cout << fvid << " ";
        //         }
        //         std::cout << std::endl;
        //     }
        // }
        // PrototypeSaveSeparatrices(std::vector<std::vector<size_t>>{p.path}, "test2");
        p.m.Update();
        // PrototypeSaveSeparatrices(std::vector<std::vector<size_t>>{p.path}, "test3");
        
        // i++;
    }
    // bool res = final_paths.size() > 0;
    // while (fixFeatures());
    Smooth(nullptr);
    // }
    return final_paths.size() > 0;
}

bool SemiGlobalSimplifier::PerformOperation(const Operation& op, vMesh* m) {
    bool res = false;
    auto Flip = [&] () {
        if (!op.isValid()) return;
        // std::unordered_set<size_t> F;
        auto edge = op.vids;
        std::vector<size_t> quads = mu->GetIntersection(m->getVertex(edge[0]).N_Fids, m->getVertex(edge[1]).N_Fids);
        if (quads.size() != 2) return;

        std::vector<size_t> qVerts1, qVerts2;
        for (auto id: quads) {
            auto& q = m->getFace(id);
            int idx = std::distance(q.Vids.begin(), std::find(q.Vids.begin(), q.Vids.end(), edge[0]));
            if (q.Vids.at((idx+1)%4) == edge[1]) {
                qVerts1.push_back(q.Vids.at((idx+2)%4));
                qVerts1.push_back(q.Vids.at((idx+3)%4));
            } else if (q.Vids.at((idx+3)%4) == edge[1]) {
                qVerts2.push_back(q.Vids.at((idx+1)%4));
                qVerts2.push_back(q.Vids.at((idx+2)%4));
            }
        }
        if (qVerts1.size() != 2 || qVerts2.size() != 2) return;
        Face q1, q2;
        if (op.clockwise) {
            q1 = Face({edge[0], qVerts2[0], qVerts2[1], qVerts1[1]});
            q2 = Face({edge[1], qVerts1[0], qVerts1[1], qVerts2[1]});
        } else {
            q1 = Face({edge[0], qVerts2[0], qVerts1[0], qVerts1[1]});
            q2 = Face({edge[1], qVerts1[0], qVerts2[0], qVerts2[1]});
        }


        m->setFace(quads[0]); m->setFace(quads[1]);
        // F.insert(quads[0]); F.insert(quads[1]);
        m->fmap[quads[0]].Vids = q1.Vids; m->fmap[quads[1]].Vids = q2.Vids;
        for (auto id: quads) {
            auto& q = m->getFace(id);  
            for (auto vid: q.Vids) {
                m->setVertex(vid);
                // auto& v = m->getVertex(vid);
                mu->UpdateContents(m->getVertex(vid).N_Fids, quads);
                // for (auto fid: v.N_Fids) {
                //     auto& f = m->getFace(fid);
                //     if (std::find(f.Vids.begin(), f.Vids.end(), vid) == f.Vids.end()) mu->UpdateContents(v.N_Fids, std::vector<size_t>{fid});
                // }
            }
        }

        for (auto fid: quads) {
            auto& f = m->getFace(fid);
            for (auto vid: f.Vids) {
                auto& v = m->getVertex(vid);
                if (std::find(v.N_Fids.begin(), v.N_Fids.end(), fid) == v.N_Fids.end()) mu->AddContents(v.N_Fids, std::vector<size_t>{fid});
            }
        }
        
        // std::unordered_set<size_t> V;
        // for (auto fid: F) {
        //     auto& f = m->getFace(fid);
        //     for (auto vid: f.Vids) {
        //         auto& v = m->getVertex(vid);
        //         if (v.N_Fids.empty()) continue;
        //     }
        // }

        // addVerticesToSmooth(op.vids[0], V);
        // addVerticesToSmooth(op.vids[1], V);
        // smooth(V);
        res = true;
    };

    auto Collapse = [&] () {
        if (!op.isValid()) return;
        // std::unordered_set<size_t> F;
        auto un = mu->GetIntersection(m->getVertex(op.vids[0]).N_Fids, m->getVertex(op.vids[1]).N_Fids);
        if (un.size() != 1) return;

        auto& quad = m->getFace(un.at(0));
        int featureCount = 0;
        for (auto qvid: quad.Vids) {
            auto& qv = m->getVertex(qvid);
            // if (qv.type == FEATURE && m->planes(qv).size() > 2) return;
            if (m->getVertex(qvid).type == FEATURE) featureCount++;
            if (featureCount > 2) return;
        }
        m->setFace(un.at(0));
        // F.insert(un.at(0));
        for (auto vid: quad.Vids) {
            m->setVertex(vid);
        }
        mu->AddContents(m->vmap[op.vids[0]].N_Fids, m->vmap[op.vids[1]].N_Fids);
        // auto coords = (m->vmap[op.vids[0]].xyz() + m->vmap[op.vids[1]].xyz())*0.5;
        auto coords = op.coords;
        // std::cout << "vids[0] boundary or feature? " << (m->vmap[op.vids[0]].isBoundary || m->vmap[op.vids[0]].type == FEATURE ? "yes" : "no") << std::endl;
        // std::cout << "vids[1] boundary or feature? " << (m->vmap[op.vids[1]].isBoundary || m->vmap[op.vids[1]].type == FEATURE ? "yes" : "no") << std::endl;
        // std::cout << "op.coords length: " << glm::length(op.coords) << std::endl;
        if (glm::length(op.coords) == 0) {
            if (!m->vmap[op.vids[0]].isBoundary && m->vmap[op.vids[0]].type != FEATURE && 
                (m->vmap[op.vids[1]].isBoundary || m->vmap[op.vids[1]].type == FEATURE)) {
                    coords = m->vmap[op.vids[1]].xyz();
                    // m->vmap[op.vids[0]].isBoundary = m->vmap[op.vids[1]].isBoundary;
                    // m->vmap[op.vids[0]].type = m->vmap[op.vids[1]].type;
            } else if ((!m->vmap[op.vids[1]].isBoundary && m->vmap[op.vids[1]].type != FEATURE && 
                (m->vmap[op.vids[0]].isBoundary || m->vmap[op.vids[0]].type == FEATURE))) {
                    coords = m->vmap[op.vids[0]].xyz();
            } else {
                coords = (m->vmap[op.vids[0]].xyz() + m->vmap[op.vids[1]].xyz())*0.5;
            }
        //     coords = m->vmap[op.vids[1]].isBoundary || m->vmap[op.vids[1]].type == FEATURE ? m->vmap[op.vids[1]].xyz() : (m->vmap[op.vids[0]].xyz() + m->vmap[op.vids[1]].xyz())*0.5;
        }
        if (!m->vmap[op.vids[0]].isBoundary && m->vmap[op.vids[0]].type != FEATURE && 
            (m->vmap[op.vids[1]].isBoundary || m->vmap[op.vids[1]].type == FEATURE)) {
                m->vmap[op.vids[0]].isBoundary = m->vmap[op.vids[1]].isBoundary;
                m->vmap[op.vids[0]].type = m->vmap[op.vids[1]].type;
        }
        // glm::dvec3 coords = glm::length(op.coords) == 0 ? (m->vmap[op.vids[0]].xyz() + m->vmap[op.vids[1]].xyz())*0.5 : op.coords;
        // std::cout << "setting coords: " << coords.x << " " << coords.y << " " << coords.z << std::endl;
        // if (m->vmap[op.vids[0]].isBoundary && m->getIdealValence(m->vmap[op.vids[0]].id) != 2) coords = m->vmap[op.vids[0]].xyz();
        // if (m->vmap[op.vids[1]].isBoundary && m->getIdealValence(m->vmap[op.vids[1]].id) != 2, false) coords = m->vmap[op.vids[1]].xyz();
        // if (m->vmap[op.vids[0]].type == FEATURE && m->getIdealValence(m->vmap[op.vids[0]].id) != 4) coords = m->vmap[op.vids[0]].xyz();
        // if (m->vmap[op.vids[1]].type == FEATURE && m->getIdealValence(m->vmap[op.vids[1]].id) != 4, false) coords = m->vmap[op.vids[1]].xyz();

        m->vmap[op.vids[0]].xyz(coords);
        // std::cout << "updating source vertex " << op.vids[1] << "s faces" << std::endl;
        for (auto fid: m->getVertex(op.vids[1]).N_Fids) {
            if (fid == un.at(0)) continue;
            // std::cout << "changing face " << fid << std::endl;
            m->setFace(fid);
            // F.insert(fid);
            auto& f = m->getFace(fid);
            // for (auto fvid: f.Vids) std::cout << fvid << " "; std::cout << std::endl;
            int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), op.vids[1]));
            f.Vids.at(idx) = op.vids[0];
            // for (auto fvid: f.Vids) std::cout << fvid << " "; std::cout << std::endl;
        }
        for (auto vid: quad.Vids) {
            auto& v = m->getVertex(vid);
            mu->UpdateContents(v.N_Fids, std::vector<size_t>{un.at(0)});
        }
        auto& v1 = m->getVertex(op.vids[0]);
        auto& v2 = m->getVertex(op.vids[1]);
        // if (!v1.isBoundary && v1.type != FEATURE && (v2.isBoundary || v2.type == FEATURE)) {
        //     v1.xyz(v2.xyz()); v1.isBoundary = v2.isBoundary; v1.type = v2.type;
        // } 
        v2.N_Fids.clear();
        m->fmap[un.at(0)].Vids.clear();

        // std::unordered_set<size_t> V;
        // for (auto fid: F) {
        //     auto& f = m->getFace(fid);
        //     for (auto vid: f.Vids) {
        //         auto& v = m->getVertex(vid);
        //         if (v.N_Fids.empty()) continue;
        //     }
        // }
        // addVerticesToSmooth(op.vids[0], V);
        // smooth(V);
        res = true;
    };

    auto Split = [&] () {
        if (!op.isValid()) return;
        // std::unordered_set<size_t> F;
        auto edge = op.vids;
        std::vector<size_t> quads = mu->GetIntersection(m->getVertex(edge[0]).N_Fids, m->getVertex(edge[1]).N_Fids);
        if (quads.size() != 2) return;

        double threshold_shape = 0.0;
        double threshold_size = 0.0;
        double avg_area = 0.0;
        auto& v1 = m->getVertex(edge[0]); auto& v2 = m->getVertex(edge[1]);
        size_t hash = std::hash<size_t>{}(v1.id) ^ std::hash<size_t>{}(v2.id);
        auto& nv = m->AddVertex((v1.xyz() + v2.xyz())*0.5, hash);
        nv.type = v1.type == FEATURE && v2.type == FEATURE ? FEATURE : REGULAR;
        nv.isBoundary = v1.isBoundary && v2.isBoundary ? true : false;
        mu->AddContents(nv.N_Fids, quads);
        std::vector<size_t> nVids = {edge[0], 0, nv.id, 0};
        m->setVertex(edge[0]);
        mu->UpdateContents(m->vmap[edge[0]].N_Fids, quads);
        for (auto id: quads) {
            m->setFace(id);
            // F.insert(id);
            auto& q = m->getFace(id);
            if (q.threshold_shape > threshold_shape) threshold_shape = q.threshold_shape;
            if (q.threshold_size > threshold_size) threshold_size = q.threshold_size;
            if (q.avg_area > avg_area) avg_area = q.avg_area;
            int idx = std::distance(q.Vids.begin(), std::find(q.Vids.begin(), q.Vids.end(), edge[0]));
            q.Vids.at(idx) = nv.id;
            if (q.Vids.at((idx+1)%4) == edge[1]) {
                nVids[3] = q.Vids.at((idx+3)%4);
                m->setVertex(nVids[3]);
            } else if (q.Vids.at((idx+3)%4) == edge[1]) {
                nVids[1] = q.Vids.at((idx+1)%4);
                m->setVertex(nVids[1]);
            }
        }
        auto& nF = m->AddFace(nVids, threshold_shape, threshold_size, avg_area);
        // F.insert(nF.id);
        for (auto vid: nF.Vids) {
            auto& v = m->getVertex(vid);
            if (std::find(v.N_Fids.begin(), v.N_Fids.end(), nF.id) == v.N_Fids.end()) mu->AddContents(v.N_Fids, std::vector<size_t>{nF.id});
        }
        // std::unordered_set<size_t> V;
        // for (auto fid: F) {
        //     auto& f = m->getFace(fid);
        //     for (auto vid: f.Vids) {
        //         auto& v = m->getVertex(vid);
        //         if (v.N_Fids.empty()) continue;
        //     }
        // }

        // addVerticesToSmooth(op.vids[0], V);
        // addVerticesToSmooth(op.vids[1], V);
        // smooth(V);
        res = true;
    };

    auto Rotate = [&] () {
        auto& v = m->getVertex(op.vid);
        // if (v.isBoundary || v.type == FEATURE) return;
        // std::unordered_set<size_t> F;
        vInfo info_v(mesh, op.vid, m);
        std::vector<size_t> fids = info_v.fids();
        std::vector<std::vector<size_t>> nVids(fids.size());
        for (int i = 0; i < fids.size(); i++) {
            auto& f1 = m->getFace(fids.at(i));
            auto& f2 = m->getFace(fids.at((i+1)%fids.size()));
            int idx1 = std::distance(f1.Vids.begin(), std::find(f1.Vids.begin(), f1.Vids.end(), op.vid));
            int idx2 = std::distance(f2.Vids.begin(), std::find(f2.Vids.begin(), f2.Vids.end(), op.vid));
            nVids.at(i) = {op.vid, f1.Vids.at((idx1+2)%4), f1.Vids.at((idx1+3)%4), f2.Vids.at((idx2+2)%4)};
        }
        for (int i = 0; i < fids.size(); i++) {
            m->setFace(fids.at(i));
            // F.insert(fids.at(i));
            auto& f = m->getFace(fids.at(i));
            f.Vids = nVids.at(i);
        }
        for (int i = 0; i < fids.size(); i++) {
            auto& f = m->getFace(fids.at(i));
            for (auto vid: f.Vids) {
                m->setVertex(vid);
                auto& v = m->getVertex(vid);
                if (std::find(v.N_Fids.begin(), v.N_Fids.end(), fids.at(i)) == v.N_Fids.end()) mu->AddContents(v.N_Fids, std::vector<size_t>{fids.at(i)});
                for (auto nfid: v.N_Fids) {
                    auto& nf = m->getFace(nfid);
                    if (std::find(nf.Vids.begin(), nf.Vids.end(), vid) == nf.Vids.end()) mu->UpdateContents(v.N_Fids, std::vector<size_t>{nfid});
                }
            }
        }
        
        // std::unordered_set<size_t> V;
        // for (auto fid: F) {
        //     auto& f = m->getFace(fid);
        //     for (auto vid: f.Vids) {
        //         auto& v = m->getVertex(vid);
        //         if (v.N_Fids.empty()) continue;
        //     }
        // }

        // addVerticesToSmooth(op.vid, V);
        // smooth(V);
        res = true;
    };

    std::unordered_map<std::string, std::function<void()>> op_map = {
        {"Flip", Flip},
        {"Collapse", Collapse},
        {"Split", Split},
        {"Rotate", Rotate}
    };

    op_map[op.name]();
    return res;
}



