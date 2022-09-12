#include "ThreeFivePair.h"

ThreeFivePair::ThreeFivePair() {}

ThreeFivePair::ThreeFivePair(Mesh& mesh_, MeshUtil& mu_, Smoother& smoother_, size_t threeId_, size_t fiveId_) : mesh(mesh_), mu(mu_), smoother(smoother_) {
    threeId = threeId_;
    fiveId = fiveId_;

    auto& three = mesh.V.at(threeId);
    auto& five = mesh.V.at(fiveId);
    for (auto eid: three.N_Eids) {
        auto & e = mesh.E.at(eid);
        if (e.Vids.at(0) == fiveId || e.Vids.at(1) == fiveId) {
            edgeId = e.id;
            break;
        }
    }
    AddContents(pairEdges, GetDifference(three.N_Eids, std::vector<size_t>{edgeId}));
    AddContents(pairEdges, GetDifference(five.N_Eids, std::vector<size_t>{edgeId}));
}

ThreeFivePair::~ThreeFivePair() {}

void ThreeFivePair::CheckValidity() {
    if (mesh.V.size() == 0 || mesh.F.size() == 0 || mesh.C.size() == 0) {
        std::cout << "No mesh to use for Simplification Operation." << std::endl;
        exit(0);
    }
}

bool ThreeFivePair::IsOperationValid() {
    if (mesh.V.at(threeId).N_Vids.size() != 3 || mesh.V.at(fiveId).N_Vids.size() != 5) return false;
    return true;
}

std::vector<size_t> ThreeFivePair::GetDifference(std::vector<size_t>& a, std::vector<size_t>& b) {
    return mu.GetDifference(a, b);
}

std::vector<size_t> ThreeFivePair::GetUnion(std::vector<size_t>& a, std::vector<size_t>& b) {
    return mu.GetUnion(a, b);
}

std::vector<size_t> ThreeFivePair::GetIntersection(std::vector<size_t>& a, std::vector<size_t>& b) {
    return mu.GetIntersection(a, b);
}

void ThreeFivePair::AddContents(std::vector<size_t>& a, std::vector<size_t>& b) {
    mu.AddContents(a, b);
}

void ThreeFivePair::UpdateContents(std::vector<size_t>& a, std::vector<size_t>& b) {
    mu.UpdateContents(a, b);
}

int ThreeFivePair::Move(size_t dest, bool skipCheck) {
    CheckValidity();
    auto& three = mesh.V.at(threeId);
    auto& five = mesh.V.at(fiveId);
    
    // if ((dest == three.id || dest == five.id) && five.N_Vids.size() == 6 && skipCheck) {
    if ((dest == three.id || dest == five.id) && five.N_Vids.size() == 6 && !skipCheck) {
        SplitSixSingularity();
        resolvedSingularity = fiveId;
        return -1;
    }
    if (!skipCheck && !IsOperationValid()) return -1;

    size_t vn = mesh.V.at(dest).N_Vids.size();
    // std::cout << "three: " << three.id << " " << three.N_Vids.size() << std::endl;
    // std::cout << "five: " << five.id << " " << five.N_Vids.size() << std::endl;
    auto& pairEdge = mesh.E.at(edgeId);
    // std::cout << "three nvids: " << std::endl;
    // for (auto nvid: three.N_Vids) {
    //     std::cout << nvid << " ";
    // }
    // std::cout << std::endl;
    // std::cout << "five nvids: " << std::endl;
    // for (auto nvid: five.N_Vids) {
    //     std::cout << nvid << " ";
    // }
    // std::cout << std::endl;
    
    if (std::find(three.N_Vids.begin(), three.N_Vids.end(), dest) == three.N_Vids.end()
    && std::find(five.N_Vids.begin(), five.N_Vids.end(), dest) == five.N_Vids.end()) return -1;
    for (auto eid: three.N_Eids) {
        if (eid == edgeId) continue;
        auto& e = mesh.E.at(eid);
        if (std::find(e.Vids.begin(), e.Vids.end(), dest) == e.Vids.end()) continue;
        size_t fid = GetIntersection(pairEdge.N_Fids, e.N_Fids).at(0);
        size_t vid = e.Vids.at(0) == dest ? e.Vids.at(0) : e.Vids.at(1);
        // size_t vid = dest;
        auto& f = mesh.F.at(fid);
        int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), threeId));
        if (f.Vids.at((idx+1)%f.Vids.size()) == vid && mesh.V.at(f.Vids.at((idx+2)%f.Vids.size())).type != FEATURE && !mesh.V.at(f.Vids.at((idx+2)%f.Vids.size())).isBoundary) {
            // std::cout << "Moving Upper Left" << std::endl;
            // SetResolvedSingularity(dest, refMap["Upper"+std::to_string(mesh.V.at(dest).N_Vids.size())]);
            // if (vn == 3) {
            //     SetResolvedSingularity(dest, refMap["Upper3"]);
            // } else if (vn == 5) {
            //     SetResolvedSingularity(dest, refMap["Upper5"]);
            // }
            // std::cout << "Upper left: " << f.Vids.at(0) << " " << f.Vids.at(1) << " " << f.Vids.at(2) << " " << f.Vids.at(3) << std::endl;
            // std::cout << "threeId: " << threeId << " vid: " << vid << std::endl;
            MoveUpperLeft();
            // std::cout << "Moved Upper Left" << std::endl;
            return 0;
        } else if (f.Vids.at((idx+3)%f.Vids.size()) == vid  && mesh.V.at(f.Vids.at((idx+2)%f.Vids.size())).type != FEATURE && !mesh.V.at(f.Vids.at((idx+2)%f.Vids.size())).isBoundary) {
            // std::cout << "Moving Upper Right" << std::endl;
            // SetResolvedSingularity(dest, refMap["Upper"+std::to_string(mesh.V.at(dest).N_Vids.size())]);
            // if (vn == 3) {
            //     SetResolvedSingularity(dest, refMap["Upper3"]);
            // } else if (vn == 5) {
            //     SetResolvedSingularity(dest, refMap["Upper5"]);
            // }
            // std::cout << "Upper right: " << f.Vids.at(0) << " " << f.Vids.at(1) << " " << f.Vids.at(2) << " " << f.Vids.at(3) << std::endl;
            // std::cout << "threeId: " << threeId << " vid: " << vid << std::endl;
            MoveUpperRight();
            // std::cout << "Moved Upper Right" << std::endl;
            return 0;
        }
    }
    
    for (auto eid: five.N_Eids) {
        if (eid == edgeId) continue;
        auto& e = mesh.E.at(eid);
        if (std::find(e.Vids.begin(), e.Vids.end(), dest) == e.Vids.end()) continue;
        if (GetIntersection(pairEdge.N_Fids, e.N_Fids).size() == 0) continue;
        size_t fid = GetDifference(e.N_Fids, pairEdge.N_Fids).at(0);
        size_t vid = e.Vids.at(0) == dest ? e.Vids.at(0) : e.Vids.at(1);
        // size_t vid = dest;
        auto& f = mesh.F.at(fid);
        int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), fiveId));
        if (f.Vids.at((idx+1)%f.Vids.size()) == vid && mesh.V.at(f.Vids.at((idx+2)%f.Vids.size())).type != FEATURE && !mesh.V.at(f.Vids.at((idx+2)%f.Vids.size())).isBoundary) {
            // std::cout << "Moving Left" << std::endl;
            // if (vn == 3) {
            //     SetResolvedSingularity(dest, refMap["Middle3"]);
            // }
            // std::cout << "left: " << f.Vids.at(0) << " " << f.Vids.at(1) << " " << f.Vids.at(2) << " " << f.Vids.at(3) << std::endl;
            // std::cout << "fiveId: " << fiveId << " vid: " << vid << std::endl;
            MoveLeft();
            // std::cout << "Moved Left" << std::endl;
            // if (vn == 5) {
            //     SetResolvedSingularity(dest, refMap["Middle5"]);
            // }
            return 0;
        } else if (f.Vids.at((idx+3)%f.Vids.size()) == vid && mesh.V.at(f.Vids.at((idx+2)%f.Vids.size())).type != FEATURE && !mesh.V.at(f.Vids.at((idx+2)%f.Vids.size())).isBoundary) {
            // std::cout << "Moving Right" << std::endl;
            // if (vn == 3) {
            //     SetResolvedSingularity(dest, refMap["Middle3"]);
            // }
            // std::cout << "right: " << f.Vids.at(0) << " " << f.Vids.at(1) << " " << f.Vids.at(2) << " " << f.Vids.at(3) << std::endl;
            // std::cout << "fiveId: " << fiveId << " vid: " << vid << std::endl;
            MoveRight();
            // std::cout << "Moved Right" << std::endl;
            // if (vn == 5) {
            //     SetResolvedSingularity(dest, refMap["Middle5"]);
            // }
            return 0;
        }
    }

    for (auto eid: five.N_Eids) {
        if (eid == edgeId) continue;
        auto& e = mesh.E.at(eid);
        if (GetIntersection(pairEdge.N_Fids, e.N_Fids).size() == 0) continue;
        size_t fid = GetDifference(e.N_Fids, pairEdge.N_Fids).at(0);
        auto& f = mesh.F.at(fid);
        size_t eid2 = GetDifference(GetIntersection(f.Eids, five.N_Eids), std::vector<size_t>{e.id}).at(0);
        auto& e2 = mesh.E.at(eid2);
        if (std::find(e2.Vids.begin(), e2.Vids.end(), dest) == e2.Vids.end()) continue;
        size_t vid = e2.Vids.at(0) == dest ? e2.Vids.at(0) : e2.Vids.at(1);
        // size_t vid = dest;
        int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), fiveId));
        if (f.Vids.at((idx+3)%f.Vids.size()) == vid && mesh.V.at(f.Vids.at((idx+1)%f.Vids.size())).type != FEATURE && !mesh.V.at(f.Vids.at((idx+1)%f.Vids.size())).isBoundary) {
            // std::cout << "Moving Lower Left" << std::endl;
            // std::cout << "Lower left: " << f.Vids.at(0) << " " << f.Vids.at(1) << " " << f.Vids.at(2) << " " << f.Vids.at(3) << std::endl;
            // std::cout << "fiveId: " << fiveId << " vid: " << vid << std::endl;
            MoveLowerLeft();
            // std::cout << "Moved Lower Left" << std::endl;
            // if (vn == 3) {
            //     SetResolvedSingularity(dest, refMap["Lower3"]);
            // } else if (vn == 5) {
            //     SetResolvedSingularity(dest, refMap["Lower5"]);
            // }
            // SetResolvedSingularity(dest, refMap["Lower"+std::to_string(mesh.V.at(dest).N_Vids.size())]);
            return 0;
        } else if (f.Vids.at((idx+1)%f.Vids.size()) == vid && mesh.V.at(f.Vids.at((idx+1)%f.Vids.size())).type != FEATURE && !mesh.V.at(f.Vids.at((idx+1)%f.Vids.size())).isBoundary) {
            // std::cout << "Moving Lower Right" << std::endl;
            // std::cout << "Lower right: " << f.Vids.at(0) << " " << f.Vids.at(1) << " " << f.Vids.at(2) << " " << f.Vids.at(3) << std::endl;
            // std::cout << "fiveId: " << fiveId << " vid: " << vid << std::endl;
            MoveLowerRight();
            // std::cout << "Moved Lower Right" << std::endl;
            // if (vn == 3) {
            //     SetResolvedSingularity(dest, refMap["Lower3"]);
            // } else if (vn == 5) {
            //     SetResolvedSingularity(dest, refMap["Lower5"]);
            // }
            // SetResolvedSingularity(dest, refMap["Lower"+std::to_string(mesh.V.at(dest).N_Vids.size())]);
            return 0;
        }
    }
    return -1;
}

void ThreeFivePair::SetResolvedSingularity(size_t dest, int ref) {
    auto& three = mesh.V.at(threeId);
    auto& five = mesh.V.at(fiveId);
    auto& v = mesh.V.at(dest);
    auto& e = mesh.E.at(edgeId);
    // if (v.N_Vids.size() == 4) return;
    switch(ref) {
        case 1:
            resolvedSingularity = GetDifference(v.N_Vids, std::vector<size_t>{threeId, fiveId}).at(0);
            break;
        case 2:
            resolvedSingularity = threeId;
            break;
        case 3:
            for (auto fid: e.N_Fids) {
                auto& f = mesh.F.at(fid);
                if (std::find(f.Vids.begin(), f.Vids.end(), v.id) != f.Vids.end()) {
                    resolvedSingularity = GetDifference(std::vector<size_t>(f.Vids.begin(), f.Vids.end()), std::vector<size_t>{threeId, fiveId, dest}).at(0);
                    break;
                }
            }
            break;
        case 4:
            resolvedSingularity = fiveId;
            break;
        case 5:
            resolvedSingularity = threeId;
            break;
        case 6:
            resolvedSingularity = fiveId;
            break;
    }

}

int ThreeFivePair::GetResolvedSingularity() {
    return resolvedSingularity;
}

void ThreeFivePair::MoveUpperLeft() {
    // CheckValidity();
    // if (!IsOperationValid()) return;

    auto& three = mesh.V.at(threeId);
    auto& pairEdge = mesh.E.at(edgeId);
    for (auto eid: three.N_Eids) {
        if (eid == edgeId) continue;
        auto& e = mesh.E.at(eid);
        size_t fid = GetIntersection(pairEdge.N_Fids, e.N_Fids).at(0);
        size_t vid = e.Vids.at(0) == threeId ? e.Vids.at(1) : e.Vids.at(0);
        auto& f = mesh.F.at(fid);
        int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), threeId));
        if (f.Vids.at((idx+1)%f.Vids.size()) == vid) {
            // size_t newFiveId = f.Vids.at(idx);
            fiveId = threeId;
            // if (mesh.V.at(vid).N_Vids.size() == 3) {
            //     threeId = GetDifference(mesh.V.at(vid).N_Vids, std::vector<size_t>{threeId, f.Vids.at((idx+2)%f.Vids.size())}).at(0);
            // } else {
            //     threeId = vid;
            // }
            std::shared_ptr<SimplificationOperation> dc = std::make_shared<DiagonalCollapse>(mesh, mu, smoother, f.id, idx, (idx+2)%f.Vids.size());
            dc->PerformOperation();
            threeId = vid;
            // fiveId = f.Vids.at((idx+2)%f.Vids.size());
            for (auto veid: mesh.V.at(threeId).N_Eids) {
                auto& ve = mesh.E.at(veid);
                if ((ve.Vids.at(0) == threeId && ve.Vids.at(1) == fiveId) || (ve.Vids.at(1) == threeId && ve.Vids.at(0) == fiveId)) {
                    edgeId = veid;
                }
            }
            break;
        }
    }
}

void ThreeFivePair::MoveUpperRight() {
    // CheckValidity();
    // if (!IsOperationValid()) return;
    
    auto& three = mesh.V.at(threeId);
    auto& pairEdge = mesh.E.at(edgeId);
    for (auto eid: three.N_Eids) {
        if (eid == edgeId) continue;
        auto& e = mesh.E.at(eid);
        size_t fid = GetIntersection(pairEdge.N_Fids, e.N_Fids).at(0);
        size_t vid = e.Vids.at(0) == threeId ? e.Vids.at(1) : e.Vids.at(0);
        auto& f = mesh.F.at(fid);
        int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), threeId));
        if (f.Vids.at((idx+3)%f.Vids.size()) == vid) {
            // size_t newFiveId = f.Vids.at(idx);
            fiveId = threeId;
            // if (mesh.V.at(vid).N_Vids.size() == 3) {
            //     threeId = GetDifference(mesh.V.at(vid).N_Vids, std::vector<size_t>{threeId, f.Vids.at((idx+2)%f.Vids.size())}).at(0);
            // } else {
            //     threeId = vid;
            // }
            std::shared_ptr<SimplificationOperation> dc = std::make_shared<DiagonalCollapse>(mesh, mu, smoother, f.id, idx, (idx+2)%f.Vids.size());
            dc->PerformOperation();
            threeId = vid;
            // fiveId = f.Vids.at((idx+2)%f.Vids.size());
            // fiveId = newFiveId;
            // std::cout << "three: " << mesh.V.at(threeId).N_Vids.size() << " five: " << mesh.V.at(fiveId).N_Vids.size() << std::endl;
            for (auto veid: mesh.V.at(threeId).N_Eids) {
                auto& ve = mesh.E.at(veid);
                if ((ve.Vids.at(0) == threeId && ve.Vids.at(1) == fiveId) || (ve.Vids.at(1) == threeId && ve.Vids.at(0) == fiveId)) {
                    edgeId = veid;
                }
            }
            break;
        }
    }
}

void ThreeFivePair::MoveLeft() {
    // CheckValidity();
    // if (!IsOperationValid()) return;
    
    auto& five = mesh.V.at(fiveId);
    auto& pairEdge = mesh.E.at(edgeId);
    for (auto eid: five.N_Eids) {
        if (eid == edgeId) continue;
        auto& e = mesh.E.at(eid);
        if (GetIntersection(pairEdge.N_Fids, e.N_Fids).size() == 0) continue;
        size_t fid = GetDifference(e.N_Fids, pairEdge.N_Fids).at(0);
        size_t vid = e.Vids.at(0) == fiveId ? e.Vids.at(1) : e.Vids.at(0);
        auto& f = mesh.F.at(fid);
        int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), fiveId));
        if (f.Vids.at((idx+1)%f.Vids.size()) == vid) {
            threeId = vid;
            fiveId = f.Vids.at((idx+2)%f.Vids.size());
            std::shared_ptr<SimplificationOperation> s = std::make_shared<EdgeRotation>(mesh, mu, smoother, e.id, false);
            s->PerformOperation();
            // threeId = vid;
            // fiveId = f.Vids.at((idx+2)%f.Vids.size());
            for (auto veid: mesh.V.at(threeId).N_Eids) {
                auto& ve = mesh.E.at(veid);
                if ((ve.Vids.at(0) == threeId && ve.Vids.at(1) == fiveId) || (ve.Vids.at(1) == threeId && ve.Vids.at(0) == fiveId)) {
                    edgeId = veid;
                }
            }
            break;
        }
    }
}

void ThreeFivePair::MoveRight() {
    // CheckValidity();
    // if (!IsOperationValid()) return;

    auto& five = mesh.V.at(fiveId);
    auto& pairEdge = mesh.E.at(edgeId);
    for (auto eid: five.N_Eids) {
        if (eid == edgeId) continue;
        auto& e = mesh.E.at(eid);
        if (GetIntersection(pairEdge.N_Fids, e.N_Fids).size() == 0) continue;
        size_t fid = GetDifference(e.N_Fids, pairEdge.N_Fids).at(0);
        size_t vid = e.Vids.at(0) == fiveId ? e.Vids.at(1) : e.Vids.at(0);
        auto& f = mesh.F.at(fid);
        int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), fiveId));
        if (f.Vids.at((idx+3)%f.Vids.size()) == vid) {
            threeId = vid;
            fiveId = f.Vids.at((idx+2)%f.Vids.size());
            std::shared_ptr<SimplificationOperation> s = std::make_shared<EdgeRotation>(mesh, mu, smoother, e.id, true);
            s->PerformOperation();
            for (auto veid: mesh.V.at(threeId).N_Eids) {
                auto& ve = mesh.E.at(veid);
                if ((ve.Vids.at(0) == threeId && ve.Vids.at(1) == fiveId) || (ve.Vids.at(1) == threeId && ve.Vids.at(0) == fiveId)) {
                    edgeId = veid;
                }
            }
            break;
        }
    }
}

void ThreeFivePair::MoveLowerLeft() {
    // CheckValidity();
    // if (!IsOperationValid()) return;
    
    auto& five = mesh.V.at(fiveId);
    auto& pairEdge = mesh.E.at(edgeId);
    for (auto eid: five.N_Eids) {
        if (eid == edgeId) continue;
        auto& e = mesh.E.at(eid);
        if (GetIntersection(pairEdge.N_Fids, e.N_Fids).size() == 0) continue;
        size_t fid = GetDifference(e.N_Fids, pairEdge.N_Fids).at(0);
        auto& f = mesh.F.at(fid);
        size_t eid2 = GetDifference(GetIntersection(f.Eids, five.N_Eids), std::vector<size_t>{e.id}).at(0);
        auto& e2 = mesh.E.at(eid2);
        size_t vid = e2.Vids.at(0) == fiveId ? e2.Vids.at(1) : e2.Vids.at(0);
        int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), fiveId));
        if (f.Vids.at((idx+3)%f.Vids.size()) == vid) {
            std::vector<size_t> verticesTosplit;
            std::vector<size_t> verticesToChange;
            verticesTosplit.push_back(vid);
            pairEdge.Vids.at(0) == fiveId ? verticesTosplit.push_back(pairEdge.Vids.at(1)) : verticesTosplit.push_back(pairEdge.Vids.at(0));
            e.Vids.at(0) == fiveId ? verticesToChange.push_back(e.Vids.at(1)) : verticesToChange.push_back(e.Vids.at(0));
            std::shared_ptr<SimplificationOperation> qs = std::make_shared<QuadSplit>(mesh, mu, smoother, fiveId, verticesTosplit, verticesToChange);
            qs->PerformOperation();
            threeId = mesh.V.size() - 1;
            fiveId = vid;
            for (auto veid: mesh.V.at(threeId).N_Eids) {
                auto& ve = mesh.E.at(veid);
                if ((ve.Vids.at(0) == threeId && ve.Vids.at(1) == fiveId) || (ve.Vids.at(1) == threeId && ve.Vids.at(0) == fiveId)) {
                    edgeId = veid;
                }
            }
            if (mesh.V.at(fiveId).N_Vids.size() == 6) {
                SplitSixSingularity();
            }
            break;
        }
    }
}

void ThreeFivePair::MoveLowerRight() {
    // CheckValidity();
    // if (!IsOperationValid()) return;
    
    auto& five = mesh.V.at(fiveId);
    auto& pairEdge = mesh.E.at(edgeId);
    for (auto eid: five.N_Eids) {
        if (eid == edgeId) continue;
        auto& e = mesh.E.at(eid);
        if (GetIntersection(pairEdge.N_Fids, e.N_Fids).size() == 0) continue;
        size_t fid = GetDifference(e.N_Fids, pairEdge.N_Fids).at(0);
        auto& f = mesh.F.at(fid);
        size_t eid2 = GetDifference(GetIntersection(f.Eids, five.N_Eids), std::vector<size_t>{e.id}).at(0);
        auto& e2 = mesh.E.at(eid2);
        size_t vid = e2.Vids.at(0) == fiveId ? e2.Vids.at(1) : e2.Vids.at(0);
        int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), fiveId));
        if (f.Vids.at((idx+1)%f.Vids.size()) == vid) {
            std::vector<size_t> verticesTosplit;
            std::vector<size_t> verticesToChange;
            verticesTosplit.push_back(vid);
            pairEdge.Vids.at(0) == fiveId ? verticesTosplit.push_back(pairEdge.Vids.at(1)) : verticesTosplit.push_back(pairEdge.Vids.at(0));
            e.Vids.at(0) == fiveId ? verticesToChange.push_back(e.Vids.at(1)) : verticesToChange.push_back(e.Vids.at(0));
            std::shared_ptr<SimplificationOperation> qs = std::make_shared<QuadSplit>(mesh, mu, smoother, fiveId, verticesTosplit, verticesToChange);
            qs->PerformOperation();
            threeId = mesh.V.size() - 1;
            fiveId = vid;
            for (auto veid: mesh.V.at(threeId).N_Eids) {
                auto& ve = mesh.E.at(veid);
                if ((ve.Vids.at(0) == threeId && ve.Vids.at(1) == fiveId) || (ve.Vids.at(1) == threeId && ve.Vids.at(0) == fiveId)) {
                    edgeId = veid;
                }
            }
            if (mesh.V.at(fiveId).N_Vids.size() == 6) {
                SplitSixSingularity();
            }
            break;
        }
    }
}

void ThreeFivePair::SplitSixSingularity() {
    std::cout << "Fixing six singularity" << std::endl;
    auto& v = mesh.V.at(fiveId);
    int d = (v.N_Vids.size() / 2) + 1;
    size_t mainV = threeId;
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
            if (f.Vids.at((idx+1)%f.Vids.size()) == ev) {
                if (j == 0 || j == d-1) {
                    verticesToSplit.push_back(ev);
                } else {
                    verticesToChange.push_back(ev);
                }
                startE = mu.GetDifference(mu.GetIntersection(f.Eids, v.N_Eids), std::vector<size_t>{e.id}).at(0);
            }
        }
    }
    std::shared_ptr<SimplificationOperation> qs = std::make_shared<QuadSplit>(mesh, mu, smoother, v.id, verticesToSplit, verticesToChange);
    qs->PerformOperation();
    threeId = mesh.V.at(mesh.V.size()-1).id;
    fiveId = verticesToSplit.at(1);
}

std::vector<size_t> ThreeFivePair::GetPairIds() {
    return std::vector<size_t>{threeId, fiveId};
}

