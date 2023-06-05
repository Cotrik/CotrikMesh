#include "DiagonalThreeFivePair.h"

DiagonalThreeFivePair::DiagonalThreeFivePair() {}

DiagonalThreeFivePair::DiagonalThreeFivePair(Mesh& mesh_, MeshUtil& mu_, Smoother& smoother_, size_t threeId_, size_t fiveId_) : SingularityPair(mesh_, mu_, smoother_, threeId_, fiveId_) {
    // threeId = threeId_;
    // fiveId = fiveId_;

    auto& three = mesh->V.at(threeId);
    auto& five = mesh->V.at(fiveId);
    faceId = GetIntersection(three.N_Fids, five.N_Fids).at(0);
    
    AddContents(pairEdges, three.N_Eids);
    AddContents(pairEdges, five.N_Eids);
}

DiagonalThreeFivePair::~DiagonalThreeFivePair() {}

void DiagonalThreeFivePair::SetMembers(Mesh& mesh_, MeshUtil& mu_, Smoother& smoother_, size_t threeId_, size_t fiveId_) {
    mesh = &mesh_;
    mu = &mu_;
    smoother = &smoother_;
    threeId = threeId_;
    fiveId = fiveId_;

    auto& three = mesh->V.at(threeId);
    auto& five = mesh->V.at(fiveId);
    faceId = GetIntersection(three.N_Fids, five.N_Fids).at(0);
    
    AddContents(pairEdges, three.N_Eids);
    AddContents(pairEdges, five.N_Eids);
}

/*
void DiagonalThreeFivePair::CheckValidity() {
    if (mesh == NULL) {
        std::cout << "No mesh to use for Three Five Pair." << std::endl;
        exit(0);
    }
    if (mu == NULL) {
        std::cout << "MeshUtil is not initialized for Three Five Pair." << std::endl;
        exit(0);
    }
    if (smoother == NULL) {
        std::cout << "Smoother is not intitalized for Three Five Pair." << std::endl;
        exit(0);
    }
    if (mesh->V.size() == 0 || mesh->F.size() == 0 || mesh->C.size() == 0) {
        std::cout << "No mesh to use for Three Five Pair." << std::endl;
        exit(0);
    }
}

bool DiagonalThreeFivePair::IsOperationValid() {
    if (mesh->V.at(threeId).N_Vids.size() != 3 || mesh->V.at(fiveId).N_Vids.size() != 5) return false;
    if (mesh->V.at(threeId).N_Eids.size() != 3 || mesh->V.at(fiveId).N_Eids.size() != 5) return false;
    if (mesh->V.at(threeId).N_Fids.size() != 3 || mesh->V.at(fiveId).N_Fids.size() != 5) return false;
    return true;
}

std::vector<size_t> DiagonalThreeFivePair::GetDifference(std::vector<size_t> a, std::vector<size_t> b) {
    return mu->GetDifference(a, b);
}

std::vector<size_t> DiagonalThreeFivePair::GetUnion(std::vector<size_t> a, std::vector<size_t> b) {
    return mu->GetUnion(a, b);
}

std::vector<size_t> DiagonalThreeFivePair::GetIntersection(std::vector<size_t> a, std::vector<size_t> b) {
    return mu->GetIntersection(a, b);
}

void DiagonalThreeFivePair::AddContents(std::vector<size_t>& a, std::vector<size_t> b) {
    mu->AddContents(a, b);
}

void DiagonalThreeFivePair::UpdateContents(std::vector<size_t>& a, std::vector<size_t> b) {
    mu->UpdateContents(a, b);
}*/

int DiagonalThreeFivePair::Move(size_t dest, double& delta, bool skipCheck) {
    CheckValidity();
    // if (!skipCheck && !IsOperationValid()) return -1;
    // std::cout << "In Move(): skipCheck: " << skipCheck << std::endl;
    auto& three = mesh->V.at(threeId);
    auto& five = mesh->V.at(fiveId);
    if (three.N_Vids.empty() || five.N_Vids.empty()) return -1;

    // std::cout << "dest: " << dest << std::endl;
    // std::cout << "three id: " << three.id << " three nfids: " << three.N_Fids.size() << " five id: " << five.id << " five nfids: " << five.N_Fids.size() << std::endl;
    // std::cout << "three id: " << three.id << " three nvids: " << three.N_Vids.size() << " five id: " << five.id << " five nvids: " << five.N_Vids.size() << std::endl;
    // std::cout << "three id: " << three.id << " three neids: " << three.N_Eids.size() << " five id: " << five.id << " five neids: " << five.N_Eids.size() << std::endl;
    // for (auto nvid: three.N_Vids) {
    //     auto& nv = mesh->V.at(nvid);
    //     std::cout << "nv id: " << nv.id << " nv nvids: " << nv.N_Vids.size() << std::endl;
    //     std::cout << "nv id: " << nv.id << " nv neids: " << nv.N_Eids.size() << std::endl;
    //     std::cout << "nv id: " << nv.id << " nv nfids: " << nv.N_Fids.size() << std::endl;
    // }
    // for (auto nvid: five.N_Vids) {
    //     auto& nv = mesh->V.at(nvid);
    //     std::cout << "nv id: " << nv.id << " nv nvids: " << nv.N_Vids.size() << std::endl;
    //     std::cout << "nv id: " << nv.id << " nv neids: " << nv.N_Eids.size() << std::endl;
    //     std::cout << "nv id: " << nv.id << " nv nfids: " << nv.N_Fids.size() << std::endl;
    // }
    std::vector<size_t> destFaces = mesh->V.at(dest).N_Fids;
    faceId = GetIntersection(three.N_Fids, five.N_Fids).at(0);
    {
        auto& f = mesh->F.at(faceId);
        int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), fiveId));
        size_t diagId = f.Vids.at((idx + 2) % f.Vids.size());
        if (diagId != threeId || mesh->V.at(diagId).N_Vids.size() != 3) return -1;
    }
    int newThreeId = -1;
    int newFiveId = -1;

    auto& f = mesh->F.at(faceId);
    if (five.N_Vids.size() == 5) {
        std::vector<size_t> fiveEdges = GetIntersection(f.Eids, five.N_Eids);
        std::vector<size_t> tempEdges;
        for (auto eid: fiveEdges) {
            for (auto fid: mesh->E.at(eid).N_Fids) {
                if (fid == f.id) continue;
                AddContents(tempEdges, GetDifference(GetIntersection(mesh->F.at(fid).Eids, five.N_Eids), std::vector<size_t>{eid}));
            }
        }
        AddContents(fiveEdges, tempEdges);
        size_t upperEid = GetDifference(five.N_Eids, fiveEdges).at(0);
        auto& upperE = mesh->E.at(upperEid);
        if ((upperE.Vids.at(0) != fiveId && upperE.Vids.at(0) == dest) || (upperE.Vids.at(1) != fiveId && upperE.Vids.at(1) == dest)) {
            // std::cout << "fiveId: " << fiveId << " dest: " << dest << " upperE.Vids.at(0): " << upperE.Vids.at(0) << " upperE.Vids.at(1): " << upperE.Vids.at(1) << std::endl;
            // std::cout << "threeId: " << threeId << " fiveId: " << fiveId << std::endl;
            // std::cout << "three: " << three.id << " " << three.N_Vids.size() << std::endl;
            // std::cout << "five: " << five.id << " " << five.N_Vids.size() << std::endl;
            // std::cout << "Move Up" << std::endl;
            // std::cout << "tempEdges: " << tempEdges.size() << std::endl;
            if (tempEdges.size() != 2) return -1;
            std::shared_ptr<SimplificationOperation> vs = std::make_shared<VertexSplit>(*mesh, *mu, *smoother, five.id, tempEdges);
            vs->PerformOperation();
            threeId = fiveId;
            fiveId = dest;
            delta += 2;
            CheckValences(fiveId, skipCheck, destFaces);
            // std::cout << "After performing operations" << std::endl;
            // std::cout << "three id: " << mesh->V.at(threeId).id << " three nfids: " << mesh->V.at(threeId).N_Fids.size() << " five id: " << mesh->V.at(fiveId).id << " five nfids: " << mesh->V.at(fiveId).N_Fids.size() << std::endl;
            // std::cout << "three id: " << mesh->V.at(threeId).id << " three nvids: " << mesh->V.at(threeId).N_Vids.size() << " five id: " << mesh->V.at(fiveId).id << " five nvids: " << mesh->V.at(fiveId).N_Vids.size() << std::endl;
            // std::cout << "three id: " << mesh->V.at(threeId).id << " three neids: " << mesh->V.at(threeId).N_Eids.size() << " five id: " << mesh->V.at(fiveId).id << " five neids: " << mesh->V.at(fiveId).N_Eids.size() << std::endl;
            // for (auto nvid: mesh->V.at(threeId).N_Vids) {
            //     auto& nv = mesh->V.at(nvid);
            //     std::cout << "nv id: " << nv.id << " nv nvids: " << nv.N_Vids.size() << std::endl;
            //     std::cout << "nv id: " << nv.id << " nv neids: " << nv.N_Eids.size() << std::endl;
            //     std::cout << "nv id: " << nv.id << " nv nfids: " << nv.N_Fids.size() << std::endl;
            // }
            // for (auto nvid: mesh->V.at(fiveId).N_Vids) {
            //     auto& nv = mesh->V.at(nvid);
            //     std::cout << "nv id: " << nv.id << " nv nvids: " << nv.N_Vids.size() << std::endl;
            //     std::cout << "nv id: " << nv.id << " nv neids: " << nv.N_Eids.size() << std::endl;
            //     std::cout << "nv id: " << nv.id << " nv nfids: " << nv.N_Fids.size() << std::endl;
            // }
            return 1;
        }
    }

    std::vector<size_t> downEdges = {GetDifference(three.N_Eids, GetIntersection(f.Eids, three.N_Eids)).at(0)};
    // std::cout << "downEdges: " << downEdges.size() << std::endl;
    newThreeId = mesh->E.at(downEdges.at(0)).Vids.at(0) == threeId ? mesh->E.at(downEdges.at(0)).Vids.at(1) : mesh->E.at(downEdges.at(0)).Vids.at(0);
    // AddContents(downEdges, GetIntersection(f.Eids, five.N_Eids));
    for (auto downEdgeId: downEdges) {
        auto& downEdge = mesh->E.at(downEdgeId);
        if ((downEdge.Vids.at(0) != threeId && downEdge.Vids.at(0) == dest) || (downEdge.Vids.at(1) != threeId && downEdge.Vids.at(1) == dest)) {
            // std::cout << "Move Down" << std::endl;
            int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), three.id));
            int t = (idx+1)%f.Vids.size();
            int s = (idx+3)%f.Vids.size();
            // std::cout << "Inside Diagonal pair: " << f.Vids.at(t) << " " << mesh->V.at(f.Vids.at(t)).type << " " << f.Vids.at(s) << " " << mesh->V.at(f.Vids.at(s)).type << std::endl;
            // std::cout << "Inside Diagonal pair: " << f.Vids.at(t) << " " << mesh->V.at(f.Vids.at(t)).isBoundary << " " << f.Vids.at(s) << " " << mesh->V.at(f.Vids.at(s)).isBoundary << std::endl;
            // if ((mesh->V.at(f.Vids.at(t)).type != FEATURE && mesh->V.at(f.Vids.at(s)).type == FEATURE) || (!mesh->V.at(f.Vids.at(t)).isBoundary && mesh->V.at(f.Vids.at(s)).isBoundary) || mu->IsSharpFeature(s)) {
            //     int tmp = s;
            //     s = t;
            //     t = tmp;
            // }
            /*if (skipCheck && mesh->V.at(f.Vids.at(t)).N_Vids.size() == 3 && !mesh->V.at(f.Vids.at(t)).isBoundary && mesh->V.at(f.Vids.at(s)).N_Vids.size() == 4) {
                Move(f.Vids.at(s), delta, skipCheck);
                auto& tThree = mesh->V.at(threeId);
                auto& tFive = mesh->V.at(fiveId);
                auto& tF = mesh->F.at(GetIntersection(tThree.N_Fids, tFive.N_Fids).at(0));
                size_t tEid = GetDifference(tThree.N_Eids, GetIntersection(tThree.N_Eids, tF.Eids)).at(0);
                size_t tDest = mesh->E.at(tEid).Vids.at(0) == threeId ? mesh->E.at(tEid).Vids.at(1) : mesh->E.at(tEid).Vids.at(0);
                int tIdx = std::distance(tF.Vids.begin(), std::find(tF.Vids.begin(), tF.Vids.end(), tThree.id));
                t = (tIdx+1)%tF.Vids.size();
                s = (tIdx+3)%tF.Vids.size();
                newThreeId = tDest;
                newFiveId = tF.Vids.at(t);
                std::shared_ptr<SimplificationOperation> dc = std::make_shared<DiagonalCollapse>(*mesh, *mu, *smoother, tF.id, t, s);
                dc->PerformOperation();
                threeId = newThreeId;
                fiveId = newFiveId;
                delta += -2;
                CheckValences(fiveId, skipCheck, destFaces);
                // Move(tDest, delta, skipCheck);
                return 1;
            } else if (skipCheck && mesh->V.at(f.Vids.at(s)).N_Vids.size() == 3 && !mesh->V.at(f.Vids.at(s)).isBoundary && mesh->V.at(f.Vids.at(t)).N_Vids.size() == 4) {
                Move(f.Vids.at(t), delta, skipCheck);
                auto& tThree = mesh->V.at(threeId);
                auto& tFive = mesh->V.at(fiveId);
                auto& tF = mesh->F.at(GetIntersection(tThree.N_Fids, tFive.N_Fids).at(0));
                size_t tEid = GetDifference(tThree.N_Eids, GetIntersection(tThree.N_Eids, tF.Eids)).at(0);
                size_t tDest = mesh->E.at(tEid).Vids.at(0) == threeId ? mesh->E.at(tEid).Vids.at(1) : mesh->E.at(tEid).Vids.at(0);
                int tIdx = std::distance(tF.Vids.begin(), std::find(tF.Vids.begin(), tF.Vids.end(), tThree.id));
                t = (tIdx+1)%tF.Vids.size();
                s = (tIdx+3)%tF.Vids.size();
                newThreeId = tDest;
                newFiveId = tF.Vids.at(t);
                std::shared_ptr<SimplificationOperation> dc = std::make_shared<DiagonalCollapse>(*mesh, *mu, *smoother, tF.id, t, s);
                dc->PerformOperation();
                threeId = newThreeId;
                fiveId = newFiveId;
                delta += -2;
                CheckValences(fiveId, skipCheck, destFaces);
                // Move(tDest, delta, skipCheck);
                return 1;
            }*/
            newFiveId = f.Vids.at(t);
            std::shared_ptr<SimplificationOperation> dc = std::make_shared<DiagonalCollapse>(*mesh, *mu, *smoother, f.id, t, s);
            dc->PerformOperation();
            threeId = newThreeId;
            fiveId = newFiveId;
            delta += -2;
            CheckValences(fiveId, skipCheck, destFaces);
            
            // std::cout << "three id: " << mesh->V.at(threeId).id << " three nfids: " << mesh->V.at(threeId).N_Fids.size() << " five id: " << mesh->V.at(fiveId).id << " five nfids: " << mesh->V.at(fiveId).N_Fids.size() << std::endl;
            // std::cout << "three id: " << mesh->V.at(threeId).id << " three nvids: " << mesh->V.at(threeId).N_Vids.size() << " five id: " << mesh->V.at(fiveId).id << " five nvids: " << mesh->V.at(fiveId).N_Vids.size() << std::endl;
            // std::cout << "three id: " << mesh->V.at(threeId).id << " three neids: " << mesh->V.at(threeId).N_Eids.size() << " five id: " << mesh->V.at(fiveId).id << " five neids: " << mesh->V.at(fiveId).N_Eids.size() << std::endl;
            // for (auto nvid: mesh->V.at(threeId).N_Vids) {
            //     auto& nv = mesh->V.at(nvid);
            //     std::cout << "nv id: " << nv.id << " nv nvids: " << nv.N_Vids.size() << std::endl;
            //     std::cout << "nv id: " << nv.id << " nv neids: " << nv.N_Eids.size() << std::endl;
            //     std::cout << "nv id: " << nv.id << " nv nfids: " << nv.N_Fids.size() << std::endl;
            // }
            // for (auto nvid: mesh->V.at(fiveId).N_Vids) {
            //     auto& nv = mesh->V.at(nvid);
            //     std::cout << "nv id: " << nv.id << " nv nvids: " << nv.N_Vids.size() << std::endl;
            //     std::cout << "nv id: " << nv.id << " nv neids: " << nv.N_Eids.size() << std::endl;
            //     std::cout << "nv id: " << nv.id << " nv nfids: " << nv.N_Fids.size() << std::endl;
            // }
            return 1;
        }
    }
    

    std::vector<size_t> fEdges = GetIntersection(f.Eids, five.N_Eids);
    bool opPerformed = false;
    for (auto eid: fEdges) {
        auto& e = mesh->E.at(eid);
        if (e.N_Fids.size() != 2) continue;
        for (auto fid: e.N_Fids) {
            auto& ef = mesh->F.at(fid);
            int idx = std::distance(ef.Vids.begin(), std::find(ef.Vids.begin(), ef.Vids.end(), five.id));
            if ((fid != f.id && ef.Vids.at((idx+1)%ef.Vids.size()) == dest && (dest != e.Vids.at(0) && dest != e.Vids.at(1))) || (fid == f.id && ef.Vids.at((idx+1)%ef.Vids.size()) == dest && (dest == e.Vids.at(0) || dest == e.Vids.at(1)))) {
                // std::cout << "Moving Left" << std::endl;
                int nFid = e.N_Fids.at(0) == f.id ? e.N_Fids.at(1) : e.N_Fids.at(0);
                auto& nEf = mesh->F.at(nFid);
                int fiveIdx = std::distance(nEf.Vids.begin(), std::find(nEf.Vids.begin(), nEf.Vids.end(), five.id));
                /*if (skipCheck && dest != e.Vids.at(0) && dest != e.Vids.at(1) && 
                    mesh->V.at(nEf.Vids.at((fiveIdx+3)%nEf.Vids.size())).N_Vids.size() == 3 && mesh->V.at(nEf.Vids.at((fiveIdx+1)%nEf.Vids.size())).N_Vids.size() == 4 &&
                    !mesh->V.at(nEf.Vids.at((fiveIdx+3)%nEf.Vids.size())).isBoundary) {
                    // std::cout << "Encountered 3 valence vertex while moving Left" << std::endl;
                    std::vector<size_t> tFs;
                    for (auto tEid: fEdges) {
                        AddContents(tFs, GetIntersection(mesh->F.at(mesh->E.at(tEid).N_Fids.at(0)).Eids, five.N_Eids));
                        AddContents(tFs, GetIntersection(mesh->F.at(mesh->E.at(tEid).N_Fids.at(1)).Eids, five.N_Eids));
                    }
                    auto tDiffs = GetDifference(five.N_Eids, tFs);
                    if (tDiffs.empty() || five.isBoundary) return 1;
                    size_t tDiffEid = tDiffs.at(0);
                    // std::cout << "After getting contents" << std::endl;
                    size_t tDest = mesh->E.at(tDiffEid).Vids.at(0) == five.id ? mesh->E.at(tDiffEid).Vids.at(1) : mesh->E.at(tDiffEid).Vids.at(0);
                    Move(tDest, delta, skipCheck);
                    auto& tF = mesh->F.at(GetIntersection(mesh->V.at(threeId).N_Fids, mesh->V.at(fiveId).N_Fids).at(0));
                    for (auto tFvid: tF.Vids) {
                        auto& tFv = mesh->V.at(tFvid);
                        if (std::find(tFv.N_Vids.begin(), tFv.N_Vids.end(), dest) != tFv.N_Vids.end()) {
                            Move(tFvid, delta, skipCheck);
                            return 1;
                        }
                    }
                }
                
                if (skipCheck && (dest == e.Vids.at(0) || dest == e.Vids.at(1)) && 
                    ((mesh->V.at(nEf.Vids.at((fiveIdx+1)%nEf.Vids.size())).N_Vids.size() == 3 && mesh->V.at(nEf.Vids.at((fiveIdx+3)%nEf.Vids.size())).N_Vids.size() == 4) ||
                    (mesh->V.at(nEf.Vids.at((fiveIdx+3)%nEf.Vids.size())).N_Vids.size() == 5 && mesh->V.at(nEf.Vids.at((fiveIdx+1)%nEf.Vids.size())).N_Vids.size() == 4)) &&
                    !mesh->V.at(nEf.Vids.at((fiveIdx+1)%nEf.Vids.size())).isBoundary && !mesh->V.at(nEf.Vids.at((fiveIdx+3)%nEf.Vids.size())).isBoundary) {
                    // std::cout << "Encountered 3 valence vertex while moving Left and dest is close to three or 5" << std::endl;
                    size_t tEid = GetDifference(three.N_Eids, GetIntersection(three.N_Eids, f.Eids)).at(0);
                    size_t tDest = mesh->E.at(tEid).Vids.at(0) == three.id ? mesh->E.at(tEid).Vids.at(1) : mesh->E.at(tEid).Vids.at(0);
                    // Move(tDest, delta, skipCheck);
                    int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), three.id));
                    int t = (idx+1)%f.Vids.size();
                    int s = (idx+3)%f.Vids.size();
                    newFiveId = f.Vids.at(t);
                    newThreeId = tDest;
                    std::shared_ptr<SimplificationOperation> dc = std::make_shared<DiagonalCollapse>(*mesh, *mu, *smoother, f.id, t, s);
                    dc->PerformOperation();
                    threeId = newThreeId;
                    fiveId = newFiveId;
                    
                    return 1;
                }*/
                fiveId = nEf.Vids.at((fiveIdx+1)%nEf.Vids.size());
                threeId = nEf.Vids.at((fiveIdx+3)%nEf.Vids.size());
                
                std::shared_ptr<SimplificationOperation> s = std::make_shared<EdgeRotation>(*mesh, *mu, *smoother, e.id, false);
                s->PerformOperation();
                opPerformed = true;
                CheckValences(fiveId, skipCheck, destFaces);
                // std::cout << "three id: " << mesh->V.at(threeId).id << " three nfids: " << mesh->V.at(threeId).N_Fids.size() << " five id: " << mesh->V.at(fiveId).id << " five nfids: " << mesh->V.at(fiveId).N_Fids.size() << std::endl;
                // std::cout << "three id: " << mesh->V.at(threeId).id << " three nvids: " << mesh->V.at(threeId).N_Vids.size() << " five id: " << mesh->V.at(fiveId).id << " five nvids: " << mesh->V.at(fiveId).N_Vids.size() << std::endl;
                // std::cout << "three id: " << mesh->V.at(threeId).id << " three neids: " << mesh->V.at(threeId).N_Eids.size() << " five id: " << mesh->V.at(fiveId).id << " five neids: " << mesh->V.at(fiveId).N_Eids.size() << std::endl;
                // for (auto nvid: mesh->V.at(threeId).N_Vids) {
                //     auto& nv = mesh->V.at(nvid);
                //     std::cout << "nv id: " << nv.id << " nv nvids: " << nv.N_Vids.size() << std::endl;
                //     std::cout << "nv id: " << nv.id << " nv neids: " << nv.N_Eids.size() << std::endl;
                //     std::cout << "nv id: " << nv.id << " nv nfids: " << nv.N_Fids.size() << std::endl;
                // }
                // for (auto nvid: mesh->V.at(fiveId).N_Vids) {
                //     auto& nv = mesh->V.at(nvid);
                //     std::cout << "nv id: " << nv.id << " nv nvids: " << nv.N_Vids.size() << std::endl;
                //     std::cout << "nv id: " << nv.id << " nv neids: " << nv.N_Eids.size() << std::endl;
                //     std::cout << "nv id: " << nv.id << " nv nfids: " << nv.N_Fids.size() << std::endl;
                // }
                return 1;
            } else if ((fid != f.id && ef.Vids.at((idx+3)%ef.Vids.size()) == dest && (dest != e.Vids.at(0) && dest != e.Vids.at(1))) || (fid == f.id && ef.Vids.at((idx+3)%ef.Vids.size()) == dest && (dest == e.Vids.at(0) || dest == e.Vids.at(1)))) {
                // std::cout << "Moving Right" << std::endl;
                int nFid = e.N_Fids.at(0) == f.id ? e.N_Fids.at(1) : e.N_Fids.at(0);
                auto& nEf = mesh->F.at(nFid);
                int fiveIdx = std::distance(nEf.Vids.begin(), std::find(nEf.Vids.begin(), nEf.Vids.end(), five.id));
                /*if (skipCheck && dest != e.Vids.at(0) && dest != e.Vids.at(1) && 
                    mesh->V.at(nEf.Vids.at((fiveIdx+1)%nEf.Vids.size())).N_Vids.size() == 3 && mesh->V.at(nEf.Vids.at((fiveIdx+3)%nEf.Vids.size())).N_Vids.size() == 4 && 
                    !mesh->V.at(nEf.Vids.at((fiveIdx+1)%nEf.Vids.size())).isBoundary) {
                    // std::cout << "Encountered 3 valence vertex while moving Right" << std::endl;
                    // std::cout << "three: ";for (auto threenvid: mesh->V.at(threeId).N_Vids) std::cout << threenvid << " ";std::cout << std::endl;
                    // std::cout << "five: ";for (auto fivenvid: mesh->V.at(fiveId).N_Vids) std::cout << fivenvid << " ";std::cout << std::endl;
                    // std::cout << "face edges: ";for (auto feid: f.Eids) std::cout << feid << " ";std::cout << std::endl;
                    // std::cout << "three edges: ";for (auto threeneid: mesh->V.at(threeId).N_Eids) std::cout << threeneid << " ";std::cout << std::endl;
                    // std::cout << "five edges: ";for (auto fiveneid: mesh->V.at(fiveId).N_Eids) std::cout << fiveneid << " ";std::cout << std::endl;
                    // std::cout << "fEdges: ";for (auto fEid: fEdges) std::cout << fEid << " ";std::cout << std::endl;

                    std::vector<size_t> tFs;
                    for (auto tEid: fEdges) {
                        AddContents(tFs, GetIntersection(mesh->F.at(mesh->E.at(tEid).N_Fids.at(0)).Eids, five.N_Eids));
                        AddContents(tFs, GetIntersection(mesh->F.at(mesh->E.at(tEid).N_Fids.at(1)).Eids, five.N_Eids));
                    }
                    auto tDiffs = GetDifference(five.N_Eids, tFs);
                    if (tDiffs.empty() || five.isBoundary) return 1;
                    size_t tDiffEid = tDiffs.at(0);
                    
                    // std::cout << "tDiffEid: " << tDiffEid << std::endl;
                    size_t tDest = mesh->E.at(tDiffEid).Vids.at(0) == five.id ? mesh->E.at(tDiffEid).Vids.at(1) : mesh->E.at(tDiffEid).Vids.at(0);
                    // std::cout << "dest: " << dest << std::endl;
                    // std::cout << "tDest: " << tDest << std::endl;
                    Move(tDest, delta, skipCheck);
                    auto& tF = mesh->F.at(GetIntersection(mesh->V.at(threeId).N_Fids, mesh->V.at(fiveId).N_Fids).at(0));
                    for (auto tFvid: tF.Vids) {
                        auto& tFv = mesh->V.at(tFvid);
                        if (std::find(tFv.N_Vids.begin(), tFv.N_Vids.end(), dest) != tFv.N_Vids.end()) {
                            Move(tFvid, delta, skipCheck);
                            return 1;
                        }
                    }
                    // std::cout << "three: ";for (auto threenvid: mesh->V.at(threeId).N_Vids) std::cout << threenvid << " ";std::cout << std::endl;
                    // std::cout << "five: ";for (auto fivenvid: mesh->V.at(fiveId).N_Vids) std::cout << fivenvid << " ";std::cout << std::endl;
                    // std::cout << "Moved Up" << std::endl;
                    // auto& tThree = mesh->V.at(threeId);
                    // std::cout << "tThree: " << tThree.id << std::endl;
                    // std::cout << "tThree nvids: " << tThree.N_Vids.size() << std::endl;
                    // std::cout << "tThree neids: " << tThree.N_Eids.size() << std::endl;
                    // std::cout << "tThree nfids: " << tThree.N_Fids.size() << std::endl;
                    // for (auto threeFid: tThree.N_Fids) {
                    //     auto& tThreeF = mesh->F.at(threeFid);
                    //     std::cout << "f: " << tThreeF.id; for (auto vid: tThreeF.Vids) std::cout << " " << vid; std::cout << std::endl;
                    //     int tThreeIdx = std::distance(tThreeF.Vids.begin(), std::find(tThreeF.Vids.begin(), tThreeF.Vids.end(), threeId));
                    //     std::cout << "tThreeIdx: " << tThreeIdx << std::endl;
                    //     std::cout << "three: " << tThreeF.Vids.at((tThreeIdx)%tThreeF.Vids.size()) << std::endl;
                    //     std::cout << "dest: " << dest << " ";
                    //     std::cout << "face vid: " << tThreeF.Vids.at((tThreeIdx+2)%tThreeF.Vids.size()) << " " << (tThreeF.Vids.at((tThreeIdx+2)%tThreeF.Vids.size()) == dest) << std::endl;
                    //     if (tThreeF.Vids.at((tThreeIdx+2)%tThreeF.Vids.size()) == dest) {
                    //         Move(tThreeF.Vids.at((tThreeIdx+3)%tThreeF.Vids.size()), delta, skipCheck);
                    //         std::cout << "Moved Right" << std::endl;
                    //         return 1;
                    //     }
                    // }
                }
                // std::cout << "e: " << e.Vids.at(0) << " " << e.Vids.at(1) << std::endl;
                // std::cout << "dest: " << dest << std::endl;
                // std::cout << "fiveIdx+1: " << mesh->V.at(nEf.Vids.at((fiveIdx+1)%nEf.Vids.size())).N_Vids.size() << std::endl;
                // std::cout << "fiveIdx+3: " << mesh->V.at(nEf.Vids.at((fiveIdx+3)%nEf.Vids.size())).N_Vids.size() << std::endl;
                if (skipCheck && (dest == e.Vids.at(0) || dest == e.Vids.at(1)) && 
                    ((mesh->V.at(nEf.Vids.at((fiveIdx+3)%nEf.Vids.size())).N_Vids.size() == 3 && mesh->V.at(nEf.Vids.at((fiveIdx+1)%nEf.Vids.size())).N_Vids.size() == 4) || 
                    (mesh->V.at(nEf.Vids.at((fiveIdx+1)%nEf.Vids.size())).N_Vids.size() == 5 && mesh->V.at(nEf.Vids.at((fiveIdx+3)%nEf.Vids.size())).N_Vids.size() == 4)) &&
                    !mesh->V.at(nEf.Vids.at((fiveIdx+3)%nEf.Vids.size())).isBoundary && !mesh->V.at(nEf.Vids.at((fiveIdx+1)%nEf.Vids.size())).isBoundary) {
                    // std::cout << "Encountered 3 valence vertex while moving Right and dest is close to three or 5" << std::endl;
                    size_t tEid = GetDifference(three.N_Eids, GetIntersection(three.N_Eids, f.Eids)).at(0);
                    size_t tDest = mesh->E.at(tEid).Vids.at(0) == three.id ? mesh->E.at(tEid).Vids.at(1) : mesh->E.at(tEid).Vids.at(0);
                    // Move(tDest, delta, skipCheck);
                    int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), three.id));
                    int t = (idx+1)%f.Vids.size();
                    int s = (idx+3)%f.Vids.size();
                    newFiveId = f.Vids.at(t);
                    newThreeId = tDest;
                    std::shared_ptr<SimplificationOperation> dc = std::make_shared<DiagonalCollapse>(*mesh, *mu, *smoother, f.id, t, s);
                    dc->PerformOperation();
                    threeId = newThreeId;
                    fiveId = newFiveId;
                    return 1;
                }*/
                fiveId = nEf.Vids.at((fiveIdx+3)%nEf.Vids.size());
                threeId = nEf.Vids.at((fiveIdx+1)%nEf.Vids.size());
                
                std::shared_ptr<SimplificationOperation> s = std::make_shared<EdgeRotation>(*mesh, *mu, *smoother, e.id, true);
                s->PerformOperation();
                opPerformed = true;
                CheckValences(fiveId, skipCheck, destFaces);
                // std::cout << "three id: " << mesh->V.at(threeId).id << " three nfids: " << mesh->V.at(threeId).N_Fids.size() << " five id: " << mesh->V.at(fiveId).id << " five nfids: " << mesh->V.at(fiveId).N_Fids.size() << std::endl;
                // std::cout << "three id: " << mesh->V.at(threeId).id << " three nvids: " << mesh->V.at(threeId).N_Vids.size() << " five id: " << mesh->V.at(fiveId).id << " five nvids: " << mesh->V.at(fiveId).N_Vids.size() << std::endl;
                // std::cout << "three id: " << mesh->V.at(threeId).id << " three neids: " << mesh->V.at(threeId).N_Eids.size() << " five id: " << mesh->V.at(fiveId).id << " five neids: " << mesh->V.at(fiveId).N_Eids.size() << std::endl;
                // for (auto nvid: mesh->V.at(threeId).N_Vids) {
                //     auto& nv = mesh->V.at(nvid);
                //     std::cout << "nv id: " << nv.id << " nv nvids: " << nv.N_Vids.size() << std::endl;
                //     std::cout << "nv id: " << nv.id << " nv neids: " << nv.N_Eids.size() << std::endl;
                //     std::cout << "nv id: " << nv.id << " nv nfids: " << nv.N_Fids.size() << std::endl;
                // }
                // for (auto nvid: mesh->V.at(fiveId).N_Vids) {
                //     auto& nv = mesh->V.at(nvid);
                //     std::cout << "nv id: " << nv.id << " nv nvids: " << nv.N_Vids.size() << std::endl;
                //     std::cout << "nv id: " << nv.id << " nv neids: " << nv.N_Eids.size() << std::endl;
                //     std::cout << "nv id: " << nv.id << " nv nfids: " << nv.N_Fids.size() << std::endl;
                // }
                return 1;
            }
        }
        // if (opPerformed) break;
    }

    return -1;
}

void DiagonalThreeFivePair::SetResolvedSingularity(size_t dest, int ref) {
    
}

int DiagonalThreeFivePair::GetResolvedSingularity() {
    return resolvedSingularity;
}
/*
void DiagonalThreeFivePair::MoveUp() {

}

void DiagonalThreeFivePair::MoveDown() {

}

void DiagonalThreeFivePair::MoveLeft() {

}

void DiagonalThreeFivePair::MoveRight() {
    
}*/

void DiagonalThreeFivePair::CheckValences(int vid, bool skipCheck, std::vector<size_t> faces) {
    if (skipCheck) return;
    // std::cout << "skipCheck: " << skipCheck << " CheckValences" << std::endl;
    if (!mesh->V.at(vid).isBoundary && mesh->V.at(vid).N_Vids.size() == 6) {
        // std::cout << "Split 6-singularity" << std::endl;
        SplitSixSingularity(vid);
    } else if (mesh->V.at(vid).isBoundary && mesh->V.at(vid).N_Vids.size() == 4) {
        // std::cout << "v: " << mesh->V.at(vid).N_Vids.size() << " is boundary" << std::endl; 
        std::shared_ptr<SimplificationOperation> s = std::make_shared<VertexSplit>(*mesh, *mu, *smoother, vid);
        s->PerformOperation();
    } else {
        // std::cout << "v: " << mesh->V.at(vid).N_Vids.size() << std::endl; 
        // for (auto fid: mesh->V.at(vid).N_Fids) {
        for (auto fid: faces) {
            auto& f = mesh->F.at(fid);
            if (f.Vids.empty()) continue;
            int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), vid));
            auto& t = mesh->V.at(f.Vids.at((idx+1)%f.Vids.size()));
            auto& d = mesh->V.at(f.Vids.at((idx+2)%f.Vids.size()));
            auto& s = mesh->V.at(f.Vids.at((idx+3)%f.Vids.size()));
            if (((t.isBoundary && t.idealValence != 2) || (t.type == FEATURE && t.idealValence != 4)) || 
                ((d.isBoundary && d.idealValence != 2) || (d.type == FEATURE && d.idealValence != 4)) ||
                ((s.isBoundary && s.idealValence != 2) || (s.type == FEATURE && s.idealValence != 4))) continue;
            // if ((!t.isBoundary && !mu->IsSharpFeature(t.id) && t.N_Vids.size() == 3) && (!s.isBoundary && !mu->IsSharpFeature(s.id) && s.N_Vids.size() == 3)) {
            if (t.N_Vids.size() == 3 && s.N_Vids.size() == 3) {
                // std::cout << "Collapsing diagonal, threeId: " << mesh->V.at(threeId).N_Vids.size() << std::endl;
                std::shared_ptr<SimplificationOperation> dc2 = std::make_shared<DiagonalCollapse>(*mesh, *mu, *smoother, f.id, (idx+1)%f.Vids.size(), (idx+3)%f.Vids.size());
                dc2->PerformOperation();
                // std::cout << "After diagonal collapse" << std::endl;
                return;
            } 
        }
    }
}

void DiagonalThreeFivePair::SplitSixSingularity(int vid) {
    auto& v = mesh->V.at(vid);
    for (auto fid: v.N_Fids) {
        auto& f = mesh->F.at(fid);
        int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
        auto& t = mesh->V.at(f.Vids.at((idx+1)%f.Vids.size()));
        auto& s = mesh->V.at(f.Vids.at((idx+3)%f.Vids.size()));
        if (mesh->V.at(f.Vids.at((idx+2)%f.Vids.size())).N_Vids.size() == 3) {
            std::vector<size_t> edgesToSplit;
            for (auto eid: mu->GetIntersection(f.Eids, v.N_Eids)) {
                auto& e = mesh->E.at(eid);
                mu->AddContents(edgesToSplit, mu->GetDifference(mu->GetIntersection(mesh->F.at(mu->GetDifference(e.N_Fids, std::vector<size_t>{fid})[0]).Eids, v.N_Eids), std::vector<size_t>{eid}));
            }
            std::shared_ptr<SimplificationOperation> s = std::make_shared<VertexSplit>(*mesh, *mu, *smoother, v.id, edgesToSplit);
            s->PerformOperation();
            return;
        } else if ((!mu->IsSharpFeature(t.id) && t.N_Vids.size() == 3) && (!mu->IsSharpFeature(s.id) && s.N_Vids.size() == 3)) {
            std::shared_ptr<SimplificationOperation> dc2 = std::make_shared<DiagonalCollapse>(*mesh, *mu, *smoother, f.id, (idx+1)%f.Vids.size(), (idx+3)%f.Vids.size());
            dc2->PerformOperation();
            return;
        }
    }
}

std::vector<size_t> DiagonalThreeFivePair::GetPairIds() {
    CheckValidity();
    return std::vector<size_t>{threeId, fiveId};
}

