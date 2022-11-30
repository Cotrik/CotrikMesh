#include "DiagonalThreeFivePair.h"

DiagonalThreeFivePair::DiagonalThreeFivePair() {}

DiagonalThreeFivePair::DiagonalThreeFivePair(Mesh& mesh_, MeshUtil& mu_, Smoother& smoother_, size_t threeId_, size_t fiveId_) : mesh(&mesh_), mu(&mu_), smoother(&smoother_) {
    threeId = threeId_;
    fiveId = fiveId_;

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
}

int DiagonalThreeFivePair::Move(size_t dest, bool skipCheck) {
    CheckValidity();
    if (!IsOperationValid()) return -1;
    auto& three = mesh->V.at(threeId);
    auto& five = mesh->V.at(fiveId);
    // std::cout << "dest: " << dest << std::endl;
    // std::cout << "three id: " << three.id << " three nfids: " << three.N_Fids.size() << " five id: " << five.id << " five nfids: " << five.N_Fids.size() << std::endl;
    // std::cout << "three id: " << three.id << " three nvids: " << three.N_Vids.size() << " five id: " << five.id << " five nvids: " << five.N_Vids.size() << std::endl;
    // std::cout << "three id: " << three.id << " three neids: " << three.N_Eids.size() << " five id: " << five.id << " five neids: " << five.N_Eids.size() << std::endl;
    faceId = GetIntersection(three.N_Fids, five.N_Fids).at(0);
    int newThreeId = -1;
    int newFiveId = -1;

    auto& f = mesh->F.at(faceId);
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
    if (upperE.Vids.at(0) == dest || upperE.Vids.at(1) == dest) {
        // std::cout << "Move Up" << std::endl;
        std::shared_ptr<SimplificationOperation> vs = std::make_shared<VertexSplit>(*mesh, *mu, *smoother, five.id, tempEdges);
        vs->PerformOperation();
        threeId = five.id;
        fiveId = dest;
        // std::cout << "three id: " << mesh->V.at(threeId).id << " three nfids: " << mesh->V.at(threeId).N_Fids.size() << " five id: " << mesh->V.at(fiveId).id << " five nfids: " << mesh->V.at(fiveId).N_Fids.size() << std::endl;
        // std::cout << "three id: " << mesh->V.at(threeId).id << " three nvids: " << mesh->V.at(threeId).N_Vids.size() << " five id: " << mesh->V.at(fiveId).id << " five nvids: " << mesh->V.at(fiveId).N_Vids.size() << std::endl;
        // std::cout << "three id: " << mesh->V.at(threeId).id << " three neids: " << mesh->V.at(threeId).N_Eids.size() << " five id: " << mesh->V.at(fiveId).id << " five neids: " << mesh->V.at(fiveId).N_Eids.size() << std::endl;
        return 1;
    }

    std::vector<size_t> downEdges = {GetDifference(three.N_Eids, GetIntersection(f.Eids, three.N_Eids)).at(0)};
    newThreeId = mesh->E.at(downEdges.at(0)).Vids.at(0) == threeId ? mesh->E.at(downEdges.at(0)).Vids.at(1) : mesh->E.at(downEdges.at(0)).Vids.at(0);
    // AddContents(downEdges, GetIntersection(f.Eids, five.N_Eids));
    for (auto downEdgeId: downEdges) {
        auto& downEdge = mesh->E.at(downEdgeId);
        if (downEdge.Vids.at(0) == dest || downEdge.Vids.at(1) == dest) {
            // std::cout << "Move Down" << std::endl;
            int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), three.id));
            newFiveId = f.Vids.at((idx+1)%f.Vids.size());
            std::shared_ptr<SimplificationOperation> dc = std::make_shared<DiagonalCollapse>(*mesh, *mu, *smoother, f.id, (idx+1)%f.Vids.size(), (idx+3)%f.Vids.size());
            dc->PerformOperation();
            threeId = newThreeId;
            fiveId = newFiveId;
            // std::cout << "three id: " << mesh->V.at(threeId).id << " three nfids: " << mesh->V.at(threeId).N_Fids.size() << " five id: " << mesh->V.at(fiveId).id << " five nfids: " << mesh->V.at(fiveId).N_Fids.size() << std::endl;
            // std::cout << "three id: " << mesh->V.at(threeId).id << " three nvids: " << mesh->V.at(threeId).N_Vids.size() << " five id: " << mesh->V.at(fiveId).id << " five nvids: " << mesh->V.at(fiveId).N_Vids.size() << std::endl;
            // std::cout << "three id: " << mesh->V.at(threeId).id << " three neids: " << mesh->V.at(threeId).N_Eids.size() << " five id: " << mesh->V.at(fiveId).id << " five neids: " << mesh->V.at(fiveId).N_Eids.size() << std::endl;
            return 1;
        }
    }
    

    std::vector<size_t> fEdges = GetIntersection(f.Eids, five.N_Eids);
    bool opPerformed = false;
    for (auto eid: fEdges) {
        auto& e = mesh->E.at(eid);
        for (auto fid: e.N_Fids) {
            auto& ef = mesh->F.at(fid);
            int idx = std::distance(ef.Vids.begin(), std::find(ef.Vids.begin(), ef.Vids.end(), five.id));
            if ((fid != f.id && ef.Vids.at((idx+1)%ef.Vids.size()) == dest && (dest != e.Vids.at(0) && dest != e.Vids.at(1))) || (fid == f.id && ef.Vids.at((idx+1)%ef.Vids.size()) == dest && (dest == e.Vids.at(0) || dest == e.Vids.at(1)))) {
                // std::cout << "Moving Left" << std::endl;
                int nFid = e.N_Fids.at(0) == f.id ? e.N_Fids.at(1) : e.N_Fids.at(0);
                auto& nEf = mesh->F.at(nFid);
                int fiveIdx = std::distance(nEf.Vids.begin(), std::find(nEf.Vids.begin(), nEf.Vids.end(), five.id));
                fiveId = nEf.Vids.at((fiveIdx+1)%nEf.Vids.size());
                threeId = nEf.Vids.at((fiveIdx+3)%nEf.Vids.size());
                
                std::shared_ptr<SimplificationOperation> s = std::make_shared<EdgeRotation>(*mesh, *mu, *smoother, e.id, false);
                s->PerformOperation();
                opPerformed = true;
                // std::cout << "three id: " << mesh->V.at(threeId).id << " three nfids: " << mesh->V.at(threeId).N_Fids.size() << " five id: " << mesh->V.at(fiveId).id << " five nfids: " << mesh->V.at(fiveId).N_Fids.size() << std::endl;
                // std::cout << "three id: " << mesh->V.at(threeId).id << " three nvids: " << mesh->V.at(threeId).N_Vids.size() << " five id: " << mesh->V.at(fiveId).id << " five nvids: " << mesh->V.at(fiveId).N_Vids.size() << std::endl;
                // std::cout << "three id: " << mesh->V.at(threeId).id << " three neids: " << mesh->V.at(threeId).N_Eids.size() << " five id: " << mesh->V.at(fiveId).id << " five neids: " << mesh->V.at(fiveId).N_Eids.size() << std::endl;
    
                return 1;
            } else if ((fid != f.id && ef.Vids.at((idx+3)%ef.Vids.size()) == dest && (dest != e.Vids.at(0) && dest != e.Vids.at(1))) || (fid == f.id && ef.Vids.at((idx+3)%ef.Vids.size()) == dest && (dest == e.Vids.at(0) || dest == e.Vids.at(1)))) {
                // std::cout << "Moving Right" << std::endl;
                int nFid = e.N_Fids.at(0) == f.id ? e.N_Fids.at(1) : e.N_Fids.at(0);
                auto& nEf = mesh->F.at(nFid);
                int fiveIdx = std::distance(nEf.Vids.begin(), std::find(nEf.Vids.begin(), nEf.Vids.end(), five.id));
                fiveId = nEf.Vids.at((fiveIdx+3)%nEf.Vids.size());
                threeId = nEf.Vids.at((fiveIdx+1)%nEf.Vids.size());
                
                std::shared_ptr<SimplificationOperation> s = std::make_shared<EdgeRotation>(*mesh, *mu, *smoother, e.id, true);
                s->PerformOperation();
                opPerformed = true;
                // std::cout << "three id: " << mesh->V.at(threeId).id << " three nfids: " << mesh->V.at(threeId).N_Fids.size() << " five id: " << mesh->V.at(fiveId).id << " five nfids: " << mesh->V.at(fiveId).N_Fids.size() << std::endl;
                // std::cout << "three id: " << mesh->V.at(threeId).id << " three nvids: " << mesh->V.at(threeId).N_Vids.size() << " five id: " << mesh->V.at(fiveId).id << " five nvids: " << mesh->V.at(fiveId).N_Vids.size() << std::endl;
                // std::cout << "three id: " << mesh->V.at(threeId).id << " three neids: " << mesh->V.at(threeId).N_Eids.size() << " five id: " << mesh->V.at(fiveId).id << " five neids: " << mesh->V.at(fiveId).N_Eids.size() << std::endl;
    
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

void DiagonalThreeFivePair::MoveUp() {

}

void DiagonalThreeFivePair::MoveDown() {

}

void DiagonalThreeFivePair::MoveLeft() {

}

void DiagonalThreeFivePair::MoveRight() {
    
}


void DiagonalThreeFivePair::SplitSixSingularity() {
    
}

std::vector<size_t> DiagonalThreeFivePair::GetPairIds() {
    CheckValidity();
    return std::vector<size_t>{threeId, fiveId};
}

