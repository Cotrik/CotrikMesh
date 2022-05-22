#include "VertexRotation.h"
#include <queue>

VertexRotation::VertexRotation() {}
VertexRotation::VertexRotation(Mesh& mesh_, MeshUtil& mu_, Smoother& smoother_, size_t vid_) : SimplificationOperation(mesh_, mu_, smoother_) {
    vid = vid_;
}
VertexRotation::~VertexRotation() {}

bool VertexRotation::IsOperationValid() {
    CheckValidity();
    if (mesh.V.at(vid).N_Vids.empty() || mesh.V.at(vid).N_Eids.empty() || mesh.V.at(vid).N_Fids.empty()) return false;
    // if (mesh.V.at(vid).N_Vids.size() != 4) return false;
    if (mesh.V.at(vid).N_Vids.size() == 2) {
        FixDoublet(vid);
        return false;
    }
    return true;
}

void VertexRotation::SetRanking(glm::dvec3 d) {
    CheckValidity();
}

void VertexRotation::PerformOperation() {
    CheckValidity();
    if (!IsOperationValid()) return;

    std::cout << "Performing Vertex Rotation for " << vid << std::endl;

    Vertex& v = mesh.V.at(vid);
    std::vector<bool> isVisited(v.N_Fids.size(), false);
    std::vector<std::vector<size_t>> newVs(v.N_Fids.size());
    std::vector<std::vector<size_t>> newEs(v.N_Fids.size());
    std::vector<Edge> newEdges;
    std::queue<size_t> q;
    q.push(0);
    while(!q.empty()) {
        size_t fidx = q.front();
        q.pop();
        isVisited.at(fidx) = true;
        Face& f = mesh.F.at(v.N_Fids.at(fidx));
        size_t vidx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), vid));
        std::vector<size_t> newV = {vid, f.Vids.at((vidx+2)%f.Vids.size()), f.Vids.at((vidx+3)%f.Vids.size())};
        size_t eid1, eid2, eid3, eid4, eid5;
        for (auto eid: f.Eids) {
            Edge& e = mesh.E.at(eid);
            if ((e.Vids.at(0) == newV.at(0) && e.Vids.at(1) == f.Vids.at((vidx+1)%f.Vids.size())) || (e.Vids.at(1) == newV.at(0) && e.Vids.at(0) == f.Vids.at((vidx+1)%f.Vids.size()))) {
                eid1 = eid;
            }
            if ((e.Vids.at(0) == newV.at(1) && e.Vids.at(1) == f.Vids.at((vidx+1)%f.Vids.size())) || (e.Vids.at(1) == newV.at(1) && e.Vids.at(0) == f.Vids.at((vidx+1)%f.Vids.size()))) {
                eid2 = eid;
            }
            if ((e.Vids.at(0) == newV.at(0) && e.Vids.at(1) == newV.at(2)) || (e.Vids.at(1) == newV.at(0) && e.Vids.at(0) == newV.at(2))) {
                eid3 = eid;
            }
            if ((e.Vids.at(0) == newV.at(1) && e.Vids.at(1) == newV.at(2)) || (e.Vids.at(1) == newV.at(1) && e.Vids.at(0) == newV.at(2))) {
                eid4 = eid;
            }
        }
        Edge& e1 = mesh.E.at(eid1);
        Edge& e2 = mesh.E.at(eid2);
        Edge& e3 = mesh.E.at(eid3);
        Edge& e4 = mesh.E.at(eid4);
        size_t fid2 = e3.N_Fids.at(0) == f.id ? fid2 = e3.N_Fids.at(1) : fid2 = e3.N_Fids.at(0);
        Face& f2 = mesh.F.at(fid2);
        size_t vidx2 = std::distance(f2.Vids.begin(), std::find(f2.Vids.begin(), f2.Vids.end(), newV.at(2)));
        newV.push_back(f2.Vids.at((vidx2+1)%f2.Vids.size()));
        for (auto eid: f2.Eids) {
            Edge& e = mesh.E.at(eid);
            if ((e.Vids.at(0) == newV.at(2) && e.Vids.at(1) == newV.at(3)) || (e.Vids.at(1) == newV.at(2) && e.Vids.at(0) == newV.at(3))) {
                eid5 = eid;
            }
        }
        Edge& e5 = mesh.E.at(eid5);
        UpdateContents(mesh.V.at(f.Vids.at((vidx+1)%f.Vids.size())).N_Vids, std::vector<size_t>{vid});
        UpdateContents(mesh.V.at(f.Vids.at((vidx+1)%f.Vids.size())).N_Eids, std::vector<size_t>{e1.id});
        UpdateContents(mesh.V.at(f.Vids.at((vidx+1)%f.Vids.size())).N_Fids, std::vector<size_t>{f.id});

        UpdateContents(v.N_Vids, std::vector<size_t>{f.Vids.at((vidx+1)%f.Vids.size())});
        AddContents(v.N_Vids, std::vector<size_t>{newV.at(1)});

        Edge newEdge(std::vector<size_t>{vid, newV.at(1)});
        newEdge.id = e1.id;
        newEdges.push_back(newEdge);

        AddContents(mesh.V.at(newV.at(1)).N_Vids, std::vector<size_t>{vid});                               
        AddContents(mesh.V.at(newV.at(1)).N_Eids, std::vector<size_t>{e1.id});

        AddContents(mesh.V.at(newV.at(3)).N_Fids, std::vector<size_t>{f.id});

        UpdateContents(e5.N_Fids, std::vector<size_t>{f2.id});
        AddContents(e5.N_Fids, std::vector<size_t>{f.id});

        newVs.at(fidx) = newV;
        newEs.at(fidx) = {e1.id, e3.id, e4.id, e5.id};                             
        
        size_t nfidx = std::distance(v.N_Fids.begin(), std::find(v.N_Fids.begin(), v.N_Fids.end(), f2.id));
        if (!isVisited.at(nfidx)) {
            q.push(nfidx);
        }
    }
    std::vector<size_t> verticesToCheck;
    for (int i = 0; i < v.N_Fids.size(); i++) {
        mesh.F.at(v.N_Fids.at(i)).Vids = newVs.at(i);
        mesh.F.at(v.N_Fids.at(i)).Eids = newEs.at(i);
        AddContents(verticesToCheck, newVs.at(i));
    }
    for (auto& e: newEdges) {
        mesh.E.at(e.id).Vids = e.Vids;
    }
    for (auto id: verticesToCheck) {
        FixDoublet(id);
    }
}
