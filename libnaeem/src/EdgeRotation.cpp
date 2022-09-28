#include "EdgeRotation.h"

EdgeRotation::EdgeRotation() {}
EdgeRotation::EdgeRotation(Mesh& mesh_, MeshUtil& mu_, Smoother& smoother_, size_t eid_, bool clockwise_) : SimplificationOperation(mesh_, mu_, smoother_) {
    eid = eid_;
    clockwise = clockwise_;
}
EdgeRotation::~EdgeRotation() {}

bool EdgeRotation::IsOperationValid() {
    CheckValidity();
    Edge& e = mesh->E.at(eid);
    if (e.Vids.empty() || e.N_Fids.empty()) return false;
    for (auto fid: e.N_Fids) {
        Face& f = mesh->F.at(fid);
        for (auto vid: f.Vids) {
            if (mesh->V.at(vid).N_Vids.size() == 2) return false;
        }
    }
    return true;
}

void EdgeRotation::SetRanking(glm::dvec3 d) {
    CheckValidity();
}

void EdgeRotation::PerformOperation() {
    CheckValidity();
    if (!IsOperationValid()) return;

    // std::cout << "Performing Edge Rotation for " << eid << std::endl;

    Edge& e = mesh->E.at(eid);
    size_t v1 = e.Vids.at(0);
    size_t v2 = e.Vids.at(1);
    std::vector<size_t> vertices = GetVertices(e, v1, v2);
    int offset_a = 1;
    int offset_b = 1;
    int offset_c = 0;
    int offset_d = 1;
    if (clockwise) {
        offset_a = 2;
        offset_b = 0;
        offset_c = 3;
        offset_d = 2;
    }
    size_t newV_1 = vertices.at((std::distance(vertices.begin(), std::find(vertices.begin(), vertices.end(), v1)) + offset_a)%vertices.size());
    size_t newV_2 = vertices.at((std::distance(vertices.begin(), std::find(vertices.begin(), vertices.end(), v2)) + offset_a)%vertices.size());
    size_t nextV_1 = vertices.at((std::distance(vertices.begin(), std::find(vertices.begin(), vertices.end(), v1)) + offset_d)%vertices.size());
    size_t nextV_2 = vertices.at((std::distance(vertices.begin(), std::find(vertices.begin(), vertices.end(), v2)) + offset_d)%vertices.size());
    size_t edgeToUpdate1;
    size_t edgeToUpdate2;
    for (auto fid: e.N_Fids) {
        auto& f = mesh->F.at(fid);
        for (auto feid: f.Eids) {
            Edge& fe = mesh->E.at(feid);
            size_t v_a = v1;
            size_t v_b = v2;
            if (clockwise) {
                v_a = v2;
                v_b = v1;
            }
            if (std::find(fe.Vids.begin(), fe.Vids.end(), v_a) != fe.Vids.end() && std::find(fe.Vids.begin(), fe.Vids.end(), nextV_1) != fe.Vids.end()) {
                edgeToUpdate1 = fe.id;
            }
            if (std::find(fe.Vids.begin(), fe.Vids.end(), v_b) != fe.Vids.end() && std::find(fe.Vids.begin(), fe.Vids.end(), nextV_2) != fe.Vids.end()) {
                edgeToUpdate2 = fe.id;
            }
        }    
    }
    for (auto fid: e.N_Fids) {
        auto& f = mesh->F.at(fid);
        int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v1));
        if (f.Vids.at((idx+1)%f.Vids.size()) == v2) {
            UpdateContents(mesh->V.at(f.Vids.at((idx+offset_b)%f.Vids.size())).N_Fids, std::vector<size_t>{fid});
            AddContents(mesh->V.at(newV_1).N_Fids, std::vector<size_t>{fid});
            f.Vids.at((idx+offset_b)%f.Vids.size()) = newV_1;
            for (auto nfid: mesh->V.at(newV_1).N_Fids) {
                AddContents(mesh->F.at(nfid).N_Fids, std::vector<size_t>{fid});
            }
            UpdateContents(f.Eids, std::vector<size_t>{edgeToUpdate2});
            UpdateContents(mesh->E.at(edgeToUpdate2).N_Fids, std::vector<size_t>{f.id});
            AddContents(f.Eids, std::vector<size_t>{edgeToUpdate1});
            AddContents(mesh->E.at(edgeToUpdate1).N_Fids, std::vector<size_t>{f.id});
        } else if (f.Vids.at((idx+3)%f.Vids.size()) == v2) {
            UpdateContents(mesh->V.at(f.Vids.at((idx+offset_c)%f.Vids.size())).N_Fids, std::vector<size_t>{fid});
            AddContents(mesh->V.at(newV_2).N_Fids, std::vector<size_t>{fid});
            f.Vids.at((idx+offset_c)%f.Vids.size()) = newV_2;
            for (auto nfid: mesh->V.at(newV_2).N_Fids) {
                AddContents(mesh->F.at(nfid).N_Fids, std::vector<size_t>{fid});
            }
            UpdateContents(f.Eids, std::vector<size_t>{edgeToUpdate1});
            UpdateContents(mesh->E.at(edgeToUpdate1).N_Fids, std::vector<size_t>{f.id});
            AddContents(f.Eids, std::vector<size_t>{edgeToUpdate2});
            AddContents(mesh->E.at(edgeToUpdate2).N_Fids, std::vector<size_t>{f.id});
        }
    }
    UpdateContents(mesh->V.at(e.Vids.at(0)).N_Vids, std::vector<size_t>{e.Vids.at(1)});
    UpdateContents(mesh->V.at(e.Vids.at(1)).N_Vids, std::vector<size_t>{e.Vids.at(0)});
    UpdateContents(mesh->V.at(e.Vids.at(0)).N_Eids, std::vector<size_t>{e.id});
    UpdateContents(mesh->V.at(e.Vids.at(1)).N_Eids, std::vector<size_t>{e.id});
    AddContents(mesh->V.at(newV_1).N_Vids, std::vector<size_t>{newV_2});
    AddContents(mesh->V.at(newV_2).N_Vids, std::vector<size_t>{newV_1});
    AddContents(mesh->V.at(newV_1).N_Eids, std::vector<size_t>{e.id});
    AddContents(mesh->V.at(newV_2).N_Eids, std::vector<size_t>{e.id});
    e.Vids.at(0) = newV_1;
    e.Vids.at(1) = newV_2;
    // std::cout << "Setting Singularities" << std::endl;
    for (auto fid: e.N_Fids) {
        auto& f = mesh->F.at(fid);
        for (auto fvid: f.Vids) {
            SetSingularity(fvid);
            auto& fv = mesh->V.at(fvid);
            // ToSmooth.push_back(fvid);
            AddContents(ToSmooth, std::vector<size_t>{fvid});
            AddContents(ToSmooth, fv.N_Vids);
            // ToSmooth.insert(ToSmooth.end(), fv.N_Vids.begin(), fv.N_Vids.end());
        }
    }
    for (auto fid: e.N_Fids) {
        Face& f = mesh->F.at(fid);
        for (auto fvid: f.Vids) {
            FixDoublet(fvid);
        }
    }
    // std::cout << "Finished Edge Rotation" << std::endl;
    Smooth();
}

std::vector<size_t> EdgeRotation::GetVertices(Edge& e, size_t v1, size_t v2) {
    std::vector<size_t> vertices;
    for (auto fid: e.N_Fids) {
        auto& f = mesh->F.at(fid);
        int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v1));
        if (f.Vids.at((idx+1)%f.Vids.size()) == v2) {
            vertices.push_back(f.Vids.at((idx+1)%f.Vids.size()));
            vertices.push_back(f.Vids.at((idx+2)%f.Vids.size()));
            vertices.push_back(f.Vids.at((idx+3)%f.Vids.size()));
        } else {
            vertices.push_back(f.Vids.at(idx));
            vertices.push_back(f.Vids.at((idx+1)%f.Vids.size()));
            vertices.push_back(f.Vids.at((idx+2)%f.Vids.size()));
        }
    }
    return vertices;
}