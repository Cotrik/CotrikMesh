/*
 * EdgeRotateSimplifier.cpp
 *
 *  Created on: Jan 20, 2019
 *      Author: cotrik
 */

#include "EdgeRotateSimplifier.h"

EdgeRotateSimplifier::EdgeRotateSimplifier(Mesh& mesh) : Simplifier(mesh) {
    // TODO Auto-generated constructor stub

}

EdgeRotateSimplifier::~EdgeRotateSimplifier() {
    // TODO Auto-generated destructor stub
}

void EdgeRotateSimplifier::Combine(std::vector<size_t>& com, std::vector<std::vector<size_t>> &res, int n, int k, int start) {
    if (k == com.size()) {
        res.push_back(com);
        return;
    }
    for (int i = start; i < n; ++i) {
        com.push_back(i);
        Combine(com, res, n, k, i + 1);
        com.pop_back();
    }
}

std::vector<std::vector<size_t>> EdgeRotateSimplifier::Combine(int n, int k) {
    std::vector<std::vector<size_t>> res;
    std::vector<size_t> com;
    Combine(com, res, n, k, 0);
    return res;
}

std::map<size_t, std::set<size_t>> EdgeRotateSimplifier::GetPatchid_Fids() {
    std::map<size_t, std::set<size_t>> patchid_fids;
    for (auto& v : mesh.V)
        if (v.patch_id != MAXID) patchid_fids[v.patch_id].insert(v.N_Fids.begin(), v.N_Fids.end());
    return patchid_fids;
}

std::map<size_t, std::set<size_t>> EdgeRotateSimplifier::GetPatchid_Vids(const std::map<size_t, std::set<size_t>>& patchid_fids) {
    std::map<size_t, std::set<size_t>> patchid_vids;
    for (auto& item : patchid_fids) {
        for (auto fid : item.second) {
            auto& f = mesh.F.at(fid);
            patchid_vids[item.first].insert(f.Vids.begin(), f.Vids.end());
        }
    }
    return patchid_vids;
}

std::set<size_t> EdgeRotateSimplifier::GetRotateFids() {
    std::set<size_t> res;
    auto patchid_fids = GetPatchid_Fids();
    auto patchid_vids = GetPatchid_Vids(patchid_fids);

    for (auto& item : patchid_vids) {
        for (auto vid : item.second) {
            auto& v = mesh.V.at(vid);
            if (v.type == CORNER && !v.isSpecial) {
                auto neighbor_f_count = 0;
                auto& patch_fids = patchid_fids[item.first];
                for (auto fid : v.N_Fids)
                    if (patch_fids.find(fid) != patch_fids.end()) ++neighbor_f_count;
                if (neighbor_f_count != 2) continue;
                for (auto fid : v.N_Fids)
                    if (patch_fids.find(fid) != patch_fids.end()) res.insert(fid);
            } else if (v.type == FEATURE) {
                auto neighbor_f_count = 0;
                auto& patch_fids = patchid_fids[item.first];
                for (auto fid : v.N_Fids)
                    if (patch_fids.find(fid) != patch_fids.end()) ++neighbor_f_count;
                if (neighbor_f_count <= 2) continue;
                for (auto fid : v.N_Fids)
                    if (patch_fids.find(fid) != patch_fids.end()) res.insert(fid);
            }
        }
    }

    return res;
}

std::vector<size_t> EdgeRotateSimplifier::GetNeighborFids(const Vertex& v, const std::set<size_t>& patch_fids) {
    std::vector<size_t> fids;
    for (auto fid : v.N_Fids)
        if (patch_fids.find(fid) != patch_fids.end()) fids.push_back(fid);
    return fids;
}

std::set<size_t> EdgeRotateSimplifier::GetRotateEids(const Vertex& v, const std::vector<size_t>& fids) {
    std::set<size_t> res;
    auto pairs = Util::combine(fids.size(), 2);
    for (auto& p : pairs) {
        auto& f0 = mesh.F.at(fids.at(p[0]));
        auto& f1 = mesh.F.at(fids.at(p[1]));
        std::set<size_t> f0eids(f0.Eids.begin(), f0.Eids.end());
        std::set<size_t> f1eids(f1.Eids.begin(), f1.Eids.end());
        auto intersect_eids = Util::get_intersect(f0eids, f1eids);
        if (intersect_eids.size() == 1) {
            auto& e = mesh.E.at(*intersect_eids.begin());
            auto vid = (v.id == e.Vids[0]) ? e.Vids[1] : e.Vids[0];
            if (mesh.V.at(vid).type == CORNER) continue;
            res.insert(e.id);
        }
    }
    return res;
}

double EdgeRotateSimplifier::GetAngle(const Vertex& v, const Vertex& v0, const Vertex& v1) {
    auto d0 = v0.xyz() - v.xyz();
    auto d1 = v1.xyz() - v.xyz();
    auto cosangle = glm::dot(glm::normalize(d0), glm::normalize(d1));
    return acos(cosangle) * 180.0 / PI;
}

std::vector<size_t> EdgeRotateSimplifier::GetNeighborVids(const Vertex& v, size_t fid) {
    std::vector<size_t> res;
    auto& f = mesh.F.at(fid);
    for (auto eid : f.Eids) {
        auto& e = mesh.E.at(eid);
        if (e.Vids[0] == v.id || e.Vids[1] == v.id) {
            auto vid = e.Vids[0] == v.id ? e.Vids[1] : e.Vids[0];
            res.push_back(vid);
        }
    }
    return res;
}

bool EdgeRotateSimplifier::IsConvex(const Vertex& v, const std::vector<size_t>& fids) {
    double total_angle = 0;
    for (auto fid : fids) {
        auto vids = get_neighbor_vids(v, fid);
        total_angle += GetAngle(v, mesh.V.at(vids[0]), mesh.V.at(vids[1]));
    }
    return total_angle < 135;
}

std::vector<size_t> EdgeRotateSimplifier::GetBoundaryEids(const Face& f0, const Face& f1, const Edge& exclude_e) {
    std::set<size_t> eids_set(f0.Eids.begin(), f0.Eids.end());
    eids_set.insert(f1.Eids.begin(), f1.Eids.end());
    eids_set.erase(exclude_e.id);
    std::vector<size_t> eids;
    std::copy(eids_set.begin(), eids_set.end(), std::back_inserter(eids));
    return eids;
}

std::vector<size_t> EdgeRotateSimplifier::GetRotateVids(const std::vector<size_t>& boundary_eids, size_t diag_vid0, size_t diag_vid1) {
    auto linkVids = GetLinkVidsFromEids(boundary_eids);
    linkVids.pop_back();
    while (linkVids.front() != diag_vid0 && linkVids.front() != diag_vid1) {
        linkVids.push_back(linkVids.front());
        linkVids.erase(linkVids.begin());
    }
    return linkVids;
}

std::set<size_t> EdgeRotateSimplifier::GetRotateEids() {
    std::set<size_t> res;
    auto patchid_fids = get_patchid_fids();
    auto patchid_vids = get_patchid_vids(patchid_fids);
    for (auto& item : patchid_vids) {
        for (auto vid : item.second) {
            auto& v = mesh.V.at(vid);
            if (v.type == CORNER && !v.isSpecial) {
                auto fids = GetNeighborFids(v, patchid_fids[item.first]);
                if (fids.size() < 2) continue;
                if (!is_convex(v, fids) && fids.size() < 4) continue;
                auto eids = GetRotateEids(v, fids);
                res.insert(eids.begin(), eids.end());
            } else if (v.type == FEATURE) {
                auto fids = GetNeighborFids(v, patchid_fids[item.first]);
                if (fids.size() < 3 || fids.size() == 4) continue;
                auto eids = GetRotateEids(v, fids);
                res.insert(eids.begin(), eids.end());
            }
        }
    }

    return res;
}

size_t EdgeRotateSimplifier::GetDiagnalVid(size_t vid, size_t fid) {
    auto& f = mesh.F.at(fid);
    for (int i = 0; i < 4; ++i) {
        if (f.Vids[i] == vid) return f.Vids.at((i + 2) % 4);
    }
    MeshFileWriter writer(mesh, "error.vtk");
    writer.WriteFile();
    std::cerr << "ERROR get_diagnal_vid\n";
    return MAXID;
}

void EdgeRotateSimplifier::Rotate(std::set<size_t>& canceledFids) {
    auto patchid_fids = GetPatchid_Fids();
    auto patchid_vids = GetPatchid_Vids(patchid_fids);
    for (auto& item : patchid_vids) {
        for (auto vid : item.second) {
            auto& v = mesh.V.at(vid);
            if ((v.type == CORNER || v.isCorner) && !v.isSpecial) {
                // v.type = CORNER;
                auto fids = GetNeighborFids(v, patchid_fids[item.first]);
                if (fids.size() < 2) continue;
                if (!is_convex(v, fids) && fids.size() < 4) continue;
                auto eids = GetRotateEids(v, fids);
                auto& e = mesh.E.at(*eids.begin());
                Rotate(e, v, canceledFids);
                return;
            } else if (v.type == FEATURE) {
                auto fids = GetNeighborFids(v, patchid_fids[item.first]);
                if (fids.size() < 3/* || fids.size() == 4*/) continue;
                auto eids = GetRotateEids(v, fids);
                auto& e = mesh.E.at(*eids.begin());
                Rotate(e, v, canceledFids);
                return;
            }
        }
    }
}

void EdgeRotateSimplifier::Rotate(const Edge& e, const Vertex& v, std::set<size_t>& canceledFids) {
    auto& f0 = mesh.F.at(e.N_Fids[0]);
    auto& f1 = mesh.F.at(e.N_Fids[1]);
    auto vid0 = e.Vids[0] == v.id ? e.Vids[1] : e.Vids[0];
    canceledFids.insert(f0.id);
    canceledFids.insert(f1.id);
    auto diag_vid0 = GetDiagnalVid(vid0, f0.id);
    auto diag_vid1 = GetDiagnalVid(v.id, f1.id);
    if (mesh.V.at(diag_vid0).type >= FEATURE) {
        diag_vid0 = GetDiagnalVid(v.id, f0.id);
        diag_vid1 = GetDiagnalVid(vid0, f1.id);
    }
    auto eids = GetBoundaryEids(f0, f1, e);
    auto linkVids = GetRotateVids(eids, diag_vid0, diag_vid1);
    insert_rotate_faces(linkVids);
}

void EdgeRotateSimplifier::Run() {
    std::set<size_t> canceledFids;
    auto count = 0;
    while (true) {
        if (canceledFids.empty() && REMOVE_DOUBLET) {
            remove_doublet(canceledFids);
            if (!canceledFids.empty()) {
                std::cout << "remove_doublet" << std::endl;
                update(canceledFids);
                canceledFids.clear();
                init();
                continue;
            }
        }
        Rotate(canceledFids);
        if (canceledFids.empty()) break;
        std::cout << "rotate " << count++ << std::endl;
        update(canceledFids);
        canceledFids.clear();
        init();
        break;
    }
}

void EdgeRotateSimplifier::Run(std::set<size_t>& canceledFids) {
    auto patchid_fids = GetPatchid_Fids();
    auto patchid_vids = GetPatchid_Vids(patchid_fids);
    for (auto& item : patchid_vids) {
        for (auto vid : item.second) {
            auto& v = mesh.V.at(vid);
            if ((v.type == CORNER || v.isCorner) && !v.isSpecial) {
                // v.type = CORNER;
                auto fids = GetNeighborFids(v, patchid_fids[item.first]);
                //if (fids.size() < 2) continue;
                //if (!is_convex(v, fids) && fids.size() < 4) continue;
                if (v.idealValence >= fids.size()) continue;
                auto eids = GetRotateEids(v, fids);
                if (eids.empty()) continue;
                auto& e = mesh.E.at(*eids.begin());
                Rotate(e, v, canceledFids);
                return;
            } else if (v.type == FEATURE) {
                auto fids = GetNeighborFids(v, patchid_fids[item.first]);
                if (fids.size() < 3/* || fids.size() == 4*/) continue;
                auto eids = GetRotateEids(v, fids);
                if (eids.empty()) continue;
                auto& e = mesh.E.at(*eids.begin());
                Rotate(e, v, canceledFids);
                return;
            }
        }
    }
}
