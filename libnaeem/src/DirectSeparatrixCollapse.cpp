#include "DirectSeparatrixCollapse.h"

DirectSeparatrixCollapse::DirectSeparatrixCollapse() : SimplificationOperation() {}
DirectSeparatrixCollapse::DirectSeparatrixCollapse(Mesh& mesh_, MeshUtil& mu_, Smoother& smoother_, size_t cid_, std::vector<size_t> s1_, std::vector<size_t> s2_, bool looseCollapse_) : SimplificationOperation(mesh_, mu_, smoother_) {
    cid = cid_;
    s1 = s1_;
    s2 = s2_;
    looseCollapse = looseCollapse_;
}

DirectSeparatrixCollapse::~DirectSeparatrixCollapse() {}

void DirectSeparatrixCollapse::SetRanking(glm::dvec3 d) {
    CheckValidity();
    if (!IsOperationValid()) {
        ranking = -1;
        return;
    }
    auto& centerV = mesh.V.at(cid);

    double min = 0.0;
    double max = 0.0;
    std::vector<size_t> sepVertices = {s1.at(0), s1.at(1), s2.at(0), s2.at(1)};
    std::vector<size_t> maxVertices;
    for (auto id: s1) {
        auto& v = mesh.V.at(id);
        auto& f = mesh.F.at(GetDifference(v.N_Fids, centerV.N_Fids).at(0));
        int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), id));
        // min += mu.GetVertexEnergy(f.Vids.at((idx + 2) % f.Vids.size()));

        sepVertices.push_back(f.Vids.at((idx + 2) % f.Vids.size()));
        auto& v1 = mesh.V.at(f.Vids.at((idx + 2) % f.Vids.size()));
        max += v1.N_Vids.size() == 4 ? v1.N_Vids.size() * 4 : v1.N_Vids.size();  

        std::vector<size_t> maxV = GetDifference(mesh.V.at(f.Vids.at((idx + 1) % f.Vids.size())).N_Vids, sepVertices);
        if (maxV.size() > 0) std::move(maxV.begin(), maxV.begin()+1, std::back_inserter(maxVertices));
        maxV.clear();
        maxV = GetDifference(mesh.V.at(f.Vids.at((idx + 3) % f.Vids.size())).N_Vids, sepVertices);
        if (maxV.size() > 0) std::move(maxV.begin(), maxV.begin()+1, std::back_inserter(maxVertices));
    }

    for (auto id: maxVertices) {
        auto& v = mesh.V.at(id);
        max += v.N_Vids.size() == 4 ? v.N_Vids.size() * 2 : v.N_Vids.size();
    }

    min += mesh.V.at(s2.at(0)).N_Vids.size() + mesh.V.at(s2.at(1)).N_Vids.size();

    ranking = max / min;

    // if (glm::length(d) > 0) {
    //     ranking /= GetDistance(d);
    // }

    // double max = mu.GetVertexEnergy(s2.at(0)) + mu.GetVertexEnergy(s1.at(1));

    // double normalized_area = 0.0;
    // std::vector<size_t> fids;
    // AddContents(fids, mesh.V.at(s1.at(0)).N_Fids);
    // AddContents(fids, mesh.V.at(s1.at(1)).N_Fids);
    // for (auto fid: fids) {
    //     normalized_area += mu.GetFaceArea(fid);
    // }
    // normalized_area /= mu.GetMeshArea();

    // ranking = (min / max) * normalized_area;
    // ranking = 1;
} 

bool DirectSeparatrixCollapse::IsOperationValid() {
    CheckValidity();

    bool isValid = true;
    auto& v = mesh.V.at(cid);
    if (v.N_Fids.size() == 0) isValid = false;
    if (mesh.V.at(s1.at(0)).N_Fids.size() != 3 || mesh.V.at(s1.at(1)).N_Fids.size() != 3) isValid = false;
    if (looseCollapse) {
        if (mesh.V.at(s2.at(0)).N_Fids.size() == 4 && mesh.V.at(s2.at(1)).N_Fids.size() == 4) isValid = false;
        if ((mesh.V.at(s2.at(0)).type == FEATURE || mesh.V.at(s2.at(0)).isBoundary) && mesh.V.at(s2.at(0)).N_Fids.size() == 4) isValid = false;
        if ((mesh.V.at(s2.at(1)).type == FEATURE || mesh.V.at(s2.at(1)).isBoundary) && mesh.V.at(s2.at(1)).N_Fids.size() == 4) isValid = false;
    } else {
        if (mesh.V.at(s2.at(0)).N_Fids.size() == 4 || mesh.V.at(s2.at(1)).N_Fids.size() == 4) isValid = false;
    }
    // for (auto id: s2) {
    //     auto& s = mesh.V.at(id);
    //     if (s.type != FEATURE) continue;
    //     for (auto eid: v.N_Eids) {
    //         auto& e = mesh.E.at(eid);
    //         if (e.Vids.at(0) == id || e.Vids.at(1) == 1) {
    //             int count = 0;
    //             for (auto fid: e.N_Fids) {
    //                 auto& f = mesh.F.at(fid);
    //                 for (auto vid: f.Vids) {
    //                     if (mesh.V.at(vid).type == FEATURE) count += 1;
    //                 }
    //             }
    //             if (count == 4) isValid = false;
    //         }
    //     }
    // }
    return isValid;
}

void DirectSeparatrixCollapse::PerformOperation() {
    CheckValidity();
    
    // centerV is the center vertex on which two neighboring 3-singularities are collapsed 
    auto& centerV = mesh.V.at(cid);    
    // for (auto nvid: centerV.N_Vids) SetUpdateElements(nvid);

    if (!IsOperationValid()) return;

    // std::cout << "Performing Direct Separatrix Collapse: " << ranking << std::endl;

    // collect vertices to update connected operations
    for (auto id: s1) {
        auto& v = mesh.V.at(id);
        std::vector<size_t> s1F = GetDifference(v.N_Fids, centerV.N_Fids);
        if (s1F.empty()) continue;
        auto& f = mesh.F.at(s1F.at(0));
        int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), id));
        SetUpdateElements(f.Vids.at((idx + 3) % f.Vids.size()));
    }
    for (auto id: s2) {
        SetUpdateElements(id);
    }
    if (ranking < 0) return;

    std::vector<size_t> prev_v = centerV.N_Vids;
    std::vector<size_t> prev_e = centerV.N_Eids;
    std::vector<size_t> prev_f = centerV.N_Fids;
    centerV.N_Vids.clear();
    centerV.N_Eids.clear();
    centerV.N_Fids.clear();

    // the edges connecting n-singularities to the center vertex
    std::vector<size_t> c_e;
    for (auto eid: prev_e) {
        auto& e = mesh.E.at(eid);
        int idx = (std::distance(e.Vids.begin(), std::find(e.Vids.begin(), e.Vids.end(), centerV.id)) + 1) % e.Vids.size();
        if (e.Vids.at(idx) == s1.at(0) || e.Vids.at(idx) == s1.at(1)) continue;
        c_e.push_back(e.id);
    }
    
    // perform collapsing operations using the extracted edges and update neighboring info
    for (auto id: c_e) {
        auto& e = mesh.E.at(id);
        auto& faceToKeep = mesh.F.at(e.N_Fids.at(0)); // one face is kept and its vertices are updated
        auto& faceToRemove = mesh.F.at(e.N_Fids.at(1)); // one face is removed as part of simplification
        
        // centerV.N_Fids.push_back(faceToKeep.id); // add faces kept to the center vertex neighborhood
        AddContents(centerV.N_Fids, std::vector<size_t>{faceToKeep.id});

        // edges kept in the above two faces to form new connections
        std::vector<size_t> edgesToKeep;
        AddContents(edgesToKeep, GetDifference(faceToKeep.Eids, prev_e));
        AddContents(edgesToKeep, GetDifference(faceToRemove.Eids, prev_e));
        
        std::vector<size_t> newVids; // initialize new vertices for the faces kept during simplification
        newVids.push_back(centerV.id);
        
        // updated neighboring info using extracted edges
        for (auto eid: edgesToKeep) {
            auto& edgeToKeep = mesh.E.at(eid);
            
            // edge of the face kept that contains a 3-singularity
            int idx = -1;
            edgeToKeep.Vids.at(0) == s1.at(0) ? idx = 0 : 1;
            edgeToKeep.Vids.at(0) == s1.at(1) ? idx = 0 : 1;
            edgeToKeep.Vids.at(1) == s1.at(0) ? idx = 1 : 1;
            edgeToKeep.Vids.at(1) == s1.at(1) ? idx = 1 : 1;
            if (idx < 0) {
                edgeToKeep.N_Fids.at(0) == faceToRemove.id ? edgeToKeep.N_Fids.at(0) = faceToKeep.id : 1;
                edgeToKeep.N_Fids.size() > 1 && edgeToKeep.N_Fids.at(1) == faceToRemove.id ? edgeToKeep.N_Fids.at(1) = faceToKeep.id : 1;
                continue;
            }
            
            // opposite vertex of the edge connecting to 3-singularity
            int idx2 = (idx+1)%edgeToKeep.Vids.size();
            auto& etk_v = mesh.V.at(edgeToKeep.Vids.at(idx2));
            
            // updates related to center vertex
            edgeToKeep.Vids.at(idx) = centerV.id;
            AddContents(centerV.N_Vids, std::vector<size_t>{etk_v.id});
            AddContents(centerV.N_Eids, std::vector<size_t>{edgeToKeep.id});
            // centerV.N_Vids.push_back(etk_v.id);
            // centerV.N_Eids.push_back(edgeToKeep.id);
            for (int i = 0; i < etk_v.N_Vids.size(); i++) {
                if (etk_v.N_Vids.at(i) == s1.at(0) || etk_v.N_Vids.at(i) == s1.at(1)) etk_v.N_Vids.at(i) = centerV.id;
            }

            // updates related to opposite vertex connecting to 3-singularity
            UpdateContents(etk_v.N_Fids, std::vector<size_t>{faceToRemove.id});
            AddContents(etk_v.N_Fids, std::vector<size_t>{faceToKeep.id});
            std::vector<size_t> etk_f = GetDifference(etk_v.N_Fids, std::vector<size_t>{faceToKeep.id});
            for (auto nfid: etk_f) {
                auto& nf = mesh.F.at(nfid);
                UpdateContents(nf.N_Fids, std::vector<size_t>{faceToRemove.id});
                AddContents(nf.N_Fids, std::vector<size_t>{faceToKeep.id});
            }

            // updates related to face kept during simplification process
            if (std::find(faceToKeep.Eids.begin(), faceToKeep.Eids.end(), edgeToKeep.id) != faceToKeep.Eids.end()) {
                edgeToKeep.N_Fids.at(0) == faceToKeep.id ? AddContents(centerV.N_Fids, std::vector<size_t>{edgeToKeep.N_Fids.at(1)}) : AddContents(centerV.N_Fids, std::vector<size_t>{edgeToKeep.N_Fids.at(0)}); 
                
                int nv_idx = std::distance(faceToKeep.Vids.begin(), std::find(faceToKeep.Vids.begin(), faceToKeep.Vids.end(), edgeToKeep.Vids.at(idx2)));
                if (s1.at(0) == faceToKeep.Vids.at((nv_idx+1)%faceToKeep.Vids.size()) || s1.at(1) == faceToKeep.Vids.at((nv_idx+1)%faceToKeep.Vids.size())) {
                    newVids.insert(newVids.begin(), faceToKeep.Vids.at(nv_idx));
                } else if (s1.at(0) == faceToKeep.Vids.at((nv_idx+3)%faceToKeep.Vids.size()) || s1.at(1) == faceToKeep.Vids.at((nv_idx+3)%faceToKeep.Vids.size())) {
                    newVids.push_back(faceToKeep.Vids.at(nv_idx));
                }
            } else {
                edgeToKeep.N_Fids.at(0) == faceToRemove.id ? AddContents(centerV.N_Fids, std::vector<size_t>{edgeToKeep.N_Fids.at(1)}) : AddContents(centerV.N_Fids, std::vector<size_t>{edgeToKeep.N_Fids.at(0)});
                edgeToKeep.N_Fids.at(0) == faceToRemove.id ? edgeToKeep.N_Fids.at(0) = faceToKeep.id : edgeToKeep.N_Fids.at(1) = faceToKeep.id;
                
                int nv_idx = std::distance(faceToRemove.Vids.begin(), std::find(faceToRemove.Vids.begin(), faceToRemove.Vids.end(), edgeToKeep.Vids.at(idx2)));
                if (s1.at(0) == faceToRemove.Vids.at((nv_idx+1)%faceToRemove.Vids.size()) || s1.at(1) == faceToRemove.Vids.at((nv_idx+1)%faceToRemove.Vids.size())) {
                    newVids.insert(newVids.begin(), faceToRemove.Vids.at(nv_idx));
                } else if (s1.at(0) == faceToRemove.Vids.at((nv_idx+3)%faceToRemove.Vids.size()) || s1.at(1) == faceToRemove.Vids.at((nv_idx+3)%faceToRemove.Vids.size())) {
                    newVids.push_back(faceToRemove.Vids.at(nv_idx));
                }
            }
        }
        std::find(faceToKeep.Vids.begin(), faceToKeep.Vids.end(), s2.at(0)) != faceToKeep.Vids.end() ? newVids.push_back(s2.at(0)) : newVids.push_back(s2.at(1));
        faceToKeep.Vids = newVids;
        faceToKeep.Eids = edgesToKeep;

        // updates related to removed face and neighboring elements
        faceToRemove.Vids.clear();
        faceToRemove.Eids.clear();
        faceToRemove.N_Fids.clear();
        UpdateContents(mesh.V.at(s2.at(0)).N_Fids, std::vector<size_t>{faceToRemove.id});
        UpdateContents(mesh.V.at(s2.at(1)).N_Fids, std::vector<size_t>{faceToRemove.id});
        for (auto fid: centerV.N_Fids) {
            auto& f = mesh.F.at(fid);
            for (int i = 0; i < f.Vids.size(); i++) {
                if (f.Vids.at(i) == s1.at(0) || f.Vids.at(i) == s1.at(1)) f.Vids.at(i) = centerV.id;
            }
            UpdateContents(f.N_Fids, std::vector<size_t>{faceToRemove.id});
            AddContents(f.N_Fids, GetDifference(centerV.N_Fids, std::vector<size_t>{f.id}));
        }
    }

    // clear 3-singularities
    mesh.V.at(s1.at(0)).N_Vids.clear();
    mesh.V.at(s1.at(0)).N_Eids.clear();
    mesh.V.at(s1.at(0)).N_Fids.clear();
    mesh.V.at(s1.at(1)).N_Vids.clear();
    mesh.V.at(s1.at(1)).N_Eids.clear();
    mesh.V.at(s1.at(1)).N_Fids.clear();
    for (auto eid: prev_e) {
        mesh.E.at(eid).Vids.clear();
        mesh.E.at(eid).N_Fids.clear();
    }

    
    // update n-singularities
    auto& s2_1 = mesh.V.at(s2.at(0));
    auto& s2_2 = mesh.V.at(s2.at(1));
    UpdateContents(s2_1.N_Vids, std::vector<size_t>{centerV.id});
    UpdateContents(s2_2.N_Vids, std::vector<size_t>{centerV.id});
    UpdateContents(s2_1.N_Eids, prev_e);
    UpdateContents(s2_2.N_Eids, prev_e);

    SetSingularity(centerV.id);
    SetSingularity(s2_1.id);
    SetSingularity(s2_2.id);
    
    // collect vertices to smooth
    // ToSmooth.push_back(centerV.id);
    // ToSmooth.push_back(s2_1.id);
    // ToSmooth.push_back(s2_2.id);
    AddContents(ToSmooth, std::vector<size_t>{centerV.id, s2_1.id, s2_2.id});
    AddContents(ToSmooth, s2_1.N_Vids);
    AddContents(ToSmooth, s2_2.N_Vids);
    // ToSmooth.insert(ToSmooth.end(), s2_1.N_Vids.begin(), s2_1.N_Vids.end());
    // ToSmooth.insert(ToSmooth.end(), s2_2.N_Vids.begin(), s2_2.N_Vids.end());
    Smooth();
}

void DirectSeparatrixCollapse::SetUpdateElements(size_t vid) {
    auto& v = mesh.V.at(vid);
    for (auto nvid: v.N_Vids) {
        if (mesh.V.at(nvid).N_Vids.size() == 4 && mesh.V.at(nvid).type != FEATURE) AddContents(toUpdate, std::vector<size_t>{nvid});
    }
    for (auto fid: v.N_Fids) {
        auto& f = mesh.F.at(fid);
        int index = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
        auto& diagV = mesh.V.at(f.Vids.at((index+2)%f.Vids.size()));
        for (auto nvid: diagV.N_Vids) {
            if (mesh.V.at(nvid).N_Vids.size() == 4 && mesh.V.at(nvid).type != FEATURE) AddContents(toUpdate, std::vector<size_t>{nvid});
        }
    }
}

glm::dvec3 DirectSeparatrixCollapse::GetLocation() {
    return mesh.V.at(cid).xyz();
}

double DirectSeparatrixCollapse::GetDistance(glm::dvec3 a) {
    return glm::length(mesh.V.at(cid).xyz() - a);
}