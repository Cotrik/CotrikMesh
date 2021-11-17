#include "DirectSeparatrixCollapse.h"

DirectSeparatrixCollapse::DirectSeparatrixCollapse() : SimplificationOperation() {}
DirectSeparatrixCollapse::DirectSeparatrixCollapse(Mesh& mesh_, MeshUtil& mu_, size_t cid_, std::vector<size_t> s1_, std::vector<size_t> s2_) : SimplificationOperation(mesh_, mu_) {
    cid = cid_;
    s1 = s1_;
    s2 = s2_;
}

DirectSeparatrixCollapse::~DirectSeparatrixCollapse() {}

void DirectSeparatrixCollapse::SetRanking() {
    CheckValidity();

    auto& centerV = mesh.V.at(cid);

    double min = 0.0;
    for (auto id: s1) {
        auto& v = mesh.V.at(id);
        auto& f = mesh.F.at(GetDifference(v.N_Fids, centerV.N_Fids).at(0));
        size_t idx = 0;
        for (idx < f.Vids.size(); idx++;) {
            if (f.Vids.at(idx) == id) break;
        }
        min += mu.GetVertexEnergy(f.Vids.at((idx + 2) % f.Vids.size()));
    }

    double max = mu.GetVertexEnergy(s2.at(0)) + mu.GetVertexEnergy(s1.at(1));

    double normalized_area = 0.0;
    std::vector<size_t> fids;
    AddContents(fids, mesh.V.at(s1.at(0)).N_Fids);
    AddContents(fids, mesh.V.at(s1.at(1)).N_Fids);
    for (auto fid: fids) {
        normalized_area += mu.GetFaceArea(fid);
    }
    normalized_area /= mu.GetMeshArea();

    ranking = min / (max * normalized_area);
} 

bool DirectSeparatrixCollapse::IsOperationValid() {
    CheckValidity();

    if (mesh.V.at(cid).N_Fids.size() == 0) return false;
    if (mesh.V.at(s1.at(0)).N_Fids.size() != 3 || mesh.V.at(s1.at(1)).N_Fids.size() != 3) return false;

    return true;
}

void DirectSeparatrixCollapse::PerformOperation() {
    CheckValidity();

    if (!IsOperationValid()) return;

    // centerV is the center vertex on which two neighboring 3-singularities are collapsed 
    auto& centerV = mesh.V.at(cid);
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
        
        centerV.N_Fids.push_back(faceToKeep.id); // add faces kept to the center vertex neighborhood
        
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
                edgeToKeep.N_Fids.at(1) == faceToRemove.id ? edgeToKeep.N_Fids.at(1) = faceToKeep.id : 1;
                continue;
            }
            
            // opposite vertex of the edge connecting to 3-singularity
            int idx2 = (idx+1)%edgeToKeep.Vids.size();
            auto& etk_v = mesh.V.at(edgeToKeep.Vids.at(idx2));
            
            // updates related to center vertex
            edgeToKeep.Vids.at(idx) = centerV.id;
            centerV.N_Vids.push_back(etk_v.id);
            centerV.N_Eids.push_back(edgeToKeep.id);
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
                edgeToKeep.N_Fids.at(0) == faceToKeep.id ? centerV.N_Fids.push_back(edgeToKeep.N_Fids.at(1)) : centerV.N_Fids.push_back(edgeToKeep.N_Fids.at(0)); 
                
                int nv_idx = std::distance(faceToKeep.Vids.begin(), std::find(faceToKeep.Vids.begin(), faceToKeep.Vids.end(), edgeToKeep.Vids.at(idx2)));
                if (s1.at(0) == faceToKeep.Vids.at((nv_idx+1)%faceToKeep.Vids.size()) || s1.at(1) == faceToKeep.Vids.at((nv_idx+1)%faceToKeep.Vids.size())) {
                    newVids.insert(newVids.begin(), faceToKeep.Vids.at(nv_idx));
                } else if (s1.at(0) == faceToKeep.Vids.at((nv_idx+3)%faceToKeep.Vids.size()) || s1.at(1) == faceToKeep.Vids.at((nv_idx+3)%faceToKeep.Vids.size())) {
                    newVids.push_back(faceToKeep.Vids.at(nv_idx));
                }
            } else {
                edgeToKeep.N_Fids.at(0) == faceToRemove.id ? centerV.N_Fids.push_back(edgeToKeep.N_Fids.at(1)) : centerV.N_Fids.push_back(edgeToKeep.N_Fids.at(0));
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

    
    // update n-singularities
    auto& s2_1 = mesh.V.at(s2.at(0));
    auto& s2_2 = mesh.V.at(s2.at(1));
    UpdateContents(s2_1.N_Vids, std::vector<size_t>{centerV.id});
    UpdateContents(s2_2.N_Vids, std::vector<size_t>{centerV.id});
    UpdateContents(s2_1.N_Eids, prev_e);
    UpdateContents(s2_2.N_Eids, prev_e);

    // collect vertices to smooth
    smoothV.push_back(centerV.id);
    smoothV.push_back(s2_1.id);
    smoothV.push_back(s2_2.id);
    smoothV.insert(smoothV.end(), s2_1.N_Vids.begin(), s2_1.N_Vids.end());
    smoothV.insert(smoothV.end(), s2_2.N_Vids.begin(), s2_2.N_Vids.end());
}