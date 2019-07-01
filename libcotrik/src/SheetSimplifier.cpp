/*
 * SheetSimplifier.cpp
 *
 *  Created on: Jan 20, 2019
 *      Author: cotrik
 */

#include "SheetSimplifier.h"

SheetSimplifier::SheetSimplifier(Mesh& mesh) : Simplifier(mesh) {
    // TODO Auto-generated constructor stub

}

SheetSimplifier::~SheetSimplifier() {
    // TODO Auto-generated destructor stub
}

void SheetSimplifier::Run(std::set<size_t>& canceledFids) {
    BaseComplexQuad baseComplex(mesh);
    baseComplex.Build();

    BaseComplexSheetQuad baseComplexSheets(baseComplex);
    baseComplexSheets.Extract();

    //auto dualMesh = Refine(mesh, 0);
    auto key_edgeId = get_key_edgeId(mesh);
    auto key_faceId = get_key_faceId(mesh);
    for (int sheetId = 0; sheetId < baseComplexSheets.sheets_componentEdgeIds.size(); ++sheetId) {
        auto& componentEdgeIds = baseComplexSheets.sheets_componentEdgeIds.at(sheetId);
        bool multiple_edges = false;
        for (auto component_eid : componentEdgeIds) {
            auto& component_e = baseComplex.componentE.at(component_eid);
            if (component_e.eids_link.size() > 1) {
                multiple_edges = true;
                break;
            }
        }
        if (multiple_edges) continue;
        bool has_interior_singularities = false;
        for (auto component_eid : componentEdgeIds) {
            auto& component_e = baseComplex.componentE.at(component_eid);
            for (auto vid : component_e.vids_link) {
                auto& v = mesh.V.at(vid);
                if (!v.isBoundary && v.N_Fids.size() != 4) {
                    has_interior_singularities = true;
                    break;
                }
            }
        }

        if (multiple_edges && !has_interior_singularities) continue;

        std::map<size_t, size_t> canceledFaceIds;
        std::set<size_t> canceledEdgeIds = GetCanceledEdgeIds(baseComplexSheets, canceledFaceIds, sheetId);

        if (!CanCollapseWithFeaturePreserved(baseComplexSheets, canceledFaceIds, sheetId)) continue;
        if (!can_collapse_vids_with_feature_preserved(canceledEdgeIds)) continue;
        {
            //std::string fname = std::string("before_global_collapsing.vtk");
            //MeshFileWriter writer(mesh, fname.c_str());
            //writer.WriteFile();
            std::cout << "collapse sheet " << sheetId << "\n";
            collapse_with_feature_preserved(key_edgeId, key_faceId, canceledFaceIds, canceledEdgeIds);
            for (auto& item : canceledFaceIds)
                canceledFids.insert(item.first);
            break;
        }
    }
}


std::set<size_t> SheetSimplifier::GetAllParallelEdgeIds(const size_t eid) {
    std::set<size_t> res;
    res.insert(eid);
    std::queue<size_t> q;
    q.push(eid);
    while (!q.empty()) {
        auto n = q.size();
        for (auto i = 0; i < n; ++i) {
            auto t = q.front();
            q.pop();
            auto& e = mesh.E.at(t);
            for (auto peid : e.parallelEids)
                if (res.find(peid) == res.end()) {
                    res.insert(peid);
                    q.push(peid);
                }
        }
    }
    return res;
}

std::set<size_t> SheetSimplifier::GetCanceledEdgeIds(const BaseComplexSheetQuad& baseComplexSheets, std::map<size_t, size_t>& canceledFaceIds, size_t sheetId) {
    std::set<size_t> canceledEdgeIds;
    for (auto baseComplexEdgeId : baseComplexSheets.sheets_componentEdgeIds[sheetId]) {
        for (auto edgeId : baseComplexSheets.baseComplex.componentE[baseComplexEdgeId].eids_link) {
            canceledEdgeIds = GetAllParallelEdgeIds(edgeId);
            for (auto& eid : canceledEdgeIds)
                for (auto n_fid : mesh.E.at(eid).N_Fids)
                    ++canceledFaceIds[n_fid];
            return canceledEdgeIds;
        }
    }
    return canceledEdgeIds;
}

bool SheetSimplifier::CanCollapseWithFeaturePreserved(const BaseComplexSheetQuad& baseComplexSheets, std::map<size_t, size_t>& canceledFaceIds, size_t sheetId) {
    const auto& baseComplexFIds = baseComplexSheets.sheets_componentFaceIds[sheetId];
    for (auto baseComplexFId : baseComplexFIds) {
        int count = 0;
        for (auto baseComplexVId : baseComplexSheets.baseComplex.componentF[baseComplexFId].Vids) {
            // if (baseComplexVId == MAXID) return false;
            auto vid = baseComplexSheets.baseComplex.Vids.at(baseComplexVId);
            auto& v = baseComplexSheets.baseComplex.mesh.V.at(vid);
            if (v.type == FEATURE || v.type == CORNER) ++count;
        }
        if (count > 3) return false;
    }

    std::set<size_t> baseComplexFaceIds_set(baseComplexFIds.begin(), baseComplexFIds.end());
    std::set<size_t> canceledEdgeIds;
    std::set<size_t> baseComplexVIds;
    for (auto baseComplexFId : baseComplexFIds)
        for (auto baseComplexVId : baseComplexSheets.baseComplex.componentF[baseComplexFId].Vids)
            baseComplexVIds.insert(baseComplexVId);

    for (auto baseComplexVId : baseComplexVIds) {
        int count = 0;
        const auto& bCFids = baseComplexSheets.baseComplex.componentV[baseComplexVId].N_Fids;
        for (auto bCFid : bCFids)
            if (baseComplexFaceIds_set.find(bCFid) != baseComplexFaceIds_set.end()) {
                ++count;
            }
        if (count >= 4)
            return false;
    }

    const auto& mesh = baseComplexSheets.baseComplex.mesh;
    for (auto& item : canceledFaceIds) {
        if (item.second >= 3) {
            auto& f = mesh.F.at(item.first);
            for (auto vid : f.Vids)
                if (vid >= mesh.V.size()) return false;
            size_t count = 0;
            for (auto vid : f.Vids) {
                auto& v = mesh.V.at(vid);
                if (v.type == CORNER || v.type == FEATURE) ++count;
            }
            if (count > 1) return false;
        }
    }
    for (auto edgeId : canceledEdgeIds) {
        auto& e = mesh.E[edgeId];
        for (auto vid : e.Vids)
            if (vid >= mesh.V.size()) return false;
        size_t count = 0;
        for (auto vid : e.Vids) {
            auto& v = mesh.V.at(vid);
            if (v.type == CORNER || v.type == FEATURE) ++count;
        }
        if (count > 1) return false;
    }

    return true;
}

void SheetSimplifier::CollapseWithFeaturePreserved(std::unordered_map<size_t, size_t>& key_edgeId, std::unordered_map<std::string, size_t>& key_faceId,
    std::map<size_t, size_t>& canceledFaceIds, std::set<size_t>& canceledEdgeIds) {
    for (auto& item : canceledFaceIds) {
        if (item.second >= 4) {
            auto& f = mesh.F.at(item.first);
            auto key = get_facekey(f);
            auto centerVid = f.Vids.front();
            size_t featureType = 0;
            for (auto vid : f.Vids) {
                auto& v = mesh.V.at(vid);
                if (v.type == FEATURE && featureType == 0) {
                    centerVid = vid;
                    featureType = FEATURE;
                }
                if (v.type == CORNER) centerVid = vid;
            }
            for (auto vid : f.Vids) {
                auto& v = mesh.V.at(vid);
                for (auto n_fid : v.N_Fids) {
                    auto& n_f = mesh.F.at(n_fid);
                    for (auto& n_vid : n_f.Vids)
                        if (n_vid == vid) n_vid = centerVid;
                }
            }
            for (auto eid : f.Eids)
                if (canceledEdgeIds.find(eid) != canceledEdgeIds.end()) canceledEdgeIds.erase(eid);
        }
    }
    for (auto edgeId : canceledEdgeIds) {
        auto& e = mesh.E[edgeId];
        auto& v0 = mesh.V[e.Vids[0]];
        auto& v1 = mesh.V[e.Vids[1]];
        auto key = (e.Vids[0] << 32) | e.Vids[1];
        if (key_edgeId.find(key) == key_edgeId.end()) std::cout << "Edge search Error !" << std::endl;
        //auto centerVid = mesh.V.size() + key_edgeId[key];
        auto centerVid = e.Vids.front();
        size_t featureType = 0;
        for (auto vid : e.Vids) {
            auto& v = mesh.V.at(vid);
            if (v.type != MAXID && v.type > featureType) {
                featureType = v.type;
                centerVid = vid;
            }
            if (v.isCorner) {
                // v.type = CORNER;
                featureType = v.type;
                centerVid = vid;
            }
        }
        for (auto vid : e.Vids) {
            auto& v = mesh.V.at(vid);
            for (auto n_fid : v.N_Fids) {
                auto& n_f = mesh.F.at(n_fid);
                for (auto& n_vid : n_f.Vids)
                    if (n_vid == vid) n_vid = centerVid;
            }
        }
    }
}
