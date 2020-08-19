/*
 * SingleSheetSimplifier.cpp
 *
 *  Created on: Jan 20, 2019
 *      Author: cotrik
 */

#include "SingleSheetSimplifier.h"

SingleSheetSimplifier::SingleSheetSimplifier(Mesh& mesh) : SheetSimplifier(mesh) {
    // TODO Auto-generated constructor stub

}

SingleSheetSimplifier::~SingleSheetSimplifier() {
    // TODO Auto-generated destructor stub
}

void SingleSheetSimplifier::Run(std::set<size_t>& canceledFids) {
    auto key_edgeId = get_key_edgeId(mesh);
    auto key_faceId = get_key_faceId(mesh);
    BaseComplexQuad baseComplex(mesh);
    baseComplex.Build();

    BaseComplexSheetQuad baseComplexSheets(baseComplex);
    baseComplexSheets.Extract();
    for (const auto& e : mesh.E) {
    //for (int i = mesh.E.size(); --i >=0;) {
        //const auto& e = mesh.E.at(i);
        auto& v0 = mesh.V.at(e.Vids[0]);
        auto& v1 = mesh.V.at(e.Vids[1]);
        if (!v0.isCorner && !v1.isCorner && (v0.N_Fids.size() == 3 || v1.N_Fids.size() == 3)/*(v0.isSingularity || v1.isSingularity)*/) {
            std::map<size_t, size_t> canceledFaceIds;
            std::vector<size_t> linkEids(1, e.id);
            std::set<size_t> canceledEdgeIds = GetCanceledEdgeIds(linkEids, canceledFaceIds);
            if (!CanCollapseWithFeaturePreserved(baseComplexSheets, canceledFaceIds)) continue;
            if (!can_collapse_vids_with_feature_preserved(canceledEdgeIds)) continue;
            {
                std::cout << "collapse sheet edge " << e.id << "(" << v0.id << "," << v1.id << ")" << "\n";
                collapse_with_feature_preserved(key_edgeId, key_faceId, canceledFaceIds, canceledEdgeIds);
                for (auto& item : canceledFaceIds)
                    canceledFids.insert(item.first);
                return;
            }
        }
    }
}

void SingleSheetSimplifier::ExtractAndCollapse(std::set<size_t>& canceledFids) {
    std::vector<std::set<size_t>> canceledEdgesIds;
    std::vector<std::map<size_t, size_t>> canceledFacesIds;
    GetSheetsProvisions(canceledEdgesIds, canceledFacesIds);
    std::vector<std::set<size_t>> canceledVerticesIds;
    for (int i = 0; i < canceledEdgesIds.size(); i++) {
        std::set<size_t> canceledEdgeIds = canceledEdgesIds.at(i);
        std::set<size_t> canceledVIds;
        for (auto eid: canceledEdgeIds) {
            canceledVIds.insert(mesh.E.at(eid).Vids.at(0));
            canceledVIds.insert(mesh.E.at(eid).Vids.at(1));
        }
        canceledVerticesIds.push_back(canceledVIds);
    }
    std::map<size_t, double> sheetsRank;
    double aggRanks = 0;
    std::vector<double> ranks;
    for (int i = 0; i < canceledEdgesIds.size(); i++) {
        double agg_center = 0;
        double agg_length = 0;
        std::set<size_t> canceledEdgeIds = canceledEdgesIds.at(i);
        for (int j = 0; j < canceledEdgeIds.size(); j++) {
            int index1 = j;
            int index2 = (j + 1) % canceledEdgeIds.size();

            std::set<size_t>::iterator it1 = canceledEdgeIds.begin();
            std::advance(it1, index1);

            std::set<size_t>::iterator it2 = canceledEdgeIds.begin();
            std::advance(it2, index2);

            auto& e1 = mesh.E.at(*it1);
            auto& e2 = mesh.E.at(*it2);

            auto& v0 = mesh.V.at(e1.Vids.at(0));
            auto& v1 = mesh.V.at(e1.Vids.at(1));
            auto& v2 = mesh.V.at(e2.Vids.at(0));
            auto& v3 = mesh.V.at(e2.Vids.at(1));

            glm::dvec3 center1((v0.x + v1.x) / 2, (v0.y + v1.y) / 2, (v0.z + v1.z) / 2);
            glm::dvec3 center2((v2.x + v3.x) / 2, (v2.y + v3.y) / 2, (v2.z + v3.z) / 2);
            double t1 = glm::length(glm::dvec3(center1.x - center2.x, center1.y - center2.y, center1.z - center2.z));
            double t2 = glm::length(glm::dvec3(v0.x - v1.x, v0.y - v1.y, v0.z - v1.z));
            agg_center += t1;
            agg_length += t2;
        }
        ranks.push_back(agg_center * agg_length);
        aggRanks += (agg_center * agg_length);
    }

    std::vector<size_t> sheetsPos;
    for (int i = 0; i < ranks.size(); i++) {
        std::vector<double>::iterator max_index = std::max_element(ranks.begin(), ranks.end());
        sheetsPos.push_back((size_t) std::distance(ranks.begin(), max_index));
        *max_index = -1;
    }
    std::vector<std::set<size_t>> finalCanceledEdgesIds;
    std::vector<std::map<size_t, size_t>> finalCanceledFacesIds;
    for (auto sId: sheetsPos) {
        if (sId == -1) {
            continue;
        }
        finalCanceledEdgesIds.push_back(canceledEdgesIds.at(sId));
        finalCanceledFacesIds.push_back(canceledFacesIds.at(sId));
        std::set<size_t> canceledEdgeIds1 = canceledEdgesIds.at(sId);
        for (int j = 0; j < canceledEdgesIds.size(); j++) {
            if (j == sId) {
                continue;
            }
            std::set<size_t> canceledEdgeIds2 = canceledEdgesIds.at(j);
            bool disjoint = true;
            for (auto e1: canceledEdgeIds1) {
                if (canceledEdgeIds2.find(e1) != canceledEdgeIds2.end()) {
                    disjoint = false;
                    break;
                }
            }
            if (!disjoint) {
                auto it = std::find(sheetsPos.begin(), sheetsPos.end(), j);
                if (it != sheetsPos.end()) {
                    sheetsPos.at(std::distance(sheetsPos.begin(), it)) = -1;
                }
            }   
        }
        
        std::set<size_t> canceledVIds1 = canceledVerticesIds.at(sId);
        for (int j = 0; j < canceledVerticesIds.size(); j++) {
            if (j == sId) {
                continue;
            }
            std::set<size_t> canceledVIds2 = canceledVerticesIds.at(j);
            bool disjoint = true;
            for (auto v1: canceledVIds1) {
                if (canceledVIds2.find(v1) != canceledVIds2.end()) {
                    disjoint = false;
                    break;
                }
            }
            if (!disjoint) {
                auto it = std::find(sheetsPos.begin(), sheetsPos.end(), j);
                if (it != sheetsPos.end()) {
                    sheetsPos.at(std::distance(sheetsPos.begin(), it)) = -1;
                }
            }   
        }
    }
    CollapseSelectedSheets(canceledFids, finalCanceledEdgesIds, finalCanceledFacesIds);
}

void SingleSheetSimplifier::GetSheetsProvisions(std::vector<std::set<size_t>>& canceledEdgesIds, std::vector<std::map<size_t, size_t>>& canceledFacesIds) {
    BaseComplexQuad baseComplex(mesh);
    baseComplex.Build();

    BaseComplexSheetQuad baseComplexSheets(baseComplex);
    baseComplexSheets.Extract();

    for (const auto& e : mesh.E) {
        auto& v0 = mesh.V.at(e.Vids[0]);
        auto& v1 = mesh.V.at(e.Vids[1]);
        if (!v0.isCorner && !v1.isCorner && (v0.N_Fids.size() == 3 || v1.N_Fids.size() == 3)/*(v0.isSingularity || v1.isSingularity)*/) {
            std::map<size_t, size_t> canceledFaceIds;
            std::vector<size_t> linkEids(1, e.id);
            std::set<size_t> canceledEdgeIds = GetCanceledEdgeIds(linkEids, canceledFaceIds);
            if (!CanCollapseWithFeaturePreserved(baseComplexSheets, canceledFaceIds)) continue;
            if (!can_collapse_vids_with_feature_preserved(canceledEdgeIds)) continue;
            {
                canceledEdgesIds.push_back(canceledEdgeIds);
                canceledFacesIds.push_back(canceledFaceIds);
            }
        }
    }
}

void SingleSheetSimplifier::CollapseSelectedSheets(std::set<size_t>& canceledFids, std::vector<std::set<size_t>>& canceledEdgesIds, std::vector<std::map<size_t, size_t>>& canceledFacesIds) {
    auto key_edgeId = get_key_edgeId(mesh);
    auto key_faceId = get_key_faceId(mesh);

    size_t n = canceledEdgesIds.size();
    for (int i = 0; i < n; i++) {
        std::set<size_t> canceledEdgeIds = canceledEdgesIds.at(i);
        std::map<size_t, size_t> canceledFaceIds = canceledFacesIds.at(i);
        collapse_with_feature_preserved(key_edgeId, key_faceId, canceledFaceIds, canceledEdgeIds);
        for (auto& item : canceledFaceIds)
            canceledFids.insert(item.first);
    }
}

std::set<size_t> SingleSheetSimplifier::GetCanceledEdgeIds(const std::vector<size_t>& linkEids, std::map<size_t, size_t>& canceledFaceIds) {
    std::set<size_t> canceledEdgeIds;
    for (auto edgeId : linkEids) {
        canceledEdgeIds = GetAllParallelEdgeIds(edgeId);
        for (auto& eid : canceledEdgeIds)
            for (auto n_fid : mesh.E.at(eid).N_Fids)
                ++canceledFaceIds[n_fid];
        return canceledEdgeIds;
    }
    return canceledEdgeIds;
}

static bool find(const std::vector<size_t>& ids, size_t target_id) {
    bool res = false;
    for (auto id : ids) if (id == target_id) return true;
    return res;
}

bool SingleSheetSimplifier::CanCollapseWithFeaturePreserved(const BaseComplexSheetQuad& baseComplexSheets, std::map<size_t, size_t>& canceledFaceIds) {
    //const auto& baseComplexFIds = baseComplexSheets.sheets_componentFaceIds[sheetId];
    std::vector<size_t> baseComplexFIds;
    for (auto& item : canceledFaceIds) {
        auto fid = item.first;
        auto& f = mesh.F.at(fid);
        if (!find(baseComplexFIds, f.componentFid))
            baseComplexFIds.push_back(f.componentFid);
    }
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
