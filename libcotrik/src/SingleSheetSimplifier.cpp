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
