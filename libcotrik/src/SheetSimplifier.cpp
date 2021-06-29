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

void SheetSimplifier::GetChordCollapseOps(BaseComplexQuad& baseComplex, std::multiset<SimplificationOperation, bool(*)(SimplificationOperation, SimplificationOperation)>& SimplificationOps) {

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
            CollapseChordWithFeaturePreserved(key_edgeId, key_faceId, canceledFaceIds, canceledEdgeIds, SimplificationOps);
        }
        // break;
    }
} 

void SheetSimplifier::CollapseChordWithFeaturePreserved(std::unordered_map<size_t, size_t>& key_edgeId, std::unordered_map<std::string, size_t>& key_faceId,
    std::map<size_t, size_t>& canceledFaceIds, std::set<size_t>& canceledEdgeIds, std::multiset<SimplificationOperation, bool(*)(SimplificationOperation, SimplificationOperation)>& SimplificationOps) {
    
    SimplificationOperation Op;
    Op.type = "Chord_Collapse";
    
    for (auto edgeId : canceledEdgeIds) {
		auto& e = mesh.E[edgeId];
		auto& v0 = mesh.V[e.Vids[0]];
		auto& v1 = mesh.V[e.Vids[1]];
		auto key = (e.Vids[0] << 32) | e.Vids[1];
		if (key_edgeId.find(key) == key_edgeId.end()) std::cout << "Edge search Error !" << std::endl;
		auto centerVid = e.Vids.front();
		bool collapseToMidPoint = true;
		size_t featureType = 0;
		for (auto vid : e.Vids) {
			auto& v = mesh.V.at(vid);
			if (v.type != MAXID && v.type > featureType) {
				featureType = v.type;
				centerVid = vid;
			}
			if (v.isCorner) {
				v.type = CORNER;
				featureType = v.type;
				centerVid = vid;
			}
		}
        size_t source_vid = centerVid == v0.id ? v1.id : v0.id;
        size_t target_vid = centerVid == v0.id ? v0.id : v1.id;
        CollapseVertexToTarget(source_vid, target_vid, Op);
    }

    // Op.profitability /= mesh.totalArea;
    Op.profitability /= Op.n;
    // Op.profitability = 1;
    SimplificationOps.insert(Op);
}

void SheetSimplifier::ExtractAndCollapse(std::set<size_t>& canceledFids) {
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
    // std::vector<double>::iterator max_index = std::max_element(ranks.begin(), ranks.end());
    // double max_rank = (double) std::distance(ranks.begin(), max_index) + 1;
    for (int i = 0; i < ranks.size(); i++) {
        // std::vector<double>::iterator index = ranks.begin() + i;
        std::vector<double>::iterator index = std::max_element(ranks.begin(), ranks.end());
        // std::vector<double>::iterator index = std::min_element(ranks.begin(), ranks.end());
        sheetsPos.push_back((size_t) std::distance(ranks.begin(), index));
        // sheetsPos.push_back(i);
        *index = -1;
        // *index = max_rank;
		// ranks.erase(ranks.begin() + (size_t) std::distance(ranks.begin(), index));
		// i = 0;
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

void SheetSimplifier::GetSheetsProvisions(std::vector<std::set<size_t>>& canceledEdgesIds, std::vector<std::map<size_t, size_t>>& canceledFacesIds) {

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
            canceledEdgesIds.push_back(canceledEdgeIds);
            canceledFacesIds.push_back(canceledFaceIds);
        }
    }
}

void SheetSimplifier::CollapseSelectedSheets(std::set<size_t>& canceledFids, std::vector<std::set<size_t>>& canceledEdgesIds, std::vector<std::map<size_t, size_t>>& canceledFacesIds) {
    auto key_edgeId = get_key_edgeId(mesh);
    auto key_faceId = get_key_faceId(mesh);

    size_t n = canceledEdgesIds.size();
    for (int i = 0; i < n; i++) {
        std::set<size_t> canceledEdgeIds = canceledEdgesIds.at(i);
        std::map<size_t, size_t> canceledFaceIds = canceledFacesIds.at(i);
        // collapse_with_feature_preserved(key_edgeId, key_faceId, canceledFaceIds, canceledEdgeIds);
        CollapseWithFeaturePreserved(key_edgeId, key_faceId, canceledFaceIds, canceledEdgeIds);
        for (auto& item : canceledFaceIds)
            canceledFids.insert(item.first);
        // break;
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
