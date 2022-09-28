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

void SheetSimplifier::GetChordCollapseOps(BaseComplexQuad& baseComplex, std::multiset<SimplificationOperationStruct, bool(*)(SimplificationOperationStruct, SimplificationOperationStruct)>& SimplificationOps) {

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
    std::map<size_t, size_t>& canceledFaceIds, std::set<size_t>& canceledEdgeIds, std::multiset<SimplificationOperationStruct, bool(*)(SimplificationOperationStruct, SimplificationOperationStruct)>& SimplificationOps) {
    
    SimplificationOperationStruct Op;
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
    /*int min_index = 0;
    double current_ranking = 1.0;
    for (int i = 0; i < canceledEdgesIds.size(); i++) {
        double r = CalculateRanking(canceledEdgesIds.at(i));
        if (r < current_ranking) {
            current_ranking = r;
            min_index = i;
        }
    }
    // std::cout << "min ranking: " << current_ranking << std::endl;
    if (current_ranking > 0.55) return;
    auto key_edgeId = get_key_edgeId(mesh);
    auto key_faceId = get_key_faceId(mesh);
    std::set<size_t> canceledEdgeIds = canceledEdgesIds.at(min_index);
    std::map<size_t, size_t> canceledFaceIds = canceledFacesIds.at(min_index);
    collapse_with_feature_preserved(key_edgeId, key_faceId, canceledFaceIds, canceledEdgeIds);
    for (auto& item : canceledFaceIds)
        canceledFids.insert(item.first);
    return;*/
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

    std::vector<std::set<size_t>> canceledChords;

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

        // bool skip = false;
        // for (auto eid: canceledEdgeIds) {
        //     auto& e = mesh.E.at(eid);
        //     // if (mesh.V.at(e.Vids.at(0)).isBoundary && !mesh.V.at(e.Vids.at(1)).isBoundary) skip = true;
        //     // if (!mesh.V.at(e.Vids.at(0)).isBoundary && mesh.V.at(e.Vids.at(1)).isBoundary) skip = true;
        //     // if (mesh.V.at(e.Vids.at(0)).type == FEATURE && !mesh.V.at(e.Vids.at(1)).type == FEATURE) skip = true;
        //     // if (!mesh.V.at(e.Vids.at(0)).type == FEATURE && mesh.V.at(e.Vids.at(1)).type == FEATURE) skip = true;
        //     if (mesh.V.at(e.Vids.at(0)).type == FEATURE || mesh.V.at(e.Vids.at(1)).type == FEATURE) skip = true;
        //     if (mesh.V.at(e.Vids.at(0)).isBoundary || mesh.V.at(e.Vids.at(1)).isBoundary) skip = true;
        // }
        // if (skip) continue;
        if (!CanCollapseWithFeaturePreserved(baseComplexSheets, canceledFaceIds, sheetId)) continue;
        if (!can_collapse_vids_with_feature_preserved(canceledEdgeIds)) continue;
        {
            canceledEdgesIds.push_back(canceledEdgeIds);
            canceledFacesIds.push_back(canceledFaceIds);
        }
    }
    /*std::cout << "writng " << canceledEdgesIds.size() << " chords" << std::endl;
    std::cout << "Writing output file" << std::endl;
    std::string outputf = "chord_collapse.vtk";
    std::ofstream ofs2(outputf.c_str());
    ofs2 << "# vtk DataFile Version 3.0\n"
        << outputf.c_str() << "\n"
        << "ASCII\n\n"
        << "DATASET UNSTRUCTURED_GRID\n";
    ofs2 << "POINTS " << mesh.V.size() << " double\n";
    std::vector<std::vector<size_t>> c_indices;
    int colorValue = 0;
    int ncolors = 15;
    std::vector<int> colors;

    for (auto c: canceledEdgesIds) {
        for (auto id: c) {
            c_indices.push_back(mesh.E.at(id).Vids);
            colors.push_back(colorValue%ncolors);
        }
        colorValue += 1;
    }
    // std::vector<size_t> c_indices = {12, 296};
    // std::cout << c_indices.size() << std::endl;
    for (size_t i = 0; i < mesh.V.size(); i++) {
        ofs2 << std::fixed << std::setprecision(7) <<  mesh.V.at(i).x << " " <<  mesh.V.at(i).y << " " <<  mesh.V.at(i).z << "\n";
    }
    ofs2 << "CELLS " << c_indices.size() << " " << 3 * c_indices.size() << std::endl;
    for (size_t i = 0; i < c_indices.size(); i++) {
        ofs2 << "2 " << c_indices.at(i)[0] << " " << c_indices.at(i)[1] << std::endl;
    }
    ofs2 << "CELL_TYPES " << c_indices.size() << "\n";
    for (size_t i = 0; i < c_indices.size(); i++) {
        ofs2 << "3" << std::endl;
    }

    ofs2 << "CELL_DATA " << c_indices.size() << "\n";
    ofs2 << "SCALARS fixed int\n";
    ofs2 << "LOOKUP_TABLE default\n";
    for (auto c: colors) {
        ofs2 << c << "\n";
    }*/
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

double SheetSimplifier::CalculateRanking(std::set<size_t>& canceledEdgeIds) {
    double alpha_q = 0.05;
    double alpha_d = 0.45;
    double alpha_v = 0.5;
    double E_q = 0.0;
    double E_d = 0.0;
    double E_v = 0.0;

    double it = 0.0;
    double max_qem = -1.0;
    double max_distance = -1.0;
    double max_val = -1.0;
    double avg_val = 0.0;
    for (auto eid: canceledEdgeIds) {
        Edge& e = mesh.E.at(eid);
        glm::dvec4 newV = CalculateQEM(e.Vids.at(0), e.Vids.at(1));
        if (newV[3] > max_qem) max_qem = newV[3];
        double d = CalculateAreaDistance(e.Vids.at(0), e.Vids.at(1));
        if (d > max_distance) max_distance = d;
        double val = CalculateValenceTerm(e.Vids.at(0), e.Vids.at(1));
        if (val > max_val) max_val = val;
        avg_val += val;
        it += 1.0;
    }
    E_q = max_qem;
    E_d = max_distance;
    E_v = max_val + (avg_val / it);
    double r = (alpha_q * (1 - exp(-E_q))) + (alpha_d * (1 - exp(-E_d))) + (alpha_v * (1 - exp(-E_v)));
    // std::cout << "ranking: " << r << std::endl;
    return r;
}

glm::dvec4 SheetSimplifier::CalculateQEM(size_t v1_id, size_t v2_id) {
    auto& v1 = mesh.V.at(v1_id);
    auto& v2 = mesh.V.at(v2_id);

    // std::cout << "Vertex v1: " << v1.id << " " << v1.x << " " << v1.y << " " << v1.z << std::endl;
    // std::cout << "Vertex v2: " << v2.id << " " << v2.x << " " << v2.y << " " << v2.z << std::endl;

    glm::dmat4 Q1(0.0);
    // std::cout << "Calculating Q1" << std::endl;
    for (auto fid: v1.N_Fids) {
        auto& f = mesh.F.at(fid);
        // std::cout << "neighbor face: " << f.id << " " << f.Vids.at(0) << " " << f.Vids.at(1) << " " << f.Vids.at(2) << " " << f.Vids.at(3) << std::endl;
        int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v1.id));
        auto& v_a = mesh.V.at(f.Vids.at((idx+1)%f.Vids.size()));
        auto& v_b = mesh.V.at(f.Vids.at((idx+3)%f.Vids.size()));

        // std::cout << "Vertex v_a: " << v_a.id << " " << v_a.x << " " << v_a.y << " " << v_a.z << std::endl;
        // std::cout << "Vertex v_b: " << v_b.id << " " << v_b.x << " " << v_b.y << " " << v_b.z << std::endl;
        
        glm::dvec3 A = v_a.xyz() - v1.xyz();
        glm::dvec3 B = v_b.xyz() - v1.xyz();
        // std::cout << "A: " << A.x << " " << A.y << " " << A.z << std::endl;
        // std::cout << "B: " << B.x << " " << B.y << " " << B.z << std::endl;
        
        glm::dvec3 n = glm::normalize(glm::cross(A, B));
        // std::cout << "normal: " << n.x << " " << n.y << " " << n.z << std::endl;

        glm::dvec4 p(n, glm::dot(n,v1.xyz()));
        // std::cout << "plane: " << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << std::endl;
        
        glm::dmat4 Q;
        Q[0] = glm::dvec4(p[0]*p[0], p[0]*p[1], p[0]*p[2], p[0]*p[3]);
        Q[1] = glm::dvec4(p[1]*p[0], p[1]*p[1], p[1]*p[2], p[1]*p[3]);
        Q[2] = glm::dvec4(p[2]*p[0], p[2]*p[1], p[2]*p[2], p[2]*p[3]);
        Q[3] = glm::dvec4(p[3]*p[0], p[3]*p[1], p[3]*p[2], p[3]*p[3]);

        // std::cout << "Q: " << std::endl;
        // std::cout << Q[0][0] << " " << Q[1][0] << " " << Q[2][0] << " " << Q[3][0] << std::endl;
        // std::cout << Q[0][1] << " " << Q[1][1] << " " << Q[2][1] << " " << Q[3][1] << std::endl;
        // std::cout << Q[0][2] << " " << Q[1][2] << " " << Q[2][2] << " " << Q[3][2] << std::endl;
        // std::cout << Q[0][3] << " " << Q[1][3] << " " << Q[2][3] << " " << Q[3][3] << std::endl;

        Q1 += Q;

        // std::cout << "Q1: " << std::endl;
        // std::cout << Q1[0][0] << " " << Q1[1][0] << " " << Q1[2][0] << " " << Q1[3][0] << std::endl;
        // std::cout << Q1[0][1] << " " << Q1[1][1] << " " << Q1[2][1] << " " << Q1[3][1] << std::endl;
        // std::cout << Q1[0][2] << " " << Q1[1][2] << " " << Q1[2][2] << " " << Q1[3][2] << std::endl;
        // std::cout << Q1[0][3] << " " << Q1[1][3] << " " << Q1[2][3] << " " << Q1[3][3] << std::endl;
        // std::cout << "***********************************" << std::endl;
    }
    // std::cout << "Q1: " << std::endl;
    // std::cout << Q1[0][0] << " " << Q1[1][0] << " " << Q1[2][0] << " " << Q1[3][0] << std::endl;
    // std::cout << Q1[0][1] << " " << Q1[1][1] << " " << Q1[2][1] << " " << Q1[3][1] << std::endl;
    // std::cout << Q1[0][2] << " " << Q1[1][2] << " " << Q1[2][2] << " " << Q1[3][2] << std::endl;
    // std::cout << Q1[0][3] << " " << Q1[1][3] << " " << Q1[2][3] << " " << Q1[3][3] << std::endl;


    glm::dmat4 Q2(0.0);
    // std::cout << "Calculating Q2" << std::endl;
    for (auto fid: v2.N_Fids) {
        auto& f = mesh.F.at(fid);
        // std::cout << "neighbor face: " << f.id << " " << f.Vids.at(0) << " " << f.Vids.at(1) << " " << f.Vids.at(2) << " " << f.Vids.at(3) << std::endl;
        int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v2.id));
        auto& v_a = mesh.V.at(f.Vids.at((idx+1)%f.Vids.size()));
        auto& v_b = mesh.V.at(f.Vids.at((idx+3)%f.Vids.size()));

        // std::cout << "Vertex v_a: " << v_a.id << " " << v_a.x << " " << v_a.y << " " << v_a.z << std::endl;
        // std::cout << "Vertex v_b: " << v_b.id << " " << v_b.x << " " << v_b.y << " " << v_b.z << std::endl;
        
        glm::dvec3 A = v_a.xyz() - v2.xyz();
        glm::dvec3 B = v_b.xyz() - v2.xyz();
        // std::cout << "A: " << A.x << " " << A.y << " " << A.z << std::endl;
        // std::cout << "B: " << B.x << " " << B.y << " " << B.z << std::endl;
        
        glm::dvec3 n = glm::normalize(glm::cross(A, B));
        // std::cout << "normal: " << n.x << " " << n.y << " " << n.z << std::endl;
        glm::dvec4 p(n, glm::dot(n,v2.xyz()));
        // std::cout << "plane: " << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << std::endl;
        
        glm::dmat4 Q;
        Q[0] = glm::dvec4(p[0]*p[0], p[0]*p[1], p[0]*p[2], p[0]*p[3]);
        Q[1] = glm::dvec4(p[1]*p[0], p[1]*p[1], p[1]*p[2], p[1]*p[3]);
        Q[2] = glm::dvec4(p[2]*p[0], p[2]*p[1], p[2]*p[2], p[2]*p[3]);
        Q[3] = glm::dvec4(p[3]*p[0], p[3]*p[1], p[3]*p[2], p[3]*p[3]);

        // std::cout << "Q: " << std::endl;
        // std::cout << Q[0][0] << " " << Q[1][0] << " " << Q[2][0] << " " << Q[3][0] << std::endl;
        // std::cout << Q[0][1] << " " << Q[1][1] << " " << Q[2][1] << " " << Q[3][1] << std::endl;
        // std::cout << Q[0][2] << " " << Q[1][2] << " " << Q[2][2] << " " << Q[3][2] << std::endl;
        // std::cout << Q[0][3] << " " << Q[1][3] << " " << Q[2][3] << " " << Q[3][3] << std::endl;

        Q2 += Q;
        // std::cout << "Q2: " << std::endl;
        // std::cout << Q2[0][0] << " " << Q2[1][0] << " " << Q2[2][0] << " " << Q2[3][0] << std::endl;
        // std::cout << Q2[0][1] << " " << Q2[1][1] << " " << Q2[2][1] << " " << Q2[3][1] << std::endl;
        // std::cout << Q2[0][2] << " " << Q2[1][2] << " " << Q2[2][2] << " " << Q2[3][2] << std::endl;
        // std::cout << Q2[0][3] << " " << Q2[1][3] << " " << Q2[2][3] << " " << Q2[3][3] << std::endl;
        // std::cout << "***********************************" << std::endl;
    }

    // std::cout << "Q2: " << std::endl;
    // std::cout << Q2[0][0] << " " << Q2[1][0] << " " << Q2[2][0] << " " << Q2[3][0] << std::endl;
    // std::cout << Q2[0][1] << " " << Q2[1][1] << " " << Q2[2][1] << " " << Q2[3][1] << std::endl;
    // std::cout << Q2[0][2] << " " << Q2[1][2] << " " << Q2[2][2] << " " << Q2[3][2] << std::endl;
    // std::cout << Q2[0][3] << " " << Q2[1][3] << " " << Q2[2][3] << " " << Q2[3][3] << std::endl;

    glm::dmat4 Q_ = Q1+Q2;
    // std::cout << "Q_: " << std::endl;
    // std::cout << Q_[0][0] << " " << Q_[1][0] << " " << Q_[2][0] << " " << Q_[3][0] << std::endl;
    // std::cout << Q_[0][1] << " " << Q_[1][1] << " " << Q_[2][1] << " " << Q_[3][1] << std::endl;
    // std::cout << Q_[0][2] << " " << Q_[1][2] << " " << Q_[2][2] << " " << Q_[3][2] << std::endl;
    // std::cout << Q_[0][3] << " " << Q_[1][3] << " " << Q_[2][3] << " " << Q_[3][3] << std::endl;

    glm::dmat4 Q_v = Q1+Q2;
    Q_v[0][3] = 0.0;
    Q_v[1][3] = 0.0;
    Q_v[2][3] = 0.0;
    Q_v[3][3] = 1.0;
    // std::cout << "Q_v: " << std::endl;
    // std::cout << Q_v[0][0] << " " << Q_v[1][0] << " " << Q_v[2][0] << " " << Q_v[3][0] << std::endl;
    // std::cout << Q_v[0][1] << " " << Q_v[1][1] << " " << Q_v[2][1] << " " << Q_v[3][1] << std::endl;
    // std::cout << Q_v[0][2] << " " << Q_v[1][2] << " " << Q_v[2][2] << " " << Q_v[3][2] << std::endl;
    // std::cout << Q_v[0][3] << " " << Q_v[1][3] << " " << Q_v[2][3] << " " << Q_v[3][3] << std::endl;

    glm::dvec4 v_in(0.0, 0.0, 0.0, 1.0);

    // std::cout << "determinant: " << glm::determinant(Q_v) << std::endl; 
    double determinant = glm::determinant(Q_v);
    glm::dvec4 v_(0.5 * (v1.xyz() + v2.xyz()), 1.0);
    if (fabs(determinant) > 0.0) {
        glm::dmat4 inv_Q_v = glm::inverse(Q_v);
        // std::cout << "inv_Q_v: " << std::endl;
        // std::cout << inv_Q_v[0][0] << " " << inv_Q_v[1][0] << " " << inv_Q_v[2][0] << " " << inv_Q_v[3][0] << std::endl;
        // std::cout << inv_Q_v[0][1] << " " << inv_Q_v[1][1] << " " << inv_Q_v[2][1] << " " << inv_Q_v[3][1] << std::endl;
        // std::cout << inv_Q_v[0][2] << " " << inv_Q_v[1][2] << " " << inv_Q_v[2][2] << " " << inv_Q_v[3][2] << std::endl;
        // std::cout << inv_Q_v[0][3] << " " << inv_Q_v[1][3] << " " << inv_Q_v[2][3] << " " << inv_Q_v[3][3] << std::endl;

        v_ = inv_Q_v*v_in;
    }

    glm::dvec4 v_int = Q_*v_;
    // std::cout << "v_: " << v_[0] << " " << v_[1] << " " << v_[2] << " " << v_[3] << std::endl;
    // std::cout << "v_int: " << v_int[0] << " " << v_int[1] << " " << v_int[2] << " " << v_int[3] << std::endl;

    v_[3] = glm::dot(v_, v_int);
    return v_;
}

double SheetSimplifier::CalculateAreaDistance(size_t v1_id, size_t v2_id) {
    return glm::distance(mesh.V.at(v1_id).xyz(), mesh.V.at(v2_id).xyz());
}

double SheetSimplifier::CalculateValenceTerm(size_t v1_id, size_t v2_id) {
    auto& v1 = mesh.V.at(v1_id);
    auto& v2 = mesh.V.at(v2_id);
    int idealVal = v1.isBoundary && v2.isBoundary ? 3 : 4;
    int target_v = v1.N_Vids.size() + v2.N_Vids.size() - idealVal;

    int val_v1 = v1.N_Vids.size() >= idealVal ? v1.N_Vids.size() - idealVal : idealVal - v1.N_Vids.size();
    int val_v2 = v2.N_Vids.size() >= idealVal ? v2.N_Vids.size() - idealVal : idealVal - v2.N_Vids.size();
    int target_valence = target_v >= idealVal ? target_v - idealVal : idealVal - target_v;
    double w_val = 0;
    w_val += target_valence > val_v1 ? target_valence - val_v1 : 0.0;
    w_val += target_valence > val_v2 ? target_valence - val_v2 : 0.0;

    return w_val;
}
