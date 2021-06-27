/*
* Simplifier.cpp
*
*  Created on: Dec 31, 2018
*      Author: cotrik
*/

#include "Simplifier.h"
#include <unordered_set>
#include <iostream>

const double PI = 3.1415926535;

std::vector<size_t> Simplifier::userCorners = {};
std::vector<size_t> Simplifier::canceledCorners = {};
double Simplifier::angle = 160;
int Simplifier::maxValence = 6;
int Simplifier::minValence = 3;
int Simplifier::iters = 10000;
int Simplifier::smoothIters = 20;
int Simplifier::resolution = 3;
bool Simplifier::featurePreserved = true;
bool Simplifier::COLLAPSE = true;
bool Simplifier::SPLIT = false;
bool Simplifier::CONFORMAL = true;
bool Simplifier::GLOBAL = true;
bool Simplifier::ROTATE = true;
bool Simplifier::COLLAPSE_DIAGNAL = true;
bool Simplifier::REMOVE_DOUBLET = true;
bool Simplifier::SHEET_SPLIT = true;
bool Simplifier::HALF = true;
bool Simplifier::TRIP = false;
bool Simplifier::checkCorner = true;
bool Simplifier::writeFile = false;

Simplifier::Simplifier(Mesh& mesh) : mesh(mesh) {

}

Simplifier::~Simplifier() {}

std::string Simplifier::get_facekey(const Face& f) {
	std::set<size_t> vids(f.Vids.begin(), f.Vids.end());
	std::string s;
	for (auto vid : vids)
		s += std::to_string(vid) + "@";
	return s;
}

std::set<size_t> Simplifier::get_allParallelEdgeIds(const size_t eid) {
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

std::map<size_t, size_t> Simplifier::get_canceledFaceIds(const std::set<size_t>& canceledEdgeIds) {
	std::map<size_t, size_t> canceledFaceIds;
	for (auto& eid : canceledEdgeIds)
		for (auto n_fid : mesh.E.at(eid).N_Fids)
			++canceledFaceIds[n_fid];
	return canceledFaceIds;
}

std::set<size_t> Simplifier::get_canceledEdgeIds(const BaseComplexSheetQuad& baseComplexSheets, std::map<size_t, size_t>& canceledFaceIds, size_t sheetId) {
	std::set<size_t> canceledEdgeIds;
	for (auto baseComplexEdgeId : baseComplexSheets.sheets_componentEdgeIds[sheetId]) {
		for (auto edgeId : baseComplexSheets.baseComplex.componentE[baseComplexEdgeId].eids_link) {
			//auto& e = baseComplexSheets.baseComplex.mesh.E[edgeId];
			//for (auto n_fid : e.N_Fids)
			//	++canceledFaceIds[n_fid];
			//canceledEdgeIds.insert(edgeId);
			canceledEdgeIds = get_allParallelEdgeIds(edgeId);
			for (auto& eid : canceledEdgeIds)
				for (auto n_fid : mesh.E.at(eid).N_Fids)
					++canceledFaceIds[n_fid];
			return canceledEdgeIds;
		}
	}
	return canceledEdgeIds;
}

bool Simplifier::can_collapse(const BaseComplexSheetQuad& baseComplexSheets, std::map<size_t, size_t>& canceledFaceIds, size_t sheetId) {
	const auto& baseComplexFIds = baseComplexSheets.sheets_componentFaceIds[sheetId];
	std::set<size_t> baseComplexFaceIds_set(baseComplexFIds.begin(), baseComplexFIds.end());
	std::set<size_t> canceledEdgeIds;
	std::set<size_t> baseComplexVIds;
	for (auto baseComplexFId : baseComplexFIds)
		for (auto baseComplexVId : baseComplexSheets.baseComplex.componentF[baseComplexFId].Vids)
			baseComplexVIds.insert(baseComplexVId);

	for (auto baseComplexVId : baseComplexVIds) {
		bool all_in = true;
		const auto& bCFids = baseComplexSheets.baseComplex.componentV[baseComplexVId].N_Fids;
		for (auto bCFid : bCFids)
			if (baseComplexFaceIds_set.find(bCFid) == baseComplexFaceIds_set.end()) {
				all_in = false;
				break;
			}
		if (all_in)
			return false;
	}

	const auto& mesh = baseComplexSheets.baseComplex.mesh;
	for (auto& item : canceledFaceIds) {
		if (item.second >= 4) {
			auto& f = mesh.F.at(item.first);
			for (auto vid : f.Vids)
				if (vid >= mesh.V.size()) return false;
		}
	}
	for (auto edgeId : canceledEdgeIds) {
		auto& e = mesh.E[edgeId];
		for (auto vid : e.Vids)
			if (vid >= mesh.V.size()) return false;
	}

	return true;
}

void Simplifier::collapse(std::unordered_map<size_t, size_t>& key_edgeId, std::unordered_map<std::string, size_t>& key_faceId,
	std::map<size_t, size_t>& canceledFaceIds, std::set<size_t>& canceledEdgeIds) {
	for (auto& item : canceledFaceIds) {
		if (item.second >= 4) {
			auto& f = mesh.F.at(item.first);
			auto key = get_facekey(f);
			auto centerVid = mesh.V.size() + mesh.E.size() + key_faceId[key];
			for (auto vid : f.Vids) {
				// if (vid >= mesh.V.size()) continue;
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
		auto centerVid = mesh.V.size() + key_edgeId[key];
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

bool Simplifier::can_collapse_with_feature_preserved(const BaseComplexSheetQuad& baseComplexSheets, std::map<size_t, size_t>& canceledFaceIds, size_t sheetId) {
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

void Simplifier::collapse_with_feature_preserved(std::unordered_map<size_t, size_t>& key_edgeId, std::unordered_map<std::string, size_t>& key_faceId,
	std::map<size_t, size_t>& canceledFaceIds, std::set<size_t>& canceledEdgeIds) {
	/*for (auto& item : canceledFaceIds) {
		if (item.second >= 4) {
			// continue;
			std::cout << "collapsing from canceled face ids" << std::endl;
			// std::cout << "face id: " << item.first << " #second: " << item.second << std::endl; 
			auto& f = mesh.F.at(item.first);
			auto key = get_facekey(f);
			auto centerVid = f.Vids.front();
			//auto centerVid = mesh.V.size() + mesh.E.size() + key_faceId[key];
			//Vertex centerV = 0.25*(mesh.V[f.Vids[0]] + mesh.V[f.Vids[1]] + mesh.V[f.Vids[2]] + mesh.V[f.Vids[3]]);
			//centerV.id = mesh.V.size();
			//centerV.type = mesh.V.at(f.Vids[0]).type;
			//centerV.label = mesh.V.at(f.Vids[0]).label;
			//centerV.labels = mesh.V.at(f.Vids[0]).labels;
			//centerV.patch_id = mesh.V.at(f.Vids[0]).patch_id;
			//centerV.patch_ids = mesh.V.at(f.Vids[0]).patch_ids;
			size_t featureType = 0;
			for (auto vid : f.Vids) {
				auto& v = mesh.V.at(vid);
				if (v.type == FEATURE && featureType == 0) {
					centerVid = vid;
				}
				if (v.type == CORNER) {
					centerVid = vid;
				}
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
	}*/
	for (auto edgeId : canceledEdgeIds) {
		auto& e = mesh.E[edgeId];
		auto& v0 = mesh.V[e.Vids[0]];
		auto& v1 = mesh.V[e.Vids[1]];
		// std::cout << "edge " << e.id << ": " << v0.id << " " << v1.id << std::endl; 
		auto key = (e.Vids[0] << 32) | e.Vids[1];
		if (key_edgeId.find(key) == key_edgeId.end()) std::cout << "Edge search Error !" << std::endl;
		//auto centerVid = mesh.V.size() + key_edgeId[key];
		auto centerVid = e.Vids.front();
		bool collapseToMidPoint = true;
		size_t featureType = 0;
		for (auto vid : e.Vids) {
			auto& v = mesh.V.at(vid);
			if (v.type != MAXID && v.type > featureType) {
				featureType = v.type;
				centerVid = vid;
				collapseToMidPoint = false;
			}
			if (v.isCorner) {
				v.type = CORNER;
				featureType = v.type;
				centerVid = vid;
				collapseToMidPoint = false;
			}
		}
		for (auto vid : e.Vids) {
			auto& v = mesh.V.at(vid);
			
			// if (vid != centerVid) {
				// auto& centerV = mesh.V.at(centerVid);
				// std::set<size_t> newE;
				// std::set<size_t> newF;
				// newE.insert(centerV.N_Eids.begin(), centerV.N_Eids.end());
				// newE.insert(v.N_Eids.begin(), v.N_Eids.end());
				// newF.insert(centerV.N_Fids.begin(), centerV.N_Fids.end());
				// newF.insert(v.N_Fids.begin(), v.N_Fids.end());

				// centerV.N_Eids.clear();
				// centerV.N_Fids.clear();
				// centerV.N_Eids.insert(centerV.N_Eids.begin(), newE.begin(), newE.end());
				// centerV.N_Fids.insert(centerV.N_Fids.begin(), newF.begin(), newF.end());
			// }

			for (auto n_fid : v.N_Fids) {
				auto& n_f = mesh.F.at(n_fid);
				for (auto& n_vid : n_f.Vids)
					if (n_vid == vid) n_vid = centerVid;
			}
			// for (auto n_eid: v.N_Eids) {
			// 	auto& n_e = mesh.E.at(n_eid);
			// 	for (auto& n_vid: n_e.Vids)
			// 		if (n_vid == vid) n_vid = centerVid;
			// }
			
		}
		if (collapseToMidPoint) {
			mesh.V.at(centerVid).x = (v1.x + v0.x) / 2;
			mesh.V.at(centerVid).y = (v1.y + v0.y) / 2;
			mesh.V.at(centerVid).z = (v1.z + v0.z) / 2;
		}
	}
}


bool Simplifier::can_collapse_vids_with_feature_preserved(const std::vector<size_t>& vids) {
	const auto& v0 = mesh.V.at(vids[0]);
	const auto& v1 = mesh.V.at(vids[1]);

	int count = 0;
	if (v0.isCorner) ++count;
	if (v1.isCorner) ++count;
	if (count > 1) return false;

	std::set<size_t> labels;
	if (!v0.isCorner && v0.label != MAXID) labels.insert(v0.label);
	if (!v1.isCorner && v1.label != MAXID) labels.insert(v1.label);



	//if (v0.N_Fids.size() + v1.N_Fids.size() - 4 > Simplifier::maxValence) return false;
	//if (v0.N_Fids.size() + v1.N_Fids.size() - 4 < Simplifier::minValence) return false;
	if (!v0.isBoundary && !v1.isBoundary && v0.N_Fids.size() + v1.N_Fids.size() - 4 > Simplifier::maxValence)
	    return false;
	if (!v0.isBoundary && !v1.isBoundary && v0.N_Fids.size() + v1.N_Fids.size() - 4 < 2/*Simplifier::minValence*/)
	    return false;

//    if (v0.isBoundary && v1.isBoundary && v0.N_Fids.size() + v1.N_Fids.size() - 2 > Simplifier::maxValence)
//        return false;
//    if (v0.isBoundary && v1.isBoundary && v0.N_Fids.size() + v1.N_Fids.size() - 2 < 2/*Simplifier::minValence*/)
//        return false;
    if (v0.type == FEATURE && v1.type == FEATURE && v0.N_Fids.size() + v1.N_Fids.size() != 4/*Simplifier::minValence*/)
        return false;
    if (v0.type == CORNER && v1.type == FEATURE && v0.N_Fids.size() + v1.N_Fids.size() != 2 + v0.idealValence/*Simplifier::minValence*/)
        return false;
    if (v1.type == CORNER && v0.type == FEATURE && v0.N_Fids.size() + v1.N_Fids.size() != 2 + v1.idealValence/*Simplifier::minValence*/)
        return false;

	if (v0.isBoundary && !v1.isBoundary && v0.N_Fids.size() + v1.N_Fids.size() - 4 > Simplifier::maxValence)
	    return false;
	if (!v0.isBoundary && v1.isBoundary && v0.N_Fids.size() + v1.N_Fids.size() - 4 > Simplifier::maxValence)
	    return false;
	if (v0.isBoundary && !v1.isBoundary && v0.N_Fids.size() + v1.N_Fids.size() - 4 < Simplifier::minValence)
	    return false;
	if (!v0.isBoundary && v1.isBoundary && v0.N_Fids.size() + v1.N_Fids.size() - 4 < Simplifier::minValence)
	    return false;

	if ((v0.idealValence >= 3 || v1.idealValence >= 3) && v0.isBoundary && v1.isBoundary &&
	        v0.N_Fids.size() + v1.N_Fids.size() - 2 < Simplifier::minValence)
	    return false;



	//if (labels.size() == 1 && v0.type == CORNER && v1.type == FEATURE) return false;
	//if (labels.size() == 1 && v0.type == FEATURE && v1.type == CORNER) return false;
    if (labels.size() == 1 && v0.type == CORNER && v1.type == FEATURE && v0.labels.find(v1.label) == v0.labels.end())
        return false;
    if (labels.size() == 1 && v0.type == FEATURE && v1.type == CORNER && v1.labels.find(v0.label) == v1.labels.end())
        return false;
	if (labels.size() > 2)
	    return false;
	if (labels.size() == 2 && (v0.isCorner || v1.isCorner))
	    return false;



	//if ((v0.type == FEATURE || v1.type == FEATURE) && v0.N_Fids.size() + v1.N_Fids.size() - 4 <= Simplifier::minValence) return false;
//	if (!v0.isBoundary && !v1.isBoundary && (v0.type == FEATURE || v1.type == FEATURE)
//		&& v0.N_Fids.size() + v1.N_Fids.size() - 4 <= Simplifier::minValence) return false;
//	if (!v0.isBoundary && v1.isBoundary && (v0.type == FEATURE || v1.type == FEATURE)
//		&& v0.N_Fids.size() + v1.N_Fids.size() - 2 <= Simplifier::minValence) return false;
//	if (v0.isBoundary && !v1.isBoundary && (v0.type == FEATURE || v1.type == FEATURE)
//		&& v0.N_Fids.size() + v1.N_Fids.size() - 2 <= Simplifier::minValence) return false;
	return true;
}

bool Simplifier::can_collapse_vids_with_feature_preserved(const std::set<size_t>& eids) {
	for (auto eid : eids) {
		auto& e = mesh.E.at(eid);
		if (!can_collapse_vids_with_feature_preserved(e.Vids)) return false;
	}
	return true;
}

bool Simplifier::collapse_sheet(std::set<size_t>& canceledFids, size_t parallel_edge_id) {
	auto canceledEdgeIds = get_allParallelEdgeIds(parallel_edge_id);
	std::map<size_t, size_t> canceledFaceIds = get_canceledFaceIds(canceledEdgeIds);
	if (can_collapse_vids_with_feature_preserved(canceledEdgeIds)) {
		auto key_edgeId = get_key_edgeId(mesh);
		auto key_faceId = get_key_faceId(mesh);
		collapse_with_feature_preserved(key_edgeId, key_faceId, canceledFaceIds, canceledEdgeIds);
		for (auto& item : canceledFaceIds)
			canceledFids.insert(item.first);
		return true;
	}
	return false;
}

void Simplifier::global_simplify(std::set<size_t>& canceledFids) {
	//{
	//	BaseComplexQuad baseComplex(mesh);
	//	baseComplex.Build();

	//	update(baseComplex);
	//	init();
	//}
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
		// if (multiple_edges) continue;
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
		std::set<size_t> canceledEdgeIds = get_canceledEdgeIds(baseComplexSheets, canceledFaceIds, sheetId);

		if (!can_collapse_with_feature_preserved(baseComplexSheets, canceledFaceIds, sheetId)) continue;
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

void Simplifier::global_simplify1(std::set<size_t>& canceledFids) {
	{
		BaseComplexQuad baseComplex(mesh);
		baseComplex.Build();

		update(baseComplex);
		init();
	}
	global_simplify(canceledFids);
}

std::unordered_map<size_t, size_t> Simplifier::get_key_edgeId(const Mesh& mesh) {
	std::unordered_map<size_t, size_t> key_edgeId;
	for (size_t i = 0; i < mesh.E.size(); ++i) {
		const auto& e = mesh.E.at(i);
		key_edgeId[(e.Vids[0] << 32) | e.Vids[1]] = i;
		key_edgeId[(e.Vids[1] << 32) | e.Vids[0]] = i;
	}
	return key_edgeId;
}

std::unordered_map<std::string, size_t> Simplifier::get_key_faceId(const Mesh& mesh) {
	std::unordered_map<std::string, size_t> key_faceId;
	for (size_t i = 0; i < mesh.F.size(); ++i) {
		const auto& f = mesh.F.at(i);
		std::set<size_t> vids(f.Vids.begin(), f.Vids.end());
		std::string s;
		for (auto vid : vids)
			s += std::to_string(vid) + "@";
		key_faceId[s] = i;
	}
	return key_faceId;
}

const int QuadRefine[4][4] = {
	0, 4, 8, 7,
	1, 5, 8, 4,
	2, 6, 8, 5,
	3, 7, 8, 6
};

const int TriRefine[3][4] = {
	0, 3, 6, 5,
	1, 4, 6, 3,
	2, 5, 6, 4
};

Mesh Simplifier::Refine(const Mesh& hex_mesh, int clockwise) {
	const Mesh& new_mesh = hex_mesh;
	////////////////////////////////////////////////////////////////////////////
	// add vertices
	std::vector<Vertex> new_vertex(new_mesh.V.size() + new_mesh.E.size() + new_mesh.C.size());
	for (size_t i = 0; i < new_mesh.V.size(); i++)
		new_vertex.at(i) = new_mesh.V.at(i);
	size_t offset = new_mesh.V.size();
	for (size_t i = 0; i < new_mesh.E.size(); i++) {
		const Edge& e = new_mesh.E.at(i);
		const Vertex& v0 = new_mesh.V.at(e.Vids[0]);
		const Vertex& v1 = new_mesh.V.at(e.Vids[1]);
		new_vertex.at(offset + i) = 0.5 * (v0.xyz() + v1.xyz());
	}
	offset = new_mesh.V.size() + new_mesh.E.size();
	size_t numOfTri = 0, numOfQuad = 0;
	for (size_t i = 0; i < new_mesh.C.size(); i++) {
		const auto& f = new_mesh.C.at(i);
		if (f.Vids.size() == 4) {
			const Vertex& v0 = new_mesh.V.at(f.Vids[0]);
			const Vertex& v1 = new_mesh.V.at(f.Vids[1]);
			const Vertex& v2 = new_mesh.V.at(f.Vids[2]);
			const Vertex& v3 = new_mesh.V.at(f.Vids[3]);
			new_vertex.at(offset + i) = 0.25 * (v0.xyz() + v1.xyz() + v2.xyz() + v3.xyz());
			++numOfQuad;
		} else  if (f.Vids.size() == 3) {
			const Vertex& v0 = new_mesh.V.at(f.Vids[0]);
			const Vertex& v1 = new_mesh.V.at(f.Vids[1]);
			const Vertex& v2 = new_mesh.V.at(f.Vids[2]);
			new_vertex.at(offset + i) = 0.3333333 * (v0.xyz() + v1.xyz() + v2.xyz());
			++numOfTri;
		}
	}
	auto key_edgeId = get_key_edgeId(new_mesh);
	//auto key_faceId = get_key_faceId(new_mesh);
	Cell cell(4);
	std::vector<Cell> new_cells(numOfTri * 3 + numOfQuad * 4, cell);
	int count = 0;
	for (size_t i = 0; i < new_mesh.C.size(); i++) {
		unsigned long v_index[9];
		const auto & f = new_mesh.C.at(i);
		for (auto j = 0; j < f.Vids.size(); j++)
			v_index[j] = f.Vids.at(j);
		//        if (clockwise != 0) {
		//            std::swap(v_index[1], v_index[3]);
		//            std::swap(v_index[5], v_index[7]);
		//        }
		if (f.Vids.size() == 4) {
			for (unsigned long j = 0; j < 4; j++) {
				const Edge e({ f.Vids.at(QuadEdge[j][0]), f.Vids.at(QuadEdge[j][1]) });
				unsigned long long key = (e.Vids[0] << 32) | e.Vids[1];
				if (key_edgeId.find(key) == key_edgeId.end()) std::cout << "Edge search Error !" << std::endl;
				auto e_index = key_edgeId[key];
				v_index[4 + j] = new_mesh.V.size() + e_index;
			}
			v_index[8] = new_mesh.V.size() + new_mesh.E.size() + i;
			for (int k = 0; k < 4; k++, count++)
				for (int j = 0; j < 4; j++)
					new_cells[count].Vids[j] = v_index[QuadRefine[k][j]];
		} else if (f.Vids.size() == 3) {
			for (unsigned long j = 0; j < 3; j++) {
				const Edge e({ f.Vids.at(TriEdge[j][0]), f.Vids.at(TriEdge[j][1]) });
				unsigned long long key = (e.Vids[0] << 32) | e.Vids[1];
				if (key_edgeId.find(key) == key_edgeId.end()) std::cout << "Edge search Error !" << std::endl;
				auto e_index = key_edgeId[key];
				v_index[3 + j] = new_mesh.V.size() + e_index;
			}
			v_index[6] = new_mesh.V.size() + new_mesh.E.size() + i;
			for (int k = 0; k < 3; k++, count++)
				for (int j = 0; j < 4; j++)
					new_cells[count].Vids[j] = v_index[TriRefine[k][j]];
		}
	}
	Mesh mesh(new_vertex, new_cells, QUAD);
	return mesh;
}

////////////////////////////////////////////////////////////////////////
void Simplifier::get_parallel_edgeids(size_t start_edge_id, size_t start_face_id,
	std::set<size_t>& parallel_edgeids, std::set<size_t>& parallel_faceids) {
	if (parallel_edgeids.find(start_edge_id) != parallel_edgeids.end()) return;
	parallel_edgeids.insert(start_edge_id);
	parallel_faceids.insert(start_face_id);
	size_t next_edge_id;
	auto& start_edge = mesh.E.at(start_edge_id);
	auto& v0 = mesh.V.at(start_edge.Vids[0]);
	auto& v1 = mesh.V.at(start_edge.Vids[1]);
	for (auto edgeid : mesh.F.at(start_face_id).Eids) {
		if (edgeid == start_edge_id) continue;
		auto& e = mesh.E.at(edgeid);
		auto& v_0 = mesh.V.at(e.Vids[0]);
		auto& v_1 = mesh.V.at(e.Vids[1]);
		if (v_0.id != v0.id && v_0.id != v1.id && v_1.id != v0.id && v_1.id != v1.id) {
			next_edge_id = edgeid;
			break;
		}
	}
	auto& next_edge = mesh.E.at(next_edge_id);
	size_t next_face_id = next_edge.N_Fids[0] == start_face_id ? next_edge.N_Fids[1] : next_edge.N_Fids[0];
	get_parallel_edgeids(next_edge_id, next_face_id, parallel_edgeids, parallel_faceids);
}

size_t Simplifier::get_faceid(size_t vid, size_t exclude_vid) {
	auto& v = mesh.V.at(vid);
	size_t res;
	for (auto fid : v.N_Fids) {
		auto& f = mesh.F.at(fid);
		bool found_exclude_vid = false;
		for (auto fvid : f.Vids)
			if (fvid == exclude_vid) {
				found_exclude_vid = true;
				break;
			}
		if (!found_exclude_vid) {
			res = fid;
			break;
		}
	}
	return res;
}

size_t Simplifier::get_diagnal_vid(size_t vid, size_t fid) {
	auto& f = mesh.F.at(fid);
	for (int i = 0; i < 4; ++i) {
		if (f.Vids[i] == vid) return f.Vids.at((i + 2) % 4);
	}
	std::cerr << "ERROR get_diagnal_vid\n";
	return MAXID;
}

size_t Simplifier::get_diagnal_vid(size_t vid, const std::vector<size_t>& fids) {
	std::set<size_t> fvids;
	for (auto fid : fids) {
		auto& f = mesh.F.at(fid);
		fvids.insert(f.Vids.begin(), f.Vids.end());
	}
	auto& v = mesh.V.at(vid);
	for (auto fid : v.N_Fids) {
		if (fid == fids[0] || fid == fids[1]) continue;
		auto& f = mesh.F.at(fid);
		bool found = false;
		for (auto fvid : f.Vids) {
			if (fvid == vid) continue;
			if (fvids.find(fvid) != fvids.end()) {
				found = true;
				break;
			}
		}
		if (!found) return get_diagnal_vid(vid, fid);
	}
	std::cerr << "ERROR get_diagnal_vid\n";
	return MAXID;
}

std::vector<size_t> Simplifier::get_collapse_vids(size_t vid, size_t eid) {
	std::vector<size_t> p;
	auto& v = mesh.V.at(vid);
	auto& e = mesh.E.at(eid);
	for (auto fid : e.N_Fids) {
		auto& f = mesh.F.at(fid);
		for (auto feid : f.Eids) {
			auto& fe = mesh.E.at(feid);
			std::set<int> vids(fe.Vids.begin(), fe.Vids.end());
			vids.insert(e.Vids.begin(), e.Vids.end());
			if (vids.size() == 4) {
				for (auto nvid : v.N_Vids) {
					if (nvid == fe.Vids[0] || nvid == fe.Vids[1]) {
						p.push_back(nvid);
						break;
					}
				}
				break;
			}
		}
	}
	return p;
}

std::vector<size_t> Simplifier::get_split_vids(size_t vid, size_t eid) {
	std::vector<size_t> p;
	auto& v = mesh.V.at(vid);
	auto& e = mesh.E.at(eid);
	for (auto peid : e.parallelEids) {
		auto& pe = mesh.E.at(peid);
		for (auto pevid : pe.Vids) {
			auto& pev = mesh.V.at(pevid);
			for (auto nvid : pev.N_Vids)
				if (nvid == vid) {
					p.push_back(pevid);
					break;
				}
		}
	}
	return p;
}

bool Simplifier::can_collapse_vids(const std::vector<size_t>& vids, size_t target_vid) {
	return ((mesh.V.at(vids[0]).N_Fids.size() + mesh.V.at(vids[1]).N_Fids.size() < 10) &&
		(mesh.V.at(vids[0]).N_Fids.size() + mesh.V.at(vids[1]).N_Fids.size() > 6));
}

void Simplifier::collapse_vids(const std::vector<size_t>& vids, size_t target_vid) {
	for (auto vid : vids) {
		auto& v = mesh.V.at(vid);
		for (auto fid : v.N_Fids) {
			auto& f = mesh.F.at(fid);
			for (auto& fvid : f.Vids) {
				if (fvid == vid) {
					fvid = target_vid;
					break;
				}
			}
		}
	}
}

bool Simplifier::can_collapse(const std::vector<size_t>& linkVids, const std::vector<size_t>& linkEids) {
	std::set<size_t> vids;
	for (size_t i = 0; i < linkEids.size(); ++i) {
		auto vid = linkVids[i];
		auto eid = linkEids[i];
		auto p = get_collapse_vids(vid, eid);
		if (!can_collapse_vids(p, vid)) return false;
		vids.insert(p.begin(), p.end());
	}
	auto vid = linkVids.back();
	auto eid = linkEids.back();
	auto p = get_collapse_vids(vid, eid);
	if (!can_collapse_vids(p, vid)) return false;
	vids.insert(p.begin(), p.end());

	if (vids.size() < 2 * linkVids.size()) return false; // tangent;
	return true;
}

void Simplifier::collapse(const std::vector<size_t>& linkVids, const std::vector<size_t>& linkEids) {
	for (size_t i = 0; i < linkEids.size(); ++i) {
		auto vid = linkVids[i];
		auto eid = linkEids[i];
		auto p = get_collapse_vids(vid, eid);
		collapse_vids(p, vid);
	}
	auto vid = linkVids.back();
	auto eid = linkEids.back();
	auto p = get_collapse_vids(vid, eid);
	collapse_vids(p, vid);
}

bool Simplifier::can_collapse_vids_with_feature_preserved(const std::vector<size_t>& vids, size_t target_vid) {
	if (vids.empty()) {
		MeshFileWriter writer(mesh, "err.vtk");
		writer.WriteFile();
	}
	const auto& v0 = mesh.V.at(vids[0]);
	const auto& v1 = mesh.V.at(vids[1]);
	const auto& v = mesh.V.at(target_vid);

	int count = 0;
	if (v0.isCorner) ++count;
	if (v1.isCorner) ++count;
	if (v.isCorner) ++count;

	std::set<size_t> labels;
	if (!v0.isCorner && v0.label != MAXID) labels.insert(v0.label);
	if (!v1.isCorner && v1.label != MAXID) labels.insert(v1.label);
	if (!v.isCorner && v.label != MAXID) labels.insert(v.label);
//	if (v0.N_Fids.size() + v1.N_Fids.size() - 4 > Simplifier::maxValence) return false;
//	if (v0.N_Fids.size() + v1.N_Fids.size() - 4 < Simplifier::minValence) return false;
    if (!v0.isBoundary && !v1.isBoundary && v0.N_Fids.size() + v1.N_Fids.size() - 4 > Simplifier::maxValence)
        return false;
    if (!v0.isBoundary && !v1.isBoundary && v0.N_Fids.size() + v1.N_Fids.size() - 4 < 2/*Simplifier::minValence*/)
        return false;

    if (v0.isBoundary && v1.isBoundary && v0.N_Fids.size() + v1.N_Fids.size() - 2 > Simplifier::maxValence)
        return false;
    if (v0.isBoundary && v1.isBoundary && v0.N_Fids.size() + v1.N_Fids.size() - 2 < 2/*Simplifier::minValence*/)
        return false;

    if (v0.isBoundary && !v1.isBoundary && v0.N_Fids.size() + v1.N_Fids.size() - 4 > Simplifier::maxValence)
        return false;
    if (!v0.isBoundary && v1.isBoundary && v0.N_Fids.size() + v1.N_Fids.size() - 4 > Simplifier::maxValence)
        return false;
    if (v0.isBoundary && !v1.isBoundary && v0.N_Fids.size() + v1.N_Fids.size() - 4 < Simplifier::minValence)
        return false;
    if (!v0.isBoundary && v1.isBoundary && v0.N_Fids.size() + v1.N_Fids.size() - 4 < Simplifier::minValence)
        return false;

    if ((v0.idealValence >= 3 || v1.idealValence >= 3) && v0.isBoundary && v1.isBoundary &&
            v0.N_Fids.size() + v1.N_Fids.size() - 2 < Simplifier::minValence)
        return false;

	if (count > 1)
	    return false;
	if (labels.size() == 1 && v.type == REGULAR && v0.type == CORNER && v1.type == FEATURE)
	    return false;
	if (labels.size() == 1 && v.type == REGULAR && v0.type == FEATURE && v1.type == CORNER)
	    return false;
	if (labels.size() > 2)
	    return false;
	if (v.isSpecial)
	    return false;
	if (labels.size() == 2 && v.isSpecial)
	    return false;
	if (labels.size() == 2 && !v.isCorner)
	    return false;
	if (labels.size() == 2 && (v0.isCorner || v1.isCorner))
	    return false;
	if (labels.size() == 2 && (/*v.isCorner || */v.labels.find(v0.label) == v.labels.end() || v.labels.find(v1.label) == v.labels.end()))
	    return false;
	//if ((v0.type == FEATURE || v1.type == FEATURE || v.type == FEATURE) && v0.N_Fids.size() + v1.N_Fids.size() - 4 <= Simplifier::minValence)
	//    return false;
	  if (!v0.isBoundary && !v1.isBoundary && (v0.type == FEATURE || v1.type == FEATURE)
	      && v0.N_Fids.size() + v1.N_Fids.size() - 4 <= Simplifier::minValence) return false;
	  if (!v0.isBoundary && v1.isBoundary && (v0.type == FEATURE || v1.type == FEATURE)
	      && v0.N_Fids.size() + v1.N_Fids.size() - 2 <= Simplifier::minValence) return false;
	  if (v0.isBoundary && !v1.isBoundary && (v0.type == FEATURE || v1.type == FEATURE)
	      && v0.N_Fids.size() + v1.N_Fids.size() - 2 <= Simplifier::minValence) return false;
	return true;
}

void Simplifier::collapse_vids_with_feature_preserved(std::vector<size_t>& vids, size_t target_vid) {
//	auto& v0 = mesh.V.at(vids[0]);
//	auto& v1 = mesh.V.at(vids[1]);
//	auto& v = mesh.V.at(target_vid);
//	if (v0.type == MAXID) v0.type = REGULAR;
//	if (v1.type == MAXID) v1.type = REGULAR;
//	if (v.type == MAXID) v.type = REGULAR;
//	if (/*v0.type <= FEATURE && v1.type <= FEATURE && */v.type == CORNER) {
//		// v.type = CORNER;
//		v0.type = REGULAR;
//		v0.label = MAXID;
//		v1.type = REGULAR;
//		v1.label = MAXID;
//	} else if (v0.type == CORNER/* && v1.type <= FEATURE && v.type <= FEATURE*/) {
//		target_vid = vids[0];
//		//mesh.V.at(target_vid).type = CORNER;
//		v.type = REGULAR;
//		v.label = MAXID;
//		v1.type = REGULAR;
//		v1.label = MAXID;
//	} else if (/*v0.type <= FEATURE && */v1.type == CORNER/* && v.type <= FEATURE*/) {
//		target_vid = vids[1];
//		// mesh.V.at(target_vid).type = CORNER;
//		v.type = REGULAR;
//		v.label = MAXID;
//		v0.type = REGULAR;
//		v0.label = MAXID;
//	} else if (v0.type == REGULAR && v1.type == REGULAR && v.type == REGULAR) {
//		//v.type = REGULAR;
//	} else if (v0.type == REGULAR && v1.type == REGULAR && v.type == REGULAR) {
//		//v.type = REGULAR;
//	} else if (v0.type == REGULAR && v1.type == REGULAR && v.type == FEATURE) {
//		//v.type = FEATURE;
//	} else if (v0.type == REGULAR && v1.type == FEATURE && v.type == FEATURE) {
//		//v.type = FEATURE;
//		v1.type = REGULAR;
//		v1.label = MAXID;
//	} else if (v0.type == FEATURE && v1.type == REGULAR && v.type == FEATURE) {
//		//v.type = FEATURE;
//		v0.type = REGULAR;
//		v0.label = MAXID;
//	} else if (v0.type == REGULAR && v1.type == FEATURE && v.type == REGULAR) {
//		target_vid = vids[1];
//		// mesh.V.at(target_vid).type = FEATURE;
//	} else if (v0.type == FEATURE && v1.type == REGULAR && v.type == REGULAR) {
//		target_vid = vids[0];
//		// mesh.V.at(target_vid).type = FEATURE;
//	} else if (v0.type <= FEATURE && v1.type <= FEATURE && v.type == FEATURE) {
//		// v.type = FEATURE;
//		v0.type = REGULAR;
//		v0.label = MAXID;
//		v1.type = REGULAR;
//		v1.label = MAXID;
//	} else {
//	    std::cerr << "###################################### Err in " << __LINE__ << std::endl;
//	}

      auto& v0 = mesh.V.at(vids[0]);
      auto& v1 = mesh.V.at(vids[1]);
      auto& v = mesh.V.at(target_vid);
	//   if (v0.type == CORNER || v1.type == CORNER) {
	// 	  std::cout << "how is this a corner?" << std::endl;
	// 	  std::cout << vids[0] << " " << vids[1] << " " << target_vid << std::endl;
	//   }
      if (v0.type == CORNER && v0.isBoundary) std::swap(vids[0], target_vid);
      else if (v1.type == CORNER && v1.isBoundary) std::swap(vids[1], target_vid);
	for (auto vid : vids) {
		Collapse(vid, target_vid);
	}
}

std::set<size_t> Simplifier::get_regionVids(const std::vector<size_t>& linkVids) {
	std::set<size_t> res;
	for (auto vid : linkVids)
		for (auto& nvid : mesh.V.at(vid).N_Vids)
			res.insert(nvid);
	return res;
}

bool Simplifier::can_collapse_with_feature_preserved(const std::vector<size_t>& linkVids, const std::vector<size_t>& linkEids) {
	std::set<size_t> vids;
	for (size_t i = 0; i < linkEids.size(); ++i) {
		auto vid = linkVids[i];
		auto eid = linkEids[i];
		auto p = get_collapse_vids(vid, eid);
		if (!can_collapse_vids_with_feature_preserved(p, vid)) return false;
		vids.insert(p.begin(), p.end());
	}
	auto vid = linkVids.back();
	auto eid = linkEids.back();
	auto p = get_collapse_vids(vid, eid);
	if (!can_collapse_vids_with_feature_preserved(p, vid)) return false;
	vids.insert(p.begin(), p.end());
	if (vids.size() < 2 * linkVids.size()) return false; // tangent;
	return true;
}

bool Simplifier::can_collapse_with_feature_preserved(const std::vector<size_t>& linkVids, const std::vector<size_t>& linkEids,
	size_t v_front_fvid, size_t v_back_fvid) {
	if (!can_collapse_with_feature_preserved(linkVids, linkEids)) return false;
	if (Simplifier::CONFORMAL) {
		auto& v_front_fv = mesh.V.at(v_front_fvid);
		auto& v_back_fv = mesh.V.at(v_back_fvid);
		if (v_front_fv.type == FEATURE && v_front_fv.N_Fids.size() == 4) return false;
		if (v_back_fv.type == FEATURE && v_back_fv.N_Fids.size() == 4) return false;
		//if (v_front_fv.isSingularity || v_back_fv.isSingularity) return false;
		if (v_front_fv.type == CORNER || v_back_fv.type == CORNER) return false;
		{
			auto& front_v = mesh.V.at(linkVids.front());
			auto& back_v = mesh.V.at(linkVids.back());
			//if (front_v.isSingularity || back_v.isSingularity) return false;
			if (front_v.type == CORNER || back_v.type == CORNER) return false;
		}
	}

	return true;
}

bool Simplifier::can_collapse_with_feature_preserved(const std::vector<size_t>& linkVids, const std::vector<size_t>& linkEids,
    size_t v_front_fvid) {
    if (!can_collapse_with_feature_preserved(linkVids, linkEids)) return false;
    if (Simplifier::CONFORMAL) {
        auto& v_front_fv = mesh.V.at(v_front_fvid);
        if (v_front_fv.type == FEATURE && v_front_fv.N_Fids.size() == 4) return false;
        if (v_front_fv.type == CORNER) return false;
        {
            auto& front_v = mesh.V.at(linkVids.front());
            auto& back_v = mesh.V.at(linkVids.back());
            if (front_v.type == CORNER || back_v.type == CORNER) return false;
        }
    }

    return true;
}

void Simplifier::collapse_with_feature_preserved(const std::vector<size_t>& linkVids, const std::vector<size_t>& linkEids) {
	for (size_t i = 0; i < linkEids.size(); ++i) {
		auto vid = linkVids[i];
		auto eid = linkEids[i];
		auto p = get_collapse_vids(vid, eid);
		// std::cout << "p size: " << p.size() << std::endl;
		collapse_vids_with_feature_preserved(p, vid);
	}
	auto vid = linkVids.back();
	auto eid = linkEids.back();
	auto p = get_collapse_vids(vid, eid);
	// std::cout << "p size: " << p.size() << std::endl;
	collapse_vids_with_feature_preserved(p, vid);
}

void Simplifier::WriteSharpEdgesVtk(const char* filename, const Mesh& m_mesh, const std::vector<FeatureLine>& featureLines) {
	const std::vector<Vertex>& V = m_mesh.V;
	const std::vector<Edge>& E = m_mesh.E;

	std::ofstream ofs(filename);
	ofs << "# vtk DataFile Version 2.0" << std::endl
		<< filename << std::endl
		<< "ASCII" << std::endl << std::endl
		<< "DATASET UNSTRUCTURED_GRID" << std::endl;
	ofs << "POINTS " << V.size() << " float" << std::endl;
	for (size_t i = 0; i < V.size(); i++)
		ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << std::endl;
	size_t numOfSharpVertices = 0;
	for (size_t i = 0; i < featureLines.size(); i++) {
		const FeatureLine& fl = featureLines.at(i);
		numOfSharpVertices += fl.Vids.size();
	}

	ofs << "CELLS " << featureLines.size() << " " << numOfSharpVertices + featureLines.size() << std::endl;
	for (size_t i = 0; i < featureLines.size(); i++) {
		const FeatureLine& fl = featureLines.at(i);
		ofs << fl.Vids.size();
		for (size_t j = 0; j < fl.Vids.size(); j++) {
			const size_t vid = fl.Vids.at(j);
			ofs << " " << vid;
		}
		ofs << std::endl;
	}

	ofs << "CELL_TYPES " << featureLines.size() << std::endl;
	for (size_t i = 0; i < featureLines.size(); i++) {
		ofs << 4 << std::endl;
	}

	ofs << "CELL_DATA " << featureLines.size() << std::endl
		<< "SCALARS " << " Feature" << " int 1" << std::endl
		<< "LOOKUP_TABLE default" << std::endl;
	for (size_t i = 0; i < featureLines.size(); i++) {
		ofs << i << std::endl;
	}
}

void Simplifier::get_feature() {
	// cos10 = 0.984807753; cos15 = 0.965925826; cos20 = 0.939692621; cos25 = 0.906307787; cos30 = 0.866025404;
	const double coslocalangle = cos((180.0 - angle) * PI / 180.0);//cos(angle);
	std::cout << "coslocalangle = " << coslocalangle << std::endl;
	mesh.SetCosAngleThreshold(coslocalangle);
	mesh.GetNormalOfSurfaceFaces();
	mesh.GetNormalOfSurfaceVertices();
	Patches patches(mesh);
	//patches.SetGlobalCosAngle(coslocalangle);
	patches.Extract();
	{
		MeshFileWriter sharpEdgesFileWriter(mesh, "SharpEdges.vtk");
		sharpEdgesFileWriter.WriteSharpEdgesVtk();
	}
	for (auto& v : mesh.V) {
		int count = 0;
		for (auto eid : v.N_Eids)
			if (mesh.E.at(eid).isSharpFeature) ++count;
		if (count == 1) {
			v.type = CORNER;
			v.isCorner = true;
			v.isSpecial = true;
			v.label = MAXID;
		}
		if (v.isCorner) {
			v.type = CORNER;
		}
	}
	// std::vector<size_t> userCorners = { 296, 312, 441, 426 };
	for (auto vid : userCorners) {
		mesh.V.at(vid).type = CORNER;
		mesh.V.at(vid).isCorner = true;
		mesh.V.at(vid).label = MAXID;
	}
	//std::vector<size_t> canceledCorners = { 310, 314, 443 };
	for (auto vid : canceledCorners) {
		mesh.V.at(vid).type = FEATURE;
		mesh.V.at(vid).isCorner = false;
	}
	{
		std::vector<size_t> cornerIds;
		for (const auto& v : mesh.V)
			if (v.isCorner) cornerIds.push_back(v.id);
		MeshFileWriter writer(mesh, "Corners.vtk");
		writer.WriteVerticesVtk(cornerIds);
	}
	{
		std::vector<size_t> faceids;
		for (auto& f : mesh.F)
			if (f.isBoundary) faceids.push_back(f.id);
		MeshFileWriter facesFileWriter(mesh, "Patches.vtk");
		facesFileWriter.WriteFacesVtk(faceids);
	}
	{
		// std::cout << "Before feature lines" << std::endl;
		mesh.LabelSharpEdges(true);
		// for (auto& e : mesh.E) e.isSharpFeature = copy[e.id];
		std::cout << "num of sharp edges: " << mesh.numOfSharpEdges << std::endl;
		std::vector<FeatureLine> featureLines(mesh.numOfSharpEdges, FeatureLine(mesh));
		for (size_t i = 0; i < mesh.numOfSharpEdges; i++) {
			featureLines.at(i).Extract(i);
		}
		WriteSharpEdgesVtk("FeatureLines.vtk", mesh, featureLines);
	}
	std::set<size_t> sharpEdgeVids;
	for (auto& e : mesh.E)
		if (e.isSharpFeature) sharpEdgeVids.insert(e.Vids.begin(), e.Vids.end());
	for (auto& v : mesh.V)
		if (v.isCorner) v.label = mesh.numOfSharpEdges;
		else if (sharpEdgeVids.find(v.id) == sharpEdgeVids.end()) v.label = MAXID;
	{
		MeshFileWriter writer(mesh, "VertexFeature0.vtk");
		writer.WriteVertexFeatureVtk();
	}

	for (auto& v : mesh.V) {
	    if (v.isCorner && !v.isSpecial) {	
	        if (is_convex(v, v.N_Fids)) v.isConvex = true;
	        v.idealValence = get_ideal_valence(v, v.N_Fids);
	    }
//	    if (v.type == FEATURE) {
//	        v.idealValence = 2;
//	        v.type = CORNER;
//	        v.isCorner = true;
//	        v.labels.insert(v.label);
//	        v.label = MAXID;
//	    }
	}
	origMesh = mesh;
}

void Simplifier::update(std::set<size_t>& canceledFids) {
	std::vector<size_t> FaceIds;
	FaceIds.reserve(mesh.F.size());
	for (auto& f : mesh.F)
		if (canceledFids.find(f.id) == canceledFids.end())
			FaceIds.push_back(f.id);
	std::vector<Vertex> newV(mesh.V.size());
//	std::vector<Face> newF(FaceIds.size());
	std::vector<Cell> newC(FaceIds.size());
	for (size_t i = 0; i < mesh.V.size(); ++i) {
	    auto& v = mesh.V.at(i);
	    auto& newv = newV.at(i);
		newv.id = i;
		newv = v.xyz();
		newv.type = v.type;
		newv.isCorner = v.isCorner;
		newv.isConvex = v.isConvex;
		newv.label = v.label;
		newv.patch_id = v.patch_id;
		newv.isSpecial = v.isSpecial;
		newv.labels = v.labels;
		newv.patch_ids = v.patch_ids;
		newv.idealValence = v.idealValence;
		newv.prescribed_length = v.prescribed_length;
		newv.smoothLocal = v.smoothLocal;
		// std::cout << v.smoothLocal << " ";
		// std::cout << newv.prescribed_length << " ";
	}
	// std::cout << "END" << std::endl;
//	for (size_t i = 0; i < FaceIds.size(); ++i) {
//		newF.at(i).id = i;
//		newF.at(i).Vids = mesh.F.at(FaceIds[i]).Vids;
//	}
	for (size_t i = 0; i < FaceIds.size(); ++i) {
		newC.at(i).id = i;
		newC.at(i).Vids = mesh.F.at(FaceIds[i]).Vids;
		newC.at(i).cellType = VTK_QUAD;
	}
	mesh.V.clear();
	mesh.E.clear();
	mesh.F.clear();
	mesh.C.clear();

	mesh.V = newV;
//	mesh.F = newF;
	mesh.C = newC;
	canceledFids.clear();
}

void Simplifier::update(const BaseComplexQuad& baseComplexQuad) {
	std::vector<Vertex> newV(baseComplexQuad.componentV.size());
	std::vector<Face> newF(baseComplexQuad.componentF.size());
	std::vector<Cell> newC(baseComplexQuad.componentF.size());
	for (auto& componentV : baseComplexQuad.componentV) {
		auto vid = baseComplexQuad.Vids.at(componentV.id);
		auto& newv = newV.at(componentV.id);
		auto& v = mesh.V.at(vid);
		newv.id = componentV.id;
		newv = v.xyz();
		newv.type = v.type;
		newv.isCorner = v.isCorner;
		newv.label = v.label;
		newv.patch_id = v.patch_id;
		newv.isSpecial = v.isSpecial;
		newv.labels = v.labels;
		newv.patch_ids = v.patch_ids;
	}
	for (auto& component_f : baseComplexQuad.componentF) {
		newF.at(component_f.id).id = component_f.id;
		newF.at(component_f.id).Vids = component_f.Vids;
	}
	for (auto& component_f : baseComplexQuad.componentF) {
		newC.at(component_f.id).id = component_f.id;
		newC.at(component_f.id).Vids = component_f.Vids;
		newC.at(component_f.id).cellType = VTK_QUAD;
	}
	mesh.V.clear();
	mesh.E.clear();
	mesh.F.clear();
	mesh.C.clear();

	mesh.V = newV;
	mesh.F = newF;
	mesh.C = newC;
}

size_t Simplifier::get_diagnal_vid(size_t vid, const std::vector<size_t>& fids, size_t end_vid) {
	size_t res = MAXID;
	std::set<size_t> fvids;
	for (auto fid : fids) {
		auto& f = mesh.F.at(fid);
		fvids.insert(f.Vids.begin(), f.Vids.end());
	}
	auto& v = mesh.V.at(vid);
	for (auto fid : v.N_Fids) {
		if (fid == fids[0] || fid == fids[1]) continue;
		auto& f = mesh.F.at(fid);
		bool found = false;
		for (auto fvid : f.Vids) {
			if (fvid == vid) continue;
			auto iter = fvids.find(fvid);
			if (iter != fvids.end() && *iter == end_vid) {
				found = true;
				res = get_diagnal_vid(vid, fid);
				break;
			}
		}
		if (found) return res;
	}
	std::cerr << "ERROR get_diagnal_vid\n";
	return res;
}

size_t Simplifier::get_id(size_t singular_vid, size_t target_vid1, size_t target_vid2) {
	auto& v_front = mesh.V.at(singular_vid);
	for (auto nvid : v_front.N_Vids) {
		auto& nv = mesh.V.at(nvid);
		int count = 0;
		for (auto vid : nv.N_Vids) {
			if (vid == target_vid1) ++count;
			if (vid == target_vid2) ++count;
		}
		if (count == 2) return nvid;
	}
	std::cerr << "ERROR get_id\n";
	return MAXID;
}

std::vector<size_t> Simplifier::get_insert_vids(const std::vector<size_t>& linkVids, const std::vector<size_t>& linkVids1) {
	std::vector<size_t> res;
	for (auto i = 0; i < linkVids.size(); ++i) {
		auto& v0 = mesh.V.at(linkVids.at(i));
		auto& v1 = mesh.V.at(linkVids1.at(i));
		Vertex v = 0.5 * (v0 + v1);
		v.id = mesh.V.size();
		if (v0.type >= FEATURE && v1.type >= FEATURE) {
			v.type = FEATURE;
			if (v0.type == FEATURE && v1.type == FEATURE) {
				v.label = v0.label;
				v.patch_ids = v0.patch_ids;
			} else if (v0.type > FEATURE && v1.type == FEATURE) {
				v.label = v1.label;
				v.patch_ids = v1.patch_ids;
			} else if (v0.type = FEATURE && v1.type > FEATURE) {
				v.label = v0.label;
				v.patch_ids = v0.patch_ids;
			} else if (v0.type > FEATURE && v1.type > FEATURE) {
				auto intersetion = Util::get_intersect(v0.labels, v1.labels);
				v.label = *intersetion.begin();
				v.patch_ids = Util::get_intersect(v0.patch_ids, v1.patch_ids);
			}
		} else if (v0.type >= FEATURE && v1.type == REGULAR) {
			v.type = REGULAR;
			v.patch_id = v1.patch_id;
		} else if (v1.type >= FEATURE && v0.type == REGULAR) {
			v.type = REGULAR;
			v.patch_id = v0.patch_id;
		} else if (v0.type == REGULAR && v1.type == REGULAR) {
			v.type = REGULAR;
			v.patch_id = v0.patch_id;
		}
		mesh.V.push_back(v);
		res.push_back(v.id);
	}
	return res;
}

std::vector<size_t> Simplifier::get_insert_fids(const std::vector<size_t>& linkVids1, const std::vector<size_t>& linkVids2) {
	std::vector<size_t> res;
	Face f;
	for (auto i = 1; i < linkVids1.size(); ++i) {
		auto vid0 = linkVids1[i - 1];
		auto vid1 = linkVids2[i - 1];
		auto vid2 = linkVids2[i];
		auto vid3 = linkVids1[i];
		f.Vids = { vid0, vid1, vid2, vid3 };
		f.id = mesh.F.size();
		mesh.F.push_back(f);
		res.push_back(f.id);
	}
	return res;
}

bool Simplifier::is_neighbor(const Vertex& v, size_t vid) {
	for (auto nvid : v.N_Vids)
		if (nvid == vid) return false;
	return true;
}

bool Simplifier::can_split_with_feature_preserved(size_t vid, size_t eid) {
	std::set<size_t> vids;
	auto& v = mesh.V.at(vid);
	auto& e = mesh.E.at(eid);
	for (auto fid : e.N_Fids) {
		auto& f = mesh.F.at(fid);
		vids.insert(f.Vids.begin(), f.Vids.end());
	}
	for (auto nvid : v.N_Vids) {
		if (vids.find(nvid) != vids.end()) continue;
		auto& nv = mesh.V.at(nvid);
		if (nv.type >= FEATURE) return false;
		//if (nv.type >= FEATURE && v.type >= FEATURE/* && nv.label == v.label*/) return false;
	}
	return true;
}

bool Simplifier::can_split_with_feature_preserved(const std::vector<size_t>& linkVids, const std::vector<size_t>& linkEids,
	size_t v_front_fvid, size_t v_back_fvid) {
		{
			auto vid = linkVids.front();
			auto eid = linkEids.front();
			if (!can_split_with_feature_preserved(vid, eid)) return false;
		}
		{
			auto vid = linkVids.back();
			auto eid = linkEids.back();
			//if (mesh.V.at(vid).type >= FEATURE) return false;
			if (!can_split_with_feature_preserved(vid, eid)) return false;
		}
		if (Simplifier::CONFORMAL) {
			auto& v_front_fv = mesh.V.at(v_front_fvid);
			auto& v_back_fv = mesh.V.at(v_back_fvid);
			if (v_front_fv.type == FEATURE && v_front_fv.N_Fids.size() == 4) return false;
			if (v_back_fv.type == FEATURE && v_back_fv.N_Fids.size() == 4) return false;
		}
		return true;
}
bool Simplifier::split_with_feature_preserved(const std::vector<size_t>& linkVids, const std::vector<size_t>& linkEids,
	size_t v_front_fvid, size_t v_back_fvid) {
	if (!can_split_with_feature_preserved(linkVids, linkEids, v_front_fvid, v_back_fvid)) return false;

	std::vector<size_t> vids1, vids2;
	{
		auto vid = linkVids.front();
		auto eid = linkEids.front();
		auto p = get_split_vids(vid, eid);
		vids1.push_back(p.front());
		vids2.push_back(p.back());
	}
	for (size_t i = 1; i < linkEids.size(); ++i) {
		auto vid = linkVids[i];
		auto eid = linkEids[i];
		auto p = get_split_vids(vid, eid);
		auto vid1 = vids1.back();
		auto& v1 = mesh.V.at(vid1);
		std::set<size_t> nvids1(v1.N_Vids.begin(), v1.N_Vids.end());
		if (nvids1.find(p.front()) != nvids1.end()) {
			vids1.push_back(p.front());
			vids2.push_back(p.back());
		} else if (nvids1.find(p.back()) != nvids1.end()) {
			vids1.push_back(p.back());
			vids2.push_back(p.front());
		} else {
			std::set<size_t> canceledFids;
			for (auto vid : linkVids) {
				auto& v = mesh.V.at(vid);
				canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
			}
			MeshFileWriter writer(mesh, "canceledFids.vtk");
			writer.WriteFacesVtk(canceledFids);
			std::cerr << "Err in split_with_feature_preserved!\n";
		}
		if (i > 1 && is_neighbor(mesh.V.at(linkVids.front()), p.front())) return false;
		if (i > 1 && is_neighbor(mesh.V.at(linkVids.front()), p.back())) return false;
		if (i < linkEids.size() - 1 && is_neighbor(mesh.V.at(linkVids.back()), p.front())) return false;
		if (i < linkEids.size() - 1 && is_neighbor(mesh.V.at(linkVids.back()), p.back())) return false;
	}
	{
		auto vid = linkVids.back();
		auto eid = linkEids.back();
		auto p = get_split_vids(vid, eid);
		auto vid1 = vids1.back();
		auto& v1 = mesh.V.at(vid1);
		std::set<size_t> nvids1(v1.N_Vids.begin(), v1.N_Vids.end());
		if (nvids1.find(p.front()) != nvids1.end()) {
			vids1.push_back(p.front());
			vids2.push_back(p.back());
		} else if (nvids1.find(p.back()) != nvids1.end()) {
			vids1.push_back(p.back());
			vids2.push_back(p.front());
		} else
			std::cerr << "Err in split_with_feature_preserved!\n";
	}

	auto insertVids1 = get_insert_vids(linkVids, vids1);
	auto insertVids2 = get_insert_vids(linkVids, vids2);

	auto v_front_eid = linkEids.front();
	auto v_back_eid = linkEids.back();
	auto& v_front_fids = mesh.E.at(v_front_eid).N_Fids;
	auto& v_back_fids = mesh.E.at(v_back_eid).N_Fids;
	auto diag_vid1 = get_diagnal_vid(linkVids.front(), v_front_fids, vids1.front());
	auto diag_vid2 = get_diagnal_vid(linkVids.front(), v_front_fids, vids2.front());
	if (diag_vid1 == MAXID) return false;
	if (diag_vid2 == MAXID) return false;
	vids1.insert(vids1.begin(), diag_vid1);
	vids2.insert(vids2.begin(), diag_vid2);
	diag_vid1 = get_diagnal_vid(linkVids.back(), v_back_fids, vids1.back());
	diag_vid2 = get_diagnal_vid(linkVids.back(), v_back_fids, vids2.back());
	if (diag_vid1 == MAXID) return false;
	if (diag_vid2 == MAXID) return false;
	vids1.push_back(diag_vid1);
	vids2.push_back(diag_vid2);
	{
		auto id = get_id(linkVids.front(), vids1.front(), v_front_fvid);
		vids1.insert(vids1.begin(), id);
		id = get_id(linkVids.front(), vids2.front(), v_front_fvid);
		vids2.insert(vids2.begin(), id);

		id = get_id(linkVids.back(), vids1.back(), v_back_fvid);
		vids1.push_back(id);
		id = get_id(linkVids.back(), vids2.back(), v_back_fvid);
		vids2.push_back(id);
	}

	//auto vids = linkVids;
	//vids.insert(vids.begin(), v_front_fvid);
	//vids.push_back(v_back_fvid);

	std::vector<size_t> linkVids1, linkVids2, linkVids3, linkVids4, linkVids5;
	std::copy(vids1.begin() + 1, vids1.begin() + vids1.size() - 1, std::back_inserter(linkVids1));

	linkVids2.push_back(vids1.front());
	std::copy(insertVids1.begin(), insertVids1.end(), std::back_inserter(linkVids2));
	linkVids2.push_back(vids1.back());

	linkVids3.push_back(v_front_fvid);
	std::copy(linkVids.begin(), linkVids.end(), std::back_inserter(linkVids3));
	linkVids3.push_back(v_back_fvid);

	linkVids4.push_back(vids2.front());
	std::copy(insertVids2.begin(), insertVids2.end(), std::back_inserter(linkVids4));
	linkVids4.push_back(vids2.back());

	std::copy(vids2.begin() + 1, vids2.begin() + vids2.size() - 1, std::back_inserter(linkVids5));

	for (auto i = 0; i < linkVids1.size(); ++i) {
		if (linkVids1[i] == linkVids2[i] && linkVids2[i] == linkVids3[i] &&
			linkVids3[i] == linkVids4[i] && linkVids4[i] == linkVids5[i]) {
			std::cerr << "adsfa\n";
			return false;
		}
	}
	{
		std::set<size_t> s(linkVids1.begin(), linkVids1.end());
		s.insert(linkVids2.begin(), linkVids2.end());
		s.insert(linkVids3.begin(), linkVids3.end());
		s.insert(linkVids4.begin(), linkVids4.end());
		s.insert(linkVids5.begin(), linkVids5.end());
		if (s.size() != 5 * linkVids1.size()) {
			std::cerr << "Error afd;\n";
			std::set<size_t> canceledFids;
			for (auto vid : linkVids) {
				auto& v = mesh.V.at(vid);
				canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
			}
			MeshFileWriter writer(mesh, "canceledFids.vtk");
			writer.WriteFacesVtk(canceledFids);
			return false;
		}
	}
	auto fids = get_insert_fids(linkVids1, linkVids2);
	fids = get_insert_fids(linkVids2, linkVids3);
	fids = get_insert_fids(linkVids3, linkVids4);
	fids = get_insert_fids(linkVids4, linkVids5);

	return true;
}

void Simplifier::strict_simplify(BaseComplexQuad& baseComplex, std::set<size_t>& canceledFids) {
	size_t id = 0;
	for (const auto& link : baseComplex.separatedVertexIdsLink) {
		const auto& linkEids = baseComplex.separatedEdgeIdsLink.at(id);
		auto& v_front = mesh.V.at(link.front());
		auto& v_back = mesh.V.at(link.back());
		if (v_front.isBoundary || v_back.isBoundary) {
			++id;
			continue;
		}
		if ((v_front.N_Fids.size() <= 5 && v_back.N_Fids.size() <= 5) && v_front.N_Fids.size() != v_back.N_Fids.size()) {
			; //ofs << 0 << std::endl;
		} else if (COLLAPSE && v_front.N_Fids.size() == 3 && v_back.N_Fids.size() == 3) {
			// ofs << 1 << std::endl;
			auto v_front_fid = get_faceid(v_front.id, link[1]);
			auto v_front_fvid = get_diagnal_vid(v_front.id, v_front_fid);
			auto v_back_fid = get_faceid(v_back.id, link[link.size() - 2]);
			auto v_back_fvid = get_diagnal_vid(v_back.id, v_back_fid);
			bool condition = false;
			if (mesh.V.at(v_front_fvid).N_Fids.size() == 5 && mesh.V.at(v_back_fvid).N_Fids.size() == 5) {
				if (can_collapse_with_feature_preserved(link, linkEids, v_front_fvid, v_back_fvid)) {
					for (auto vid : link) {
						auto& v = mesh.V.at(vid);
						canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
					}
					collapse_with_feature_preserved(link, linkEids);
					break;
				}
			}
		} else if (Simplifier::SPLIT && v_front.N_Fids.size() == 5 && v_back.N_Fids.size() == 5) {
			//ofs << 2 << std::endl;
			auto v_front_eid = baseComplex.separatedEdgeIdsLink.at(id).front();
			auto v_back_eid = baseComplex.separatedEdgeIdsLink.at(id).back();
			auto& v_front_fids = mesh.E.at(v_front_eid).N_Fids;
			auto& v_back_fids = mesh.E.at(v_back_eid).N_Fids;

			auto v_front_fvid = get_diagnal_vid(v_front.id, v_front_fids);
			auto v_back_fvid = get_diagnal_vid(v_back.id, v_back_fids);
			bool condition = false;
			if (mesh.V.at(v_front_fvid).N_Fids.size() == 3 && mesh.V.at(v_back_fvid).N_Fids.size() == 3) {
				if (!split_with_feature_preserved(link, linkEids, v_front_fvid, v_back_fvid)) {
					++id;
					continue;
				}
				std::cout << "STRICT SPLIT OPERATION" << std::endl;
				for (auto vid : link) {
					auto& v = mesh.V.at(vid);
					canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
				}
				break;
			}
		} else {
			; //ofs << 3 << std::endl;
		}
		++id;
	}

}

void Simplifier::strict_simplify_reverse(BaseComplexQuad& baseComplex, std::set<size_t>& canceledFids) {
    for (int id = baseComplex.separatedVertexIdsLink.size(); --id >= 0;) {
        const auto& link = baseComplex.separatedVertexIdsLink.at(id);
        const auto& linkEids = baseComplex.separatedEdgeIdsLink.at(id);
        auto& v_front = mesh.V.at(link.front());
        auto& v_back = mesh.V.at(link.back());
        if (v_front.isBoundary || v_back.isBoundary) continue;
        if ((v_front.N_Fids.size() <= 5 && v_back.N_Fids.size() <= 5) && v_front.N_Fids.size() != v_back.N_Fids.size()) {
            ; //ofs << 0 << std::endl;
        } else if (COLLAPSE && v_front.N_Fids.size() == 3 && v_back.N_Fids.size() == 3) {
            // ofs << 1 << std::endl;
            auto v_front_fid = get_faceid(v_front.id, link[1]);
            auto v_front_fvid = get_diagnal_vid(v_front.id, v_front_fid);
            auto v_back_fid = get_faceid(v_back.id, link[link.size() - 2]);
            auto v_back_fvid = get_diagnal_vid(v_back.id, v_back_fid);
            bool condition = false;
            if (mesh.V.at(v_front_fvid).N_Fids.size() == 5 && mesh.V.at(v_back_fvid).N_Fids.size() == 5) {
                if (can_collapse_with_feature_preserved(link, linkEids, v_front_fvid, v_back_fvid)) {
                    for (auto vid : link) {
                        auto& v = mesh.V.at(vid);
                        canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
                    }
                    collapse_with_feature_preserved(link, linkEids);
                    break;
                }
            }
        } else if (Simplifier::SPLIT && v_front.N_Fids.size() == 5 && v_back.N_Fids.size() == 5) {
            //ofs << 2 << std::endl;
            auto v_front_eid = baseComplex.separatedEdgeIdsLink.at(id).front();
            auto v_back_eid = baseComplex.separatedEdgeIdsLink.at(id).back();
            auto& v_front_fids = mesh.E.at(v_front_eid).N_Fids;
            auto& v_back_fids = mesh.E.at(v_back_eid).N_Fids;

            auto v_front_fvid = get_diagnal_vid(v_front.id, v_front_fids);
            auto v_back_fvid = get_diagnal_vid(v_back.id, v_back_fids);
            bool condition = false;
            if (mesh.V.at(v_front_fvid).N_Fids.size() == 3 && mesh.V.at(v_back_fvid).N_Fids.size() == 3) {
                if (!split_with_feature_preserved(link, linkEids, v_front_fvid, v_back_fvid)) {
                    ++id;
                    continue;
                }
                for (auto vid : link) {
                    auto& v = mesh.V.at(vid);
                    canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
                }
                break;
            }
        } else {
            ; //ofs << 3 << std::endl;
        }
    }
}

void Simplifier::loose_simplify(BaseComplexQuad& baseComplex, std::set<size_t>& canceledFids) {
	size_t id = -1;
	for (const auto& link : baseComplex.separatedVertexIdsLink) {
		const auto& linkEids = baseComplex.separatedEdgeIdsLink.at(++id);
		auto& v_front = mesh.V.at(link.front());
		auto& v_back = mesh.V.at(link.back());
		if (v_front.isBoundary || v_back.isBoundary) continue;
		if ((v_front.N_Fids.size() <= 5 && v_back.N_Fids.size() <= 5) && v_front.N_Fids.size() != v_back.N_Fids.size()) {
			; //ofs << 0 << std::endl;
		} else if (COLLAPSE && v_front.N_Fids.size() == 3 && v_back.N_Fids.size() == 3) {
			// ofs << 1 << std::endl;
			auto v_front_fid = get_faceid(v_front.id, link[1]);
			auto v_front_fvid = get_diagnal_vid(v_front.id, v_front_fid);
			auto v_back_fid = get_faceid(v_back.id, link[link.size() - 2]);
			auto v_back_fvid = get_diagnal_vid(v_back.id, v_back_fid);
			if (mesh.V.at(v_front_fvid).N_Fids.size() > Simplifier::minValence && mesh.V.at(v_back_fvid).N_Fids.size() > Simplifier::minValence) {
				if (can_collapse_with_feature_preserved(link, linkEids, v_front_fvid, v_back_fvid)) {
					for (auto vid : link) {
						auto& v = mesh.V.at(vid);
						canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
					}
					collapse_with_feature_preserved(link, linkEids);
					break;
				}
			}
		} else if (Simplifier::SPLIT && v_front.N_Fids.size() == 5 && v_back.N_Fids.size() == 5) {
			//ofs << 2 << std::endl;
			auto v_front_eid = baseComplex.separatedEdgeIdsLink.at(id).front();
			auto v_back_eid = baseComplex.separatedEdgeIdsLink.at(id).back();
			auto& v_front_fids = mesh.E.at(v_front_eid).N_Fids;
			auto& v_back_fids = mesh.E.at(v_back_eid).N_Fids;

			auto v_front_fvid = get_diagnal_vid(v_front.id, v_front_fids);
			auto v_back_fvid = get_diagnal_vid(v_back.id, v_back_fids);
			bool condition = false;
			if (mesh.V.at(v_front_fvid).N_Fids.size() < Simplifier::maxValence && mesh.V.at(v_back_fvid).N_Fids.size() < Simplifier::maxValence) {
				if (split_with_feature_preserved(link, linkEids, v_front_fvid, v_back_fvid)) {
				std::cout << "LOOSE SPLIT OPERATION" << std::endl;
					for (auto vid : link) {
						auto& v = mesh.V.at(vid);
						canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
					}
					break;
				}
			}
		} else {
			; //ofs << 3 << std::endl;
		}
	}
}

void Simplifier::three_connections_collapse(BaseComplexQuad& baseComplex, std::set<size_t>& canceledFids, bool looseSimplify) {
	size_t id = 0;
	if (looseSimplify) {
		id = -1;
	}
	struct collapsableThreeLink {
		std::vector<size_t> target;
		std::vector<size_t> link;
		std::vector<std::vector<size_t>> collapse;
	};
	std::vector<collapsableThreeLink> threeLinks;

	size_t numElements = 0;
	for (const auto& link : baseComplex.separatedVertexIdsLink) {
		if (looseSimplify) {
			++id;
		}
		const auto& linkEids = baseComplex.separatedEdgeIdsLink.at(id);
		auto& v_front = mesh.V.at(link.front());
		auto& v_back = mesh.V.at(link.back());
		if (v_front.isBoundary || v_back.isBoundary) {
			if (!looseSimplify) {
				++id;
			}
			continue;
		}
		if ((v_front.N_Fids.size() <= 5 && v_back.N_Fids.size() <= 5) && v_front.N_Fids.size() != v_back.N_Fids.size()) {
			; //ofs << 0 << std::endl;
		} else if (COLLAPSE && v_front.N_Fids.size() == 3 && v_back.N_Fids.size() == 3) {
			// ofs << 1 << std::endl;
			auto v_front_fid = get_faceid(v_front.id, link[1]);
			auto v_front_fvid = get_diagnal_vid(v_front.id, v_front_fid);
			auto v_back_fid = get_faceid(v_back.id, link[link.size() - 2]);
			auto v_back_fvid = get_diagnal_vid(v_back.id, v_back_fid);
			bool condition = mesh.V.at(v_front_fvid).N_Fids.size() == 5 && mesh.V.at(v_back_fvid).N_Fids.size() == 5;
			if (looseSimplify) {
				condition = mesh.V.at(v_front_fvid).N_Fids.size() > Simplifier::minValence && mesh.V.at(v_back_fvid).N_Fids.size() > Simplifier::minValence;
			}
			if (condition) {
				if (can_collapse_with_feature_preserved(link, linkEids, v_front_fvid, v_back_fvid)) {
					// for (auto vid : link) {
					// 	auto& v = mesh.V.at(vid);
					// 	canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
					// }
					// collapse_with_feature_preserved(link, linkEids);
					collapsableThreeLink l;
					l.link = link;
					for (size_t i = 0; i < linkEids.size(); ++i) {
						auto vid = link[i];
						auto eid = linkEids[i];
						auto p = get_collapse_vids(vid, eid);
						auto& v0 = mesh.V.at(p[0]);
						auto& v1 = mesh.V.at(p[1]);
						auto& v = mesh.V.at(vid);
						if (v0.type == CORNER && v0.isBoundary) std::swap(p[0], vid);
						else if (v1.type == CORNER && v1.isBoundary) std::swap(p[1], vid);
						l.target.push_back(vid);
						l.collapse.push_back(p);
					// // 	collapse_vids_with_feature_preserved(p, vid);
					}
					auto vid = link.back();
					auto eid = linkEids.back();
					auto p = get_collapse_vids(vid, eid);
					auto& v0 = mesh.V.at(p[0]);
					auto& v1 = mesh.V.at(p[1]);
					auto& v = mesh.V.at(vid);
					if (v0.type == CORNER && v0.isBoundary) std::swap(p[0], vid);
					else if (v1.type == CORNER && v1.isBoundary) std::swap(p[1], vid);
					l.target.push_back(vid);
					l.collapse.push_back(p);

					threeLinks.push_back(l);
					// collapseVids.insert(std::pair<size_t, std::vector<size_t>>(vid, p));
					// break;
				}
			}
		} else if (Simplifier::SPLIT && v_front.N_Fids.size() == 5 && v_back.N_Fids.size() == 5) {
			; //ofs << 2 << std::endl;
		} else {
			; //ofs << 3 << std::endl;
		}
		if (!looseSimplify) {
			++id;
		}
	}
	std::vector<double> ranks;
	for (auto l : threeLinks) {
		double rank = 0;
		for (int i = 0; i < l.target.size(); i++) {
			auto& v1 = mesh.V.at(l.target.at(i));
			for (auto value: l.collapse.at(i)) {
				auto& v2 = mesh.V.at(value);
				rank += glm::length(glm::dvec3(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z));
			}
		}
		ranks.push_back(rank);
	}
	// std::cout << "RANKS: " << ranks.size() << std::endl;
	std::vector<size_t> targetVidsPos;
	// std::vector<double>::iterator max_index = std::max_element(ranks.begin(), ranks.end());
    // double max_rank = (double) std::distance(ranks.begin(), max_index) + 1;
	for (int i = 0; i < ranks.size(); i++) {
        std::vector<double>::iterator index = std::max_element(ranks.begin(), ranks.end());
        // std::vector<double>::iterator index = std::min_element(ranks.begin(), ranks.end());
        targetVidsPos.push_back((size_t) std::distance(ranks.begin(), index));
        // targetVidsPos.push_back(i);
        // *index = max_rank;
        *index = -1;
		// ranks.erase(index);
		// i = 0;
    }
	// std::cout << "TARGET VIDS POS: " << targetVidsPos.size() << std::endl;
	std::vector<collapsableThreeLink> finalThreeLinks;
	for (int i = 0; i < targetVidsPos.size(); i++) {
		if (targetVidsPos.at(i) == -1) {
			continue;
		}
		collapsableThreeLink l = threeLinks.at(targetVidsPos.at(i));
		finalThreeLinks.push_back(l);
		for (int j = 0; j < threeLinks.size(); j++) {
			if (j == targetVidsPos.at(i)) {
				continue;
			}
			bool disjoint = true;
			collapsableThreeLink l2 = threeLinks.at(j);
			for (auto id: l.target) {
				if (std::find(l2.target.begin(), l2.target.end(), id) != l2.target.end()) {
					disjoint = false;
					break;
				}
				if (!disjoint) {
					break;
				}
				for (auto vec: l2.collapse) {
					if (std::find(vec.begin(), vec.end(), id) != vec.end()) {
						disjoint = false;
						break;
					}
				}
				if (!disjoint) {
					break;
				}
			}
			for (auto vec: l.collapse) {
				for (auto id: l2.target) {
					if (std::find(vec.begin(), vec.end(), id) != vec.end()) {
						disjoint = false;
						break;
					}
				}
				if (!disjoint) {
					break;
				}
				for (auto vec2: l2.collapse) {
					for (auto id: vec2) {
						if (std::find(vec.begin(), vec.end(), id) != vec.end()) {
							disjoint = false;
							break;
						}
					}
					if (!disjoint) {
						break;
					}
				}
				if (!disjoint) {
					break;
				}
			}
			if (!disjoint) {
				auto it = std::find(targetVidsPos.begin(), targetVidsPos.end(), j);
                if (it != targetVidsPos.end()) {
                    targetVidsPos.at(std::distance(targetVidsPos.begin(), it)) = -1;
                }
			}
		}
	}
	// std::cout << "FINAL THREE LINKS: " << finalThreeLinks.size() << std::endl;
	for (auto l: finalThreeLinks) {
		for (auto vid : l.link) {
			auto& v = mesh.V.at(vid);
			canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
		}
		for (int i = 0; i < l.target.size(); i++) {
			// auto& v = mesh.V.at(l.target.at(i));
			// canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
			for (auto vid: l.collapse.at(i)) {
				Collapse(vid, l.target.at(i));
			}
			// collapse_vids_with_feature_preserved(l.collapse.at(i), l.target.at(i));
		}
		// break;
	}
	// std::cout << finalThreeLinks.size() << std::endl;
	// for (auto l: finalThreeLinks) {
	// 	std::cout << "--------------LINK-------------------" << std::endl;
	// 	for (int i = 0; i < l.target.size(); i++) {
	// 		std::cout << l.target.at(i) << ":";
	// 		for (auto value: l.collapse.at(i)) {
	// 			std::cout << " " << value;
	// 		}
	// 		std::cout << std::endl;
	// 	}
	// 	std::cout << "--------------------------------------" << std::endl;
	// }
	// std::vector<size_t> target_indices;
	// std::vector<size_t> collapse_indices;
	// for (auto l: finalThreeLinks) {
	// 	for (int i = 0; i < l.target.size(); i++) {
	// 		target_indices.push_back(l.target.at(i));
	// 		for (auto value: l.collapse.at(i)) {
	// 			collapse_indices.push_back(value);
	// 		}
	// 	}
	// }
	// std::ofstream ofs("target_vids.vtk");
	// ofs << "# vtk DataFile Version 3.0\n"
	// 	<< "output.vtk\n"
	// 	<< "ASCII\n\n"
	// 	<< "DATASET UNSTRUCTURED_GRID\n";
	// ofs << "POINTS " << mesh.V.size() << " double\n";
	
	// for (size_t i = 0; i < mesh.V.size(); i++) {
	// 	ofs << mesh.V.at(i).x << " " <<  mesh.V.at(i).y << " " <<  mesh.V.at(i).z << "\n";
	// }
	// ofs << "CELLS " << target_indices.size() << " " << 2 * target_indices.size() << std::endl;
	// for (size_t i = 0; i < target_indices.size(); i++) {
	// 	ofs << "1 " << target_indices.at(i) << std::endl;
	// }
	// ofs << "CELL_TYPES " << target_indices.size() << "\n";
	// for (size_t i = 0; i < target_indices.size(); i++) {
	// 	ofs << "1" << std::endl;
	// }

	// std::ofstream ofs2("collapse_vids.vtk");
	// ofs2 << "# vtk DataFile Version 3.0\n"
	// 	<< "output.vtk\n"
	// 	<< "ASCII\n\n"
	// 	<< "DATASET UNSTRUCTURED_GRID\n";
	// ofs2 << "POINTS " << mesh.V.size() << " double\n";
	
	// for (size_t i = 0; i < mesh.V.size(); i++) {
	// 	ofs2 << mesh.V.at(i).x << " " <<  mesh.V.at(i).y << " " <<  mesh.V.at(i).z << "\n";
	// }
	// ofs2 << "CELLS " << collapse_indices.size() << " " << 2 * collapse_indices.size() << std::endl;
	// for (size_t i = 0; i < collapse_indices.size(); i++) {
	// 	ofs2 << "1 " << collapse_indices.at(i) << std::endl;
	// }
	// ofs2 << "CELL_TYPES " << collapse_indices.size() << "\n";
	// for (size_t i = 0; i < collapse_indices.size(); i++) {
	// 	ofs2 << "1" << std::endl;
	// }
}

void Simplifier::GetSeparatrixCollapseOps(BaseComplexQuad& baseComplex, bool looseSimplify, std::multiset<SimplificationOperation, bool(*)(SimplificationOperation, SimplificationOperation)>& SimplificationOps) {
	size_t id = 0;
	if (looseSimplify) id = -1;
	for (const auto& link : baseComplex.separatedVertexIdsLink) {
		if (looseSimplify) ++id;
		const auto& linkEids = baseComplex.separatedEdgeIdsLink.at(id);
		auto& v_front = mesh.V.at(link.front());
		auto& v_back = mesh.V.at(link.back());
		if (v_front.isBoundary || v_back.isBoundary) {
			if (!looseSimplify) ++id;
			continue;
		}
		if ((v_front.N_Fids.size() <= 5 && v_back.N_Fids.size() <= 5) && v_front.N_Fids.size() != v_back.N_Fids.size()) {;}
		else if (v_front.N_Fids.size() == 3 && v_back.N_Fids.size() == 3) {
			auto v_front_fid = get_faceid(v_front.id, link[1]);
			auto v_front_fvid = get_diagnal_vid(v_front.id, v_front_fid);
			auto v_back_fid = get_faceid(v_back.id, link[link.size() - 2]);
			auto v_back_fvid = get_diagnal_vid(v_back.id, v_back_fid);
			bool condition = mesh.V.at(v_front_fvid).N_Fids.size() == 5 && mesh.V.at(v_back_fvid).N_Fids.size() == 5;
			if (looseSimplify) {
				condition = mesh.V.at(v_front_fvid).N_Fids.size() > Simplifier::minValence && mesh.V.at(v_back_fvid).N_Fids.size() > Simplifier::minValence;
			}
			if (condition && can_collapse_with_feature_preserved(link, linkEids, v_front_fvid, v_back_fvid)) {
				SimplificationOperation Op;
				Op.type = looseSimplify ? "Loose_Separatrix_Collapse" : "Strict_Separatrix_Collapse";
				for (size_t i = 0; i < linkEids.size(); ++i) {
					auto vid = link[i];
					auto eid = linkEids[i];
					auto p = get_collapse_vids(vid, eid);
					CollapseVertexToTarget(p[0], vid, Op);
					CollapseVertexToTarget(p[1], vid, Op);
					Op.updatedVertexPos.push_back(0.5 * (mesh.V.at(p[0]).xyz() + mesh.V.at(p[1]).xyz()));
					Op.updateVertexIds.push_back(vid);
				}
				auto vid = link.back();
				auto eid = linkEids.back();
				auto p = get_collapse_vids(vid, eid);
				CollapseVertexToTarget(p[0], vid, Op);
				CollapseVertexToTarget(p[1], vid, Op);
				Op.updatedVertexPos.push_back(0.5 * (mesh.V.at(p[0]).xyz() + mesh.V.at(p[1]).xyz()));
				Op.updateVertexIds.push_back(vid);
				// Op.profitability /= mesh.totalArea;
				Op.profitability /= Op.n;
				// Op.profitability = 1;
				SimplificationOps.insert(Op);
			}
		}
		if (!looseSimplify) ++id;
	}
}

void Simplifier::GetHalfSeparatrixOps(BaseComplexQuad& baseComplex, std::multiset<SimplificationOperation, bool(*)(SimplificationOperation, SimplificationOperation)>& SimplificationOps) {
	size_t id = -1;
	for (auto link : baseComplex.separatedVertexIdsLink) {
        const auto& linkEids = baseComplex.separatedEdgeIdsLink.at(++id);
        auto& v_front = mesh.V.at(link.front());
        auto& v_back = mesh.V.at(link.back());
        if (!(v_front.isBoundary ^ v_back.isBoundary)) continue;
        if (v_front.N_Fids.size() == 3 && v_back.N_Fids.size() == 2) {
            bool onthesameline = true;
            for (auto nvid : v_back.N_Vids) {
                auto& nv = mesh.V.at(nvid);
                if (nv.type == CORNER) {
                    onthesameline = false;
                    break;
                }
            }
            if (!onthesameline) continue;
	        auto v_front_fid = get_faceid(v_front.id, link[1]);
            auto v_front_fvid = get_diagnal_vid(v_front.id, v_front_fid);
            auto& v_front_fv = mesh.V.at(v_front_fvid);
            if (v_front_fvid >= mesh.V.size()) {
                MeshFileWriter writer(mesh, "error.vtk");
                writer.WriteFile();
            }
            if (v_front_fv.N_Fids.size() >= Simplifier::minValence + 1) {
                if (can_collapse_with_feature_preserved(link, linkEids, v_front_fvid)) {
    				SimplificationOperation Op;
					Op.type = "Half_Separatrix_Collapse";
					for (size_t i = 0; i < linkEids.size(); ++i) {
						auto vid = link[i];
						auto eid = linkEids[i];
						auto p = get_collapse_vids(vid, eid);
						CollapseVerticesToTargetWithFeaturePreserved(p, vid, Op);
						Op.updatedVertexPos.push_back(0.5 * (mesh.V.at(p[0]).xyz() + mesh.V.at(p[1]).xyz()));
						Op.updateVertexIds.push_back(vid);
					}
					auto vid = link.back();
					auto eid = linkEids.back();
					auto p = get_collapse_vids(vid, eid);
					CollapseVerticesToTargetWithFeaturePreserved(p, vid, Op);
					Op.updatedVertexPos.push_back(0.5 * (mesh.V.at(p[0]).xyz() + mesh.V.at(p[1]).xyz()));
					Op.updateVertexIds.push_back(vid);
					// Op.profitability /= mesh.totalArea;
					Op.profitability /= Op.n;
					// Op.profitability = 1;
					SimplificationOps.insert(Op);
                }
            }
        } else if (!v_back.isBoundary && v_back.N_Fids.size() == 3 && v_front.isBoundary && v_front.N_Fids.size() == 2) {
            std::reverse(link.begin(), link.end());
            bool onthesameline = true;
            for (auto nvid : v_front.N_Vids) {
                auto& nv = mesh.V.at(nvid);
                if (nv.type == CORNER) {
                    onthesameline = false;
                    break;
                }
            }
            if (!onthesameline) continue;
	        auto v_front_fid = get_faceid(v_back.id, link[1]);
            auto v_front_fvid = get_diagnal_vid(v_back.id, v_front_fid);
            auto& v_front_fv = mesh.V.at(v_front_fvid);
            if (v_front_fv.isBoundary/* && v_front_fv.N_Fids.size() <= v_front_fv.idealValence*/) continue;
            if (v_front_fvid >= mesh.V.size()) {
                MeshFileWriter writer(mesh, "error.vtk");
                writer.WriteFile();
            }
            if (mesh.V.at(v_front_fvid).N_Fids.size() >= Simplifier::minValence + 1) {
                if (can_collapse_with_feature_preserved(link, linkEids, v_front_fvid)) {
    				SimplificationOperation Op;
					Op.type = "Half_Separatrix_Collapse";
					for (size_t i = 0; i < linkEids.size(); ++i) {
						auto vid = link[i];
						auto eid = linkEids[i];
						auto p = get_collapse_vids(vid, eid);
						CollapseVerticesToTargetWithFeaturePreserved(p, vid, Op);
						Op.updatedVertexPos.push_back(0.5 * (mesh.V.at(p[0]).xyz() + mesh.V.at(p[1]).xyz()));
						Op.updateVertexIds.push_back(vid);
					}
					auto vid = link.back();
					auto eid = linkEids.back();
					auto p = get_collapse_vids(vid, eid);
					CollapseVerticesToTargetWithFeaturePreserved(p, vid, Op);
					Op.updatedVertexPos.push_back(0.5 * (mesh.V.at(p[0]).xyz() + mesh.V.at(p[1]).xyz()));
					Op.updateVertexIds.push_back(vid);
					// Op.profitability /= mesh.totalArea;
					Op.profitability /= Op.n;
					// Op.profitability = 1;
					SimplificationOps.insert(Op);
                }
            }
        }
    }
}

void Simplifier::CollapseVertexToTarget(size_t source_vid, size_t target_vid, SimplificationOperation& Op) {
	auto& source = mesh.V.at(source_vid);
	auto& target = mesh.V.at(target_vid);
	Op.profitability += glm::distance(source.xyz(), target.xyz());
	Op.n += 1;
	for (auto fid: source.N_Fids) {
		if (std::find(target.N_Fids.begin(), target.N_Fids.end(), fid) == target.N_Fids.end()) {
			bool skip = false;
			for (auto& Op_face: Op.newFaces) {
				if (Op_face.id == fid) {
					for (int i = 0; i < Op_face.Vids.size(); i++) {
						Op_face.Vids.at(i) = Op_face.Vids.at(i) == source.id ? target.id : Op_face.Vids.at(i);
					}
					skip = true;
				}
			}
			if (skip) continue;
			Face& n_f = mesh.F.at(fid);
			Face newF;
			newF.id = n_f.id;
			for (int i = 0; i < n_f.Vids.size(); i++) {
				n_f.Vids.at(i) == source.id ? newF.Vids.push_back(target.id) : newF.Vids.push_back(n_f.Vids.at(i));
			}
			Op.newFaces.push_back(newF);
			Op.canceledFids.insert(n_f.id);
			// Op.profitability += mesh.GetQuadFaceArea(mesh.F.at(fid).Vids);
			// Op.processedFids.insert(n_f.id);
			// Op.processedFids.insert(n_f.N_Fids.begin(), n_f.N_Fids.end());
		} else {
			Op.canceledFids.insert(fid);
			// Op.processedFids.insert(fid);
			// Op.profitability += mesh.GetQuadFaceArea(mesh.F.at(fid).Vids);
		}
	}
}

void Simplifier::CollapseVerticesToTargetWithFeaturePreserved(std::vector<size_t>& source_vids, size_t target_vid, SimplificationOperation& Op) {
	auto& v0 = mesh.V.at(source_vids[0]);
	auto& v1 = mesh.V.at(source_vids[1]);
	
	if (v0.type == CORNER && v0.isBoundary) std::swap(source_vids[0], target_vid);
	else if (v1.type == CORNER && v1.isBoundary) std::swap(source_vids[1], target_vid);
	for (auto source_vid : source_vids) {
		CollapseVertexToTarget(source_vid, target_vid, Op);
	}
}

void Simplifier::five_connections_split(BaseComplexQuad& baseComplex, std::set<size_t>& canceledFids, bool looseSimplify) {
	size_t id = 0;
	if (looseSimplify) {
		id = -1;
	}
	for (const auto& link : baseComplex.separatedVertexIdsLink) {
		if (looseSimplify) {
			++id;
		}
		const auto& linkEids = baseComplex.separatedEdgeIdsLink.at(id);
		auto& v_front = mesh.V.at(link.front());
		auto& v_back = mesh.V.at(link.back());
		if (v_front.isBoundary || v_back.isBoundary) {
			if (!looseSimplify) {
				++id;
			}
			continue;
		}
		if ((v_front.N_Fids.size() <= 5 && v_back.N_Fids.size() <= 5) && v_front.N_Fids.size() != v_back.N_Fids.size()) {
			; //ofs << 0 << std::endl;
		} else if (COLLAPSE && v_front.N_Fids.size() == 3 && v_back.N_Fids.size() == 3) {
			; // ofs << 1 << std::endl;
		} else if (Simplifier::SPLIT && v_front.N_Fids.size() == 5 && v_back.N_Fids.size() == 5) {
			//ofs << 2 << std::endl;
			auto v_front_eid = baseComplex.separatedEdgeIdsLink.at(id).front();
			auto v_back_eid = baseComplex.separatedEdgeIdsLink.at(id).back();
			auto& v_front_fids = mesh.E.at(v_front_eid).N_Fids;
			auto& v_back_fids = mesh.E.at(v_back_eid).N_Fids;

			auto v_front_fvid = get_diagnal_vid(v_front.id, v_front_fids);
			auto v_back_fvid = get_diagnal_vid(v_back.id, v_back_fids);
			bool condition = false;
			if (looseSimplify) {
				if (mesh.V.at(v_front_fvid).N_Fids.size() < Simplifier::maxValence && mesh.V.at(v_back_fvid).N_Fids.size() < Simplifier::maxValence) {
					if (split_with_feature_preserved(link, linkEids, v_front_fvid, v_back_fvid)) {
						for (auto vid : link) {
							auto& v = mesh.V.at(vid);
							canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
						}
						break;
					}
				}
			} else {
				if (mesh.V.at(v_front_fvid).N_Fids.size() == 3 && mesh.V.at(v_back_fvid).N_Fids.size() == 3) {
					if (!split_with_feature_preserved(link, linkEids, v_front_fvid, v_back_fvid)) {
						++id;
						continue;
					}
					for (auto vid : link) {
						auto& v = mesh.V.at(vid);
						canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
					}
					break;
				}
			}
		} else {
			; //ofs << 3 << std::endl;
		}
		if (!looseSimplify) {
			++id;
		}
	}
}

void Simplifier::loose_simplify_reverse(BaseComplexQuad& baseComplex, std::set<size_t>& canceledFids) {
    for (int id = baseComplex.separatedVertexIdsLink.size(); --id >= 0;) {
        const auto& link = baseComplex.separatedVertexIdsLink.at(id);
        const auto& linkEids = baseComplex.separatedEdgeIdsLink.at(id);
        auto& v_front = mesh.V.at(link.front());
        auto& v_back = mesh.V.at(link.back());
        if (v_front.isBoundary || v_back.isBoundary) continue;
        if ((v_front.N_Fids.size() <= 5 && v_back.N_Fids.size() <= 5) && v_front.N_Fids.size() != v_back.N_Fids.size()) {
            ; //ofs << 0 << std::endl;
        } else if (COLLAPSE && v_front.N_Fids.size() == 3 && v_back.N_Fids.size() == 3) {
            // ofs << 1 << std::endl;
            auto v_front_fid = get_faceid(v_front.id, link[1]);
            auto v_front_fvid = get_diagnal_vid(v_front.id, v_front_fid);
            auto v_back_fid = get_faceid(v_back.id, link[link.size() - 2]);
            auto v_back_fvid = get_diagnal_vid(v_back.id, v_back_fid);
            if (mesh.V.at(v_front_fvid).N_Fids.size() > Simplifier::minValence && mesh.V.at(v_back_fvid).N_Fids.size() > Simplifier::minValence) {
                if (can_collapse_with_feature_preserved(link, linkEids, v_front_fvid, v_back_fvid)) {
                    for (auto vid : link) {
                        auto& v = mesh.V.at(vid);
                        canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
                    }
                    collapse_with_feature_preserved(link, linkEids);
                    break;
                }
            }
        } else if (Simplifier::SPLIT && v_front.N_Fids.size() == 5 && v_back.N_Fids.size() == 5) {
            //ofs << 2 << std::endl;
            auto v_front_eid = baseComplex.separatedEdgeIdsLink.at(id).front();
            auto v_back_eid = baseComplex.separatedEdgeIdsLink.at(id).back();
            auto& v_front_fids = mesh.E.at(v_front_eid).N_Fids;
            auto& v_back_fids = mesh.E.at(v_back_eid).N_Fids;

            auto v_front_fvid = get_diagnal_vid(v_front.id, v_front_fids);
            auto v_back_fvid = get_diagnal_vid(v_back.id, v_back_fids);
            bool condition = false;
            if (mesh.V.at(v_front_fvid).N_Fids.size() < Simplifier::maxValence && mesh.V.at(v_back_fvid).N_Fids.size() < Simplifier::maxValence) {
                if (split_with_feature_preserved(link, linkEids, v_front_fvid, v_back_fvid)) {
                    for (auto vid : link) {
                        auto& v = mesh.V.at(vid);
                        canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
                    }
                    break;
                }
            }
        } else {
            ; //ofs << 3 << std::endl;
        }
    }
}
#include <ctime>
#include <stdlib.h>     /* srand, rand */
int Simplifier::get_id(std::unordered_set<int>& ids, int n) {
    while (true) {
        int id = rand() % n;
        if (ids.count(id) == 0) {
            ids.insert(id);
            return id;
        }
    }
    return INT_MIN;
}
void Simplifier::loose_simplify_random(BaseComplexQuad& baseComplex, std::set<size_t>& canceledFids) {
    std::unordered_set<int> ids;
    srand (time(NULL));
    while (ids.size() < baseComplex.separatedVertexIdsLink.size()) {
    //for (int id = baseComplex.separatedVertexIdsLink.size(); --id >= 0;) {
        int id = get_id(ids, (int)baseComplex.separatedVertexIdsLink.size());
        const auto& link = baseComplex.separatedVertexIdsLink.at(id);
        const auto& linkEids = baseComplex.separatedEdgeIdsLink.at(id);
        auto& v_front = mesh.V.at(link.front());
        auto& v_back = mesh.V.at(link.back());
        if (v_front.isBoundary || v_back.isBoundary) continue;
        if ((v_front.N_Fids.size() <= 5 && v_back.N_Fids.size() <= 5) && v_front.N_Fids.size() != v_back.N_Fids.size()) {
            ; //ofs << 0 << std::endl;
        } else if (COLLAPSE && v_front.N_Fids.size() == 3 && v_back.N_Fids.size() == 3) {
            // ofs << 1 << std::endl;
            auto v_front_fid = get_faceid(v_front.id, link[1]);
            auto v_front_fvid = get_diagnal_vid(v_front.id, v_front_fid);
            auto v_back_fid = get_faceid(v_back.id, link[link.size() - 2]);
            auto v_back_fvid = get_diagnal_vid(v_back.id, v_back_fid);
            if (mesh.V.at(v_front_fvid).N_Fids.size() > Simplifier::minValence && mesh.V.at(v_back_fvid).N_Fids.size() > Simplifier::minValence) {
                if (can_collapse_with_feature_preserved(link, linkEids, v_front_fvid, v_back_fvid)) {
                    for (auto vid : link) {
                        auto& v = mesh.V.at(vid);
                        canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
                    }
                    collapse_with_feature_preserved(link, linkEids);
                    break;
                }
            }
        } else if (Simplifier::SPLIT && v_front.N_Fids.size() == 5 && v_back.N_Fids.size() == 5) {
            //ofs << 2 << std::endl;
            auto v_front_eid = baseComplex.separatedEdgeIdsLink.at(id).front();
            auto v_back_eid = baseComplex.separatedEdgeIdsLink.at(id).back();
            auto& v_front_fids = mesh.E.at(v_front_eid).N_Fids;
            auto& v_back_fids = mesh.E.at(v_back_eid).N_Fids;

            auto v_front_fvid = get_diagnal_vid(v_front.id, v_front_fids);
            auto v_back_fvid = get_diagnal_vid(v_back.id, v_back_fids);
            bool condition = false;
            if (mesh.V.at(v_front_fvid).N_Fids.size() < Simplifier::maxValence && mesh.V.at(v_back_fvid).N_Fids.size() < Simplifier::maxValence) {
                if (split_with_feature_preserved(link, linkEids, v_front_fvid, v_back_fvid)) {
                    for (auto vid : link) {
                        auto& v = mesh.V.at(vid);
                        canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
                    }
                    break;
                }
            }
        } else {
            ; //ofs << 3 << std::endl;
        }
    }
}

void Simplifier::half_separatrix_collapse(BaseComplexQuad& baseComplex, std::set<size_t>& canceledFids) {
	struct collapsableHalfSeparatrix {
		std::vector<size_t> target;
		std::vector<std::vector<size_t>> collapse;
	};
	std::vector<collapsableHalfSeparatrix> halfSeparatrices;
    size_t id = -1;
	for (auto link : baseComplex.separatedVertexIdsLink) {
        const auto& linkEids = baseComplex.separatedEdgeIdsLink.at(++id);
        auto& v_front = mesh.V.at(link.front());
        auto& v_back = mesh.V.at(link.back());
        if (!(v_front.isBoundary ^ v_back.isBoundary)) continue;
        //if (v_front.isBoundary || !v_back.isBoundary) continue;
        if (v_front.N_Fids.size() == 3 && v_back.N_Fids.size() == 2) {
            bool onthesameline = true;
            for (auto nvid : v_back.N_Vids) {
                auto& nv = mesh.V.at(nvid);
                if (nv.type == CORNER) {
                    onthesameline = false;
                    break;
                }
            }
            if (!onthesameline) continue;
			// indices.push_back(v_front.id);
			// indices.push_back(v_back.id);
			// numElements += 1;
			// continue;
            auto v_front_fid = get_faceid(v_front.id, link[1]);
            auto v_front_fvid = get_diagnal_vid(v_front.id, v_front_fid);
            auto& v_front_fv = mesh.V.at(v_front_fvid);
            //if (v_front_fv.isBoundary/* && v_front_fv.N_Fids.size() <= v_front_fv.idealValence*/) continue;
            if (v_front_fvid >= mesh.V.size()) {
                MeshFileWriter writer(mesh, "error.vtk");
                writer.WriteFile();
            }
            if (v_front_fv.N_Fids.size() >= Simplifier::minValence + 1) {
                if (can_collapse_with_feature_preserved(link, linkEids, v_front_fvid)) {
                    // for (auto vid : link) {
                    //     auto& v = mesh.V.at(vid);
                    //     canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
                    // }
                    // collapse_with_feature_preserved(link, linkEids);
                    // break;
					collapsableHalfSeparatrix l;
					for (size_t i = 0; i < linkEids.size(); ++i) {
						auto vid = link[i];
						auto eid = linkEids[i];
						auto p = get_collapse_vids(vid, eid);
						l.target.push_back(vid);
						l.collapse.push_back(p);
					// // 	collapse_vids_with_feature_preserved(p, vid);
					}
					auto vid = link.back();
					auto eid = linkEids.back();
					auto p = get_collapse_vids(vid, eid);
					l.target.push_back(vid);
					l.collapse.push_back(p);

					halfSeparatrices.push_back(l);
                }
            }
        }
        else if (!v_back.isBoundary && v_back.N_Fids.size() == 3 && v_front.isBoundary && v_front.N_Fids.size() == 2) {
            std::reverse(link.begin(), link.end());
            bool onthesameline = true;
            for (auto nvid : v_front.N_Vids) {
                auto& nv = mesh.V.at(nvid);
                if (nv.type == CORNER) {
                    onthesameline = false;
                    break;
                }
            }
            if (!onthesameline) continue;
			// indices.push_back(v_front.id);
			// indices.push_back(v_back.id);
			// numElements += 1;
			// continue;
            auto v_front_fid = get_faceid(v_back.id, link[1]);
            auto v_front_fvid = get_diagnal_vid(v_back.id, v_front_fid);
            auto& v_front_fv = mesh.V.at(v_front_fvid);
            if (v_front_fv.isBoundary/* && v_front_fv.N_Fids.size() <= v_front_fv.idealValence*/) continue;
            if (v_front_fvid >= mesh.V.size()) {
                MeshFileWriter writer(mesh, "error.vtk");
                writer.WriteFile();
            }
            if (mesh.V.at(v_front_fvid).N_Fids.size() >= Simplifier::minValence + 1) {
                if (can_collapse_with_feature_preserved(link, linkEids, v_front_fvid)) {
                    // for (auto vid : link) {
                    //     auto& v = mesh.V.at(vid);
                    //     canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
                    // }
                    // collapse_with_feature_preserved(link, linkEids);
                    // break;
					collapsableHalfSeparatrix l;
					for (size_t i = 0; i < linkEids.size(); ++i) {
						auto vid = link[i];
						auto eid = linkEids[i];
						auto p = get_collapse_vids(vid, eid);
						l.target.push_back(vid);
						l.collapse.push_back(p);
					// // 	collapse_vids_with_feature_preserved(p, vid);
					}
					auto vid = link.back();
					auto eid = linkEids.back();
					auto p = get_collapse_vids(vid, eid);
					l.target.push_back(vid);
					l.collapse.push_back(p);

					halfSeparatrices.push_back(l);
                }
            }
        }
    }
	std::vector<double> ranks;
	for (auto l : halfSeparatrices) {
		double rank = 0;
		for (int i = 0; i < l.target.size(); i++) {
			auto& v1 = mesh.V.at(l.target.at(i));
			for (auto value: l.collapse.at(i)) {
				auto& v2 = mesh.V.at(value);
				rank += glm::length(glm::dvec3(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z));
			}
		}
		ranks.push_back(rank);
	}
	std::vector<size_t> targetVidsPos;
	// std::vector<double>::iterator max_index = std::max_element(ranks.begin(), ranks.end());
    // double max_rank = (double) std::distance(ranks.begin(), max_index) + 1;
	for (int i = 0; i < ranks.size(); i++) {
        std::vector<double>::iterator index = std::max_element(ranks.begin(), ranks.end());
        // std::vector<double>::iterator index = std::min_element(ranks.begin(), ranks.end());
        targetVidsPos.push_back((size_t) std::distance(ranks.begin(), index));
        // targetVidsPos.push_back(i);
        *index = -1;
        // *index = max_rank;
		// ranks.erase(ranks.begin() + (size_t) std::distance(ranks.begin(), index));
		// i = 0;
    }
	std::vector<collapsableHalfSeparatrix> finalHalfSeparatrices;
	for (int i = 0; i < targetVidsPos.size(); i++) {
		if (targetVidsPos.at(i) == -1) {
			continue;
		}
		collapsableHalfSeparatrix l = halfSeparatrices.at(targetVidsPos.at(i));
		finalHalfSeparatrices.push_back(l);
		for (int j = 0; j < halfSeparatrices.size(); j++) {
			if (j == targetVidsPos.at(i)) {
				continue;
			}
			bool disjoint = true;
			collapsableHalfSeparatrix l2 = halfSeparatrices.at(j);
			for (auto id: l.target) {
				if (std::find(l2.target.begin(), l2.target.end(), id) != l2.target.end()) {
					disjoint = false;
					break;
				}
				if (!disjoint) {
					break;
				}
				for (auto vec: l2.collapse) {
					if (std::find(vec.begin(), vec.end(), id) != vec.end()) {
						disjoint = false;
						break;
					}
				}
				if (!disjoint) {
					break;
				}
			}
			for (auto vec: l.collapse) {
				for (auto id: l2.target) {
					if (std::find(vec.begin(), vec.end(), id) != vec.end()) {
						disjoint = false;
						break;
					}
				}
				if (!disjoint) {
					break;
				}
				for (auto vec2: l2.collapse) {
					for (auto id: vec2) {
						if (std::find(vec.begin(), vec.end(), id) != vec.end()) {
							disjoint = false;
							break;
						}
					}
					if (!disjoint) {
						break;
					}
				}
				if (!disjoint) {
					break;
				}
			}
			if (!disjoint) {
				auto it = std::find(targetVidsPos.begin(), targetVidsPos.end(), j);
                if (it != targetVidsPos.end()) {
                    targetVidsPos.at(std::distance(targetVidsPos.begin(), it)) = -1;
                }
			}
		}
	}
	for (auto l: finalHalfSeparatrices) {
		for (int i = 0; i < l.target.size(); i++) {
			auto& v = mesh.V.at(l.target.at(i));
			canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
			collapse_vids_with_feature_preserved(l.collapse.at(i), l.target.at(i));
		}
		// break;
	}
}

void Simplifier::half_simplify(BaseComplexQuad& baseComplex, std::set<size_t>& canceledFids) {
    size_t id = -1;
    for (auto link : baseComplex.separatedVertexIdsLink) {
        const auto& linkEids = baseComplex.separatedEdgeIdsLink.at(++id);
        auto& v_front = mesh.V.at(link.front());
        auto& v_back = mesh.V.at(link.back());
        if (!(v_front.isBoundary ^ v_back.isBoundary)) continue;
        //if (v_front.isBoundary || !v_back.isBoundary) continue;
        if (v_front.N_Fids.size() == 3 && v_back.N_Fids.size() == 2) {
            bool onthesameline = true;
            for (auto nvid : v_back.N_Vids) {
                auto& nv = mesh.V.at(nvid);
                if (nv.type == CORNER) {
                    onthesameline = false;
                    break;
                }
            }
            if (!onthesameline) continue;
            auto v_front_fid = get_faceid(v_front.id, link[1]);
            auto v_front_fvid = get_diagnal_vid(v_front.id, v_front_fid);
            auto& v_front_fv = mesh.V.at(v_front_fvid);
            //if (v_front_fv.isBoundary/* && v_front_fv.N_Fids.size() <= v_front_fv.idealValence*/) continue;
            if (v_front_fvid >= mesh.V.size()) {
                MeshFileWriter writer(mesh, "error.vtk");
                writer.WriteFile();
            }
            if (v_front_fv.N_Fids.size() >= Simplifier::minValence + 1) {
                if (can_collapse_with_feature_preserved(link, linkEids, v_front_fvid)) {
                    for (auto vid : link) {
                        auto& v = mesh.V.at(vid);
                        canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
                    }
                    collapse_with_feature_preserved(link, linkEids);
                    break;
                }
            }
        }
        else if (!v_back.isBoundary && v_back.N_Fids.size() == 3 && v_front.isBoundary && v_front.N_Fids.size() == 2) {
            std::reverse(link.begin(), link.end());
            bool onthesameline = true;
            for (auto nvid : v_front.N_Vids) {
                auto& nv = mesh.V.at(nvid);
                if (nv.type == CORNER) {
                    onthesameline = false;
                    break;
                }
            }
            if (!onthesameline) continue;
            auto v_front_fid = get_faceid(v_back.id, link[1]);
            auto v_front_fvid = get_diagnal_vid(v_back.id, v_front_fid);
            auto& v_front_fv = mesh.V.at(v_front_fvid);
            if (v_front_fv.isBoundary/* && v_front_fv.N_Fids.size() <= v_front_fv.idealValence*/) continue;
            if (v_front_fvid >= mesh.V.size()) {
                MeshFileWriter writer(mesh, "error.vtk");
                writer.WriteFile();
            }
            if (mesh.V.at(v_front_fvid).N_Fids.size() >= Simplifier::minValence + 1) {
                if (can_collapse_with_feature_preserved(link, linkEids, v_front_fvid)) {
                    for (auto vid : link) {
                        auto& v = mesh.V.at(vid);
                        canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
                    }
                    collapse_with_feature_preserved(link, linkEids);
                    break;
                }
            }
        }
    }
}

// void Simplifier::half_simplify(BaseComplexQuad& baseComplex, std::set<size_t>& canceledFids) {
// 	std::vector<size_t> indices;
// 	size_t numElements = 0;
//     size_t id = -1;
//     for (auto link : baseComplex.separatedVertexIdsLink) {
//         const auto& linkEids = baseComplex.separatedEdgeIdsLink.at(++id);
//         auto& v_front = mesh.V.at(link.front());
//         auto& v_back = mesh.V.at(link.back());
//         if (!(v_front.isBoundary ^ v_back.isBoundary)) continue;
//         //if (v_front.isBoundary || !v_back.isBoundary) continue;
//         if (v_front.N_Fids.size() == 3 && v_back.N_Fids.size() == 2) {
//             bool onthesameline = true;
//             for (auto nvid : v_back.N_Vids) {
//                 auto& nv = mesh.V.at(nvid);
//                 if (nv.type == CORNER) {
//                     onthesameline = false;
//                     break;
//                 }
//             }
//             if (!onthesameline) continue;
// 			// indices.push_back(v_front.id);
// 			// indices.push_back(v_back.id);
// 			// numElements += 1;
// 			// continue;
//             auto v_front_fid = get_faceid(v_front.id, link[1]);
//             auto v_front_fvid = get_diagnal_vid(v_front.id, v_front_fid);
//             auto& v_front_fv = mesh.V.at(v_front_fvid);
//             //if (v_front_fv.isBoundary/* && v_front_fv.N_Fids.size() <= v_front_fv.idealValence*/) continue;
//             if (v_front_fvid >= mesh.V.size()) {
//                 MeshFileWriter writer(mesh, "error.vtk");
//                 writer.WriteFile();
//             }
//             if (v_front_fv.N_Fids.size() >= Simplifier::minValence + 1) {
//                 if (can_collapse_with_feature_preserved(link, linkEids, v_front_fvid)) {
//                     for (auto vid : link) {
//                         auto& v = mesh.V.at(vid);
//                         canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
//                     }
//                     collapse_with_feature_preserved(link, linkEids);
//                     break;
//                 }
//             }
//         }
//         else if (!v_back.isBoundary && v_back.N_Fids.size() == 3 && v_front.isBoundary && v_front.N_Fids.size() == 2) {
//             std::reverse(link.begin(), link.end());
//             bool onthesameline = true;
//             for (auto nvid : v_front.N_Vids) {
//                 auto& nv = mesh.V.at(nvid);
//                 if (nv.type == CORNER) {
//                     onthesameline = false;
//                     break;
//                 }
//             }
//             if (!onthesameline) continue;
// 			// indices.push_back(v_front.id);
// 			// indices.push_back(v_back.id);
// 			// numElements += 1;
// 			// continue;
//             auto v_front_fid = get_faceid(v_back.id, link[1]);
//             auto v_front_fvid = get_diagnal_vid(v_back.id, v_front_fid);
//             auto& v_front_fv = mesh.V.at(v_front_fvid);
//             if (v_front_fv.isBoundary/* && v_front_fv.N_Fids.size() <= v_front_fv.idealValence*/) continue;
//             if (v_front_fvid >= mesh.V.size()) {
//                 MeshFileWriter writer(mesh, "error.vtk");
//                 writer.WriteFile();
//             }
//             if (mesh.V.at(v_front_fvid).N_Fids.size() >= Simplifier::minValence + 1) {
//                 if (can_collapse_with_feature_preserved(link, linkEids, v_front_fvid)) {
//                     for (auto vid : link) {
//                         auto& v = mesh.V.at(vid);
//                         canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
//                     }
//                     collapse_with_feature_preserved(link, linkEids);
//                     break;
//                 }
//             }
//         }
//     }
// 	/*std::ofstream ofs("half_separatrix_connections.vtk");
//     ofs << "# vtk DataFile Version 3.0\n"
//         << "output.vtk\n"
//         << "ASCII\n\n"
//         << "DATASET UNSTRUCTURED_GRID\n";
//     ofs << "POINTS " << mesh.V.size() << " double\n";
    
//     for (size_t i = 0; i < mesh.V.size(); i++) {
//         ofs << mesh.V.at(i).x << " " <<  mesh.V.at(i).y << " " <<  mesh.V.at(i).z << "\n";
//     }
//     ofs << "CELLS " << numElements << " " << 3 * numElements << std::endl;
//     for (size_t i = 0; i < indices.size(); i+=2) {
//         ofs << "2 " << indices.at(i) << " " << indices.at(i+1) << std::endl;
//     }
//     ofs << "CELL_TYPES " << numElements << "\n";
//     for (size_t i = 0; i < numElements; i++) {
//         ofs << "3" << std::endl;
//     }
// 	return;*/
// }

void Simplifier::sheet_simplify(BaseComplexQuad& baseComplex, std::set<size_t>& canceledFids) {
	size_t id = 0;
	for (const auto& link : baseComplex.separatedVertexIdsLink) {
		const auto& linkEids = baseComplex.separatedEdgeIdsLink.at(id);
		auto& v_front = mesh.V.at(link.front());
		auto& v_back = mesh.V.at(link.back());
		if (v_front.isBoundary || v_back.isBoundary) {
			++id;
			continue;
		}
		if (Simplifier::GLOBAL && ((v_front.N_Fids.size() < 4 && v_back.N_Fids.size() > 4) || (v_back.N_Fids.size() < 4 && v_front.N_Fids.size() > 4))) {
			for (auto eid : linkEids)
				if (collapse_sheet(canceledFids, eid)) return;
		}
		++id;
	}
}

void Simplifier::simplify(BaseComplexQuad& baseComplex, std::set<size_t>& canceledFids) {
	for (auto t = 0; t < 2; ++t) {
		size_t id = 0;
		for (const auto& link : baseComplex.separatedVertexIdsLink) {
			const auto& linkEids = baseComplex.separatedEdgeIdsLink.at(id);
			auto& v_front = mesh.V.at(link.front());
			auto& v_back = mesh.V.at(link.back());
			if ((v_front.N_Fids.size() <= 5 && v_back.N_Fids.size() <= 5) && v_front.N_Fids.size() != v_back.N_Fids.size()) {
				; //ofs << 0 << std::endl;
			} else if (COLLAPSE && v_front.N_Fids.size() == 3 && v_back.N_Fids.size() == 3) {
				// ofs << 1 << std::endl;
				auto v_front_fid = get_faceid(v_front.id, link[1]);
				auto v_front_fvid = get_diagnal_vid(v_front.id, v_front_fid);
				auto v_back_fid = get_faceid(v_back.id, link[link.size() - 2]);
				auto v_back_fvid = get_diagnal_vid(v_back.id, v_back_fid);
				bool condition = false;
				if ((t == 0 && mesh.V.at(v_front_fvid).N_Fids.size() == 5 && mesh.V.at(v_back_fvid).N_Fids.size() == 5)
					|| (t == 1 && (mesh.V.at(v_front_fvid).N_Fids.size() > Simplifier::minValence && mesh.V.at(v_back_fvid).N_Fids.size() > Simplifier::minValence)))
					condition = true;
				if (condition) {
					if (can_collapse_with_feature_preserved(link, linkEids)) {
						for (auto vid : link) {
							auto& v = mesh.V.at(vid);
							canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
						}
						collapse_with_feature_preserved(link, linkEids);
						break;
					}
				}
			} else if (Simplifier::SPLIT && v_front.N_Fids.size() == 5 && v_back.N_Fids.size() == 5) {
				//ofs << 2 << std::endl;
				auto v_front_eid = baseComplex.separatedEdgeIdsLink.at(id).front();
				auto v_back_eid = baseComplex.separatedEdgeIdsLink.at(id).back();
				auto& v_front_fids = mesh.E.at(v_front_eid).N_Fids;
				auto& v_back_fids = mesh.E.at(v_back_eid).N_Fids;

				auto v_front_fvid = get_diagnal_vid(v_front.id, v_front_fids);
				auto v_back_fvid = get_diagnal_vid(v_back.id, v_back_fids);
				bool condition = false;
				if ((t == 0 && mesh.V.at(v_front_fvid).N_Fids.size() == 3 && mesh.V.at(v_back_fvid).N_Fids.size() == 3)
					|| (t == 1 && (mesh.V.at(v_front_fvid).N_Fids.size() < Simplifier::maxValence && mesh.V.at(v_back_fvid).N_Fids.size() < Simplifier::maxValence)))
					condition = true;
				if (condition) {
					if (!split_with_feature_preserved(link, linkEids, v_front_fvid, v_back_fvid)) {
						++id;
						continue;
					}
					for (auto vid : link) {
						auto& v = mesh.V.at(vid);
						canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
					}
					break;
				}
			} else {
				; //ofs << 3 << std::endl;
			}
			++id;
		}
		if (!canceledFids.empty()) break;
	}
}

void Simplifier::smooth_project() {
	std::map<size_t, std::set<size_t>> origLabel_vids;
	std::map<size_t, std::set<size_t>> origPatch_vids;
	std::map<size_t, std::set<size_t>> origLabel_eids;
	std::map<size_t, std::set<size_t>> origSharpEdgeVid_NVids;

	std::map<size_t, std::set<size_t>> label_vids;
	std::map<size_t, std::set<size_t>> label_eids;
	std::map<size_t, std::set<size_t>> sharpEdgeVid_NVids;
	Vertex vertex;
	// std::cout << "origMesh.V: " << origMesh.V.size() << std::endl;

	std::vector<Vertex> centerVertices;
	centerVertices.insert(centerVertices.begin(), origMesh.V.begin(), origMesh.V.end());

	for (auto& f : origMesh.F) {
		glm::dvec3 center(0, 0, 0);
		for (auto nvid : origMesh.F.at(f.id).Vids)
			center += origMesh.V.at(nvid).xyz();
		center *= 0.25;
		vertex = center;
		// vertex.id = origMesh.V.size();
		vertex.id = centerVertices.size();
		vertex.patch_id = f.label;
		// origMesh.V.push_back(vertex);
		centerVertices.push_back(vertex);
		//origPatch_vids[vertex.patch_id].insert(vertex.id);
	}
	for (auto& e : origMesh.E) {
		glm::dvec3 center(0, 0, 0);
		for (auto nvid : e.Vids)
			center += origMesh.V.at(nvid).xyz();
		center *= 0.5;
		vertex = center;
		// vertex.id = origMesh.V.size();
		vertex.id = centerVertices.size();

		if (!e.isSharpFeature) {
			vertex.patch_id = origMesh.V.at(e.Vids[0]).type == REGULAR ? origMesh.V.at(e.Vids[0]).patch_id : origMesh.V.at(e.Vids[1]).patch_id;
		} else {
			continue;
			vertex.label = origMesh.V.at(e.Vids[0]).label == MAXID ? origMesh.V.at(e.Vids[1]).label : origMesh.V.at(e.Vids[0]).label;
			vertex.type = FEATURE;
		}
		vertex.isBoundary = e.isBoundary;
		// origMesh.V.push_back(vertex);
		centerVertices.push_back(vertex);
		//origPatch_vids[vertex.patch_id].insert(vertex.id);
	}
	std::cout << "mesh: " << mesh.V.size() << std::endl;
	for (auto& v : centerVertices)
		if (v.label != MAXID) origLabel_vids[v.label].insert(v.id);
	for (auto& v : centerVertices) {
		if (v.label != MAXID) {
			origLabel_vids[v.label].insert(v.id);
			for (auto& nvid : v.N_Vids) {
				auto& nv = centerVertices.at(nvid);
				auto& lineVids = origLabel_vids[v.label];
				if (nv.isCorner || lineVids.find(nvid) != lineVids.end()) origSharpEdgeVid_NVids[v.id].insert(nvid);
			}
		} else if (v.label == MAXID && !v.isCorner) {
			origPatch_vids[v.patch_id].insert(v.id);
		}
	}
	std::cout << "origPatch_vids: " << origPatch_vids[0].size() << std::endl;

	for (auto& v : mesh.V)
		if (v.label != MAXID) label_vids[v.label].insert(v.id);
	for (auto& v : mesh.V) {
		if (v.label != MAXID) {
			label_vids[v.label].insert(v.id);
			for (auto& nvid : v.N_Vids) {
				auto& nv = mesh.V.at(nvid);
				auto& lineVids = label_vids[v.label];
				if (nv.isCorner || lineVids.find(nvid) != lineVids.end()) sharpEdgeVid_NVids[v.id].insert(nvid);
			}
		}
	}

	// smooth and project
	int iters = smoothIters;
	int iter = 0;
	while (iters--) {
		std::cout << "smooth iter = " << iter++ << std::endl;
		// std::cout << "sharpEdgeVid_NVids: " << sharpEdgeVid_NVids.size() << std::endl;
		for (auto& item : sharpEdgeVid_NVids) {
			auto& v = mesh.V.at(item.first);
			if (v.type != FEATURE/*v.isCorner*/) continue;
			glm::dvec3 center(0, 0, 0);
			for (auto nvid : item.second)
				center += mesh.V.at(nvid).xyz();
			center /= item.second.size();
			const auto& origLineVids = origLabel_vids[v.label];
			size_t closest_origLineVid = *origLineVids.begin();
			double closest_distance = 100000000.0;
			for (auto vid : origLineVids) {
				auto& origv = centerVertices.at(vid);
				auto distance = glm::length(origv.xyz() - center);
				if (distance < closest_distance) {
					closest_origLineVid = vid;
					closest_distance = distance;
				}
			}
			v = centerVertices.at(closest_origLineVid).xyz();
		}
		// std::cout << "mesh.V: " << mesh.V.size() << std::endl;
		for (auto& v : mesh.V) {
			if (v.type >= FEATURE) continue;
			glm::dvec3 center(0, 0, 0);
			if (iters < 3)
				for (auto nvid : v.N_Vids)
					center += mesh.V.at(nvid).xyz();
			else
				for (auto nfid : v.N_Fids)
					for (auto nvid : mesh.F.at(nfid).Vids)
						center += 0.25 * mesh.V.at(nvid).xyz();
			center /= v.N_Fids.size();
			v = center;
			// // if (iters < 3) continue;

			auto& patchVids = origPatch_vids[v.patch_id];
			size_t closest_origVid = *patchVids.begin();
			double closest_distance = 100000000.0;
			int i = 0;
			for (auto patchVid : patchVids) {
				auto& origv = centerVertices.at(patchVid);
				auto distance = glm::length(origv.xyz() - center);
				if (distance < closest_distance) {
					closest_origVid = origv.id;
					closest_distance = distance;
				}
				i += 1;
			}
			// std::cout << "num Patch Vids iterations: " << i << std::endl;

			// v = centerVertices.at(closest_origVid).xyz();
		}
	}
}

static void refineVertexInFaces(Mesh& mesh, std::vector<Vertex>& refinedV, int resolution = 3) {
    Vertex vertex;
    for (auto& f : mesh.F) {
        for (double u = 0; u < resolution; ++u)
            for (double v = 0; v < resolution; ++v) {
                auto base = 1.0 / (resolution + 1);
                auto v01 = (1.0 + u) * base * mesh.V.at(f.Vids[0]).xyz() + (resolution - u) * base * mesh.V.at(f.Vids[1]).xyz();
                auto v32 = (1.0 + u) * base * mesh.V.at(f.Vids[3]).xyz() + (resolution - u) * base * mesh.V.at(f.Vids[2]).xyz();
                
				vertex = (1.0 + v) * base * v01 + (resolution - v) * base * v32;
                // vertex.id = mesh.V.size();
                vertex.id = refinedV.size();
                vertex.patch_id = f.label;
                // mesh.V.push_back(vertex);
                refinedV.push_back(vertex);
            }
    }
}

static void refineVertexInEdges(Mesh& mesh, std::vector<Vertex>& refinedV, int resolution = 3) {
    auto& origMesh = mesh;
    Vertex vertex;
    for (auto& e : origMesh.E) {
        for (double u = 0; u < resolution; ++u) {
            auto base = 1.0 / (resolution + 1);
            vertex = (1.0 + u) * base * origMesh.V.at(e.Vids[0]).xyz() + (resolution - u) * base * origMesh.V.at(e.Vids[1]).xyz();
            // vertex.id = origMesh.V.size();
            vertex.id = refinedV.size();

            if (!e.isSharpFeature) {
                vertex.patch_id = origMesh.V.at(e.Vids[0]).type == REGULAR ? origMesh.V.at(e.Vids[0]).patch_id : origMesh.V.at(e.Vids[1]).patch_id;
                vertex.label = MAXID;
            } else {
                if (origMesh.V.at(e.Vids[0]).label == origMesh.V.at(e.Vids[1]).label)
                    vertex.label = origMesh.V.at(e.Vids[0]).label;
                else vertex.label = origMesh.V.at(e.Vids[0]).isCorner ? origMesh.V.at(e.Vids[1]).label : origMesh.V.at(e.Vids[0]).label;
                vertex.type = FEATURE;
            }
			vertex.isBoundary = e.isBoundary;
            // origMesh.V.push_back(vertex);
            refinedV.push_back(vertex);
        }
    }
}

std::vector<size_t> get_pair(const Mesh& mesh, const Vertex& vi, const Edge& e) {
    std::vector<size_t> res;
    for (auto nfid : e.N_Fids) {
        auto& nf = mesh.F.at(nfid);
        for (auto nfvid : nf.Vids)
            if (std::find(vi.N_Vids.begin(), vi.N_Vids.end(), nfvid) != vi.N_Vids.end()
                    && std::find(e.Vids.begin(), e.Vids.end(), nfvid) == e.Vids.end()) {
                res.push_back(nfvid);
                break;
            }
    }
    return res;
}

double Simplifier::laplacian_positive_cotan_weight(const Vertex& vi, const Edge& e) {
    auto vj_id = e.Vids[0] == vi.id ? e.Vids[1] : e.Vids[0];
    auto p = get_pair(mesh, vi, e);
    auto vj_next = p.front();
    auto vj_prev = p.back();

    auto pi      = vi.xyz();
    auto pj      = mesh.V.at(vj_id).xyz();
    auto pj_prev = mesh.V.at(vj_next).xyz();
    auto pj_next = mesh.V.at(vj_prev).xyz();

    auto e1 = glm::length(pi - pj);
    auto e2 = glm::length(pi - pj_prev);
    auto e3 = glm::length(pj_prev - pj);
    // NOTE: cos(alpha) = (a^2.b^2  - c^2) / (2.a.b)
    // with a, b, c the lengths of of the sides of the triangle and (a, b) forming the angle alpha.
    auto cos_alpha = fabs((e3 * e3 + e2 * e2 - e1 * e1) / (2.0 * e3 * e2));

    auto e4 = glm::length(pi - pj_next);
    auto e5 = glm::length(pj_next - pj);
    auto cos_beta = fabs((e4 * e4 + e5 * e5 - e1 * e1) / (2.0 * e4 * e5));

    // NOTE: cot(x) = cos(x)/sin(x)
    // and recall cos(x)^2 + sin(x)^2 = 1
    // then sin(x) = sqrt(1-cos(x))
    auto cotan1 = cos_alpha / sqrt(1.0f - cos_alpha * cos_alpha);
    auto cotan2 = cos_beta  / sqrt(1.0f - cos_beta  * cos_beta );

    // wij = (cot(alpha) + cot(beta))
    float wij = (cotan1 + cotan2) / 2.0f;

    if (isnan(wij)) wij = 0.0f;
    // compute the cotangent value close to 0.0f. as cotan approaches infinity close to zero we clamp higher values
    auto eps = 1e-6;
    auto cotan_max = cos(eps) / sin(eps);
    if (wij >= cotan_max) wij = cotan_max;
    return wij;
}

void Simplifier::smooth_project(int resolution) {
//	{
//        MeshFileWriter writer(origMesh, "Orig.vtk");
//        writer.WriteVertexFeatureVtk();
//	}
	std::vector<Vertex> refinedV;
	refinedV.insert(refinedV.begin(), origMesh.V.begin(), origMesh.V.end());
	refineVertexInFaces(origMesh, refinedV, resolution);
	refineVertexInEdges(origMesh, refinedV, resolution);
//	{
//        MeshFileWriter writer(origMesh, "Refine.vtk");
//        writer.WriteVertexFeatureVtk();
//	}

    std::map<size_t, std::set<size_t>> origLabel_vids;
    std::map<size_t, std::set<size_t>> origPatch_vids;
    std::map<size_t, std::set<size_t>> origLabel_eids;
    std::map<size_t, std::set<size_t>> origSharpEdgeVid_NVids;

    std::map<size_t, std::set<size_t>> label_vids;
    std::map<size_t, std::set<size_t>> label_eids;
    std::map<size_t, std::set<size_t>> sharpEdgeVid_NVids;
	// for (auto& v : origMesh.V)
	for (auto& v : refinedV)
		if (v.label != MAXID) origLabel_vids[v.label].insert(v.id);
	// for (auto& v : origMesh.V) {
	for (auto& v : refinedV) {
		if (v.label != MAXID) {
			origLabel_vids[v.label].insert(v.id);
			for (auto& nvid : v.N_Vids) {
				// auto& nv = origMesh.V.at(nvid);
				auto& nv = refinedV.at(nvid);
				auto& lineVids = origLabel_vids[v.label];
				if (nv.isCorner || lineVids.find(nvid) != lineVids.end()) origSharpEdgeVid_NVids[v.id].insert(nvid);
			}
		} else if (v.label == MAXID && !v.isCorner) {
			origPatch_vids[v.patch_id].insert(v.id);
		}
	}
	for (auto& v : mesh.V)
		if (v.label != MAXID) label_vids[v.label].insert(v.id);
	for (auto& v : mesh.V) {
		if (v.label != MAXID) {
			label_vids[v.label].insert(v.id);
			for (auto& nvid : v.N_Vids) {
				auto& nv = mesh.V.at(nvid);
				auto& lineVids = label_vids[v.label];
				if (nv.isCorner || lineVids.find(nvid) != lineVids.end()) sharpEdgeVid_NVids[v.id].insert(nvid);
			}
		}
	}
	// return;
	// smooth and project
	int iters = smoothIters;
	int iter = 0;
	while (iters--) {
		std::cout << "smooth refine iter inside loop center = " << iter++ << std::endl;
		/*int it = 0;
		while (it < 100) {
			std::vector<glm::dvec3> centers(mesh.V.size());
			for (int i = 0; i < mesh.V.size(); i++) {
				auto& v = mesh.V.at(i);
				if (v.isBoundary) continue;
				glm::dvec3 center(0, 0, 0);
				for (auto nvid: v.N_Vids)
					center += mesh.V.at(nvid).xyz();
				center /= v.N_Vids.size();
				// auto w = 0.0;
				// for (auto neid : v.N_Eids) {
				// 	auto& e = mesh.E.at(neid);
				// 	auto wij = laplacian_positive_cotan_weight(v, e);
				// 	auto nvid = e.Vids[0] == v.id ? e.Vids[1] : e.Vids[0];
				// 	center += wij * mesh.V.at(nvid).xyz();
				// 	w += wij;
				// }
				// center /= w;
				// v = center;
				centers.at(i) = center;
				// centers.push_back(center);
			}
			for (int i = 0; i < mesh.V.size(); i++) {
				auto& v = mesh.V.at(i);
				if (v.isBoundary) continue;
				// std::cout << "center: " << centers.at(i).x << " " << centers.at(i).y << " " << centers.at(i).z << std::endl;
				v = centers.at(i);
			}
			it++;
		}*/
		// continue;
		// for (auto& item : sharpEdgeVid_NVids) {
		// 	auto& v = mesh.V.at(item.first);
		// 	// if (!mesh.smoothGlobal && !v.smoothLocal) continue;
		// 	if (v.type != FEATURE/*v.isCorner*/) continue;
		// 	// glm::dvec3 center(0, 0, 0);
		// 	// for (auto nvid : item.second)
		// 	// 	center += mesh.V.at(nvid).xyz();
		// 	// center /= item.second.size();
		// 	// v = center;
		// 	const auto& origLineVids = origLabel_vids[v.label];
		// 	// size_t closest_origLineVid = *origLineVids.begin();
		// 	//std::cout << "closest_origLineVid = " << closest_origLineVid << std::endl;
		// 	double closest_distance = 100000000.0;
		// 	glm::dvec3 curr_pos = v.xyz();
		// 	for (auto vid : origLineVids) {
		// 		// auto& origv = origMesh.V.at(vid);
		// 		auto& origv = refinedV.at(vid);
		// 		auto distance = glm::length(origv.xyz() - curr_pos);
		// 		if (distance < closest_distance) {
		// 			// closest_origLineVid = vid;
		// 			closest_distance = distance;
		// 			v = origv.xyz();
		// 		}
		// 	}
		// 	// v = origMesh.V.at(closest_origLineVid).xyz();
		// 	// v = refinedV.at(closest_origLineVid).xyz();
        //     //std::cout << "closest_origLineVid = " << closest_origLineVid << std::endl;
        //     //std::cout << "----------------------------" << std::endl;
		// }
		for (auto& v : mesh.V) {
			// if (!mesh.smoothGlobal && !v.smoothLocal) continue;
			glm::dvec3 center(0, 0, 0);
//			if (iters > 10 || iters < 10)
			// {
			int n = 0;
			for (auto nvid : v.N_Vids) {
				// if (v.isBoundary && !mesh.V.at(nvid).isBoundary) continue;
				center += mesh.V.at(nvid).xyz();
				n += 1;
			}
			center /= v.N_Vids.size();
			// center /= n;
			// }
//			else {
//				for (auto nfid : v.N_Fids)
//					for (auto nvid : mesh.F.at(nfid).Vids)
//						center += 0.25 * mesh.V.at(nvid).xyz();
//				center /= v.N_Fids.size();
//			}
			// auto w = 0.0;
			// for (auto neid : v.N_Eids) {
			// 	auto& e = mesh.E.at(neid);
			// 	auto wij = laplacian_positive_cotan_weight(v, e);
			// 	auto nvid = e.Vids[0] == v.id ? e.Vids[1] : e.Vids[0];
			// 	center += wij * mesh.V.at(nvid).xyz();
			// 	w += wij;
			// }
			// center /= w;
			// if (isnan(center.x) || isnan(center.y) || isnan(center.z)) center = glm::dvec3(0, 0, 0);
			if (!v.isCorner) v = center;
			if (!v.isBoundary) continue;
			// if (iters < 3) continue;

			auto& patchVids = origPatch_vids[v.patch_id];
			// size_t closest_origVid = *patchVids.begin();
			size_t closest_origVid = refinedV.at(0).id;
			double closest_distance = 100000000.0;
			// glm::dvec3 curr_pos = v.xyz();
			// for (auto patchVid : patchVids) {
			for (auto& origv: refinedV) {
				// if (v.isBoundary != origv.isBoundary) continue;
				if (!origv.isBoundary) continue;
				// auto& origv = origMesh.V.at(patchVid);
				// auto& origv = refinedV.at(patchVid);
				auto distance = glm::length(origv.xyz() - v.xyz());
				if (distance < closest_distance) {
					closest_origVid = origv.id;
					closest_distance = distance;
					// v = origv.xyz();
				}
			}

			// v = origMesh.V.at(closest_origVid).xyz();
			v = refinedV.at(closest_origVid).xyz();
		}
	}
}
void Simplifier::smooth_project1(int resolution) {
    std::vector<Vertex> refinedV;
	refinedV.insert(refinedV.begin(), origMesh.V.begin(), origMesh.V.end());
	refineVertexInFaces(origMesh, refinedV, resolution);
	refineVertexInEdges(origMesh, refinedV, resolution);

    std::map<size_t, std::set<size_t>> origLabel_vids;
    std::map<size_t, std::set<size_t>> origPatch_vids;
    std::map<size_t, std::set<size_t>> origLabel_eids;
    std::map<size_t, std::set<size_t>> origSharpEdgeVid_NVids;

    std::map<size_t, std::set<size_t>> label_vids;
    std::map<size_t, std::set<size_t>> label_eids;
    std::map<size_t, std::set<size_t>> sharpEdgeVid_NVids;
    for (auto& v : origMesh.V)
        if (v.label != MAXID) origLabel_vids[v.label].insert(v.id);
    for (auto& v : origMesh.V) {
        if (v.label != MAXID) {
            origLabel_vids[v.label].insert(v.id);
            for (auto& nvid : v.N_Vids) {
                auto& nv = origMesh.V.at(nvid);
                auto& lineVids = origLabel_vids[v.label];
                if (nv.isCorner || lineVids.find(nvid) != lineVids.end()) origSharpEdgeVid_NVids[v.id].insert(nvid);
            }
        } else if (v.label == MAXID && !v.isCorner) {
            origPatch_vids[v.patch_id].insert(v.id);
        }
    }
    for (auto& v : mesh.V)
        if (v.label != MAXID) label_vids[v.label].insert(v.id);
    for (auto& v : mesh.V) {
        if (v.label != MAXID) {
            label_vids[v.label].insert(v.id);
            for (auto& nvid : v.N_Vids) {
                auto& nv = mesh.V.at(nvid);
                auto& lineVids = label_vids[v.label];
                if (nv.isCorner || lineVids.find(nvid) != lineVids.end()) sharpEdgeVid_NVids[v.id].insert(nvid);
            }
        }
    }

    // smooth and project
    int iters = smoothIters;
    int iter = 0;
    while (iters--) {
        std::cout << "smooth iter = " << iter++ << std::endl;
        for (auto& item : sharpEdgeVid_NVids) {
            auto& v = mesh.V.at(item.first);
            if (v.type != FEATURE/*v.isCorner*/) continue;
            glm::dvec3 center(0, 0, 0);
            for (auto nvid : item.second)
                center += mesh.V.at(nvid).xyz();
            center /= item.second.size();
            const auto& origLineVids = origLabel_vids[v.label];
            size_t closest_origLineVid = *origLineVids.begin();
            //std::cout << "closest_origLineVid = " << closest_origLineVid << std::endl;
            double closest_distance = 100000000.0;
            for (auto vid : origLineVids) {
                auto& origv = origMesh.V.at(vid);
                auto distance = glm::length(origv.xyz() - center);
                if (distance < closest_distance) {
                    closest_origLineVid = vid;
                    closest_distance = distance;
                }
            }
            v = origMesh.V.at(closest_origLineVid).xyz();
            //std::cout << "closest_origLineVid = " << closest_origLineVid << std::endl;
            //std::cout << "----------------------------" << std::endl;
        }
        for (auto& v : mesh.V) {
            if (v.type >= FEATURE) continue;
            glm::dvec3 center(0, 0, 0);
          if (iters > 10)
          {
              for (auto nvid : v.N_Vids)
                  center += mesh.V.at(nvid).xyz();
              center /= v.N_Vids.size();
          }
//          else {
//              for (auto nfid : v.N_Fids)
//                  for (auto nvid : mesh.F.at(nfid).Vids)
//                      center += 0.25 * mesh.V.at(nvid).xyz();
//              center /= v.N_Fids.size();
//          }
          else {
            auto w = 0.0;
            for (auto neid : v.N_Eids) {
                auto& e = mesh.E.at(neid);
                auto wij = laplacian_positive_cotan_weight(v, e);
                auto nvid = e.Vids[0] == v.id ? e.Vids[1] : e.Vids[0];
                center += wij * mesh.V.at(nvid).xyz();
                w += wij;
            }
            center /= w;
          }
            v = center;
            // if (iters < 3) continue;

//            auto& patchVids = origPatch_vids[v.patch_id];
//            size_t closest_origVid = *patchVids.begin();
//            double closest_distance = 100000000.0;
//            for (auto patchVid : patchVids) {
//                auto& origv = origMesh.V.at(patchVid);
//                auto distance = glm::length(origv.xyz() - center);
//                if (distance < closest_distance) {
//                    closest_origVid = origv.id;
//                    closest_distance = distance;
//                }
//            }
//
//            v = origMesh.V.at(closest_origVid).xyz();
        }
    }
}
Mesh Simplifier::RefineWithFeaturePreserved(const Mesh& hex_mesh, int clockwise) {
	const Mesh& new_mesh = hex_mesh;
	////////////////////////////////////////////////////////////////////////////
	// add vertices
	std::vector<Vertex> new_vertex(new_mesh.V.size() + new_mesh.E.size() + new_mesh.C.size());
	for (size_t i = 0; i < new_mesh.V.size(); i++) {
		auto& v = new_vertex.at(i);
		auto& newv = new_mesh.V.at(i);
		v = newv;
		v.id = i;
		v.type = newv.type;
		v.isCorner = newv.isCorner;
		v.label = newv.label;
		v.patch_id = newv.patch_id;
		v.isSpecial = newv.isSpecial;
		v.labels = newv.labels;
	}
	size_t offset = new_mesh.V.size();
	for (size_t i = 0; i < new_mesh.E.size(); i++) {
		const Edge& e = new_mesh.E.at(i);
		const Vertex& v0 = new_mesh.V.at(e.Vids[0]);
		const Vertex& v1 = new_mesh.V.at(e.Vids[1]);
		auto& v = new_vertex.at(offset + i);
		v = 0.5 * (v0.xyz() + v1.xyz());
		v.id = offset + i;
		std::set<size_t> labels = Util::get_intersect(v0.labels, v1.labels);
		std::set<size_t> patch_ids = Util::get_intersect(v0.patch_ids, v1.patch_ids);
		if (labels.size() == 1 && v0.type >= FEATURE && v1.type >= FEATURE) {
			auto label = *labels.begin();
			v.label = label;
			v.type = FEATURE;
			v.labels.insert(label);
		}

		//if (v0.type == FEATURE && v1.type == CORNER) {
		//	v.label = v0.label;
		//	v.type = FEATURE;
		//	v.labels.insert(v0.label);
		//} else if (v1.type == FEATURE && v0.type == CORNER) {
		//	v.label = v1.label;
		//	v.type = FEATURE;
		//	v.labels.insert(v1.label);
		//}

		if (patch_ids.size() == 1) {
			auto patch_id = *patch_ids.begin();
			v.patch_id = patch_id;
		}
		v.patch_ids = patch_ids;
	}
	offset = new_mesh.V.size() + new_mesh.E.size();
	size_t numOfTri = 0, numOfQuad = 0;
	for (size_t i = 0; i < new_mesh.F.size(); i++) {
		const auto& f = new_mesh.F.at(i);
		const Vertex& v0 = new_mesh.V.at(f.Vids[0]);
		const Vertex& v1 = new_mesh.V.at(f.Vids[1]);
		const Vertex& v2 = new_mesh.V.at(f.Vids[2]);
		const Vertex& v3 = new_mesh.V.at(f.Vids[3]);

		auto& v = new_vertex.at(offset + i);
		v = 0.25 * (v0.xyz() + v1.xyz() + v2.xyz() + v3.xyz());
		v.id = offset + i;
		++numOfQuad;

		std::set<size_t> patch_ids1 = Util::get_intersect(v0.patch_ids, v1.patch_ids);
		std::set<size_t> patch_ids2 = Util::get_intersect(v2.patch_ids, v3.patch_ids);
		std::set<size_t> patch_ids = Util::get_intersect(patch_ids1, patch_ids2);
		if (patch_ids.size() == 1) {
			auto patch_id = *patch_ids.begin();
			v.patch_id = patch_id;
		}
		v.patch_ids = patch_ids;
	}
	auto key_edgeId = get_key_edgeId(new_mesh);
	//auto key_faceId = get_key_faceId(new_mesh);
	Cell cell(4);
	std::vector<Cell> new_cells(numOfTri * 3 + numOfQuad * 4, cell);
	int count = 0;
	for (size_t i = 0; i < new_mesh.F.size(); i++) {
		unsigned long v_index[9];
		const auto & f = new_mesh.F.at(i);
		for (auto j = 0; j < f.Vids.size(); j++)
			v_index[j] = f.Vids.at(j);
		for (unsigned long j = 0; j < 4; j++) {
			const Edge e({ f.Vids.at(QuadEdge[j][0]), f.Vids.at(QuadEdge[j][1]) });
			unsigned long long key = (e.Vids[0] << 32) | e.Vids[1];
			if (key_edgeId.find(key) == key_edgeId.end()) std::cout << "Edge search Error !" << std::endl;
			auto e_index = key_edgeId[key];
			v_index[4 + j] = new_mesh.V.size() + e_index;
		}
		v_index[8] = new_mesh.V.size() + new_mesh.E.size() + i;
		for (int k = 0; k < 4; k++, count++)
			for (int j = 0; j < 4; j++)
				new_cells[count].Vids[j] = v_index[QuadRefine[k][j]];
	}
	Mesh mesh(new_vertex, new_cells, QUAD);
	return mesh;
}

Mesh Simplifier::RefineWithFeaturePreserved2(const Mesh& hex_mesh, int clockwise) {
	const Mesh& new_mesh = hex_mesh;
	////////////////////////////////////////////////////////////////////////////
	// add vertices
	std::vector<Vertex> new_vertex(new_mesh.V.size() + new_mesh.E.size() + new_mesh.C.size());
	for (size_t i = 0; i < new_mesh.V.size(); i++)
		new_vertex.at(i) = new_mesh.V.at(i);
	size_t offset = new_mesh.V.size();
	for (size_t i = 0; i < new_mesh.E.size(); i++) {
		const Edge& e = new_mesh.E.at(i);
		const Vertex& v0 = new_mesh.V.at(e.Vids[0]);
		const Vertex& v1 = new_mesh.V.at(e.Vids[1]);
		new_vertex.at(offset + i) = 0.5 * (v0.xyz() + v1.xyz());
	}
	offset = new_mesh.V.size() + new_mesh.E.size();
	size_t numOfTri = 0, numOfQuad = 0;
	for (size_t i = 0; i < new_mesh.C.size(); i++) {
		const auto& f = new_mesh.C.at(i);
		if (f.Vids.size() == 4) {
			const Vertex& v0 = new_mesh.V.at(f.Vids[0]);
			const Vertex& v1 = new_mesh.V.at(f.Vids[1]);
			const Vertex& v2 = new_mesh.V.at(f.Vids[2]);
			const Vertex& v3 = new_mesh.V.at(f.Vids[3]);
			new_vertex.at(offset + i) = 0.25 * (v0.xyz() + v1.xyz() + v2.xyz() + v3.xyz());
			++numOfQuad;
		} else  if (f.Vids.size() == 3) {
			const Vertex& v0 = new_mesh.V.at(f.Vids[0]);
			const Vertex& v1 = new_mesh.V.at(f.Vids[1]);
			const Vertex& v2 = new_mesh.V.at(f.Vids[2]);
			new_vertex.at(offset + i) = 0.3333333 * (v0.xyz() + v1.xyz() + v2.xyz());
			++numOfTri;
		}
	}
	auto key_edgeId = get_key_edgeId(new_mesh);
	//auto key_faceId = get_key_faceId(new_mesh);
	Cell cell(4);
	std::vector<Cell> new_cells(numOfTri * 3 + numOfQuad * 4, cell);
	int count = 0;
	for (size_t i = 0; i < new_mesh.C.size(); i++) {
		unsigned long v_index[9];
		const auto & f = new_mesh.C.at(i);
		for (auto j = 0; j < f.Vids.size(); j++)
			v_index[j] = f.Vids.at(j);
		//        if (clockwise != 0) {
		//            std::swap(v_index[1], v_index[3]);
		//            std::swap(v_index[5], v_index[7]);
		//        }
		if (f.Vids.size() == 4) {
			for (unsigned long j = 0; j < 4; j++) {
				const Edge e({ f.Vids.at(QuadEdge[j][0]), f.Vids.at(QuadEdge[j][1]) });
				unsigned long long key = (e.Vids[0] << 32) | e.Vids[1];
				if (key_edgeId.find(key) == key_edgeId.end()) std::cout << "Edge search Error !" << std::endl;
				auto e_index = key_edgeId[key];
				v_index[4 + j] = new_mesh.V.size() + e_index;
			}
			v_index[8] = new_mesh.V.size() + new_mesh.E.size() + i;
			for (int k = 0; k < 4; k++, count++)
				for (int j = 0; j < 4; j++)
					new_cells[count].Vids[j] = v_index[QuadRefine[k][j]];
		} else if (f.Vids.size() == 3) {
			for (unsigned long j = 0; j < 3; j++) {
				const Edge e({ f.Vids.at(TriEdge[j][0]), f.Vids.at(TriEdge[j][1]) });
				unsigned long long key = (e.Vids[0] << 32) | e.Vids[1];
				if (key_edgeId.find(key) == key_edgeId.end()) std::cout << "Edge search Error !" << std::endl;
				auto e_index = key_edgeId[key];
				v_index[3 + j] = new_mesh.V.size() + e_index;
			}
			v_index[6] = new_mesh.V.size() + new_mesh.E.size() + i;
			for (int k = 0; k < 3; k++, count++)
				for (int j = 0; j < 4; j++)
					new_cells[count].Vids[j] = v_index[TriRefine[k][j]];
		}
	}
	Mesh mesh(new_vertex, new_cells, QUAD);
	return mesh;
}

std::vector<size_t> Simplifier::get_ids(const std::string str) {
	std::stringstream ss(str);
	std::vector<size_t> res;
	size_t id;
	while (ss >> id) res.push_back(id);
	return res;
}

std::map<size_t, std::set<size_t>> Simplifier::get_patchid_fids() {
	std::map<size_t, std::set<size_t>> patchid_fids;
	for (auto& v : mesh.V)
		if (v.patch_id != MAXID) patchid_fids[v.patch_id].insert(v.N_Fids.begin(), v.N_Fids.end());
	return patchid_fids;
}

std::map<size_t, std::set<size_t>> Simplifier::get_patchid_vids(const std::map<size_t, std::set<size_t>>& patchid_fids) {
	std::map<size_t, std::set<size_t>> patchid_vids;
	for (auto& item : patchid_fids) {
		for (auto fid : item.second) {
			auto& f = mesh.F.at(fid);
			patchid_vids[item.first].insert(f.Vids.begin(), f.Vids.end());
		}
	}
	return patchid_vids;
}

std::set<size_t> Simplifier::get_rotate_fids() {
	std::set<size_t> res;
	auto patchid_fids = get_patchid_fids();
	auto patchid_vids = get_patchid_vids(patchid_fids);

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

std::vector<size_t> Simplifier::get_neighbor_fids(const Vertex& v, const std::set<size_t>& patch_fids) {
	std::vector<size_t> fids;
	for (auto fid : v.N_Fids)
		if (patch_fids.find(fid) != patch_fids.end()) fids.push_back(fid);
	return fids;
}

std::set<size_t> Simplifier::get_rotate_eids(const Vertex& v, const std::vector<size_t>& fids) {
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
			res.insert(e.id);
		}
	}
	return res;
}

std::set<size_t> Simplifier::get_rotate_eids_() {
	std::set<size_t> res;
	auto patchid_fids = get_patchid_fids();
	auto patchid_vids = get_patchid_vids(patchid_fids);
	auto convex_corners = get_convex_corners();
	for (auto& item : patchid_vids) {
		for (auto vid : item.second) {
			auto& v = mesh.V.at(vid);
			if (v.type == CORNER && !v.isSpecial) {
				auto fids = get_neighbor_fids(v, patchid_fids[item.first]);
				if (fids.size() < 2) continue;
				if (convex_corners.find(v.id) == convex_corners.end() && fids.size() < 4) continue;
				auto eids = get_rotate_eids(v, fids);
				res.insert(eids.begin(), eids.end());
			} else if (v.type == FEATURE) {
				auto fids = get_neighbor_fids(v, patchid_fids[item.first]);
				if (fids.size() < 3 || fids.size() == 4) continue;
				auto eids = get_rotate_eids(v, fids);
				res.insert(eids.begin(), eids.end());
			}
		}
	}

	return res;
}

double Simplifier::get_angle(const Vertex& v, const Vertex& v0, const Vertex& v1) {
	auto d0 = v0.xyz() - v.xyz();
	auto d1 = v1.xyz() - v.xyz();
	auto cosangle = glm::dot(glm::normalize(d0), glm::normalize(d1));
	return acos(cosangle) * 180.0 / PI;
}

std::vector<size_t> Simplifier::get_neighbor_vids(const Vertex& v, size_t fid) {
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

double Simplifier::get_angle(const Vertex& v, const std::vector<size_t>& fids) {
    double total_angle = 0;
    for (auto fid : fids) {
        auto vids = get_neighbor_vids(v, fid);
        total_angle += get_angle(v, mesh.V.at(vids[0]), mesh.V.at(vids[1]));
    }
    return total_angle;
}

bool Simplifier::is_convex(const Vertex& v, const std::vector<size_t>& fids) {
    return get_angle(v, fids) < 135;
}

bool Simplifier::is_concave(const Vertex& v, const std::vector<size_t>& fids) {
    return get_angle(v, fids) > 225;
}

size_t Simplifier::get_ideal_valence(const Vertex& v, const std::vector<size_t>& fids) {
    size_t angle = (size_t) get_angle(v, fids);
    auto remain = angle % 90;
    return angle / 90 + int(remain > 45);
}

std::set<size_t> Simplifier::get_rotate_eids() {
	std::set<size_t> res;
	auto patchid_fids = get_patchid_fids();
	auto patchid_vids = get_patchid_vids(patchid_fids);
	for (auto& item : patchid_vids) {
		for (auto vid : item.second) {
			auto& v = mesh.V.at(vid);
			if (v.type == CORNER && !v.isSpecial) {
				auto fids = get_neighbor_fids(v, patchid_fids[item.first]);
				if (fids.size() < 2) continue;
				if (!is_convex(v, fids) && fids.size() < 4) continue;
				auto eids = get_rotate_eids(v, fids);
				res.insert(eids.begin(), eids.end());
			} else if (v.type == FEATURE) {
				auto fids = get_neighbor_fids(v, patchid_fids[item.first]);
				if (fids.size() < 3 || fids.size() == 4) continue;
				auto eids = get_rotate_eids(v, fids);
				res.insert(eids.begin(), eids.end());
			}
		}
	}

	return res;
}

void Simplifier::rotate(const Edge& e, const Vertex& v, std::set<size_t>& canceledFids) {
	auto& f0 = mesh.F.at(e.N_Fids[0]);
	auto& f1 = mesh.F.at(e.N_Fids[1]);
	auto vid0 = e.Vids[0] == v.id ? e.Vids[1] : e.Vids[0];
	canceledFids.insert(f0.id);
	canceledFids.insert(f1.id);
	auto diag_vid0 = get_diagnal_vid(vid0, f0.id);
	auto diag_vid1 = get_diagnal_vid(v.id, f1.id);
	if (mesh.V.at(diag_vid0).type >= FEATURE) {
		diag_vid0 = get_diagnal_vid(v.id, f0.id);
		diag_vid1 = get_diagnal_vid(vid0, f1.id);
	}
	auto eids = get_boundary_eids(f0, f1, e);
	auto linkVids = get_rotate_vids(eids, diag_vid0, diag_vid1);
	insert_rotate_faces(linkVids);
}

void Simplifier::rotate(std::set<size_t>& canceledFids) {
	auto patchid_fids = get_patchid_fids();
	auto patchid_vids = get_patchid_vids(patchid_fids);
	for (auto& item : patchid_vids) {
		for (auto vid : item.second) {
			auto& v = mesh.V.at(vid);
			if ((v.type == CORNER || v.isCorner) && !v.isSpecial) {
			    v.type = CORNER;
				auto fids = get_neighbor_fids(v, patchid_fids[item.first]);
				if (fids.size() < 2) continue;
				if (!is_convex(v, fids) && fids.size() < 4) continue;
				auto eids = get_rotate_eids(v, fids);
				auto& e = mesh.E.at(*eids.begin());
				rotate(e, v, canceledFids);
				return;
			} else if (v.type == FEATURE) {
				auto fids = get_neighbor_fids(v, patchid_fids[item.first]);
				if (fids.size() < 3 || fids.size() == 4) continue;
				auto eids = get_rotate_eids(v, fids);
				auto& e = mesh.E.at(*eids.begin());
				rotate(e, v, canceledFids);
				return;
			}
		}
	}

}

//void Simplifier::rotate(std::set<size_t>& canceledFids) {
//    auto patchid_fids = get_patchid_fids();
//    auto patchid_vids = get_patchid_vids(patchid_fids);
//    for (auto& item : patchid_vids) {
//        for (auto vid : item.second) {
//            auto& v = mesh.V.at(vid);
//            if ((v.type == CORNER || v.isCorner) && !v.isSpecial) {
//                v.type = CORNER;
//                auto fids = get_neighbor_fids(v, patchid_fids[item.first]);
//                if (fids.size() < 2) continue;
//                if (!is_convex(v, fids) && fids.size() < 4) continue;
//                auto eids = get_rotate_eids(v, fids);
//                auto& e = mesh.E.at(*eids.begin());
//                rotate(e, v, canceledFids);
//                return;
//            }
//        }
//    }
//
//    for (auto& item : patchid_vids) {
//        for (auto vid : item.second) {
//            auto& v = mesh.V.at(vid);
//            if (v.type == FEATURE) {
//                auto fids = get_neighbor_fids(v, patchid_fids[item.first]);
//                if (fids.size() < 3) continue;
//                auto eids = get_rotate_eids(v, fids);
//                auto& e = mesh.E.at(*eids.begin());
//                rotate(e, v, canceledFids);
//                return;
//            }
//        }
//    }
//
//}

std::vector<size_t> Simplifier::GetLinkVids(const std::vector<size_t>& vids) {
	auto Vids = vids;
	std::vector<size_t> ringVids;
	ringVids.push_back(Vids.back());
	Vids.pop_back();
	while (ringVids.size() < vids.size()) {
		auto last_vid = ringVids.back();
		int i = 0;
		for (auto next_vid : Vids) {
			auto& last_v = mesh.V.at(last_vid);
			for (auto nvid : last_v.N_Vids)
				if (nvid == next_vid) {
					ringVids.push_back(next_vid);
					Vids.erase(Vids.begin() + i);
					break;
				}
			++i;
		}
	}
	return ringVids;
}

std::vector<size_t> Simplifier::GetLinkEVids(const std::vector<size_t>& eids) {
	std::set<size_t> vids;
	for (auto eid : eids) {
		const auto& e = mesh.E.at(eid);
		vids.insert(e.Vids.begin(), e.Vids.end());
	}
	std::vector<size_t> res;
	for (auto vid : vids) res.push_back(vid);
	return res;
}

std::vector<size_t> Simplifier::GetLinkVidsFromEids(const std::vector<size_t>& eids) {
	std::set<size_t> eids_set(eids.begin(), eids.end());
	const auto vids = GetLinkEVids(eids);
	auto Vids = vids;
	std::vector<size_t> ringVids;
	ringVids.push_back(Vids.back());
	Vids.pop_back();
	while (ringVids.size() <= vids.size()) {
		auto last_vid = ringVids.back();
		int i = 0;
		auto& last_v = mesh.V.at(last_vid);
		auto next_eid = MAXID;
		for (auto neid : last_v.N_Eids) {
			if (eids_set.find(neid) == eids_set.end()) continue;
			next_eid = neid;
			break;
		}
		if (next_eid == MAXID) {
			static bool flag = true;
			if (flag) {
				MeshFileWriter writer(mesh, "err.vtk");
				writer.WriteFile();
				flag = false;
				MeshFileWriter writer_(mesh, "errVids.vtk");
				writer_.WriteVerticesVtk(vids);
			}
			std::cerr << "Err in GetLinkVidsFromEids\n";
			return {};
			break;
		}
		auto& ne = mesh.E.at(next_eid);
		auto next_vid = ne.Vids[0] == last_vid ? ne.Vids[1] : ne.Vids[0];
		ringVids.push_back(next_vid);
		eids_set.erase(next_eid);
		++i;
	}
	return ringVids;
}

std::vector<size_t> Simplifier::get_boundary_eids(const Face& f0, const Face& f1, const Edge& exclude_e) {
	std::set<size_t> eids_set(f0.Eids.begin(), f0.Eids.end());
	eids_set.insert(f1.Eids.begin(), f1.Eids.end());
	eids_set.erase(exclude_e.id);
	std::vector<size_t> eids;
	std::copy(eids_set.begin(), eids_set.end(), std::back_inserter(eids));
	return eids;
}

std::vector<size_t> Simplifier::get_rotate_vids(const std::vector<size_t>& boundary_eids, size_t diag_vid0, size_t diag_vid1) {
	auto linkVids = GetLinkVidsFromEids(boundary_eids);
	linkVids.pop_back();
	while (linkVids.front() != diag_vid0 && linkVids.front() != diag_vid1) {
		linkVids.push_back(linkVids.front());
		linkVids.erase(linkVids.begin());
	}
	return linkVids;
}

void Simplifier::insert_rotate_faces(const std::vector<size_t>& rotate_vids) {
	std::vector<size_t> vids0, vids1;
	std::copy(rotate_vids.begin(), rotate_vids.begin() + 4, std::back_inserter(vids0));
	std::copy(rotate_vids.begin() + 3, rotate_vids.end(), std::back_inserter(vids1));
	vids1.push_back(rotate_vids.front());
	Face newf;
	newf.id = mesh.F.size();
	newf.Vids = vids0;
	mesh.F.push_back(newf);
	newf.id = mesh.F.size();
	newf.Vids = vids1;
	mesh.F.push_back(newf);
}

bool Simplifier::can_rotate1(const Edge& e) {
	auto& v0 = mesh.V.at(e.Vids[0]);
	auto& v1 = mesh.V.at(e.Vids[1]);
	if (v0.type == CORNER && v1.type == CORNER) return false;
	if (v0.type == CORNER && v1.type == FEATURE) return false;
	if (v1.type == CORNER && v0.type == FEATURE) return false;
	if (v1.type == FEATURE && v0.type == FEATURE && v0.label == v1.label) return false;
	return true;
}

bool Simplifier::can_rotate(const Edge& e) {
	auto& v0 = mesh.V.at(e.Vids[0]);
	auto& v1 = mesh.V.at(e.Vids[1]);
	if (v0.type == CORNER && v1.type == CORNER) return false;
	if (v0.type == CORNER && v1.type == FEATURE && v0.labels.find(v1.label) != v0.labels.end()) return false;
	if (v1.type == CORNER && v0.type == FEATURE && v1.labels.find(v0.label) != v1.labels.end()) return false;
	if (v1.type == FEATURE && v0.type == FEATURE && v0.label == v1.label) return false;
	return true;
}

void Simplifier::rotate(const Face& f0, const Face& f1, const Edge& e, std::set<size_t>& canceledFids) {
	canceledFids.insert(f0.id);
	canceledFids.insert(f1.id);
	auto diag_vid0 = get_diagnal_vid(e.Vids[0], f0.id);
	auto diag_vid1 = get_diagnal_vid(e.Vids[1], f1.id);
	auto eids = get_boundary_eids(f0, f1, e);
	auto linkVids = get_rotate_vids(eids, diag_vid0, diag_vid1);
	insert_rotate_faces(linkVids);
}

bool Simplifier::rotate(const Face& f0, const Face& f1, std::set<size_t>& canceledFids) {
	std::set<size_t> f0eids(f0.Eids.begin(), f0.Eids.end());
	std::set<size_t> f1eids(f1.Eids.begin(), f1.Eids.end());
	auto intersect_eids = Util::get_intersect(f0eids, f1eids);
	if (intersect_eids.size() != 1) {
		// std::cerr << "Err in auto eids = Util::get_intersect(f0eids, f1eids);\n";
		return false;
	}
	auto& e = mesh.E.at(*intersect_eids.begin());
	if (!can_rotate(e)) return false;
	rotate(f0, f1, e, canceledFids);
	return true;
}

bool Simplifier::rotate1(const std::vector<size_t>& fids, std::set<size_t>& canceledFids) {
	auto& f0 = mesh.F.at(fids.at(0));
	auto& f1 = mesh.F.at(fids.at(1));
	if (!rotate(f0, f1, canceledFids) && fids.size() == 3) {
		auto& f2 = mesh.F.at(fids.at(2));
		if (!rotate(f0, f2, canceledFids)) return false;
	} else if (fids.size() == 4) {
		auto& f2 = mesh.F.at(fids.at(2));
		auto& f3 = mesh.F.at(fids.at(3));
		if (!rotate(f0, f2, canceledFids) &&
			!rotate(f0, f3, canceledFids) &&
			!rotate(f1, f2, canceledFids)) return false;
	}
	return true;
}

bool Simplifier::rotate(const std::vector<size_t>& fids, std::set<size_t>& canceledFids) {
	auto pairs = Util::combine(fids.size(), 2);
	for (auto& p : pairs) {
		auto& f0 = mesh.F.at(fids.at(p[0]));
		auto& f1 = mesh.F.at(fids.at(p[1]));
		if (!rotate(f0, f1, canceledFids)) continue;
		return true;
	}
	return false;
}

void Simplifier::insert_rotate_fids(std::set<size_t>& canceledFids) {
	auto patchid_fids = get_patchid_fids();
	auto patchid_vids = get_patchid_vids(patchid_fids);

	for (auto& item : patchid_vids) {
		for (auto vid : item.second) {
			auto& v = mesh.V.at(vid);
			if (v.type == CORNER && !v.isSpecial) {
				auto fids = get_neighbor_fids(v, patchid_fids[item.first]);
				if (fids.size() < 2) continue;
				if (!rotate(fids, canceledFids)) continue;
				return;
			} else if (v.type == FEATURE) {
				auto fids = get_neighbor_fids(v, patchid_fids[item.first]);
				if (fids.size() < 3) continue;
				if (!rotate(fids, canceledFids)) continue;
				return;
			}
		}
	}
}

void Simplifier::remove_doublet(std::set<size_t>& canceledFids) {
	for (auto& v : mesh.V) {
		if (v.N_Fids.size() != 2 || v.isBoundary) continue;
		auto& v0 = mesh.V.at(v.N_Vids[0]);
		auto& v1 = mesh.V.at(v.N_Vids[1]);
		//if (v0.type >= FEATURE && v1.type >= FEATURE && v.type < FEATURE) continue;
		//if (v0.type > FEATURE && v1.type >= FEATURE && v.type < FEATURE) continue;
		//if (v0.type >= FEATURE && v1.type > FEATURE && v.type < FEATURE) continue;
		//if (v0.type == FEATURE && v1.type == FEATURE && v.type < FEATURE) continue;
		auto& f0 = mesh.F.at(v.N_Fids.front());
		auto& f1 = mesh.F.at(v.N_Fids.back());
		std::set<size_t> eids_set(f0.Eids.begin(), f0.Eids.end());
		eids_set.insert(f1.Eids.begin(), f1.Eids.end());
		for (auto eid : v.N_Eids)
			eids_set.erase(eid);
		std::vector<size_t> eids;
		std::copy(eids_set.begin(), eids_set.end(), std::back_inserter(eids));
		auto linkVids = GetLinkVidsFromEids(eids);
		if (linkVids.empty()) continue;
		linkVids.pop_back();
		Face newf;
		newf.id = mesh.F.size();
		newf.Vids = linkVids;
		mesh.F.push_back(newf);
		canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
		return;
	}

	for (auto& v : mesh.V) {
		if (v.N_Fids.size() != 1 || !v.isBoundary) continue;
		auto& v0 = mesh.V.at(v.N_Vids[0]);
		auto& v1 = mesh.V.at(v.N_Vids[1]);
		//if (v0.type != FEATURE || v1.type != FEATURE || v0.label != v.label || v1.label != v.label) continue;
		auto target_vid = MAXID;
		if (v0.type == FEATURE && v1.type == FEATURE && v.type == FEATURE && (v0.label && v.label && v1.label == v.label)) {
			target_vid = v.id;
		} else if (v0.type == FEATURE && v1.type == FEATURE && v.type == CORNER &&
			(v.labels.find(v0.label) != v.labels.end() && v.labels.find(v1.label) != v.labels.end())) {
			target_vid = v.id;
		} else if (v0.type == CORNER && v1.type == FEATURE && v.type == FEATURE &&
			(v0.labels.find(v.label) != v0.labels.end() && v1.label == v.label)) {
			target_vid = v0.id;
		} else if (v0.type == FEATURE && v1.type == CORNER && v.type == FEATURE &&
			(v1.labels.find(v.label) != v1.labels.end() && v0.label == v.label)) {
			target_vid = v1.id;
		}
		if (target_vid == MAXID) continue;
		auto vids = { v.id, v0.id, v1.id };
		for (auto vid : vids)
			for (auto nfid : mesh.V.at(vid).N_Fids) {
				auto& nf = mesh.F.at(nfid);
				for (auto& fvid : nf.Vids)
					if (fvid == vid) {
						fvid = target_vid;
						break;
					}
			}
		canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
		return;
	}
}

void Simplifier::collapse_diagnal1(std::set<size_t>& canceledFids) {
	for (auto& f : mesh.F) {
		auto& v0 = mesh.V.at(f.Vids[0]);
		auto& v1 = mesh.V.at(f.Vids[1]);
		auto& v2 = mesh.V.at(f.Vids[2]);
		auto& v3 = mesh.V.at(f.Vids[3]);
		if (v0.type == CORNER || v1.type == CORNER || v2.type == CORNER || v3.type == CORNER) continue;
		//if (v0.type == FEATURE && v0.N_Fids.size() == 4) continue;
		//if (v1.type == FEATURE && v1.N_Fids.size() == 4) continue;
		//if (v2.type == FEATURE && v2.N_Fids.size() == 4) continue;
		//if (v3.type == FEATURE && v3.N_Fids.size() == 4) continue;

		if (v0.type == FEATURE && v1.type == FEATURE && v0.label == v1.label) continue;
		if (v1.type == FEATURE && v2.type == FEATURE && v1.label == v2.label) continue;
		if (v2.type == FEATURE && v3.type == FEATURE && v2.label == v3.label) continue;
		if (v3.type == FEATURE && v0.type == FEATURE && v3.label == v0.label) continue;

		if ((v0.N_Fids.size() == 3 && v1.N_Fids.size() >= 4 && v2.N_Fids.size() == 3 && v3.N_Fids.size() >= 4) ||
			(v0.N_Fids.size() >= 4 && v1.N_Fids.size() == 3 && v2.N_Fids.size() >= 4 && v3.N_Fids.size() == 3)) {
			canceledFids.insert(f.id);
			auto target_vid = v0.N_Fids.size() == 3 ? v0.id : v1.id;

			//for (auto vid : f.Vids)
			{
				auto vid = v0.N_Fids.size() == 3 ? v2.id : v3.id;
				auto& v = mesh.V.at(vid);
				for (auto nfid : v.N_Fids) {
					auto& nf = mesh.F.at(nfid);
					for (auto& fvid : nf.Vids) {
						if (fvid == vid) {
							fvid = target_vid;
							break;
						}
					}
				}
			}
			//return;
		}
	}
}

void Simplifier::collapse_diagnal(std::set<size_t>& canceledFids) {
	for (auto& f : mesh.F) {
		auto& v0 = mesh.V.at(f.Vids[0]);
		auto& v1 = mesh.V.at(f.Vids[1]);
		auto& v2 = mesh.V.at(f.Vids[2]);
		auto& v3 = mesh.V.at(f.Vids[3]);
		if (v0.type == CORNER || v1.type == CORNER || v2.type == CORNER || v3.type == CORNER) continue;
		//if (v0.type == FEATURE && v0.N_Fids.size() == 4) continue;
		//if (v1.type == FEATURE && v1.N_Fids.size() == 4) continue;
		//if (v2.type == FEATURE && v2.N_Fids.size() == 4) continue;
		//if (v3.type == FEATURE && v3.N_Fids.size() == 4) continue;

		if (v0.type == FEATURE && v1.type == FEATURE && v0.label == v1.label) continue;
		if (v1.type == FEATURE && v2.type == FEATURE && v1.label == v2.label) continue;
		if (v2.type == FEATURE && v3.type == FEATURE && v2.label == v3.label) continue;
		if (v3.type == FEATURE && v0.type == FEATURE && v3.label == v0.label) continue;

		if ((v0.N_Fids.size() == 3 && v1.N_Fids.size() >= 5 && v2.N_Fids.size() == 3 && v3.N_Fids.size() >= 5) ||
			(v0.N_Fids.size() >= 5 && v1.N_Fids.size() == 3 && v2.N_Fids.size() >= 5 && v3.N_Fids.size() == 3)) {
			canceledFids.insert(f.id);
			auto target_vid = v0.N_Fids.size() == 3 ? v0.id : v1.id;

			//for (auto vid : f.Vids)
			{
				auto vid = v0.N_Fids.size() == 3 ? v2.id : v3.id;
				auto& v = mesh.V.at(vid);
				for (auto nfid : v.N_Fids) {
					auto& nf = mesh.F.at(nfid);
					for (auto& fvid : nf.Vids) {
						if (fvid == vid) {
							fvid = target_vid;
							break;
						}
					}
				}
			}
			//return;
		}
	}
	if (canceledFids.empty()) collapse_diagnal1(canceledFids);
}

void Simplifier::init() {
	mesh.CompressWithFeaturePreserved();
	mesh.RemoveUselessVertices();
	mesh.BuildAllConnectivities();
	mesh.ExtractBoundary();
	mesh.ExtractSingularities();
	mesh.BuildParallelE();
	mesh.unifyOrientation();
	mesh.GetQuadMeshArea();
	// mesh.SetOneRingNeighborhood();
	mesh.V.resize(mesh.V.size());
    for (auto& v : mesh.V)
        if (v.isCorner) {
            v.type = CORNER;
        }

}

bool Simplifier::hasSingularities() const {
    int numOfCorners = 0, numOfSingularities = 0;
    for (auto& v : mesh.V)
        if (v.isCorner) ++numOfCorners;
    for (auto& v : mesh.V)
        if (v.isSingularity) return ++numOfSingularities;
    return numOfSingularities > numOfCorners;
}

void Simplifier::align_feature() {
	std::set<size_t> canceledFids;
	auto count = 0;
	while (true) {
		if (canceledFids.empty() && Simplifier::REMOVE_DOUBLET) {
			remove_doublet(canceledFids);
			if (!canceledFids.empty()) {
				std::cout << "remove_doublet" << std::endl;
				update(canceledFids);
				canceledFids.clear();
				init();
				continue;
			}
		}
		rotate(canceledFids);
		if (canceledFids.empty()) break;
        std::cout << "rotate " << count++ << std::endl;
		update(canceledFids);
		canceledFids.clear();
		init();
	}
}

bool Simplifier::simplify(int& iter) {
	std::set<size_t> canceledFids;
	init();
	if (!hasSingularities()) return false;

	if (iter == 0 && featurePreserved) {
		get_feature();
		auto eids = get_rotate_eids();
		{
			MeshFileWriter writer(mesh, "rotate_eids.vtk");
			writer.WriteEdgesVtk(eids);
		}
		align_feature();
		{
			std::cout << "writing rotate.vtk " << std::endl;
			MeshFileWriter writer(mesh, "rotate.vtk");
			writer.WriteFile();
		}
	}
	if (iter > 0) {
		align_feature();
	}
	BaseComplexQuad baseComplex(mesh);
	baseComplex.ExtractSingularVandE();
	baseComplex.BuildE();
	strict_simplify(baseComplex, canceledFids);
	if (canceledFids.empty()) {
		std::cout << "loose_simplify\n";
		loose_simplify(baseComplex, canceledFids);
	}
	if (canceledFids.empty() && Simplifier::REMOVE_DOUBLET) {
		std::cout << "remove_doublet" << std::endl;
		remove_doublet(canceledFids);
	}
	if (canceledFids.empty() && Simplifier::GLOBAL) {
		update(canceledFids);
		init();
		global_simplify(canceledFids);
		std::cout << "Finished global_simplify\n";
	}	
	//if (canceledFids.empty() && Simplifier::GLOBAL) {
	//	std::cout << "sheet_simplify\n";
	//	sheet_simplify(baseComplex, canceledFids);
	//}
//	if (canceledFids.empty() && Simplifier::ROTATE) {
//		std::cout << "rotate_edges" << std::endl;
//		rotate(canceledFids);
//	}
	if (canceledFids.empty() && Simplifier::COLLAPSE_DIAGNAL) {
		std::cout << "collapse_diagnal" << std::endl;
		collapse_diagnal(canceledFids);
		//if (!canceledFids.empty()) {
		//	collapse_diagnal(canceledFids);
		//	update(canceledFids);
		//	MeshFileWriter writer(mesh, "collapse_diagnal.vtk");
		//	writer.WriteFile();
		//	return true;
		//}
	}
	//if (canceledFids.empty() && Simplifier::GLOBAL) {
	//	update(canceledFids);
	//	init();
	//	global_simplify1(canceledFids);
	//	std::cout << "Finished global_simplify1\n";
	//}
	if (canceledFids.empty()) return false;
	update(canceledFids);
	// if (iter > 100) 
	{
	    auto num = std::to_string(iter);
	    while (num.size() < 3) num.insert(num.begin(), '0');
		std::string fname = std::string("iter") + num + ".vtk";
		MeshFileWriter writer(mesh, fname.c_str());
		writer.WriteFile();
	}
	std::cout << "iter = " << iter++ << std::endl;
	return true;
}


size_t permutation[3][2] = { { 0, 1 },{ 1, 2 },{ 2, 0 } };

const double PI2 = 3.1415926 * 2;
double Simplifier::GetAngle(const Vertex& v, const Face& c) {
	size_t vid1 = 0;
	size_t vid2 = 0;
	for (int i = 0; i < 3; i++) {
		if (v.id != c.Vids[permutation[i][0]] && v.id != c.Vids[permutation[i][1]]) {
			vid1 = c.Vids[permutation[i][0]];
			vid2 = c.Vids[permutation[i][1]];
			break;
		}
	}
	glm::dvec3 v1 = mesh.V.at(vid1).xyz() - v.xyz();
	glm::dvec3 v2 = mesh.V.at(vid2).xyz() - v.xyz();

	return acos(glm::dot(v1, v2) / (glm::length(v1) * glm::length(v2)));
}

std::set<size_t> Simplifier::get_convex_corners() {
	std::set<size_t> res;
	for (auto& v : mesh.V) {
		if (!v.isCorner || v.isSpecial) continue;
		double gaussianCurvature = 0;
		for (auto fid : v.N_Fids) {
			const Face& f = mesh.F.at(fid);
			gaussianCurvature += GetAngle(v, f);
		}
		gaussianCurvature = PI2 - gaussianCurvature;
		if (gaussianCurvature > PI) res.insert(v.id);
	}
	return res;
}

void Simplifier::Collapse(size_t vid, size_t target_vid) {
    auto& v = mesh.V.at(vid);
    for (auto nfid : v.N_Fids) {
        auto& nf = mesh.F.at(nfid);
        for (auto& fvid : nf.Vids) {
            if (fvid == vid) {
                fvid = target_vid;
                break;
            }
        }
    }
}
