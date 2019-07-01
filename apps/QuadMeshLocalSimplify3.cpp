/*
 * QuadMeshLocalSimplify3.cpp
 *
 *  Created on: Dec 27, 2018
 *      Author: cotrik
 */
#include "QuadMeshLocalSimplify3.h"
int maxValence = 5;
int minValence = 3;
int smoothIters = 20;
bool featurePreserved = true;
double angle = 160;
const double PI = 3.1415926535;
Mesh origMesh;	
std::vector<size_t> userCorners;
std::vector<size_t> canceledCorners;
// std::vector<size_t> userCorners = { 296, 312, 441, 426 }; // fandisk
// std::vector<size_t> canceledCorners = { 310, 314, 443 };  // fandisk
bool COLLAPSE = false;
bool SPLIT = true;
bool conformal = true;
bool global = true;
bool ROTATE = true;
bool COLLAPSE_DIAGNAL = true;
bool REMOVE_DOUBLET = true;
std::string get_facekey(const Face& f) {
	std::set<size_t> vids(f.Vids.begin(), f.Vids.end());
	std::string s;
	for (auto vid : vids)
		s += std::to_string(vid) + "@";
	return s;
}

std::set<size_t> get_canceledEdgeIds(const BaseComplexSheetQuad& baseComplexSheets, std::map<size_t, size_t>& canceledFaceIds, size_t sheetId) {
	std::set<size_t> canceledEdgeIds;
	for (auto baseComplexEdgeId : baseComplexSheets.sheets_componentEdgeIds[sheetId]) {
		for (auto edgeId : baseComplexSheets.baseComplex.componentE[baseComplexEdgeId].eids_link) {
			auto& e = baseComplexSheets.baseComplex.mesh.E[edgeId];
			for (auto n_fid : e.N_Fids)
				++canceledFaceIds[n_fid];
			canceledEdgeIds.insert(edgeId);
		}
	}
	return canceledEdgeIds;
}

bool can_collapse(const BaseComplexSheetQuad& baseComplexSheets, std::map<size_t, size_t>& canceledFaceIds, size_t sheetId) {
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

void collapse(Mesh& mesh, std::unordered_map<size_t, size_t>& key_edgeId, std::unordered_map<std::string, size_t>& key_faceId,
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

bool can_collapse_with_feature_preserved(const BaseComplexSheetQuad& baseComplexSheets, std::map<size_t, size_t>& canceledFaceIds, size_t sheetId) {
	const auto& baseComplexFIds = baseComplexSheets.sheets_componentFaceIds[sheetId];
	for (auto baseComplexFId : baseComplexFIds) {
		int count = 0;
		for (auto baseComplexVId : baseComplexSheets.baseComplex.componentF[baseComplexFId].Vids) {
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

void collapse_with_feature_preserved(Mesh& mesh, std::unordered_map<size_t, size_t>& key_edgeId, std::unordered_map<std::string, size_t>& key_faceId,
	std::map<size_t, size_t>& canceledFaceIds, std::set<size_t>& canceledEdgeIds) {
	for (auto& item : canceledFaceIds) {
		if (item.second >= 4) {
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
				if (v.type == FEATURE && featureType == 0) centerVid = vid;
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


bool can_collapse_vids_with_feature_preserved(Mesh& mesh, const std::vector<size_t>& vids) {
	const auto& v0 = mesh.V.at(vids[0]);
	const auto& v1 = mesh.V.at(vids[1]);

	int count = 0;
	if (v0.isCorner) ++count;
	if (v1.isCorner) ++count;

	std::set<size_t> labels;
	if (!v0.isCorner && v0.label != MAXID) labels.insert(v0.label);
	if (!v1.isCorner && v1.label != MAXID) labels.insert(v1.label);
	if (v0.N_Fids.size() + v1.N_Fids.size() - 4 > maxValence) return false;
	if (v0.N_Fids.size() + v1.N_Fids.size() - 4 < minValence) return false;
	if (count > 1) return false;
	if (labels.size() == 1 && v0.type == CORNER && v1.type == FEATURE) return false;
	if (labels.size() == 1 && v0.type == FEATURE && v1.type == CORNER) return false;
	if (labels.size() > 2) return false;
	if (labels.size() == 2 && (v0.isCorner || v1.isCorner)) return false;
	if ((v0.type == FEATURE || v1.type == FEATURE) && v0.N_Fids.size() + v1.N_Fids.size() - 4 <= minValence) return false;
	return true;
}

bool can_collapse_vids_with_feature_preserved(Mesh& mesh, const std::set<size_t>& eids) {
	for (auto eid : eids) {
		auto& e = mesh.E.at(eid);
		if (!can_collapse_vids_with_feature_preserved(mesh, e.Vids)) return false;
	}
	return true;
}

void global_simplify(Mesh& mesh, std::set<size_t>& canceledFids) {
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

		std::map<size_t, size_t> canceledFaceIds;
		std::set<size_t> canceledEdgeIds = get_canceledEdgeIds(baseComplexSheets, canceledFaceIds, sheetId);

		if (!can_collapse_with_feature_preserved(baseComplexSheets, canceledFaceIds, sheetId)) continue;
		if (!can_collapse_vids_with_feature_preserved(mesh, canceledEdgeIds)) continue;
		{

			std::cout << "collapse sheet " << sheetId << "\n";
			collapse_with_feature_preserved(mesh, key_edgeId, key_faceId, canceledFaceIds, canceledEdgeIds);
			for (auto& item : canceledFaceIds)
				canceledFids.insert(item.first);
			break;
		}
	}
}

std::unordered_map<size_t, size_t> get_key_edgeId(const Mesh& mesh) {
	std::unordered_map<size_t, size_t> key_edgeId;
	for (size_t i = 0; i < mesh.E.size(); ++i) {
		const auto& e = mesh.E.at(i);
		key_edgeId[(e.Vids[0] << 32) | e.Vids[1]] = i;
		key_edgeId[(e.Vids[1] << 32) | e.Vids[0]] = i;
	}
	return key_edgeId;
}

std::unordered_map<std::string, size_t> get_key_faceId(const Mesh& mesh) {
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

Mesh Refine(const Mesh& hex_mesh, int clockwise) {
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
void get_parallel_edgeids(const Mesh& mesh, size_t start_edge_id, size_t start_face_id,
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
    get_parallel_edgeids(mesh, next_edge_id, next_face_id, parallel_edgeids, parallel_faceids);
}

size_t get_faceid(const Mesh& mesh, size_t vid, size_t exclude_vid) {
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

size_t get_diagnal_vid(const Mesh& mesh, size_t vid, size_t fid) {
    auto& f = mesh.F.at(fid);
    for (int i = 0; i < 4; ++i) {
        if (f.Vids[i] == vid) return f.Vids.at((i + 2) % 4);
    }
    std::cerr << "ERROR get_diagnal_vid\n";
    return MAXID;
}

size_t get_diagnal_vid(const Mesh& mesh, size_t vid, const std::vector<size_t>& fids) {
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
		if (!found) return get_diagnal_vid(mesh, vid, fid);
	}
	std::cerr << "ERROR get_diagnal_vid\n";
	return MAXID;
}

std::vector<size_t> get_collapse_vids(Mesh& mesh, size_t vid, size_t eid) {
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

std::vector<size_t> get_split_vids(Mesh& mesh, size_t vid, size_t eid) {
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

bool can_collapse_vids(Mesh& mesh, const std::vector<size_t>& vids, size_t target_vid) {
    return ((mesh.V.at(vids[0]).N_Fids.size() + mesh.V.at(vids[1]).N_Fids.size() < 10) &&
            (mesh.V.at(vids[0]).N_Fids.size() + mesh.V.at(vids[1]).N_Fids.size() > 6));
}

void collapse_vids(Mesh& mesh, const std::vector<size_t>& vids, size_t target_vid) {
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

bool can_collapse(Mesh& mesh, const std::vector<size_t>& linkVids, const std::vector<size_t>& linkEids) {
    std::set<size_t> vids;
    for (size_t i = 0; i < linkEids.size(); ++i) {
        auto vid = linkVids[i];
        auto eid = linkEids[i];
        auto p = get_collapse_vids(mesh, vid, eid);
        if (!can_collapse_vids(mesh, p, vid)) return false;
        vids.insert(p.begin(), p.end());
    }
    auto vid = linkVids.back();
    auto eid = linkEids.back();
    auto p = get_collapse_vids(mesh, vid, eid);
    if (!can_collapse_vids(mesh, p, vid)) return false;
    vids.insert(p.begin(), p.end());

    if (vids.size() < 2 * linkVids.size()) return false; // tangent;
    return true;
}

void collapse(Mesh& mesh, const std::vector<size_t>& linkVids, const std::vector<size_t>& linkEids) {
    for (size_t i = 0; i < linkEids.size(); ++i) {
        auto vid = linkVids[i];
        auto eid = linkEids[i];
        auto p = get_collapse_vids(mesh, vid, eid);
        collapse_vids(mesh, p, vid);
    }
    auto vid = linkVids.back();
    auto eid = linkEids.back();
    auto p = get_collapse_vids(mesh, vid, eid);
    collapse_vids(mesh, p, vid);
}

bool can_collapse_vids_with_feature_preserved(Mesh& mesh, const std::vector<size_t>& vids, size_t target_vid) {
	const auto& v0 = mesh.V.at(vids[0]);
	const auto& v1 = mesh.V.at(vids[1]);
	const auto& v = mesh.V.at(target_vid);

    int count  = 0;
    if (v0.isCorner) ++count;
    if (v1.isCorner) ++count;
    if (v.isCorner) ++count;

	std::set<size_t> labels;
	if (!v0.isCorner && v0.label != MAXID) labels.insert(v0.label);
	if (!v1.isCorner && v1.label != MAXID) labels.insert(v1.label);
	if (!v.isCorner && v.label != MAXID) labels.insert(v.label);	
	if (v0.N_Fids.size() + v1.N_Fids.size() - 4 > maxValence) return false;
	if (v0.N_Fids.size() + v1.N_Fids.size() - 4 < minValence) return false;
	if (count > 1) return false;
	if (labels.size() == 1 && v.type == REGULAR && v0.type == CORNER && v1.type == FEATURE) return false;
	if (labels.size() == 1 && v.type == REGULAR && v0.type == FEATURE && v1.type == CORNER) return false;
	if (labels.size() > 2) return false;
	if (v.isSpecial) return false;
	if (labels.size() == 2 && v.isSpecial) return false;
	if (labels.size() == 2 && !v.isCorner) return false;
	if (labels.size() == 2 && (v0.isCorner || v1.isCorner)) return false;
	if (labels.size() == 2 && (/*v.isCorner || */v.labels.find(v0.label) == v.labels.end() || v.labels.find(v1.label) == v.labels.end())) return false;
	if ((v0.type == FEATURE || v1.type == FEATURE || v.type == FEATURE) && v0.N_Fids.size() + v1.N_Fids.size() - 4 <= minValence) return false;
	return true;
}

void collapse_vids_with_feature_preserved(Mesh& mesh, const std::vector<size_t>& vids, size_t target_vid) {
	auto& v0 = mesh.V.at(vids[0]);
	auto& v1 = mesh.V.at(vids[1]);
	auto& v = mesh.V.at(target_vid);
	if (v0.type == MAXID) v0.type = REGULAR;
	if (v1.type == MAXID) v1.type = REGULAR;
	if (v.type == MAXID) v.type = REGULAR;
	if (/*v0.type <= FEATURE && v1.type <= FEATURE && */v.type == CORNER) {
		// v.type = CORNER;
		v0.type = REGULAR;
		v0.label = MAXID;
		v1.type = REGULAR;
		v1.label = MAXID;
	} else if (v0.type == CORNER/* && v1.type <= FEATURE && v.type <= FEATURE*/) {
		target_vid = vids[0];
		//mesh.V.at(target_vid).type = CORNER; 
		v.type = REGULAR;
		v.label = MAXID;
		v1.type = REGULAR;
		v1.label = MAXID;
	} else if (/*v0.type <= FEATURE && */v1.type == CORNER/* && v.type <= FEATURE*/) {
		target_vid = vids[1];
		// mesh.V.at(target_vid).type = CORNER; 
		v.type = REGULAR;
		v.label = MAXID;
		v0.type = REGULAR;
		v0.label = MAXID;
	} else if (v0.type == REGULAR && v1.type == REGULAR && v.type == REGULAR) {
		//v.type = REGULAR; 
	} else if (v0.type == REGULAR && v1.type == REGULAR && v.type == REGULAR) { 
		//v.type = REGULAR; 
	} else if (v0.type == REGULAR && v1.type == REGULAR && v.type == FEATURE) { 
		//v.type = FEATURE; 
	} else if (v0.type == REGULAR && v1.type == FEATURE && v.type == FEATURE) { 
		//v.type = FEATURE; 
		v1.type = REGULAR;
		v1.label = MAXID;
	} else if (v0.type == FEATURE && v1.type == REGULAR && v.type == FEATURE) { 
		//v.type = FEATURE; 
		v0.type = REGULAR;
		v0.label = MAXID;
	} else if (v0.type == REGULAR && v1.type == FEATURE && v.type == REGULAR) { 
		target_vid = vids[1]; 
		// mesh.V.at(target_vid).type = FEATURE; 
	} else if (v0.type == FEATURE && v1.type == REGULAR && v.type == REGULAR) { 
		target_vid = vids[0]; 
		// mesh.V.at(target_vid).type = FEATURE; 
	} else if (v0.type <= FEATURE && v1.type <= FEATURE && v.type == FEATURE) { 
		// v.type = FEATURE; 
		v0.type = REGULAR;
		v0.label = MAXID;
		v1.type = REGULAR;
		v1.label = MAXID;
	} 
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

std::set<size_t> get_regionVids(const Mesh& mesh, const std::vector<size_t>& linkVids) {
	std::set<size_t> res;
	for (auto vid : linkVids)
		for (auto& nvid : mesh.V.at(vid).N_Vids)
			res.insert(nvid);
	return res;
}

bool can_collapse_with_feature_preserved(Mesh& mesh, const std::vector<size_t>& linkVids, const std::vector<size_t>& linkEids) {
    std::set<size_t> vids;
    for (size_t i = 0; i < linkEids.size(); ++i) {
        auto vid = linkVids[i];
        auto eid = linkEids[i];
        auto p = get_collapse_vids(mesh, vid, eid);
        if (!can_collapse_vids_with_feature_preserved(mesh, p, vid)) return false;
        vids.insert(p.begin(), p.end());
    }
    auto vid = linkVids.back();
    auto eid = linkEids.back();
    auto p = get_collapse_vids(mesh, vid, eid);
    if (!can_collapse_vids_with_feature_preserved(mesh, p, vid)) return false;
    vids.insert(p.begin(), p.end());

    if (vids.size() < 2 * linkVids.size()) return false; // tangent;
    return true;
}

bool can_collapse_with_feature_preserved(Mesh& mesh, const std::vector<size_t>& linkVids, const std::vector<size_t>& linkEids,
	size_t v_front_fvid, size_t v_back_fvid) {
	if (!can_collapse_with_feature_preserved(mesh, linkVids, linkEids)) return false;
	if (conformal) {
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

void collapse_with_feature_preserved(Mesh& mesh, const std::vector<size_t>& linkVids, const std::vector<size_t>& linkEids) {
    for (size_t i = 0; i < linkEids.size(); ++i) {
        auto vid = linkVids[i];
        auto eid = linkEids[i];
        auto p = get_collapse_vids(mesh, vid, eid);
        collapse_vids_with_feature_preserved(mesh, p, vid);
    }
    auto vid = linkVids.back();
    auto eid = linkEids.back();
    auto p = get_collapse_vids(mesh, vid, eid);
    collapse_vids_with_feature_preserved(mesh, p, vid);
}

void WriteSharpEdgesVtk(const char* filename, const Mesh& m_mesh, const std::vector<FeatureLine>& featureLines) {
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
void get_feature(Mesh& mesh) {
	// cos10 = 0.984807753; cos15 = 0.965925826; cos20 = 0.939692621; cos25 = 0.906307787; cos30 = 0.866025404;
	const double coslocalangle = cos((180.0 - angle) * PI / 180.0);//cos(angle);
	std::cout << "coslocalangle = " << coslocalangle << std::endl;
	mesh.SetCosAngleThreshold(coslocalangle);
	mesh.GetNormalOfSurfaceFaces();
	mesh.GetNormalOfSurfaceVertices();
	Patches patches(mesh);
	patches.SetGlobalCosAngle(coslocalangle);
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
	mesh.LabelSharpEdges(true);
	{
		mesh.LabelSharpEdges(true);
		// for (auto& e : mesh.E) e.isSharpFeature = copy[e.id];
		std::vector<FeatureLine> featureLines(mesh.numOfSharpEdges, FeatureLine(mesh));
		for (size_t i = 0; i < mesh.numOfSharpEdges; i++)
			featureLines.at(i).Extract(i);
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

	origMesh = mesh;
}

bool remove_one_vertex_with_sharp_feature_preserved(Mesh& mesh, std::set<size_t>& canceledFids) {
	for (auto& v : mesh.V) {
		if (v.N_Fids.size() != 2 || v.type >= FEATURE) continue;
		auto& f = mesh.F.at(v.N_Fids.back());
		std::set<size_t> fvids(f.Vids.begin(), f.Vids.end());
		fvids.erase(v.id);
		fvids.erase(v.N_Vids.front());
		fvids.erase(v.N_Vids.back());
		auto targetvid = *fvids.begin();
		v.id = targetvid;
		canceledFids.insert(v.N_Fids.back());
		break;
	}
	return !canceledFids.empty();
}

void update(Mesh& mesh, std::set<size_t>& canceledFids) {
    std::vector<size_t> FaceIds;
    FaceIds.reserve(mesh.F.size());
    for (auto& f : mesh.F)
        if (canceledFids.find(f.id) == canceledFids.end())
            FaceIds.push_back(f.id);
    std::vector<Vertex> newV(mesh.V.size());
    std::vector<Face> newF(FaceIds.size());
    std::vector<Cell> newC(FaceIds.size());
    for (size_t i = 0; i < mesh.V.size(); ++i) {
        newV.at(i).id = i;
        newV.at(i) = mesh.V.at(i).xyz();
        newV.at(i).type = mesh.V.at(i).type;
        newV.at(i).isCorner = mesh.V.at(i).isCorner;
		newV.at(i).label = mesh.V.at(i).label;
		newV.at(i).patch_id = mesh.V.at(i).patch_id;
		newV.at(i).isSpecial = mesh.V.at(i).isSpecial;
		newV.at(i).labels = mesh.V.at(i).labels;
		newV.at(i).patch_ids = mesh.V.at(i).patch_ids;
    }
    for (size_t i = 0; i < FaceIds.size(); ++i) {
        newF.at(i).id = i;
        newF.at(i).Vids = mesh.F.at(FaceIds[i]).Vids;
    }
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
    mesh.F = newF;
    mesh.C = newC;
    canceledFids.clear();
}

void KillHighValenceSingularities_(Mesh& mesh) {
	size_t numOfSingularV = 0;
	for (const auto& v : mesh.V)
		if (v.N_Fids.size() == 6) ++numOfSingularV;

	while (numOfSingularV--) {
		bool found_valence6 = false;
		for (const auto& v : mesh.V) {
			if (v.N_Fids.size() != 6) continue;
			if (v.label != MAXID) continue;
			auto centerV = v;
			found_valence6 = true;
			std::set<size_t> neighborFVids;
			for (auto nfid : v.N_Fids) {
				auto& nf = mesh.F.at(nfid);
				neighborFVids.insert(nf.Vids.begin(), nf.Vids.end());
			}
			neighborFVids.erase(v.id);
			std::vector<size_t> linkVids;
			for (auto nfvid : neighborFVids)
				if (mesh.V.at(nfvid).N_Fids.size() == 4) {
					linkVids.push_back(nfvid);
					neighborFVids.erase(nfvid);
					break;
				}

			auto start_vid = linkVids.back();
			while (!neighborFVids.empty()) {
				auto& start_v = mesh.V.at(start_vid);
				for (auto nvid : start_v.N_Vids) {
					if (neighborFVids.find(nvid) != neighborFVids.end()) {
						linkVids.push_back(nvid);
						neighborFVids.erase(nvid);
						start_vid = nvid;
						break;
					}
				}
			}

			Vertex new_v;
			std::vector<size_t> ids = { 2, 3, 4, 8, 9, 10 };
			auto xyz = centerV.xyz();
			for (auto id : ids) {
				new_v = 0.5 * (xyz + mesh.V.at(linkVids[id]).xyz());
				new_v.id = mesh.V.size();
				mesh.V.push_back(new_v);
				linkVids.push_back(new_v.id);
			}
			const std::vector<std::vector<size_t>> new_Fvids = {
				{ linkVids[0], linkVids[1], linkVids[2], linkVids[12] },
			{ linkVids[2], linkVids[3], linkVids[13], linkVids[12] },
			{ linkVids[3], linkVids[4], linkVids[14], linkVids[13] },
			{ linkVids[4], linkVids[5], linkVids[6], linkVids[14] },
			{ linkVids[6], linkVids[7], linkVids[8], linkVids[15] },
			{ linkVids[8], linkVids[9], linkVids[16], linkVids[15] },
			{ linkVids[9], linkVids[10], linkVids[17], linkVids[16] },
			{ linkVids[10], linkVids[11], linkVids[0], linkVids[17] },
			{ centerV.id, linkVids[17], linkVids[0], linkVids[12] },
			{ centerV.id, linkVids[12], linkVids[13], linkVids[14] },
			{ centerV.id, linkVids[14], linkVids[6], linkVids[15] },
			{ centerV.id, linkVids[15], linkVids[16], linkVids[17] }
			};
			Face new_f;
			for (auto& vids : new_Fvids) {
				new_f.Vids = vids;
				new_f.id = mesh.F.size();
				mesh.F.push_back(new_f);
			}

			std::set<size_t> canceledFids(centerV.N_Fids.begin(), centerV.N_Fids.end());
			update(mesh, canceledFids);
			mesh.CompressWithFeaturePreserved();
			mesh.BuildAllConnectivities();
			mesh.ExtractBoundary();
			mesh.ExtractSingularities();
			break;
		}
		if (!found_valence6) break;
	}
}

bool KillHighValenceSingularities(Mesh& mesh) {
	bool found_valence6 = false;
	for (const auto& v : mesh.V) {
		if (v.N_Fids.size() != 6) continue;
		if (v.label != MAXID) continue;
		auto centerV = v;
		std::set<size_t> neighborFVids;
		for (auto nfid : v.N_Fids) {
			auto& nf = mesh.F.at(nfid);
			neighborFVids.insert(nf.Vids.begin(), nf.Vids.end());
		}
		neighborFVids.erase(v.id);
		std::vector<size_t> linkVids;
		for (auto nfvid : neighborFVids)
			if (mesh.V.at(nfvid).N_Fids.size() == 4) {
				linkVids.push_back(nfvid);
				neighborFVids.erase(nfvid);
				break;
			}

		auto start_vid = linkVids.back();
		while (!neighborFVids.empty()) {
			auto& start_v = mesh.V.at(start_vid);
			for (auto nvid : start_v.N_Vids) {
				if (neighborFVids.find(nvid) != neighborFVids.end()) {
					linkVids.push_back(nvid);
					neighborFVids.erase(nvid);
					start_vid = nvid;
					break;
				}
			}
		}
		// vector<size_t> check_ids = { 0, 2, 4, 6, 8, 10 };
		// vector<size_t> check_ids = { 0, 3, 6, 9 };
		std::vector<size_t> check_ids = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
		bool will_introduce_valence6 = false;
		for (auto id : check_ids) {
			if (mesh.V.at(linkVids[id]).N_Fids.size() > 4) {
				will_introduce_valence6 = true;
				break;
			}
		}
		if (will_introduce_valence6) continue;
		// if (mesh.V.at(linkVids[0]).N_Fids.size() > 4 || mesh.V.at(linkVids[6]).N_Fids.size() > 4) continue;

		Vertex new_v;
		std::vector<size_t> ids = { 2, 3, 4, 8, 9, 10 };
		auto xyz = centerV.xyz();
		for (auto id : ids) {
			new_v = 0.5 * (xyz + mesh.V.at(linkVids[id]).xyz());
			new_v.id = mesh.V.size();
			mesh.V.push_back(new_v);
			linkVids.push_back(new_v.id);
		}
		const std::vector<std::vector<size_t>> new_Fvids = {
			{ linkVids[0], linkVids[1], linkVids[2], linkVids[12] },
		{ linkVids[2], linkVids[3], linkVids[13], linkVids[12] },
		{ linkVids[3], linkVids[4], linkVids[14], linkVids[13] },
		{ linkVids[4], linkVids[5], linkVids[6], linkVids[14] },
		{ linkVids[6], linkVids[7], linkVids[8], linkVids[15] },
		{ linkVids[8], linkVids[9], linkVids[16], linkVids[15] },
		{ linkVids[9], linkVids[10], linkVids[17], linkVids[16] },
		{ linkVids[10], linkVids[11], linkVids[0], linkVids[17] },
		{ centerV.id, linkVids[17], linkVids[0], linkVids[12] },
		{ centerV.id, linkVids[12], linkVids[13], linkVids[14] },
		{ centerV.id, linkVids[14], linkVids[6], linkVids[15] },
		{ centerV.id, linkVids[15], linkVids[16], linkVids[17] }
		};
		Face new_f;
		for (auto& vids : new_Fvids) {
			new_f.Vids = vids;
			new_f.id = mesh.F.size();
			mesh.F.push_back(new_f);
		}

		std::set<size_t> canceledFids(centerV.N_Fids.begin(), centerV.N_Fids.end());
		update(mesh, canceledFids);
		found_valence6 = true;
		break;
	}

	return found_valence6;
}

size_t get_diagnal_vid(const Mesh& mesh, size_t vid, const std::vector<size_t>& fids, size_t end_vid) {
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
				res = get_diagnal_vid(mesh, vid, fid);
				break;
			}
		}
		if (found) return res;
	}
	std::cerr << "ERROR get_diagnal_vid\n";
	return res;
}

size_t get_id(const Mesh& mesh, size_t singular_vid, size_t target_vid1, size_t target_vid2) {
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

std::vector<size_t> get_insert_vids(Mesh& mesh, const std::vector<size_t>& linkVids, const std::vector<size_t>& linkVids1) {
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

std::vector<size_t> get_insert_fids(Mesh& mesh, const std::vector<size_t>& linkVids1, const std::vector<size_t>& linkVids2) {
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

bool is_neighbor(const Vertex& v, size_t vid) {
	for (auto nvid : v.N_Vids)
		if (nvid == vid) return false;
	return true;
}

bool can_split_with_feature_preserved(Mesh& mesh, size_t vid, size_t eid) {
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

bool can_split_with_feature_preserved(Mesh& mesh, const std::vector<size_t>& linkVids, const std::vector<size_t>& linkEids,
	size_t v_front_fvid, size_t v_back_fvid) {
	{
		auto vid = linkVids.front();
		auto eid = linkEids.front();
		if (!can_split_with_feature_preserved(mesh, vid, eid)) return false;
	}
	{
		auto vid = linkVids.back();
		auto eid = linkEids.back();
		//if (mesh.V.at(vid).type >= FEATURE) return false;
		if (!can_split_with_feature_preserved(mesh, vid, eid)) return false;
	}
	if (conformal) {
		auto& v_front_fv = mesh.V.at(v_front_fvid);
		auto& v_back_fv = mesh.V.at(v_back_fvid);
		if (v_front_fv.type == FEATURE && v_front_fv.N_Fids.size() == 4) return false;
		if (v_back_fv.type == FEATURE && v_back_fv.N_Fids.size() == 4) return false;
	}
	return true;
}
bool split_with_feature_preserved(Mesh& mesh, const std::vector<size_t>& linkVids, const std::vector<size_t>& linkEids,
	size_t v_front_fvid, size_t v_back_fvid) {
	if (!can_split_with_feature_preserved(mesh, linkVids, linkEids, v_front_fvid, v_back_fvid)) return false;

	std::vector<size_t> vids1, vids2;
	{
		auto vid = linkVids.front();
		auto eid = linkEids.front();
		auto p = get_split_vids(mesh, vid, eid);
		vids1.push_back(p.front());
		vids2.push_back(p.back());
	}
	for (size_t i = 1; i < linkEids.size(); ++i) {
		auto vid = linkVids[i];
		auto eid = linkEids[i];
		auto p = get_split_vids(mesh, vid, eid);
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
		auto p = get_split_vids(mesh, vid, eid);
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

	auto insertVids1 = get_insert_vids(mesh, linkVids, vids1);
	auto insertVids2 = get_insert_vids(mesh, linkVids, vids2);

	auto v_front_eid = linkEids.front();
	auto v_back_eid = linkEids.back();
	auto& v_front_fids = mesh.E.at(v_front_eid).N_Fids;
	auto& v_back_fids = mesh.E.at(v_back_eid).N_Fids;
	auto diag_vid1 = get_diagnal_vid(mesh, linkVids.front(), v_front_fids, vids1.front());
	auto diag_vid2 = get_diagnal_vid(mesh, linkVids.front(), v_front_fids, vids2.front());
	if (diag_vid1 == MAXID) return false;
	if (diag_vid2 == MAXID) return false;
	vids1.insert(vids1.begin(), diag_vid1);
	vids2.insert(vids2.begin(), diag_vid2);
	diag_vid1 = get_diagnal_vid(mesh, linkVids.back(), v_back_fids, vids1.back());
	diag_vid2 = get_diagnal_vid(mesh, linkVids.back(), v_back_fids, vids2.back());
	if (diag_vid1 == MAXID) return false;
	if (diag_vid2 == MAXID) return false;
	vids1.push_back(diag_vid1);
	vids2.push_back(diag_vid2);
	{
		auto id = get_id(mesh, linkVids.front(), vids1.front(), v_front_fvid);
		vids1.insert(vids1.begin(), id);
		id = get_id(mesh, linkVids.front(), vids2.front(), v_front_fvid);
		vids2.insert(vids2.begin(), id);

		id = get_id(mesh, linkVids.back(), vids1.back(), v_back_fvid);
		vids1.push_back(id);
		id = get_id(mesh, linkVids.back(), vids2.back(), v_back_fvid);
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
	auto fids = get_insert_fids(mesh, linkVids1, linkVids2);
	fids = get_insert_fids(mesh, linkVids2, linkVids3);
	fids = get_insert_fids(mesh, linkVids3, linkVids4);
	fids = get_insert_fids(mesh, linkVids4, linkVids5);

	return true;
}

void strict_simplify(Mesh& mesh, BaseComplexQuad& baseComplex, std::set<size_t>& canceledFids) {
	size_t id = 0;
	for (const auto& link : baseComplex.separatedVertexIdsLink) {
		const auto& linkEids = baseComplex.separatedEdgeIdsLink.at(id);
		auto& v_front = mesh.V.at(link.front());
		auto& v_back = mesh.V.at(link.back());
		if ((v_front.N_Fids.size() <= 5 && v_back.N_Fids.size() <= 5) && v_front.N_Fids.size() != v_back.N_Fids.size()) {
			; //ofs << 0 << std::endl;
		} else if (COLLAPSE && v_front.N_Fids.size() == 3 && v_back.N_Fids.size() == 3) {
			// ofs << 1 << std::endl;
			auto v_front_fid = get_faceid(mesh, v_front.id, link[1]);
			auto v_front_fvid = get_diagnal_vid(mesh, v_front.id, v_front_fid);
			auto v_back_fid = get_faceid(mesh, v_back.id, link[link.size() - 2]);
			auto v_back_fvid = get_diagnal_vid(mesh, v_back.id, v_back_fid);
			bool condition = false;
			if (mesh.V.at(v_front_fvid).N_Fids.size() == 5 && mesh.V.at(v_back_fvid).N_Fids.size() == 5) {
				if (can_collapse_with_feature_preserved(mesh, link, linkEids, v_front_fvid, v_back_fvid)) {
					for (auto vid : link) {
						auto& v = mesh.V.at(vid);
						canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
					}
					collapse_with_feature_preserved(mesh, link, linkEids);
					break;
				}
			}
		} else if (SPLIT && v_front.N_Fids.size() == 5 && v_back.N_Fids.size() == 5) {
			//ofs << 2 << std::endl;
			auto v_front_eid = baseComplex.separatedEdgeIdsLink.at(id).front();
			auto v_back_eid = baseComplex.separatedEdgeIdsLink.at(id).back();
			auto& v_front_fids = mesh.E.at(v_front_eid).N_Fids;
			auto& v_back_fids = mesh.E.at(v_back_eid).N_Fids;

			auto v_front_fvid = get_diagnal_vid(mesh, v_front.id, v_front_fids);
			auto v_back_fvid = get_diagnal_vid(mesh, v_back.id, v_back_fids);
			bool condition = false;
			if (mesh.V.at(v_front_fvid).N_Fids.size() == 3 && mesh.V.at(v_back_fvid).N_Fids.size() == 3) {
				if (!split_with_feature_preserved(mesh, link, linkEids, v_front_fvid, v_back_fvid)) {
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

}

void loose_simplify(Mesh& mesh, BaseComplexQuad& baseComplex, std::set<size_t>& canceledFids) {
	size_t id = 0;
	for (const auto& link : baseComplex.separatedVertexIdsLink) {
		const auto& linkEids = baseComplex.separatedEdgeIdsLink.at(id);
		auto& v_front = mesh.V.at(link.front());
		auto& v_back = mesh.V.at(link.back());
		if ((v_front.N_Fids.size() <= 5 && v_back.N_Fids.size() <= 5) && v_front.N_Fids.size() != v_back.N_Fids.size()) {
			; //ofs << 0 << std::endl;
		} else if (COLLAPSE && v_front.N_Fids.size() == 3 && v_back.N_Fids.size() == 3) {
			// ofs << 1 << std::endl;
			auto v_front_fid = get_faceid(mesh, v_front.id, link[1]);
			auto v_front_fvid = get_diagnal_vid(mesh, v_front.id, v_front_fid);
			auto v_back_fid = get_faceid(mesh, v_back.id, link[link.size() - 2]);
			auto v_back_fvid = get_diagnal_vid(mesh, v_back.id, v_back_fid);
			if (mesh.V.at(v_front_fvid).N_Fids.size() > minValence && mesh.V.at(v_back_fvid).N_Fids.size() > minValence) {
				if (can_collapse_with_feature_preserved(mesh, link, linkEids, v_front_fvid, v_back_fvid)) {
					for (auto vid : link) {
						auto& v = mesh.V.at(vid);
						canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
					}
					collapse_with_feature_preserved(mesh, link, linkEids);
					break;
				}
			}
		} else if (SPLIT && v_front.N_Fids.size() == 5 && v_back.N_Fids.size() == 5) {
			//ofs << 2 << std::endl;
			auto v_front_eid = baseComplex.separatedEdgeIdsLink.at(id).front();
			auto v_back_eid = baseComplex.separatedEdgeIdsLink.at(id).back();
			auto& v_front_fids = mesh.E.at(v_front_eid).N_Fids;
			auto& v_back_fids = mesh.E.at(v_back_eid).N_Fids;

			auto v_front_fvid = get_diagnal_vid(mesh, v_front.id, v_front_fids);
			auto v_back_fvid = get_diagnal_vid(mesh, v_back.id, v_back_fids);
			bool condition = false;
			if (mesh.V.at(v_front_fvid).N_Fids.size() < maxValence && mesh.V.at(v_back_fvid).N_Fids.size() < maxValence) {
				if (split_with_feature_preserved(mesh, link, linkEids, v_front_fvid, v_back_fvid)) {
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
		++id;
	}
}

void simplify(Mesh& mesh, BaseComplexQuad& baseComplex, std::set<size_t>& canceledFids) {
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
				auto v_front_fid = get_faceid(mesh, v_front.id, link[1]);
				auto v_front_fvid = get_diagnal_vid(mesh, v_front.id, v_front_fid);
				auto v_back_fid = get_faceid(mesh, v_back.id, link[link.size() - 2]);
				auto v_back_fvid = get_diagnal_vid(mesh, v_back.id, v_back_fid);
				bool condition = false;
				if ((t == 0 && mesh.V.at(v_front_fvid).N_Fids.size() == 5 && mesh.V.at(v_back_fvid).N_Fids.size() == 5)
					|| (t == 1 && (mesh.V.at(v_front_fvid).N_Fids.size() > minValence && mesh.V.at(v_back_fvid).N_Fids.size() > minValence)))
					condition = true;
				if (condition) {
					if (can_collapse_with_feature_preserved(mesh, link, linkEids)) {
						for (auto vid : link) {
							auto& v = mesh.V.at(vid);
							canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
						}
						collapse_with_feature_preserved(mesh, link, linkEids);
						break;
					}
				}
			} else if (SPLIT && v_front.N_Fids.size() == 5 && v_back.N_Fids.size() == 5) {
				//ofs << 2 << std::endl;
				auto v_front_eid = baseComplex.separatedEdgeIdsLink.at(id).front();
				auto v_back_eid = baseComplex.separatedEdgeIdsLink.at(id).back();
				auto& v_front_fids = mesh.E.at(v_front_eid).N_Fids;
				auto& v_back_fids = mesh.E.at(v_back_eid).N_Fids;

				auto v_front_fvid = get_diagnal_vid(mesh, v_front.id, v_front_fids);
				auto v_back_fvid = get_diagnal_vid(mesh, v_back.id, v_back_fids);
				bool condition = false;
				if ((t == 0 && mesh.V.at(v_front_fvid).N_Fids.size() == 3 && mesh.V.at(v_back_fvid).N_Fids.size() == 3)
					|| (t == 1 && (mesh.V.at(v_front_fvid).N_Fids.size() < maxValence && mesh.V.at(v_back_fvid).N_Fids.size() < maxValence)))
					condition = true;
				if (condition) {
					if (!split_with_feature_preserved(mesh, link, linkEids, v_front_fvid, v_back_fvid)) {
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

void collapse_(Mesh& mesh, BaseComplexQuad& baseComplex, std::set<size_t>& canceledFids) {
	auto copyMesh = mesh;
	std::set<size_t> total_regionVids;
	for (auto t = 0; t < 2; ++t) {
		size_t id = 0;
		for (const auto& link : baseComplex.separatedVertexIdsLink) {
			auto& v_front = mesh.V.at(link.front());
			auto& v_back = mesh.V.at(link.back());
			if ((v_front.N_Fids.size() <= 5 && v_back.N_Fids.size() <= 5) && v_front.N_Fids.size() != v_back.N_Fids.size()) {
				; //ofs << 0 << std::endl;
			} else if (v_front.N_Fids.size() == 3 && v_back.N_Fids.size() == 3) {
				// ofs << 1 << std::endl;
				auto v_front_fid = get_faceid(mesh, v_front.id, link[1]);
				auto v_front_fvid = get_diagnal_vid(mesh, v_front.id, v_front_fid);

				auto v_back_fid = get_faceid(mesh, v_back.id, link[link.size() - 2]);
				auto v_back_fvid = get_diagnal_vid(mesh, v_back.id, v_back_fid);

				bool condition = false;
				if ((t == 0 && mesh.V.at(v_front_fvid).N_Fids.size() == 5 && mesh.V.at(v_back_fvid).N_Fids.size() == 5)
					|| (t == 1 && (mesh.V.at(v_front_fvid).N_Fids.size() > 3 && mesh.V.at(v_back_fvid).N_Fids.size() > 3)))
					condition = true;
				if (condition) {
					if (can_collapse_with_feature_preserved(mesh, baseComplex.separatedVertexIdsLink.at(id),
						baseComplex.separatedEdgeIdsLink.at(id))) {
						auto regionVids = get_regionVids(copyMesh, baseComplex.separatedVertexIdsLink.at(id));
						bool found = false;
						for (auto rvid : regionVids)
							if (total_regionVids.find(rvid) != total_regionVids.end()) {
								found = true;
								break;
							}
						if (found) continue;
						total_regionVids.insert(regionVids.begin(), regionVids.end());
						for (auto vid : link) {
							auto& v = copyMesh.V.at(vid);
							canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
						}
						collapse_with_feature_preserved(mesh, baseComplex.separatedVertexIdsLink.at(id),
							baseComplex.separatedEdgeIdsLink.at(id));

						// break;
					}
				}
			} else if (v_front.N_Fids.size() == 5 && v_back.N_Fids.size() == 5) {
				; //ofs << 2 << std::endl;
			} else {
				; //ofs << 3 << std::endl;
			}
			++id;
		}
		if (!canceledFids.empty())  break;
	}
}
void smooth_project(Mesh& origMesh, Mesh& mesh) {
	std::map<size_t, std::set<size_t>> origLabel_vids;
	std::map<size_t, std::set<size_t>> origPatch_vids;
	std::map<size_t, std::set<size_t>> origLabel_eids;
	std::map<size_t, std::set<size_t>> origSharpEdgeVid_NVids;

	std::map<size_t, std::set<size_t>> label_vids;
	std::map<size_t, std::set<size_t>> label_eids;
	std::map<size_t, std::set<size_t>> sharpEdgeVid_NVids;
	Vertex vertex;
	for (auto& f : origMesh.F) {
		glm::dvec3 center(0, 0, 0);
		for (auto nvid : origMesh.F.at(f.id).Vids)
			center += origMesh.V.at(nvid).xyz();
		center *= 0.25;
		vertex = center;
		vertex.id = origMesh.V.size();
		vertex.patch_id = f.label;
		origMesh.V.push_back(vertex);
		//origPatch_vids[vertex.patch_id].insert(vertex.id);
	}
	for (auto& e : origMesh.E) {
		glm::dvec3 center(0, 0, 0);
		for (auto nvid : e.Vids)
			center += origMesh.V.at(nvid).xyz();
		center *= 0.5;
		vertex = center;
		vertex.id = origMesh.V.size();

		if (!e.isSharpFeature) {
			vertex.patch_id = origMesh.V.at(e.Vids[0]).type == REGULAR ? origMesh.V.at(e.Vids[0]).patch_id : origMesh.V.at(e.Vids[1]).patch_id;
		} else {
			continue;
			vertex.label = origMesh.V.at(e.Vids[0]).label == MAXID ? origMesh.V.at(e.Vids[1]).label : origMesh.V.at(e.Vids[0]).label;
			vertex.type = FEATURE;
		}
		origMesh.V.push_back(vertex);
		//origPatch_vids[vertex.patch_id].insert(vertex.id);
	}
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
			if (v.isCorner) continue;
			glm::dvec3 center(0, 0, 0);
			for (auto nvid : item.second)
				center += mesh.V.at(nvid).xyz();
			center /= item.second.size();
			const auto& origLineVids = origLabel_vids[v.label];
			size_t closest_origLineVid = *origLineVids.begin();
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
		}
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
			// if (iters < 3) continue;

			auto& patchVids = origPatch_vids[v.patch_id];
			size_t closest_origVid = *patchVids.begin();
			double closest_distance = 100000000.0;
			for (auto patchVid : patchVids) {
				auto& origv = origMesh.V.at(patchVid);
				auto distance = glm::length(origv.xyz() - center);
				if (distance < closest_distance) {
					closest_origVid = origv.id;
					closest_distance = distance;
				}
			}

			v = origMesh.V.at(closest_origVid).xyz();
		}
	}
}

Mesh RefineWithFeaturePreserved(const Mesh& hex_mesh, int clockwise) {
	const Mesh& new_mesh = hex_mesh;
	////////////////////////////////////////////////////////////////////////////
	// add vertices
	std::vector<Vertex> new_vertex(new_mesh.V.size() + new_mesh.E.size() + new_mesh.C.size());
	for (size_t i = 0; i < new_mesh.V.size(); i++) {
		auto& v = new_vertex.at(i);
		auto& newv= new_mesh.V.at(i);
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

Mesh RefineWithFeaturePreserved2(const Mesh& hex_mesh, int clockwise) {
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

std::vector<size_t> get_ids(const std::string str) {
	std::stringstream ss(str);
	std::vector<size_t> res;
	size_t id;
	while (ss >> id) res.push_back(id);
	return res;
}

std::map<size_t, std::set<size_t>> get_patchid_fids(const Mesh& mesh) {
	std::map<size_t, std::set<size_t>> patchid_fids;
	for (auto& v : mesh.V)
		if (v.patch_id != MAXID) patchid_fids[v.patch_id].insert(v.N_Fids.begin(), v.N_Fids.end());
	return patchid_fids;
}

std::map<size_t, std::set<size_t>> get_patchid_vids(const Mesh& mesh, const std::map<size_t, std::set<size_t>>& patchid_fids) {
	std::map<size_t, std::set<size_t>> patchid_vids;
	for (auto& item : patchid_fids) {
		for (auto fid : item.second) {
			auto& f = mesh.F.at(fid);
			patchid_vids[item.first].insert(f.Vids.begin(), f.Vids.end());
		}
	}
	return patchid_vids;
}

std::set<size_t> get_rotate_fids(const Mesh& mesh) {
	std::set<size_t> res;
	auto patchid_fids = get_patchid_fids(mesh);
	auto patchid_vids = get_patchid_vids(mesh, patchid_fids);

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

std::vector<size_t> get_neighbor_fids(const Vertex& v, const std::set<size_t>& patch_fids) {
	std::vector<size_t> fids;
	for (auto fid : v.N_Fids)
		if (patch_fids.find(fid) != patch_fids.end()) fids.push_back(fid);
	return fids;
}

std::set<size_t> get_rotate_eids(const Mesh& mesh, const Vertex& v, const std::vector<size_t>& fids) {
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
std::set<size_t> get_rotate_eids_(const Mesh& mesh) {
	std::set<size_t> res;
	auto patchid_fids = get_patchid_fids(mesh);
	auto patchid_vids = get_patchid_vids(mesh, patchid_fids);
	auto convex_corners = get_convex_corners(mesh);
	for (auto& item : patchid_vids) {
		for (auto vid : item.second) {
			auto& v = mesh.V.at(vid);
			if (v.type == CORNER && !v.isSpecial) {
				auto fids = get_neighbor_fids(v, patchid_fids[item.first]);
				if (fids.size() < 2) continue;
				if (convex_corners.find(v.id) == convex_corners.end() && fids.size() < 4) continue;
				auto eids = get_rotate_eids(mesh, v, fids);
				res.insert(eids.begin(), eids.end());
			} else if (v.type == FEATURE) {
				auto fids = get_neighbor_fids(v, patchid_fids[item.first]);
				if (fids.size() < 3 || fids.size() == 4) continue;
				auto eids = get_rotate_eids(mesh, v, fids);
				res.insert(eids.begin(), eids.end());
			 }
		}
	}

	return res;
}

double get_angle(const Vertex& v, const Vertex& v0, const Vertex& v1) {
	auto d0 = v0.xyz() - v.xyz();
	auto d1 = v1.xyz() - v.xyz();
	auto cosangle = glm::dot(glm::normalize(d0), glm::normalize(d1));
	return acos(cosangle) * 180.0 / PI;
}

std::vector<size_t> get_neighbor_vids(const Mesh& mesh, const Vertex& v, size_t fid) {
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

bool is_convex(const Mesh& mesh, const Vertex& v, const std::vector<size_t>& fids) {
	double total_angle = 0;
	for (auto fid : fids) {
		auto vids = get_neighbor_vids(mesh, v, fid);
		total_angle += get_angle(v, mesh.V.at(vids[0]), mesh.V.at(vids[1]));
	}
	return total_angle < 135;
}

std::set<size_t> get_rotate_eids(const Mesh& mesh) {
	std::set<size_t> res;
	auto patchid_fids = get_patchid_fids(mesh);
	auto patchid_vids = get_patchid_vids(mesh, patchid_fids);
	for (auto& item : patchid_vids) {
		for (auto vid : item.second) {
			auto& v = mesh.V.at(vid);
			if (v.type == CORNER && !v.isSpecial) {
				auto fids = get_neighbor_fids(v, patchid_fids[item.first]);
				if (fids.size() < 2) continue;
				if (!is_convex(mesh, v, fids) && fids.size() < 4) continue;
				auto eids = get_rotate_eids(mesh, v, fids);
				res.insert(eids.begin(), eids.end());
			} else if (v.type == FEATURE) {
				auto fids = get_neighbor_fids(v, patchid_fids[item.first]);
				if (fids.size() < 3 || fids.size() == 4) continue;
				auto eids = get_rotate_eids(mesh, v, fids);
				res.insert(eids.begin(), eids.end());
			}
		}
	}

	return res;
}

void rotate(Mesh& mesh, const Edge& e, const Vertex& v, std::set<size_t>& canceledFids) {
	auto& f0 = mesh.F.at(e.N_Fids[0]);
	auto& f1 = mesh.F.at(e.N_Fids[1]);
	auto vid0 = e.Vids[0] == v.id ? e.Vids[1] : e.Vids[0];
	canceledFids.insert(f0.id);
	canceledFids.insert(f1.id);
	auto diag_vid0 = get_diagnal_vid(mesh, vid0, f0.id);
	auto diag_vid1 = get_diagnal_vid(mesh, v.id, f1.id);
	if (mesh.V.at(diag_vid0).type >= FEATURE) {
		diag_vid0 = get_diagnal_vid(mesh, v.id, f0.id);
		diag_vid1 = get_diagnal_vid(mesh, vid0, f1.id);
	}
	auto eids = get_boundary_eids(f0, f1, e);
	auto linkVids = get_rotate_vids(mesh, eids, diag_vid0, diag_vid1);
	insert_rotate_faces(mesh, linkVids);
}

void rotate(Mesh& mesh, std::set<size_t>& canceledFids) {
	auto patchid_fids = get_patchid_fids(mesh);
	auto patchid_vids = get_patchid_vids(mesh, patchid_fids);
	for (auto& item : patchid_vids) {
		for (auto vid : item.second) {
			auto& v = mesh.V.at(vid);
			if (v.type == CORNER && !v.isSpecial) {
				auto fids = get_neighbor_fids(v, patchid_fids[item.first]);
				if (fids.size() < 2) continue;
				if (!is_convex(mesh, v, fids) && fids.size() < 4) continue;
				auto eids = get_rotate_eids(mesh, v, fids);
				auto& e = mesh.E.at(*eids.begin());
				rotate(mesh, e, v, canceledFids);
				return;
			} else if (v.type == FEATURE) {
				auto fids = get_neighbor_fids(v, patchid_fids[item.first]);
				if (fids.size() < 3 || fids.size() == 4) continue;
				auto eids = get_rotate_eids(mesh, v, fids);
				auto& e = mesh.E.at(*eids.begin());
				rotate(mesh, e, v, canceledFids);
				return;
			}
		}
	}

}

std::vector<size_t> GetLinkVids(const Mesh& mesh, const std::vector<size_t>& vids) {
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

std::vector<size_t> GetLinkEVids(const Mesh& mesh, const std::vector<size_t>& eids) {
	std::set<size_t> vids;
	for (auto eid : eids) {
		const auto& e = mesh.E.at(eid);
		vids.insert(e.Vids.begin(), e.Vids.end());
	}
	std::vector<size_t> res;
	for (auto vid : vids) res.push_back(vid);
	return res;
}

std::vector<size_t> GetLinkVidsFromEids(const Mesh& mesh, const std::vector<size_t>& eids) {
	std::set<size_t> eids_set(eids.begin(), eids.end());
	const auto vids = GetLinkEVids(mesh, eids);
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
			return { };
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

std::vector<size_t> get_boundary_eids(const Face& f0, const Face& f1, const Edge& exclude_e) {
	std::set<size_t> eids_set(f0.Eids.begin(), f0.Eids.end());
	eids_set.insert(f1.Eids.begin(), f1.Eids.end());
	eids_set.erase(exclude_e.id);
	std::vector<size_t> eids;
	std::copy(eids_set.begin(), eids_set.end(), std::back_inserter(eids));
	return eids;
}

std::vector<size_t> get_rotate_vids(const Mesh& mesh, const std::vector<size_t>& boundary_eids, size_t diag_vid0, size_t diag_vid1) {
	auto linkVids = GetLinkVidsFromEids(mesh, boundary_eids);
	linkVids.pop_back();
	while (linkVids.front() != diag_vid0 && linkVids.front() != diag_vid1) {
		linkVids.push_back(linkVids.front());
		linkVids.erase(linkVids.begin());
	}
	return linkVids;
}

void insert_rotate_faces(Mesh& mesh, const std::vector<size_t>& rotate_vids) {
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

bool can_rotate1(const Mesh& mesh, const Edge& e) {
	auto& v0 = mesh.V.at(e.Vids[0]);
	auto& v1 = mesh.V.at(e.Vids[1]);
	if (v0.type == CORNER && v1.type == CORNER) return false;
	if (v0.type == CORNER && v1.type == FEATURE) return false;
	if (v1.type == CORNER && v0.type == FEATURE) return false;
	if (v1.type == FEATURE && v0.type == FEATURE && v0.label == v1.label) return false;
	return true;
}

bool can_rotate(const Mesh& mesh, const Edge& e) {
	auto& v0 = mesh.V.at(e.Vids[0]);
	auto& v1 = mesh.V.at(e.Vids[1]);
	if (v0.type == CORNER && v1.type == CORNER) return false;
	if (v0.type == CORNER && v1.type == FEATURE && v0.labels.find(v1.label) != v0.labels.end()) return false;
	if (v1.type == CORNER && v0.type == FEATURE && v1.labels.find(v0.label) != v1.labels.end()) return false;
	if (v1.type == FEATURE && v0.type == FEATURE && v0.label == v1.label) return false;
	return true;
}

void rotate(Mesh& mesh, const Face& f0, const Face& f1, const Edge& e, std::set<size_t>& canceledFids) {
	canceledFids.insert(f0.id);
	canceledFids.insert(f1.id);
	auto diag_vid0 = get_diagnal_vid(mesh, e.Vids[0], f0.id);
	auto diag_vid1 = get_diagnal_vid(mesh, e.Vids[1], f1.id);
	auto eids = get_boundary_eids(f0, f1, e);
	auto linkVids = get_rotate_vids(mesh, eids, diag_vid0, diag_vid1);
	insert_rotate_faces(mesh, linkVids);
}

bool rotate(Mesh& mesh, const Face& f0, const Face& f1, std::set<size_t>& canceledFids) {
	std::set<size_t> f0eids(f0.Eids.begin(), f0.Eids.end());
	std::set<size_t> f1eids(f1.Eids.begin(), f1.Eids.end());
	auto intersect_eids = Util::get_intersect(f0eids, f1eids);
	if (intersect_eids.size() != 1) {
		// std::cerr << "Err in auto eids = Util::get_intersect(f0eids, f1eids);\n";
		return false;
	}
	auto& e = mesh.E.at(*intersect_eids.begin());
	if (!can_rotate(mesh, e)) return false;
	rotate(mesh, f0, f1, e, canceledFids);
	return true;
}

bool rotate1(Mesh& mesh, const std::vector<size_t>& fids, std::set<size_t>& canceledFids) {
	auto& f0 = mesh.F.at(fids.at(0));
	auto& f1 = mesh.F.at(fids.at(1));
	if (!rotate(mesh, f0, f1, canceledFids) && fids.size() == 3) {
		auto& f2 = mesh.F.at(fids.at(2));
		if (!rotate(mesh, f0, f2, canceledFids)) return false;
	} else if (fids.size() == 4) {
		auto& f2 = mesh.F.at(fids.at(2));
		auto& f3 = mesh.F.at(fids.at(3));
		if (!rotate(mesh, f0, f2, canceledFids) &&
			!rotate(mesh, f0, f3, canceledFids) &&
			!rotate(mesh, f1, f2, canceledFids)) return false;
	}
	return true;
}

bool rotate(Mesh& mesh, const std::vector<size_t>& fids, std::set<size_t>& canceledFids) {
	auto pairs = Util::combine(fids.size(), 2);
	for (auto& p : pairs) {
		auto& f0 = mesh.F.at(fids.at(p[0]));
		auto& f1 = mesh.F.at(fids.at(p[1]));
		if (!rotate(mesh, f0, f1, canceledFids)) continue;
		return true;
	}
	return false;
}

void insert_rotate_fids(Mesh& mesh, std::set<size_t>& canceledFids) {
	auto patchid_fids = get_patchid_fids(mesh);
	auto patchid_vids = get_patchid_vids(mesh, patchid_fids);

	for (auto& item : patchid_vids) {
		for (auto vid : item.second) {
			auto& v = mesh.V.at(vid);
			if (v.type == CORNER && !v.isSpecial) {
				auto fids = get_neighbor_fids(v, patchid_fids[item.first]);
				if (fids.size() < 2) continue;
				if (!rotate(mesh, fids, canceledFids)) continue;
				return;
			} else if (v.type == FEATURE) {
				auto fids = get_neighbor_fids(v, patchid_fids[item.first]);
				if (fids.size() < 3) continue;
				if (!rotate(mesh, fids, canceledFids)) continue;
				return;
			}
		}
	}
}

void remove_doublet(Mesh& mesh, std::set<size_t>& canceledFids) {
	for (auto& v : mesh.V) {
		if (v.N_Fids.size() != 2) continue;
		auto& v0 = mesh.V.at(v.N_Vids[0]);
		auto& v1 = mesh.V.at(v.N_Vids[1]);
		if (v0.type >= FEATURE && v1.type >= FEATURE && v.type < FEATURE) continue;
		auto& f0 = mesh.F.at(v.N_Fids.front());
		auto& f1 = mesh.F.at(v.N_Fids.back());
		std::set<size_t> eids_set(f0.Eids.begin(), f0.Eids.end());
		eids_set.insert(f1.Eids.begin(), f1.Eids.end());
		for (auto eid : v.N_Eids)
			eids_set.erase(eid);
		std::vector<size_t> eids;
		std::copy(eids_set.begin(), eids_set.end(), std::back_inserter(eids));
		auto linkVids = GetLinkVidsFromEids(mesh, eids);
		if (linkVids.empty()) continue;
		linkVids.pop_back();
		Face newf;
		newf.id = mesh.F.size();
		newf.Vids = linkVids;
		mesh.F.push_back(newf);
		canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
		return;
	}
}

void collapse_diagnal1(Mesh& mesh, std::set<size_t>& canceledFids) {
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

void collapse_diagnal(Mesh& mesh, std::set<size_t>& canceledFids) {
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
	if (canceledFids.empty()) collapse_diagnal1(mesh, canceledFids);
}

bool simplify(Mesh& mesh, int& iter) {
	std::set<size_t> canceledFids;
	mesh.CompressWithFeaturePreserved();
	mesh.BuildAllConnectivities();
	mesh.ExtractBoundary();
	mesh.ExtractSingularities();
	mesh.BuildParallelE();
	if (iter == 0 && featurePreserved) {
		get_feature(mesh);
		auto eids = get_rotate_eids(mesh);
		{
			MeshFileWriter writer(mesh, "rotate_eids.vtk");
			writer.WriteEdgesVtk(eids);
		}
		auto count = 0;
		while (true) {
			std::cout << "rotate " << count++ << std::endl;
			if (canceledFids.empty() && REMOVE_DOUBLET) {
				remove_doublet(mesh, canceledFids);
				if (!canceledFids.empty()) {
					std::cout << "remove_doublet" << std::endl;
					update(mesh, canceledFids);
					canceledFids.clear();
					mesh.CompressWithFeaturePreserved();
					mesh.BuildAllConnectivities();
					mesh.ExtractBoundary();
					mesh.ExtractSingularities();
					mesh.BuildParallelE();
					continue;
				}
			}
			rotate(mesh, canceledFids);
			if (canceledFids.empty()) break;
			update(mesh, canceledFids);
			canceledFids.clear();
			mesh.CompressWithFeaturePreserved();
			mesh.BuildAllConnectivities();
			mesh.ExtractBoundary();
			mesh.ExtractSingularities();
			mesh.BuildParallelE();
		}
		{
			std::cout << "writing rotate.vtk "<< std::endl;
			MeshFileWriter writer(mesh, "rotate.vtk");
			writer.WriteFile();
		}
	}
	BaseComplexQuad baseComplex(mesh);
	baseComplex.ExtractSingularVandE();
	baseComplex.BuildE();
	strict_simplify(mesh, baseComplex, canceledFids);
	if (canceledFids.empty()) {
		std::cout << "loose_simplify\n";
		loose_simplify(mesh, baseComplex, canceledFids);
	}
	std::cout << "iter = " << iter++ << std::endl;
	if (canceledFids.empty() && global) {
		update(mesh, canceledFids);
		mesh.CompressWithFeaturePreserved();
		mesh.BuildAllConnectivities();
		mesh.ExtractBoundary();
		mesh.ExtractSingularities();
		mesh.BuildParallelE();
		global_simplify(mesh, canceledFids);
		std::cout << "Finished global_simplify\n";
	}
	if (canceledFids.empty() && REMOVE_DOUBLET) {
		std::cout << "remove_doublet" << std::endl;
		remove_doublet(mesh, canceledFids);
	}
	if (canceledFids.empty() && ROTATE) {
		std::cout << "rotate_edges" << std::endl;
		rotate(mesh, canceledFids);
	}
	if (canceledFids.empty() && COLLAPSE_DIAGNAL) {
		std::cout << "collapse_diagnal" << std::endl;
		collapse_diagnal(mesh, canceledFids);
		if (!canceledFids.empty()) {
			collapse_diagnal(mesh, canceledFids);
			update(mesh, canceledFids);
			MeshFileWriter writer(mesh, "collapse_diagnal.vtk");
			writer.WriteFile();
			return true;
		}
	}
	if (canceledFids.empty()) return false;
	update(mesh, canceledFids);
	//std::string fname = std::string("iter") + std::to_string(iter) + ".vtk";
	//MeshFileWriter writer(mesh, fname.c_str());
	//writer.WriteFile();
	return true;
}


size_t permutation[3][2] = { { 0, 1 }, { 1, 2 }, { 2, 0 }};

const double PI2 = 3.1415926 * 2;
double GetAngle(const Mesh& mesh, const Vertex& v, const Face& c) {
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

std::set<size_t> get_convex_corners(const Mesh& mesh) {
	std::set<size_t> res;
	for (auto& v : mesh.V) {
		if (!v.isCorner || v.isSpecial) continue;
		double gaussianCurvature = 0;
		for (auto fid : v.N_Fids) {
			const Face& f = mesh.F.at(fid);
			gaussianCurvature += GetAngle(mesh, v, f);
		}
		gaussianCurvature = PI2 - gaussianCurvature;
		if (gaussianCurvature > PI) res.insert(v.id);
	}
	return res;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: QuadMeshLocalSimplify2 quad.vtk simplified.vtk iters=<1> maxValence=<3> maxValence=<5> " 
			<< "featurePreserved=true smoothIters=20 angle=<160> userCorners=\"\" canceledCorners=\"\"" <<
			" collapse=false split=true conformal=true global=true rotate=true remove_doublet=true collapse_diagnal=true "<< std::endl;
        return -1;
    }
    ArgumentManager argumentManager(argc, argv);
    std::string input = argv[1];
    std::string output = argv[2];
    int iters = 10000;
    auto strIters = argumentManager.get("iters");
    if (!strIters.empty()) iters = std::stoi(strIters);
	auto strmaxValence = argumentManager.get("maxValence");
	if (!strmaxValence.empty()) maxValence = std::stoi(strmaxValence);
	auto strminValence = argumentManager.get("minValence");
	if (!strminValence.empty()) minValence = std::stoi(strminValence);
	auto strsmoothIters = argumentManager.get("smoothIters");
	if (!strsmoothIters.empty()) smoothIters = std::stoi(strsmoothIters);

	auto strfeaturePreserved = argumentManager.get("featurePreserved");
	if (!strfeaturePreserved.empty()) featurePreserved = strfeaturePreserved == "false" ? false : true;

	auto strcollapse = argumentManager.get("collapse");
	if (!strcollapse.empty()) COLLAPSE = strcollapse == "false" ? false : true;

	auto strsplit = argumentManager.get("split");
	if (!strsplit.empty()) SPLIT = strsplit == "false" ? false : true;

	auto strrotate = argumentManager.get("rotate");
	if (!strrotate.empty()) ROTATE = strrotate == "false" ? false : true;

	auto strremove_doublet = argumentManager.get("remove_doublet");
	if (!strremove_doublet.empty()) REMOVE_DOUBLET = strremove_doublet == "false" ? false : true;

	auto strcollapse_diagnal = argumentManager.get("collapse_diagnal");
	if (!strcollapse_diagnal.empty()) COLLAPSE_DIAGNAL = strcollapse_diagnal == "false" ? false : true;

	auto strconformal = argumentManager.get("conformal");
	if (!strconformal.empty()) conformal = strconformal == "false" ? false : true;

	auto strglobal = argumentManager.get("global");
	if (!strglobal.empty()) global = strglobal == "false" ? false : true;

	auto strAngle = argumentManager.get("angle");
	if (!strAngle.empty()) angle = std::stod(strAngle);

	auto struserCorners = argumentManager.get("userCorners");
	if (!struserCorners.empty()) userCorners = get_ids(struserCorners);

	auto strcanceledCorners = argumentManager.get("canceledCorners");
	if (!strcanceledCorners.empty()) canceledCorners = get_ids(strcanceledCorners);
	
    std::cout << "---------------------------------------" << std::endl;
    std::cout << "input  = " << input << std::endl;
    std::cout << "output = " << output << std::endl;
	std::cout << "iters = " << iters << std::endl;
	std::cout << "maxValence = " << maxValence << std::endl;
	std::cout << "minValence = " << minValence << std::endl;
	std::cout << "smoothIters = " << smoothIters << std::endl;
	std::cout << "featurePreserved = " << featurePreserved << std::endl;
	std::cout << "angle = " << angle << std::endl;
	std::cout << "collapse = " << COLLAPSE << std::endl;
	std::cout << "split = " << SPLIT << std::endl;
	std::cout << "conformal = " << conformal << std::endl;
	std::cout << "global = " << global << std::endl;
	std::cout << "rotate = " << ROTATE << std::endl;
	std::cout << "remove_doublet = " << REMOVE_DOUBLET << std::endl;
	std::cout << "collapse_diagnal = " << COLLAPSE_DIAGNAL << std::endl;
    std::cout << "---------------------------------------" << std::endl;
    int iter = 0;
    MeshFileReader reader(input.c_str());
    Mesh& mesh = (Mesh&) reader.GetMesh();
	auto maxValence_copy = maxValence;
	if (maxValence_copy > 5) maxValence = 5;
	while (maxValence <= maxValence_copy) {
		while (iters--) {
			if (!simplify(mesh, iter)) break;
		}
		++maxValence;
	}
	{
		MeshFileWriter writer(mesh, "VertexFeature.vtk");
		writer.WriteVertexFeatureVtk();
	}
	{
		auto fids = get_rotate_fids(mesh);
		MeshFileWriter writer(mesh, "RotateFaces.vtk");
		writer.WriteFacesVtk(fids);
	}
	smooth_project(origMesh, mesh); 
	{
		MeshFileWriter writer(mesh, output.c_str());
		writer.WriteFile();
		std::cout << "V = " << mesh.V.size() << std::endl;
		std::cout << "E = " << mesh.E.size() << std::endl;
		std::cout << "F = " << mesh.F.size() << std::endl;
	}
	//{
	//	auto refinedMesh = RefineWithFeaturePreserved(mesh, 0);
	//	refinedMesh.BuildAllConnectivities();
	//	MeshFileWriter writer(refinedMesh, "refine.vtk");
	//	writer.WriteFile();
	//	MeshFileWriter writer0(refinedMesh, "refine_VertexFeature.vtk");
	//	writer0.WriteVertexFeatureVtk();
	//	smooth_project(origMesh, refinedMesh);
	//	MeshFileWriter writer1(refinedMesh, "refine_smooth.vtk");
	//	writer1.WriteFile();
	//	//writer1.WriteVertexFeatureVtk();
	//}
    return 0;
}
