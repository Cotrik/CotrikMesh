/*
* MapLoops.cpp
*
*  Created on: Nov 25, 2018
*      Author: cotrik
*/

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "BaseComplex.h"
#include "BaseComplexQuad.h"
#include "BaseComplexSheetQuad.h"
#include "RefinedDual.h"
#include "ArgumentManager.h"
#include <iostream>
#include <fstream>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <stack>
#include <algorithm>

struct VertexUV : public Vertex {
	size_t father;
	glm::dvec2 uv;
};
void read(const char* filename, std::vector<Vertex>& V, std::vector<Face>& F) {
	std::ifstream ifs(filename);
	std::string line;
	while (getline(ifs, line)) {
		std::stringstream ss(line);
		std::string keyword;
		Vertex v;
		ss >> keyword;
		if (keyword == "Vertex") {
			ss >> v.id >> v.x >> v.y >> v.z;
			--v.id;
			V.push_back(v);
		} else if (keyword == "Face") {
			Face f;
			ss >> f.id;
			--f.id;
			size_t vid;
			while (ss >> vid)
				f.Vids.push_back(vid - 1);

			F.push_back(f);
		} else if (keyword == "Edge") {
			break;
		}
	}
}
void read(const char* filename, std::vector<VertexUV>& V, std::vector<Face>& F) {
	std::ifstream ifs(filename);
	std::string line;
	while (getline(ifs, line)) {
		std::stringstream ss(line);
		std::string keyword;
		VertexUV v;
		ss >> keyword;
		if (keyword == "Vertex") {
			ss >> v.id >> v.x >> v.y >> v.z;
			--v.id;
			std::string father; // , uv, vv;
			ss >> father;// >> uv >> vv;

			for (auto& c : father)
				if (c == '(' || c == ')') c = ' ';
			std::istringstream ss_father(father);
			ss_father >> father >> v.father;
			--v.father;

			V.push_back(v);
		} else if (keyword == "Face") {
			Face f;
			ss >> f.id;
			--f.id;
			size_t vid;
			while (ss >> vid)
				f.Vids.push_back(vid - 1);

			F.push_back(f);
		} else if (keyword == "Edge") {
			break;
		}
	}
}

void generate_quad_mesh(const std::vector<Vertex>& V, const std::vector<Face>& F, Mesh& quadMesh) {
	std::vector<Cell> C;
	Cell c;
	for (auto& f : F) {
		c.id = f.id;
		c.Vids = f.Vids;
		C.push_back(c);
	}

	quadMesh.V = V;
	quadMesh.C = C;
	quadMesh.m_cellType = QUAD;
	//quadMesh.BuildAllConnectivities();
	//quadMesh.ExtractBoundary();
}

void generate_quad_mesh(const std::vector<VertexUV>& V, const std::vector<Face>& F, Mesh& quadMesh) {
	std::vector<Vertex> VV;
	Vertex vv;
	for (auto& v : V) {
		vv = v;
		vv.id = v.id;
		VV.push_back(vv);
	}
	std::vector<Cell> C;
	Cell c;
	for (auto& f : F) {
		c.id = f.id;
		c.Vids = f.Vids;
		C.push_back(c);
	}

	quadMesh.V = VV;
	quadMesh.C = C;
	quadMesh.m_cellType = QUAD;
	quadMesh.BuildAllConnectivities();
	quadMesh.ExtractBoundary();
}

void write(const std::vector<Vertex>& V, const std::vector<Face>& F, const char* filename) {
	std::vector<Vertex> triV;
	Vertex x;
	for (auto& v : V) {
		x.id = v.id;
		x.x = v.x;
		x.y = v.y;
		x.z = v.z;
		triV.push_back(x);
	}
	std::vector<Cell> C;
	Cell c;
	for (auto& f : F) {
		c.id = f.id;
		c.Vids = f.Vids;
		C.push_back(c);
	}
	MeshFileWriter writer(triV, C, filename, F.front().Vids.size() == 3 ? TRIANGLE : QUAD);
	writer.WriteFile();
}

std::set<size_t> get_boundary_vids(const Mesh& uvtriMesh) {
	std::set<size_t> boundary_vids;
	for (auto& e : uvtriMesh.E)
		if (e.N_Fids.size() == 1) {
			boundary_vids.insert(e.Vids[0]);
			boundary_vids.insert(e.Vids[1]);
		}
	return boundary_vids;
}

std::set<size_t> get_boundary_eids(const Mesh& uvtriMesh) {
	std::set<size_t> edgeids;
	for (auto& e : uvtriMesh.E)
		if (e.N_Fids.size() == 1)
			edgeids.insert(e.id);
	return edgeids;
}


static std::unordered_map<size_t, size_t> get_key_edgeId(const Mesh& mesh) {
	std::unordered_map<size_t, size_t> key_edgeId;
	for (size_t i = 0; i < mesh.E.size(); ++i) {
		const auto& e = mesh.E.at(i);
		key_edgeId[(e.Vids[0] << 32) | e.Vids[1]] = i;
		key_edgeId[(e.Vids[1] << 32) | e.Vids[0]] = i;
	}
	return key_edgeId;
}

static std::unordered_map<size_t, size_t> get_key_edgeId(const RefinedDualQuad& mesh) {
	std::unordered_map<size_t, size_t> key_edgeId;
	for (size_t i = 0; i < mesh.E.size(); ++i) {
		const auto& e = mesh.E.at(i);
		key_edgeId[(e.Vids[0] << 32) | e.Vids[1]] = i;
		key_edgeId[(e.Vids[1] << 32) | e.Vids[0]] = i;
	}
	return key_edgeId;
}

static std::string make_facekey(const Face& f) {
	std::set<size_t> vids(f.Vids.begin(), f.Vids.end());
	std::string s;
	for (auto vid : vids)
		s += std::to_string(vid) + "@";
	return s;
}

static std::string make_facekey(size_t vid0, size_t vid1, size_t vid2, size_t vid3) {
	std::set<size_t> vids = { vid0, vid1, vid2, vid3 };
	std::string s;
	for (auto vid : vids)
		s += std::to_string(vid) + "@";
	return s;
}

static std::unordered_map<std::string, size_t> get_key_faceId(const Mesh& mesh) {
	std::unordered_map<std::string, size_t> key_faceId;
	for (size_t i = 0; i < mesh.F.size(); ++i) {
		const auto& f = mesh.F.at(i);
		key_faceId[make_facekey(f)] = i;
	}
	return key_faceId;
}

static std::unordered_map<std::string, size_t> get_key_faceId(const RefinedDualQuad& mesh) {
	std::unordered_map<std::string, size_t> key_faceId;
	for (size_t i = 0; i < mesh.F.size(); ++i) {
		const auto& f = mesh.F.at(i);
		key_faceId[make_facekey(f)] = i;
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

void WriteDualVTK(const RefinedDualQuad& dual, const std::vector<std::unordered_set<size_t>> all_dualEdgeIds, const char *filename) {
	size_t numOfLines = 0;
	for (const auto& dualEdgeIds : all_dualEdgeIds)
		numOfLines += dualEdgeIds.size();

	const auto& V = dual.V;
	const auto& E = dual.E;
	std::ofstream ofs(filename);
	ofs << "# vtk DataFile Version 2.0\n"
		<< filename << "\n"
		<< "ASCII\n\n"
		<< "DATASET POLYDATA\n";
	ofs << "POINTS " << V.size() << " double" << "\n";
	for (const auto& v : V)
		ofs << v.x << " " << v.y << " " << v.z << "\n";

	ofs << "Lines " << numOfLines << " " << 3 * numOfLines << "\n";
	for (const auto& dualEdgeIds : all_dualEdgeIds)
		for (const auto edge_id : dualEdgeIds) {
			const auto& edge = E.at(edge_id);
			ofs << edge.Vids.size() << " " << edge.Vids[0] << " " << edge.Vids[1] << "\n";
		}

	ofs << "CELL_DATA " << numOfLines << "\n"
		<< "SCALARS " << "id" << " int 1\n"
		<< "LOOKUP_TABLE default\n";
	size_t sheet_id = 0;
	for (const auto& dualEdgeIds : all_dualEdgeIds) {
		for (const auto id : dualEdgeIds)
			ofs << sheet_id << "\n";
		++sheet_id;
	}
}

void WriteDualVerticesVTK(const RefinedDualQuad& dual, const std::vector<std::unordered_set<size_t>> all_dualVIds, const char *filename) {
	size_t numOfLines = 0;
	for (const auto& dualEdgeIds : all_dualVIds)
		numOfLines += dualEdgeIds.size();

	const auto& V = dual.V;
	const auto& E = dual.E;
	std::ofstream ofs(filename);
	ofs << "# vtk DataFile Version 2.0\n"
		<< filename << "\n"
		<< "ASCII\n\n"
		<< "DATASET POLYDATA\n";
	ofs << "POINTS " << V.size() << " double" << "\n";
	for (const auto& v : V)
		ofs << v.x << " " << v.y << " " << v.z << "\n";

	ofs << "VERTICES " << numOfLines << " " << 2 * numOfLines << "\n";
	for (const auto& dualVIds : all_dualVIds)
		for (const auto id : dualVIds)
		ofs << "1 " << id << "\n";

	ofs << "CELL_DATA " << numOfLines << "\n"
		<< "SCALARS " << "id" << " int 1\n"
		<< "LOOKUP_TABLE default\n";
	size_t sheet_id = 0;
	for (const auto& dualVIds : all_dualVIds) {
		for (const auto id : dualVIds)
			ofs << sheet_id << "\n";
		++sheet_id;
	}
}

int main(int argc, char* argv[]) {
	if (argc < 3) {
		std::cout << "Usage: MapLoops orig.quad.m open.quad.m open.uv.quad.m output.quad.vtk\n";
		return -1;
	}
	std::vector<Vertex> origV;
	std::vector<Face> origF;
	read(argv[1], origV, origF);
	Mesh origMesh;
	generate_quad_mesh(origV, origF, origMesh);

	std::vector<VertexUV> openV;
	std::vector<Face> openF;
	read(argv[2], openV, openF);
	Mesh openMesh;
	generate_quad_mesh(openV, openF, openMesh);

	std::vector<Vertex> uvV;
	std::vector<Face> uvF;
	read(argv[3], uvV, uvF);
	Mesh uvMesh;
	generate_quad_mesh(uvV, uvF, uvMesh);
	
	origMesh.BuildAllConnectivities();
	origMesh.ExtractBoundary();
	origMesh.ExtractSingularities();
	origMesh.BuildParallelE();
	origMesh.BuildConsecutiveE();
	origMesh.BuildOrthogonalE();

	openMesh.BuildAllConnectivities();
	openMesh.ExtractBoundary();
	openMesh.ExtractSingularities();
	openMesh.BuildParallelE();
	openMesh.BuildConsecutiveE();
	openMesh.BuildOrthogonalE();

	uvMesh.BuildAllConnectivities();
	uvMesh.ExtractBoundary();
	uvMesh.ExtractSingularities();
	uvMesh.BuildParallelE();
	uvMesh.BuildConsecutiveE();
	uvMesh.BuildOrthogonalE();

	BaseComplexQuad baseComplex(origMesh);
	baseComplex.Build();
	BaseComplexSheetQuad baseComplexSheets(baseComplex);
	baseComplexSheets.Extract();
	//baseComplexSheets.WriteDualVTK("ChordDual.vtk");


	//auto open_boundary_eids = get_boundary_eids(openMesh);
	//auto open_boundary_vids = get_boundary_vids(openMesh);
	std::cout << "building key_edgeId...\n";
	auto key_edgeId = get_key_edgeId(origMesh);
	std::cout << "building key_faceId...\n";
	auto key_faceId = get_key_faceId(origMesh);
	std::cout << "building dual...\n";
	RefinedDualQuad dual(origMesh);
	dual.Build();
	auto dual_key_edgeId = get_key_edgeId(dual);
	auto dual_key_faceId = get_key_faceId(dual);

	std::cout << "building open_key_edgeId...\n";
	auto open_key_edgeId = get_key_edgeId(openMesh);
	std::cout << "building open_key_faceId...\n";
	auto open_key_faceId = get_key_faceId(openMesh);
	std::cout << "building open_dual...\n";
	RefinedDualQuad open_dual(openMesh);
	open_dual.Build();
	auto open_dual_key_edgeId = get_key_edgeId(open_dual);
	auto open_dual_key_faceId = get_key_faceId(open_dual);

	std::cout << "building uv_key_edgeId...\n";
	auto uv_key_edgeId = get_key_edgeId(uvMesh);
	std::cout << "building open_key_faceId...\n";
	auto uv_key_faceId = get_key_faceId(uvMesh);
	std::cout << "building open_dual...\n";
	RefinedDualQuad uv_dual(uvMesh);
	uv_dual.Build();
	auto uv_dual_key_edgeId = get_key_edgeId(open_dual);
	auto uv_dual_key_faceId = get_key_faceId(open_dual);

	std::cout << "building origVid_openVids...\n";
	std::unordered_map<size_t, std::set<size_t>> origVid_DualVids;
	for (auto& open_v : openMesh.V)
		origVid_DualVids[openV.at(open_v.id).father].insert(open_v.id);

	std::cout << "building origDualVid_openDualVids...\n";
	std::unordered_map<size_t, std::set<size_t>> origDualVid_openDualVids;
	for (auto& open_v : open_dual.V) {
		if (open_v.id < openMesh.V.size()) origDualVid_openDualVids[openV.at(open_v.id).father].insert(open_v.id);
		else if (open_v.id < openMesh.V.size() + openMesh.E.size()) {
			auto edgeid = open_v.id - openMesh.V.size();
			auto& e = openMesh.E.at(edgeid);
			auto key = (openV.at(e.Vids[0]).father << 32) | openV.at(e.Vids[1]).father;
			auto orig_eid = key_edgeId[key];
			auto orig_vid = origMesh.V.size() + orig_eid;
			origDualVid_openDualVids[orig_vid].insert(open_v.id);
		} else if (open_v.id < openMesh.V.size() + openMesh.E.size() + openMesh.F.size()) {
			auto faceid = open_v.id - openMesh.V.size() - openMesh.E.size();
			//auto& face = openMesh.F.at(faceid);
			//auto orig_face_key = make_facekey(openV.at(face.Vids[0]).father, openV.at(face.Vids[1]).father, 
			//	                              openV.at(face.Vids[2]).father, openV.at(face.Vids[3]).father);
			//auto orig_faceid = key_faceId[orig_face_key];
			auto orig_vid = origMesh.V.size() + origMesh.E.size() + faceid;
			origDualVid_openDualVids[orig_vid].insert(open_v.id);
		}
	}

	std::cout << "building origDualEid_openDualEid...\n";
	std::unordered_map<size_t, std::set<size_t>> origDualEid_openDualEid;
	for (auto& orig_e : dual.E) {
		if (origDualVid_openDualVids[orig_e.Vids[0]].size() == 1 && origDualVid_openDualVids[orig_e.Vids[1]].size() == 1) {
			auto open_vid0 = *origDualVid_openDualVids[orig_e.Vids[0]].begin();
			auto open_vid1 = *origDualVid_openDualVids[orig_e.Vids[1]].begin();
			auto open_key = (open_vid0 << 32) | open_vid1;
			auto open_eid = open_dual_key_edgeId[open_key];
			origDualEid_openDualEid[orig_e.id].insert(open_eid);
		}
		//if (origDualVid_openDualVids[orig_e.Vids[0]].size() == 1 && origDualVid_openDualVids[orig_e.Vids[1]].size() == 1) {
		//	auto open_vid0 = *origDualVid_openDualVids[orig_e.Vids[0]].begin();
		//	auto open_vid1 = *origDualVid_openDualVids[orig_e.Vids[1]].begin();
		//	if (open_vid0 < openMesh.V.size() && open_vid1 < openMesh.V.size()) {
		//		auto open_key = (open_vid0 << 32) | open_vid1;
		//		auto open_eid = open_dual_key_edgeId[open_key];
		//		origDualEid_openDualEid[orig_e.id].insert(open_eid);
		//	} /*else if (open_vid0 > openMesh.V.size() + openMesh.E.size()) {
		//		auto open_fid = open_vid1 - openMesh.V.size() - openMesh.E.size();
		//		auto& open_f = openMesh.F.at(open_fid);
		//		auto open_eid = open_vid0 - openMesh.V.size();
		//		origDualEid_openDualEid[orig_e.id].insert(open_eid);
		//	} else if (open_vid1 > openMesh.V.size() + openMesh.E.size()) {
		//		auto open_fid = open_vid1 - openMesh.V.size() - openMesh.E.size();
		//		auto& open_f = openMesh.F.at(open_fid);
		//		auto open_eid = open_vid0 - openMesh.V.size();
		//		origDualEid_openDualEid[orig_e.id].insert(open_eid);
		//	}*/
		//}
		//if (orig_e.id == 24084) {
		//	auto& open_vid0s = origDualVid_openDualVids[orig_e.Vids[0]];
		//	auto& open_vid1s = origDualVid_openDualVids[orig_e.Vids[1]];
		//	continue;
		//}
		//else if (origDualVid_openDualVids[orig_e.Vids[0]].size() == 1 && origDualVid_openDualVids[orig_e.Vids[1]].size() == 2) {
		//	auto open_vid0 = *origDualVid_openDualVids[orig_e.Vids[0]].begin();
		//	if (open_vid0 > openMesh.V.size() + openMesh.E.size()) {
		//		auto open_fid0 = open_vid0 - openMesh.V.size() - openMesh.E.size();
		//		auto& open_f0 = openMesh.F.at(open_fid0);
		//		
		//		for (auto vid : open_f0.Vids) {
		//			bool found = false;
		//			for (auto open_vid1 : origDualVid_openDualVids[orig_e.Vids[1]])
		//				if (vid == open_vid1) {
		//					auto open_eid = open_vid1 - openMesh.V.size();
		//					origDualEid_openDualEid[orig_e.id].insert(open_eid);
		//					found = true;
		//					break;
		//				}
		//			if (found) break;
		//		}
		//	}
		//} else if (origDualVid_openDualVids[orig_e.Vids[1]].size() == 1 && origDualVid_openDualVids[orig_e.Vids[0]].size() == 2) {
		//	auto open_vid0 = *origDualVid_openDualVids[orig_e.Vids[1]].begin();
		//	if (open_vid0 > openMesh.V.size() + openMesh.E.size()) {
		//		auto open_fid0 = open_vid0 - openMesh.V.size() - openMesh.E.size();
		//		auto& open_f0 = openMesh.F.at(open_fid0);

		//		for (auto vid : open_f0.Vids) {
		//			bool found = false;
		//			for (auto open_vid1 : origDualVid_openDualVids[orig_e.Vids[0]])
		//				if (vid == open_vid1) {
		//					auto open_eid = open_vid1 - openMesh.V.size();
		//					origDualEid_openDualEid[orig_e.id].insert(open_eid);
		//					found = true;
		//					break;
		//				}
		//			if (found) break;
		//		}
		//	}
		//} 
		//else if (origDualVid_openDualVids[orig_e.Vids[0]].size() == 1 && origDualVid_openDualVids[orig_e.Vids[1]].size() == 2) {
		//	auto open_vid0 = *origDualVid_openDualVids[orig_e.Vids[0]].begin();
		//	auto open_vid1_0 = *origDualVid_openDualVids[orig_e.Vids[1]].begin();
		//	auto open_vid1_1 = *origDualVid_openDualVids[orig_e.Vids[1]].rbegin();
		//	if (open_vid0 > openMesh.V.size() + openMesh.E.size()) {
		//		auto open_fid0 = open_vid0 - openMesh.V.size() - openMesh.E.size();
		//		auto& open_f0 = openMesh.F.at(open_fid0);

		//		for (auto vid : open_f0.Vids)
		//			if (vid == open_vid1_0) {
		//				auto open_eid = open_vid1_0 - openMesh.V.size();
		//				origDualEid_openDualEid[orig_e.id].insert(open_eid);
		//				break;
		//			} else if (vid == open_vid1_1) {
		//				auto open_eid = open_vid1_1 - openMesh.V.size();
		//				origDualEid_openDualEid[orig_e.id].insert(open_eid);
		//				break;
		//			}
		//	}
		//}
		//else if (origDualVid_openDualVids[orig_e.Vids[1]].size() == 1 && origDualVid_openDualVids[orig_e.Vids[0]].size() == 2) {
		//	auto open_vid0 = *origDualVid_openDualVids[orig_e.Vids[1]].begin();
		//	auto open_vid1_0 = *origDualVid_openDualVids[orig_e.Vids[0]].begin();
		//	auto open_vid1_1 = *origDualVid_openDualVids[orig_e.Vids[0]].rbegin();
		//	if (open_vid0 > openMesh.V.size() + openMesh.E.size()) {
		//		auto open_fid0 = open_vid0 - openMesh.V.size() - openMesh.E.size();
		//		auto& open_f0 = openMesh.F.at(open_fid0);

		//		for (auto vid : open_f0.Vids)
		//			if (vid == open_vid1_0) {
		//				auto open_eid = open_vid1_0 - openMesh.V.size();
		//				origDualEid_openDualEid[orig_e.id].insert(open_eid);
		//				break;
		//			} else if (vid == open_vid1_1) {
		//				auto open_eid = open_vid1_1 - openMesh.V.size();
		//				origDualEid_openDualEid[orig_e.id].insert(open_eid);
		//				break;
		//			}
		//	}
		//}
	}

	std::cout << "building all_dualEdgeIds...\n";
	std::vector<std::unordered_set<size_t>> all_dualEdgeIds;
	for (size_t i = 0; i < baseComplexSheets.sheets_componentFaceIds.size(); ++i) {
		std::unordered_set<size_t> dualEdgeIds = baseComplexSheets.GetDualEdgeIds(i);
		std::unordered_set<size_t> open_dualEdgeIds;
		for (auto dualEdgeId : dualEdgeIds)
			if (origDualEid_openDualEid[dualEdgeId].size() == 1) 
				open_dualEdgeIds.insert(*origDualEid_openDualEid[dualEdgeId].begin());
		all_dualEdgeIds.push_back(open_dualEdgeIds);
	}

	std::cout << "writing all_dualEdgeIds...\n";
	//WriteDualVTK(open_dual, all_dualEdgeIds, argv[3]);
	WriteDualVTK(uv_dual, all_dualEdgeIds, argv[4]);

	//std::vector<std::unordered_set<size_t>> all_dualVIds;
	//for (size_t i = 0; i < baseComplexSheets.sheets_componentFaceIds.size(); ++i) {
	//	std::unordered_set<size_t> dualEdgeIds = baseComplexSheets.GetDualEdgeIds(i);
	//	std::unordered_set<size_t> open_dualVIds;
	//	for (auto dualEdgeId : dualEdgeIds) {
	//		auto& dualEdge = dual.E.at(dualEdgeId);
	//		auto& v0 = origDualVid_openDualVids[dualEdge.Vids[0]];
	//		auto& v1 = origDualVid_openDualVids[dualEdge.Vids[1]];
	//		open_dualVIds.insert(v0.begin(), v0.end());
	//		open_dualVIds.insert(v1.begin(), v1.end());
	//	}
	//	all_dualVIds.push_back(open_dualVIds);
	//}
	//std::cout << "writing all_dualVIds...\n";
	//WriteDualVerticesVTK(open_dual, all_dualVIds, argv[3]);
	return 0;
}
