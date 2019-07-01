/*
* ExtractBaseHomology.cpp
*
*  Created on: Oct 20, 2018
*      Author: cotrik
*/

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "ArgumentManager.h"
#include <iostream>
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

void read(const char* filename, std::vector<VertexUV>& V, std::vector<Face>& F, std::vector<Edge>& E) {
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
			Edge e(2);
			ss >> e.Vids[0] >> e.Vids[1];
			--e.Vids[0];
			--e.Vids[1];
			E.push_back(e);
		}
	}
}

void write(const char* filename, const std::vector<VertexUV>& V, const std::vector<Face>& F) {
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

void generate_tri_mesh(const std::vector<VertexUV>& V, const std::vector<Face>& F, Mesh& triMesh) {
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

	triMesh.V = triV;
	triMesh.C = C;
	triMesh.m_cellType = TRIANGLE;
	triMesh.BuildAllConnectivities();
	triMesh.ExtractBoundary();
}

void generate_tri_mesh(const std::vector<Vertex>& V, const std::vector<Face>& F, Mesh& triMesh) {
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

	triMesh.V = triV;
	triMesh.C = C;
	triMesh.m_cellType = TRIANGLE;
	triMesh.BuildAllConnectivities();
	triMesh.ExtractBoundary();
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
	quadMesh.BuildAllConnectivities();
	quadMesh.ExtractBoundary();
}

void write(const char* filename, const std::vector<Vertex>& V, const std::vector<Face>& F) {
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

std::vector<std::vector<size_t>> GetBaseEdgeIds(const Mesh& mesh, const size_t numOfBase) {
	std::unordered_map<size_t, size_t> key_edgeId;
	for (auto& e : mesh.E) {
		key_edgeId[(e.Vids[0] << 32) | e.Vids[1]] = e.id;
		key_edgeId[(e.Vids[1] << 32) | e.Vids[0]] = e.id;
	}

	std::vector<std::vector<size_t>> baseEdgeIds;
	for (auto i = 0; i < numOfBase; ++i) {
		std::vector<VertexUV> V;
		std::vector<Face> F;
		std::vector<Edge> E;
		std::string filename = std::string("base_") + std::to_string(i) + ".m";
		read(filename.c_str(), V, F, E);
		std::vector<size_t> eids;
		for (auto& e : E)
			eids.push_back(key_edgeId[(e.Vids[0] << 32) | e.Vids[1]]);
		baseEdgeIds.push_back(eids);
	}
	return baseEdgeIds;
}

std::vector<std::vector<size_t>> GetBaseLinkVids(const Mesh& mesh, const std::vector<std::vector<size_t>>& baseEdgeIds) {
	std::vector<std::vector<size_t>> baseLinkVids;
	for (auto& eids_ : baseEdgeIds) {
		auto eids = eids_;
		std::vector<size_t> ringEids;
		ringEids.push_back(eids.back());
		eids.pop_back();
		while (ringEids.size() < eids_.size()) {
			auto last_eid = ringEids.back();
			int i = 0;
			for (auto next_eid : eids) {
				auto& last_e = mesh.E.at(last_eid);
				std::set<size_t> vids(last_e.Vids.begin(), last_e.Vids.end());
				auto& next_e = mesh.E.at(next_eid);
				vids.insert(next_e.Vids.begin(), next_e.Vids.end());
				if (vids.size() == 3) {
					ringEids.push_back(next_eid);
					eids.erase(eids.begin() + i);
					break;
				}
				++i;
			}
		}
		std::vector<size_t> ringVids;
		auto& front_e = mesh.E.at(ringEids.front());
		ringVids.push_back(front_e.Vids[0]);
		for (auto eid : ringEids) {
			auto last_vid = ringVids.back();
			auto& e = mesh.E.at(eid);
			auto next_vid = e.Vids[0] == last_vid ? e.Vids[1] : e.Vids[0];
			ringVids.push_back(next_vid);
		}
		baseLinkVids.push_back(ringVids);
	}
	return baseLinkVids;
}

void Clean(std::vector<std::vector<size_t>>& baseLinkVids) {
	for (auto& link : baseLinkVids) {
		for (int i = 1; i < link.size() - 1; ++i) {
			auto& vid_curr = link[i];
			auto& vid_prev = link[i - 1];
			auto& vid_next = link[i + 1];
		}
	}
}
int main(int argc, char* argv[]) {
	if (argc < 4) {
		std::cout << "Usage: ExtractBaseHomology input.vtk output.vtk num_of_base\n";
		return -1;
	}
	ArgumentManager argumentManager(argc, argv);
	MeshFileReader reader(argv[1]);
	auto mesh = reader.GetMesh();
	mesh.BuildAllConnectivities();
	mesh.ExtractBoundary();
	mesh.ExtractSingularities();
	mesh.BuildParallelE();
	mesh.BuildConsecutiveE();
	mesh.BuildOrthogonalE();

	size_t numOfBase = std::stoi(argv[3]);
	
	std::vector<std::vector<size_t>> baseEdgeIds = GetBaseEdgeIds(mesh, numOfBase);
	std::vector<std::vector<size_t>> baseLinkVids = GetBaseLinkVids(mesh, baseEdgeIds);
	Clean(baseLinkVids);

	MeshFileWriter writer(mesh, argv[2]);
	writer.WriteLinksVtk(baseLinkVids);
	return 0;
}
