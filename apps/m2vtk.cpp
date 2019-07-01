/*
* ExtractQuadMesh.cpp
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
			//std::string father, uv, vv;
			//ss >> father >> uv >> vv;

			//for (auto& c : father)
			//	if (c == '(' || c == ')') c = ' ';
			//std::istringstream ss_father(father);
			//Vertex u_v;
			//ss_father >> father >> v.father;
			//--v.father;

			//for (auto& c : uv)
			//	if (c == '(' || c == ')') c = ' ';
			//std::istringstream ss_uv(uv);
			//ss_uv >> uv >> v.uv.x;

			//for (auto& c : vv)
			//	if (c == '(' || c == ')') c = ' ';
			//std::istringstream ss_vv(vv);
			//ss_vv >> v.uv.y;

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

int main(int argc, char* argv[]) {
	if (argc < 3) {
		std::cout << "Usage: m2vtk input.m output.[vtk|off|mesh]\n";
		return -1;
	}

	std::vector<VertexUV> V;
	std::vector<Face> F;
	read(argv[1], V, F);
	write(argv[2], V, F);
	return 0;
}
