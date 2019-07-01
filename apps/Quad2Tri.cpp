/*
* Quad2Tri.cpp
*
*  Created on: Oct 20, 2018
*      Author: cotrik
*/

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
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

Mesh quad2tri(const Mesh& mesh);
void write(const Mesh& mesh, const char* filename) {
	std::ofstream ofs(filename);
	ofs << "# vtk DataFile Version 4.0" << std::endl
		<< filename << std::endl
		<< "ASCII" << std::endl << std::endl
		<< "DATASET POLYDATA" << std::endl;
	ofs << "POINTS " << mesh.V.size() << " double" << std::endl;
	for (auto& v : mesh.V)
		ofs << v.x << " " << v.y << " " << v.z << std::endl;
	ofs << "POLYGONS " << mesh.C.size() * 2 << " " << mesh.C.size() * 8 << "\n";
	for (auto& c : mesh.C) {
		ofs << "3 " << c.Vids[0] << " " << c.Vids[1] << " " << c.Vids[2] << std::endl;
		ofs << "3 " << c.Vids[0] << " " << c.Vids[2] << " " << c.Vids[3] << std::endl;
	}
}

int main(int argc, char* argv[]) {
	if (argc < 3) {
		std::cout << "Usage: Quad2Tri quad.[vtk|off|mesh|obj] tri.vtk\n";
		return -1;
	}

	MeshFileReader reader(argv[1]);
	auto mesh = (Mesh&)reader.GetMesh();
	mesh.BuildAllConnectivities();
	write(mesh, argv[2]);
//	auto triMesh = quad2tri(mesh);
//	MeshFileWriter writer(triMesh, argv[2]);
//	writer.WriteFile();
	return 0;
}

glm::dvec3 get_center(const Mesh& mesh, const Face& f) {
    glm::dvec3 res;
    for (auto vid : f.Vids)
        res += mesh.V.at(vid).xyz();
    res /= f.Vids.size();
    return res;
}

Mesh quad2tri(const Mesh& mesh) {
    Mesh triMesh;
    triMesh.m_cellType = TRIANGLE;
    triMesh.V = mesh.V;
    triMesh.V.resize(mesh.V.size() + mesh.F.size());
    triMesh.C.resize(4 * mesh.F.size());
    for (auto& f : mesh.F) {
        auto& triV = triMesh.V.at(mesh.V.size() + f.id);
        triV = get_center(mesh, f);
        triV.id = mesh.V.size() + f.id;
    }

    for (auto& f : mesh.F) {
        auto& centerV = triMesh.V.at(mesh.V.size() + f.id);
        size_t fid = 0;
        for (auto eid : f.Eids) {
            auto& e = mesh.E.at(eid);
            auto& triF = triMesh.C.at(f.id * 4 + fid);
            triF.id = f.id * 4 + fid++;
            triF.Vids = {centerV.id, e.Vids[0], e.Vids[1]};
        }
    }

    return triMesh;
}
