/*
* Tri2Quad.cpp
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

void write(const Mesh& mesh, const char* filename) {
	std::ofstream ofs(filename);
	ofs << "# vtk DataFile Version 4.0" << std::endl
		<< filename << std::endl
		<< "ASCII" << std::endl << std::endl
		<< "DATASET POLYDATA" << std::endl;
	ofs << "POINTS " << mesh.V.size() << " double" << std::endl;
	for (auto& v : mesh.V)
		ofs << v.x << " " << v.y << " " << v.z << std::endl;
	ofs << "POLYGONS " << mesh.C.size() / 2 << " " << mesh.C.size() / 2 * 5 << "\n";
	for (size_t i = 0; i < mesh.C.size(); i += 2) {
		auto& c0 = mesh.C.at(i);
		auto& c1 = mesh.C.at(i + 1);
		ofs << "4 " << c0.Vids[0] << " " << c0.Vids[1] << " " << c0.Vids[2] << " " << c1.Vids[2] << std::endl;
	}
}

int main(int argc, char* argv[]) {
	if (argc < 3) {
		std::cout << "Usage: Quad2Tri quad.[vtk|off|mesh|obj] tri.vtk\n";
		return -1;
	}

	MeshFileReader reader(argv[1]);
	auto mesh = (Mesh&)reader.GetMesh();
	//mesh.RemoveUselessVertices();
	write(mesh, argv[2]);
	return 0;
}
