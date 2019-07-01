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
	for (auto& v : mesh.V)
		ofs << "Vertex " << v.id + 1 << " " << v.x << " " << v.y << " " << v.z << "\n";
	size_t fid = 0;
	for (auto& f : mesh.C) { 
		ofs << "Face " << ++fid;
		for (auto vid : f.Vids)
			ofs << " " << vid + 1;
		ofs << std::endl;
	}

}
int main(int argc, char* argv[]) {
	if (argc < 3) {
		std::cout << "Usage: m2vtk input.[vtk|off|mesh|obj] output.m\n";
		return -1;
	}

	MeshFileReader reader(argv[1]);
	auto mesh = (Mesh&)reader.GetMesh();
	write(mesh, argv[2]);
	return 0;
}
