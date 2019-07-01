/*
 * ExtractDualMesh.cpp
 *
 *  Created on: Nove 27, 2018
 *      Author: cotrik
 */

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "ArgumentManager.h"
#include "DualMesh.h"
#include <iostream>
#include <fstream>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <stack>
#include <algorithm>


int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: ExtractDualMesh <file> <out.vtk> removedNonQuad=false";
        return -1;
    }
    ArgumentManager argumentManager(argc, argv);
	bool removedNonQuad = false;
	std::string strremovedNonQuad = argumentManager.get("removedNonQuad");
	if (!strremovedNonQuad.empty()) removedNonQuad = strremovedNonQuad != "true" ? false : true;
	MeshFileReader reader(argv[1]);
	auto mesh = reader.GetMesh();
	mesh.BuildAllConnectivities();
	mesh.ExtractBoundary();
	mesh.ExtractSingularities();
	mesh.BuildParallelE();
	mesh.BuildConsecutiveE();
	mesh.BuildOrthogonalE();

	DualMesh dualMesh;
	dualMesh.Build(mesh);
	if (dualMesh.m_cellType == POLYGON)
		dualMesh.BuildAllConnectivities();
	else if (dualMesh.m_cellType == POLYHEDRA)
		dualMesh.BuildConnection(mesh);
	//dualMesh.ExtractBoundary();
	//dualMesh.ExtractSingularities();
	//dualMesh.BuildParallelE();
	//dualMesh.BuildConsecutiveE();
	//dualMesh.BuildOrthogonalE();
	//dualMesh.FixOrientation(mesh);
	//dualMesh.FixOrientation();

	//std::set<size_t> canceled_fids;
	//for (auto& f : dualMesh.F)
	//	if (f.Vids.size() != 4) canceled_fids.insert(f.id);
	//std::cout << "canceled_fids.size() = " << canceled_fids.size() << std::endl;
	//std::vector<size_t> fids;
	//for (auto& f : dualMesh.F)
	//	if (canceled_fids.find(f.id) == canceled_fids.end()) fids.push_back(f.id);
	//std::cout << "fids.size() = " << fids.size() << std::endl;

	MeshFileWriter writer(dualMesh, argv[2]);
	if (!removedNonQuad) writer.WriteFile();
	//else writer.WriteFacesVtk(fids);
    return 0;
}

