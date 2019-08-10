/*
 * ExtractBaseMesh.cpp
 *
 *  Created on: July 4, 2019
 *      Author: cotrik
 */

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "EdgeLines.h"
#include "BaseComplexQuad.h"
#include "BaseComplexSheetQuad.h"
#include "BaseComplexChord.h"
#include "BaseComplexEditor.h"
#include "SingularityGraph.h"
#include "ArgumentManager.h"

#include <iostream>

int main(int argc, char* argv[]) {
	if (argc < 2) {
		std::cout << "Usage: ExtractBaseMesh input.vtk output.vtk" << std::endl;
		return -1;
	}
	ArgumentManager argumentManager(argc, argv);
	std::string filename = argv[1];
	std::string output_filename = argv[2];
	std::cout << "---------------------------------------" << std::endl;
	std::cout << "input_filename = " << filename << std::endl;
	std::cout << "output_filename = " << output_filename << std::endl;
	std::cout << "---------------------------------------" << std::endl;

	MeshFileReader reader(filename.c_str());
	Mesh& mesh = (Mesh&)reader.GetMesh();
	mesh.RemoveUselessVertices();
	mesh.BuildAllConnectivities();
	mesh.ExtractBoundary();
	mesh.ExtractSingularities();
	mesh.BuildParallelE();
	mesh.BuildConsecutiveE();
	mesh.BuildOrthogonalE();

	if (mesh.m_cellType == QUAD) {
		BaseComplexQuad baseComplex(mesh);
		baseComplex.Build();
		baseComplex.WriteBaseComplexQuadVTK(output_filename.c_str());
	} else if (mesh.m_cellType == HEXAHEDRA) {
		BaseComplex baseComplex(mesh);
		baseComplex.Build();
		baseComplex.WriteBaseComplexHexVTK(output_filename.c_str());
	}
	return 0;
}


