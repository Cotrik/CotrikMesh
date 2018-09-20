/*
 * CleanStructure.cpp
 *
 *  Created on: Aug 28, 2018
 *      Author: cotrik
 */

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "BaseComplexCleaner.h"
#include "ArgumentManager.h"

#include <iostream>
#include <fstream>
#include <algorithm>

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: CleanStructure hex.vtk output=out.vtk" << std::endl;
        return -1;
    }
    ArgumentManager argumentManager(argc, argv);
    std::string filename = argv[1];
    std::string output = "out.vtk";
    if (!argumentManager.get("output").empty()) output = argumentManager.get("output");
    std::cout << "---------------------------------------" << std::endl;
    std::cout << "filename = " << filename << std::endl;
    std::cout << "output = " << output << std::endl;
    std::cout << "---------------------------------------" << std::endl;

    MeshFileReader reader(filename.c_str());
    Mesh& mesh = (Mesh&)reader.GetMesh();
    mesh.BuildAllConnectivities();
    //mesh.RemoveUselessVertices();
    mesh.ExtractBoundary();
    mesh.ExtractSingularities();
    //"Info: cos10 = 0.984807753; cos15 = 0.965925826; cos20 = 0.939692621; cos25 = 0.906307787; cos30 = 0.866025404\n\n";
    const double cosangle = 0.866025404;
    mesh.SetCosAngleThreshold(cosangle);
    mesh.LabelSurface();
    mesh.LabelSharpEdges(true);
    // For extracting singularity Graph
    mesh.BuildParallelE();
    mesh.BuildConsecutiveE();
    mesh.BuildOrthogonalE();

    BaseComplex baseComplex(mesh);
    baseComplex.Build();

    BaseComplexCleaner baseComplexCleaner(baseComplex);
    baseComplexCleaner.Run();

    return 0;
}


