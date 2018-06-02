/*
 * HexRefine.cpp
 *
 *  Created on: Dec 5, 2016
 *      Author: cotrik
 */
#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "RefinedDual.h"
#include "ArgumentManager.h"

int main(int argc, char* argv[])
{
    if (argc < 2) {
        std::cout << "Usage: HexRefine hex.vtk refine.vtk" << std::endl;
        return -1;
    }
    ArgumentManager argumentManager(argc, argv);
    std::string input = argv[1];
    std::string output = argv[2];
    std::cout << "---------------------------------------" << std::endl;
    std::cout << "input  = " << input << std::endl;
    std::cout << "output = " << output << std::endl;
    std::cout << "---------------------------------------" << std::endl;

    MeshFileReader reader(input.c_str());
    Mesh& mesh = (Mesh&) reader.GetMesh();
    mesh.BuildAllConnectivities();
    mesh.ExtractBoundary();
    mesh.ExtractSingularities();
    mesh.BuildParallelE();
    mesh.BuildConsecutiveE();
    mesh.BuildOrthogonalE();

    auto refinedMesh = GetRefineMesh3(mesh);
    MeshFileWriter writer(refinedMesh, output.c_str());
    writer.WriteFile();
}


