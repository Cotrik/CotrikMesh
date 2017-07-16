/*
 * ExtractBaseComplex.cpp
 *
 *  Created on: Jun 1, 2017
 *      Author: cotrik
 */

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "BaseComplex.h"
#include "ArgumentManager.h"

#include <iostream>

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cout << "Usage: ExtractBaseComplex hex.vtk" << std::endl;
        return -1;
    }
    ArgumentManager argumentManager(argc, argv);
    std::string filename = argv[1];
    std::cout << "---------------------------------------" << std::endl;
    std::cout << "filename = " << filename << std::endl;
    std::cout << "---------------------------------------" << std::endl;

    MeshFileReader reader(filename.c_str());
    Mesh& mesh = (Mesh&)reader.GetMesh();
    mesh.BuildAllConnectivities();
    mesh.ExtractBoundary();
    mesh.ExtractSingularities();

    BaseComplex baseComplex(mesh);
    baseComplex.Build();
    baseComplex.WriteBaseComplex_VTK("BaseComplex.vtk");
    baseComplex.WriteSingularV_VTK("singularV.vtk");
    baseComplex.WriteSingularE_VTK("singularE.vtk");
    baseComplex.WriteSingularities_VTK("singularities.vtk");
    baseComplex.WriteBaseComplexSeparatedFacePatchesVTK("BaseComplexSeparatedFacePatches.vtk");
    baseComplex.WriteComponentEdge_NeighborComponentFaces_VTK("ComponentEdge_NeighborComponentFaces.vtk");
    return 0;
}


