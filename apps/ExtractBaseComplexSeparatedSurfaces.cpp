/*
 * ExtractBaseComplex.cpp
 *
 *  Created on: Jun 1, 2017
 *      Author: cotrik
 */

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "EdgeLines.h"
#include "BaseComplex.h"
#include "BaseComplexQuad.h"
#include "BaseComplexSheet.h"
#include "BaseComplexSheetQuad.h"
#include "BaseComplexChord.h"
#include "BaseComplexEditor.h"
#include "SingularityGraph.h"
#include "ArgumentManager.h"

#include <iostream>

int main(int argc, char* argv[]){
    if (argc < 2) {
        std::cout << "Usage: ExtractBaseComplexSeparatedSurfaces hex.vtk" << std::endl;
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
    baseComplex.WriteBaseComplexHexVTK("BaseComplexHex.vtk");
    baseComplex.WriteBaseComplex_VTK("BaseComplex.vtk");
    baseComplex.WriteSingularV_VTK("singularV.vtk");
    baseComplex.WriteSingularE_VTK("singularE.vtk");
    baseComplex.WriteSingularities_VTK("singularities.vtk");
	baseComplex.WriteBaseComplexSeparatedFacePatchesVTK("BaseComplexSeparatedFacePatches.vtk");
	baseComplex.WriteBaseComplexSeparatedSurfacesVTK("BaseComplexSeparatedSurfaces.vtk");
    return 0;
}
