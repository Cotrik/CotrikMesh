/*
 * ExtractSingularities.cpp
 *
 *  Created on: Sep 9, 2018
 *      Author: davim
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

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: ExtractSingularities hex.vtk" << std::endl;
        return -1;
    }
    ArgumentManager argumentManager(argc, argv);
    std::string filename = argv[1];
    std::cout << "---------------------------------------" << std::endl;
    std::cout << "filename = " << filename << std::endl;
    std::cout << "---------------------------------------" << std::endl;

    MeshFileReader reader(filename.c_str());
    Mesh& mesh = (Mesh&)reader.GetMesh();
    mesh.RemoveUselessVertices();
    mesh.BuildAllConnectivities();
    mesh.ExtractBoundary();
    mesh.ExtractSingularities();
    //"Info: cos10 = 0.984807753; cos15 = 0.965925826; cos20 = 0.939692621; cos25 = 0.906307787; cos30 = 0.866025404\n\n";
    //const double cosangle = 0.866025404;
    //mesh.SetCosAngleThreshold(cosangle);
    //mesh.LabelSurface();
    //mesh.LabelSharpEdges(true);
    std::cout << "mesh cellType = " << mesh.m_cellType << "\n";
    if (mesh.m_cellType == HEXAHEDRA) {
        // For extracting singularity Graph
        mesh.BuildParallelE();
        mesh.BuildConsecutiveE();
        mesh.BuildOrthogonalE();

        BaseComplex baseComplex(mesh);
        baseComplex.Build();

        std::cout << "#Singularites = " << baseComplex.SingularityI.E.size() << "\n";

        baseComplex.WriteBaseComplexHexVTK("BaseComplexHex.vtk");
        baseComplex.WriteBaseComplex_VTK("BaseComplex.vtk");
        baseComplex.WriteSingularV_VTK("singularV.vtk");
        baseComplex.WriteSingularE_VTK("singularE.vtk");
        baseComplex.WriteSingularities_VTK("singularities.vtk");
        baseComplex.WriteBaseComplexComponentsVTK("BaseComplexComponents.vtk");
        baseComplex.WriteBaseComplex_ColorFacesVTK("BaseComplexColorFaces.vtk");
        baseComplex.WriteBaseComplex_ColorEdgesVTK("BaseComplexColorEdges.vtk");
        baseComplex.WriteBaseComplex_ColorVerticesVTK("BaseComplexColorVertices.vtk");

    //    baseComplex.WriteBaseComplexAllComponentsVTK("ComponentCells");

        BaseComplexSheet baseComplexSheet(baseComplex);
        baseComplexSheet.Extract();

        std::cout << "#Sheets = " << baseComplexSheet.sheets_componentCellIds.size() << "\n";
        std::cout << "#Singularities = " << baseComplex.SingularityI.E.size() << "\n";
    } else if (mesh.m_cellType == QUAD) {
        std::vector<size_t> singularVids;
        for (auto& v : mesh.V)
            if (v.isSingularity) singularVids.push_back(v.id);
        MeshFileWriter writer(mesh, "singular_V.vtk");
        writer.WriteVerticesVtk(singularVids);
        std::cout << "#Singularities = " << singularVids.size() << "\n";
    }

    return 0;
}
