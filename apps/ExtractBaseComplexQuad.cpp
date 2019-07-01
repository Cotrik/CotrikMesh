/*
 * ExtractBaseComplex.cpp
 *
 *  Created on: Jun 1, 2017
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

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cout << "Usage: ExtractBaseComplexQuad quad.vtk" << std::endl;
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
    // For extracting singularity Graph
    mesh.BuildParallelE();
    mesh.BuildConsecutiveE();
    mesh.BuildOrthogonalE();

//    EdgeLines edgeLines(mesh);
//    edgeLines.Build();
//    for (auto& edgeLine : edgeLines.edgeLines) {
//        for (size_t i = 0; i < edgeLine.Eids.size(); ++i) {
//            auto& edge = mesh.E.at(edgeLine.Eids.at(i));
//            auto vid0 = edgeLine.Vids.at(i + 0);
//            auto vid1 = edgeLine.Vids.at(i + 1);
//            if (edge.Vids[0] != vid0 && edge.Vids[1] != vid1 && edge.Vids[0] == vid1 && edge.Vids[1] == vid0) {
//                std::swap(edge.Vids[0], edge.Vids[1]);
//                //std::cout << "std::swap(edge.Vids[0], edge.Vids[1]);\n";
//            }
//        }
//    }

    BaseComplexQuad baseComplex(mesh);
    baseComplex.Build();
    baseComplex.WriteBaseComplexQuadVTK("BaseComplexQuad.vtk");
//    baseComplex.WriteBaseComplex_VTK("BaseComplex.vtk");
//    baseComplex.WriteSingularV_VTK("singularV.vtk");
//    baseComplex.WriteSingularE_VTK("singularE.vtk");
//    baseComplex.WriteSingularities_VTK("singularities.vtk");
    baseComplex.WriteBaseComplexSeparatedEdgeLinksVTK("BaseComplexSeparatedEdgeLinks.vtk");
//    baseComplex.WriteComponentEdge_NeighborComponentFaces_VTK("ComponentEdge_NeighborComponentFaces.vtk");
//    baseComplex.WriteSingularEdge_NeighborSeparatedComponentFacePatches_VTK("SingularEdge_NeighborSeparatedComponentFacePatches.vtk");
//    baseComplex.WriteAllSingularEdge_NeighborSeparatedComponentFacePatches_VTK("SingularFaces");
//    baseComplex.WriteBaseComplexComponentsWithoutSingularitiesVTK("BaseComplexComponentsWithoutSingularities.vtk");
//    baseComplex.WriteBaseComplexComponentsWithSingularitiesVTK("BaseComplexComponentsWithSingularities.vtk");
//    baseComplex.WriteBaseComplexAllComponentsVTK("ComponentFaces");
    baseComplex.WriteBaseComplexComponentsVTK("BaseComplexComponents.vtk");
//    baseComplex.WriteBaseComplex_ColorFacesVTK("BaseComplexColorFaces.vtk");
//    baseComplex.WriteBaseComplex_ColorEdgesVTK("BaseComplexColorEdges.vtk");
//    baseComplex.WriteBaseComplex_ColorVerticesVTK("BaseComplexColorVertices.vtk");

    BaseComplexSheetQuad baseComplexSheets(baseComplex);
    baseComplexSheets.Extract();
////    baseComplexSheets.ExtractSets();
////    baseComplexSheets.WriteAllSheetsCellsVTK("SheetCells");
//    baseComplexSheets.WriteAllSheetsFacesVTK("SheetFaces");
//    baseComplexSheets.WriteAllSheetsEdgesVTK("SheetEdges");
////    baseComplexSheets.WriteAllSheetsFacesAndEdgesVTK("SheetFacesAndEdges");
	baseComplexSheets.WriteAllSheetsFacesInOneVTK("SheetFaces.vtk");
//    baseComplexSheets.WriteAllSheetsFacesDualVTK("SheetDual");
	baseComplexSheets.WriteDualVTK("ChordDual.vtk");
//	baseComplexSheets.WriteDualLinksVTK("ChordDualLinks.vtk");
//    baseComplexSheets.ExtractSheetDecompositions();
//    baseComplexSheets.WriteSheetDecompositionsFile(NULL);
//
//    baseComplexSheets.ExtractSheetConnectivities();
//    baseComplexSheets.WriteSheetsConnectivitiesMatrixVTK("SheetsConnectivities.vtk");
//    baseComplexSheets.WriteSheetsConnectivitiesMatrixMat("SheetsConnectivities.mat");
//    baseComplexSheets.ExtractMainSheetConnectivities();
//    baseComplexSheets.WriteSheetsConnectivitiesMatrixVTK("DominantSheetsConnectivities.vtk");
//    baseComplexSheets.WriteSheetsConnectivitiesMatrixMat("DominantSheetsConnectivities.mat");
//    baseComplexSheets.ComputeComplexity();

//    BaseComplexChord baseComplexChords(baseComplex);
//    baseComplexChords.Extract();
//    baseComplexChords.WriteAllChordsCellsVTK("ChordCells");
//    baseComplexChords.WriteAllChordsFacesVTK("ChordFaces");
//    baseComplexChords.WriteAllChordsEdgesVTK("ChordEdges");
//    baseComplexChords.WriteAllChordsFacesAndEdgesVTK("ChordFacesAndEdges");
//    baseComplexChords.WriteAllChordsCurvesVTK("ChordCurves");
    return 0;
}


