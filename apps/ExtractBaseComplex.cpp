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
#include "BaseComplexSheet.h"
#include "BaseComplexChord.h"
#include "BaseComplexEditor.h"
#include "SingularityGraph.h"
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
    //"Info: cos10 = 0.984807753; cos15 = 0.965925826; cos20 = 0.939692621; cos25 = 0.906307787; cos30 = 0.866025404\n\n";
    const double cosangle = 0.866025404;
    mesh.SetCosAngleThreshold(cosangle);
    mesh.LabelSurface();
    mesh.LabelSharpEdges(true);
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

    BaseComplex baseComplex(mesh);
    baseComplex.Build();
    baseComplex.WriteBaseComplexHexVTK("BaseComplexHex.vtk");
    baseComplex.WriteBaseComplex_VTK("BaseComplex.vtk");
    baseComplex.WriteSingularV_VTK("singularV.vtk");
    baseComplex.WriteSingularE_VTK("singularE.vtk");
    baseComplex.WriteSingularities_VTK("singularities.vtk");
    baseComplex.WriteBaseComplexSeparatedFacePatchesVTK("BaseComplexSeparatedFacePatches.vtk");
    baseComplex.WriteComponentEdge_NeighborComponentFaces_VTK("ComponentEdge_NeighborComponentFaces.vtk");
    baseComplex.WriteSingularEdge_NeighborSeparatedComponentFacePatches_VTK("SingularEdge_NeighborSeparatedComponentFacePatches.vtk");
    baseComplex.WriteBaseComplexComponentsWithoutSingularitiesVTK("BaseComplexComponentsWithoutSingularities.vtk");
    baseComplex.WriteBaseComplexComponentsWithSingularitiesVTK("BaseComplexComponentsWithSingularities.vtk");

    BaseComplexSheet baseComplexSheets(baseComplex);
    baseComplexSheets.Extract();
    baseComplexSheets.WriteAllSheetsCellsVTK("SheetCells");
    baseComplexSheets.WriteAllSheetsFacesVTK("SheetFaces");
    baseComplexSheets.WriteAllSheetsEdgesVTK("SheetEdges");
    baseComplexSheets.WriteAllSheetsFacesAndEdgesVTK("SheetFacesAndEdges");
    baseComplexSheets.WriteAllSheetsCellsDualVTK("SheetDual");

    BaseComplexChord baseComplexChords(baseComplex);
    baseComplexChords.Extract();
    baseComplexChords.WriteAllChordsCellsVTK("ChordCells");
    baseComplexChords.WriteAllChordsFacesVTK("ChordFaces");
    baseComplexChords.WriteAllChordsEdgesVTK("ChordEdges");
    baseComplexChords.WriteAllChordsFacesAndEdgesVTK("ChordFacesAndEdges");
    baseComplexChords.WriteAllChordsCurvesVTK("ChordCurves");
    //baseComplexChords.WriteAllChordsFramesVTK("ChordFrames");

//    BaseComplexEditor baseComplexEditor(mesh, baseComplex);
//    baseComplexEditor.Run();
//
//    SingularityGraph singularityGraph(baseComplex);
//    singularityGraph.Build();
//    singularityGraph.WriteMatrixVTK("matrix.directlyLinkedSingularEdge.vtk", singularityGraph.directlyLinkedSingularEdgeIds);
//    singularityGraph.WriteMatrixVTK("matrix.linkedByOneComponentEdgeSingularEdge.vtk", singularityGraph.linkedByOneComponentEdgeSingularEdgeIds);
//    singularityGraph.WriteMatrixVTK("matrix.orthogonalDirectionSingularEdgeIds.vtk", singularityGraph.orthogonalDirectionSingularEdgeIds);
//    singularityGraph.WriteMatrixVTK("matrix.onTheSameFacesPatchSingularEdgeIds.vtk", singularityGraph.onTheSameFacesPatchSingularEdgeIds);
//    singularityGraph.WriteMatrixMat("matrix.directlyLinkedSingularEdge.mat", singularityGraph.directlyLinkedSingularEdgeIds);
//    singularityGraph.WriteMatrixMat("matrix.linkedByOneComponentEdgeSingularEdge.mat", singularityGraph.linkedByOneComponentEdgeSingularEdgeIds);
//    singularityGraph.WriteMatrixMat("matrix.orthogonalDirectionSingularEdgeIds.mat", singularityGraph.orthogonalDirectionSingularEdgeIds);
//    singularityGraph.WriteMatrixMat("matrix.onTheSameFacesPatchSingularEdgeIds.mat", singularityGraph.onTheSameFacesPatchSingularEdgeIds);
    return 0;
}


