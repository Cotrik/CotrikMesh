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

void OutputSheetDecompositonsDual(const char* mesh_filename, const char* dual_filename, const int scalar);
void OutputSheetDecompositonsDual_New(const char* mesh_filename, const char* dual_filename, const int scalar);

int main(int argc, char* argv[]){
    if (argc < 2) {
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
//    baseComplex.WriteAllSingularEdge_NeighborSeparatedComponentFacePatches_VTK("SingularFaces");
    baseComplex.WriteBaseComplexComponentsWithoutSingularitiesVTK("BaseComplexComponentsWithoutSingularities.vtk");
    baseComplex.WriteBaseComplexComponentsWithSingularitiesVTK("BaseComplexComponentsWithSingularities.vtk");
//    baseComplex.WriteBaseComplexAllComponentsVTK("ComponentCells");
//    baseComplex.WriteBaseComplexAllComponentsEdgesAndFacesVTK("ComponentEdgesAndFaces");
    baseComplex.WriteBaseComplexComponentsVTK("BaseComplexComponents.vtk");
    baseComplex.WriteBaseComplex_ColorFacesVTK("BaseComplexColorFaces.vtk");
    baseComplex.WriteBaseComplex_ColorEdgesVTK("BaseComplexColorEdges.vtk");
    baseComplex.WriteBaseComplex_ColorVerticesVTK("BaseComplexColorVertices.vtk");

    BaseComplexSheet baseComplexSheets(baseComplex);
    baseComplexSheets.Extract();
    baseComplexSheets.ExtractSets();
    baseComplexSheets.WriteAllSheetsCellsVTK("SheetCells");
//    baseComplexSheets.WriteAllSheetsFacesVTK("SheetFaces");
//    baseComplexSheets.WriteAllSheetsEdgesVTK("SheetEdges");
    baseComplexSheets.WriteAllSheetsFacesAndEdgesVTK("TestSheetFacesAndEdges");
    baseComplexSheets.WriteAllSheetsCellsDualVTK("SheetDual");
    baseComplexSheets.ExtractSheetDecompositions();

    for (auto sheet_id = 0; sheet_id < baseComplexSheets.GetNumOfSheets(); ++sheet_id) {
        std::string mesh_filename = "SheetDual" + std::to_string(sheet_id) + ".vtk";
        std::string dual_filename = "QuadDual" + std::to_string(sheet_id) + ".vtk";
        OutputSheetDecompositonsDual_New(mesh_filename.c_str(), dual_filename.c_str(), sheet_id);
    }

//    baseComplexSheets.ExtractSheetConnectivities();
//    baseComplexSheets.WriteSheetsConnectivitiesMatrixVTK("SheetsConnectivities.vtk");
//    baseComplexSheets.WriteSheetsConnectivitiesMatrixMat("SheetsConnectivities.mat");
//    baseComplexSheets.ExtractMainSheetConnectivities();
//    baseComplexSheets.WriteSheetsConnectivitiesMatrixVTK("DominantSheetsConnectivities.vtk");
//    baseComplexSheets.WriteSheetsConnectivitiesMatrixMat("DominantSheetsConnectivities.mat");
//    baseComplexSheets.ComputeComplexity();

    BaseComplexChord baseComplexChords(baseComplex);
    baseComplexChords.Extract();
//    baseComplexChords.WriteAllChordsCellsVTK("ChordCells");
//    baseComplexChords.WriteAllChordsFacesVTK("ChordFaces");
//    baseComplexChords.WriteAllChordsEdgesVTK("ChordEdges");
//    baseComplexChords.WriteAllChordsFacesAndEdgesVTK("ChordFacesAndEdges");
//    baseComplexChords.WriteAllChordsCurvesVTK("ChordCurves");
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


void OutputSheetDecompositonsDual(const char* mesh_filename, const char* dual_filename, const int scalar) {
    MeshFileReader reader(mesh_filename);
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

    BaseComplexQuad baseComplex(mesh);
    baseComplex.Build();
    BaseComplexSheetQuad baseComplexSheets(baseComplex);
    baseComplexSheets.Extract();
    baseComplexSheets.ExtractSheetDecompositions();
    baseComplexSheets.WriteSheetDecompositionsDuaVTK(dual_filename, scalar);

    baseComplexSheets.ExtractSheetConnectivities();
//    baseComplexSheets.WriteSheetsConnectivitiesMatrixVTK("SheetsConnectivities.vtk");
//    baseComplexSheets.WriteSheetsConnectivitiesMatrixMat("SheetsConnectivities.mat");
    baseComplexSheets.ExtractMainSheetConnectivities();
//    baseComplexSheets.WriteSheetsConnectivitiesMatrixVTK("DominantSheetsConnectivities.vtk");
//    baseComplexSheets.WriteSheetsConnectivitiesMatrixMat("DominantSheetsConnectivities.mat");
    baseComplexSheets.ComputeComplexity();
}

void OutputSheetDecompositonsDual_New(const char* mesh_filename, const char* dual_filename, const int scalar) {
    MeshFileReader reader(mesh_filename);
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

    BaseComplexQuad baseComplex(mesh);
    baseComplex.Build();
    BaseComplexSheetQuad baseComplexSheets(baseComplex);
    baseComplexSheets.Extract();
    baseComplexSheets.ExtractSheetDecompositionsAll();
    baseComplexSheets.ExtractSheetConnectivities();

    auto sheetType = baseComplexSheets.GetSheetType();
    if (sheetType == NON_SIMPLE) {
        std::cout << "sheet " << scalar << " is NON_SIMPLE\n";
        baseComplexSheets.WriteSelfIntersectingEdges((std::string("SelfIntersectingEdges") + std::to_string(scalar) + ".vtk").c_str());
    }

    int sheet_id = 0;
    float min_complexity = INT_MAX;
    int min_complexity_id = 0;
    std::vector<pair<int, float>> representativeSheetId_complexities(baseComplexSheets.Get_sheets_coverSheetIds().size());
    for (auto& sheetIds : baseComplexSheets.Get_sheets_coverSheetIds()) {
        baseComplexSheets.ExtractMainSheetConnectivities(sheet_id);
        float complexity = baseComplexSheets.ComputeComplexityUnbalancedMatrix(sheet_id);
        representativeSheetId_complexities.push_back(std::make_pair(sheet_id, complexity));
        if (complexity < min_complexity) {
            min_complexity = complexity;
            min_complexity_id = sheet_id;
        }
        ++sheet_id;
    }
    std::cout << "min_main_chords_complexity = " << min_complexity << "\n";
    baseComplexSheets.WriteSheetDecompositionsDuaVTK(dual_filename, scalar, min_complexity_id);

    {
        std::ofstream ofs("chord_decompositions.txt");
        for (auto& sheetIds : baseComplexSheets.Get_sheets_coverSheetIds()) {
            for (auto sheetId : sheetIds) ofs << sheetId << " ";
            ofs << "\n";
        }
    }
    baseComplexSheets.ExtractMainSheetConnectivities(0);
    baseComplexSheets.ComputeComplexityUnbalancedMatrix(0);
}
