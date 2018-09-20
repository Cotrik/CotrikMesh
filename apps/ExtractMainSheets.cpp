/*
 * ExtractMainSheets.cpp
 *
 *  Created on: March 25, 2018
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

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: ExtractMainSheets hex.vtk" << std::endl;
        return -1;
    }
    ArgumentManager argumentManager(argc, argv);
    std::string filename = argv[1];
    bool bfs = false;
    if (argumentManager.get("bfs") == "true") bfs = true;
    std::cout << "---------------------------------------" << std::endl;
    std::cout << "filename = " << filename << std::endl;
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
    BaseComplexSheet baseComplexSheets(baseComplex);
    baseComplexSheets.Extract();
    baseComplexSheets.ExtractSheetDecompositionsAll();
    baseComplexSheets.ExtractSheetConnectivities();
    std::cout << "#Sheets = " <<  baseComplexSheets.sheets_componentCellIds.size() << "\n";

    //    BaseComplexSheet baseComplexSheets(baseComplex);
    //    baseComplexSheets.Extract();
    //    baseComplexSheets.ExtractSheetDecompositions();
    //    baseComplexSheets.ExtractSheetConnectivities();
    //    baseComplexSheets.ExtractMainSheets();

//    BaseComplexSheet baseComplexSheets(baseComplex);
//    baseComplexSheets.Extract();
//    baseComplexSheets.ExtractSheetDecompositions();
//    baseComplexSheets.ExtractSheetConnectivities();
//    baseComplexSheets.ExtractMainSheets();

//    baseComplexSheets.WriteSheetsConnectivitiesMatrixVTK("SheetsConnectivities.vtk");
//    baseComplexSheets.WriteSheetsConnectivitiesMatrixMat("SheetsConnectivities.mat");
//    baseComplexSheets.ExtractMainSheetConnectivities(0);
//    baseComplexSheets.WriteSheetsConnectivitiesMatrixVTK("DominantSheetsConnectivities.vtk");
//    baseComplexSheets.WriteSheetsConnectivitiesMatrixMat("DominantSheetsConnectivities.mat");

    return 0;
}


