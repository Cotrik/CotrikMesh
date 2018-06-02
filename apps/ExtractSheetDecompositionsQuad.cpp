/*
 * ExtractSheetDecompositionsQuad.cpp
 *
 *  Created on: Jun 11, 2018
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
        std::cout << "Usage: ExtractSheetDecompositionsQuad quad.vtk bfs=false" << std::endl;
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

    BaseComplexQuad baseComplex(mesh);
    baseComplex.Build();

    BaseComplexSheetQuad baseComplexSheets(baseComplex);
    baseComplexSheets.Extract();
    //baseComplexSheets.ExtractSets();
    baseComplexSheets.ExtractSheetDecompositions(bfs);
    baseComplexSheets.ExtractSheetConnectivities();
    int sheet_id = 0;
    for (auto& sheetIds : baseComplexSheets.Get_sheets_coverSheetIds()){
        baseComplexSheets.ExtractMainSheetConnectivities(sheet_id);
        baseComplexSheets.ComputeComplexityDrChen(sheet_id++);
    }
    {
        std::ofstream ofs("sheet_decompositions.txt");
        for (auto& sheetIds : baseComplexSheets.Get_sheets_coverSheetIds()) {
            for (auto sheetId : sheetIds) ofs << sheetId << " ";
            ofs << "\n";
        }
    }
    baseComplexSheets.ExtractMainSheetConnectivities(0);
    baseComplexSheets.ComputeComplexityDrChen(0);

    return 0;
}


