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
#include <fstream>
#include <algorithm>

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: ExtractSheetDecompositions hex.vtk verify=false" << std::endl;
        return -1;
    }
    ArgumentManager argumentManager(argc, argv);
    std::string filename = argv[1];
    bool verify = false;
    if (argumentManager.get("verify") == "true") verify = true;
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

    BaseComplexSheet baseComplexSheets(baseComplex);
    baseComplexSheets.Extract();
    // baseComplexSheets.ExtractSets();
    // baseComplexSheets.ExtractSheetDecompositions(bfs);
    baseComplexSheets.ExtractSheetDecompositionsAll();
    baseComplexSheets.VerifySheetDecompositions();
    if (verify)  baseComplexSheets.Verify();
    else {
    	std::vector<size_t> xx;
		for (int i = 0; i < baseComplexSheets.sheets_componentCellIds.size(); ++i) xx.push_back(i);
		baseComplexSheets.sheets_coverSheetIds.push_back(xx);
    }
    baseComplexSheets.ExtractSheetConnectivities();
    baseComplexSheets.WriteSheetsConnectivitiesMatrixVTK("SheetsConnectivities.vtk");
    baseComplexSheets.WriteSheetsConnectivitiesMatrixMat("SheetsConnectivities.mat");
    baseComplexSheets.ExtractMainSheetConnectivities();

    int sheet_id = 0;
    for (auto& sheetIds : baseComplexSheets.sheets_coverSheetIds){
        baseComplexSheets.ExtractMainSheetConnectivities(sheet_id);
        baseComplexSheets.ComputeComplexityUnbalancedMatrix(sheet_id++);
    }
    {

    	auto& vs = baseComplexSheets.socs;
    	  sort(vs.begin(), vs.end(),
    	       [&](const sheetIds_overlaps_complexity & a, const sheetIds_overlaps_complexity& b) {
    		  if (a.complexity < b.complexity) {
    			  return true;
    		  } else if (a.complexity > b.complexity) {
    			  return false;
    		  } else {
				   if (a.sheetIds.size() < b.sheetIds.size()) return true;
				   else if (a.sheetIds.size() == b.sheetIds.size()) {
					   std::vector<size_t> as, bs;
					   for (auto x : a.sheetIds) as.push_back(x);
					   for (auto x : b.sheetIds) bs.push_back(x);
					   for (int i = 0; i < as.size(); ++i) {
					       if (as.at(i) < bs.at(i)) return true;
					       else if (as.at(i) > bs.at(i)) return false;
					   }
					   return false;
				   }
				   return false;
			  }
    		  return false;
    	  });
    	  if (verify) {
    		  std::ofstream ofs("all.txt");
    		  ofs << "complexity\toverlaps\tsheets\n";
			  for (auto& v : vs) {
				  ofs << v.complexity << "/" << v.complexity*v.sheetIds.size() << "/" << v.sheetIds.size() << "\t";
				  ofs << v.overlaps << "\t";
				for (auto i : v.sheetIds) ofs << i << " ";
				ofs << "\n";
			  }
    	  } else {
    		  std::ofstream ofs("ours.txt");
    		  ofs << "complexity\toverlaps\tsheets\n";
			  for (auto& v : vs) {
				  ofs << v.complexity << "/" << v.complexity*v.sheetIds.size() << "/" << v.sheetIds.size() << "\t";
				  ofs << v.overlaps << "\t";
				for (auto i : v.sheetIds) ofs << i << " ";
				ofs << "\n";
			  }
    	  }
    	  baseComplexSheets.sheets_coverSheetIds.clear();
//    	  for (auto it = vs.rbegin(); it != vs.rend(); ++it) {
//    		  std::vector<size_t> x;
//    		  for (auto xx : it->sheetIds) x.push_back(xx);
//    		  baseComplexSheets.sheets_coverSheetIds.push_back(x);
//    	  }
          for (auto it = vs.rbegin(); it != vs.rend(); ++it) {
              std::vector<size_t> x;
              for (auto xx : it->sorted_sheetIds) x.push_back(xx);
              baseComplexSheets.sheets_coverSheetIds.push_back(x);
          }
        {
            std::ofstream ofs("sheet_decompositions.txt");
            for (auto& sheetIds : baseComplexSheets.sheets_coverSheetIds) {
                for (auto sheetId : sheetIds)
                    ofs << sheetId << " ";
                ofs << "\n";
            }
        }
//    	  vs.clear();
    	  sheet_id = 0;
		  for (auto& sheetIds : baseComplexSheets.sheets_coverSheetIds){
			  baseComplexSheets.ExtractMainSheetConnectivities(sheet_id);
			  baseComplexSheets.ComputeComplexityUnbalancedMatrix(sheet_id++, true);
		  }
    }

//    baseComplexSheets.ExtractMainSheetConnectivities(0);
//    baseComplexSheets.ComputeComplexityUnbalancedMatrix(0);
//    baseComplexSheets.WriteSheetsConnectivitiesMatrixVTK("DominantSheetsConnectivities.vtk");
//    baseComplexSheets.WriteSheetsConnectivitiesMatrixMat("DominantSheetsConnectivities.mat");
//    baseComplexSheets.ComputeComplexity();
//    baseComplexSheets.WriteAllDominantSheetsConnectivitiesMatrixMat("DominantSheetsConnectivities");
    return 0;
}


