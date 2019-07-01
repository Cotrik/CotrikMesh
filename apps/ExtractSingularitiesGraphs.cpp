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
        std::cout << "Usage: ExtractSingularitiesGraphs hex.vtk" << std::endl;
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
    auto sg = baseComplex.SingularityI.ExtractSubConnectedGraphs();

    std::cout << "#Singularites = " << baseComplex.SingularityI.E.size() << "\n";
    std::cout << "#Singularity Graphs = " << sg.size() << "\n";
    int x = 0;
    for (auto& graph : sg) {
        double total_index = 0;
        for (auto& svid : graph.first) {
            auto& sv = baseComplex.SingularityI.V.at(svid);
            const auto& v = mesh.V.at(sv.id_mesh);
            const auto valence = v.N_Cids.size();
            if (v.isBoundary) total_index += 0.5 - valence * 0.125;
            else total_index += 1.0 - valence * 0.125;
            // std::cout << total_index << std::endl;
        }
        for (auto& seid : graph.second) {
            auto& se = baseComplex.SingularityI.E.at(seid);
            const auto& e = mesh.E.at(se.es_link.front());
            const auto valence = e.N_Cids.size();
            if (e.isBoundary) total_index -= 0.5 - valence * 0.25;
            else total_index -= 1.0 - valence * 0.25;
            // std::cout << total_index << std::endl;
        }
        std::cout << "TotalIndex = " << total_index << "\n";
        baseComplex.WriteSingularities_VTK((std::string("singularities") + std::to_string(x) + ".vtk").c_str(), graph);

        std::set<size_t> set_cids;
        std::map<size_t, size_t> cid_count;
        for (auto& seid : graph.second) {
            auto& se = baseComplex.SingularityI.E.at(seid);
            for (auto eid : se.es_link) {
                const auto& e = mesh.E.at(eid);
                set_cids.insert(e.N_Cids.begin(), e.N_Cids.end());
                for (auto cid : e.N_Cids)
                    ++cid_count[cid];
            }
        }
        std::ofstream ofs(std::string("cids") + std::to_string(x) + ".txt");
        for (auto& item : cid_count)
            ofs << item.first << " " << item.second << "\n";

        std::set<size_t> removable_seids;
        for (auto& seid : graph.second) {
            auto& se = baseComplex.SingularityI.E.at(seid);
            bool removable = true;
            for (auto eid : se.es_link) {
                const auto& e = mesh.E.at(eid);
                for (auto cid : e.N_Cids)
                    if (cid_count[cid] <= 1) {
                        removable = false;
                        break;
                    }
                if (!removable)
                    break;
            }
            if (removable) {
                for (auto eid : se.es_link) {
                    const auto& e = mesh.E.at(eid);
                    for (auto cid : e.N_Cids)
                        --cid_count[cid];
                }
                removable_seids.insert(seid);
            }
        }
        for (auto& seid : removable_seids)
        graph.second.erase(seid);
        baseComplex.WriteSingularities_VTK((std::string("singularitiesReduced") + std::to_string(x) + ".vtk").c_str(), graph);
        //std::vector<size_t> cids(set_cids.begin(), set_cids.end());
        std::vector<Cell> cells;
        for (auto cid : set_cids)
            cells.push_back(mesh.C.at(cid));
        MeshFileWriter writer(mesh.V, cells, (std::string("singularitiesHex") + std::to_string(x++) + ".vtk").c_str());
        writer.WriteFile();
    }
//    baseComplex.WriteBaseComplexHexVTK("BaseComplexHex.vtk");
//    baseComplex.WriteBaseComplex_VTK("BaseComplex.vtk");
//    baseComplex.WriteSingularV_VTK("singularV.vtk");
//    baseComplex.WriteSingularE_VTK("singularE.vtk");
//    baseComplex.WriteSingularities_VTK("singularities.vtk");
//    baseComplex.WriteBaseComplexComponentsVTK("BaseComplexComponents.vtk");
//    baseComplex.WriteBaseComplex_ColorFacesVTK("BaseComplexColorFaces.vtk");
//    baseComplex.WriteBaseComplex_ColorEdgesVTK("BaseComplexColorEdges.vtk");
//    baseComplex.WriteBaseComplex_ColorVerticesVTK("BaseComplexColorVertices.vtk");
//    baseComplex.WriteBaseComplexAllComponentsVTK("ComponentCells");

//    BaseComplexSheet baseComplexSheet(baseComplex);
//    baseComplexSheet.Extract();
//
//    std::cout << "#Sheets = " << baseComplexSheet.sheets_componentCellIds.size() << "\n";
    return 0;
}
