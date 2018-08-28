/*
 * BaseComplexSheet.cpp
 *
 *  Created on: Oct 5, 2017
 *      Author: cotrik
 */

#include "BaseComplexSheet.h"
#include "MeshFileWriter.h"
#include "Dual.h"
#include "RefinedDual.h"
#include <algorithm>
#include <sstream>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <stack>
#include <queue>
#include <iterator>
#include <Eigen/Dense>

const std::string red_color = "\033[1;31m";
const std::string end_color = "\033[0m";

Mesh GetRefineMesh(const Mesh& hex_mesh, int clockwise = 0)
{
    const Mesh& new_mesh = hex_mesh;
    ////////////////////////////////////////////////////////////////////////////
    // add vertices
    std::vector<Vertex> new_vertex(new_mesh.V.size() + new_mesh.E.size() + new_mesh.F.size() + new_mesh.C.size());
    for (size_t i = 0; i < new_mesh.V.size(); i++)
        new_vertex.at(i) = new_mesh.V.at(i);
    size_t offset = new_mesh.V.size();
    for (size_t i = 0; i < new_mesh.E.size(); i++) {
        const Edge& e = new_mesh.E.at(i);
        const Vertex& v0 = new_mesh.V.at(e.Vids[0]);
        const Vertex& v1 = new_mesh.V.at(e.Vids[1]);
        new_vertex.at(offset + i) = 0.5f * (v0.xyz() + v1.xyz());
    }
    offset = new_mesh.V.size() + new_mesh.E.size();
    for (size_t i = 0; i < new_mesh.F.size(); i++) {
        const Face& f = new_mesh.F.at(i);
        const Vertex& v0 = new_mesh.V.at(f.Vids[0]);
        const Vertex& v1 = new_mesh.V.at(f.Vids[2]);
        new_vertex.at(offset + i) = 0.5f * (v0.xyz() + v1.xyz());
    }
    offset = new_mesh.V.size() + new_mesh.E.size() + new_mesh.F.size();
    for (size_t i = 0; i < new_mesh.C.size(); i++) {
        const Cell& c = new_mesh.C.at(i);
        const Vertex& v0 = new_mesh.V.at(c.Vids[0]);
        const Vertex& v1 = new_mesh.V.at(c.Vids[6]);
        new_vertex.at(offset + i) = 0.5f * (v0.xyz() + v1.xyz());
    }
    //new_mesh.V = new_vertex;
    /////////////////////////////////////////////////////////////////
    // add cells
    const unsigned int HexEdge[12][2] = { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 }, { 4, 5 }, { 5, 6 }, { 6, 7 }, { 7, 4 }, { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 }, };

    const int HexRefine[8][8] =
    {
        11, 20, 10, 3, 22, 26, 25, 19,
        20, 9, 2, 10, 26, 23, 18, 25,
        22, 26, 25, 19, 15, 21, 14, 7,
        26, 23, 18, 25, 21, 13, 6, 14,
        0, 8, 20, 11, 16, 24, 26, 22,
        8, 1, 9, 20, 24, 17, 23, 26,
        16, 24, 26, 22, 4, 12, 21, 15,
        24, 17, 23, 26, 12, 5, 13, 21
    };

    Cell cell(8);
    std::vector<Cell> new_cells(8 * new_mesh.C.size(), cell);
    int count = 0;
    for (size_t i = 0; i < new_mesh.C.size(); i++) {
        unsigned long v_index[27];
        const Cell& c = new_mesh.C.at(i);
        for (auto j = 0; j < 8; j++)
            v_index[j] = c.Vids.at(j);
        if (clockwise != 0) {
            std::swap(v_index[1], v_index[3]);
            std::swap(v_index[5], v_index[7]);
        }
        for (unsigned long j = 0; j < 12; j++) {
            const Edge e({c.Vids.at(HexEdge[j][0]), c.Vids.at(HexEdge[j][1])});
            auto it = std::find(new_mesh.E.begin(), new_mesh.E.end(), e);
            if (it == new_mesh.E.end()) std::cout << "Edge search Error !" << std::endl;
            const unsigned long e_index = std::distance(new_mesh.E.begin(), it);
            v_index[8 + j] = new_mesh.V.size() + e_index;
        }
        auto face_qual = [](const Face& a, const Face& b) {
            std::unordered_set<size_t> s;
            s.insert(a.Vids.begin(), a.Vids.end());
            s.insert(b.Vids.begin(), b.Vids.end());
            return s.size() == 4;
        };
        for (unsigned long j = 0; j < 6; j++) {
            Face f({c.Vids.at(HexFaces[j][0]), c.Vids.at(HexFaces[j][1]), c.Vids.at(HexFaces[j][2]), c.Vids.at(HexFaces[j][3])});
            for (auto f_index = 0; f_index < new_mesh.F.size(); ++f_index) {
                auto& a = new_mesh.F.at(f_index);
                std::unordered_set<size_t> s;
                s.insert(a.Vids.begin(), a.Vids.end());
                s.insert(f.Vids.begin(), f.Vids.end());
                if (s.size() == 4) {
                    v_index[20 + j] = new_mesh.V.size() + new_mesh.E.size() + f_index;
                    break;
                }
            }
//            auto it = std::find_if(new_mesh.F.begin(), new_mesh.F.end(), f);
//            if (it == new_mesh.F.end()) std::cout << "Face search Error !" << std::endl;
//            const unsigned long f_index = std::distance(new_mesh.F.begin(), it);
//            v_index[20 + j] = new_mesh.V.size() + new_mesh.E.size() + f_index;
        }
        v_index[26] = new_mesh.V.size() + new_mesh.E.size() + new_mesh.F.size() + i;
        //Cell new_cell(8, 0);
        for (int k = 0; k < 8; k++, count++)
            for (int j = 0; j < 8; j++)
                new_cells[count].Vids[j] = v_index[HexRefine[k][j]];
    }
    Mesh mesh(new_vertex, new_cells, HEXAHEDRA);
    return mesh;
}

std::set<std::set<size_t>> shortestCombination(const std::set<size_t>& filter, const std::vector<std::set<size_t>>& listOfSets) {
    size_t size = listOfSets.size();
    if (size > 20) {
        std::cerr << "Too many combinations\n";
        return {};
    }
    int combinations = 1 << size;
    std::vector<std::set<std::set<size_t>>> possibleSolutions;
    for (int l = 0; l < combinations; l++) {
        std::set<std::set<size_t>> combination;
        for (int j = 0; j < size; j++)
            if (((l >> j) & 1) != 0) combination.insert(listOfSets.at(j));
        possibleSolutions.push_back(combination);
    }

    for (auto& possibleSolution : possibleSolutions) {
        std::set<size_t> s;
        for (auto& ss : possibleSolution)
            for (auto sss : ss)
                s.insert(sss);
        if (filter == s) return possibleSolution;
    }
    return {};
}
void FindSetCombination(std::vector<std::set<size_t>>& input, std::set<size_t>& target, std::vector<std::set<size_t>>& output) {
    std::set<int> full;
    for (auto it : input)
        full.insert(it.begin(), it.end());

    if (!includes(full.begin(), full.end(), target.begin(), target.end())) return;

    for (size_t i = input.size() - 1; i > 0; --i) {
        std::vector<bool> vec(input.size(), false);
        std::fill(vec.begin() + i, vec.end(), true);
        std::set<int> comb;

        do {
            for (size_t j = 0; j < vec.size(); ++j)
                if (vec[j]) comb.insert(input[j].begin(), input[j].end());
            if (includes(comb.begin(), comb.end(), target.begin(), target.end())) {
                for (size_t j = 0; j < vec.size(); ++j)
                    if (vec[j]) output.push_back(input[j]);
                return;
            }
            comb.clear();

        } while (next_permutation(vec.begin(), vec.end()));
    }
}

std::vector<std::vector<std::set<size_t>>> FindSetCombination(std::vector<std::set<size_t>>& input, std::set<size_t>& target) {

    std::set<int> full;
    for (auto it : input)
        full.insert(it.begin(), it.end());

    if (!includes(full.begin(), full.end(), target.begin(), target.end())) return {};
    std::vector<std::vector<std::set<size_t>>> res;
    std::vector<std::set<size_t>> output;
    for (size_t i = input.size() - 1; i > 0; --i) {
        std::vector<bool> vec(input.size(), false);
        std::fill(vec.begin() + i, vec.end(), true);
        std::set<int> comb;

        do {
            for (size_t j = 0; j < vec.size(); ++j)
                if (vec[j]) comb.insert(input[j].begin(), input[j].end());
            if (includes(comb.begin(), comb.end(), target.begin(), target.end())) {
                for (size_t j = 0; j < vec.size(); ++j)
                    if (vec[j]) output.push_back(input[j]);
                res.push_back(output);
            }
            comb.clear();

        } while (next_permutation(vec.begin(), vec.end()));
    }
}

BaseComplexSheet::BaseComplexSheet(BaseComplex& baseComplex)
: baseComplex(baseComplex)
{
    // TODO Auto-generated constructor stub

}

BaseComplexSheet::~BaseComplexSheet()
{
    // TODO Auto-generated destructor stub
}
/*
         3____________________2
         /|                 /|
        / |                / |
       /  |               /  |
   0  /___|______________/ 1 |
      |   |              |   |
      |   |              |   |
      |   |              |   |
      |   |______________|___|
      |   / 7            |  / 6
      |  /               | /
      | /                |/
      |/_________________/
    4                   5
*/
const size_t parallel_edge_ids[3][4][2] = {
        {{0, 4}, {1, 5}, {2, 6}, {3, 7}},
        {{0, 1}, {3, 2}, {7, 6}, {4, 5}},
        {{0, 3}, {1, 2}, {5, 6}, {1, 7}}
};
void BaseComplexSheet::Extract()
{
    std::vector<bool> visited(baseComplex.componentE.size(), false);
    for (auto& componentEdge : baseComplex.componentE) {
        if (visited.at(componentEdge.id)) continue;
        std::vector<size_t> sheetComponentEdgeIds;
        std::vector<size_t> sheetComponentFaceIds;
        std::vector<size_t> sheetComponentCellIds;
        GetParallelComponents(componentEdge, sheetComponentEdgeIds, sheetComponentFaceIds, sheetComponentCellIds);
        sheets_componentEdgeIds.push_back(sheetComponentEdgeIds);
        sheets_componentFaceIds.push_back(sheetComponentFaceIds);
        sheets_componentCellIds.push_back(sheetComponentCellIds);

        for (auto id : sheetComponentEdgeIds) visited.at(id) = true;
    }
}

void BaseComplexSheet::ExtractSets()
{
//    std::vector<std::set<size_t>> input;
//    std::set<size_t> target;
//    std::vector<std::set<size_t>> output;
//    for (auto& i : sheets_componentCellIds)
//        input.push_back(std::set<size_t>(i.begin(), i.end()));
//    for (auto i = 0; i < baseComplex.componentC.size(); target.insert(i++));
//    FindSetCombination(input, target, output);
//    std::cout << "*** sheets set size " << output.size() << " ***" << std::endl;
//    std::set<size_t> ids;
//    for (auto i : output)
//        for (auto j = 0; j < input.size(); ++j)
//            if (i == input[j]) {
//                ids.insert(j);
//                break;
//            }
//    std::cout << "***************************" << std::endl;
//    for (auto i : ids)
//        std::cout << i << " ";
//    std::cout << std::endl;
//    std::cout << "***************************" << std::endl;
//
//    Mesh mesh = GetRefineMesh(baseComplex.mesh);
//    MeshFileWriter writer(mesh, "refined.vtk");
//    writer.WriteFile();
}

void BaseComplexSheet::ExtractSheetDecompositionsAll() {
    ExtractSheetDecompositions(false);
    std::cout << "\n---- sheetDecomposition bfs = false! ----" << "\n";
    VerifySheetDecompositions();
    auto representativeSheetSets = Get_sheets_coverSheetIds();
    {
        BaseComplexSheet baseComplexSheets1(baseComplex);
        baseComplexSheets1.Extract();
        baseComplexSheets1.ExtractSets();
        baseComplexSheets1.ExtractSheetDecompositions(true);
        std::cout << "\n---- sheetDecomposition bfs = true! ----" << "\n";
        baseComplexSheets1.VerifySheetDecompositions();
        std::copy(baseComplexSheets1.Get_sheets_coverSheetIds().begin(), baseComplexSheets1.Get_sheets_coverSheetIds().end(), back_inserter(representativeSheetSets));
    }
    std::vector<bool> componentFaceBoundary(baseComplex.componentF.size(), false);
    for (auto& bF : baseComplex.componentF) {
        componentFaceBoundary[bF.id] = bF.isBoundary;
        bF.isBoundary = false;
    }
    {
        BaseComplexSheet baseComplexSheets1(baseComplex);
        baseComplexSheets1.Extract();
        baseComplexSheets1.ExtractSets();
        baseComplexSheets1.ExtractSheetDecompositions(false);
        std::cout << "\n---- bF.isBoundary = false; sheetDecomposition bfs = false! ----" << "\n";
        baseComplexSheets1.VerifySheetDecompositions();
        std::copy(baseComplexSheets1.Get_sheets_coverSheetIds().begin(), baseComplexSheets1.Get_sheets_coverSheetIds().end(), back_inserter(representativeSheetSets));
    }
    {
        BaseComplexSheet baseComplexSheets1(baseComplex);
        baseComplexSheets1.Extract();
        baseComplexSheets1.ExtractSets();
        baseComplexSheets1.ExtractSheetDecompositions(true);
        std::cout << "\n---- bF.isBoundary = false; sheetDecomposition bfs = true! ----" << "\n";
        baseComplexSheets1.VerifySheetDecompositions();
        std::copy(baseComplexSheets1.Get_sheets_coverSheetIds().begin(), baseComplexSheets1.Get_sheets_coverSheetIds().end(), back_inserter(representativeSheetSets));
    }
    for (auto& bF : baseComplex.componentF)
        bF.isBoundary = componentFaceBoundary[bF.id];
    sheets_coverSheetIds = representativeSheetSets;
    RemoveSheetSetsRedundancy();
    std::sort(sheets_coverSheetIds.begin(), sheets_coverSheetIds.end(), [&](const std::vector<size_t>& a, const std::vector<size_t>& b) {return a.size() < b.size();});
    std::cout << "****** SheetDecompositions after RemoveSheetSetsRedundancy******\n";
    for (auto& sheetIds : sheets_coverSheetIds) {
        for (auto sheetId : sheetIds) std::cout << sheetId << " ";
        std::cout << "\n";
//        break;
    }
}
void BaseComplexSheet::ExtractSheetDecompositions(const bool bfs) {
//    sheets_hashComponentCellIds.resize(sheets_componentCellIds.size());
//    for (size_t sheetId = 0; sheetId < sheets_componentCellIds.size(); ++sheetId) {
//        auto& cellIds = sheets_componentCellIds.at(sheetId);
//        sheets_hashComponentCellIds[sheetId].insert(cellIds.begin(), cellIds.end());
//    }
    cellComponent_sheetIds.resize(baseComplex.componentC.size());
    for (size_t sheetId = 0; sheetId < sheets_componentCellIds.size(); ++sheetId) {
        for (auto component : sheets_componentCellIds[sheetId])
            if (cellComponent_sheetIds[component].find(sheetId) == cellComponent_sheetIds[component].end())
                cellComponent_sheetIds[component].insert(sheetId);
    }
    for (size_t sheetId = 0; sheetId < sheets_componentCellIds.size(); ++sheetId) {
        auto sheetBoundaryFaceComponentIds = GetSheetBoundaryFaceComponentIds(sheetId);
        for (auto faceComponentId : sheetBoundaryFaceComponentIds)
            faceComponent_neighborSheetIds[faceComponentId].insert(sheetId);
        sheets_BoundaryFaceComponentIds.push_back(sheetBoundaryFaceComponentIds);
    }
    for (size_t sheetId = 0; sheetId < sheets_componentCellIds.size(); ++sheetId)
        sheets_coverSheetIds.push_back(!bfs ? GetCoverSheetIds(sheetId) : GetCoverSheetIdsBFS(sheetId));

    std::sort(sheets_coverSheetIds.begin(), sheets_coverSheetIds.end(), [&](const std::vector<size_t>& a, const std::vector<size_t>& b) {return a.size() < b.size();});
//    std::cout << "****** SheetDecompositions ******\n";
//    for (auto& sheetIds : sheets_coverSheetIds) {
//        for (auto sheetId : sheetIds) std::cout << sheetId << " ";
//        std::cout << "\n";
////        break;
//    }
    RemoveSheetSetsRedundancy();
    std::sort(sheets_coverSheetIds.begin(), sheets_coverSheetIds.end(), [&](const std::vector<size_t>& a, const std::vector<size_t>& b) {return a.size() < b.size();});
//    std::cout << "****** SheetDecompositions after RemoveSheetSetsRedundancy******\n";
//    for (auto& sheetIds : sheets_coverSheetIds) {
//        for (auto sheetId : sheetIds) std::cout << sheetId << " ";
//        std::cout << "\n";
////        break;
//    }

    std::ofstream ofs("sheet_decompositions.txt");
    for (auto& sheetIds : sheets_coverSheetIds) {
        for (auto sheetId : sheetIds) ofs << sheetId << " ";
        ofs << "\n";
    }
    WriteSheetNeighborsSheetIdsJS("sheet_neighborsheetids.js");
}

static void combine(std::vector<size_t>& com, std::vector<std::vector<size_t> > &res, int n, int k, int start) {
    if (k == com.size()) {
        res.push_back(com);
        return;
    }
    for (int i = start; i < n; ++i) {
        com.push_back(i);
        combine(com, res, n, k, i + 1);
        com.pop_back();
    }
}

static std::vector<std::vector<size_t>> combine(int n, int k) {
    std::vector<std::vector<size_t>> res;
    std::vector<size_t> com;
    combine(com, res, n, k, 0);
    return res;
}

void BaseComplexSheet::ExtractSheetConnectivities() {
    size_t n = sheets_componentCellIds.size();
    sheets_connectivities.resize(n, std::vector<size_t>(n, 0));

    if (n < 1) return;
    std::vector<std::unordered_set<size_t>> sheets_hashcomponentEdgeIds;  // each sheet consists a list of base-complex componentEdge ids;
    std::vector<std::unordered_set<size_t>> sheets_hashcomponentFaceIds;  // each sheet consists a list of base-complex componentFace ids;
    std::vector<std::unordered_set<size_t>> sheets_hashcomponentCellIds;  // each sheet consists a list of base-complex componentCell ids;

    for (auto & sheet_componentEdgeIds : sheets_componentEdgeIds)
        sheets_hashcomponentEdgeIds.push_back(std::unordered_set<size_t> (sheet_componentEdgeIds.begin(), sheet_componentEdgeIds.end()));
    for (auto & sheet_componentFaceIds : sheets_componentFaceIds)
        sheets_hashcomponentFaceIds.push_back(std::unordered_set<size_t> (sheet_componentFaceIds.begin(), sheet_componentFaceIds.end()));
    for (auto & sheet_componentCellIds : sheets_componentCellIds)
        sheets_hashcomponentCellIds.push_back(std::unordered_set<size_t> (sheet_componentCellIds.begin(), sheet_componentCellIds.end()));

//    // Get neigbor connectivities
//    for (auto& sheetBoundaryFaceComponentIds : sheets_BoundaryFaceComponentIds) {
//        std::unordered_set<size_t> neighborSheetIds;
//        for (auto sheetBoundaryFaceComponentId : sheetBoundaryFaceComponentIds) {
//            if (baseComplex.componentF.at(sheetBoundaryFaceComponentId).isBoundary) continue;
//            const auto & sheetIds = faceComponent_neighborSheetIds[sheetBoundaryFaceComponentId];
//            if (sheetIds.size() <= 1) continue;
//            auto& neighborComponentCellIds = baseComplex.componentF.at(sheetBoundaryFaceComponentId).N_Cids;
//            auto id1 = neighborComponentCellIds[0];
//            auto id2 = neighborComponentCellIds[1];
//            auto combinations = combine(sheetIds.size(), 2);
//            for (auto combination : combinations) {
//                auto sheetId1 = combination[0];
//                auto sheetId2 = combination[1];
//                if ((sheets_hashcomponentCellIds[sheetId1].find(id1) != sheets_hashcomponentCellIds[sheetId1].end() &&
//                        sheets_hashcomponentCellIds[sheetId2].find(id2) != sheets_hashcomponentCellIds[sheetId2].end()) ||
//                        (sheets_hashcomponentCellIds[sheetId1].find(id2) != sheets_hashcomponentCellIds[sheetId1].end() &&
//                        sheets_hashcomponentCellIds[sheetId2].find(id1) != sheets_hashcomponentCellIds[sheetId2].end())) {
//                    sheets_connectivities[sheetId1][sheetId2] = SheetsConnectivity_NEIGHBOR;
//                    sheets_connectivities[sheetId2][sheetId1] = SheetsConnectivity_NEIGHBOR;
//                }
//            }
//        }
//    }
    // Get intersections
    for (const auto& sheetIds : cellComponent_sheetIds)
        if (sheetIds.size() > 1) {
            auto combinations = combine(sheetIds.size(), 2);
            std::vector<size_t> sheet_ids(sheetIds.begin(), sheetIds.end());
            for (auto combination : combinations) {
                auto sheetid1 = sheet_ids[combination[0]];
                auto sheetid2 = sheet_ids[combination[1]];
                sheets_connectivities[sheetid1][sheetid2] = SheetsConnectivity_INTERSECT;
                sheets_connectivities[sheetid2][sheetid1] = SheetsConnectivity_INTERSECT;
            }
        }

    // Get neigbor connectivities
    size_t sheetid = 0;
    for (auto& sheetBoundaryFaceComponentIds : sheets_BoundaryFaceComponentIds) {
        std::unordered_set<size_t> neighborSheetIds;
        for (auto sheetBoundaryFaceComponentId : sheetBoundaryFaceComponentIds) {
            const auto & sheetIds = faceComponent_neighborSheetIds[sheetBoundaryFaceComponentId];
            if (sheetIds.size() <= 1) continue;
            if (baseComplex.componentF.at(sheetBoundaryFaceComponentId).isBoundary) continue;

            for (auto neighborSheetId : faceComponent_neighborSheetIds[sheetBoundaryFaceComponentId])
                if (neighborSheetId != sheetid) neighborSheetIds.insert(neighborSheetId);
        }
        for (auto neighborSheetId : neighborSheetIds) {
            if (sheets_connectivities[sheetid][neighborSheetId] == SheetsConnectivity_UNKNOWN) {
                sheets_connectivities[sheetid][neighborSheetId] = SheetsConnectivity_NEIGHBOR;
                sheets_connectivities[neighborSheetId][sheetid] = SheetsConnectivity_NEIGHBOR;
            }
            else if (sheets_connectivities[sheetid][neighborSheetId] == SheetsConnectivity_INTERSECT) {
                sheets_connectivities[sheetid][neighborSheetId] = SheetsConnectivity_NEIGHBOR_AND_INTERSECT;
                sheets_connectivities[neighborSheetId][sheetid] = SheetsConnectivity_NEIGHBOR_AND_INTERSECT;
            }
        }
        ++sheetid;
    }



    sheets_connectivities_float.resize(n, std::vector<float>(n, 0.0f));
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            if (sheets_connectivities[i][j] == SheetsConnectivity_NEIGHBOR) sheets_connectivities_float[i][j] = -1.0f;
            else if (sheets_connectivities[i][j] == SheetsConnectivity_INTERSECT) sheets_connectivities_float[i][j] = 1.0f;
            else if (sheets_connectivities[i][j] == SheetsConnectivity_NEIGHBOR_AND_INTERSECT) sheets_connectivities_float[i][j] = 0.5f;

    all_sheets_connectivities = sheets_connectivities;
}

//void BaseComplexSheet::ExtractSheetConnectivities() {
//    size_t n = sheets_componentCellIds.size();
//    sheets_connectivities.resize(n, std::vector<size_t>(n, 0));
//
//    // Get neigbor connectivities
//    size_t sheetid = 0;
//    for (auto& sheetBoundaryFaceComponentIds : sheets_BoundaryFaceComponentIds) {
//        std::unordered_set<size_t> neighborSheetIds;
//        for (auto sheetBoundaryFaceComponentId : sheetBoundaryFaceComponentIds) {
//            const auto & sheetIds = faceComponent_neighborSheetIds[sheetBoundaryFaceComponentId];
//            if (sheetIds.size() <= 1) continue;
//            if (baseComplex.componentF.at(sheetBoundaryFaceComponentId).isBoundary) {
////                const auto & sheetIds = faceComponent_neighborSheetIds[sheetBoundaryFaceComponentId];
////                if (sheetIds.size() > 1) {
////                    auto combinations = combine(sheetIds.size(), 2);
////                    for (auto combination : combinations)
////                        if (sheets_connectivities[combination[0]][combination[1]] == 0) {
////                            sheets_connectivities[combination[0]][combination[1]] = SheetsConnectivity_INTERSECT;
////                            sheets_connectivities[combination[1]][combination[0]] = SheetsConnectivity_INTERSECT;
////                        }
////                }
//                continue;
//            } else
//            for (auto neighborSheetId : faceComponent_neighborSheetIds[sheetBoundaryFaceComponentId])
//                neighborSheetIds.insert(neighborSheetId);
//        }
//        for (auto neighborSheetId : neighborSheetIds) {
//            sheets_connectivities[sheetid][neighborSheetId] = SheetsConnectivity_NEIGHBOR;
//            sheets_connectivities[neighborSheetId][sheetid] = SheetsConnectivity_NEIGHBOR;
//        }
//        ++sheetid;
//    }
//
//    // Get intersections
//    for (const auto& sheetIds : cellComponent_sheetIds)
//        if (sheetIds.size() > 1) {
//            auto combinations = combine(sheetIds.size(), 2);
//            for (auto combination : combinations) {
//                if (sheets_connectivities[combination[0]][combination[1]] == SheetsConnectivity_UNKNOWN) {
//                    sheets_connectivities[combination[0]][combination[1]] = SheetsConnectivity_INTERSECT;
//                    sheets_connectivities[combination[1]][combination[0]] = SheetsConnectivity_INTERSECT;
//                } else if (sheets_connectivities[combination[0]][combination[1]] == SheetsConnectivity_NEIGHBOR) {
//                    sheets_connectivities[combination[0]][combination[1]] = SheetsConnectivity_NEIGHBOR_AND_INTERSECT;
//                    sheets_connectivities[combination[1]][combination[0]] = SheetsConnectivity_NEIGHBOR_AND_INTERSECT;
//                }
//            }
//        }
//
//    sheets_connectivities_float.resize(n, std::vector<float>(n, 0.0f));
//    for (size_t i = 0; i < n; ++i)
//        for (size_t j = 0; j < n; ++j)
//            if (sheets_connectivities[i][j] == SheetsConnectivity_NEIGHBOR) sheets_connectivities_float[i][j] = 0.33f;
//            else if (sheets_connectivities[i][j] == SheetsConnectivity_INTERSECT) sheets_connectivities_float[i][j] = 0.67f;
//            else if (sheets_connectivities[i][j] == SheetsConnectivity_NEIGHBOR_AND_INTERSECT) sheets_connectivities_float[i][j] = 1.0f;
//}

void BaseComplexSheet::ExtractMainSheetConnectivities(int main_sheets_id) {
    //static bool flag = false;
    if (!m_flag) {
        all_sheets_connectivities = sheets_connectivities;
        m_flag = true;
    }
    size_t n = sheets_componentCellIds.size();
    std::vector<std::vector<size_t>> new_sheets_connectivities;
    std::vector<size_t> sheetIds;
    std::vector<bool> active(n, false);
    for (auto id : sheets_coverSheetIds[main_sheets_id])
        active[id] = true;
    for (size_t i = 0; i < n; ++i)
        if (active[i]) new_sheets_connectivities.push_back(all_sheets_connectivities[i]);

    for (auto& row : new_sheets_connectivities) {
        std::vector<size_t> new_sheet_connectivities;
        for (size_t i = 0; i < n; ++i)
            if (active[i]) new_sheet_connectivities.push_back(row[i]);
        row = new_sheet_connectivities;
    }
    sheets_connectivities = new_sheets_connectivities;

//    for (size_t i = 0; i < n; ++i) {
//        auto& row = sheets_connectivities[i];
//        if (!active[i]) {
//            row.clear();
//            row.resize(n, 0.0f);
//        } else {
//            for (size_t j = 0; j < n; ++j)
//                if (!active[j]) row[j] = 0.0f;
//        }
//    }

    n = sheets_connectivities.size();
    sheets_connectivities_float.resize(n);
    for (auto& row : sheets_connectivities_float) {
        row.clear();
        row.resize(n, 0.0f);
    }
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            if (sheets_connectivities[i][j] == SheetsConnectivity_NEIGHBOR) sheets_connectivities_float[i][j] = -1.0f;
            else if (sheets_connectivities[i][j] == SheetsConnectivity_INTERSECT) sheets_connectivities_float[i][j] = 1.0f;
            else if (sheets_connectivities[i][j] == SheetsConnectivity_NEIGHBOR_AND_INTERSECT) sheets_connectivities_float[i][j] = 0.5f;
}

void BaseComplexSheet::WriteSheetsConnectivitiesMatrixVTK(const char *filename) const {
    MeshFileWriter writer(baseComplex.mesh, filename);
    writer.WriteMatrixVTK(sheets_connectivities);
}

void BaseComplexSheet::RemoveSheetSetsRedundancy() {
    for (auto& coverSheetIds : sheets_coverSheetIds)
        RemoveSheetSetsRedundancy(coverSheetIds);
    std::unordered_set<std::string> hash_sheets_coverSheetIds;
    std::vector<std::vector<size_t>> res;
    for (auto& coverSheetIds : sheets_coverSheetIds) {
        std::set<size_t> s(coverSheetIds.begin(), coverSheetIds.end());
        std::string key;
        for (auto v : s)
            key += "@" + std::to_string(v);
        auto oldsize = hash_sheets_coverSheetIds.size();
        hash_sheets_coverSheetIds.insert(key);
        if (hash_sheets_coverSheetIds.size() != oldsize)
            res.push_back(coverSheetIds);
    }
    sheets_coverSheetIds = res;
}

void BaseComplexSheet::WriteSheetsConnectivitiesMatrixMat(const char* filename) const {
    std::ofstream ofs(filename);
    const size_t n = sheets_connectivities_float.size();
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j)
            ofs << sheets_connectivities_float[i][j] << "\t";
        ofs << "\n";
    }
}

void BaseComplexSheet::WriteAllDominantSheetsConnectivitiesMatrixMat(const char* filename_prefix) {
    int id = 0;
    for (auto& sheet_coverSheetIds : sheets_coverSheetIds) {
        std::string filename = filename_prefix + std::to_string(id) + ".mat";
        ExtractMainSheetConnectivities(id++);
        WriteSheetsConnectivitiesMatrixMat(filename.c_str());
    }
}

void BaseComplexSheet::WriteDominantSheetsConnectivitiesMatrixMat(const char* filename) const {
    std::ofstream ofs(filename);
    const size_t n = sheets_connectivities_float.size();
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j)
            ofs << sheets_connectivities_float[i][j] << "\t";
        ofs << "\n";
    }
}

void BaseComplexSheet::RemoveSheetSetsRedundancy(std::vector<size_t>& coverSheetIds) {
    auto sheetIds = coverSheetIds;
    std::sort(sheetIds.begin(), sheetIds.end(), [&](const size_t& a, const size_t& b) {
        return sheets_componentCellIds[a].size() > sheets_componentCellIds[b].size();});
    std::vector<size_t> componentCellId_count(baseComplex.componentC.size(), 0);
    for (auto sheetId : sheetIds)
        for (auto componentCellId : sheets_componentCellIds[sheetId])
            ++componentCellId_count[componentCellId];
    while (true) {
        bool sheetIdsChanged = false;
        for (auto iter = sheetIds.begin(); iter != sheetIds.end();) {
            auto sheetId = *iter;
            bool canRemove = true;
            for (auto componentCellId : sheets_componentCellIds[sheetId])
                if (componentCellId_count[componentCellId] <= 1) {
                    canRemove = false;
                    break;
                }
            if (canRemove) {
                for (auto componentCellId : sheets_componentCellIds[sheetId])
                    --componentCellId_count[componentCellId];
                iter = sheetIds.erase(iter);
                sheetIdsChanged = true;
            } else ++iter;
        }
        if (!sheetIdsChanged) break;
    }
    std::vector<size_t> new_sheetIds;
    new_sheetIds.reserve(sheetIds.size());
    for (auto coverSheetId : coverSheetIds)
        for (auto sheetId : sheetIds) {
            if (coverSheetId == sheetId) {
                new_sheetIds.push_back(coverSheetId);
                break;
            }
        }

    coverSheetIds = new_sheetIds;
}

std::unordered_set<size_t> BaseComplexSheet::GetSheetBoundaryFaceComponentIds(size_t sheetId) const {
    std::unordered_set<size_t> sheetBoundaryFaceComponentIds;
    std::unordered_set<size_t> sheetCellComponentIds(sheets_componentCellIds[sheetId].begin(), sheets_componentCellIds[sheetId].end());
    std::unordered_set<size_t> sheetFaceComponentIds;
    for (auto componentId : sheetCellComponentIds) {
        const auto& component = baseComplex.componentC.at(componentId);
        sheetFaceComponentIds.insert(component.Fids.begin(), component.Fids.end());
    }
    for (auto faceComponentId : sheetFaceComponentIds) {
        const auto& faceComponent = baseComplex.componentF.at(faceComponentId);
        if (faceComponent.isBoundary) {
            sheetBoundaryFaceComponentIds.insert(faceComponentId);
            continue;
        }
        for (auto neighborCellComponentId : faceComponent.N_Cids)
            if (sheetCellComponentIds.find(neighborCellComponentId) == sheetCellComponentIds.end())
                sheetBoundaryFaceComponentIds.insert(faceComponentId);
    }
    return sheetBoundaryFaceComponentIds;
}

bool isSheetIdExistedInResult(size_t sheetId, const std::vector<size_t>& res) {
    bool existed = false;
    for (auto id : res)
        if (id == sheetId) {
            existed = true;
            break;
        }
    return existed;
}

bool BaseComplexSheet::IsSheetRedundant(size_t sheetId, const std::vector<bool>& componentCovered) const{
    bool coverAll = true;
    for (size_t componentId : sheets_componentCellIds.at(sheetId))
        if (!componentCovered[componentId]) {
            coverAll = false;
            break;
        }
    return coverAll;
}

bool hasCoveredAllComponents(size_t sheetId, const std::vector<bool>& componentCovered, const BaseComplex& baseComplex, bool& coveredAllComponent) {
    bool coverAll = true;
    for (size_t componentId = 0; componentId < baseComplex.componentC.size(); ++componentId)
        if (!componentCovered[componentId]) {
            coverAll = false;
            break;
        }
    if (coverAll) coveredAllComponent = true;
    else coveredAllComponent = false;
    return coveredAllComponent;
}

bool BaseComplexSheet::HasCoveredSheetComponents(size_t sheetId, const std::vector<bool>& componentCovered) const {
    bool coveredSheetComponent = true;
    for (auto component : sheets_componentCellIds[sheetId])
        if (!componentCovered[component]) {
            coveredSheetComponent = false;
            break;
        }
    return coveredSheetComponent;
}

//std::vector<size_t> BaseComplexSheet::GetCoverSheetIds(size_t beginSheetId) const {
//    std::vector<size_t> res;
//    std::unordered_set<size_t> resSet;
//    std::queue<size_t> q;
//    q.push(beginSheetId);
//    std::vector<size_t> componentCovered(baseComplex.componentC.size(), false);
//    bool coveredAllComponent = false;
//    while (!q.empty()) {
//        size_t n = q.size();
//        for (size_t i = 0; i < n; ++i) {
//            auto sheetId = q.front();
//            q.pop();
//            if (isSheetIdExistedInResult(sheetId, res)) continue;
//            res.push_back(sheetId);
//            resSet.insert(sheetId);
//
//            for (auto component : sheets_componentCellIds[sheetId])
//                componentCovered[component] = true;
//            if (hasCoveredAllComponents(sheetId, componentCovered, baseComplex, coveredAllComponent)) break;
//
//            std::unordered_set<size_t> neighborSheetIds = GetNeighborSheetIds(sheetId);
//            q.push(GetMinComponentIntersectionNeighborSheetId(sheetId, neighborSheetIds, resSet));
////            auto candidateSheetIds = GetMinComponentIntersectionNeighborSheetIds(sheetId, neighborSheetIds, resSet);
////            for (auto candidateSheetId : candidateSheetIds)
////                q.push(candidateSheetId);
//        }
//        if (coveredAllComponent) break;
//    }
//    if (!coveredAllComponent)
//        std::cerr << red_color << "Decomposition incorrect in GetCoverSheetIds for begining with " << beginSheetId << "\n" << end_color;
//    return res;
//}

//std::vector<size_t> BaseComplexSheet::GetCoverSheetIds(size_t beginSheetId) const {
//    std::vector<size_t> res;
//    std::unordered_set<size_t> resSet;
//    std::queue<size_t> q;
//    std::stack<size_t> st;
//    q.push(beginSheetId);
//    std::vector<size_t> componentCovered(baseComplex.componentC.size(), false);
//    bool coveredAllComponent = false;
//    while (!coveredAllComponent) {
//        while (!q.empty()) {
//            size_t n = q.size();
//            for (size_t i = 0; i < n; ++i) {
//                auto sheetId = q.front();
//                q.pop();
//                if (isSheetIdExistedInResult(sheetId, res)) continue;
//                res.push_back(sheetId);
//                resSet.insert(sheetId);
//
//                for (auto component : sheets_componentCellIds[sheetId])
//                    componentCovered[component] = true;
//                if (hasCoveredAllComponents(sheetId, componentCovered, baseComplex, coveredAllComponent)) break;
//
//                std::unordered_set<size_t> neighborSheetIds = GetNeighborSheetIds(sheetId);
//                auto target = GetMinComponentIntersectionNeighborSheetId(sheetId, neighborSheetIds, resSet);
//                q.push(target);
//                for (auto neighborSheetId : neighborSheetIds)
//                    if (neighborSheetId != target) st.push(neighborSheetId);
//            }
//            if (coveredAllComponent) break;
//        }
//        if (!coveredAllComponent) {
//            q.push(st.top());
//            st.pop();
//        }
//    }
//    if (!coveredAllComponent)
//        std::cerr << red_color << "Decomposition incorrect in GetCoverSheetIds for begining with " << beginSheetId << "\n" << end_color;
//    return res;
//}

std::vector<size_t> BaseComplexSheet::GetCoverSheetIds(size_t beginSheetId) const {
    std::vector<size_t> res;
    std::unordered_set<size_t> resSet;
    std::queue<size_t> q;
    std::stack<std::pair<size_t, std::unordered_set<size_t>>> st; // sheetid_neigboringSheetIds
    q.push(beginSheetId);
    std::vector<bool> componentCovered(baseComplex.componentC.size(), false);
    bool coveredAllComponent = false;
    while (!coveredAllComponent) {
        while (!q.empty()) {
            size_t n = q.size();
            for (size_t i = 0; i < n; ++i) {
                auto sheetId = q.front();
                q.pop();
                if (isSheetIdExistedInResult(sheetId, res)) continue;
                if (IsSheetRedundant(sheetId, componentCovered)) continue;
                res.push_back(sheetId);
                resSet.insert(sheetId);

                for (auto component : sheets_componentCellIds[sheetId])
                    componentCovered[component] = true;
                if (hasCoveredAllComponents(sheetId, componentCovered, baseComplex, coveredAllComponent)) break;

                std::unordered_set<size_t> neighborSheetIds = GetNeighborSheetIds(sheetId);
                auto target = GetMinComponentIntersectionNeighborSheetId(sheetId, neighborSheetIds, resSet);
                q.push(target);
                neighborSheetIds.erase(target);
                if (!neighborSheetIds.empty())
                	st.push(std::make_pair(sheetId, neighborSheetIds));
            }
            if (coveredAllComponent) break;
        }
        if (!coveredAllComponent) {
        	auto& t = st.top();
            auto target = GetMinComponentIntersectionNeighborSheetId(t.first, t.second, resSet);
            q.push(target);
            t.second.erase(target);
            if (t.second.empty())
            	st.pop();
        }
    }
    if (!coveredAllComponent)
        std::cerr << red_color << "Decomposition incorrect in GetCoverSheetIds for begining with " << beginSheetId << "\n" << end_color;
    return res;
}

std::vector<size_t> BaseComplexSheet::GetCoverSheetIdsBFS(size_t beginSheetId) const {
    std::vector<size_t> res;
    std::unordered_set<size_t> resSet;
    std::queue<size_t> q;
    std::stack<std::pair<size_t, std::unordered_set<size_t>>> st; // sheetid_neigboringSheetIds
    q.push(beginSheetId);
    std::vector<bool> componentCovered(baseComplex.componentC.size(), false);
    bool coveredAllComponent = false;
    while (!coveredAllComponent) {
        while (!q.empty()) {
            size_t n = q.size();
            for (size_t i = 0; i < n; ++i) {
                auto sheetId = q.front();
                q.pop();
                if (isSheetIdExistedInResult(sheetId, res)) continue;
                if (HasCoveredSheetComponents(sheetId, componentCovered)) continue;
                res.push_back(sheetId);
                resSet.insert(sheetId);

                for (auto component : sheets_componentCellIds[sheetId])
                    componentCovered[component] = true;
                if (hasCoveredAllComponents(sheetId, componentCovered, baseComplex, coveredAllComponent)) break;

                std::unordered_set<size_t> neighborSheetIds = GetNeighborSheetIds(sheetId);
                auto candidateSheetIds = GetMinComponentIntersectionNeighborSheetIds(sheetId, neighborSheetIds, resSet);
                for (auto candidateSheetId : candidateSheetIds) {
                    if (!isSheetIdExistedInResult(candidateSheetId, res)) {
                        q.push(candidateSheetId);
                    }
                    neighborSheetIds.erase(candidateSheetId);
                }
                if (!neighborSheetIds.empty())
                	st.push(std::make_pair(sheetId, neighborSheetIds));
            }
            if (!coveredAllComponent) {
            	auto& t = st.top();
                std::unordered_set<size_t>& neighborSheetIds = t.second;
                auto candidateSheetIds = GetMinComponentIntersectionNeighborSheetIds(t.first, neighborSheetIds, resSet);
                for (auto candidateSheetId : candidateSheetIds) {
                    if (!isSheetIdExistedInResult(candidateSheetId, res)) {
                        q.push(candidateSheetId);
                    }
                    neighborSheetIds.erase(candidateSheetId);
                }
                if (neighborSheetIds.empty())
                	st.pop();
            }
        }
    }
    if (!coveredAllComponent)
        std::cerr << red_color << "Decomposition incorrect in GetCoverSheetIdsBFS for begining with " << beginSheetId << "\n" << end_color;
    return res;
}

std::unordered_set<size_t> BaseComplexSheet::GetNeighborSheetIds(size_t sheetId) const {
    std::unordered_set<size_t> neighborSheetIds;
    for (auto& boundaryFaceComponentId : sheets_BoundaryFaceComponentIds[sheetId]) {
        auto& boundaryFaceComponent = baseComplex.componentF.at(boundaryFaceComponentId);
        for (auto neiborghComponentId : boundaryFaceComponent.N_Cids) {
            for (auto neighborSheetId : cellComponent_sheetIds[neiborghComponentId]) {
                if (neighborSheetId == sheetId) continue;
                neighborSheetIds.insert(neighborSheetId);
            }
        }
    }
    return neighborSheetIds;
}

size_t BaseComplexSheet::GetMinComponentIntersectionNeighborSheetId(size_t sheetId,
        const std::unordered_set<size_t>& neighborSheetIds, const std::unordered_set<size_t>& resSet) const {
    size_t minComponentIntersectionNums = MAXID;
    size_t minComponentIntersectionNeighborSheetId = *neighborSheetIds.begin();
    std::unordered_set<size_t> sheetComponentIds(sheets_componentCellIds[sheetId].begin(), sheets_componentCellIds[sheetId].end());
    for (auto neighborSheetId : neighborSheetIds) {
        if (resSet.find(neighborSheetId) != resSet.end()) continue;
        std::unordered_set<size_t> neighborSheetComponentIds(sheets_componentCellIds[neighborSheetId].begin(), sheets_componentCellIds[neighborSheetId].end());
        size_t intersectionNum = 0;
        for (auto id : sheetComponentIds)
            if (neighborSheetComponentIds.find(id) != neighborSheetComponentIds.end()) ++intersectionNum;
        if (intersectionNum < minComponentIntersectionNums ||
                (intersectionNum == minComponentIntersectionNums &&
                        sheets_componentCellIds[neighborSheetId].size() > sheets_componentCellIds[minComponentIntersectionNums].size())) {
            minComponentIntersectionNums = intersectionNum;
            minComponentIntersectionNeighborSheetId = neighborSheetId;
        }
    }
    return minComponentIntersectionNeighborSheetId;
}

std::vector<size_t> BaseComplexSheet::GetMinComponentIntersectionNeighborSheetIds(size_t sheetId,
        const std::unordered_set<size_t>& neighborSheetIds, const std::unordered_set<size_t>& resSet) const {
    std::vector<size_t> res;
    size_t minComponentIntersectionNums = MAXID;
    size_t minComponentIntersectionNeighborSheetId = *neighborSheetIds.begin();
    std::unordered_set<size_t> sheetComponentIds(sheets_componentCellIds[sheetId].begin(), sheets_componentCellIds[sheetId].end());
    //const auto& sheetComponentIds = sheets_hashComponentCellIds.at(sheetId);
    for (auto neighborSheetId : neighborSheetIds) {
        if (resSet.find(neighborSheetId) != resSet.end()) continue;
        std::unordered_set<size_t> neighborSheetComponentIds(sheets_componentCellIds[neighborSheetId].begin(), sheets_componentCellIds[neighborSheetId].end());
        //const auto& neighborSheetComponentIds = sheets_hashComponentCellIds.at(neighborSheetId);
        size_t intersectionNum = 0;
        for (auto id : sheetComponentIds)
            if (neighborSheetComponentIds.find(id) != neighborSheetComponentIds.end()) ++intersectionNum;
        if (intersectionNum < minComponentIntersectionNums ||
                (intersectionNum == minComponentIntersectionNums &&
                        sheets_componentCellIds[neighborSheetId].size() > sheets_componentCellIds[minComponentIntersectionNums].size())) {
            minComponentIntersectionNums = intersectionNum;
            minComponentIntersectionNeighborSheetId = neighborSheetId;
        }
    }
    for (auto neighborSheetId : neighborSheetIds) {
        if (resSet.find(neighborSheetId) != resSet.end()) continue;
        std::unordered_set<size_t> neighborSheetComponentIds(sheets_componentCellIds[neighborSheetId].begin(), sheets_componentCellIds[neighborSheetId].end());
        //const auto& neighborSheetComponentIds = sheets_hashComponentCellIds.at(neighborSheetId);
        size_t intersectionNum = 0;
        for (auto id : sheetComponentIds)
            if (neighborSheetComponentIds.find(id) != neighborSheetComponentIds.end()) ++intersectionNum;
        if (intersectionNum == minComponentIntersectionNums) res.push_back(neighborSheetId);
    }
    return res;
}

void BaseComplexSheet::WriteSheetNeighborsSheetIdsJS(const char* filename) {
    std::ofstream ofs(filename);
    ofs << "var data = [";
    for (size_t sheetId = 0; sheetId < sheets_componentCellIds.size(); ++sheetId) {
        //if (sheetId != 0) ofs << ",";
        std::unordered_set<size_t> neighborSheetIds = GetNeighborSheetIds(sheetId);
        std::unordered_set<size_t> sheetComponentIds(sheets_componentCellIds[sheetId].begin(), sheets_componentCellIds[sheetId].end());
        for (auto neighborSheetId : neighborSheetIds) {
            std::unordered_set<size_t> neighborSheetComponentIds(sheets_componentCellIds[neighborSheetId].begin(), sheets_componentCellIds[neighborSheetId].end());
            size_t intersectionNum = 0;
            for (auto id : sheetComponentIds)
                if (neighborSheetComponentIds.find(id) != neighborSheetComponentIds.end()) ++intersectionNum;
            ofs << ",[\'s" << sheetId << "\',\'s" << neighborSheetId << "\'," << intersectionNum << "]\n";
        }
    }
    ofs << "];";
}

//void BaseComplexSheet::ExtractSets()
//{
//    std::vector<std::set<size_t>> input;
//    std::set<size_t> target;
//    for (auto& i : sheets_componentCellIds)
//        input.push_back(std::set<size_t>(i.begin(), i.end()));
//    for (auto i = 0; i < baseComplex.componentC.size(); target.insert(i++));
//    auto output_set = FindSetCombination(input, target);
//    for (auto& output : output_set) {
//        std::set<size_t> ids;
//        for (auto i : output)
//            for (auto j = 0; j < input.size(); ++j)
//                if (i == input[j]) {
//                    ids.insert(j);
//                    break;
//                }
//        std::cout << "output size " << output.size() << " *************" << std::endl;
//        for (auto i : ids)
//            std::cout << i << " ";
//        std::cout << std::endl;
//        std::cout << "***************************" << std::endl;
//    }
//}

void BaseComplexSheet::GetParallelComponents(const ComponentEdge & componentEdge,
        std::vector<size_t>& sheetComponentEdgeIds, std::vector<size_t>& sheetComponentFaceIds, std::vector<size_t>& sheetComponentCellIds)
{
    std::vector<bool> edge_visited(baseComplex.componentE.size(), false);
    std::vector<bool> face_visited(baseComplex.componentF.size(), false);
    std::vector<bool> cell_visited(baseComplex.componentC.size(), false);
    edge_visited.at(componentEdge.id) = true;
    sheetComponentEdgeIds.push_back(componentEdge.id);
    std::queue<size_t> q;
    q.push(componentEdge.id);
    while (!q.empty()) {
        size_t n = q.size();
        for (int i = 0; i < n; ++i) {
            const size_t componentEdgeId = q.front();
            q.pop();

            //if (!edge_visited.at(componentEdgeId)) sheetComponentEdgeIds.push_back(componentEdgeId);
            for (auto neighborComponentFaceId : baseComplex.componentE.at(componentEdgeId).N_Fids)
                if (!face_visited.at(neighborComponentFaceId)) {
                    face_visited.at(neighborComponentFaceId) = true;
                    sheetComponentFaceIds.push_back(neighborComponentFaceId);
                }
            for (auto neighborComponentCellId : baseComplex.componentE.at(componentEdgeId).N_Cids)
                if (!cell_visited.at(neighborComponentCellId)) {
                    cell_visited.at(neighborComponentCellId) = true;
                    sheetComponentCellIds.push_back(neighborComponentCellId);
                }

            const auto nextComponentEdgeIds = GetParallelComponentEdgeIds(baseComplex.componentE.at(componentEdgeId));
            for (auto nextComponentEdgeId : nextComponentEdgeIds)
                if (!edge_visited.at(nextComponentEdgeId)) {
                    edge_visited.at(nextComponentEdgeId) = true;
                    sheetComponentEdgeIds.push_back(nextComponentEdgeId);
                    q.push(nextComponentEdgeId);
                }
        }
    }
}
/*
         3____________________2
         /|                 /|
        / |                / |
       /  |               /  |
   0  /___|______________/ 1 |
      |   |              |   |
      |   |              |   |
      |   |              |   |
      |   |______________|___|
      |   / 7            |  / 6
      |  /               | /
      | /                |/
      |/_________________/
    4                   5
*/
std::vector<size_t> BaseComplexSheet::GetParallelComponentEdgeIds(const ComponentEdge & componentEdge) {
    std::vector<size_t> res;
    for (auto neighborComponentFaceId : baseComplex.componentE.at(componentEdge.id).N_Fids) {
        for (auto componentFaceEdgeId : baseComplex.componentF.at(neighborComponentFaceId).Eids) {
            auto& componentFaceEdge = baseComplex.componentE.at(componentFaceEdgeId);
            if (componentFaceEdge.vids_link.front() != componentEdge.vids_link.front() &&
                componentFaceEdge.vids_link.front() != componentEdge.vids_link.back() &&
                componentFaceEdge.vids_link.back() != componentEdge.vids_link.front() &&
                componentFaceEdge.vids_link.back() != componentEdge.vids_link.back()) {
                res.push_back(componentFaceEdgeId);
            }
        }
    }
    return res;
}

void BaseComplexSheet::WriteSheetsEdgesVTK(const char *filename) const
{

}

void BaseComplexSheet::WriteSheetsFacesVTK(const char *filename) const
{

}

void BaseComplexSheet::WriteSheetsCellsVTK(const char *filename) const
{
    const auto& mesh = baseComplex.GetMesh();
    const auto& V = mesh.V;
    const auto& E = mesh.E;
    std::ofstream ofs(filename);
    ofs << "# vtk DataFile Version 2.0\n"
        << filename << "\n"
        << "ASCII\n\n"
        << "DATASET UNSTRUCTURED_GRID\n";
    ofs << "POINTS " << V.size() << " float" << "\n";
    for (size_t i = 0; i < V.size(); i++)
        ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << "\n";

    size_t cells_num = 0;
    for (const auto& sheet : sheets_componentCellIds)
        for (const auto componentid : sheet)
            cells_num += baseComplex.componentC.at(componentid).cids_patch.size();

    ofs << "CELLS " << cells_num << " " << 9 * cells_num << "\n";
    for (const auto& sheet : sheets_componentCellIds)
        for (const auto componentid : sheet)
            for (const auto cell_id : baseComplex.componentC.at(componentid).cids_patch) {
                    const Cell& cell = mesh.C.at(cell_id);
                    ofs << cell.Vids.size();
                    for (auto vid : cell.Vids)
                        ofs << " " << vid;
                    ofs << "\n";
            }

    ofs << "CELL_TYPES " << cells_num << "\n";
    for (size_t i = 0; i < cells_num; i++)
        ofs << "12\n";

    ofs << "CELL_DATA " << cells_num << "\n"
        << "SCALARS " << "sheet" << " int 1\n"
        << "LOOKUP_TABLE default\n";
    size_t sheet_id = 0;
    for (const auto& sheet : sheets_componentCellIds) {
        for (const auto componentid : sheet)
            for (const auto cell_id : baseComplex.componentC.at(componentid).cids_patch)
                ofs << sheet_id << "\n";
        ++sheet_id;
    }
}

void BaseComplexSheet::WriteAllSheetsCellsDualVTK(const char *filename_prefix) const
{
    RefinedDual dual(baseComplex.mesh);
    dual.Build();

    for (int i = 0; i < sheets_componentCellIds.size(); ++i) {
        std::string filename = std::string(filename_prefix) + std::to_string(i) + ".vtk";
        WriteSheetCellsDualVTK(filename.c_str(), i, dual);
    }
}

void BaseComplexSheet::WriteAllSheetsCellsVTK(const char *filename_prefix) const
{
    for (int i = 0; i < sheets_componentCellIds.size(); ++i) {
        std::string filename = std::string(filename_prefix) + std::to_string(i) + ".vtk";
        WriteSheetCellsVTK(filename.c_str(), i);
    }
}

void BaseComplexSheet::WriteAllSheetsFacesVTK(const char *filename_prefix) const
{
    for (int i = 0; i < sheets_componentFaceIds.size(); ++i) {
        std::string filename = std::string(filename_prefix) + std::to_string(i) + ".vtk";
        WriteSheetFacesVTK(filename.c_str(), i);
    }
}

void BaseComplexSheet::WriteAllSheetsEdgesVTK(const char *filename_prefix) const
{
    for (int i = 0; i < sheets_componentEdgeIds.size(); ++i) {
        std::string filename = std::string(filename_prefix) + std::to_string(i) + ".vtk";
        WriteSheetEdgesVTK(filename.c_str(), i);
    }
}

void BaseComplexSheet::WriteAllSheetsFacesAndEdgesVTK(const char *filename_prefix) const
{
    for (int i = 0; i < sheets_componentEdgeIds.size(); ++i) {
        std::string filename = std::string(filename_prefix) + std::to_string(i) + ".vtk";
        WriteSheetFacesAndEdgesVTK(filename.c_str(), i);
    }
}

std::unordered_set<size_t> BaseComplexSheet::GetParallelEdgeIds(const size_t sheet_id) const {
    std::unordered_set<size_t> allEdgeIds;
    for (const auto componentid : sheets_componentCellIds.at(sheet_id))
        for (const auto cellid : baseComplex.componentC.at(componentid).cids_patch) {
            const auto& cell = baseComplex.mesh.C.at(cellid);
            for (const auto edgeid : cell.Eids)
                if (allEdgeIds.find(edgeid) == allEdgeIds.end()) allEdgeIds.insert(edgeid);
        }

    std::unordered_set<size_t> parallelEdgeIds;
    size_t edges_num = 0;
    for (const auto componentid : sheets_componentEdgeIds.at(sheet_id)) {
        //for (const auto edgeid : baseComplex.componentE.at(componentid).eids_link)
        size_t mid = baseComplex.componentE.at(componentid).eids_link.size() / 2;
        size_t edgeid = baseComplex.componentE.at(componentid).eids_link.at(mid);
            parallelEdgeIds.insert(edgeid);
            break;
//        parallelEdgeIds.insert(baseComplex.componentE.at(componentid).eids_link.front());
//        if (baseComplex.componentE.at(componentid).eids_link.size() > 1)
//        parallelEdgeIds.insert(baseComplex.componentE.at(componentid).eids_link.back());
    }
    size_t old = parallelEdgeIds.size();
    while (true) {
        for (auto edgeid : parallelEdgeIds) {
            const auto & edge = baseComplex.mesh.E.at(edgeid);
            for (auto parallel_edge : edge.parallelEids)
                if (allEdgeIds.find(parallel_edge) != allEdgeIds.end() &&  parallelEdgeIds.find(parallel_edge) == parallelEdgeIds.end())
                    parallelEdgeIds.insert(parallel_edge);
        }
        if (parallelEdgeIds.size() == old) break;
        old = parallelEdgeIds.size();
    }
    return parallelEdgeIds;
}

std::unordered_set<size_t> BaseComplexSheet::GetParallelSingularEdgeIds(const size_t sheet_id) const {
    std::unordered_set<size_t> parallelSingularEdgeIds;
    for (const auto componentEdgeId : sheets_componentEdgeIds.at(sheet_id))
        for (const auto edgeid : baseComplex.componentE.at(componentEdgeId).eids_link) {
            const auto& edge = baseComplex.mesh.E.at(edgeid);
            if (edge.isSingularity) parallelSingularEdgeIds.insert(edgeid);
        }
    return parallelSingularEdgeIds;
}

void BaseComplexSheet::WriteSheetCellsDualVTK(const char *filename, const size_t sheet_id) const
{
    auto parallelEdgeIds = GetParallelEdgeIds(sheet_id);
    Dual dual(baseComplex.mesh);
    dual.Build();

    std::unordered_set<size_t> dualFaceIds;
    for (const auto componentid : sheets_componentCellIds.at(sheet_id))
        for (const auto cellid : baseComplex.componentC.at(componentid).cids_patch) {
            const auto& cell = baseComplex.mesh.C.at(cellid);
            auto cellParallelEdgeIds = getParallelEdgeIds(baseComplex.mesh, cell);
            for (auto i = 0; i < 3; ++i) {
                bool canFindFourEdges = true;
                for (auto parallelEdgeId : cellParallelEdgeIds.at(i))
                    if (parallelEdgeIds.find(parallelEdgeId) == parallelEdgeIds.end()) {
                        canFindFourEdges = false;
                        break;
                    }
                if (canFindFourEdges) {
                    dualFaceIds.insert(cell.id * 3 + i);
                }
            }
        }

    const auto& mesh = baseComplex.GetMesh();
    const auto& V = dual.V;
    const auto& F = dual.F;
    std::ofstream ofs(filename);
    ofs << "# vtk DataFile Version 2.0\n"
        << filename << "\n"
        << "ASCII\n\n"
        << "DATASET POLYDATA\n";
    ofs << "POINTS " << V.size() << " double" << "\n";
    for (const auto& v : V)
        ofs << v.x << " " << v.y << " " << v.z << "\n";

    ofs << "Polygons " << dualFaceIds.size() << " " << 5 * dualFaceIds.size() << "\n";
    for (const auto face_id : dualFaceIds) {
        const auto& face = F.at(face_id);
        ofs << face.Vids.size();
        for (auto vid : face.Vids)
            ofs << " " << vid;
        ofs << "\n";
    }

    ofs << "CELL_DATA " << dualFaceIds.size() << "\n"
        << "SCALARS " << "sheet" << " int 1\n"
        << "LOOKUP_TABLE default\n";
    for (const auto id : dualFaceIds)
        ofs << sheet_id << "\n";
}

std::unordered_set<size_t> BaseComplexSheet::GetDualFaceIds(const size_t sheet_id) const {
    auto parallelEdgeIds = GetParallelEdgeIds(sheet_id);
    std::unordered_set<size_t> dualFaceIds;
    for (const auto componentid : sheets_componentCellIds.at(sheet_id))
        for (const auto cellid : baseComplex.componentC.at(componentid).cids_patch) {
            const auto& cell = baseComplex.mesh.C.at(cellid);
            auto cellParallelEdgeIds = getParallelEdgeIds(baseComplex.mesh, cell);
            for (auto i = 0; i < 3; ++i) {
                bool canFindFourEdges = true;
                for (auto parallelEdgeId : cellParallelEdgeIds.at(i))
                    if (parallelEdgeIds.find(parallelEdgeId) == parallelEdgeIds.end()) {
                        canFindFourEdges = false;
                        break;
                    }
                if (canFindFourEdges)
                    for (auto j = 0; j < 4; ++j)
                        dualFaceIds.insert(baseComplex.mesh.F.size() * 4 + cell.id * 12 + i * 4 + j);
            }
        }
    return dualFaceIds;
}

void BaseComplexSheet::WriteSheetCellsDualVTK(const char *filename, const size_t sheet_id, const RefinedDual& dual) const {
    std::unordered_set<size_t> dualFaceIds = GetDualFaceIds(sheet_id);
    std::unordered_set<size_t> parallelSingularEdgeIds = GetParallelSingularEdgeIds(sheet_id);
    WriteSheetCellsDualVTK(filename, sheet_id, dual, dualFaceIds);
    //WriteSheetCellsDualVTK(filename, sheet_id, dual, dualFaceIds, parallelSingularEdgeIds);
}

void BaseComplexSheet::WriteSheetCellsDualVTK(const char *filename, const size_t sheet_id, const RefinedDual& dual, const std::unordered_set<size_t>& dualFaceIds) const
{
    const auto& mesh = baseComplex.GetMesh();
    const auto& V = dual.V;
    const auto& F = dual.F;
    std::ofstream ofs(filename);
    ofs << "# vtk DataFile Version 2.0\n"
        << filename << "\n"
        << "ASCII\n\n"
        << "DATASET POLYDATA\n";
    ofs << "POINTS " << V.size() << " double" << "\n";
    for (const auto& v : V)
        ofs << v.x << " " << v.y << " " << v.z << "\n";

    ofs << "Polygons " << dualFaceIds.size() << " " << 5 * dualFaceIds.size() << "\n";
    for (const auto face_id : dualFaceIds) {
        const auto& face = F.at(face_id);
        ofs << face.Vids.size();
        for (auto vid : face.Vids)
            ofs << " " << vid;
        ofs << "\n";
    }

    ofs << "CELL_DATA " << dualFaceIds.size() << "\n"
        << "SCALARS " << "sheet" << " int 1\n"
        << "LOOKUP_TABLE default\n";
    for (const auto id : dualFaceIds)
        ofs << sheet_id << "\n";
}

void BaseComplexSheet::WriteSheetCellsDualVTK(const char *filename, const size_t sheet_id, const RefinedDual& dual,
        const std::unordered_set<size_t>& dualFaceIds, const std::unordered_set<size_t>& parallelSingularEdgeIds) const
{
    const auto& mesh = baseComplex.GetMesh();
    const auto& V = dual.V;
    const auto& F = dual.F;
    std::ofstream ofs(filename);
    ofs << "# vtk DataFile Version 2.0\n"
        << filename << "\n"
        << "ASCII\n\n"
        << "DATASET POLYDATA\n";
    ofs << "POINTS " << V.size() << " double" << "\n";
    for (const auto& v : V)
        ofs << v.x << " " << v.y << " " << v.z << "\n";

    ofs << "Lines " << parallelSingularEdgeIds.size() << " " << 3 * parallelSingularEdgeIds.size() << "\n";
    for (const auto edge_id : parallelSingularEdgeIds) {
            const Edge& cell = mesh.E.at(edge_id);
            ofs << cell.Vids.size();
            for (auto vid : cell.Vids)
                ofs << " " << vid;
            ofs << "\n";
    }

    ofs << "Polygons " << dualFaceIds.size() << " " << 5 * dualFaceIds.size() << "\n";
    for (const auto face_id : dualFaceIds) {
        const auto& face = F.at(face_id);
        ofs << face.Vids.size();
        for (auto vid : face.Vids)
            ofs << " " << vid;
        ofs << "\n";
    }

    ofs << "CELL_DATA " << parallelSingularEdgeIds.size() + dualFaceIds.size() << "\n"
        << "SCALARS " << "sheet" << " int 1\n"
        << "LOOKUP_TABLE default\n";
    for (const auto edge_id : parallelSingularEdgeIds)
        ofs << sheet_id << "\n";
    for (const auto id : dualFaceIds)
        ofs << sheet_id << "\n";
}

void BaseComplexSheet::WriteSheetCellsVTK(const char *filename, const size_t sheet_id) const
{
    const auto& mesh = baseComplex.GetMesh();
    const auto& V = mesh.V;
    const auto& E = mesh.E;
    std::ofstream ofs(filename);
    ofs << "# vtk DataFile Version 2.0\n"
        << filename << "\n"
        << "ASCII\n\n"
        << "DATASET UNSTRUCTURED_GRID\n";
    ofs << "POINTS " << V.size() << " float" << "\n";
    for (size_t i = 0; i < V.size(); i++)
        ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << "\n";

    size_t cells_num = 0;
    for (const auto componentid : sheets_componentCellIds.at(sheet_id))
        cells_num += baseComplex.componentC.at(componentid).cids_patch.size();

    ofs << "CELLS " << cells_num << " " << 9 * cells_num << "\n";
    for (const auto componentid : sheets_componentCellIds.at(sheet_id))
        for (const auto cell_id : baseComplex.componentC.at(componentid).cids_patch) {
                const Cell& cell = mesh.C.at(cell_id);
                ofs << cell.Vids.size();
                for (auto vid : cell.Vids)
                    ofs << " " << vid;
                ofs << "\n";
        }

    ofs << "CELL_TYPES " << cells_num << "\n";
    for (size_t i = 0; i < cells_num; i++)
        ofs << "12\n";

    ofs << "CELL_DATA " << cells_num << "\n";
    ofs << "SCALARS " << "component" << " int 1\n"
        << "LOOKUP_TABLE default\n";
    for (const auto componentid : sheets_componentCellIds.at(sheet_id))
        for (const auto cell_id : baseComplex.componentC.at(componentid).cids_patch)
            ofs << baseComplex.componentC.at(componentid).color % 8 << "\n";
    ofs << "SCALARS " << "raw_component_id" << " int 1\n"
        << "LOOKUP_TABLE default\n";
    for (const auto componentid : sheets_componentCellIds.at(sheet_id))
        for (const auto cell_id : baseComplex.componentC.at(componentid).cids_patch)
            ofs << componentid << "\n";
    ofs << "SCALARS " << "component_id" << " int 1\n"
        << "LOOKUP_TABLE default\n";
    int sheet_component_id = 0;
    for (const auto componentid : sheets_componentCellIds.at(sheet_id)) {
        for (const auto cell_id : baseComplex.componentC.at(componentid).cids_patch)
            ofs << sheet_component_id << "\n";
        ++sheet_component_id;
    }
    ofs << "SCALARS " << "sheet" << " int 1\n"
        << "LOOKUP_TABLE default\n";
    for (const auto componentid : sheets_componentCellIds.at(sheet_id))
        for (const auto cell_id : baseComplex.componentC.at(componentid).cids_patch)
            ofs << sheet_id << "\n";
}

void BaseComplexSheet::WriteSheetFacesVTK(const char *filename, const size_t sheet_id) const
{
    const auto& mesh = baseComplex.GetMesh();
    const auto& V = mesh.V;
    std::ofstream ofs(filename);
    ofs << "# vtk DataFile Version 2.0\n"
        << filename << "\n"
        << "ASCII\n\n"
        << "DATASET POLYDATA\n";
    ofs << "POINTS " << V.size() << " double" << "\n";
    for (size_t i = 0; i < V.size(); i++)
        ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << "\n";

    size_t faces_num = 0;
    for (const auto componentid : sheets_componentFaceIds.at(sheet_id))
        faces_num += baseComplex.componentF.at(componentid).fids_patch.size();

    ofs << "Polygons " << faces_num << " " << 5 * faces_num << "\n";
    for (const auto componentid : sheets_componentFaceIds.at(sheet_id))
        for (const auto cell_id : baseComplex.componentF.at(componentid).fids_patch) {
                const Face& cell = mesh.F.at(cell_id);
                ofs << cell.Vids.size();
                for (auto vid : cell.Vids)
                    ofs << " " << vid;
                ofs << "\n";
        }

    ofs << "CELL_DATA " << faces_num << "\n"
        << "SCALARS " << "sheet" << " int 1\n"
        << "LOOKUP_TABLE default\n";
    for (const auto componentid : sheets_componentFaceIds.at(sheet_id))
        for (const auto cell_id : baseComplex.componentF.at(componentid).fids_patch)
            ofs << sheet_id << "\n";
}

void BaseComplexSheet::WriteSheetEdgesVTK(const char *filename, const size_t sheet_id) const
{
    const auto& mesh = baseComplex.GetMesh();
    const auto& V = mesh.V;
    std::ofstream ofs(filename);
    ofs << "# vtk DataFile Version 2.0\n"
        << filename << "\n"
        << "ASCII\n\n"
        << "DATASET POLYDATA\n";
    ofs << "POINTS " << V.size() << " double" << "\n";
    for (size_t i = 0; i < V.size(); i++)
        ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << "\n";

    size_t edges_num = 0;
    for (const auto componentid : sheets_componentEdgeIds.at(sheet_id))
        edges_num += baseComplex.componentE.at(componentid).eids_link.size();

    ofs << "Lines " << edges_num << " " << 3 * edges_num << "\n";
    for (const auto componentid : sheets_componentEdgeIds.at(sheet_id))
        for (const auto cell_id : baseComplex.componentE.at(componentid).eids_link) {
                const Edge& cell = mesh.E.at(cell_id);
                ofs << cell.Vids.size();
                for (auto vid : cell.Vids)
                    ofs << " " << vid;
                ofs << "\n";
        }

    ofs << "CELL_DATA " << edges_num << "\n"
        << "SCALARS " << "sheet" << " int 1\n"
        << "LOOKUP_TABLE default\n";
    for (const auto componentid : sheets_componentEdgeIds.at(sheet_id))
        for (const auto edge_id : baseComplex.componentE.at(componentid).eids_link)
            ofs << sheet_id << "\n";
    ofs << "SCALARS " << "valence" << " int 1\n"
        << "LOOKUP_TABLE default\n";
    for (const auto componentid : sheets_componentEdgeIds.at(sheet_id))
        for (const auto edge_id : baseComplex.componentE.at(componentid).eids_link)
            ofs << baseComplex.mesh.E.at(edge_id).N_Cids.size() << "\n";
}

void BaseComplexSheet::WriteSheetFacesAndEdgesVTK(const char *filename, const size_t sheet_id) const
{
    const auto& mesh = baseComplex.GetMesh();
    const auto& V = mesh.V;
    std::ofstream ofs(filename);
    ofs << "# vtk DataFile Version 2.0\n"
        << filename << "\n"
        << "ASCII\n\n"
        << "DATASET POLYDATA\n";
    ofs << "POINTS " << V.size() << " double" << "\n";
    for (size_t i = 0; i < V.size(); i++)
        ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << "\n";

    size_t faces_num = 0;
        for (const auto componentid : sheets_componentFaceIds.at(sheet_id))
            faces_num += baseComplex.componentF.at(componentid).fids_patch.size();
    size_t edges_num = 0;
        for (const auto componentid : sheets_componentEdgeIds.at(sheet_id))
            edges_num += baseComplex.componentE.at(componentid).eids_link.size();

    ofs << "Lines " << edges_num << " " << 3 * edges_num << "\n";
    for (const auto componentid : sheets_componentEdgeIds.at(sheet_id))
        for (const auto cell_id : baseComplex.componentE.at(componentid).eids_link) {
                const Edge& cell = mesh.E.at(cell_id);
                ofs << cell.Vids.size();
                for (auto vid : cell.Vids)
                    ofs << " " << vid;
                ofs << "\n";
        }
    ofs << "Polygons " << faces_num << " " << 5 * faces_num << "\n";
    for (const auto componentid : sheets_componentFaceIds.at(sheet_id))
        for (const auto cell_id : baseComplex.componentF.at(componentid).fids_patch) {
                const Face& cell = mesh.F.at(cell_id);
                ofs << cell.Vids.size();
                for (auto vid : cell.Vids)
                    ofs << " " << vid;
                ofs << "\n";
        }

    ofs << "CELL_DATA " << edges_num + faces_num << "\n"
        << "SCALARS " << "sheet" << " int 1\n"
        << "LOOKUP_TABLE default\n";
    for (const auto componentid : sheets_componentFaceIds.at(sheet_id))
        for (const auto cell_id : baseComplex.componentF.at(componentid).fids_patch)
            ofs << sheet_id << "\n";
    for (const auto componentid : sheets_componentEdgeIds.at(sheet_id))
        for (const auto cell_id : baseComplex.componentE.at(componentid).eids_link)
            ofs << sheet_id << "\n";
    ofs << "SCALARS " << "valence" << " int 1\n"
        << "LOOKUP_TABLE default\n";
    for (const auto componentid : sheets_componentEdgeIds.at(sheet_id))
        for (const auto edge_id : baseComplex.componentE.at(componentid).eids_link)
            ofs << baseComplex.mesh.E.at(edge_id).N_Cids.size() << "\n";
    for (const auto componentid : sheets_componentFaceIds.at(sheet_id))
        for (const auto cell_id : baseComplex.componentF.at(componentid).fids_patch)
            ofs << 0 << "\n";
}

const std::vector<std::vector<size_t>>& BaseComplexSheet::Get_sheets_coverSheetIds() const {
    return sheets_coverSheetIds;
}

const size_t BaseComplexSheet::GetNumOfSheets() const {
    return sheets_componentFaceIds.size();
}

void BaseComplexSheet::ComputeComplexity() {
    {
    const auto& main_sheet_ids = sheets_coverSheetIds.front();
    auto combinations = combine(main_sheet_ids.size(), 2);
    size_t complexity = 0;
    for (auto& p : combinations) {
        auto sheetid1 = main_sheet_ids[p[0]];
        auto sheetid2 = main_sheet_ids[p[1]];
        complexity += all_sheets_connectivities[sheetid1][sheetid2] != SheetsConnectivity_UNKNOWN;
    }
    //if (complexity == 0) complexity = 1;
    std::cout << "##############################\n";
    std::cout << "### Complexity of main sheets = " << complexity << "\n";
    std::cout << "##############################\n";

    //std::unordered_set<size_t> hash_main_sheet_ids(main_sheet_ids.begin(), main_sheet_ids.end());
    if (complexity == 0) sheets_importance = {0};
    else
    {
    sheets_importance.reserve(main_sheet_ids.size());
    for (auto main_sheet_id : main_sheet_ids) {
        size_t importance = 0;
        for (auto id : main_sheet_ids)
            importance += all_sheets_connectivities[main_sheet_id][id] != SheetsConnectivity_UNKNOWN;
        sheets_importance.push_back((float(importance))/complexity);
    }
    }
    std::cout << "##############################\n";
    std::cout << "### Complexity of each sheet : ";
    for (auto sheet_importance : sheets_importance)
        std::cout << sheet_importance << " ";
    std::cout << "\n##############################\n";

    std::ifstream ifs("sheet_complexity.txt");
    std::vector<size_t> sheets_complexities;
    std::string line;
    while (getline(ifs, line)) {
        std::stringstream ss(line);
        std::string str;
        while (ss >> str);
        sheets_complexities.push_back(std::stoi(str));
    }

    float total_complexity = 0.0f;
    int i = 0;
    for (auto main_sheet_id : main_sheet_ids)
        total_complexity += sheets_complexities[main_sheet_id] * sheets_importance[i++];

    std::cout << "### Total Complexity = " << total_complexity << "\n";
    }
    {
        const auto& main_sheet_ids = sheets_coverSheetIds.front();
        auto combinations = combine(main_sheet_ids.size(), 2);
        size_t complexity = 0;
        for (auto& row : all_sheets_connectivities)
            for (auto& ele : row)
                complexity += ele != SheetsConnectivity_UNKNOWN;
        //if (complexity == 0) complexity = 1;
        std::cout << "##############################\n";
        std::cout << "### Complexity of all sheets = " << complexity << "\n";
        std::cout << "##############################\n";

        //std::unordered_set<size_t> hash_main_sheet_ids(main_sheet_ids.begin(), main_sheet_ids.end());
//        if (complexity == 1) sheets_importance = {1};
//        else
        {
        sheets_importance.reserve(all_sheets_connectivities.size());
        for (auto& row : all_sheets_connectivities) {
            size_t importance = 0;
            for (auto& ele : row)
                importance += ele != SheetsConnectivity_UNKNOWN;
            sheets_importance.push_back((float(importance))/complexity);
        }
        }
        std::cout << "##############################\n";
        std::cout << "### Complexity of each sheet : ";
        for (auto sheet_importance : sheets_importance)
            std::cout << sheet_importance << " ";
        std::cout << "\n##############################\n";

        std::ifstream ifs("sheet_complexity.txt");
        std::vector<size_t> sheets_complexities;
        std::string line;
        while (getline(ifs, line)) {
            std::stringstream ss(line);
            std::string str;
            while (ss >> str);
            sheets_complexities.push_back(std::stoi(str));
        }

        float total_complexity = 0.0f;

        for (int i = 0; i < all_sheets_connectivities.size(); ++i)
            total_complexity += sheets_complexities[i] * sheets_importance[i];

        std::cout << "### All Total Complexity = " << total_complexity << "\n";
    }
}

void BaseComplexSheet::ComputeImportance() {

}

bool BaseComplexSheet::IsAdjacent(const int component_id1, const int component_id2) const {
    const auto& component1 = baseComplex.componentC.at(component_id1);
    const auto& component2 = baseComplex.componentC.at(component_id2);

    std::unordered_set<size_t> s(component1.Fids.begin(), component1.Fids.end());
    s.insert(component2.Fids.begin(), component2.Fids.end());

    return s.size() == 11;
}

std::unordered_set<size_t> BaseComplexSheet::GetCommonComponentCellIds(const std::unordered_set<size_t>& common_component_face_ids) const {
    std::unordered_set<size_t> common_component_cell_ids;
    for (auto common_component_face_id : common_component_face_ids) {
        auto& component = baseComplex.componentF[common_component_face_id];
        common_component_cell_ids.insert(component.N_Cids.begin(), component.N_Cids.end());
    }
    return common_component_cell_ids;
}

std::unordered_set<size_t> BaseComplexSheet::GetCommonComponentFaceIds(size_t sheetid1, size_t sheetid2) const {
    const auto& sheet1_component_face_ids = sheets_BoundaryFaceComponentIds[sheetid1];
    const auto& sheet2_component_face_ids = sheets_BoundaryFaceComponentIds[sheetid2];
    std::unordered_set<size_t> sheet1_component_face_ids_set(sheet1_component_face_ids.begin(), sheet1_component_face_ids.end());
    std::unordered_set<size_t> common_component_face_ids;
    for (auto component_face_id : sheet2_component_face_ids)
        if (sheet1_component_face_ids_set.find(component_face_id) != sheet1_component_face_ids_set.end())
            common_component_face_ids.insert(component_face_id);
    return common_component_face_ids;
}

std::unordered_set<size_t> BaseComplexSheet::GetCommonComponentCellIds(size_t sheetid1, size_t sheetid2) const {
    const auto& sheet1_component_cell_ids = sheets_componentCellIds[sheetid1];
    const auto& sheet2_component_cell_ids = sheets_componentCellIds[sheetid2];
    std::unordered_set<size_t> sheet1_component_cell_ids_set(sheet1_component_cell_ids.begin(), sheet1_component_cell_ids.end());
    std::unordered_set<size_t> common_component_cell_ids;
    for (auto component_cell_id : sheet2_component_cell_ids)
        if (sheet1_component_cell_ids_set.find(component_cell_id) != sheet1_component_cell_ids_set.end())
            common_component_cell_ids.insert(component_cell_id);
    return common_component_cell_ids;
}

int find(int x, std::vector<int>& parents) {
    return parents[x] == x ? x : find(parents[x], parents);
}

int findCircleNum(std::vector<std::vector<int>>& M) {
    if (M.empty()) return 0;
    int n = M.size();

    std::vector<int> leads(n, 0);
    for (int i = 0; i < n; i++) { leads[i] = i; }   // initialize leads for every kid as themselves

    int groups = n;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {   // avoid recalculate M[i][j], M[j][i]
            if (M[i][j]) {
                int lead1 = find(i, leads);
                int lead2 = find(j, leads);
                if (lead1 != lead2) {       // if 2 group belongs 2 different leads, merge 2 group to 1
                    leads[lead1] = lead2;
                    groups--;
                }
            }
        }
    }
    return groups;
}

size_t find(int x, std::vector<std::pair<size_t, std::unordered_set<size_t>>>& parents) {
    return parents[x].first == x ? x : find(parents[x].first, parents);
}

std::vector<std::unordered_set<size_t>> findCircleGroups(std::vector<std::vector<size_t>>& M) {
    if (M.empty()) return {};
    int n = M.size();

    std::vector<std::pair<size_t, std::unordered_set<size_t>>> leads(n);
    for (int i = 0; i < n; i++) { leads[i].first = i; }   // initialize leads for every kid as themselves

    int groups = n;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {   // avoid recalculate M[i][j], M[j][i]
            if (M[i][j]) {
                size_t lead1 = find(i, leads);
                size_t lead2 = find(j, leads);
                leads[lead1].second.insert(i);
                leads[lead2].second.insert(j);
                if (lead1 != lead2) {       // if 2 group belongs 2 different leads, merge 2 group to 1
                    leads[lead1].first = lead2;
                    leads[lead2].second.insert(leads[lead1].second.begin(), leads[lead1].second.end());
                    leads[lead1].second.clear();
                    groups--;
                }
            }
        }
    }

//    std::sort(leads.begin(), leads.end(), [&](const std::pair<size_t, std::unordered_set<size_t>>& a, const std::pair<size_t, std::unordered_set<size_t>>& b){
//        return a.first < b.first;});
//    leads.resize(std::distance(leads.begin(), std::unique(leads.begin(), leads.end())));
//    std::vector<std::unordered_set<size_t>> res;
//    for (auto& e : leads)
//        res.push_back(e.second);

    std::vector<std::unordered_set<size_t>> res;
    for (auto& e : leads)
        if (!e.second.empty()) res.push_back(e.second);
    if (res.size() > 1) {
        std::cout << "findCircleGroups numberOfgroups = " << groups << " results groups = " << res.size() << std::endl;
    }
    return res;
}


//size_t BaseComplexSheet::GetNumOfIntersections(const std::unordered_set<size_t>& common_component_cell_ids) const {
//    size_t common_component_cell_ids_size = common_component_cell_ids.size();
//    auto relation_float = common_component_cell_ids.size();
//    if (common_component_cell_ids_size > 1) {
//        std::vector<size_t> common_component_cell_ids_array(common_component_cell_ids.begin(), common_component_cell_ids.end());
//        auto combs = combine(common_component_cell_ids_size, 2);
//        for (auto& comb : combs) {
//            auto component_id1 = common_component_cell_ids_array[comb[0]];
//            auto component_id2 = common_component_cell_ids_array[comb[1]];
//            if (IsAdjacent(component_id1, component_id2)) --relation_float;
//        }
//    }
//    return relation_float;
//}

size_t BaseComplexSheet::GetNumOfIntersections(const std::unordered_set<size_t>& common_component_cell_ids) const {
    //int n = baseComplex.componentC.size();
    int n = common_component_cell_ids.size();
    std::vector<std::vector<int>> M(n, std::vector<int>(n, 0));
    size_t common_component_cell_ids_size = common_component_cell_ids.size();
    std::vector<size_t> common_component_cell_ids_array(common_component_cell_ids.begin(), common_component_cell_ids.end());
    auto combs = combine(common_component_cell_ids_size, 2);
    for (auto& comb : combs) {
        auto component_id1 = common_component_cell_ids_array[comb[0]];
        auto component_id2 = common_component_cell_ids_array[comb[1]];
        if (IsAdjacent(component_id1, component_id2))
            M[comb[0]][comb[1]] = M[comb[1]][comb[0]] = 1;
    }
    return findCircleNum(M);
}

void BaseComplexSheet::ComputeComplexityDrChen(int sheetid) {

    size_t max_num_of_intersections = 1;
    float max_num_of_hybrid = 1;

    auto sorted_main_sheet_ids = sheets_coverSheetIds.at(sheetid);
    std::sort(sorted_main_sheet_ids.begin(), sorted_main_sheet_ids.end());

    size_t n = sheets_connectivities.size();
    sheets_connectivities_float.clear();
    sheets_connectivities_float.resize(n,std::vector<float>(n, 0));
    auto combinations = combine(n, 2);
    for (auto& p : combinations) {
        auto sheetid1 = sorted_main_sheet_ids[p[0]];
        auto sheetid2 = sorted_main_sheet_ids[p[1]];
        auto relation = sheets_connectivities[p[0]][p[1]];
        if (relation == ADJACENT) {
            auto common_component_face_ids = GetCommonComponentFaceIds(sheetid1, sheetid2);
            auto common_component_size = GetCommonComponentCellIds(common_component_face_ids).size();//common_component_face_ids.size() * 2;
            auto total_component_size = sheets_componentCellIds[sheetid1].size() + sheets_componentCellIds[sheetid2].size();
            auto relation_float = (float(common_component_size))/total_component_size;
            sheets_connectivities_float[p[0]][p[1]] = relation_float;
            sheets_connectivities_float[p[1]][p[0]] = relation_float;
        } else if (relation == INTERSECTING) {
            auto relation_float = GetNumOfIntersections(GetCommonComponentCellIds(sheetid1, sheetid2));
            ++relation_float;
            if (relation_float > max_num_of_intersections) max_num_of_intersections = relation_float;
            sheets_connectivities_float[p[0]][p[1]] = relation_float;
            sheets_connectivities_float[p[1]][p[0]] = relation_float;
        } else if (relation == HYBRID) {
            auto common_component_face_ids = GetCommonComponentFaceIds(sheetid1, sheetid2);
            auto common_component_size = GetCommonComponentCellIds(common_component_face_ids).size();//common_component_face_ids.size() * 2;
            auto total_component_size = sheets_componentCellIds[sheetid1].size() + sheets_componentCellIds[sheetid2].size();
            auto relation_float = (float(common_component_size)) / total_component_size;

            auto numOfIntersections = GetNumOfIntersections(GetCommonComponentCellIds(sheetid1, sheetid2));

            relation_float = (1 + relation_float) * numOfIntersections;
            if (relation_float > max_num_of_hybrid) max_num_of_hybrid = relation_float;
            sheets_connectivities_float[p[0]][p[1]] = relation_float;
            sheets_connectivities_float[p[1]][p[0]] = relation_float;
        }
    }

    auto sheets_connectivities_real = sheets_connectivities_float;
    for (auto& p : combinations) {
        auto sheetid1 = sorted_main_sheet_ids[p[0]];
        auto sheetid2 = sorted_main_sheet_ids[p[1]];
        auto relation = sheets_connectivities[p[0]][p[1]];
        if (relation == ADJACENT) {
            ;
        } else if (relation == INTERSECTING) {
            sheets_connectivities_float[p[0]][p[1]] /= max_num_of_intersections;
            sheets_connectivities_float[p[1]][p[0]] /= max_num_of_intersections;
        } else if (relation == HYBRID) {
            sheets_connectivities_float[p[0]][p[1]] /= max_num_of_hybrid;
            sheets_connectivities_float[p[1]][p[0]] /= max_num_of_hybrid;
        }
    }

    /////////////////////////////////
    std::ifstream ifs("sheet_complexity.txt");
    std::vector<size_t> sheets_complexities;
    std::string line;
    while (getline(ifs, line)) {
        std::stringstream ss(line);
        std::string str;
        while (ss >> str);
        sheets_complexities.push_back(std::stoi(str));
    }
    if (sheets_connectivities.size() == 1) {
        std::cout << "\n\n\n ********** COMPLEXITY *********** \n";
        std::cout << "#COMPLEXITY = " << sheets_complexities[sheets_coverSheetIds[0][0]];
        std::cout << "\n\n\n ********** COMPLEXITY *********** \n";
        sheets_connectivities[0][0] = sheets_complexities[sheets_coverSheetIds[0][0]];
    }
//    float total_complexity = 0.0f;
//    for (int i = 0; i < all_sheets_connectivities.size(); ++i)
//        total_complexity += sheets_complexities[i] * sheets_importance[i];
    ////////////////////////////////////
//    int num = 0;
//    for (auto& row : sheets_connectivities_float)
//        for (auto ele : row) num += (ele != 0);

    int id = 0;
    float max_sum = 0;
    for (auto& row : sheets_connectivities_real) {
        float sum = 0.0f;
        for (auto ele: row) sum += fabs(ele);
        sum *= sheets_complexities[sorted_main_sheet_ids[id]];
        if (sum > max_sum) max_sum = sum;
        row[id++] = sum;
    }
    id = 0;
    for (auto& row : sheets_connectivities_real) {
        sheets_connectivities_float[id][id] = sheets_connectivities_real[id][id];
        //row[id++] /= max_sum;
        sheets_connectivities_float[id][id] /= max_sum;
        ++id;
    }
    struct sc {
        int s; // sheet_id
        std::vector<float> c; // connectivities
        sc(){}
        sc(const sc& o) : s(o.s), c(o.c){}
        sc& operator = (const sc& o) {
            s = o.s;
            c = o.c;
            return *this;
        }
    };
    std::vector<sc> ss(n);
    for (int i = 0; i < n; ++i) {
        ss[i].c = sheets_connectivities_real[i];
        ss[i].s = i;
    }
    std::sort(ss.begin(), ss.end(), [&](const sc& a, const sc& b){return a.c[a.s] > b.c[b.s];});

    std::cout << "\n\n\n ********** NEW *********** \n";
    for (auto& row : ss)
        std::cout << sorted_main_sheet_ids[row.s] << " ";
    std::cout << "\n *************************************** \n\n\n";
    {
        int i = 0;
        for (auto& row : ss)
            sheets_coverSheetIds.at(sheetid).at(i++) = sorted_main_sheet_ids[row.s];
    }

//    auto sheets_connectivities_float_new = sheets_connectivities_float;
//    auto sheets_connectivities_new = sheets_connectivities;
//    for (int i = 0; i < n; ++i)
//        for (int j = 0; j < n; ++j) {
//            sheets_connectivities_float_new[i][j] = sheets_connectivities_float[ss[i].s][ss[j].s];
//            sheets_connectivities_new[i][j] = sheets_connectivities[ss[i].s][ss[j].s];
//        }
//
//    sheets_connectivities_float = sheets_connectivities_float_new;
//    sheets_connectivities = sheets_connectivities_new;

    /////////////////////////////////
    {
        std::ofstream ofs("complexity.mat");
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j)
                ofs << sheets_connectivities_real[i][j] << "\t";
            ofs << "\n";
        }
    }
    {
        std::ofstream ofs("diagonal.mat");
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j)
                ofs << (i == j ? sheets_connectivities_float[i][j] : 0) << "\t";
            ofs << "\n";
        }
    }
    {
        std::ofstream ofs("adjacent.mat");
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j)
                ofs << (sheets_connectivities[i][j] == ADJACENT ? sheets_connectivities_float[i][j] : 0) << "\t";
            ofs << "\n";
        }
    }
    {
        std::ofstream ofs("intersecting.mat");
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j)
                ofs << (sheets_connectivities[i][j] == INTERSECTING ? sheets_connectivities_float[i][j] : 0) << "\t";
            ofs << "\n";
        }
    }
    {
        std::ofstream ofs("hybrid.mat");
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j)
                ofs << (sheets_connectivities[i][j] == HYBRID ? sheets_connectivities_float[i][j] : 0) << "\t";
            ofs << "\n";
        }
    }

    {
        Eigen::MatrixXf m(n,n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                m(i, j) = sheets_connectivities_real[i][j];
        //if (sheetid == 0)
        {
            std::cout << "\n\n\n ********** COMPLEXITY *********** \n";
            std::cout << "#mainsheets = " << sheets_coverSheetIds.at(sheetid).size() << " main sheets set = " << sheetid << " #COMPLEXITY = " << m.norm();
            std::cout << "\n\n\n ********** COMPLEXITY *********** \n";
        }
    }
}

void BaseComplexSheet::ComputeComplexityUnbalancedMatrix(int sheetid) {

    size_t max_num_of_intersections = 1;
    float max_num_of_hybrid = 1;

    auto sorted_main_sheet_ids = sheets_coverSheetIds.at(sheetid);
    std::sort(sorted_main_sheet_ids.begin(), sorted_main_sheet_ids.end());

    size_t n = sheets_connectivities.size();
    sheets_connectivities_float.clear();
    sheets_connectivities_float.resize(n,std::vector<float>(n, 0));
    auto combinations = combine(n, 2);
    for (auto& p : combinations) {
        auto sheetid1 = sorted_main_sheet_ids[p[0]];
        auto sheetid2 = sorted_main_sheet_ids[p[1]];
        auto relation = sheets_connectivities[p[0]][p[1]];
        if (relation == ADJACENT) {
            auto common_component_face_ids = GetCommonComponentFaceIds(sheetid1, sheetid2);
            auto common_component_size = GetCommonComponentCellIds(common_component_face_ids).size();//common_component_face_ids.size() * 2;
            auto total_component_size = sheets_componentCellIds[sheetid1].size() + sheets_componentCellIds[sheetid2].size();
            auto relation_float = (float(common_component_size))/total_component_size;
            sheets_connectivities_float[p[0]][p[1]] = relation_float;
            sheets_connectivities_float[p[1]][p[0]] = relation_float;
        } else if (relation == INTERSECTING) {
            auto num_of_intersections = GetNumOfIntersections(GetCommonComponentCellIds(sheetid1, sheetid2));
            if (num_of_intersections > max_num_of_intersections) max_num_of_intersections = num_of_intersections;
            sheets_connectivities_float[p[0]][p[1]] = num_of_intersections;
            sheets_connectivities_float[p[1]][p[0]] = num_of_intersections;
        } else if (relation == HYBRID) {
            auto common_component_face_ids = GetCommonComponentFaceIds(sheetid1, sheetid2);
            auto common_component_size = GetCommonComponentCellIds(common_component_face_ids).size();//common_component_face_ids.size() * 2;
            auto total_component_size = sheets_componentCellIds[sheetid1].size() + sheets_componentCellIds[sheetid2].size();
            auto relation_float = (float(common_component_size)) / total_component_size;

            auto numOfIntersections = GetNumOfIntersections(GetCommonComponentCellIds(sheetid1, sheetid2));

            relation_float = (1 + relation_float) * numOfIntersections;
            if (relation_float > max_num_of_hybrid) max_num_of_hybrid = relation_float;
            sheets_connectivities_float[p[0]][p[1]] = relation_float;
            sheets_connectivities_float[p[1]][p[0]] = relation_float;
        }
    }

    auto sheets_connectivities_real = sheets_connectivities_float;
    for (auto& p : combinations) {
        auto sheetid1 = sorted_main_sheet_ids[p[0]];
        auto sheetid2 = sorted_main_sheet_ids[p[1]];
        auto relation = sheets_connectivities[p[0]][p[1]];
        if (relation == ADJACENT) {
            ;
        } else if (relation == INTERSECTING) {
            sheets_connectivities_float[p[0]][p[1]] /= max_num_of_intersections;
            sheets_connectivities_float[p[1]][p[0]] /= max_num_of_intersections;
        } else if (relation == HYBRID) {
            sheets_connectivities_float[p[0]][p[1]] /= max_num_of_hybrid;
            sheets_connectivities_float[p[1]][p[0]] /= max_num_of_hybrid;
        }
    }

    /////////////////////////////////
//    std::ifstream ifs("sheet_complexity.txt");
//    //std::ifstream ifs("min_main_chords_complexity.txt");
//    std::vector<size_t> sheets_complexities;
//    std::string line;
//    while (getline(ifs, line)) {
//        std::stringstream ss(line);
//        std::string str;
//        while (ss >> str);
//        sheets_complexities.push_back(std::stof(str));
//    }
//
//    for (int i = 1; i < n; ++i)
//        for (int j = 0; j < i; ++j)
//            sheets_connectivities_real[i][j] *= sheets_complexities[sorted_main_sheet_ids[j]];
//    for (int i = 0; i < n; ++i)
//        for (int j = i + 1; j < n; ++j)
//            sheets_connectivities_real[i][j] *= sheets_complexities[sorted_main_sheet_ids[j]];
    int id = 0;
    float max_sum = 0;
    for (auto& row : sheets_connectivities_real) {
        float sum = 0.0f;
        for (auto ele: row) sum += fabs(ele);
        // sum *= sheets_complexities[sorted_main_sheet_ids[id]];
        if (sum > max_sum) max_sum = sum;
        row[id++] = sum;
    }
    id = 0;
    for (auto& row : sheets_connectivities_real) {
        sheets_connectivities_float[id][id] = sheets_connectivities_real[id][id];
        sheets_connectivities_float[id][id] /= 1.25 * max_sum;
        ++id;
    }
    struct sc {
        int s; // sheet_id
        std::vector<float> c; // connectivities
        sc(){}
        sc(const sc& o) : s(o.s), c(o.c){}
        sc& operator = (const sc& o) {
            s = o.s;
            c = o.c;
            return *this;
        }
    };
    std::vector<sc> ss(n);
    for (int i = 0; i < n; ++i) {
        ss[i].c = sheets_connectivities_real[i];
        ss[i].s = i;
    }
    std::sort(ss.begin(), ss.end(), [&](const sc& a, const sc& b){return a.c[a.s] > b.c[b.s];});


    auto sheets_connectivities_float_new = sheets_connectivities_float;
    auto sheets_connectivities_real_new = sheets_connectivities_real;
    auto sheets_connectivities_new = sheets_connectivities;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            sheets_connectivities_float_new[i][j] = sheets_connectivities_float[ss[i].s][ss[j].s];
            sheets_connectivities_real_new[i][j] = sheets_connectivities_real[ss[i].s][ss[j].s];
            sheets_connectivities_new[i][j] = sheets_connectivities[ss[i].s][ss[j].s];
        }

    sheets_connectivities_float = sheets_connectivities_float_new;
    sheets_connectivities_real = sheets_connectivities_real_new;
    sheets_connectivities = sheets_connectivities_new;


    // std::cout << "\n\n\n ********** NEW *********** \n";
    std::cout << "---------------------------- \n";
    for (auto& row : ss)
        std::cout << sorted_main_sheet_ids[row.s] << " ";
    std::cout << "\n";
    int xx = 0;
    for (auto& row : ss) {
        auto num_of_singularities = GetNumberOfSingularities(sorted_main_sheet_ids[row.s]);
        std::cout << num_of_singularities << " ";
        //sheets_connectivities_real[xx++][n] = num_of_singularities;
    }
    std::cout << "\n";
    {
        int i = 0;
        for (auto& row : ss)
            sheets_coverSheetIds.at(sheetid).at(i++) = sorted_main_sheet_ids[row.s];
    }
    WriteComplexityMat(sheets_connectivities_real, ("complexity" + std::to_string(sheetid) + ".mat").c_str());
    WriteDiagonalMat(sheets_connectivities_float, ("diagonal" + std::to_string(sheetid) + ".mat").c_str());
    WriteAdjacentMat(sheets_connectivities_float, ("adjacent" + std::to_string(sheetid) + ".mat").c_str());
    WriteIntersectingMat(sheets_connectivities_float, ("intersecting" + std::to_string(sheetid) + ".mat").c_str());
    WriteHybridMat(sheets_connectivities_float, ("hybrid" + std::to_string(sheetid) + ".mat").c_str());

    {
        Eigen::MatrixXf m(n,n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                m(i, j) = (i == j ? 0 : sheets_connectivities_real[i][j]);
                //m(i, j) = sheets_connectivities_real[i][j];
        std::cout << "#mainsheets = " << sheets_coverSheetIds.at(sheetid).size() << " main sheets set = " << sheetid << " #COMPLEXITY = " << m.norm()
                << " Average = " << m.norm() / sheets_coverSheetIds.at(sheetid).size() << "\n";

		sheetIds_overlaps_complexity soc;
		size_t overlaps = 0;
		for (auto& row : ss) {
			auto sheetid = sorted_main_sheet_ids[row.s];
			soc.sheetIds.insert(sheetid);
			overlaps += sheets_componentCellIds.at(sheetid).size();
		}
		soc.overlaps = overlaps - baseComplex.componentC.size();
		soc.complexity = m.norm() / sheets_coverSheetIds.at(sheetid).size();
		socs.push_back(soc);
	}
}

void BaseComplexSheet::WriteComplexityMat(const std::vector<std::vector<float>>& M, const char* filename) const {
    std::ofstream ofs(filename);
    for (int i = 0; i < M.size(); ++i) {
        for (int j = 0; j < M[0].size(); ++j)
            ofs << M[i][j] << "\t";
        ofs << "\n";
    }
}
void BaseComplexSheet::WriteDiagonalMat(const std::vector<std::vector<float>>& M, const char* filename) const {
    std::ofstream ofs(filename);
    for (int i = 0; i < M.size(); ++i) {
        for (int j = 0; j < M[0].size(); ++j)
            ofs << (i == j ? M[i][j] : 0) << "\t";
        ofs << "\n";
    }
}
void BaseComplexSheet::WriteAdjacentMat(const std::vector<std::vector<float>>& M, const char* filename) const {
    std::ofstream ofs(filename);
    for (int i = 0; i < M.size(); ++i) {
        for (int j = 0; j < M[0].size(); ++j)
            ofs << (sheets_connectivities[i][j] == ADJACENT ? M[i][j] : 0) << "\t";
        ofs << "\n";
    }
}
void BaseComplexSheet::WriteIntersectingMat(const std::vector<std::vector<float>>& M, const char* filename) const {
    std::ofstream ofs(filename);
    for (int i = 0; i < M.size(); ++i) {
        for (int j = 0; j < M[0].size(); ++j)
            ofs << (sheets_connectivities[i][j] == INTERSECTING ? M[i][j] : 0) << "\t";
        ofs << "\n";
    }
}
void BaseComplexSheet::WriteHybridMat(const std::vector<std::vector<float>>& M, const char* filename) const {
    std::ofstream ofs(filename);
    for (int i = 0; i < M.size(); ++i) {
        for (int j = 0; j < M[0].size(); ++j)
            ofs << (sheets_connectivities[i][j] == HYBRID ? M[i][j] : 0) << "\t";
        ofs << "\n";
    }
}

int BaseComplexSheet::GetNumberOfSingularities(int sheet_id) const {
    std::unordered_set<size_t> singularity_ids;
    for (auto component_id : sheets_componentCellIds.at(sheet_id)) {
        for (auto component_edge_id : baseComplex.componentC.at(component_id).Eids) {
            auto& component_edge = baseComplex.componentE.at(component_edge_id);
            auto& edge = baseComplex.mesh.E.at(component_edge.eids_link.front());
            if (edge.isSingularity)
                singularity_ids.insert(edge.singularEid);
        }
    }
    return singularity_ids.size();
}

//void dfs(const float alpha = 1.0f, const std::vector<std::vector<size_t>>& w, std::unordered_set<size_t>& current_sheet_set, std::vector<bool>& HB) {
//    ;
//}

void BaseComplexSheet::ExtractMainSheets() {
    auto sheet_intersecting_component_ids_groups = Get_sheet_intersecting_component_ids_groups();

}


std::vector<std::unordered_set<size_t>> BaseComplexSheet::GetIntersectionGroups(const std::unordered_set<size_t>& common_component_cell_ids) const {
    //int n = baseComplex.componentC.size();
    int n = common_component_cell_ids.size();
    std::vector<std::vector<size_t>> M(n, std::vector<size_t>(n, 0));
    size_t common_component_cell_ids_size = common_component_cell_ids.size();
    std::vector<size_t> common_component_cell_ids_array(common_component_cell_ids.begin(), common_component_cell_ids.end());
    auto combs = combine(common_component_cell_ids_size, 2);
    for (auto& comb : combs) {
        auto component_id1 = common_component_cell_ids_array[comb[0]];
        auto component_id2 = common_component_cell_ids_array[comb[1]];
        if (IsAdjacent(component_id1, component_id2))
            M[comb[0]][comb[1]] = M[comb[1]][comb[0]] = 1;
    }

    auto groups = findCircleGroups(M);
    std::vector<std::unordered_set<size_t>> res(groups.size());
    int id = 0;
    for (auto& group : groups) {
        for (auto& e : group)
            res[id].insert(common_component_cell_ids_array[e]);
        ++id;
    }
    return res;
}

std::vector<std::vector<std::vector<std::unordered_set<size_t>>>> BaseComplexSheet::Get_sheet_intersecting_component_ids_groups() const {
    // sheet_intersecting_component_ids_groups[i][j][k]
    // i, j denotes sheet i and sheet j
    // k denotes the intersection group id
    std::vector<std::vector<std::vector<std::unordered_set<size_t>>>> sheet_intersecting_component_ids_groups;
    int n = all_sheets_connectivities.size();
    sheet_intersecting_component_ids_groups.resize(n, std::vector<std::vector<std::unordered_set<size_t>>>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            if (all_sheets_connectivities[i][j] == INTERSECTING) {
                // int num_of_intersection = GetNumOfIntersections(GetCommonComponentCellIds(i, j));
                auto groups = GetIntersectionGroups(GetCommonComponentCellIds(i, j));
                sheet_intersecting_component_ids_groups[i][j] = sheet_intersecting_component_ids_groups[j][i] = groups;
            }
        }
    }
    return sheet_intersecting_component_ids_groups;
}

void BaseComplexSheet::VerifySheetDecompositions() {
    for (const auto& coverSheetIds : sheets_coverSheetIds) {
        std::vector<size_t> componentCellId_count(baseComplex.componentC.size(), 0);
        for (auto sheetId : coverSheetIds)
            for (auto componentCellId : sheets_componentCellIds.at(sheetId))
                ++componentCellId_count.at(componentCellId);
        bool correct = true;
        for (auto count : componentCellId_count)
            if (count == 0) {
                correct = false;
                break;
            }
        if (!correct) {
            std::cerr << red_color << "\n---- SheetDecomposition is incorrect! ----\n"<< end_color;
            for (auto sheetId : coverSheetIds)
                std::cerr << " " << sheetId;
            std::cerr << "\n";
        }
    }
//	std::vector<size_t> xx;
//	for (int i = 0; i < sheets_componentCellIds.size(); ++i) xx.push_back(i);
//	sheets_coverSheetIds.push_back(xx);
}

bool BaseComplexSheet::Verify() {
	std::ofstream ofs("verification.txt");
	ofs << "---- sheets_componentCellIds\n";
//	for (auto& componentCellIds : sheets_componentCellIds) {
//		for (auto componentCellId : componentCellIds)
//			ofs << " " << componentCellId;
//		ofs << "\n";
//	}
	ofs << "---- all combinations\n";
	auto n = sheets_componentCellIds.size(); // number of sheets;
	std::vector<std::vector<size_t>> combinations;
	for (int i = 1; i <= n; ++i) {
		for (int j = 1; j <= i; ++j) {
			auto combs = combine(i, j);
//			for (auto& comb : combs) {
//				for (auto c : comb)
//					ofs << " " << c;
//				ofs << "\n";
//			}
			std::copy(combs.begin(), combs.end(), std::back_inserter(combinations));
		}
	}
	ofs << "---- all valid combinations\n";
	std::vector<std::vector<size_t>> nc;
	for (auto& combs : combinations) {
		std::vector<bool> covered(baseComplex.componentC.size(), false);
		for (auto c : combs) {
			for (auto componentCellId: sheets_componentCellIds.at(c))
				covered[componentCellId] = true;
		}
		bool coverAll = true;
		for (auto c : covered)
			if (!c) {
				coverAll = false;
				break;
			}
		if (coverAll) {
			nc.push_back(combs);
//			for (auto c : combs)
//				ofs << " " << c;
//			ofs << "\n";
		}
	}
	combinations = nc;
	auto old = sheets_coverSheetIds;
	sheets_coverSheetIds = combinations;
	RemoveSheetSetsRedundancy();
	combinations = sheets_coverSheetIds;
	ofs << "---- all combinations without redundancy\n";
	for (auto& combs : combinations) {
		for (auto c : combs)
			ofs << " " << c;
		ofs << "\n";
	}
	std::vector<size_t> xx;
	for (int i = 0; i < n; ++i) xx.push_back(i);
	sheets_coverSheetIds.push_back(xx);
	//sheets_coverSheetIds = old;
	return true;
}
