/*
 * BaseComplexSheet.cpp
 *
 *  Created on: Oct 5, 2017
 *      Author: cotrik
 */

#include "BaseComplexSheetQuad.h"
#include "MeshFileWriter.h"
#include "Dual.h"
#include "RefinedDualQuad.h"
#include <algorithm>
#include <sstream>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <iterator>
#include <Eigen/Dense>

static std::set<std::set<size_t>> shortestCombination(const std::set<size_t>& filter, const std::vector<std::set<size_t>>& listOfSets) {
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
static void FindSetCombination(std::vector<std::set<size_t>>& input, std::set<size_t>& target, std::vector<std::set<size_t>>& output) {
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

static std::vector<std::vector<std::set<size_t>>> FindSetCombination(std::vector<std::set<size_t>>& input, std::set<size_t>& target) {

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

BaseComplexSheetQuad::BaseComplexSheetQuad(BaseComplexQuad& baseComplex)
: baseComplex(baseComplex)
{
    // TODO Auto-generated constructor stub

}

BaseComplexSheetQuad::~BaseComplexSheetQuad()
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
void BaseComplexSheetQuad::Extract()
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

void BaseComplexSheetQuad::ExtractSheetDecompositions(const bool bfs) {
    faceComponent_sheetIds.resize(baseComplex.componentF.size());
    for (size_t sheetId = 0; sheetId < sheets_componentFaceIds.size(); ++sheetId) {
        for (auto component : sheets_componentFaceIds[sheetId])
            if (faceComponent_sheetIds[component].find(sheetId) == faceComponent_sheetIds[component].end())
                faceComponent_sheetIds[component].insert(sheetId);
    }
    for (size_t sheetId = 0; sheetId < sheets_componentFaceIds.size(); ++sheetId) {
        auto sheetBoundaryEdgeComponentIds = GetSheetBoundaryEdgeComponentIds(sheetId);
        for (auto edgeComponentId : sheetBoundaryEdgeComponentIds)
            edgeComponent_neighborSheetIds[edgeComponentId].insert(sheetId);
        sheets_BoundaryEdgeComponentIds.push_back(sheetBoundaryEdgeComponentIds);
    }
    for (size_t sheetId = 0; sheetId < sheets_componentCellIds.size(); ++sheetId)
        sheets_coverSheetIds.push_back(!bfs ? GetCoverSheetIds(sheetId) : GetCoverSheetIdsBFS(sheetId));

    std::sort(sheets_coverSheetIds.begin(), sheets_coverSheetIds.end(), [&](const std::vector<size_t>& a, const std::vector<size_t>& b) {return a.size() < b.size();});
    std::cout << "****** SheetDecompositions ******\n";
    for (auto& sheetIds : sheets_coverSheetIds) {
        for (auto sheetId : sheetIds) std::cout << sheetId << " ";
        std::cout << "\n";
    }
    RemoveSheetSetsRedundancy();
    std::sort(sheets_coverSheetIds.begin(), sheets_coverSheetIds.end(), [&](const std::vector<size_t>& a, const std::vector<size_t>& b) {return a.size() < b.size();});
    std::cout << "****** SheetDecompositions after RemoveSheetSetsRedundancy******\n";
    for (auto& sheetIds : sheets_coverSheetIds) {
        for (auto sheetId : sheetIds) std::cout << sheetId << " ";
        std::cout << "\n";
    }
}

void BaseComplexSheetQuad::WriteSheetDecompositionsFile(const char *filename) const {
    std::ofstream ofs("sheet_decompositions.txt");
    for (auto& sheetIds : sheets_coverSheetIds) {
        for (auto sheetId : sheetIds) ofs << sheetId << " ";
        ofs << "\n";
    }
    WriteSheetNeighborsSheetIdsJS("sheet_neighborsheetids.js");
}

void BaseComplexSheetQuad::WriteSheetDecompositionsDuaVTK(const char *filename, const int scalar) const {
    if (sheets_coverSheetIds.empty()) return;
    RefinedDualQuad dual(baseComplex.mesh);
    dual.Build();
    std::vector<std::unordered_set<size_t>> dualEdgeIdsSet;
    for (auto sheet_id : sheets_coverSheetIds[0])
        dualEdgeIdsSet.push_back(GetDualEdgeIds(sheet_id));

    const auto& mesh = baseComplex.GetMesh();
    const auto& V = dual.V;
    const auto& E = dual.E;
    std::ofstream ofs(filename);
    ofs << "# vtk DataFile Version 2.0\n"
        << filename << "\n"
        << "ASCII\n\n"
        << "DATASET POLYDATA\n";
    ofs << "POINTS " << V.size() << " double" << "\n";
    for (const auto& v : V)
        ofs << v.x << " " << v.y << " " << v.z << "\n";

    size_t numOfDualEdgeIds = 0;
    for (auto& dualEdgeIds : dualEdgeIdsSet)
        numOfDualEdgeIds += dualEdgeIds.size();
    ofs << "Lines " << numOfDualEdgeIds << " " << 3 * numOfDualEdgeIds << "\n";
    for (auto& dualEdgeIds : dualEdgeIdsSet)
        for (const auto edge_id : dualEdgeIds) {
            const auto& edge = E.at(edge_id);
            ofs << edge.Vids.size() << " " << edge.Vids[0] << " " << edge.Vids[1] << "\n";
        }

    ofs << "CELL_DATA " << numOfDualEdgeIds << "\n"
        << "SCALARS " << "sheet" << " int 1\n"
        << "LOOKUP_TABLE default\n";
    for (auto& dualEdgeIds : dualEdgeIdsSet)
        for (const auto id : dualEdgeIds)
            ofs << scalar << "\n";
}

void BaseComplexSheetQuad::RemoveSheetSetsRedundancy() {
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

void BaseComplexSheetQuad::RemoveSheetSetsRedundancy(std::vector<size_t>& coverSheetIds) {
    auto sheetIds = coverSheetIds;
    std::sort(sheetIds.begin(), sheetIds.end(), [&](const size_t& a, const size_t& b) {
        return sheets_componentFaceIds[a].size() < sheets_componentFaceIds[b].size();});
    std::vector<size_t> componentFaceId_count(baseComplex.componentF.size(), 0);
    for (auto sheetId : sheetIds)
        for (auto componentFaceId : sheets_componentFaceIds[sheetId])
            ++componentFaceId_count[componentFaceId];
    while (true) {
        bool sheetIdsChanged = false;
        for (auto iter = sheetIds.begin(); iter != sheetIds.end();) {
            auto sheetId = *iter;
            bool canRemove = true;
            for (auto componentFaceId : sheets_componentFaceIds[sheetId])
                if (componentFaceId_count[componentFaceId] <= 1) {
                    canRemove = false;
                    break;
                }
            if (canRemove) {
                for (auto componentFaceId : sheets_componentFaceIds[sheetId])
                    --componentFaceId_count[componentFaceId];
                iter = sheetIds.erase(iter);
                sheetIdsChanged = false;
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

std::unordered_set<size_t> BaseComplexSheetQuad::GetSheetBoundaryFaceComponentIds(size_t sheetId) const {
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

std::unordered_set<size_t> BaseComplexSheetQuad::GetSheetBoundaryEdgeComponentIds(size_t sheetId) const {
    std::unordered_set<size_t> sheetBoundaryEdgeComponentIds;
    std::unordered_set<size_t> sheetFaceComponentIds(sheets_componentFaceIds[sheetId].begin(), sheets_componentFaceIds[sheetId].end());
    std::unordered_set<size_t> sheetEdgeComponentIds;
    for (auto componentId : sheetFaceComponentIds) {
        const auto& component = baseComplex.componentF.at(componentId);
        sheetEdgeComponentIds.insert(component.Eids.begin(), component.Eids.end());
    }
    for (auto edgeComponentId : sheetEdgeComponentIds) {
        const auto& edgeComponent = baseComplex.componentE.at(edgeComponentId);
        if (edgeComponent.isBoundary) {
            sheetBoundaryEdgeComponentIds.insert(edgeComponentId);
            continue;
        }
        for (auto neighborFaceComponentId : edgeComponent.N_Fids)
            if (sheetFaceComponentIds.find(neighborFaceComponentId) == sheetFaceComponentIds.end())
                sheetBoundaryEdgeComponentIds.insert(edgeComponentId);
    }
    return sheetBoundaryEdgeComponentIds;
}

static bool isSheetIdExistedInResult(size_t sheetId, const std::vector<size_t>& res) {
    bool existed = false;
    for (auto id : res)
        if (id == sheetId) {
            existed = true;
            break;
        }
    return existed;
}

static bool hasCoveredAllComponents(size_t sheetId, const std::vector<size_t>& componentCovered, const BaseComplex& baseComplex, bool& coveredAllComponent) {
    for (size_t componentId = 0; componentId < baseComplex.componentF.size(); ++componentId)
        if (!componentCovered[componentId]) {
            coveredAllComponent = false;
            break;
        }
    return coveredAllComponent;
}

bool BaseComplexSheetQuad::HasCoveredSheetComponents(size_t sheetId, const std::vector<size_t>& componentCovered) const {
    bool coveredSheetComponent = true;
    for (auto component : sheets_componentCellIds[sheetId])
        if (!componentCovered[component]) {
            coveredSheetComponent = false;
            break;
        }
    return coveredSheetComponent;
}

std::vector<size_t> BaseComplexSheetQuad::GetCoverSheetIds(size_t beginSheetId) const {
    std::vector<size_t> res;
    std::unordered_set<size_t> resSet;
    std::queue<size_t> q;
    q.push(beginSheetId);
    std::vector<size_t> componentCovered(baseComplex.componentF.size(), false);
    while (!q.empty()) {
        bool coveredAllComponent = true;
        size_t n = q.size();
        for (size_t i = 0; i < n; ++i) {
            auto sheetId = q.front();
            q.pop();
            if (isSheetIdExistedInResult(sheetId, res)) continue;
            res.push_back(sheetId);
            resSet.insert(sheetId);

            for (auto component : sheets_componentFaceIds[sheetId])
                componentCovered[component] = true;
            if (hasCoveredAllComponents(sheetId, componentCovered, baseComplex, coveredAllComponent)) break;

            std::unordered_set<size_t> neighborSheetIds = GetNeighborSheetIds(sheetId);
            q.push(GetMinComponentIntersectionNeighborSheetId(sheetId, neighborSheetIds, resSet));
//            auto candidateSheetIds = GetMinComponentIntersectionNeighborSheetIds(sheetId, neighborSheetIds, resSet);
//            for (auto candidateSheetId : candidateSheetIds)
//                q.push(candidateSheetId);
        }
        if (coveredAllComponent) break;
    }
    return res;
}

std::vector<size_t> BaseComplexSheetQuad::GetCoverSheetIdsBFS(size_t beginSheetId) const {
    std::vector<size_t> res;
    std::unordered_set<size_t> resSet;
    std::queue<size_t> q;
    q.push(beginSheetId);
    std::vector<size_t> componentCovered(baseComplex.componentC.size(), false);
    while (!q.empty()) {
        bool coveredAllComponent = true;
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
            for (auto candidateSheetId : candidateSheetIds)
                if (!isSheetIdExistedInResult(candidateSheetId, res))
                    q.push(candidateSheetId);
        }
        if (coveredAllComponent) break;
    }
    return res;
}

std::unordered_set<size_t> BaseComplexSheetQuad::GetNeighborSheetIds(size_t sheetId) const {
    std::unordered_set<size_t> neighborSheetIds;
    for (auto& boundaryEdgeComponentId : sheets_BoundaryEdgeComponentIds[sheetId]) {
        auto& boundaryEdgeComponent = baseComplex.componentE.at(boundaryEdgeComponentId);
        for (auto neiborghComponentId : boundaryEdgeComponent.N_Fids) {
            for (auto neighborSheetId : faceComponent_sheetIds[neiborghComponentId]) {
                if (neighborSheetId == sheetId) continue;
                neighborSheetIds.insert(neighborSheetId);
            }
        }
    }
    return neighborSheetIds;
}

size_t BaseComplexSheetQuad::GetMinComponentIntersectionNeighborSheetId(size_t sheetId,
        const std::unordered_set<size_t>& neighborSheetIds, const std::unordered_set<size_t>& resSet) const {
    size_t minComponentIntersectionNums = MAXID;
    size_t minComponentIntersectionNeighborSheetId = *neighborSheetIds.begin();
    std::unordered_set<size_t> sheetComponentIds(sheets_componentFaceIds[sheetId].begin(), sheets_componentFaceIds[sheetId].end());
    for (auto neighborSheetId : neighborSheetIds) {
        if (resSet.find(neighborSheetId) != resSet.end()) continue;
        std::unordered_set<size_t> neighborSheetComponentIds(sheets_componentFaceIds[neighborSheetId].begin(), sheets_componentFaceIds[neighborSheetId].end());
        size_t intersectionNum = 0;
        for (auto id : sheetComponentIds)
            if (neighborSheetComponentIds.find(id) != neighborSheetComponentIds.end()) ++intersectionNum;
        if (intersectionNum < minComponentIntersectionNums ||
                (intersectionNum == minComponentIntersectionNums &&
                        sheets_componentFaceIds[neighborSheetId].size() > sheets_componentFaceIds[minComponentIntersectionNums].size())) {
            minComponentIntersectionNums = intersectionNum;
            minComponentIntersectionNeighborSheetId = neighborSheetId;
        }
    }
    return minComponentIntersectionNeighborSheetId;
}

std::vector<size_t> BaseComplexSheetQuad::GetMinComponentIntersectionNeighborSheetIds(size_t sheetId,
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

void BaseComplexSheetQuad::WriteSheetNeighborsSheetIdsJS(const char* filename) const {
    std::ofstream ofs(filename);
    ofs << "var data = [";
    for (size_t sheetId = 0; sheetId < sheets_componentFaceIds.size(); ++sheetId) {
        //if (sheetId != 0) ofs << ",";
        std::unordered_set<size_t> neighborSheetIds = GetNeighborSheetIds(sheetId);
        std::unordered_set<size_t> sheetComponentIds(sheets_componentFaceIds[sheetId].begin(), sheets_componentFaceIds[sheetId].end());
        for (auto neighborSheetId : neighborSheetIds) {
            std::unordered_set<size_t> neighborSheetComponentIds(sheets_componentFaceIds[neighborSheetId].begin(), sheets_componentFaceIds[neighborSheetId].end());
            size_t intersectionNum = 0;
            for (auto id : sheetComponentIds)
                if (neighborSheetComponentIds.find(id) != neighborSheetComponentIds.end()) ++intersectionNum;
            ofs << ",[\'s" << sheetId << "\',\'s" << neighborSheetId << "\'," << intersectionNum << "]\n";
        }
    }
    ofs << "];";
}

void BaseComplexSheetQuad::ExtractSets()
{
    std::vector<std::set<size_t>> input;
    std::set<size_t> target;
    for (auto& i : sheets_componentFaceIds)
        input.push_back(std::set<size_t>(i.begin(), i.end()));
    for (auto i = 0; i < baseComplex.componentF.size(); target.insert(i++));
    auto output_set = FindSetCombination(input, target);
    for (auto& output : output_set) {
        std::set<size_t> ids;
        for (auto i : output)
            for (auto j = 0; j < input.size(); ++j)
                if (i == input[j]) {
                    ids.insert(j);
                    break;
                }
        std::cout << "output size " << output.size() << " *************" << std::endl;
        for (auto i : ids)
            std::cout << i << " ";
        std::cout << std::endl;
        std::cout << "***************************" << std::endl;
    }
}

void BaseComplexSheetQuad::GetParallelComponents(const ComponentEdge & componentEdge,
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
std::vector<size_t> BaseComplexSheetQuad::GetParallelComponentEdgeIds(const ComponentEdge & componentEdge) {
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

void BaseComplexSheetQuad::WriteSheetsEdgesVTK(const char *filename) const
{

}

void BaseComplexSheetQuad::WriteSheetsFacesVTK(const char *filename) const
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

    size_t faces_num = 0;
    for (const auto& sheet : sheets_componentFaceIds)
        for (const auto componentid : sheet)
            faces_num += baseComplex.componentF.at(componentid).fids_patch.size();

    ofs << "CELLS " << faces_num << " " << 5 * faces_num << "\n";
    for (const auto& sheet : sheets_componentFaceIds)
        for (const auto componentid : sheet)
            for (const auto cell_id : baseComplex.componentF.at(componentid).fids_patch) {
                    const Face& cell = mesh.F.at(cell_id);
                    ofs << cell.Vids.size();
                    for (auto vid : cell.Vids)
                        ofs << " " << vid;
                    ofs << "\n";
            }

    ofs << "CELL_TYPES " << faces_num << "\n";
    for (size_t i = 0; i < faces_num; i++)
        ofs << "9\n";

    ofs << "CELL_DATA " << faces_num << "\n"
        << "SCALARS " << "sheet" << " int 1\n"
        << "LOOKUP_TABLE default\n";
    size_t sheet_id = 0;
    for (const auto& sheet : sheets_componentFaceIds) {
        for (const auto componentid : sheet)
            for (const auto face_id : baseComplex.componentF.at(componentid).fids_patch)
                ofs << sheet_id << "\n";
        ++sheet_id;
    }
}

//void BaseComplexSheetQuad::WriteSheetsCellsVTK(const char *filename) const
//{
//    const auto& mesh = baseComplex.GetMesh();
//    const auto& V = mesh.V;
//    const auto& E = mesh.E;
//    std::ofstream ofs(filename);
//    ofs << "# vtk DataFile Version 2.0\n"
//        << filename << "\n"
//        << "ASCII\n\n"
//        << "DATASET UNSTRUCTURED_GRID\n";
//    ofs << "POINTS " << V.size() << " float" << "\n";
//    for (size_t i = 0; i < V.size(); i++)
//        ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << "\n";
//
//    size_t cells_num = 0;
//    for (const auto& sheet : sheets_componentCellIds)
//        for (const auto componentid : sheet)
//            cells_num += baseComplex.componentC.at(componentid).cids_patch.size();
//
//    ofs << "CELLS " << cells_num << " " << 9 * cells_num << "\n";
//    for (const auto& sheet : sheets_componentCellIds)
//        for (const auto componentid : sheet)
//            for (const auto cell_id : baseComplex.componentC.at(componentid).cids_patch) {
//                    const Cell& cell = mesh.C.at(cell_id);
//                    ofs << cell.Vids.size();
//                    for (auto vid : cell.Vids)
//                        ofs << " " << vid;
//                    ofs << "\n";
//            }
//
//    ofs << "CELL_TYPES " << cells_num << "\n";
//    for (size_t i = 0; i < cells_num; i++)
//        ofs << "12\n";
//
//    ofs << "CELL_DATA " << cells_num << "\n"
//        << "SCALARS " << "sheet" << " int 1\n"
//        << "LOOKUP_TABLE default\n";
//    size_t sheet_id = 0;
//    for (const auto& sheet : sheets_componentCellIds) {
//        for (const auto componentid : sheet)
//            for (const auto cell_id : baseComplex.componentC.at(componentid).cids_patch)
//                ofs << sheet_id << "\n";
//        ++sheet_id;
//    }
//}

void BaseComplexSheetQuad::WriteAllSheetsFacesDualVTK(const char *filename_prefix) const
{
    RefinedDualQuad dual(baseComplex.mesh);
    dual.Build();

    for (int i = 0; i < sheets_componentFaceIds.size(); ++i) {
        std::string filename = std::string(filename_prefix) + std::to_string(i) + ".vtk";
        WriteSheetFacesDualVTK(filename.c_str(), i, dual);
    }
}

//void BaseComplexSheetQuad::WriteAllSheetsCellsVTK(const char *filename_prefix) const
//{
//    for (int i = 0; i < sheets_componentCellIds.size(); ++i) {
//        std::string filename = std::string(filename_prefix) + std::to_string(i) + ".vtk";
//        WriteSheetCellsVTK(filename.c_str(), i);
//    }
//}

void BaseComplexSheetQuad::WriteAllSheetsFacesVTK(const char *filename_prefix) const
{
    for (int i = 0; i < sheets_componentFaceIds.size(); ++i) {
        std::string filename = std::string(filename_prefix) + std::to_string(i) + ".vtk";
        WriteSheetFacesVTK(filename.c_str(), i);
    }
}

void BaseComplexSheetQuad::WriteAllSheetsEdgesVTK(const char *filename_prefix) const
{
    for (int i = 0; i < sheets_componentEdgeIds.size(); ++i) {
        std::string filename = std::string(filename_prefix) + std::to_string(i) + ".vtk";
        WriteSheetEdgesVTK(filename.c_str(), i);
    }
}

void BaseComplexSheetQuad::WriteAllSheetsFacesAndEdgesVTK(const char *filename_prefix) const
{
    for (int i = 0; i < sheets_componentEdgeIds.size(); ++i) {
        std::string filename = std::string(filename_prefix) + std::to_string(i) + ".vtk";
        WriteSheetFacesAndEdgesVTK(filename.c_str(), i);
    }
}

std::unordered_set<size_t> BaseComplexSheetQuad::GetParallelEdgeIds(const size_t sheet_id) const {
    std::unordered_set<size_t> allEdgeIds;
    for (const auto componentid : sheets_componentFaceIds.at(sheet_id))
        for (const auto faceid : baseComplex.componentF.at(componentid).fids_patch) {
            const auto& face = baseComplex.mesh.F.at(faceid);
            for (const auto edgeid : face.Eids)
                if (allEdgeIds.find(edgeid) == allEdgeIds.end()) allEdgeIds.insert(edgeid);
        }

    std::unordered_set<size_t> parallelEdgeIds;
    size_t edges_num = 0;
    for (const auto componentid : sheets_componentEdgeIds.at(sheet_id)) {
        size_t mid = baseComplex.componentE.at(componentid).eids_link.size() / 2;
        size_t edgeid = baseComplex.componentE.at(componentid).eids_link.at(mid);
        parallelEdgeIds.insert(edgeid);
        break;
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

std::unordered_set<size_t> BaseComplexSheetQuad::GetParallelSingularEdgeIds(const size_t sheet_id) const {
    std::unordered_set<size_t> parallelSingularEdgeIds;
    for (const auto componentEdgeId : sheets_componentEdgeIds.at(sheet_id))
        for (const auto edgeid : baseComplex.componentE.at(componentEdgeId).eids_link) {
            const auto& edge = baseComplex.mesh.E.at(edgeid);
            if (edge.isSingularity) parallelSingularEdgeIds.insert(edgeid);
        }
    return parallelSingularEdgeIds;
}

std::unordered_set<size_t> BaseComplexSheetQuad::GetDualEdgeIds(const size_t sheet_id) const {
    auto parallelEdgeIds = GetParallelEdgeIds(sheet_id);
    std::unordered_set<size_t> dualEdgeIds;
    for (const auto componentid : sheets_componentFaceIds.at(sheet_id))
        for (const auto faceid : baseComplex.componentF.at(componentid).fids_patch) {
            const auto& face = baseComplex.mesh.F.at(faceid);
            auto faceParallelEdgeIds = getParallelEdgeIds(baseComplex.mesh, face);
            for (auto i = 0; i < 2; ++i) {
                bool canFindTwoEdges = true;
                for (auto parallelEdgeId : faceParallelEdgeIds.at(i))
                    if (parallelEdgeIds.find(parallelEdgeId) == parallelEdgeIds.end()) {
                        canFindTwoEdges = false;
                        break;
                    }
                if (canFindTwoEdges)
                    for (auto j = 0; j < 2; ++j)
                        dualEdgeIds.insert(baseComplex.mesh.E.size() * 2 + face.id * 4 + j * 2 + i);
            }
        }
    return dualEdgeIds;
}

void BaseComplexSheetQuad::WriteSheetFacesDualVTK(const char *filename, const size_t sheet_id, const RefinedDualQuad& dual) const {
    std::unordered_set<size_t> dualEdgeIds = GetDualEdgeIds(sheet_id);
    //std::unordered_set<size_t> parallelSingularEdgeIds = GetParallelSingularEdgeIds(sheet_id);
    WriteSheetFacesDualVTK(filename, sheet_id, dual, dualEdgeIds);
    //WriteSheetFacesDualVTK(filename, sheet_id, dual, dualFaceIds, parallelSingularEdgeIds);
}

void BaseComplexSheetQuad::WriteSheetFacesDualVTK(const char *filename, const size_t sheet_id, const RefinedDualQuad& dual, const std::unordered_set<size_t>& dualEdgeIds) const
{
    const auto& mesh = baseComplex.GetMesh();
    const auto& V = dual.V;
    const auto& E = dual.E;
    std::ofstream ofs(filename);
    ofs << "# vtk DataFile Version 2.0\n"
        << filename << "\n"
        << "ASCII\n\n"
        << "DATASET POLYDATA\n";
    ofs << "POINTS " << V.size() << " double" << "\n";
    for (const auto& v : V)
        ofs << v.x << " " << v.y << " " << v.z << "\n";

    ofs << "Lines " << dualEdgeIds.size() << " " << 3 * dualEdgeIds.size() << "\n";
    for (const auto edge_id : dualEdgeIds) {
        const auto& edge = E.at(edge_id);
        ofs << edge.Vids.size() << " " << edge.Vids[0] << " " << edge.Vids[1] << "\n";
    }

    ofs << "CELL_DATA " << dualEdgeIds.size() << "\n"
        << "SCALARS " << "sheet" << " int 1\n"
        << "LOOKUP_TABLE default\n";
    for (const auto id : dualEdgeIds)
        ofs << sheet_id << "\n";
}

void BaseComplexSheetQuad::WriteSheetFacesDualVTK(const char *filename, const size_t sheet_id, const RefinedDualQuad& dual,
        const std::unordered_set<size_t>& dualEdgeIds, const std::unordered_set<size_t>& singularVertexIds) const
{
    const auto& mesh = baseComplex.GetMesh();
    const auto& V = dual.V;
    const auto& E = dual.E;
    std::ofstream ofs(filename);
    ofs << "# vtk DataFile Version 2.0\n"
        << filename << "\n"
        << "ASCII\n\n"
        << "DATASET POLYDATA\n";
    ofs << "POINTS " << V.size() << " double" << "\n";
    for (const auto& v : V)
        ofs << v.x << " " << v.y << " " << v.z << "\n";

    ofs << "Vertices " << singularVertexIds.size() << " " << 2 * singularVertexIds.size() << "\n";
    for (const auto singularVertexId : singularVertexIds)
        ofs << "1 " << singularVertexId << "\n";

    ofs << "Lines " << dualEdgeIds.size() << " " << 3 * dualEdgeIds.size() << "\n";
    for (const auto edge_id : dualEdgeIds) {
        const auto& edge = E.at(edge_id);
        ofs << edge.Vids.size() << edge.Vids[0] << " " << edge.Vids[1] << "\n";
    }

    ofs << "CELL_DATA " << singularVertexIds.size() + dualEdgeIds.size() << "\n"
        << "SCALARS " << "sheet" << " int 1\n"
        << "LOOKUP_TABLE default\n";
    for (const auto singularVertexId : singularVertexIds)
        ofs << sheet_id << "\n";
    for (const auto id : dualEdgeIds)
        ofs << sheet_id << "\n";
}

//void BaseComplexSheetQuad::WriteSheetCellsVTK(const char *filename, const size_t sheet_id) const
//{
//    const auto& mesh = baseComplex.GetMesh();
//    const auto& V = mesh.V;
//    const auto& E = mesh.E;
//    std::ofstream ofs(filename);
//    ofs << "# vtk DataFile Version 2.0\n"
//        << filename << "\n"
//        << "ASCII\n\n"
//        << "DATASET UNSTRUCTURED_GRID\n";
//    ofs << "POINTS " << V.size() << " float" << "\n";
//    for (size_t i = 0; i < V.size(); i++)
//        ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << "\n";
//
//    size_t cells_num = 0;
//    for (const auto componentid : sheets_componentCellIds.at(sheet_id))
//        cells_num += baseComplex.componentC.at(componentid).cids_patch.size();
//
//    ofs << "CELLS " << cells_num << " " << 9 * cells_num << "\n";
//    for (const auto componentid : sheets_componentCellIds.at(sheet_id))
//        for (const auto cell_id : baseComplex.componentC.at(componentid).cids_patch) {
//                const Cell& cell = mesh.C.at(cell_id);
//                ofs << cell.Vids.size();
//                for (auto vid : cell.Vids)
//                    ofs << " " << vid;
//                ofs << "\n";
//        }
//
//    ofs << "CELL_TYPES " << cells_num << "\n";
//    for (size_t i = 0; i < cells_num; i++)
//        ofs << "12\n";
//
//    ofs << "CELL_DATA " << cells_num << "\n";
//    ofs << "SCALARS " << "component" << " int 1\n"
//        << "LOOKUP_TABLE default\n";
//    for (const auto componentid : sheets_componentCellIds.at(sheet_id))
//        for (const auto cell_id : baseComplex.componentC.at(componentid).cids_patch)
//            ofs << baseComplex.componentC.at(componentid).color % 8 << "\n";
//    ofs << "SCALARS " << "raw_component_id" << " int 1\n"
//        << "LOOKUP_TABLE default\n";
//    for (const auto componentid : sheets_componentCellIds.at(sheet_id))
//        for (const auto cell_id : baseComplex.componentC.at(componentid).cids_patch)
//            ofs << componentid << "\n";
//    ofs << "SCALARS " << "component_id" << " int 1\n"
//        << "LOOKUP_TABLE default\n";
//    int sheet_component_id = 0;
//    for (const auto componentid : sheets_componentCellIds.at(sheet_id)) {
//        for (const auto cell_id : baseComplex.componentC.at(componentid).cids_patch)
//            ofs << sheet_component_id << "\n";
//        ++sheet_component_id;
//    }
//    ofs << "SCALARS " << "sheet" << " int 1\n"
//        << "LOOKUP_TABLE default\n";
//    for (const auto componentid : sheets_componentCellIds.at(sheet_id))
//        for (const auto cell_id : baseComplex.componentC.at(componentid).cids_patch)
//            ofs << sheet_id << "\n";
//}

void BaseComplexSheetQuad::WriteSheetFacesVTK(const char *filename, const size_t sheet_id) const
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

void BaseComplexSheetQuad::WriteSheetEdgesVTK(const char *filename, const size_t sheet_id) const
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
        for (const auto cell_id : baseComplex.componentE.at(componentid).eids_link)
            ofs << sheet_id << "\n";
}

void BaseComplexSheetQuad::WriteSheetFacesAndEdgesVTK(const char *filename, const size_t sheet_id) const
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

enum SheetsConnectivity {
    SheetsConnectivity_UNKNOWN = 0,
    SheetsConnectivity_NEIGHBOR,
    SheetsConnectivity_INTERSECT,
    SheetsConnectivity_NEIGHBOR_AND_INTERSECT
};

void BaseComplexSheetQuad::ExtractSheetConnectivities() {
    size_t n = sheets_componentFaceIds.size();
    sheets_connectivities.resize(n, std::vector<size_t>(n, 0));

    if (n < 1) return;
    std::vector<std::unordered_set<size_t>> sheets_hashcomponentEdgeIds;  // each sheet consists a list of base-complex componentEdge ids;
    std::vector<std::unordered_set<size_t>> sheets_hashcomponentFaceIds;  // each sheet consists a list of base-complex componentFace ids;
    //std::vector<std::unordered_set<size_t>> sheets_hashcomponentCellIds;  // each sheet consists a list of base-complex componentCell ids;

    for (auto & sheet_componentEdgeIds : sheets_componentEdgeIds)
        sheets_hashcomponentEdgeIds.push_back(std::unordered_set<size_t> (sheet_componentEdgeIds.begin(), sheet_componentEdgeIds.end()));
    for (auto & sheet_componentFaceIds : sheets_componentFaceIds)
        sheets_hashcomponentFaceIds.push_back(std::unordered_set<size_t> (sheet_componentFaceIds.begin(), sheet_componentFaceIds.end()));
    // for (auto & sheet_componentCellIds : sheets_componentCellIds)
    //     sheets_hashcomponentCellIds.push_back(std::unordered_set<size_t> (sheet_componentCellIds.begin(), sheet_componentCellIds.end()));
    // Get intersections
    for (const auto& sheetIds : faceComponent_sheetIds)
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
    for (auto& sheetBoundaryEdgeComponentIds : sheets_BoundaryEdgeComponentIds) {
        std::unordered_set<size_t> neighborSheetIds;
        for (auto sheetBoundaryEdgeComponentId : sheetBoundaryEdgeComponentIds) {
            const auto & sheetIds = edgeComponent_neighborSheetIds[sheetBoundaryEdgeComponentId];
            if (sheetIds.size() <= 1) continue;
            if (baseComplex.componentE.at(sheetBoundaryEdgeComponentId).isBoundary) continue;

            for (auto neighborSheetId : edgeComponent_neighborSheetIds[sheetBoundaryEdgeComponentId])
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
}

//void BaseComplexSheetQuad::ExtractMainSheetConnectivities() {
//    all_sheets_connectivities = sheets_connectivities;
//    size_t n = sheets_componentFaceIds.size();
//    std::vector<std::vector<size_t>> new_sheets_connectivities;
//    std::vector<size_t> sheetIds;
//    std::vector<bool> active(n, false);
//    for (auto id : sheets_coverSheetIds[0])
//        active[id] = true;
//    for (size_t i = 0; i < n; ++i)
//        if (active[i]) new_sheets_connectivities.push_back(all_sheets_connectivities[i]);
//
//    for (auto& row : new_sheets_connectivities) {
//        std::vector<size_t> new_sheet_connectivities;
//        for (size_t i = 0; i < n; ++i)
//            if (active[i]) new_sheet_connectivities.push_back(row[i]);
//        row = new_sheet_connectivities;
//    }
//    sheets_connectivities = new_sheets_connectivities;
//
//    n = sheets_connectivities.size();
//    sheets_connectivities_float.resize(n);
//    for (auto& row : sheets_connectivities_float) {
//        row.clear();
//        row.resize(n, 0.0f);
//    }
//    for (size_t i = 0; i < n; ++i)
//        for (size_t j = 0; j < n; ++j)
//            if (sheets_connectivities[i][j] == SheetsConnectivity_NEIGHBOR) sheets_connectivities_float[i][j] = -1.0f;
//            else if (sheets_connectivities[i][j] == SheetsConnectivity_INTERSECT) sheets_connectivities_float[i][j] = 1.0f;
//            else if (sheets_connectivities[i][j] == SheetsConnectivity_NEIGHBOR_AND_INTERSECT) sheets_connectivities_float[i][j] = 0.5f;
//}

void BaseComplexSheetQuad::WriteSheetsConnectivitiesMatrixVTK(const char *filename) const {
    MeshFileWriter writer(baseComplex.mesh, filename);
    writer.WriteMatrixVTK(sheets_connectivities);
}

void BaseComplexSheetQuad::WriteSheetsConnectivitiesMatrixMat(const char* filename) const {
    std::ofstream ofs(filename);
    const size_t n = sheets_connectivities_float.size();
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j)
            ofs << sheets_connectivities_float[i][j] << "\t";
        ofs << "\n";
    }
}

const std::vector<std::vector<size_t>>& BaseComplexSheetQuad::Get_sheets_coverSheetIds() const {
    return sheets_coverSheetIds;
}

const size_t BaseComplexSheetQuad::GetNumOfSheets() const {
    return sheets_componentFaceIds.size();
}


void BaseComplexSheetQuad::ComputeComplexity() {
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
    std::cout << "Complexity of main sheets = " << complexity << "\n";
    std::cout << "##############################\n";

    //std::unordered_set<size_t> hash_main_sheet_ids(main_sheet_ids.begin(), main_sheet_ids.end());
    sheets_importance.reserve(main_sheet_ids.size());
    for (auto main_sheet_id : main_sheet_ids) {
        size_t importance = 0;
        for (auto id : main_sheet_ids)
            importance += all_sheets_connectivities[main_sheet_id][id] != SheetsConnectivity_UNKNOWN;
        sheets_importance.push_back(float(importance)/complexity);
    }
    std::cout << "##############################\n";
    std::cout << "Complexity of each sheet : ";
    for (auto sheet_importance : sheets_importance)
        std::cout << sheet_importance << " ";
    std::cout << "\n##############################\n";
}
bool BaseComplexSheetQuad::IsAdjacent(const int component_id1, const int component_id2) const {
    const auto& component1 = baseComplex.componentF.at(component_id1);
    const auto& component2 = baseComplex.componentF.at(component_id2);

    std::unordered_set<size_t> s(component1.Eids.begin(), component1.Eids.end());
    s.insert(component2.Eids.begin(), component2.Eids.end());

    return s.size() == 7;
}
void BaseComplexSheetQuad::ComputeComplexityDrChen(int sheetid) {
    enum Relation {
        Relation_UNKNOWN = 0,
        ADJACENT,
        INTERSECTING,
        HYBRID
    };

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
            const auto& sheet1_component_edge_ids = sheets_BoundaryEdgeComponentIds[sheetid1];
            const auto& sheet2_component_edg_ids = sheets_BoundaryEdgeComponentIds[sheetid2];
            std::unordered_set<size_t> sheet1_component_edg_ids_set(sheet1_component_edge_ids.begin(), sheet1_component_edge_ids.end());
            std::unordered_set<size_t> common_component_edg_ids;
            for (auto component_edg_id : sheet2_component_edg_ids)
                if (sheet1_component_edg_ids_set.find(component_edg_id) != sheet1_component_edg_ids_set.end())
                    common_component_edg_ids.insert(component_edg_id);

            auto common_component_size = common_component_edg_ids.size() * 2;
            auto total_component_size = sheets_componentFaceIds[sheetid1].size() + sheets_componentFaceIds[sheetid2].size();
            auto relation_float = (float(common_component_size))/total_component_size;
            sheets_connectivities_float[p[0]][p[1]] = relation_float;
            sheets_connectivities_float[p[1]][p[0]] = relation_float;
        } else if (relation == INTERSECTING) {
            const auto& sheet1_component_face_ids = sheets_componentFaceIds[sheetid1];
            const auto& sheet2_component_face_ids = sheets_componentFaceIds[sheetid2];
            std::unordered_set<size_t> sheet1_component_face_ids_set(sheet1_component_face_ids.begin(), sheet1_component_face_ids.end());
            std::unordered_set<size_t> common_component_face_ids;
            for (auto component_face_id : sheet2_component_face_ids)
                if (sheet1_component_face_ids_set.find(component_face_id) != sheet1_component_face_ids_set.end())
                    common_component_face_ids.insert(component_face_id);

            size_t common_component_face_ids_size = common_component_face_ids.size();
            auto relation_float = common_component_face_ids.size();
            if (common_component_face_ids_size > 1) {
                std::vector<size_t> common_component_face_ids_array(common_component_face_ids.begin(), common_component_face_ids.end());
                auto combs = combine(common_component_face_ids_size, 2);
                for (auto& comb : combs) {
                    auto component_id1 = common_component_face_ids_array[comb[0]];
                    auto component_id2 = common_component_face_ids_array[comb[1]];
                    if (IsAdjacent(component_id1, component_id2)) --relation_float;
                }
            }
            ++relation_float;
            if (relation_float > max_num_of_intersections) max_num_of_intersections = relation_float;

            sheets_connectivities_float[p[0]][p[1]] = relation_float;
            sheets_connectivities_float[p[1]][p[0]] = relation_float;
        } else if (relation == HYBRID) {
            const auto& sheet1_component_edge_ids = sheets_BoundaryEdgeComponentIds[sheetid1];
            const auto& sheet2_component_edg_ids = sheets_BoundaryEdgeComponentIds[sheetid2];
            std::unordered_set<size_t> sheet1_component_edg_ids_set(sheet1_component_edge_ids.begin(), sheet1_component_edge_ids.end());
            std::unordered_set<size_t> common_component_edg_ids;
            for (auto component_edg_id : sheet2_component_edg_ids)
                if (sheet1_component_edg_ids_set.find(component_edg_id) != sheet1_component_edg_ids_set.end())
                    common_component_edg_ids.insert(component_edg_id);

            auto common_component_size = common_component_edg_ids.size() * 2;
            auto total_component_size = sheets_componentFaceIds[sheetid1].size() + sheets_componentFaceIds[sheetid2].size();
            auto relation_float = (float(common_component_size))/total_component_size;

            const auto& sheet1_component_face_ids = sheets_componentFaceIds[sheetid1];
            const auto& sheet2_component_face_ids = sheets_componentFaceIds[sheetid2];
            std::unordered_set<size_t> sheet1_component_face_ids_set(sheet1_component_face_ids.begin(), sheet1_component_face_ids.end());
            std::unordered_set<size_t> common_component_face_ids;
            for (auto component_face_id : sheet2_component_face_ids)
                if (sheet1_component_face_ids_set.find(component_face_id) != sheet1_component_face_ids_set.end())
                    common_component_face_ids.insert(component_face_id);

            size_t common_component_face_ids_size = common_component_face_ids.size();
            if (common_component_face_ids_size > 1) {
                std::vector<size_t> common_component_face_ids_array(common_component_face_ids.begin(), common_component_face_ids.end());
                auto combs = combine(common_component_face_ids_size, 2);
                for (auto& comb : combs) {
                    auto component_id1 = common_component_face_ids_array[comb[0]];
                    auto component_id2 = common_component_face_ids_array[comb[1]];
                    if (IsAdjacent(component_id1, component_id2)) --common_component_face_ids_size;
                }
            }

            relation_float = (1 + relation_float) * common_component_face_ids_size;
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
//    std::vector<size_t> sheets_complexities;
//    std::string line;
//    while (getline(ifs, line)) {
//        std::stringstream ss(line);
//        std::string str;
//        while (ss >> str);
//        sheets_complexities.push_back(std::stoi(str));
//    }
//    if (sheets_connectivities.size() == 1) {
//        std::cout << "\n\n\n ********** COMPLEXITY *********** \n";
//        std::cout << "#COMPLEXITY = " << sheets_complexities[sheets_coverSheetIds[0][0]];
//        std::cout << "\n\n\n ********** COMPLEXITY *********** \n";
//        sheets_connectivities[0][0] = sheets_complexities[sheets_coverSheetIds[0][0]];
//    }

    if (sheets_connectivities.size() == 1) {
        std::cout << "\n\n\n ********** COMPLEXITY *********** \n";
        std::cout << "#COMPLEXITY = " << 1;
        std::cout << "\n\n\n ********** COMPLEXITY *********** \n";
        sheets_connectivities[0][0] = 1;
    }

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

    auto sheets_connectivities_float_new = sheets_connectivities_float;
    auto sheets_connectivities_new = sheets_connectivities;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            sheets_connectivities_float_new[i][j] = sheets_connectivities_float[ss[i].s][ss[j].s];
            sheets_connectivities_new[i][j] = sheets_connectivities[ss[i].s][ss[j].s];
        }

//    sheets_connectivities_float = sheets_connectivities_float_new;
//    sheets_connectivities = sheets_connectivities_new;

    /////////////////////////////////
    {
        std::ofstream ofs("complexity.mat");
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j)
                ofs << sheets_connectivities_float[i][j] << "\t";
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


void BaseComplexSheetQuad::ExtractMainSheetConnectivities(int main_sheets_id) {
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
