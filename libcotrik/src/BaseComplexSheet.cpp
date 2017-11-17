/*
 * BaseComplexSheet.cpp
 *
 *  Created on: Oct 5, 2017
 *      Author: cotrik
 */

#include "BaseComplexSheet.h"
#include "MeshFileWriter.h"
#include "Dual.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <queue>

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

//void BaseComplexSheet::ExtractSets()
//{
//        std::vector<size_t> roots(baseComplex.componentE.size(), 0);
//        for (auto i = 0; i < baseComplex.componentE.size(); ++i) roots[i] = i;
//        for (auto& componentEdge : baseComplex.componentE) {
//            const auto parallelComponentEdgeIds = GetParallelComponentEdgeIds(componentEdge);
//            const auto copys = parallelComponentEdgeIds;
//            for (auto& copy : copys) {
//                while (copy != roots[copy]) copy = roots[copy];
//            }
//            for (int u = 1; u < copys.size(); ++u)
//                for (int v = 0; v < u; ++v)
//                    if (copys[u] == copys[v])
//            int f = e.first;
//            int s = e.second;
//            while (f != roots[f]) f = roots[f];
//            while (s != roots[s]) s = roots[s];
//            if (roots[f] != roots[s]) {
//                --n;
//                roots[s] = f;
//            }
//        }
//        return n;
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
    for (int i = 0; i < sheets_componentCellIds.size(); ++i) {
        std::string filename = std::string(filename_prefix) + std::to_string(i) + ".vtk";
        WriteSheetCellsDualVTK(filename.c_str(), i);
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
        for (const auto cell_id : baseComplex.componentE.at(componentid).eids_link)
            ofs << sheet_id << "\n";
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
}
