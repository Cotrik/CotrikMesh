/*
 * BaseComplexChord.cpp
 *
 *  Created on: Oct 5, 2017
 *      Author: cotrik
 */

#include "BaseComplexChord.h"
#include "MeshFileWriter.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <stack>

BaseComplexChord::BaseComplexChord(BaseComplex& baseComplex)
: baseComplex(baseComplex)
{
    // TODO Auto-generated constructor stub

}

BaseComplexChord::~BaseComplexChord()
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
void BaseComplexChord::Extract()
{
    std::vector<bool> visited(baseComplex.componentF.size(), false);
    for (auto& componentFace : baseComplex.componentF) {
        if (visited.at(componentFace.id)) continue;
        std::vector<size_t> chordComponentEdgeIds;
        std::vector<size_t> chordComponentFaceIds;
        std::vector<size_t> chordComponentCellIds;
        GetParallelComponents(componentFace, chordComponentEdgeIds, chordComponentFaceIds, chordComponentCellIds);
        chords_componentEdgeIds.push_back(chordComponentEdgeIds);
        chords_componentFaceIds.push_back(chordComponentFaceIds);
        chords_componentCellIds.push_back(chordComponentCellIds);

        for (auto id : chordComponentFaceIds) visited.at(id) = true;
    }
}

void BaseComplexChord::GetParallelComponents(const ComponentFace & componentFace,
        std::vector<size_t>& chordComponentEdgeIds, std::vector<size_t>& chordComponentFaceIds, std::vector<size_t>& chordComponentCellIds)
{
    std::vector<bool> edge_visited(baseComplex.componentE.size(), false);
    std::vector<bool> face_visited(baseComplex.componentF.size(), false);
    std::vector<bool> cell_visited(baseComplex.componentC.size(), false);
    face_visited.at(componentFace.id) = true;
    chordComponentFaceIds.push_back(componentFace.id);
    std::queue<size_t> q;
    q.push(componentFace.id);
    while (!q.empty()) {
        size_t n = q.size();
        for (int i = 0; i < n; ++i) {
            const size_t componentFaceId = q.front();
            q.pop();

            for (auto neighborComponentCellId : baseComplex.componentF.at(componentFaceId).N_Cids)
                if (!cell_visited.at(neighborComponentCellId)) {
                    cell_visited.at(neighborComponentCellId) = true;
                    chordComponentCellIds.push_back(neighborComponentCellId);
                    for (auto componentEdgeId : baseComplex.componentC.at(neighborComponentCellId).Eids)
                        if (!edge_visited.at(componentEdgeId)) {
                            edge_visited.at(componentEdgeId) = true;
                            chordComponentEdgeIds.push_back(componentEdgeId);
                        }
                }
            const auto nextComponentFaceIds = GetParallelComponentFaceIds(baseComplex.componentF.at(componentFaceId));
            for (auto nextComponentFaceId : nextComponentFaceIds)
                if (!face_visited.at(nextComponentFaceId)) {
                    face_visited.at(nextComponentFaceId) = true;
                    chordComponentFaceIds.push_back(nextComponentFaceId);
                    q.push(nextComponentFaceId);
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
std::vector<size_t> BaseComplexChord::GetParallelComponentFaceIds(const ComponentFace & componentFace) const{
    std::vector<size_t> res;
    for (auto neighborComponentCellId : baseComplex.componentF.at(componentFace.id).N_Cids) {
        for (auto componentFaceId : baseComplex.componentC.at(neighborComponentCellId).Fids) {
            auto& componentCellFace = baseComplex.componentF.at(componentFaceId);
            auto vids = baseComplex.componentF.at(componentFace.id).Vids;
            std::copy(componentCellFace.Vids.begin(), componentCellFace.Vids.end(), back_inserter(vids));
            std::sort(vids.begin(), vids.end());
            auto num = std::distance(vids.begin(), std::unique(vids.begin(), vids.end()));
            if (num == 8) res.push_back(componentFaceId);
        }
    }
    return res;
}

void BaseComplexChord::WriteChordsEdgesVTK(const char *filename) const
{

}

void BaseComplexChord::WriteChordsFacesVTK(const char *filename) const
{

}

void BaseComplexChord::WriteChordsCellsVTK(const char *filename) const
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
    for (const auto& sheet : chords_componentCellIds)
        for (const auto componentid : sheet)
            cells_num += baseComplex.componentC.at(componentid).cids_patch.size();

    ofs << "CELLS " << cells_num << " " << 9 * cells_num << "\n";
    for (const auto& sheet : chords_componentCellIds)
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
        << "SCALARS " << "chord" << " int 1\n"
        << "LOOKUP_TABLE default\n";
    size_t chord_id = 0;
    for (const auto& sheet : chords_componentCellIds) {
        for (const auto componentid : sheet)
            for (const auto cell_id : baseComplex.componentC.at(componentid).cids_patch)
                ofs << chord_id << "\n";
        ++chord_id;
    }
}

void BaseComplexChord::WriteAllChordsCellsVTK(const char *filename_prefix) const
{
    for (int i = 0; i < chords_componentCellIds.size(); ++i) {
        std::string filename = std::string(filename_prefix) + std::to_string(i) + ".vtk";
        WriteChordCellsVTK(filename.c_str(), i);
    }
}

void BaseComplexChord::WriteAllChordsFacesVTK(const char *filename_prefix) const
{
    for (int i = 0; i < chords_componentFaceIds.size(); ++i) {
        std::string filename = std::string(filename_prefix) + std::to_string(i) + ".vtk";
        WriteChordFacesVTK(filename.c_str(), i);
    }
}

void BaseComplexChord::WriteAllChordsEdgesVTK(const char *filename_prefix) const
{
    for (int i = 0; i < chords_componentEdgeIds.size(); ++i) {
        std::string filename = std::string(filename_prefix) + std::to_string(i) + ".vtk";
        WriteChordEdgesVTK(filename.c_str(), i);
    }
}

void BaseComplexChord::WriteAllChordsFacesAndEdgesVTK(const char *filename_prefix) const
{
    for (int i = 0; i < chords_componentEdgeIds.size(); ++i) {
        std::string filename = std::string(filename_prefix) + std::to_string(i) + ".vtk";
        WriteChordFacesAndEdgesVTK(filename.c_str(), i);
    }
}

void BaseComplexChord::WriteAllChordsCurvesVTK(const char *filename_prefix) const
{
    for (int i = 0; i < chords_componentEdgeIds.size(); ++i) {
        std::string filename = std::string(filename_prefix) + std::to_string(i) + ".vtk";
        WriteChordCurvesVTK(filename.c_str(), i);
    }
}

void BaseComplexChord::WriteAllChordsFramesVTK(const char *filename_prefix) const
{
    for (int i = 0; i < chords_componentEdgeIds.size(); ++i) {
        std::string filename = std::string(filename_prefix) + std::to_string(i) + ".vtk";
        WriteChordFramesVTK(filename.c_str(), i);
    }
}

void BaseComplexChord::WriteChordCellsVTK(const char *filename, const size_t chord_id) const
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
    for (const auto componentid : chords_componentCellIds.at(chord_id))
        cells_num += baseComplex.componentC.at(componentid).cids_patch.size();

    ofs << "CELLS " << cells_num << " " << 9 * cells_num << "\n";
    for (const auto componentid : chords_componentCellIds.at(chord_id))
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
    for (const auto componentid : chords_componentCellIds.at(chord_id))
        for (const auto cell_id : baseComplex.componentC.at(componentid).cids_patch)
            ofs << baseComplex.componentC.at(componentid).color % 8 << "\n";
    ofs << "SCALARS " << "raw_component_id" << " int 1\n"
        << "LOOKUP_TABLE default\n";
    for (const auto componentid : chords_componentCellIds.at(chord_id))
        for (const auto cell_id : baseComplex.componentC.at(componentid).cids_patch)
            ofs << componentid << "\n";
    ofs << "SCALARS " << "component_id" << " int 1\n"
        << "LOOKUP_TABLE default\n";
    int chord_component_id = 0;
    for (const auto componentid : chords_componentCellIds.at(chord_id)) {
        for (const auto cell_id : baseComplex.componentC.at(componentid).cids_patch)
            ofs << chord_component_id << "\n";
        ++chord_component_id;
    }
    ofs << "SCALARS " << "chord" << " int 1\n"
        << "LOOKUP_TABLE default\n";
    for (const auto componentid : chords_componentCellIds.at(chord_id))
        for (const auto cell_id : baseComplex.componentC.at(componentid).cids_patch)
            ofs << chord_id << "\n";
}

void BaseComplexChord::WriteChordFacesVTK(const char *filename, const size_t chord_id) const
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
    for (const auto componentid : chords_componentFaceIds.at(chord_id))
        faces_num += baseComplex.componentF.at(componentid).fids_patch.size();

    ofs << "Polygons " << faces_num << " " << 5 * faces_num << "\n";
    for (const auto componentid : chords_componentFaceIds.at(chord_id))
        for (const auto cell_id : baseComplex.componentF.at(componentid).fids_patch) {
                const Face& cell = mesh.F.at(cell_id);
                ofs << cell.Vids.size();
                for (auto vid : cell.Vids)
                    ofs << " " << vid;
                ofs << "\n";
        }

    ofs << "CELL_DATA " << faces_num << "\n"
        << "SCALARS " << "chord" << " int 1\n"
        << "LOOKUP_TABLE default\n";
    for (const auto componentid : chords_componentFaceIds.at(chord_id))
        for (const auto cell_id : baseComplex.componentF.at(componentid).fids_patch)
            ofs << chord_id << "\n";
}

void BaseComplexChord::WriteChordEdgesVTK(const char *filename, const size_t chord_id) const
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
    for (const auto componentid : chords_componentEdgeIds.at(chord_id))
        edges_num += baseComplex.componentE.at(componentid).eids_link.size();

    ofs << "Lines " << edges_num << " " << 3 * edges_num << "\n";
    for (const auto componentid : chords_componentEdgeIds.at(chord_id))
        for (const auto cell_id : baseComplex.componentE.at(componentid).eids_link) {
                const Edge& cell = mesh.E.at(cell_id);
                ofs << cell.Vids.size();
                for (auto vid : cell.Vids)
                    ofs << " " << vid;
                ofs << "\n";
        }

    ofs << "CELL_DATA " << edges_num << "\n"
        << "SCALARS " << "chord" << " int 1\n"
        << "LOOKUP_TABLE default\n";
    for (const auto componentid : chords_componentEdgeIds.at(chord_id))
        for (const auto cell_id : baseComplex.componentE.at(componentid).eids_link)
            ofs << chord_id << "\n";
}

void BaseComplexChord::WriteChordFacesAndEdgesVTK(const char *filename, const size_t chord_id) const
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
    auto linkedChordComponentFaceIds = GetLinkedChordComponentFaceIds(chords_componentFaceIds.at(chord_id));
    size_t faces_num = 0;
        for (const auto componentid : chords_componentFaceIds.at(chord_id))
            faces_num += baseComplex.componentF.at(componentid).fids_patch.size();
    size_t edges_num = 0;
        for (const auto componentid : chords_componentEdgeIds.at(chord_id))
            edges_num += baseComplex.componentE.at(componentid).eids_link.size();

    ofs << "Lines " << edges_num << " " << 3 * edges_num << "\n";
    for (const auto componentid : chords_componentEdgeIds.at(chord_id))
        for (const auto cell_id : baseComplex.componentE.at(componentid).eids_link) {
                const Edge& cell = mesh.E.at(cell_id);
                ofs << cell.Vids.size();
                for (auto vid : cell.Vids)
                    ofs << " " << vid;
                ofs << "\n";
        }
    ofs << "Polygons " << faces_num << " " << 5 * faces_num << "\n";
    for (const auto componentid : linkedChordComponentFaceIds)
        for (const auto cell_id : baseComplex.componentF.at(componentid).fids_patch) {
                const Face& cell = mesh.F.at(cell_id);
                ofs << cell.Vids.size();
                for (auto vid : cell.Vids)
                    ofs << " " << vid;
                ofs << "\n";
        }

    ofs << "CELL_DATA " << edges_num + faces_num << "\n"
        << "SCALARS " << "chord" << " int 1\n"
        << "LOOKUP_TABLE default\n";
    for (const auto componentid : chords_componentEdgeIds.at(chord_id))
        for (const auto cell_id : baseComplex.componentE.at(componentid).eids_link)
            ofs << chord_id << "\n";
    for (const auto componentid : chords_componentFaceIds.at(chord_id))
        for (const auto cell_id : baseComplex.componentF.at(componentid).fids_patch)
            ofs << chord_id << "\n";
    ofs << "SCALARS " << "chordlink" << " int 1\n"
        << "LOOKUP_TABLE default\n";
    for (const auto componentid : chords_componentEdgeIds.at(chord_id))
        for (const auto cell_id : baseComplex.componentE.at(componentid).eids_link)
            ofs << 0 << "\n";
    int linkid = 0;
    for (const auto componentid : linkedChordComponentFaceIds) {
        for (const auto cell_id : baseComplex.componentF.at(componentid).fids_patch)
            ofs << linkid << "\n";
        ++linkid;
    }
}

//void BaseComplexChord::WriteChordCurvesVTK(const char *filename, const size_t chord_id) const
//{
//    auto linkedChordComponentFaceIds = GetLinkedChordComponentFaceIds(chords_componentFaceIds.at(chord_id));
//    auto linkedChordComponentCellIds = GetLinkedChordComponentCellIds(chords_componentCellIds.at(chord_id), linkedChordComponentFaceIds);
//    auto curveVertices = GetLinkedChordCurveVertices(linkedChordComponentFaceIds, linkedChordComponentCellIds);
//    const auto& mesh = baseComplex.GetMesh();
//    const auto& V = mesh.V;
//    std::ofstream ofs(filename);
//    ofs << "# vtk DataFile Version 2.0\n"
//        << filename << "\n"
//        << "ASCII\n\n"
//        << "DATASET POLYDATA\n";
//    ofs << "POINTS " << V.size() + curveVertices.size() << " double" << "\n";
//    for (auto& v : V)
//        ofs << v.x << " " << v.y << " " << v.z << "\n";
//    for (auto& v : curveVertices)
//        ofs << v.x << " " << v.y << " " << v.z << "\n";
//
//    size_t faces_num = 0;
//        for (const auto componentid : chords_componentFaceIds.at(chord_id))
//            faces_num += baseComplex.componentF.at(componentid).fids_patch.size();
////    size_t edges_num = 0;
////        for (const auto componentid : chords_componentEdgeIds.at(chord_id))
////            edges_num += baseComplex.componentE.at(componentid).eids_link.size();
//
//    ofs << "Lines " << 1 << " " << 1 + curveVertices.size() << "\n";
//    ofs << curveVertices.size() << " ";
//    for (auto i = 0; i < curveVertices.size(); ++i)
//        ofs << " " << V.size() + i;
//    ofs << "\n";
//
//    ofs << "Polygons " << faces_num << " " << 5 * faces_num << "\n";
//    for (const auto componentid : linkedChordComponentFaceIds)
//        for (const auto cell_id : baseComplex.componentF.at(componentid).fids_patch) {
//                const Face& cell = mesh.F.at(cell_id);
//                ofs << cell.Vids.size();
//                for (auto vid : cell.Vids)
//                    ofs << " " << vid;
//                ofs << "\n";
//        }
//
//    ofs << "CELL_DATA " << 1 + faces_num << "\n"
//        << "SCALARS " << "chord" << " int 1\n"
//        << "LOOKUP_TABLE default\n";
//    ofs << chord_id << "\n";
//    for (const auto componentid : chords_componentFaceIds.at(chord_id))
//        for (const auto cell_id : baseComplex.componentF.at(componentid).fids_patch)
//            ofs << chord_id << "\n";
//    ofs << "SCALARS " << "chordlink" << " int 1\n"
//        << "LOOKUP_TABLE default\n";
//    ofs << 0 << "\n";
//    int linkid = 0;
//    for (const auto componentid : linkedChordComponentFaceIds) {
//        for (const auto cell_id : baseComplex.componentF.at(componentid).fids_patch)
//            ofs << linkid << "\n";
//        ++linkid;
//    }
//}

void BaseComplexChord::WriteChordCurvesVTK(const char *filename, const size_t chord_id) const
{
    auto linkedChordComponentFaceIds = GetLinkedChordComponentFaceIds(chords_componentFaceIds.at(chord_id));
    auto linkedChordComponentCellIds = GetLinkedChordComponentCellIds(chords_componentCellIds.at(chord_id), linkedChordComponentFaceIds);
    auto curveVertices = GetLinkedChordCurveVertices(linkedChordComponentFaceIds, linkedChordComponentCellIds);
    std::ofstream ofs(filename);
    ofs << "# vtk DataFile Version 2.0\n"
        << filename << "\n"
        << "ASCII\n\n"
        << "DATASET POLYDATA\n";
    ofs << "POINTS " << curveVertices.size() << " double" << "\n";
    for (auto& v : curveVertices)
        ofs << v.x << " " << v.y << " " << v.z << "\n";
    ofs << "Lines " << 1 << " " << 1 + curveVertices.size() << "\n";
    ofs << curveVertices.size() << " ";
    for (auto i = 0; i < curveVertices.size(); ++i)
        ofs << " " << i;
    ofs << "\n";
    ofs << "CELL_DATA " << curveVertices.size() << "\n"
        << "SCALARS " << "chord" << " int 1\n"
        << "LOOKUP_TABLE default\n";
    for (auto i = 0; i < curveVertices.size(); ++i)
        ofs << chord_id << "\n";
    ofs << "SCALARS " << "chordlink" << " int 1\n"
        << "LOOKUP_TABLE default\n";
    ofs << 0 << "\n";
    for (auto i = 0; i < curveVertices.size(); ++i)
        ofs << i << "\n";
}

void BaseComplexChord::WriteChordFramesVTK(const char *filename, const size_t chord_id) const
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

    auto linkedChordComponentFaceIds = GetLinkedChordComponentFaceIds(chords_componentFaceIds.at(chord_id));
//    std::cout << "linkedChordComponentFaceIds : " << linkedChordComponentFaceIds.size() << "\n";
//    std::cout << "chords_componentFaceIds.at(chord_id) : " << chords_componentFaceIds.at(chord_id).size() << "\n";
    size_t faces_num = 0;
        for (const auto componentid : chords_componentFaceIds.at(chord_id))
            faces_num += baseComplex.componentF.at(componentid).fids_patch.size();
    size_t edges_num = 0;
        for (const auto componentid : chords_componentEdgeIds.at(chord_id))
            edges_num += baseComplex.componentE.at(componentid).eids_link.size();

    ofs << "Lines " << edges_num << " " << 3 * edges_num << "\n";
    for (const auto componentid : chords_componentEdgeIds.at(chord_id))
        for (const auto cell_id : baseComplex.componentE.at(componentid).eids_link) {
                const Edge& cell = mesh.E.at(cell_id);
                ofs << cell.Vids.size();
                for (auto vid : cell.Vids)
                    ofs << " " << vid;
                ofs << "\n";
        }
    ofs << "Polygons " << faces_num << " " << 5 * faces_num << "\n";
    for (const auto componentid : linkedChordComponentFaceIds)
        for (const auto cell_id : baseComplex.componentF.at(componentid).fids_patch) {
                const Face& cell = mesh.F.at(cell_id);
                ofs << cell.Vids.size();
                for (auto vid : cell.Vids)
                    ofs << " " << vid;
                ofs << "\n";
        }

    ofs << "CELL_DATA " << edges_num + faces_num << "\n";
    ofs << "SCALARS " << "chord" << " int 1\n"
        << "LOOKUP_TABLE default\n";
    for (const auto componentid : linkedChordComponentFaceIds)
        for (const auto cell_id : baseComplex.componentF.at(componentid).fids_patch)
            ofs << chord_id << "\n";
    for (const auto componentid : chords_componentEdgeIds.at(chord_id))
        for (const auto cell_id : baseComplex.componentE.at(componentid).eids_link)
            ofs << chord_id << "\n";
    ofs << "SCALARS " << "chordlink" << " int 1\n"
        << "LOOKUP_TABLE default\n";
    for (const auto componentid : chords_componentEdgeIds.at(chord_id))
        for (const auto cell_id : baseComplex.componentE.at(componentid).eids_link)
            ofs << 0 << "\n";
    int linkid = 0;
    for (const auto componentid : linkedChordComponentFaceIds) {
        for (const auto cell_id : baseComplex.componentF.at(componentid).fids_patch)
            ofs << linkid << "\n";
        ++linkid;
    }

}

std::vector<size_t> BaseComplexChord::GetLinkedChordComponentCellIds(const std::vector<size_t>& chordComponentCellIds, const std::vector<size_t>& linkedChordComponentFaceIds) const{
    std::vector<size_t> res;
    bool isLoop = true;
    size_t startComponentFaceId = linkedChordComponentFaceIds.at(0);
    for (auto chordComponentFaceId : linkedChordComponentFaceIds) {
        const auto& componentFace = baseComplex.componentF.at(chordComponentFaceId);
        if (componentFace.N_Cids.size() == 1) {
            startComponentFaceId = chordComponentFaceId;
            isLoop = false;
            break;
        }
    }

    for (int i = 1; i < linkedChordComponentFaceIds.size(); ++i) {
        auto id1 = linkedChordComponentFaceIds.at(i - 1);
        auto id2 = linkedChordComponentFaceIds.at(i);
        for (auto componentCellId : chordComponentCellIds) {
            const auto& componentCell = baseComplex.componentC.at(componentCellId);
            size_t count = 0;
            for (auto componentFaceId : componentCell.Fids) {
                if (componentFaceId == id1 || componentFaceId == id2 ) ++count;
            }
            if (count == 2) res.push_back(componentCellId);
        }
    }
    if (isLoop && linkedChordComponentFaceIds.size() != 0 && linkedChordComponentFaceIds.size() != 2) {
        auto id1 = linkedChordComponentFaceIds.front();
        auto id2 = linkedChordComponentFaceIds.back();
        for (auto componentCellId : chordComponentCellIds) {
            const auto& componentCell = baseComplex.componentC.at(componentCellId);
            size_t count = 0;
            for (auto componentFaceId : componentCell.Fids) {
                if (componentFaceId == id1 || componentFaceId == id2 ) ++count;
            }
            if (count == 2) res.push_back(componentCellId);
        }
    }
    return res;
}

std::vector<size_t> BaseComplexChord::GetLinkedChordComponentFaceIds(const std::vector<size_t>& chordComponentFaceIds) const{
    std::vector<size_t> res;
    bool isLoop = true;
    size_t startComponentFaceId = chordComponentFaceIds.at(0);
    for (auto chordComponentFaceId : chordComponentFaceIds) {
        const auto& componentFace = baseComplex.componentF.at(chordComponentFaceId);
        if (componentFace.N_Cids.size() == 1) {
            startComponentFaceId = chordComponentFaceId;
            isLoop = false;
            break;
        }
    }

    std::unordered_map<size_t, bool> ids_visited;
    for (auto id : chordComponentFaceIds)
        ids_visited[id] = false;
    std::stack<size_t> st;
    st.push(startComponentFaceId);
    //ids_visited[startComponentFaceId] = true;
    while (!st.empty()) {
        auto id = st.top();
        st.pop();
        ids_visited[id] = true;
        res.push_back(id);
        const auto& componentFace = baseComplex.componentF.at(id);
        const auto parallelComponentFaceIds = GetParallelComponentFaceIds(componentFace);
        for (auto next : parallelComponentFaceIds)
            if (!ids_visited[next]) st.push(next);
    }
    return res;
}

glm::dvec3 BaseComplexChord::GetCenter(const ComponentFace& c) const {
    glm::dvec3 sum;
    for (auto vid : c.Vids)
        sum += baseComplex.V.at(vid).xyz();
    return sum / double(c.Vids.size());
}

glm::dvec3 BaseComplexChord::GetCenter(const ComponentCell& c) const {
    glm::dvec3 sum;
    for (auto vid : c.Vids)
        sum += baseComplex.V.at(vid).xyz();
    return sum / double(c.Vids.size());
}
std::vector<glm::dvec3> BaseComplexChord::GetLinkedChordCurveVertices(const std::vector<size_t>& linkedChordComponentFaceIds, const std::vector<size_t>& linkedChordComponentCellIds) const {
    std::vector<glm::dvec3> centers;
    for (auto i = 0; i < linkedChordComponentCellIds.size(); ++i) {
        centers.push_back(GetCenter(baseComplex.componentF.at(linkedChordComponentFaceIds.at(i))));
        centers.push_back(GetCenter(baseComplex.componentC.at(linkedChordComponentCellIds.at(i))));
    }
    if (linkedChordComponentCellIds.size() == linkedChordComponentFaceIds.size()) centers.push_back(GetCenter(baseComplex.componentF.at(linkedChordComponentFaceIds.at(0))));
    else centers.push_back(GetCenter(baseComplex.componentF.at(linkedChordComponentFaceIds.back())));
    return centers;
}
