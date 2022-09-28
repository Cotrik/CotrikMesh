/*
 * Dual.cpp
 *
 *  Created on: Oct 30, 2017
 *      Author: cotrik
 */

#include "Dual.h"
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <algorithm>

Mesh GetRefineMesh1(const Mesh& hex_mesh, int clockwise = 0)
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
        new_vertex.at(offset + i) = 0.5 * (v0.xyz() + v1.xyz());
    }
    offset = new_mesh.V.size() + new_mesh.E.size();
    for (size_t i = 0; i < new_mesh.F.size(); i++) {
        const Face& f = new_mesh.F.at(i);
        const Vertex& v0 = new_mesh.V.at(f.Vids[0]);
        const Vertex& v1 = new_mesh.V.at(f.Vids[2]);
        new_vertex.at(offset + i) = 0.5 * (v0.xyz() + v1.xyz());
    }
    offset = new_mesh.V.size() + new_mesh.E.size() + new_mesh.F.size();
    for (size_t i = 0; i < new_mesh.C.size(); i++) {
        const Cell& c = new_mesh.C.at(i);
        const Vertex& v0 = new_mesh.V.at(c.Vids[0]);
        const Vertex& v1 = new_mesh.V.at(c.Vids[6]);
        new_vertex.at(offset + i) = 0.5 * (v0.xyz() + v1.xyz());
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
        }
        v_index[26] = new_mesh.V.size() + new_mesh.E.size() + new_mesh.F.size() + i;
        for (int k = 0; k < 8; k++, count++)
            for (int j = 0; j < 8; j++)
                new_cells[count].Vids[j] = v_index[HexRefine[k][j]];
    }
    Mesh mesh(new_vertex, new_cells, HEXAHEDRA);
    return mesh;
}

Dual::Dual(Mesh& mesh)
: mesh(mesh)
{
    // TODO Auto-generated constructor stub

}

Dual::~Dual() {
    // TODO Auto-generated destructor stub
}

void Dual::Build() {
    BuildV();
    BuildE();
    BuildF();
    BuildC();
}

void Dual::BuildV() {
    V.resize(mesh.E.size());
    size_t id = 0;
    for (auto& e : mesh.E) {
        const auto& v1 = mesh.V.at(e.Vids[0]);
        const auto& v2 = mesh.V.at(e.Vids[1]);
        auto v = DualVertex((v1 + v2) * 0.5);
        v.id = id;
        V[id++] = v;
    }

//    size_t surfaceSize = 0;
//    size_t innerFaceSize = 0;
//    for (size_t i = 0; i < mesh.F.size(); i++)
//        if (mesh.F.at(i).isBoundary) surfaceSize++;
//        else innerFaceSize++;
//
//    V.resize(surfaceSize + mesh.C.size());
//    for (size_t i = 0; i < mesh.C.size(); i++)
//        V.at(i).BuildFrom(mesh, i, false);
//
//    size_t cellId = mesh.C.size();
//    for (size_t i = 0; i < mesh.F.size(); i++)
//        if (mesh.F.at(i).isBoundary) V.at(cellId).BuildFrom(mesh, cellId++, true, i);
}

void Dual::BuildE() {
    E.resize(2 * mesh.F.size());
    size_t id = 0;
    for (auto& f : mesh.F) {
        const auto& e1 = mesh.E.at(f.Eids[0]);
        const auto& e2 = mesh.E.at(f.Eids[1]);
        const auto& e3 = mesh.E.at(f.Eids[2]);
        const auto& e4 = mesh.E.at(f.Eids[3]);

        const auto& e1v1 = mesh.V.at(e1.Vids[0]);
        const auto& e1v2 = mesh.V.at(e1.Vids[1]);
        const auto& e2v1 = mesh.V.at(e2.Vids[0]);
        const auto& e2v2 = mesh.V.at(e2.Vids[1]);

        std::vector<size_t> vids1(2, MAXID);
        std::vector<size_t> vids2(2, MAXID);
        if (e1v1.id == e2v1.id || e1v2.id == e2v1.id || e1v1.id == e2v2.id || e1v2.id == e2v2.id) {
            vids1[0] = e1.id;
            vids1[1] = e3.id;
            vids2[0] = e2.id;
            vids2[1] = e4.id;
        } else {
            vids1[0] = e1.id;
            vids1[1] = e2.id;
            vids2[0] = e3.id;
            vids2[1] = e4.id;
        }
        auto de1 = DualEdge(vids1);
        de1.id = id;
        E[id++] = de1;

        auto de2 = DualEdge(vids2);
        de2.id = id;
        E[id++] = de2;
    }

//    E.resize(mesh.F.size());
//    size_t innerFaceId = 0;
//    size_t surfaceId = 0;
//    for (size_t i = 0; i < mesh.F.size(); i++) {
//        E.at(i).id = mesh.F.at(i).id;
//        if (!mesh.F.at(i).isBoundary) E.at(i).Vids = mesh.F.at(i).N_Cids;
//        else {
//            E.at(i).Vids.push_back(mesh.C.size() + surfaceId++);
//            E.at(i).Vids.push_back(mesh.F.at(i).N_Cids.at(0));
//        }
//    }
}

/*
                  3__________________2
                  /|                 /|
                 / |                / |
                /  |               /  |
            0  /___|_____________1/   |
               |   |              |   |
               |   |              |   |
               |   |              |   |
               |   |______________|___|
               |   / 7            |  /6
               |  /               | /
               | /                |/
               |/_________________/
             4                    5
*/
const size_t parallelEdges[3][4][2] = {
        {{0, 4}, {3, 7}, {2, 6}, {1, 5}},
        {{0, 1}, {3, 2}, {7, 6}, {4, 5}},
        {{0, 3}, {1, 2}, {5, 6}, {4, 7}}
};

/*

            0   __________________3
               |                  |
               |                  |
               |                  |
               |                  |
               |                  |
               |                  |
               |                  |
               |__________________|
             1                    2
*/

const size_t parallelQuadEdges[2][2][2] = {
        {{0, 1}, {3, 2}},
        {{0, 3}, {1, 2}}
};

std::vector<std::vector<size_t>> getParallelEdgeIds(const Mesh& mesh, const Face& face) {
    std::vector<std::vector<size_t>> parallel_edge_ids(2, std::vector<size_t>(2, MAXID));
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            auto v1 = face.Vids.at(parallelQuadEdges[i][j][0]);
            auto v2 = face.Vids.at(parallelQuadEdges[i][j][1]);
            for (const auto edge_id : face.Eids) {
                const auto& edge = mesh.E.at(edge_id);
                auto v_1 = edge.Vids[0];
                auto v_2 = edge.Vids[1];
                if ((v1 == v_1 && v2 == v_2) || (v1 == v_2 && v2 == v_1)) {
                    parallel_edge_ids[i][j] = edge_id;
                    break;
                }
            }
        }
    }
    return parallel_edge_ids;
}

std::vector<std::vector<size_t>> getParallelEdgeIds(const Mesh& mesh, const Cell& cell) {
    std::vector<std::vector<size_t>> parallel_edge_ids(3, std::vector<size_t>(4, MAXID));
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 4; ++j) {
            auto v1 = cell.Vids.at(parallelEdges[i][j][0]);
            auto v2 = cell.Vids.at(parallelEdges[i][j][1]);
            for (const auto edge_id : cell.Eids) {
                const auto& edge = mesh.E.at(edge_id);
                auto v_1 = edge.Vids[0];
                auto v_2 = edge.Vids[1];
                if ((v1 == v_1 && v2 == v_2) || (v1 == v_2 && v2 == v_1)) {
                    parallel_edge_ids[i][j] = edge_id;
                    break;
                }
            }
        }
    }
    return parallel_edge_ids;
}

void Dual::BuildF() {
    F.resize(3 * mesh.C.size());
    size_t id = 0;
    for (auto& c : mesh.C) {
        std::vector<std::vector<size_t>> parallel_edge_ids = getParallelEdgeIds(mesh, c);
        for (int i = 0; i < 3; ++i) {
            auto f = DualFace(parallel_edge_ids[i]);
            f.id = id;
            F[id++] = f;
        }
    }
}

void Dual::BuildC() {

}
