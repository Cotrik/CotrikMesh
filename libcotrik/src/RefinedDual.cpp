/*
 * RefineDual.cpp
 *
 *  Created on: Nov 21, 2017
 *      Author: cotrik
 */

#include "RefinedDual.h"
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <iostream>
#include <algorithm>

Mesh GetRefineMesh2(const Mesh& hex_mesh, int clockwise)
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
        const Vertex& v1 = new_mesh.V.at(f.Vids[1]);
        const Vertex& v2 = new_mesh.V.at(f.Vids[2]);
        const Vertex& v3 = new_mesh.V.at(f.Vids[3]);
        new_vertex.at(offset + i) = 0.25f * (v0.xyz() + v1.xyz() + v2.xyz() + v3.xyz());
    }
    offset = new_mesh.V.size() + new_mesh.E.size() + new_mesh.F.size();
    for (size_t i = 0; i < new_mesh.C.size(); i++) {
        const Cell& c = new_mesh.C.at(i);
        const Vertex& v0 = new_mesh.V.at(c.Vids[0]);
        const Vertex& v1 = new_mesh.V.at(c.Vids[1]);
        const Vertex& v2 = new_mesh.V.at(c.Vids[2]);
        const Vertex& v3 = new_mesh.V.at(c.Vids[3]);
        const Vertex& v4 = new_mesh.V.at(c.Vids[4]);
        const Vertex& v5 = new_mesh.V.at(c.Vids[5]);
        const Vertex& v6 = new_mesh.V.at(c.Vids[6]);
        const Vertex& v7 = new_mesh.V.at(c.Vids[7]);
        new_vertex.at(offset + i) = 0.125f * (v0.xyz() + v1.xyz() + v2.xyz() + v3.xyz() + v4.xyz() + v5.xyz() + v6.xyz() + v7.xyz());
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
std::unordered_map<unsigned long long, size_t> get_key_edgeId(const Mesh& mesh) {
    std::unordered_map<unsigned long long, size_t> key_edgeId;
//#pragma omp parallel for
    for (size_t i = 0; i < mesh.E.size(); ++i) {
        const auto& e = mesh.E.at(i);
        key_edgeId[(e.Vids[0] << 32) | e.Vids[1]] = i;
        key_edgeId[(e.Vids[1] << 32) | e.Vids[0]] = i;
    }
    return key_edgeId;
}

std::unordered_map<std::string, size_t> get_key_faceId(const Mesh& mesh) {
    std::unordered_map<std::string, size_t> key_faceId;
//#pragma omp parallel for
    for (size_t i = 0; i < mesh.F.size(); ++i) {
        const auto& f = mesh.F.at(i);
        std::set<size_t> vids(f.Vids.begin(), f.Vids.end());
        std::string s;
        for (auto vid : vids)
            s += std::to_string(vid) + "@";
        key_faceId[s] = i;
    }
    return key_faceId;
}
Mesh GetRefineMesh3(const Mesh& hex_mesh, int clockwise)
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
        const Vertex& v1 = new_mesh.V.at(f.Vids[1]);
        const Vertex& v2 = new_mesh.V.at(f.Vids[2]);
        const Vertex& v3 = new_mesh.V.at(f.Vids[3]);
        new_vertex.at(offset + i) = 0.25f * (v0.xyz() + v1.xyz() + v2.xyz() + v3.xyz());
    }
    offset = new_mesh.V.size() + new_mesh.E.size() + new_mesh.F.size();
    for (size_t i = 0; i < new_mesh.C.size(); i++) {
        const Cell& c = new_mesh.C.at(i);
        const Vertex& v0 = new_mesh.V.at(c.Vids[0]);
        const Vertex& v1 = new_mesh.V.at(c.Vids[1]);
        const Vertex& v2 = new_mesh.V.at(c.Vids[2]);
        const Vertex& v3 = new_mesh.V.at(c.Vids[3]);
        const Vertex& v4 = new_mesh.V.at(c.Vids[4]);
        const Vertex& v5 = new_mesh.V.at(c.Vids[5]);
        const Vertex& v6 = new_mesh.V.at(c.Vids[6]);
        const Vertex& v7 = new_mesh.V.at(c.Vids[7]);
        new_vertex.at(offset + i) = 0.125f * (v0.xyz() + v1.xyz() + v2.xyz() + v3.xyz() + v4.xyz() + v5.xyz() + v6.xyz() + v7.xyz());
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

    auto key_edgeId = get_key_edgeId(new_mesh);
    auto key_faceId = get_key_faceId(new_mesh);
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
            unsigned long long key = (e.Vids[0] << 32) | e.Vids[1];
            if (key_edgeId.find(key) == key_edgeId.end()) std::cout << "Edge search Error !" << std::endl;
            auto e_index = key_edgeId[key];
            v_index[8 + j] = new_mesh.V.size() + e_index;
        }
//#pragma omp parallel for
        for (unsigned long j = 0; j < 6; j++) {
            std::set<size_t> vids({c.Vids.at(HexFaces[j][0]), c.Vids.at(HexFaces[j][1]), c.Vids.at(HexFaces[j][2]), c.Vids.at(HexFaces[j][3])});
            std::string s;
            for (auto vid : vids)
                s += std::to_string(vid) + "@";
            v_index[20 + j] = new_mesh.V.size() + new_mesh.E.size() + key_faceId[s];
        }
        v_index[26] = new_mesh.V.size() + new_mesh.E.size() + new_mesh.F.size() + i;
        for (int k = 0; k < 8; k++, count++)
            for (int j = 0; j < 8; j++)
                new_cells[count].Vids[j] = v_index[HexRefine[k][j]];
    }
    Mesh mesh(new_vertex, new_cells, HEXAHEDRA);
    return mesh;
}
RefinedDual::RefinedDual(Mesh& mesh)
: mesh(mesh)
{
    // TODO Auto-generated constructor stub

}

RefinedDual::~RefinedDual() {
    // TODO Auto-generated destructor stub
}

void RefinedDual::Build() {
    //refinedMesh = GetRefineMesh(mesh);
    BuildV();
    BuildE();
    BuildF();
    BuildC();
}

void RefinedDual::BuildV() {
    V.resize(mesh.V.size() + mesh.E.size() + mesh.F.size() + mesh.C.size());
    for (size_t i = 0; i < mesh.V.size(); i++)
        V.at(i) = mesh.V.at(i).xyz();
    size_t offset = mesh.V.size();
    for (size_t i = 0; i < mesh.E.size(); i++) {
        const Edge& e = mesh.E.at(i);
        const Vertex& v0 = mesh.V.at(e.Vids[0]);
        const Vertex& v1 = mesh.V.at(e.Vids[1]);
        V.at(offset + i) = 0.5f * (v0.xyz() + v1.xyz());
    }
    offset = mesh.V.size() + mesh.E.size();
    for (size_t i = 0; i < mesh.F.size(); i++) {
        const Face& f = mesh.F.at(i);
        const Vertex& v0 = mesh.V.at(f.Vids[0]);
        const Vertex& v1 = mesh.V.at(f.Vids[1]);
        const Vertex& v2 = mesh.V.at(f.Vids[2]);
        const Vertex& v3 = mesh.V.at(f.Vids[3]);
        V.at(offset + i) = 0.25f * (v0.xyz() + v1.xyz() + v2.xyz() + v3.xyz());
    }
    offset = mesh.V.size() + mesh.E.size() + mesh.F.size();
    for (size_t i = 0; i < mesh.C.size(); i++) {
        const Cell& c = mesh.C.at(i);
        const Vertex& v0 = mesh.V.at(c.Vids[0]);
        const Vertex& v1 = mesh.V.at(c.Vids[1]);
        const Vertex& v2 = mesh.V.at(c.Vids[2]);
        const Vertex& v3 = mesh.V.at(c.Vids[3]);
        const Vertex& v4 = mesh.V.at(c.Vids[4]);
        const Vertex& v5 = mesh.V.at(c.Vids[5]);
        const Vertex& v6 = mesh.V.at(c.Vids[6]);
        const Vertex& v7 = mesh.V.at(c.Vids[7]);
        V.at(offset + i) = 0.125f * (v0.xyz() + v1.xyz() + v2.xyz() + v3.xyz() + v4.xyz() + v5.xyz() + v6.xyz() + v7.xyz());
    }
    size_t id = 0;
    for (auto& v : V)
        v.id = id++;
}

void RefinedDual::BuildE() {
    E.resize(2 * mesh.E.size() + 4 * mesh.F.size() + 6 * mesh.C.size());
    size_t offset = 0;
    for (size_t i = 0; i < mesh.E.size(); i++) {
        const Edge& e = mesh.E.at(i);
        const Vertex& v0 = V.at(e.Vids[0]);
        const Vertex& v1 = V.at(e.Vids[1]);
        const Vertex& mid = V.at(mesh.V.size() + i);
        E.at(offset + i * 2 + 0).Vids = {v0.id, mid.id};
        E.at(offset + i * 2 + 1).Vids = {v1.id, mid.id};
    }
    offset = 2 * mesh.E.size();
    for (size_t i = 0; i < mesh.F.size(); i++) {
        const Face& f = mesh.F.at(i);
        for (size_t j = 0; j < f.Eids.size(); j++) {
            auto eid = f.Eids.at(j);
            const Vertex& mid = V.at(mesh.V.size() + mesh.E.size() + i);
            E.at(offset + i * f.Eids.size() + j).Vids = {mesh.V.size() + eid, mid.id};
        }
    }
    offset = 2 * mesh.E.size() + 4 * mesh.F.size();
    for (size_t i = 0; i < mesh.C.size(); i++) {
        const Cell& c = mesh.C.at(i);
        for (size_t j = 0; j < c.Fids.size(); j++) {
            auto fid = c.Fids.at(j);
            const Vertex& mid = V.at(mesh.V.size() + mesh.E.size() + mesh.F.size() + i);
            E.at(offset + i * c.Fids.size() + j).Vids = {mesh.V.size() + mesh.E.size() + fid, mid.id};
        }
    }
    size_t id = 0;
    for (auto& e : E)
        e.id = id++;
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

size_t getParallelEdgeIds(const Mesh& mesh, const Face& face, const Edge& edge) {
    for (int i = 1; i < face.Eids.size(); ++i) {
        const auto& e = mesh.E.at(face.Eids.at(i));
        if (e.Vids[0] != edge.Vids[0] && e.Vids[1] != edge.Vids[1] && e.Vids[0] != edge.Vids[1] && e.Vids[1] != edge.Vids[0]) return e.id;
    }
    return MAXID;
}
size_t getCrossEdgeIds(const Mesh& mesh, const Face& face, const Edge& edge, const Vertex& v) {
    for (auto eid : face.Eids) {
        const auto& e = mesh.E.at(eid);
        if ((e.Vids[0] == v.id || e.Vids[1] == v.id) && eid != edge.id) return eid;
    }
    return MAXID;
}
void RefinedDual::BuildF() {
    DualFace ff({MAXID, MAXID, MAXID, MAXID});
    F.resize(4 * mesh.F.size() + 12 * mesh.C.size(), ff);
    size_t offset = 0;
    for (size_t i = 0; i < mesh.F.size(); i++) {
        const auto& f = mesh.F.at(i);
        const Edge& e0 = mesh.E.at(f.Eids.at(0));
        const Edge& e1 = mesh.E.at(getParallelEdgeIds(mesh, f, e0));
        const Edge& e0v0 = getCrossEdgeIds(mesh, f, e0, mesh.V.at(e0.Vids[0]));
        const Edge& e0v1 = getCrossEdgeIds(mesh, f, e0, mesh.V.at(e0.Vids[1]));
        const Edge& e1v0 = getCrossEdgeIds(mesh, f, e1, mesh.V.at(e1.Vids[0]));
        const Edge& e1v1 = getCrossEdgeIds(mesh, f, e1, mesh.V.at(e1.Vids[1]));

        const Vertex& mid = V.at(mesh.V.size() + mesh.E.size() + i);
        F.at(offset + i * 4 + 0).Vids = {e0.Vids[0], mesh.V.size() + e0.id, mid.id, mesh.V.size() + e0v0.id};
        F.at(offset + i * 4 + 1).Vids = {e0.Vids[1], mesh.V.size() + e0.id, mid.id, mesh.V.size() + e0v1.id};
        F.at(offset + i * 4 + 2).Vids = {e1.Vids[0], mesh.V.size() + e1.id, mid.id, mesh.V.size() + e1v0.id};
        F.at(offset + i * 4 + 3).Vids = {e1.Vids[1], mesh.V.size() + e1.id, mid.id, mesh.V.size() + e1v1.id};
    }
    const unsigned int HexEdge[12][2] = { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 }, { 4, 5 }, { 5, 6 }, { 6, 7 }, { 7, 4 }, { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 }, };

    const int HexRefine[8][8] = {
        11, 20, 10, 3, 22, 26, 25, 19,
        20, 9, 2, 10, 26, 23, 18, 25,
        22, 26, 25, 19, 15, 21, 14, 7,
        26, 23, 18, 25, 21, 13, 6, 14,
        0, 8, 20, 11, 16, 24, 26, 22,
        8, 1, 9, 20, 24, 17, 23, 26,
        16, 24, 26, 22, 4, 12, 21, 15,
        24, 17, 23, 26, 12, 5, 13, 21
    };
//    const int HexRefineFace[12][4] =  {
//        8, 20, 26, 24,
//        9, 20, 26, 23,
//        10, 20, 26, 25,
//        11, 20, 26, 22,
//        12, 21, 26, 24,
//        13, 21, 26, 23,
//        14, 21, 26, 25,
//        15, 21, 26, 22,
//        16, 24, 26, 22,
//        17, 24, 26, 23,
//        18, 25, 26, 23,
//        19, 25, 26, 22
//    };
    const int HexRefineFace[12][4] = {
        {16, 24, 26, 22},
        {17, 24, 26, 23},
        {18, 25, 26, 23},
        {19, 25, 26, 22},

        {8, 20, 26, 24},
        {10, 20, 26, 25},
        {12, 21, 26, 24},
        {14, 21, 26, 25},

        {9, 20, 26, 23},
        {11, 20, 26, 22},
        {13, 21, 26, 23},
        {15, 21, 26, 22}
    };

    auto key_edgeId = get_key_edgeId(mesh);
    auto key_faceId = get_key_faceId(mesh);

    int clockwise = 0;
    size_t count = 0;
    offset = 4 * mesh.F.size();
    for (size_t i = 0; i < mesh.C.size(); i++) {
        unsigned long v_index[27];
        const Cell& c = mesh.C.at(i);
        for (auto j = 0; j < 8; j++)
            v_index[j] = c.Vids.at(j);
        if (clockwise != 0) {
            std::swap(v_index[1], v_index[3]);
            std::swap(v_index[5], v_index[7]);
        }
        for (unsigned long j = 0; j < 12; j++) {
            const Edge e({c.Vids.at(HexEdge[j][0]), c.Vids.at(HexEdge[j][1])});
            unsigned long long key = (e.Vids[0] << 32) | e.Vids[1];
            if (key_edgeId.find(key) == key_edgeId.end()) std::cout << "Edge search Error !" << std::endl;
            auto e_index = key_edgeId[key];
            v_index[8 + j] = mesh.V.size() + e_index;
        }
        for (unsigned long j = 0; j < 6; j++) {
            std::set<size_t> vids({c.Vids.at(HexFaces[j][0]), c.Vids.at(HexFaces[j][1]), c.Vids.at(HexFaces[j][2]), c.Vids.at(HexFaces[j][3])});
            std::string s;
            for (auto vid : vids)
                s += std::to_string(vid) + "@";
            v_index[20 + j] = mesh.V.size() + mesh.E.size() + key_faceId[s];
        }
        v_index[26] = mesh.V.size() + mesh.E.size() + mesh.F.size() + i;
        for (int k = 0; k < 12; k++, count++)
            for (int j = 0; j < 4; j++)
                F[offset + count].Vids[j] = v_index[HexRefineFace[k][j]];
    }
    size_t id = 0;
    for (auto& f : F)
        f.id = id++;
}

void RefinedDual::BuildC() {

}
