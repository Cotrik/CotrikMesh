/*
 * RefineDualQuad.cpp
 *
 *  Created on: Jan 12, 2018
 *      Author: cotrik
 */

#include "RefinedDualQuad.h"
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <iostream>
#include <algorithm>

static std::unordered_map<unsigned long long, size_t> get_key_edgeId(const Mesh& mesh) {
    std::unordered_map<unsigned long long, size_t> key_edgeId;
//#pragma omp parallel for
    for (size_t i = 0; i < mesh.E.size(); ++i) {
        const auto& e = mesh.E.at(i);
        key_edgeId[(e.Vids[0] << 32) | e.Vids[1]] = i;
        key_edgeId[(e.Vids[1] << 32) | e.Vids[0]] = i;
    }
    return key_edgeId;
}

static std::unordered_map<std::string, size_t> get_key_faceId(const Mesh& mesh) {
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

Mesh GetRefineQuadMesh(const Mesh& quad_mesh, int clockwise)
{
    const Mesh& new_mesh = quad_mesh;
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

    /////////////////////////////////////////////////////////////////
    /*
     * 0       7        3
     |---------|---------|
     |         |         |
     |         |         |
     |         |         |
    4|---------8---------|6
     |         |         |
     |         |         |
     |         |         |
     |---------|---------|
      1        5        2
     * */
    // add cells
    const unsigned int QuadEdge[4][2] = {{ 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 }};
    const int QuadRefine[4][4] = {{0, 4, 8, 7}, {4, 1, 5, 8}, {8, 5, 2, 6}, {7, 8, 6, 3}};
    auto key_edgeId = get_key_edgeId(new_mesh);
    auto key_faceId = get_key_faceId(new_mesh);
    Face face(4);
    std::vector<Face> new_faces(4 * new_mesh.F.size(), face);
    int count = 0;
    for (const auto& f : new_mesh.F) {
        unsigned long v_index[9];
        for (auto j = 0; j < 4; j++)
            v_index[j] = f.Vids.at(j);
        if (clockwise != 0) std::swap(v_index[1], v_index[3]);
        for (unsigned long j = 0; j < 4; j++) {
            const Edge e({f.Vids.at(QuadEdge[j][0]), f.Vids.at(QuadEdge[j][1])});
            unsigned long long key = (e.Vids[0] << 32) | e.Vids[1];
            if (key_edgeId.find(key) == key_edgeId.end()) std::cout << "Edge search Error !" << std::endl;
            auto e_index = key_edgeId[key];
            v_index[4 + j] = new_mesh.V.size() + e_index;
        }
        v_index[8] = new_mesh.V.size() + new_mesh.E.size() + f.id;
        for (int k = 0; k < 4; k++, count++)
            for (int j = 0; j < 4; j++)
                new_faces[count].Vids[j] = v_index[QuadRefine[k][j]];
    }
    Mesh mesh(new_vertex, new_faces, QUAD);
    return mesh;
}

RefinedDualQuad::RefinedDualQuad(Mesh& mesh)
: mesh(mesh)
{
    // TODO Auto-generated constructor stub

}

RefinedDualQuad::~RefinedDualQuad() {
    // TODO Auto-generated destructor stub
}

void RefinedDualQuad::Build() {
    //refinedMesh = GetRefineMesh(mesh);
    BuildV();
    BuildE();
    BuildF();
    BuildC();
}

void RefinedDualQuad::BuildV() {
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
    size_t id = 0;
    for (auto& v : V)
        v.id = id++;
}

void RefinedDualQuad::BuildE() {
    E.resize(2 * mesh.E.size() + 4 * mesh.F.size());
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
        const Vertex& mid = V.at(mesh.V.size() + mesh.E.size() + i);
        for (size_t j = 0; j < f.Eids.size(); j++) {
            auto eid = f.Eids.at(j);
            E.at(offset + i * f.Eids.size() + j).Vids = {mesh.V.size() + eid, mid.id};
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

static size_t getParallelEdgeIds(const Mesh& mesh, const Face& face, const Edge& edge) {
    for (int i = 1; i < face.Eids.size(); ++i) {
        const auto& e = mesh.E.at(face.Eids.at(i));
        if (e.Vids[0] != edge.Vids[0] && e.Vids[1] != edge.Vids[1] && e.Vids[0] != edge.Vids[1] && e.Vids[1] != edge.Vids[0]) return e.id;
    }
    return MAXID;
}
static size_t getCrossEdgeIds(const Mesh& mesh, const Face& face, const Edge& edge, const Vertex& v) {
    for (auto eid : face.Eids) {
        const auto& e = mesh.E.at(eid);
        if ((e.Vids[0] == v.id || e.Vids[1] == v.id) && eid != edge.id) return eid;
    }
    return MAXID;
}
void RefinedDualQuad::BuildF() {
    DualFace ff({MAXID, MAXID, MAXID, MAXID});
    F.resize(4 * mesh.F.size(), ff);
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
    size_t id = 0;
    for (auto& f : F)
        f.id = id++;
}

void RefinedDualQuad::BuildC() {

}
