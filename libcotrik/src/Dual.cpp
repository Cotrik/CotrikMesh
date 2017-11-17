/*
 * Dual.cpp
 *
 *  Created on: Oct 30, 2017
 *      Author: cotrik
 */

#include "Dual.h"
#include <unordered_map>
#include <unordered_set>

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
        auto v = DualVertex((v1 + v2) * 0.5f);
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
