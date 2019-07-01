/*
 * QuadTriMeshRefine.cpp
 *
 *  Created on: Nov 6, 2018
 *      Author: cotrik
 */
#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "ArgumentManager.h"

#include <unordered_map>
#include <unordered_set>
#include <set>
#include <iostream>
#include <algorithm>

Mesh Refine(const Mesh& quad_mesh, int clockwise);

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: QuadTriMeshRefine QuadTri.vtk refine.vtk" << std::endl;
        return -1;
    }
    ArgumentManager argumentManager(argc, argv);
    std::string input = argv[1];
    std::string output = argv[2];
    std::cout << "---------------------------------------" << std::endl;
    std::cout << "input  = " << input << std::endl;
    std::cout << "output = " << output << std::endl;
    std::cout << "---------------------------------------" << std::endl;

    MeshFileReader reader(input.c_str());
    Mesh& mesh = (Mesh&) reader.GetMesh();
    mesh.BuildAllConnectivities();
//    mesh.ExtractBoundary();
//    mesh.ExtractSingularities();
//    mesh.BuildParallelE();
//    mesh.BuildConsecutiveE();
//    mesh.BuildOrthogonalE();

    auto refinedMesh = Refine(mesh, 0);
    MeshFileWriter writer(refinedMesh, output.c_str());
    writer.WriteFile();
}

static std::unordered_map<unsigned long long, size_t> get_key_edgeId(const Mesh& mesh) {
    std::unordered_map<unsigned long long, size_t> key_edgeId;
    for (size_t i = 0; i < mesh.E.size(); ++i) {
        const auto& e = mesh.E.at(i);
        key_edgeId[(e.Vids[0] << 32) | e.Vids[1]] = i;
        key_edgeId[(e.Vids[1] << 32) | e.Vids[0]] = i;
    }
    return key_edgeId;
}

static std::unordered_map<std::string, size_t> get_key_faceId(const Mesh& mesh) {
    std::unordered_map<std::string, size_t> key_faceId;
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

const int QuadRefine[4][4] = {
    0, 4, 8, 7,
    1, 5, 8, 4,
    2, 6, 8, 5,
    3, 7, 8, 6
};

const int TriRefine[3][4] = {
    0, 3, 6, 5,
    1, 4, 6, 3,
    2, 5, 6, 4
};

Mesh Refine(const Mesh& hex_mesh, int clockwise) {
    const Mesh& new_mesh = hex_mesh;
    ////////////////////////////////////////////////////////////////////////////
    // add vertices
    std::vector<Vertex> new_vertex(new_mesh.V.size() + new_mesh.E.size() + new_mesh.C.size());
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
    size_t numOfTri = 0, numOfQuad = 0;
    for (size_t i = 0; i < new_mesh.C.size(); i++) {
        const auto& f = new_mesh.C.at(i);
        if (f.Vids.size() == 4) {
            const Vertex& v0 = new_mesh.V.at(f.Vids[0]);
            const Vertex& v1 = new_mesh.V.at(f.Vids[1]);
            const Vertex& v2 = new_mesh.V.at(f.Vids[2]);
            const Vertex& v3 = new_mesh.V.at(f.Vids[3]);
            new_vertex.at(offset + i) = 0.25 * (v0.xyz() + v1.xyz() + v2.xyz() + v3.xyz());
            ++numOfQuad;
        } else  if (f.Vids.size() == 3) {
            const Vertex& v0 = new_mesh.V.at(f.Vids[0]);
            const Vertex& v1 = new_mesh.V.at(f.Vids[1]);
            const Vertex& v2 = new_mesh.V.at(f.Vids[2]);
            new_vertex.at(offset + i) = 0.3333333 * (v0.xyz() + v1.xyz() + v2.xyz());
            ++numOfTri;
        }
    }
    auto key_edgeId = get_key_edgeId(new_mesh);
    //auto key_faceId = get_key_faceId(new_mesh);
    Cell cell(4);
    std::vector<Cell> new_cells(numOfTri * 3 + numOfQuad * 4, cell);
    int count = 0;
    for (size_t i = 0; i < new_mesh.C.size(); i++) {
        unsigned long v_index[9];
        const auto & f = new_mesh.C.at(i);
        for (auto j = 0; j < f.Vids.size(); j++)
            v_index[j] = f.Vids.at(j);
//        if (clockwise != 0) {
//            std::swap(v_index[1], v_index[3]);
//            std::swap(v_index[5], v_index[7]);
//        }
        if (f.Vids.size() == 4) {
            for (unsigned long j = 0; j < 4; j++) {
                const Edge e( { f.Vids.at(QuadEdge[j][0]), f.Vids.at(QuadEdge[j][1]) });
                unsigned long long key = (e.Vids[0] << 32) | e.Vids[1];
                if (key_edgeId.find(key) == key_edgeId.end()) std::cout << "Edge search Error !" << std::endl;
                auto e_index = key_edgeId[key];
                v_index[4 + j] = new_mesh.V.size() + e_index;
            }
            v_index[8] = new_mesh.V.size() + new_mesh.E.size() + i;
            for (int k = 0; k < 4; k++, count++)
                for (int j = 0; j < 4; j++)
                    new_cells[count].Vids[j] = v_index[QuadRefine[k][j]];
        } else if (f.Vids.size() == 3) {
            for (unsigned long j = 0; j < 3; j++) {
                const Edge e( { f.Vids.at(TriEdge[j][0]), f.Vids.at(TriEdge[j][1]) });
                unsigned long long key = (e.Vids[0] << 32) | e.Vids[1];
                if (key_edgeId.find(key) == key_edgeId.end()) std::cout << "Edge search Error !" << std::endl;
                auto e_index = key_edgeId[key];
                v_index[3 + j] = new_mesh.V.size() + e_index;
            }
            v_index[6] = new_mesh.V.size() + new_mesh.E.size() + i;
            for (int k = 0; k < 3; k++, count++)
                for (int j = 0; j < 4; j++)
                    new_cells[count].Vids[j] = v_index[TriRefine[k][j]];
        }
    }
    Mesh mesh(new_vertex, new_cells, QUAD);
    return mesh;
}
