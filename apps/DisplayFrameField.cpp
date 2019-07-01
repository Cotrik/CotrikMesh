/*
 * DisplayFrameField.cpp
 *
 *  Created on: May 19, 2017
 *      Author: cotrik
 */


#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "ArgumentManager.h"

#include <iostream>

int main(int argc, char* argv[])
{
    if (argc < 3)
    {
        std::cout << "Usage: ExtractPatches tet.vtk quaternion.txt frames.vtk scale=<0.05>" << std::endl;
        return -1;
    }
    ArgumentManager argumentManager(argc, argv);
    std::string tetfilename = argv[1];
    std::string filename = argv[2];
    std::string output_filename = argv[3];
    double scale = 0.05;
    const std::string strscale = argumentManager.get("scale");
    if (!strscale.empty()) scale = std::stod(strscale);
    std::cout << "---------------------------------------" << std::endl;
    std::cout << "filename = " << filename << std::endl;
    std::cout << "strscale = " << strscale << std::endl;
//    std::cout << "min_numOfFaces_per_patch = " << min_numOfFaces_per_patch << std::endl;
//    std::cout << "localangle = " << strLocalAngle << std::endl;
//    std::cout << "globalangle = " << strGlobalaAngle << std::endl;
    std::cout << "---------------------------------------" << std::endl;


    // Read quaternions from file
    std::vector<glm::dquat> quaternions;
    quaternions.reserve(50000);
    {
        std::ifstream ifs(filename.c_str());
        glm::dquat q;
        while (ifs >> q.x >> q.y >> q.z >> q.w)
            quaternions.push_back(q);
        quaternions.resize(quaternions.size());
    }

    // Construct Frame Field
    struct frame {
        glm::dvec3 x , y, z;
        frame() {}
        frame(const frame& rhs):x(rhs.x), y(rhs.y), z(rhs.z) {}
        ~frame() {}
    };
    std::vector<frame> frames(quaternions.size(), frame());
    {
        const glm::dvec3 x(1.0, 0.0, 0.0);
        const glm::dvec3 y(0.0, 1.0, 0.0);
        const glm::dvec3 z(0.0, 0.0, 1.0);

        int i = 0;
        for (const auto& q : quaternions) {
            frames[i].x = glm::rotate(q, x);
            frames[i].y = glm::rotate(q, y);
            frames[i].z = glm::rotate(q, z);
            i++;
        }
    }

    MeshFileReader reader(tetfilename.c_str());
    Mesh& mesh = (Mesh&)reader.GetMesh();

    const double s = scale;
    for (size_t i = 0; i < mesh.V.size(); ++i) {
        const Vertex& v = mesh.V.at(i);
        glm::dvec3& x = frames[i].x;
        glm::dvec3& y = frames[i].y;
        glm::dvec3& z = frames[i].z;
        x = glm::dvec3(x.x * s, x.y * s, x.z * s) + v.xyz();
        y = glm::dvec3(y.x * s, y.y * s, y.z * s) + v.xyz();
        z = glm::dvec3(z.x * s, z.y * s, z.z * s) + v.xyz();
    }
// Make Edge frames
    std::vector<Vertex> V(mesh.V.size() * 4, Vertex());
    std::vector<Edge> E(mesh.V.size() * 3, Edge());
    for (size_t i = 0; i < mesh.V.size(); ++i) {
        const Vertex& v = mesh.V.at(i);
        V[i] = mesh.V.at(i).xyz();
        V[i].id = i;
    }
    for (size_t i = 0; i < mesh.V.size(); ++i) {
        const Vertex& v = mesh.V.at(i);
        V[i] = mesh.V.at(i).xyz();
        V[i].id = i;
    }
    for (size_t i = 0; i < frames.size(); ++i) {
        V[frames.size() + i] = frames[i].x;
        V[frames.size() + i].id = frames.size() + i;
    }
    for (size_t i = 0; i < frames.size(); ++i) {
        V[2*frames.size() + i] = frames[i].y;
        V[2*frames.size() + i].id = 2*frames.size() + i;
    }
    for (size_t i = 0; i < frames.size(); ++i) {
        V[3*frames.size() + i] = frames[i].z;
        V[3*frames.size() + i].id = 3*frames.size() + i;
    }
    for (size_t i = 0; i < mesh.V.size(); ++i) {
        const Vertex& v = mesh.V.at(i);
        Edge& edge1 = E.at(i);
        Edge& edge2 = E.at(i + mesh.V.size());
        Edge& edge3 = E.at(i + mesh.V.size() * 2);
        edge1.Vids = std::vector<size_t> {i, i + mesh.V.size()};
        edge2.Vids = std::vector<size_t> {i, i + mesh.V.size() * 2};
        edge3.Vids = std::vector<size_t> {i, i + mesh.V.size() * 3};
    }

    Mesh newMesh;
    newMesh.V = V;
    newMesh.E = E;
    std::vector<size_t> Eids(3 * mesh.V.size(), 0);
    for (size_t i = 0; i < Eids.size(); ++i)
        Eids[i] = i;
    MeshFileWriter writer(newMesh, output_filename.c_str());
    writer.WriteFramesVtk();

// Make Cube frames
    // Construct Frame Field
    struct v8 {
        std::vector<glm::dvec3> v;
        v8() {v.resize(8);}
        v8(const v8& rhs):v(rhs.v) {}
        ~v8() {}
    };
    std::vector<v8> v8s(quaternions.size(), v8());
    {
        const glm::dvec3 v0(-0.5, -0.5, +0.5);
        const glm::dvec3 v1(+0.5, -0.5, +0.5);
        const glm::dvec3 v2(+0.5, -0.5, -0.5);
        const glm::dvec3 v3(-0.5, -0.5, -0.5);
        const glm::dvec3 v4(-0.5, +0.5, +0.5);
        const glm::dvec3 v5(+0.5, +0.5, +0.5);
        const glm::dvec3 v6(+0.5, +0.5, -0.5);
        const glm::dvec3 v7(-0.5, +0.5, -0.5);

        int i = 0;
        for (const auto& q : quaternions) {
            v8s[i].v[0] = glm::rotate(q, v0);
            v8s[i].v[1] = glm::rotate(q, v1);
            v8s[i].v[2] = glm::rotate(q, v2);
            v8s[i].v[3] = glm::rotate(q, v3);
            v8s[i].v[4] = glm::rotate(q, v4);
            v8s[i].v[5] = glm::rotate(q, v5);
            v8s[i].v[6] = glm::rotate(q, v6);
            v8s[i].v[7] = glm::rotate(q, v7);
            i++;
        }
    }

    //std::vector<Vertex> V(mesh.V.size() * 8, Vertex());
    V.resize(mesh.V.size() * 8);
    std::vector<Cell> C(mesh.V.size() * 1, Cell());
    double fs = s;
    for (size_t i = 0; i < mesh.V.size(); ++i) {
        const Vertex& v = mesh.V.at(i);
        for (size_t j = 0; j < 8; ++j) {
            V[i * 8 + j] = v8s[i].v[j] * fs + v.xyz();
            V[i * 8 + j].id = i * 8 + j;
        }
    }
    for (size_t i = 0; i < mesh.V.size(); ++i) {
        Cell& cell = C.at(i);
        const size_t o = 8*i;
        cell.Vids = std::vector<size_t> {o + 0, o + 1, o + 2, o + 3, o + 4, o + 5, o + 6, o + 7};
        cell.id = i;
    }
//////////////////////////////////////////////////////////////
    //
//////////////////////////////////////////////////////////////
    MeshFileWriter cubeMeshWriter(V, C, "CUBE.vtk");
    cubeMeshWriter.WriteFile();

    return 0;
}
