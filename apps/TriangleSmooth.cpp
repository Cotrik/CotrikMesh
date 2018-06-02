/*
 * TriangleSmooth.cpp
 *
 *  Created on: Dec 29, 2017
 *      Author: cotrik
 */

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include <iostream>
#include <fstream>

void Smooth(Mesh& mesh) {
    for (auto& v: mesh.V) {
        if (v.isBoundary) continue;
        glm::vec3 sum(0.0, 0.0, 0.0);
        for (auto vid : v.N_Vids)
            sum += mesh.V.at(vid).xyz();
        v = sum / float(v.N_Vids.size());
    }
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: TriangleSmooth input output iters" << std::endl;
        return EXIT_FAILURE;
    }
    MeshFileReader reader(argv[1]);
    Mesh& mesh = (Mesh&)reader.GetMesh();
    mesh.RemoveUselessVertices();
    mesh.BuildAllConnectivities();
    mesh.ExtractBoundary();
    int iters = argc >= 4 ? std::stoi(argv[3]) : 1;
    while (iters--)
        Smooth(mesh);
    MeshFileWriter writer(mesh, argv[2]);
    writer.WriteFile();
    return 0;
}
