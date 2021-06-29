/*
 * ConvertFileFormat.cpp
 *
 *  Created on: Dec 19, 2019
 *      Author: Naeem
 */

#include "MeshFileReader.h"
#include <iostream>
#include "ArgumentManager.h"
int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cout << "Usage: Fix2DMeshVertices input_file[.vtk] output_file.txt\n";
        return -1;
    }
    ArgumentManager argumentManager(argc, argv);
    MeshFileReader reader(argv[1]);
    Mesh& mesh = (Mesh&)reader.GetMesh();
    mesh.BuildAllConnectivities();
    mesh.unifyOrientation();

    std::ofstream ofs(argv[2]);
    std::vector<int> tags(mesh.V.size(), 0);
    const Vertex& v0_r = mesh.V[mesh.C.at(0).Vids[0]];
    const Vertex& v1_r = mesh.V[mesh.C.at(0).Vids[1]];
    const Vertex& v2_r = mesh.V[mesh.C.at(0).Vids[2]];
    const glm::dvec3 v10_r = v0_r.xyz() - v1_r.xyz();
    const glm::dvec3 v12_r = v2_r.xyz() - v1_r.xyz();
    const glm::dvec3 n_r = glm::normalize(glm::cross(v12_r, v10_r));
    ofs << n_r.x << " " << n_r.y << " " << n_r.z << std::endl;
    for (auto &v : mesh.V) {
        if (v.N_Vids.size() < 4) {
            tags[v.id] = 1;
        }
        if (v.N_Vids.size() == 2) {
            tags[mesh.V.at(v.N_Vids[0]).id] = 1;
            tags[mesh.V.at(v.N_Vids[1]).id] = 1;
        }
    }
    for (int i = 0; i < tags.size(); i++) {
        ofs << tags[i];
        if (i != tags.size() - 1) {
            ofs << std::endl;
        }
    }
    return 0;
}
