/*
 * FixedSmooth.cpp
 *
 *  Created on: Jan 18, 2018
 *      Author: cotrik
 */

#include "ArgumentManager.h"
#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include <iostream>
#include <fstream>

struct Id_Coordinate {
    size_t id = MAXID;
    glm::vec3 xyz;
    bool operator == (const Id_Coordinate& rhs) const {
        return id == rhs.id;
    }
};

std::vector<Id_Coordinate> GetId_Coordinates(const char* filename) {
    std::vector<Id_Coordinate> res;
    std::ifstream ifs(filename);
    Id_Coordinate idCoordinate;
    while (ifs >> idCoordinate.id >> idCoordinate.xyz.x >> idCoordinate.xyz.y >> idCoordinate.xyz.z)
        res.push_back(idCoordinate);

    return res;
}

void Smooth(Mesh& mesh, const std::vector<Id_Coordinate>& idCoordinates) {
    for (auto& idCoordinate : idCoordinates) {
        auto& v = mesh.V.at(idCoordinate.id);
        v = idCoordinate.xyz;
        v.isActive = false;
    }

    for (auto& v: mesh.V) {
        if (!v.isActive) continue;
        glm::vec3 sum(0.0, 0.0, 0.0);
        for (auto vid : v.N_Vids)
            sum += mesh.V.at(vid).xyz();
        v = sum / float(v.N_Vids.size());
    }
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: TriangleSmooth input output iters=1 fixed_points_file=fixed.txt" << std::endl;
        return EXIT_FAILURE;
    }
    ArgumentManager am(argc, argv);
    std::string str_iters = am.get("iters");
    std::string str_fixed_points_file = am.get("fixed_points_file");
    MeshFileReader reader(argv[1]);
    Mesh& mesh = (Mesh&)reader.GetMesh();
    mesh.RemoveUselessVertices();
    mesh.BuildAllConnectivities();
    mesh.ExtractBoundary();
    int iters = str_iters.empty() ? 1 : std::stoi(str_iters);
    auto idCoordinates = GetId_Coordinates(str_fixed_points_file.c_str());
    while (iters--)
        Smooth(mesh, idCoordinates);
    MeshFileWriter writer(mesh, argv[2]);
    writer.WriteFile();
    return 0;
}



