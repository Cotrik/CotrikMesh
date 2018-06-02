/*
 * SliceMesh.cpp
 *
 *  Created on: Dec 3, 2017
 *      Author: cotrik
 */

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "ArgumentManager.h"
#include <iostream>
#include <sstream>
#include <iomanip>


void printMeshSizeInfo(const Mesh& mesh) {
    glm::vec3 mn(INT_MAX, INT_MAX, INT_MAX);
    glm::vec3 mx(INT_MIN, INT_MIN, INT_MIN);
    for (const auto& v : mesh.V) {
        if (v.x < mn.x) mn.x = v.x;
        if (v.y < mn.y) mn.y = v.y;
        if (v.z < mn.z) mn.z = v.z;

        if (v.x > mx.x) mx.x = v.x;
        if (v.y > mx.y) mx.y = v.y;
        if (v.z > mx.z) mx.z = v.z;
    }
    glm::vec3 s = mx - mn;
    std::cout << "######## mesh size info ########\n";
    std::cout << "x : " << mn.x << " ~ " << mx.x << " size : " << s.x << "\n";
    std::cout << "y : " << mn.y << " ~ " << mx.y << " size : " << s.y << "\n";
    std::cout << "z : " << mn.z << " ~ " << mx.z << " size : " << s.z << "\n";
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: SliceMesh <file>\n";
        return -1;
    }
    ArgumentManager am(argc, argv);
    MeshFileReader reader(argv[1]);
    Mesh& mesh = (Mesh&)reader.GetMesh();
    printMeshSizeInfo(mesh);

    return 0;
}

