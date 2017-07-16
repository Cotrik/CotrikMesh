/*
 * PerturbInnerVertices.cpp
 *
 *  Created on: Nov 20, 2016
 *      Author: cotrik
 */

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "FrameField.h"
#include "PolyLine.h"
#include "FrameOpt.h"
#include <iostream>
#include <random>

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cout << "Usage: perturb <file> [magnitude=10]\n";
        return -1;
    }

    MeshFileReader reader(argv[1]);
    Mesh& mesh = (Mesh&)reader.GetMesh();
    mesh.BuildAllConnectivities();
    mesh.ExtractBoundary();
    std::cout << "=============================\n";
    double sumEdgeLength = 0.0;
    size_t numOfBoundaryEdges = 0;
    for (size_t i = 0; i < mesh.E.size(); i++)
    {
        //if (mesh.E.at(i).isBoundary)
        {
            const Vertex& v1 = mesh.V.at(mesh.E.at(i).Vids.at(0));
            const Vertex& v2 = mesh.V.at(mesh.E.at(i).Vids.at(1));
            mesh.E[i].length = glm::length(v1 - v2);
            sumEdgeLength += mesh.E[i].length;
            //sumEdgeLength += glm::length(v1 - v2);
            numOfBoundaryEdges++;
        }
    }
    double avgMeshEdgeLength = sumEdgeLength / mesh.E.size();
    std::cout << "Average Edge Length = " << avgMeshEdgeLength << std::endl;

    int magnitude = 10;
    if (argc == 3)
        magnitude = std::stoi(argv[2]);

    srandom(time(NULL));
    for (size_t i = 0; i < mesh.V.size(); i++){
        Vertex& v = mesh.V.at(i);
        if (v.isBoundary)
            continue;
        const double radius = ((double)rand()) / (double)RAND_MAX;  // radius in (0, 1)
        v.x += avgMeshEdgeLength * magnitude * radius;
        v.y += avgMeshEdgeLength * magnitude * radius;
        v.z += avgMeshEdgeLength * magnitude * radius;
    }

    MeshFileWriter fileWriter(mesh, "input.vtk");
    fileWriter.WriteFile();

    return 0;
}


