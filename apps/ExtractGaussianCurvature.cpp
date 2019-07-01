/*
 * ExtractGaussianCurvature.cpp
 *
 *  Created on: Feb 23, 2017
 *      Author: cotrik
 */


#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include <iostream>
#include "ArgumentManager.h"
size_t permutation[3][2] = {
        {0, 1},
        {1, 2},
        {2, 0}
};

const double PI2 = 3.1415926 * 2;
double GetAngle(Mesh& mesh, const Vertex& v, const Face& c)
{
    size_t vid1 = 0;
    size_t vid2 = 0;
    for (int i = 0; i < 3; i++) {
        if (v.id != c.Vids[permutation[i][0]] && v.id != c.Vids[permutation[i][1]]) {
            vid1 = c.Vids[permutation[i][0]];
            vid2 = c.Vids[permutation[i][1]];
            break;
        }
    }
    glm::dvec3 v1 = mesh.V.at(vid1).xyz() - v.xyz();
    glm::dvec3 v2 = mesh.V.at(vid2).xyz() - v.xyz();

    return acos(glm::dot(v1, v2) / (glm::length(v1) * glm::length(v2)));
}

int main(int argc, char* argv[])
{
    if (argc < 3)
    {
        std::cerr << "Usage: ExtractGaussianCurvature input.vtk output.vtk\n";
        return -1;
    }

    MeshFileReader reader(argv[1]);
    //reader.GetPointsScalarFields();
    Mesh& mesh = (Mesh&)reader.GetMesh();
    mesh.BuildAllConnectivities();
    //mesh.ExtractBoundary();

    std::vector<double> gaussianCurvatureOnVertices(mesh.V.size());
    for (size_t i = 0; i < mesh.V.size(); i++) {
        const Vertex& v = mesh.V.at(i);
        for (size_t j = 0; j < v.N_Fids.size(); j++) {
            const Face& f = mesh.F.at(v.N_Fids.at(j));
            gaussianCurvatureOnVertices.at(i) += GetAngle(mesh, v, f);
        }
        gaussianCurvatureOnVertices.at(i) = PI2 - gaussianCurvatureOnVertices.at(i);
    }

    MeshFileWriter writer(mesh, argv[2]);
    writer.WriteFile();
    writer.WritePointData(gaussianCurvatureOnVertices, "Gaussian");
    return 0;
}
