/*
 * MeshShrink.cpp
 *
 *  Created on: Feb 27, 2017
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

const float PI2 = 3.1415926f * 2;
float GetAngle(Mesh& mesh, const Vertex& v, const Face& c)
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
    if (argc < 6)
    {
        std::cerr << "Usage: ExtractGaussianCurvature input_tri_file iters base_edge_length_ratio curvature_base smooth_iters\n";
        return -1;
    }

    MeshFileReader reader(argv[1]);
    //reader.GetPointsScalarFields();
    Mesh& mesh = (Mesh&)reader.GetMesh();
    mesh.BuildAllConnectivities();
    mesh.ExtractBoundary();

    std::vector<bool> stop(mesh.V.size(), false);
    std::vector<int> count(mesh.V.size(), 0);
    for (size_t i = 0; i < mesh.V.size(); i++) {
        Vertex& v = mesh.V.at(i);
        v.type = REGULAR;
    }
    int iters = std::stoi(argv[2]);
    double base_edge_length_ratio = std::stof(argv[3]);
    double curvature_base = std::stof(argv[4]);
    int smooth_iters = std::stoi(argv[5]);
    int iter = 0;
    while (iter++ != iters) {
        std::cout << "---------- iter " << iter << " ----------\n";
        //mesh.ClearLabelOfSurface();
        //mesh.LabelSurface();
        //mesh.SmoothSurface(1, LAPLACE_EDGE, true, false, false);
        mesh.SmoothSurface(1, LAPLACE_EDGE, false, true, true);
    std::vector<float> gaussianCurvatureOnVertices(mesh.V.size());
    for (size_t i = 0; i < mesh.V.size(); i++) {
        const Vertex& v = mesh.V.at(i);
        for (size_t j = 0; j < v.N_Fids.size(); j++) {
            const Face& f = mesh.F.at(v.N_Fids.at(j));
            gaussianCurvatureOnVertices.at(i) += GetAngle(mesh, v, f);
        }
        gaussianCurvatureOnVertices.at(i) = PI2 - gaussianCurvatureOnVertices.at(i);
        if (gaussianCurvatureOnVertices.at(i) < -PI2)
            gaussianCurvatureOnVertices.at(i) = -PI2;
    }

//    MeshFileWriter writer(mesh, argv[2]);
//    writer.WriteFile();
//    writer.WritePointData(gaussianCurvatureOnVertices, "Gaussian");

    double sum_edge_length = 0.0;
    for (size_t i = 0; i < mesh.E.size(); i++) {
        const Edge& e = mesh.E.at(i);
        const Vertex& v1 = mesh.V.at(e.Vids.at(0));
        const Vertex& v2 = mesh.V.at(e.Vids.at(1));
        sum_edge_length += glm::length(v1.xyz() - v2.xyz());
    }
    double avg_edge_length = sum_edge_length/mesh.E.size();
    double base_edge_length = base_edge_length_ratio * avg_edge_length;
    std::cout << "base_edge_length = " << base_edge_length << std::endl;

    mesh.GetNormalOfSurfaceFaces();
    mesh.GetNormalOfSurfaceVertices();

    std::vector<Vertex> newV = mesh.V;
#pragma omp parallel for
    for (size_t i = 0; i < mesh.V.size(); i++) {
        Vertex& v = mesh.V.at(i);
        //if (v.type == FEATURE && v.type == CORNER) continue;

        glm::dvec3 new_v(0.0, 0.0, 0.0);
        new_v.x = v.x - v.normal.x * base_edge_length / (curvature_base + gaussianCurvatureOnVertices.at(i));
        new_v.y = v.y - v.normal.y * base_edge_length / (curvature_base + gaussianCurvatureOnVertices.at(i));
        new_v.z = v.z - v.normal.z * base_edge_length / (curvature_base + gaussianCurvatureOnVertices.at(i));

            if (!stop[i] && mesh.IsPointInside(new_v)) {
                //v = new_v;
                newV[i] = new_v;
            }
            else {
                count[i]++;
                if (count[i] > smooth_iters) {
                    v.type = MAXID;
                    stop[i] = true;
                }
            }
    }
    mesh.V = newV;

    std::string filename = std::string("shrink.") + std::to_string(iter) + ".vtk";
    MeshFileWriter writer3(mesh, filename.c_str());
    writer3.WriteFile();
    writer3.WritePointData(gaussianCurvatureOnVertices, "Gaussian");
    //mesh.SmoothSurface(1, LAPLACE_EDGE, true, false, true);
    }
    return 0;
}


