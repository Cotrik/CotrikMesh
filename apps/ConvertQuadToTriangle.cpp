/*
 * convert_quad_to_triangle.cpp
 *
 *  Created on: Apr 8, 2017
 *      Author: cotrik
 */

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "SmoothAlgorithm.h"
#include <iostream>
#include "ArgumentManager.h"
int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cout << "Usage: convert_quad_to_triangle input.vtk output.vtk\n";
        return -1;
    }
    ArgumentManager argumentManager(argc, argv);
    std::string input_file_name(argv[1]);

    MeshFileReader reader(input_file_name.c_str());
    Mesh& mesh = (Mesh&)reader.GetMesh();
    mesh.RemoveUselessVertices();
    mesh.BuildAllConnectivities();
    mesh.ExtractBoundary();
    //mesh.ExtractSingularities();
    //mesh.SetCosAngleThreshold(cosangle);
    //mesh.LabelSurface();
    //mesh.BuildParallelE();
    //mesh.BuildConsecutiveE();
    //mesh.BuildOrthogonalE();
    //mesh.GetNormalOfSurfaceFaces();
    //mesh.GetNormalOfSurfaceVertices();

//    SmoothAlgorithm smoother(mesh, LAPLACIAN);
//    const int iters = 1e2;
//    const double eps = 1e-4;
//    smoother.Run(iters, eps);

    std::vector<Vertex> V(mesh.V.size() + mesh.F.size());
    std::vector<Cell> C(4 * mesh.F.size());

    for (size_t i = 0; i < mesh.V.size(); i++) {
        V[i] = mesh.V[i].xyz();
        V[i].id = i;
    }
    const size_t off = mesh.V.size();
    for (size_t i = 0; i < mesh.F.size(); i++) {
        const Face& f = mesh.F.at(i);
        const glm::dvec3 v = mesh.V[f.Vids[0]].xyz() + mesh.V[f.Vids[2]].xyz();
        const glm::dvec3 cv(0.5*v.x, 0.5*v.y, 0.5*v.z);
        V[off + i] = cv;
        V[off + i].id = off + i;
    }

    for (size_t i = 0; i < mesh.F.size(); i++) {
        const Face& f = mesh.F.at(i);
        const size_t cvid = off + i;
        Cell tri1(3);
        tri1.Vids[0] = cvid;
        tri1.Vids[1] = f.Vids[0];
        tri1.Vids[2] = f.Vids[1];
        tri1.id = 4 * i + 0;

        Cell tri2(3);
        tri2.Vids[0] = cvid;
        tri2.Vids[1] = f.Vids[1];
        tri2.Vids[2] = f.Vids[2];
        tri2.id = 4 * i + 1;

        Cell tri3(3);
        tri3.Vids[0] = cvid;
        tri3.Vids[1] = f.Vids[2];
        tri3.Vids[2] = f.Vids[3];
        tri3.id = 4 * i + 2;

        Cell tri4(3);
        tri4.Vids[0] = cvid;
        tri4.Vids[1] = f.Vids[3];
        tri4.Vids[2] = f.Vids[0];
        tri4.id = 4 * i + 3;

        C[4 * i + 0] = tri1;
        C[4 * i + 1] = tri2;
        C[4 * i + 2] = tri3;
        C[4 * i + 3] = tri4;
    }

    MeshFileWriter writer(V, C, argv[2], TRIANGLE);
    writer.WriteFile();
    return 0;
}


