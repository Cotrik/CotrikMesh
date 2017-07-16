/*
 * HexSmooth.cpp
 *
 *  Created on: Mar 29, 2017
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
        std::cout << "Usage: Smooth input.vtk\n";
        return -1;
    }
    ArgumentManager argumentManager(argc, argv);
    std::string input_file_name(argv[1]);

    MeshFileReader reader(input_file_name.c_str());
    Mesh& mesh = (Mesh&)reader.GetMesh();
    mesh.RemoveUselessVertices();
    mesh.BuildAllConnectivities();
    mesh.ExtractBoundary();
    mesh.ExtractSingularities();
    //mesh.SetCosAngleThreshold(cosangle);
    mesh.LabelSurface();
    //mesh.BuildParallelE();
    //mesh.BuildConsecutiveE();
    //mesh.BuildOrthogonalE();
    mesh.GetNormalOfSurfaceFaces();
    mesh.GetNormalOfSurfaceVertices();

    //SmoothAlgorithm smoother(mesh, LAPLACIAN);
    SmoothAlgorithm smoother(mesh, SCALED_JACOBIAN);
    const int iters = 1e2;
    const double eps = 1e-4;
    smoother.Run(iters, eps);

    std::string filename = "smoothed.hex.vtk";
    MeshFileWriter writer(mesh, filename.c_str());
    writer.WriteFile();
    return 0;
}

