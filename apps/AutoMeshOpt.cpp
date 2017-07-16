/*
 * automeshopt.cpp
 *
 *  Created on: Nov 6, 2016
 *      Author: cotrik
 */

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "AutoMeshOpt.h"
#include <iostream>
#include "ArgumentManager.h"
int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cout << "Usage: AutoMeshOpt input=<input.vtk> orig=<orig.vtk> hausdorff=<0.01>\n";
        return -1;
    }
    ArgumentManager argumentManager(argc, argv);
    std::string input = "input.vtk";
    std::string orig = "orig.vtk";
    double hausdorff = 0.01;
    {
        const std::string strFilename = argumentManager.get("input");
        if (!strFilename.empty()) input = strFilename;
        const std::string strOrig = argumentManager.get("orig");
        if (!strOrig.empty()) orig = strOrig;
        const std::string strHausdorff = argumentManager.get("hausdorff");
        if (!strHausdorff.empty()) hausdorff = std::stod(strHausdorff);
        std::cout << "-----------------------------------\n";
        std::cout << "input = " << input << std::endl;
        std::cout << "orig = " << orig << std::endl;
        std::cout << "hausdorff = " << hausdorff << std::endl;
        std::cout << "-----------------------------------\n";
    }
    MeshFileReader inputReader(input.c_str());
    Mesh& mesh = (Mesh&)inputReader.GetMesh();
    mesh.RemoveUselessVertices();
    mesh.BuildAllConnectivities();
    mesh.ExtractBoundary();
    mesh.ExtractSingularities();
    mesh.SetCosAngleThreshold(0.939692621);
    mesh.LabelSurface();
    mesh.BuildParallelE();
    mesh.BuildConsecutiveE();
    mesh.BuildOrthogonalE();
    mesh.GetNormalOfSurfaceFaces();
    mesh.GetNormalOfSurfaceVertices();
    std::cout << "genus = " <<  1 - (mesh.V.size() - mesh.E.size() + mesh.F.size() - mesh.C.size()) << std::endl;
    std::cout << "#V:" << mesh.V.size() << " - #E:" << mesh.E.size() << " + #F:" << mesh.F.size() << " - #C:" << mesh.C.size() << " = " << "1 - genus" << std::endl;

    MeshFileReader origReader(orig.c_str());
    Mesh& origMesh = (Mesh&)origReader.GetMesh();
    origMesh.BuildAllConnectivities();
    origMesh.ExtractBoundary();
    origMesh.GetNormalOfSurfaceFaces();
    origMesh.GetNormalOfSurfaceVertices();

    AutoMeshOpt meshOpt(mesh, origMesh);
    meshOpt.SetHausdorffError(hausdorff);
    meshOpt.Run();

    return 0;
}
