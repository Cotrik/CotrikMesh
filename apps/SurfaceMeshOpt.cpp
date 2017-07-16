/*
 * surfacemeshopt.cpp
 *
 *  Created on: April 8, 2017
 *      Author: cotrik
 */

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "EdgeLines.h"
#include "SurfaceMeshOpt.h"
#include <iostream>
#include "ArgumentManager.h"
int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cout << "Usage: MeshOpt input=input.vtk iters=50 alpha=1000000 beta=1000000 gamma=1.0 cosangle=0.984807753 "
                  << "stepSize=0.9 anisotropy=0.5 useAvgLength=false converge=true allowBigStep=false\n";
        return -1;
    }
    ArgumentManager argumentManager(argc, argv);
    std::string filename = "input.vtk";
    size_t iters = 50;
    double alpha = 1000000.0;
    double beta = 1000000.0;
    double gamma = 1.0;
    double stepSize = 0.9;
    double anisotropy = 0.05;
    double cosangle = 0.984807753;
    double minScaledJacobian = 0.0;
    bool useAvgLength = false;
    bool allowBigStep = false;
    bool changeBoundary = false;
    bool converge = true;
    {
        const std::string strFilename = argumentManager.get("input");
        if (!strFilename.empty()) filename = strFilename;

        const std::string strIters = argumentManager.get("iters");
        if (!strIters.empty()) iters = std::stoi(strIters);

        const std::string strAlpha = argumentManager.get("alpha");
        if (!strAlpha.empty()) alpha = std::stod(strAlpha);

        const std::string strBeta = argumentManager.get("beta");
        if (!strBeta.empty()) beta = std::stod(strBeta);

        const std::string strGamma = argumentManager.get("gamma");
        if (!strGamma.empty()) gamma = std::stod(strGamma);

        const std::string strStepSize = argumentManager.get("stepSize");
        if (!strStepSize.empty()) stepSize = std::stod(strStepSize);

        const std::string strAnisotropy = argumentManager.get("anisotropy");
        if (!strAnisotropy.empty()) anisotropy = std::stod(strAnisotropy);

        const std::string strMinScaledJacobian = argumentManager.get("minScaledJacobian");
        if (!strMinScaledJacobian.empty()) minScaledJacobian = std::stod(strMinScaledJacobian);

        const std::string strCosangle = argumentManager.get("cosangle");
        if (!strCosangle.empty()) cosangle = std::stod(strCosangle);

        const std::string strUseAvgLength = argumentManager.get("useAvgLength");
        if (!strUseAvgLength.empty())
            useAvgLength = strUseAvgLength == "true" ? true : false;

        const std::string strConverge = argumentManager.get("converge");
        if (!strConverge.empty())
            converge = strConverge == "false" ? false : true;

        const std::string strAllowBigStep = argumentManager.get("allowBigStep");
        if (!strAllowBigStep.empty())
            allowBigStep = strAllowBigStep == "true" ? true : false;

        const std::string strChangeBoundary = argumentManager.get("changeBoundary");
        if (!strChangeBoundary.empty())
            changeBoundary = strChangeBoundary == "true" ? true : false;

        std::cout << "-----------------------------------\n";
        std::cout << "input = " << filename << std::endl;
        std::cout << "iters = " << iters << std::endl;
        std::cout << "alpha = " << alpha << std::endl;
        std::cout << "beta = " << beta << std::endl;
        std::cout << "gamma = " << gamma << std::endl;
        std::cout << "stepSize = " << stepSize << std::endl;
        std::cout << "anisotropy = " << anisotropy << std::endl;
        std::cout << "minScaledJacobian = " << minScaledJacobian << std::endl;
        std::cout << "cosangle = " << cosangle << std::endl;
        std::cout << "useAvgLength = " << useAvgLength << std::endl;
        std::cout << "converge = " << converge << std::endl;
        std::cout << "allowBigStep = " << allowBigStep << std::endl;
        std::cout << "changeBoundary = " << changeBoundary << std::endl;
        std::cout << "-----------------------------------\n";
    }
    MeshFileReader reader(filename.c_str());
    Mesh& mesh = (Mesh&)reader.GetMesh();
    mesh.RemoveUselessVertices();
    mesh.BuildAllConnectivities();
    mesh.ExtractBoundary();
    mesh.ExtractSingularities();
    mesh.SetCosAngleThreshold(cosangle);
    mesh.LabelSurface();
    mesh.BuildParallelE();
    mesh.BuildConsecutiveE();
    mesh.BuildOrthogonalE();
    mesh.GetNormalOfSurfaceFaces();
    mesh.GetNormalOfSurfaceVertices();
    std::cout << "genus = " <<  1 - (mesh.V.size() - mesh.E.size() + mesh.F.size() - mesh.C.size()) << std::endl;
    std::cout << "#V:" << mesh.V.size() << " - #E:" << mesh.E.size() << " + #F:" << mesh.F.size() << " - #C:" << mesh.C.size() << " = " << "1 - genus" << std::endl;


    size_t nF = 0;
    for (size_t i = 0; i < mesh.F.size(); i++)
        if (mesh.F.at(i).isBoundary) ++nF;
    std::vector<Face> F(nF);
    nF = 0;
    for (size_t i = 0; i < mesh.F.size(); i++)
        if (mesh.F.at(i).isBoundary) {
            F[nF].id = nF;
            F[nF++].Vids = mesh.F.at(i).Vids;
        }
    std::vector<Cell> C(nF);
    nF = 0;
    for (size_t i = 0; i < mesh.F.size(); i++)
        if (mesh.F.at(i).isBoundary) {
            C[nF].id = nF;
            C[nF++].Vids = mesh.F.at(i).Vids;
        }
    //Mesh surfaceMesh(mesh.V, C, QUAD);

    MeshFileWriter surfaceWriter(mesh.V, C, "surface.quad.off", QUAD);
    surfaceWriter.WriteFile();
    MeshFileReader surfaceReader("surface.quad.off");
    Mesh& surfaceMesh = (Mesh&)surfaceReader.GetMesh();

    surfaceMesh.RemoveUselessVertices();
    surfaceMesh.BuildAllConnectivities();
    surfaceMesh.ExtractBoundary();
    surfaceMesh.ExtractSingularities();
    surfaceMesh.SetCosAngleThreshold(cosangle);
    surfaceMesh.LabelSurface();
    surfaceMesh.BuildParallelE();
    surfaceMesh.BuildConsecutiveE();
    surfaceMesh.BuildOrthogonalE();
    surfaceMesh.GetNormalOfSurfaceFaces();
    surfaceMesh.GetNormalOfSurfaceVertices();

    SurfaceMeshOpt meshOpt(surfaceMesh);
    meshOpt.SetAlpha(alpha);
    meshOpt.SetBeta(beta);
    meshOpt.SetGamma(gamma);
    meshOpt.SetStepSize(stepSize);
    meshOpt.SetUseAverageTargetLength(useAvgLength);
    meshOpt.SetAnisotropy(anisotropy);
    meshOpt.SetMinScaledJacobian(minScaledJacobian);
    meshOpt.SetRecoverable(converge);
    meshOpt.SetAllowBigStep(allowBigStep);
//    meshOpt.SetChangeBoundary(changeBoundary);
    meshOpt.SetRefMesh(mesh);
    size_t iter = meshOpt.Run(iters);

    return 0;
}
