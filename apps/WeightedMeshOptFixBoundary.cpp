/*
 * main.cpp
 *
 *  Created on: Nov 6, 2016
 *      Author: cotrik
 */

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "EdgeLines.h"
#include "WeightedMeshOptFixBoundary.h"
#include <iostream>
#include "ArgumentManager.h"
#include "MeshQuality.h"
int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cout << "Usage: WeightedMeshOptFixBoundary input=input.vtk iters=50 alpha=1000000 beta=1.0 gamma=1.0 cosangle=0.965925826"
                  << "stepSize=0.9 anisotropy=0.05 useAvgLength=false converge=true allowBigStep=false changeBoundary=false minScaledJacobian=0.0\n";
        return -1;
    }
    ArgumentManager argumentManager(argc, argv);
    std::string filename = "input.vtk";
    size_t iters = 50;
    double alpha = 1000000.0;
    double beta = 1.0;
    double gamma = 1.0;
    double stepSize = 0.9;
    double anisotropy = 0.05;
    double cosangle = 0.965925826;
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
        if (!strUseAvgLength.empty()) useAvgLength = strUseAvgLength == "true" ? true : false;
        const std::string strConverge = argumentManager.get("converge");
        if (!strConverge.empty()) converge = strConverge == "false" ? false : true;
        const std::string strAllowBigStep = argumentManager.get("allowBigStep");
        if (!strAllowBigStep.empty()) allowBigStep = strAllowBigStep == "true" ? true : false;
        const std::string strChangeBoundary = argumentManager.get("changeBoundary");
        if (!strChangeBoundary.empty()) changeBoundary = strChangeBoundary == "true" ? true : false;

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
//    mesh.ExtractLayers();
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

//    double minimumScaledJacobian = 0.0;
//    std::vector<size_t> badCellIds;
//    size_t m_numOfInvertdElements = GetMinScaledJacobianVerdict(mesh, minimumScaledJacobian, badCellIds, minScaledJacobian);
//    std::cout << "iter = " << 0 << " #inverted = " << m_numOfInvertdElements << " MSJ = " << minimumScaledJacobian << "#badCells = " << badCellIds.size() << std::endl;
//    std::vector<size_t> remainIds;
//    for (size_t j = 0; j < mesh.C.size(); j++) {
//        bool found = false;
//        for (size_t i = 0; i < badCellIds.size(); i++)
//            if (badCellIds.at(i) == j) {
//                found = true;
//                break;
//            }
//        if (!found) remainIds.push_back(j);
//    }
//
//    Mesh remainMesh(mesh, remainIds);
//    remainMesh.RemoveUselessVertices();
//    remainMesh.BuildAllConnectivities();
//    remainMesh.ExtractBoundary();
//    remainMesh.ExtractSingularities();
//    remainMesh.SetCosAngleThreshold(cosangle);
//    remainMesh.LabelSurface();
//    remainMesh.BuildParallelE();
//    remainMesh.BuildConsecutiveE();
//    remainMesh.BuildOrthogonalE();
//    remainMesh.GetNormalOfSurfaceFaces();
//    remainMesh.GetNormalOfSurfaceVertices();

    WeightedMeshOptFixBoundary meshOpt(mesh);
    meshOpt.SetAlpha(alpha);
    meshOpt.SetBeta(beta);
    meshOpt.SetGamma(gamma);
    meshOpt.SetStepSize(stepSize);
    meshOpt.SetUseAverageTargetLength(useAvgLength);
    meshOpt.SetAnisotropy(anisotropy);
    meshOpt.SetMinScaledJacobian(minScaledJacobian);
    meshOpt.SetRecoverable(converge);
    meshOpt.SetAllowBigStep(allowBigStep);
    meshOpt.SetChangeBoundary(changeBoundary);
    meshOpt.SetRefMesh(mesh);
    size_t iter = meshOpt.Run(iters);

    return 0;
}

//int main(int argc, char* argv[])
//{
//    if (argc < 2)
//    {
//        std::cout << "Usage: MeshOpt <file> <alpha> <iters> <!useAvgTargetLength>\n";
//        return -1;
//    }
//    MeshFileReader reader(argv[1]);
//    Mesh& mesh = (Mesh&)reader.GetMesh();
//    mesh.BuildAllConnectivities();
//    mesh.ExtractBoundary();
//    mesh.ExtractSingularities();
//    mesh.BuildParallelE();
//    mesh.BuildConsecutiveE();
//    mesh.BuildOrthogonalE();
//    std::cout << "genus = " <<  1 - (mesh.V.size() - mesh.E.size() + mesh.F.size() - mesh.C.size()) << std::endl;
//    std::cout << "#V:" << mesh.V.size() << " - #E:" << mesh.E.size() << " + #F:" << mesh.F.size() << " - #C:" << mesh.C.size() << " = " << "1 - genus" << std::endl;
//
//    MeshOpt meshOpt(mesh, 1000000.0, 1.0, 1.0);
//    if (argc >= 3)
//        meshOpt.SetAlpha(std::stod(argv[2]));
//    size_t iters = 10;
//    if (argc >= 4)
//        iters = std::stoi(argv[3]);
//    if (argc >= 5)
//        meshOpt.SetUseAverageTargetLength(std::stoi(argv[4]) == 0 ? false : true);
//    if (argc >= 6)
//        meshOpt.SetStepSize(std::stod(argv[5]));
//    if (argc >= 7)
//        meshOpt.SetAnisotropy(std::stod(argv[6]));
//    if (argc >= 8)
//        meshOpt.SetRecoverable(std::stoi(argv[7]) == 0 ? false : true);
//    meshOpt.Run(iters);
//
//    return 0;
//}
