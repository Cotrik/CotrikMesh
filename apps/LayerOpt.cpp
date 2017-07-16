/*
 * layeropt.cpp
 *
 *  Created on: Nov 6, 2016
 *      Author: cotrik
 */

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "EdgeLines.h"
#include "MeshOpt.h"
#include <iostream>
#include "ArgumentManager.h"
int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cout << "Usage: LayerOpt input=<input.vtk> iters=<10> alpha=<1000000> gamma=<1.0> stepSize=<0.9> anisotropy=<0.05> useAvgLength=<false> converge=<true>\n";
        return -1;
    }
    ArgumentManager argumentManager(argc, argv);
    std::string filename = "input.vtk";
    size_t iters = 50;
    double alpha = 1000000.0;
    double gamma = 1.0;
    double stepSize = 0.9;
    double anisotropy = 0.05;
    bool useAvgLength = false;
    bool converge = true;
    {
        const std::string strFilename = argumentManager.get("input");
        if (!strFilename.empty())
            filename = strFilename;

        const std::string strIters = argumentManager.get("iters");
        if (!strIters.empty())
            iters = std::stoi(strIters);

        const std::string strAlpha = argumentManager.get("alpha");
        if (!strAlpha.empty())
            alpha = std::stod(strAlpha);

        const std::string strGamma = argumentManager.get("gamma");
        if (!strGamma.empty())
            gamma = std::stod(strGamma);

        const std::string strStepSize = argumentManager.get("stepSize");
        if (!strStepSize.empty())
            stepSize = std::stod(strStepSize);

        const std::string strAnisotropy = argumentManager.get("anisotropy");
        if (!strAnisotropy.empty())
            anisotropy = std::stod(strAnisotropy);

        const std::string strUseAvgLength = argumentManager.get("useAvgLength");
        if (!strUseAvgLength.empty())
            useAvgLength = strUseAvgLength == "true" ? true : false;

        const std::string strConverge = argumentManager.get("converge");
        if (!strConverge.empty())
            converge = strConverge == "false" ? false : true;

        std::cout << "-----------------------------------\n";
        std::cout << "input = " << filename << std::endl;
        std::cout << "iters = " << iters << std::endl;
        std::cout << "alpha = " << alpha << std::endl;
        std::cout << "gamma = " << gamma << std::endl;
        std::cout << "stepSize = " << stepSize << std::endl;
        std::cout << "anisotropy = " << anisotropy << std::endl;
        std::cout << "useAvgLength = " << useAvgLength << std::endl;
        std::cout << "converge = " << converge << std::endl;
        std::cout << "-----------------------------------\n";
    }
    MeshFileReader reader(filename.c_str());
    Mesh& mesh = (Mesh&)reader.GetMesh();
    mesh.BuildAllConnectivities();
    mesh.ExtractBoundary();
    const size_t numOfLayers = mesh.ExtractLayers();
    mesh.ExtractSingularities();
    mesh.BuildParallelE();
    mesh.BuildConsecutiveE();
    mesh.BuildOrthogonalE();
    mesh.GetNormalOfSurfaceFaces();
    mesh.GetNormalOfSurfaceVertices();
    std::cout << "--------------------------------\n";
    std::cout << "Mesh : " << filename << "\n";
    std::cout << "genus = " <<  1 - (mesh.V.size() - mesh.E.size() + mesh.F.size() - mesh.C.size()) << std::endl;
    std::cout << "#V:" << mesh.V.size() << " - #E:" << mesh.E.size() << " + #F:" << mesh.F.size() << " - #C:" << mesh.C.size() << " = " << "1 - genus" << std::endl;

    size_t count = 0;
    while (count++ < 1){
    //while (count++ < numOfLayers){
        std::string layerFilename = std::string("OutLayer") + std::to_string(count) +  ".vtk";
        MeshFileReader layerReader(layerFilename.c_str());
        Mesh& outLayer = (Mesh&)layerReader.GetMesh();
        //outLayer.RemoveUselessVertices();
        outLayer.BuildAllConnectivities();
        outLayer.ExtractBoundary();
        //mesh1.ExtractLayers();
        outLayer.ExtractSingularities();
        outLayer.BuildParallelE();
        outLayer.BuildConsecutiveE();
        outLayer.BuildOrthogonalE();
        outLayer.GetNormalOfSurfaceFaces();
        outLayer.GetNormalOfSurfaceVertices();
        std::cout << "--------------------------------\n";
        std::cout << "Mesh : " << layerFilename << "\n";
        std::cout << "genus = " <<  1 - (outLayer.V.size() - outLayer.E.size() + outLayer.F.size() - outLayer.C.size()) << std::endl;
        std::cout << "#V:" << outLayer.V.size() << " - #E:" << outLayer.E.size() << " + #F:" << outLayer.F.size() << " - #C:" << outLayer.C.size() << " = " << "1 - genus" << std::endl;

        for (size_t i = 0; i < outLayer.V.size(); i++)
            outLayer.V[i].isBoundary = mesh.V[mesh.layers[0].Vids[i]].isBoundary;

        MeshOpt layerOpt(outLayer);
        layerOpt.AdjustOutlayer();
        layerOpt.SetAlpha(alpha);
        layerOpt.SetGamma(gamma);
        layerOpt.SetStepSize(stepSize);
        layerOpt.SetUseAverageTargetLength(useAvgLength);
        layerOpt.SetAnisotropy(anisotropy);
        layerOpt.SetRecoverable(converge);
        layerOpt.Run(iters);

        for (size_t i = 0; i < outLayer.V.size(); i++) {
            mesh.V[mesh.layers[0].Vids[i]].x = outLayer.V[i].x;
            mesh.V[mesh.layers[0].Vids[i]].y = outLayer.V[i].y;
            mesh.V[mesh.layers[0].Vids[i]].z = outLayer.V[i].z;
        }

        std::string outLayerFilename = std::string("OutLayerOpt.") + std::to_string(count) +  ".vtk";
        MeshFileWriter outLayerWriter(outLayer, outLayerFilename.c_str());
        outLayerWriter.WriteFile();


        std::string innerLayerFilename = std::string("InnerLayer") + std::to_string(count) +  ".vtk";
        MeshFileReader innerLayerReader(innerLayerFilename.c_str());
        Mesh& innerLayer = (Mesh&)innerLayerReader.GetMesh();
        //innerLayer.RemoveUselessVertices();
        innerLayer.BuildAllConnectivities();
        innerLayer.ExtractBoundary();
        innerLayer.ExtractSingularities();
        innerLayer.BuildParallelE();
        innerLayer.BuildConsecutiveE();
        innerLayer.BuildOrthogonalE();
        std::cout << "--------------------------------\n";
        std::cout << "Mesh : " << innerLayerFilename << "\n";
        std::cout << "genus = " <<  1 - (innerLayer.V.size() - innerLayer.E.size() + innerLayer.F.size() - innerLayer.C.size()) << std::endl;
        std::cout << "#V:" << innerLayer.V.size() << " - #E:" << innerLayer.E.size() << " + #F:" << innerLayer.F.size() << " - #C:" << innerLayer.C.size() << " = " << "1 - genus" << std::endl;

        for (size_t i = 0; i < innerLayer.V.size(); i++) {
            innerLayer.V[i].x = mesh.V[mesh.innerLayers[0].Vids[i]].x;
            innerLayer.V[i].y = mesh.V[mesh.innerLayers[0].Vids[i]].y;
            innerLayer.V[i].z = mesh.V[mesh.innerLayers[0].Vids[i]].z;
        }

        std::string innerLayerBeforeFilename = std::string("InnerLayer_before.") + std::to_string(count) +  ".vtk";
        MeshFileWriter innerLayerWriter1(innerLayer, innerLayerBeforeFilename.c_str());
        innerLayerWriter1.WriteFile();

        MeshOpt innerLayerOpt(innerLayer);
        innerLayerOpt.AdjustOutlayer();
        innerLayerOpt.SetAlpha(alpha);
        innerLayerOpt.SetGamma(gamma);
        innerLayerOpt.SetStepSize(stepSize);
        innerLayerOpt.SetUseAverageTargetLength(useAvgLength);
        innerLayerOpt.SetAnisotropy(anisotropy);
        innerLayerOpt.SetRecoverable(converge);
        innerLayerOpt.Run(iters);

        for (size_t i = 0; i < innerLayer.V.size(); i++) {
            mesh.V[mesh.innerLayers[0].Vids[i]].x = innerLayer.V[i].x;
            mesh.V[mesh.innerLayers[0].Vids[i]].y = innerLayer.V[i].y;
            mesh.V[mesh.innerLayers[0].Vids[i]].z = innerLayer.V[i].z;
        }

        std::string innerLayerOptFilename = std::string("InnerLayerOpt.") + std::to_string(count) +  ".vtk";
        MeshFileWriter innerLayerWriter(innerLayer, innerLayerOptFilename.c_str());
        innerLayerWriter.WriteFile();
    }

    MeshFileWriter optwriter(mesh, "layeropt.vtk");
    optwriter.WriteFile();
    return 0;

//    for (size_t i = 0; i < mesh.layers[0].Vids.size(); i++){
//        mesh.V[mesh.layers[0].Vids[i]].isBoundary = true;
//        mesh.V[mesh.layers[0].Vids[i]].x = mesh1.V[i].x;
//        mesh.V[mesh.layers[0].Vids[i]].y = mesh1.V[i].y;
//        mesh.V[mesh.layers[0].Vids[i]].z = mesh1.V[i].z;
//    }

    MeshOpt meshOpt(mesh);
    meshOpt.SetAlpha(alpha);
    meshOpt.SetGamma(gamma);
    meshOpt.SetStepSize(stepSize);
    meshOpt.SetUseAverageTargetLength(useAvgLength);
    meshOpt.SetAnisotropy(anisotropy);
    meshOpt.SetRecoverable(converge);
    meshOpt.Run(iters);

    return 0;
}
