/*
 * MeshScale.cpp
 *
 *  Created on: Jan 7, 2017
 *      Author: cotrik
 */

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "ArgumentManager.h"
#include <iostream>

void GetArguments(ArgumentManager& argumentManager,
    std::string& orig,// = "tri.vtk";
    std::string& input,// = "quad.vtk";
    std::string& output,// = "out.vtk";
    size_t& iters, // = 1;
    bool& preserveSharpFeature, // = false
    bool& treatSharpFeatureAsRegular,// = true
    bool& treatCornerAsRegular, // = false
    double& cosangle// = 0.939692621
)
{   const std::string strOrigTriFilename = argumentManager.get("orig");
    if (!strOrigTriFilename.empty()) orig = strOrigTriFilename;

    const std::string strInputFilename = argumentManager.get("input");
    if (!strInputFilename.empty()) input = strInputFilename;

    const std::string strOutputFilename = argumentManager.get("output");
    if (!strOutputFilename.empty()) output = strOutputFilename;

    const std::string strIters = argumentManager.get("iters");
    if (!strIters.empty()) iters = std::stoi(strIters);

    const std::string strCosangle = argumentManager.get("cosangle");
    if (!strCosangle.empty()) cosangle = std::stod(strCosangle);

    const std::string strPreserveSharpFeature = argumentManager.get("preserveSharpFeature");
    if (!strPreserveSharpFeature.empty()) preserveSharpFeature = strPreserveSharpFeature == "true" ? true : false;

    const std::string strTreatSharpFeatureAsRegular = argumentManager.get("treatSharpFeatureAsRegular");
    if (!strTreatSharpFeatureAsRegular.empty()) treatSharpFeatureAsRegular = strTreatSharpFeatureAsRegular == "false" ? false : true;

    const std::string strTreatCornerAsRegular = argumentManager.get("treatCornerAsRegular");
    if (!strTreatCornerAsRegular.empty()) treatCornerAsRegular = strTreatCornerAsRegular == "true" ? true : false;

    std::cout << "-----------------------------------\n";
    std::cout << "orig = " << orig << std::endl;
    std::cout << "input = " << input << std::endl;
    std::cout << "output = " << output << std::endl;
    std::cout << "iters = " << iters << std::endl;
    std::cout << "cosangle = " << cosangle << std::endl;
    std::cout << "preserveSharpFeature = " << preserveSharpFeature << std::endl;
    std::cout << "treatSharpFeatureAsRegular = " << treatSharpFeatureAsRegular << std::endl;
    std::cout << "treatCornerAsRegular = " << treatCornerAsRegular << std::endl;
    std::cout << "-----------------------------------\n";
}

int main(int argc, char* argv[])
{
    if (argc < 3)
    {
        std::cout << "Usage: MeshScale input output scale=<0.025> auto=<true> orig=<input.vtk>\n\n";
        return -1;
    }
    ArgumentManager argumentManager(argc, argv);
    double scale = 0.025;
    std::string orig;
    bool Auto = true;
    {
        const std::string strScale = argumentManager.get("scale");
        if (!strScale.empty()) scale = std::stod(strScale);

        const std::string strAuto = argumentManager.get("auto");
        if (!strAuto.empty()) Auto = strAuto == "false" ? false : true;

        const std::string strScaleBackRef = argumentManager.get("orig");
        if (!strScaleBackRef.empty()) orig = strScaleBackRef;

        std::cout << "scale = " << scale << std::endl;
        std::cout << "auto = " << Auto << std::endl;
        std::cout << "orig = " << orig << std::endl;
    }

    if (orig.empty()) {
        MeshFileReader reader(argv[1]);
        Mesh& mesh = (Mesh&)reader.GetMesh();

        if (Auto) {
            mesh.BuildAllConnectivities();
            mesh.ExtractBoundary();
            mesh.GetAvgEdgeLength();
            scale = mesh.avgEdgeLength / scale;
        }

        mesh.m_center = Util::GetCenter(mesh.V);
        mesh.Zoom(mesh.m_center, 1.0/scale);
        std::cout << "center = (" << mesh.m_center.x << ", " << mesh.m_center.y << ", " << mesh.m_center.z << ")\n";
        std::cout << "scale = " << scale << " 1.0/scale = " << 1.0/scale << "\n";

        MeshFileWriter writer(mesh, argv[2]);
        writer.WriteFile();
    }
    else {
        MeshFileReader origReader(orig.c_str());
        Mesh& origMesh = (Mesh&)origReader.GetMesh();
        MeshFileReader reader(argv[1]);
        Mesh& mesh = (Mesh&)reader.GetMesh();

        if (Auto) {
            origMesh.BuildAllConnectivities();
            origMesh.ExtractBoundary();
            origMesh.GetAvgEdgeLength();
            scale = origMesh.avgEdgeLength / scale;
            //cout << "scale = " << scale << "\n";
        }
        mesh.m_center = Util::GetCenter(origMesh.V);
        if (Auto) mesh.Zoom(mesh.m_center, scale);
        else
            mesh.Zoom(mesh.m_center, 1.0 / scale);

        std::cout << "center = (" << mesh.m_center.x << ", " << mesh.m_center.y << ", " << mesh.m_center.z << ")\n";
        std::cout << "scale = " << scale << "\n";

        MeshFileWriter writer(mesh, argv[2]);
        writer.WriteFile();
    }


    return 0;
}



