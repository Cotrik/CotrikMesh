/*
 * ExtractFrameField.cpp
 *
 *  Created on: Nov 20, 2016
 *      Author: cotrik
 */

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "FrameField.h"
#include "PolyLine.h"
#include "FrameOpt.h"
#include "ArgumentManager.h"
#include <iostream>

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cout << "Usage: FrameOpt <file>\n";
        return -1;
    }
    MeshFileReader reader(argv[1]);
    Mesh& mesh = (Mesh&)reader.GetMesh();
    mesh.BuildAllConnectivities();
    mesh.ExtractBoundary();
    std::cout << "genus = " <<  1 - (mesh.V.size() - mesh.E.size() + mesh.F.size() - mesh.C.size()) << std::endl;
    std::cout << "#V:" << mesh.V.size() << " - #E:" << mesh.E.size() << " + #F:" << mesh.F.size() << " - #C:" << mesh.C.size() << " = " << "1 - genus" << std::endl;

    FrameField framefield;
    framefield.BuildFrom(mesh);
    std::string filename = std::string(argv[1]).substr(0, std::string(argv[1]).size() - 4) + ".framefield.vtk";
    framefield.WriteFile(/*"FrameField.vtk"*/filename.c_str());
    std::cout << "Frame field extracted" << std::endl;

    PolyLines polylines(mesh, framefield);
    polylines.Build();
    std::string filename1 = std::string(argv[1]).substr(0, std::string(argv[1]).size() - 4) + ".polyline.vtk";
    polylines.WriteFile(/*"PolyLines.vtk"*/filename1.c_str());

    return 0;
}
