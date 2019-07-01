/*
 * PolycubeProcess.cpp
 *
 *  Created on: Nov 5, 2018
 *      Author: cotrik
 */

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "PolyLine.h"
#include "FeatureLine.h"
#include "ArgumentManager.h"
#include <iostream>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <stack>
#include <algorithm>
#include <iomanip>

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cout << "Usage: PolycubeProcess origTet.vtk polycubeTet.vtk angle=<5|10|15|20|25|30>\n";
        return -1;
    }
    ArgumentManager argumentManager(argc, argv);
    std::string origTetFilename = argv[1];
    std::string polycubeTetFilename = argv[2];
    double angle = 15;
    {
        const std::string strAngle = argumentManager.get("angle");
        if (!strAngle.empty()) angle = std::stod(strAngle);
    }
    MeshFileReader origTetMeshReader(origTetFilename.c_str());
    Mesh& origTetMesh = (Mesh&) origTetMeshReader.GetMesh();
    origTetMesh.BuildAllConnectivities();
    origTetMesh.ExtractBoundary();
    // cos10 = 0.984807753; cos15 = 0.965925826; cos20 = 0.939692621; cos25 = 0.906307787; cos30 = 0.866025404;
    origTetMesh.SetCosAngleThreshold(cos(angle));
    origTetMesh.LabelSurface();
    origTetMesh.LabelSharpEdges(true);
    origTetMesh.ExtractSingularities();

    MeshFileReader polycubeTetMeshReader(polycubeTetFilename.c_str());
    Mesh& polycubeTetMesh = (Mesh&) polycubeTetMeshReader.GetMesh();
    polycubeTetMesh.BuildAllConnectivities();
    polycubeTetMesh.ExtractBoundary();
    // cos10 = 0.984807753; cos15 = 0.965925826; cos20 = 0.939692621; cos25 = 0.906307787; cos30 = 0.866025404;
    polycubeTetMesh.SetCosAngleThreshold(cos(angle));
    polycubeTetMesh.LabelSurface();
    polycubeTetMesh.LabelSharpEdges(true);
    polycubeTetMesh.ExtractSingularities();

//    MeshFileWriter sharpEdgesFileWriter(mesh, "SharpEdges.vtk");
//    sharpEdgesFileWriter.WriteSharpEdgesVtk();
//
//    MeshFileWriter cornersFileWriter(mesh, "Corners.vtk");
//    cornersFileWriter.WriteCornersVtk();

    std::vector<FeatureLine> featureLines(polycubeTetMesh.numOfSharpEdges, FeatureLine(polycubeTetMesh));
    for (size_t i = 0; i < polycubeTetMesh.numOfSharpEdges; i++)
        featureLines.at(i).Extract(i);
//    {
//        MeshFileWriter writer(polycubeTetMesh, "FeatureLines.vtk");
//        writer.WriteFeatureLinesVtk(featureLines);
//    }
    {
        std::vector<size_t> cornerVids;
        for (auto& v : polycubeTetMesh.V)
            if (v.isCorner) cornerVids.push_back(v.id);
        MeshFileWriter writer(origTetMesh, "OrigTetMeshWithFeatures.vtk");
        //writer.WriteCellsWithFeatureLinesVtk(featureLines);
        writer.WriteCellsWithFeatureVtk(featureLines, cornerVids);
    }
    return 0;
}
