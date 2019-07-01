/*
 * ExtractSharpEdges.cpp
 *
 *  Created on: Nov 17, 2017
 *      Author: cotrik
 */
#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "FrameField.h"
#include "PolyLine.h"
#include "FrameOpt.h"
#include "FeatureLine.h"
#include "MeshQuality.h"
#include "ArgumentManager.h"
#include <math.h>
#include <iostream>
#include <iomanip>

void WriteSharpEdgesVtk(const char* filename, const Mesh& m_mesh, const std::vector<FeatureLine>& featureLines);

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: ReadMesh <file> angle=<5|10|15|20|25|30>\n";
        return -1;
    }
    ArgumentManager argumentManager(argc, argv);
    double angle = 15;
    {
        const std::string strAngle = argumentManager.get("angle");
        if (!strAngle.empty()) angle = std::stod(strAngle);
    }
    MeshFileReader reader(argv[1]);
    Mesh& mesh = (Mesh&)reader.GetMesh();
    mesh.BuildAllConnectivities();
    mesh.ExtractBoundary();
    // cos10 = 0.984807753; cos15 = 0.965925826; cos20 = 0.939692621; cos25 = 0.906307787; cos30 = 0.866025404;
    mesh.SetCosAngleThreshold(cos(angle));
    mesh.LabelSurface();
    mesh.LabelSharpEdges(true);
    //mesh.ExtractLayers();
    mesh.ExtractSingularities();

//    MeshFileWriter edgesFileWriter(mesh, "Edges.vtk");
//    edgesFileWriter.WriteEdgesVtk();

    MeshFileWriter sharpEdgesFileWriter(mesh, "SharpEdges.vtk");
    sharpEdgesFileWriter.WriteSharpEdgesVtk();

//    MeshFileWriter facesFileWriter(mesh, "Faces.vtk");
//    facesFileWriter.WriteFacesVtk();

    MeshFileWriter cornersFileWriter(mesh, "Corners.vtk");
    cornersFileWriter.WriteCornersVtk();


    std::vector<FeatureLine> featureLines(mesh.numOfSharpEdges, FeatureLine(mesh));
    for (size_t i = 0; i < mesh.numOfSharpEdges; i++)
        featureLines.at(i).Extract(i);
    WriteSharpEdgesVtk("FeatureLines.vtk", mesh, featureLines);
    return 0;
}

void WriteSharpEdgesVtk(const char* filename, const Mesh& m_mesh, const std::vector<FeatureLine>& featureLines)
{
    const std::vector<Vertex>& V = m_mesh.V;
    const std::vector<Edge>& E = m_mesh.E;

    std::ofstream ofs(filename);
    ofs << "# vtk DataFile Version 2.0" << std::endl
        << filename << std::endl
        << "ASCII" << std::endl << std::endl
        << "DATASET UNSTRUCTURED_GRID" << std::endl;
    ofs << "POINTS " << V.size() << " float" << std::endl;
    ofs << std::fixed << std::setprecision(7);
    for (size_t i = 0; i < V.size(); i++)
        ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << std::endl;
    size_t numOfSharpVertices = 0;
    for (size_t i = 0; i < featureLines.size(); i++) {
        const FeatureLine& fl = featureLines.at(i);
        numOfSharpVertices += fl.Vids.size();
    }

    ofs << "CELLS " << featureLines.size() << " " << numOfSharpVertices + featureLines.size() << std::endl;
    for (size_t i = 0; i < featureLines.size(); i++) {
        const FeatureLine& fl = featureLines.at(i);
        ofs << fl.Vids.size();
        for (size_t j = 0; j < fl.Vids.size(); j++) {
            const size_t vid = fl.Vids.at(j);
            ofs << " " << vid;
        }
        ofs << std::endl;
    }

    ofs << "CELL_TYPES " << featureLines.size() << std::endl;
    for (size_t i = 0; i < featureLines.size(); i++) {
        ofs << 4 << std::endl;
    }

    ofs << "CELL_DATA " << featureLines.size() << std::endl
        << "SCALARS " << " Feature" << " int 1" << std::endl
        << "LOOKUP_TABLE default" << std::endl;
    for (size_t i = 0; i < featureLines.size(); i++) {
        ofs << i << std::endl;
    }
}
