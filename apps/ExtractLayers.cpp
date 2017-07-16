/*
 * ExtractLayers.cpp
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

void test(const Mesh& mesh);

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cout << "Usage: ExtractLayers <file>\n";
        return -1;
    }
    MeshFileReader reader(argv[1]);
    Mesh& mesh = (Mesh&)reader.GetMesh();
    mesh.BuildAllConnectivities();
    mesh.ExtractBoundary();
    mesh.ExtractLayers();
    mesh.ExtractSingularities();
    std::cout << "genus = " <<  1 - (mesh.V.size() - mesh.E.size() + mesh.F.size() - mesh.C.size()) << std::endl;
    std::cout << "#V:" << mesh.V.size() << " - #E:" << mesh.E.size() << " + #F:" << mesh.F.size() << " - #C:" << mesh.C.size() << " = " << "1 - genus" << std::endl;
    test(mesh);

    MeshFileWriter edgesFileWriter(mesh, "Edges.vtk");
    edgesFileWriter.WriteEdgesVtk();

    MeshFileWriter facesFileWriter(mesh, "Faces.vtk");
    facesFileWriter.WriteFacesVtk();

    MeshFileWriter surfaceFileWriter(mesh, "Surface.off");
    surfaceFileWriter.WriteSurfaceOff();

    return 0;
}

void test(const Mesh& mesh)
{
    const size_t numOfV = mesh.V.size();
    const size_t numOfE = mesh.E.size();
    const size_t numOfF = mesh.F.size();
    const size_t numOfC = mesh.C.size();
    std::cout << "#V:" << mesh.V.size() << " - #E:" << mesh.E.size() << " + #F:" << mesh.F.size() << " - #C:" << mesh.C.size() << " = " << "1 - genus" << std::endl;

    const Vertex& v = mesh.V.at(0);
    std::cout << "isBoundary: " << v.isBoundary << std::endl;
    std::cout << "---- Vertex 0 ----\n";
    std::cout << "N_Vids: ";
    for (size_t i = 0; i < v.N_Vids.size(); i++)
        std::cout << v.N_Vids[i] << " ";
    std::cout << "\n";

    std::cout << "N_Eids: ";
    for (size_t i = 0; i < v.N_Eids.size(); i++)
        std::cout << v.N_Eids[i] << " ";
    std::cout << "\n";

    std::cout << "N_Fids: ";
    for (size_t i = 0; i < v.N_Fids.size(); i++)
        std::cout << v.N_Fids[i] << " ";
    std::cout << "\n";

    std::cout << "N_Cids: ";
    for (size_t i = 0; i < v.N_Cids.size(); i++)
        std::cout << v.N_Cids[i] << " ";
    std::cout << "\n";

    const Edge& e = mesh.E.at(0);
    std::cout << "isBoundary: " << e.isBoundary << std::endl;
    std::cout << "---- Edge 0 ----\n";
    std::cout << "Vids: ";
    for (size_t i = 0; i < e.Vids.size(); i++)
        std::cout << e.Vids[i] << " ";
    std::cout << "\n";

    std::cout << "N_Vids: ";
    for (size_t i = 0; i < e.N_Vids.size(); i++)
        std::cout << e.N_Vids[i] << " ";
    std::cout << "\n";

    std::cout << "N_Eids: ";
    for (size_t i = 0; i < e.N_Eids.size(); i++)
        std::cout << e.N_Eids[i] << " ";
    std::cout << "\n";

    std::cout << "N_Fids: ";
    for (size_t i = 0; i < e.N_Fids.size(); i++)
        std::cout << e.N_Fids[i] << " ";
    std::cout << "\n";

    std::cout << "N_Cids: ";
    for (size_t i = 0; i < e.N_Cids.size(); i++)
        std::cout << e.N_Cids[i] << " ";
    std::cout << "\n";

    const Face& f = mesh.F.at(0);
    std::cout << "---- Face 0 ----\n";
    std::cout << "isBoundary: " << f.isBoundary << std::endl;

    std::cout << "Vids: ";
    for (size_t i = 0; i < f.Vids.size(); i++)
        std::cout << f.Vids[i] << " ";
    std::cout << "\n";

    std::cout << "Eids: ";
    for (size_t i = 0; i < f.Eids.size(); i++)
        std::cout << f.Eids[i] << " ";
    std::cout << "\n";

    std::cout << "N_Vids: ";
    for (size_t i = 0; i < f.N_Vids.size(); i++)
        std::cout << f.N_Vids[i] << " ";
    std::cout << "\n";

    std::cout << "N_Eids: ";
    for (size_t i = 0; i < f.N_Eids.size(); i++)
        std::cout << f.N_Eids[i] << " ";
    std::cout << "\n";

    std::cout << "N_Fids: ";
    for (size_t i = 0; i < f.N_Fids.size(); i++)
        std::cout << f.N_Fids[i] << " ";
    std::cout << "\n";

    std::cout << "N_Cids: ";
    for (size_t i = 0; i < f.N_Cids.size(); i++)
        std::cout << f.N_Cids[i] << " ";
    std::cout << "\n";

    const Cell& c = mesh.C.at(0);
    std::cout << "---- Cell 0 ----\n";
    std::cout << "isBoundary: " << c.isBoundary << std::endl;

    std::cout << "Vids: ";
    for (size_t i = 0; i < c.Vids.size(); i++)
        std::cout << c.Vids[i] << " ";
    std::cout << "\n";

    std::cout << "Eids: ";
    for (size_t i = 0; i < c.Eids.size(); i++)
        std::cout << c.Eids[i] << " ";
    std::cout << "\n";

    std::cout << "Fids: ";
    for (size_t i = 0; i < c.Fids.size(); i++)
        std::cout << c.Fids[i] << " ";
    std::cout << "\n";

    std::cout << "N_Vids: ";
    for (size_t i = 0; i < c.N_Vids.size(); i++)
        std::cout << c.N_Vids[i] << " ";
    std::cout << "\n";

    std::cout << "N_Eids: ";
    for (size_t i = 0; i < c.N_Eids.size(); i++)
        std::cout << c.N_Eids[i] << " ";
    std::cout << "\n";

    std::cout << "N_Fids: ";
    for (size_t i = 0; i < c.N_Fids.size(); i++)
        std::cout << c.N_Fids[i] << " ";
    std::cout << "\n";

    std::cout << "N_Cids: ";
    for (size_t i = 0; i < c.N_Cids.size(); i++)
        std::cout << c.N_Cids[i] << " ";
    std::cout << "\n";
}
