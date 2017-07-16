/*
 * ReadMesh.cpp
 *
 *  Created on: Dec 12, 2016
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
#include <iostream>
#include <iomanip>
void test(const Mesh& mesh);
void WriteSharpEdgesVtk(const char* filename, const Mesh& m_mesh, const std::vector<FeatureLine>& featureLines)
{
    const std::vector<Vertex>& V = m_mesh.V;
    const std::vector<Edge>& E = m_mesh.E;

    std::ofstream ofs(filename);
    ofs << "# vtk DataFile Version 2.0" << std::endl
        << filename << std::endl
        << "ASCII" << std::endl << std::endl
        << "DATASET POLYDATA" << std::endl;
    ofs << "POINTS " << V.size() << " float" << std::endl;
    ofs << std::fixed << setprecision(7);
    for (size_t i = 0; i < V.size(); i++)
        ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << std::endl;
    size_t numOfSharpVertices = 0;
    for (size_t i = 0; i < featureLines.size(); i++) {
        const FeatureLine& fl = featureLines.at(i);
        numOfSharpVertices += fl.Vids.size();
    }

    ofs << "LINES " << featureLines.size() << " " << numOfSharpVertices + featureLines.size() << std::endl;
    for (size_t i = 0; i < featureLines.size(); i++) {
        const FeatureLine& fl = featureLines.at(i);
        ofs << fl.Vids.size();
        for (size_t j = 0; j < fl.Vids.size(); j++) {
            const size_t vid = fl.Vids.at(j);
            ofs << " " << vid << std::endl;
        }
    }

    ofs << "CELL_DATA " << featureLines.size() << std::endl
        << "SCALARS " << " Feature" << " int 1" << std::endl
        << "LOOKUP_TABLE default" << std::endl;
    for (size_t i = 0; i < featureLines.size(); i++) {
        ofs << i << std::endl;
    }
}
int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cout << "Usage: ReadMesh <file> cosangle=<0.939692621> \n"
                "Info: cos10 = 0.984807753; cos15 = 0.965925826; cos20 = 0.939692621; cos25 = 0.906307787; cos30 = 0.866025404\n\n";
        return -1;
    }
    ArgumentManager argumentManager(argc, argv);
    double cosangle = 0.939692621;
    {
        const std::string strCosangle = argumentManager.get("cosangle");
        if (!strCosangle.empty()) cosangle = std::stod(strCosangle);
    }
    MeshFileReader reader(argv[1]);
    Mesh& mesh = (Mesh&)reader.GetMesh();
    mesh.BuildAllConnectivities();
    mesh.ExtractBoundary();
    // cos10 = 0.984807753; cos15 = 0.965925826; cos20 = 0.939692621; cos25 = 0.906307787; cos30 = 0.866025404;
    mesh.SetCosAngleThreshold(cosangle);
    mesh.LabelSurface();
    mesh.LabelSharpEdges(true);
    //mesh.ExtractLayers();
    mesh.ExtractSingularities();

    test(mesh);

    MeshFileWriter edgesFileWriter(mesh, "Edges.vtk");
    edgesFileWriter.WriteEdgesVtk();

    MeshFileWriter sharpEdgesFileWriter(mesh, "SharpEdges.vtk");
    sharpEdgesFileWriter.WriteSharpEdgesVtk();

    MeshFileWriter facesFileWriter(mesh, "Faces.vtk");
    facesFileWriter.WriteFacesVtk();

    MeshFileWriter cornersFileWriter(mesh, "Corners.vtk");
    cornersFileWriter.WriteCornersVtk();

    MeshFileWriter tempFileWriter(mesh, "Temp.vtk");
    tempFileWriter.WriteFile();
    {
    std::vector<size_t> badCellIds;
    std::vector<size_t> warningCellIds;
    std::vector<size_t> goodCellIds;
    std::vector<size_t> highCellIds;
    std::vector<size_t> excellentCellIds;
    GetQuality("Temp.vtk", badCellIds, warningCellIds, goodCellIds, highCellIds, excellentCellIds);
    std::vector<Cell> badCells(badCellIds.size());
    for (size_t i = 0; i < badCellIds.size(); i++)
        badCells.at(i) = mesh.C.at(badCellIds.at(i));
    MeshFileWriter badCellsFileWriter(mesh.V, badCells, "badCells.vtk", HEXAHEDRA);
    badCellsFileWriter.WriteFile();
    std::vector<int> badCellLabel(badCellIds.size(), -1);
    badCellsFileWriter.WriteCellData(badCellLabel);

    std::vector<Cell> warningCells(warningCellIds.size());
    for (size_t i = 0; i < warningCellIds.size(); i++)
        warningCells.at(i) = mesh.C.at(warningCellIds.at(i));
    MeshFileWriter warningCellsFileWriter(mesh.V, warningCells, "warningCells.vtk", HEXAHEDRA);
    warningCellsFileWriter.WriteFile();
    std::vector<int> warningCellLabel(warningCellIds.size(), 0);
    warningCellsFileWriter.WriteCellData(warningCellLabel);

    std::vector<Cell> goodCells(goodCellIds.size());
    for (size_t i = 0; i < goodCellIds.size(); i++)
        goodCells.at(i) = mesh.C.at(goodCellIds.at(i));
    MeshFileWriter goodCellsFileWriter(mesh.V, goodCells, "goodCells.vtk", HEXAHEDRA);
    goodCellsFileWriter.WriteFile();
    std::vector<int> goodCellLabel(goodCellIds.size(), 1);
    goodCellsFileWriter.WriteCellData(goodCellLabel);

    std::vector<Cell> highCells(highCellIds.size());
    for (size_t i = 0; i < highCellIds.size(); i++)
        highCells.at(i) = mesh.C.at(highCellIds.at(i));
    MeshFileWriter highCellsFileWriter(mesh.V, highCells, "highCells.vtk", HEXAHEDRA);
    highCellsFileWriter.WriteFile();
    std::vector<int> highCellLabel(highCellIds.size(), 2);
    highCellsFileWriter.WriteCellData(highCellLabel);

    std::vector<Cell> excellentCells(excellentCellIds.size());
    for (size_t i = 0; i < excellentCellIds.size(); i++)
        excellentCells.at(i) = mesh.C.at(excellentCellIds.at(i));
    MeshFileWriter excellentCellsFileWriter(mesh.V, excellentCells, "excellentCells.vtk", HEXAHEDRA);
    excellentCellsFileWriter.WriteFile();
    std::vector<int> excellentCellLabel(excellentCellIds.size(), 3);
    excellentCellsFileWriter.WriteCellData(excellentCellLabel);
    }

//    MeshFileWriter surfaceFileWriter(mesh, "Surface.off");
//    surfaceFileWriter.WriteSurfaceOff();

    std::vector<FeatureLine> featureLines(mesh.numOfSharpEdges, FeatureLine(mesh));
    for (size_t i = 0; i < mesh.numOfSharpEdges; i++)
        featureLines.at(i).Extract(i);
    WriteSharpEdgesVtk("FeatureLines.vtk", mesh, featureLines);
    return 0;
}

void test(const Mesh& mesh)
{
    const size_t numOfV = mesh.V.size();
    const size_t numOfE = mesh.E.size();
    const size_t numOfF = mesh.F.size();
    const size_t numOfC = mesh.C.size();

    if (mesh.m_cellType == HEXAHEDRA || mesh.m_cellType == TETRAHEDRA) {
        std::cout << "genus = " <<  1 - (mesh.V.size() - mesh.E.size() + mesh.F.size() - mesh.C.size()) << std::endl;
        std::cout << "#V:" << mesh.V.size() << " - #E:" << mesh.E.size() << " + #F:" << mesh.F.size() << " - #C:" << mesh.C.size() << " = " << "1 - genus" << std::endl;
    }
    else if (mesh.m_cellType == TRIANGLE || mesh.m_cellType == QUAD) {
        std::cout << "genus = " <<  (2 - (mesh.V.size() - mesh.E.size() + mesh.F.size()))/2 << std::endl;
        std::cout << "#V:" << mesh.V.size() << " - #E:" << mesh.E.size() << " + #F:" << mesh.F.size() << " = " << "2 - 2*genus" << std::endl;
    }

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
    if (mesh.m_cellType == HEXAHEDRA || mesh.m_cellType == TETRAHEDRA) {
    std::cout << "N_Cids: ";
    for (size_t i = 0; i < v.N_Cids.size(); i++)
        std::cout << v.N_Cids[i] << " ";
    std::cout << "\n";
    }
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
    if (mesh.m_cellType == HEXAHEDRA || mesh.m_cellType == TETRAHEDRA) {
    std::cout << "N_Cids: ";
    for (size_t i = 0; i < e.N_Cids.size(); i++)
        std::cout << e.N_Cids[i] << " ";
    std::cout << "\n";
    }
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
    if (mesh.m_cellType == HEXAHEDRA || mesh.m_cellType == TETRAHEDRA) {
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

    const glm::vec3 a(1, 2, 1);
    const glm::vec3 dir(1, 2, 3);
    Line line(a, dir);
    const glm::vec3 p(2, 3, 4);
    glm::vec3 intersection;
    const double d = line.Perpendicular(p, intersection);
    std::cout << "d = " << d << std::endl;
}


