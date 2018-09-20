/*
 * SLIM.cpp
 *
 *  Created on: Dec 15, 2016
 *      Author: cotrik
 */

#include <iostream>

#include <igl/slim.h>

#include <igl/components.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/Timer.h>

#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/harmonic.h>
#include <igl/serialize.h>
#include <igl/read_triangle_mesh.h>

#include <igl/flipped_triangles.h>
#include <igl/euler_characteristic.h>
#include <igl/barycenter.h>
#include <igl/adjacency_matrix.h>
#include <igl/is_edge_manifold.h>
#include <igl/doublearea.h>
#include <igl/cat.h>

#include "Mesh.h"
#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "ArgumentManager.h"
#include <stdlib.h>

#include <string>
#include <vector>
#include <map>

void BuildEnergyMap(std::map<std::string, igl::SLIMData::SLIM_ENERGY>& energyMap)
{
    energyMap["ARAP"]                    = igl::SLIMData::ARAP;
    energyMap["LOG_ARAP"]                = igl::SLIMData::LOG_ARAP;
    energyMap["SYMMETRIC_DIRICHLET"]     = igl::SLIMData::SYMMETRIC_DIRICHLET;
    energyMap["CONFORMAL"]               = igl::SLIMData::CONFORMAL;
    energyMap["EXP_CONFORMAL"]           = igl::SLIMData::EXP_CONFORMAL;
    energyMap["EXP_SYMMETRIC_DIRICHLET"] = igl::SLIMData::EXP_SYMMETRIC_DIRICHLET;
}

void GetVertices(const Eigen::MatrixXd& V_o, std::vector<Vertex>& V)
{
    V.resize(V_o.rows());
    for (size_t i = 0; i < V_o.rows(); i++) {
        Vertex& v = V.at(i);
        v.id = i;
        for (size_t j = 0; j < V_o.cols(); j++)
            v[j] = V_o(i, j);
    }
}
void GetCells(const Eigen::MatrixXi& F, std::vector<Cell>& C)
{
    C.resize(F.rows());
    for (size_t i = 0; i < F.rows(); i++) {
        Cell& c = C.at(i);
        c.id = i;
        c.Vids.resize(F.cols());
        for (size_t j = 0; j < F.cols(); j++)
            c.Vids[j] = F(i, j);
    }
}
void WriteVtk(const igl::SLIMData& sData, const ElementType cellType = TRIANGLE)
{
    std::vector<Vertex> V;
    GetVertices(sData.V_o, V);
    std::vector<Cell> C;
    GetCells(sData.F, C);
    static int iter = 0;
    std::string filename = std::string("iter.") + std::to_string(iter++) + ".vtk";
    MeshFileWriter writer(V, C, filename.c_str(), cellType);
    writer.WriteFile();
}
const size_t TetToHex[8] = {4, 5, 6, 7, 2, 3, 0, 1};
void GetHexCells(const Eigen::MatrixXi& F, std::vector<Cell>& C)
{
    C.resize(F.rows()/8);
    for (size_t i = 0; i < F.rows(); i+=8) {
        Cell& c = C.at(i/8);
        c.id = i/8;
        c.Vids.resize(8);
        for (size_t j = 0; j < 8; j++)
            c.Vids[TetToHex[j]] = F(i + j, 0);
    }
}
void WriteHexVtk(const igl::SLIMData& sData, const ElementType cellType = HEXAHEDRA)
{
    std::vector<Vertex> V;
    GetVertices(sData.V_o, V);
    std::vector<Cell> C;
    GetHexCells(sData.F, C);
    static int iter = 0;
    std::string filename = std::string("Hex.") + std::to_string(iter++) + ".vtk";
    MeshFileWriter writer(V, C, filename.c_str(), cellType);
    writer.WriteFile();
}

void GetVertices(const std::vector<Vertex>& V, Eigen::MatrixXd& V_o)
{
    V_o.resize(V.size(), 3);
    for (size_t i = 0; i < V_o.rows(); i++) {
        const Vertex& v = V.at(i);
        for (size_t j = 0; j < V_o.cols(); j++)
            V_o(i, j) = v[j];
    }
}

void GetCells(const std::vector<Cell>& C, Eigen::MatrixXi& F)
{
    F.resize(C.size(), C[0].Vids.size());
    for (size_t i = 0; i < F.rows(); i++) {
        const Cell& c = C.at(i);
        for (size_t j = 0; j < F.cols(); j++)
            F(i, j) = c.Vids[j];
    }
}
const size_t HexToTet[8][4] =
{
    {4, 6, 7, 3},
    {5, 7, 4, 0},
    {6, 4, 5, 1},
    {7, 5, 6, 2},
    {2, 0, 3, 7},
    {3, 1, 0, 4},
    {0, 2, 1, 5},
    {1, 3, 2, 6}
};

void ConvertHexToTet(const Eigen::MatrixXi& H, Eigen::MatrixXi& T)
{
    T.resize(H.rows() * 8, 4);
    for (size_t i = 0; i < H.rows(); i++)
        for (size_t j = 0; j < 8; j++)
            for (size_t k = 0; k < 4; k++)
                T(i * 8 + j, k) = H(i, HexToTet[j][k]);
}

size_t GetClosestBoundaryVertexId(const Vertex& v, const std::vector<Vertex>& V)
{
    size_t id = 0;
    double distance = MAXID;
    for (size_t i = 0; i < V.size(); i++) {
        const Vertex& tv = V.at(i);
        if (tv.isBoundary) {
            double d = glm::length(v - tv);
            if (d < distance) {
                distance = d;
                id = tv.id;
            }
        }
    }
    return id;
}
void GetPolycubeConstraints(const Mesh& tetOrigMesh, const Mesh& tetPolycubeMesh, const Mesh& hexPolycubeMesh,
                            Eigen::VectorXi& b, Eigen::MatrixXd& bc)
{
    int numOfVerticesOnBondary = 0;
    for (size_t i = 0; i < hexPolycubeMesh.V.size(); i++) {
        const Vertex& v = hexPolycubeMesh.V.at(i);
        if (v.isBoundary)
            numOfVerticesOnBondary++;
    }
    b.resize(numOfVerticesOnBondary);
    bc.resize(numOfVerticesOnBondary, 3);
    numOfVerticesOnBondary = 0;
    for (size_t i = 0; i < hexPolycubeMesh.V.size(); i++) {
        const Vertex& v = hexPolycubeMesh.V.at(i);
        if (v.isBoundary) {
            b(numOfVerticesOnBondary) = v.id;
            const size_t id = GetClosestBoundaryVertexId(v, tetPolycubeMesh.V);
            for (size_t j = 0; j < 3; j++)
                bc(numOfVerticesOnBondary, j) = tetOrigMesh.V[id][j];
            numOfVerticesOnBondary++;
        }
    }
}
void GetArguments(ArgumentManager& argumentManager,
    std::string& hex,// = "polycube.hex.vtk";
    std::string& tri,// = "polycube.tri.vtk";
    std::string& orig,// = "orig.tri.vtk";
    std::string& energy,// = "EXP_CONFORMAL";
    size_t& iters,// = 50;
    double& soft_const_p,// = 1e5;
    double& exp_factor// = 5.0;
)
{
    const std::string strHexPolycubeFilename = argumentManager.get("hex");
    if (!strHexPolycubeFilename.empty()) hex = strHexPolycubeFilename;

    const std::string strTriPoylycubeFilename = argumentManager.get("tri");
    if (!strTriPoylycubeFilename.empty()) tri = strTriPoylycubeFilename;

    const std::string strOrigTriFilename = argumentManager.get("orig");
    if (!strOrigTriFilename.empty()) orig = strOrigTriFilename;

    const std::string strEnergy = argumentManager.get("energy");
    if (!strEnergy.empty()) energy = strEnergy;

    const std::string strIters = argumentManager.get("iters");
    if (!strIters.empty()) iters = std::stoi(strIters);

    const std::string strSoftConstP = argumentManager.get("soft_const_p");
    if (!strSoftConstP.empty()) soft_const_p = std::stod(strSoftConstP);

    const std::string strExpFactor = argumentManager.get("exp_factor");
    if (!strExpFactor.empty()) exp_factor = std::stod(strExpFactor);

    std::cout << "-----------------------------------\n";
    std::cout << "hex = " << hex << std::endl;
    std::cout << "tri = " << tri << std::endl;
    std::cout << "orig = " << orig << std::endl;
    std::cout << "energy = " << energy << std::endl;
    std::cout << "iters = " << iters << std::endl;
    std::cout << "soft_const_p = " << soft_const_p << std::endl;
    std::cout << "exp_factor = " << exp_factor << std::endl;
    std::cout << "-----------------------------------\n";
}
int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cout << "Usage: Slim hex=<polycube.hex.vtk> tri=<polycube.tri.vtk> orig=<orig.tri.vtk> iters=<50> soft_const_p=<1e5> exp_factor=<5.0> "
                  << "energy=<ARAP|LOG_ARAP|SYMMETRIC_DIRICHLET|CONFORMAL|EXP_CONFORMAL|EXP_SYMMETRIC_DIRICHLET>\n\n";
        std::cout << "Example: Slim hex=\033[1;32mpolycube.hex.vtk\033[0m "
                  << "tri=\033[1;32mpolycube.tri.vtk\033[0m "
                  << "orig=\033[1;32morig.tri.vtk\033[0m "
                  << "iters=\033[1;32m50\033[0m "
                  << "soft_const_p=\033[1;32m1e5\033[0m "
                  << "exp_factor=\033[1;32m5.0\033[0m "
                  << "energy=\033[1;32mEXP_CONFORMAL\033[0m\n\n";
        return -1;
    }

    std::string hex = "polycube.hex.vtk";
    std::string tri = "polycube.tri.vtk";
    std::string orig = "orig.tri.vtk";
    std::string energy = "EXP_CONFORMAL";
    size_t iters = 50;
    double soft_const_p = 1e5;
    double exp_factor = 5.0;
    ArgumentManager argumentManager(argc, argv);
    GetArguments(argumentManager, hex, tri, orig, energy, iters, soft_const_p, exp_factor);

    //////////////////////////////////////
    MeshFileReader origReader(orig.c_str());
    Mesh& origMesh = (Mesh&)origReader.GetMesh();
    origMesh.BuildAllConnectivities();
    origMesh.ExtractBoundary();

    MeshFileReader triPolycubeReader(tri.c_str());
    Mesh& triPolycubeMesh = (Mesh&)triPolycubeReader.GetMesh();
    triPolycubeMesh.BuildAllConnectivities();
    triPolycubeMesh.ExtractBoundary();

    MeshFileReader hexPolycubeReader(hex.c_str());
    Mesh& hexPolycubeMesh = (Mesh&)hexPolycubeReader.GetMesh();
    hexPolycubeMesh.BuildAllConnectivities();
    hexPolycubeMesh.ExtractBoundary();

    Eigen::MatrixXd hexPolycubeV;
    Eigen::MatrixXi hexPolycubeH;
    GetVertices(hexPolycubeMesh.V, hexPolycubeV);
    GetCells(hexPolycubeMesh.C, hexPolycubeH);

    Eigen::MatrixXi hexPolycubeT;
    ConvertHexToTet(hexPolycubeH, hexPolycubeT);

    Eigen::MatrixXd V = hexPolycubeV;
    Eigen::MatrixXi F = hexPolycubeT;
    Eigen::MatrixXd V_0 = V;
    Eigen::VectorXi b;
    Eigen::MatrixXd bc;
    GetPolycubeConstraints(origMesh, triPolycubeMesh, hexPolycubeMesh, b, bc);

    igl::SLIMData sData;
    igl::Timer timer;
    sData.exp_factor = exp_factor;
    std::map<std::string, igl::SLIMData::SLIM_ENERGY> energyMap;
    BuildEnergyMap(energyMap);

    timer.start();
    slim_precompute(V, F, V_0, sData, energyMap[energy], b, bc, soft_const_p);
    std::cout << "precomputed time = " << timer.getElapsedTime() << std::endl;
    WriteVtk(sData, TETRAHEDRA);
    WriteHexVtk(sData, HEXAHEDRA);
    while (iters--) {
        timer.start();
        slim_solve(sData, 1); // 1 iter
        static int iter = 1;
        std::cout << "iter = " << iter++ << " time = " << timer.getElapsedTime() << std::endl;
        WriteVtk(sData, TETRAHEDRA);
        WriteHexVtk(sData, HEXAHEDRA);
    }
}
