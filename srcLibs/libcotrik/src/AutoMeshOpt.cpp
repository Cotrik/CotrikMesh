/*
 * AutoMeshOpt.cpp
 *
 *  Created on: Nov 25, 2016
 *      Author: cotrik
 */

#include "AutoMeshOpt.h"
#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "MeshQuality.h"

#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <queue>
#include <iostream>
#include <math.h>
#include <igl/hausdorff.h>

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
void WriteHexVtk(const igl::SLIMData& sData, const ElementType cellType = HEXAHEDRA, const char* pFilename = NULL)
{
    std::vector<Vertex> V;
    GetVertices(sData.V_o, V);
    std::vector<Cell> C;
    GetHexCells(sData.F, C);
    static int iter = 0;
    if (pFilename == NULL) {
        std::string filename = std::string("Hex.") + std::to_string(iter++) + ".vtk";
        MeshFileWriter writer(V, C, filename.c_str(), cellType);
        writer.WriteFile();
    } else {
        MeshFileWriter writer(V, C, pFilename, cellType);
        writer.WriteFile();
    }
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
void SetConstraints(const Mesh& targetMesh, Eigen::VectorXi& b, Eigen::MatrixXd& bc)
{
    int numOfVerticesOnBondary = 0;
    for (size_t i = 0; i < targetMesh.V.size(); i++) {
        const Vertex& v = targetMesh.V.at(i);
        if (v.isBoundary)
            numOfVerticesOnBondary++;
    }
    b.resize(numOfVerticesOnBondary);
    bc.resize(numOfVerticesOnBondary, 3);
    numOfVerticesOnBondary = 0;
    for (size_t i = 0; i < targetMesh.V.size(); i++) {
        const Vertex& v = targetMesh.V.at(i);
        if (v.isBoundary) {
            b(numOfVerticesOnBondary) = v.id;
            for (size_t j = 0; j < 3; j++)
                bc(numOfVerticesOnBondary, j) = targetMesh.V[i][j];
            numOfVerticesOnBondary++;
        }
    }
}
void AutoMeshOpt::InversionFreeDeformToTargetMesh()
{
    Eigen::MatrixXd inputV;
    Eigen::MatrixXi inputH;
    GetVertices(mesh.V, inputV);
    GetCells(mesh.C, inputH);

    Eigen::MatrixXi inputT;
    ConvertHexToTet(inputH, inputT);

    Eigen::MatrixXd V = inputV;
    Eigen::MatrixXi F = inputT;
    Eigen::MatrixXd V_0 = V;
    Eigen::VectorXi b;
    Eigen::MatrixXd bc;
    SetConstraints(origMesh, b, bc);

    std::string energy = "EXP_CONFORMAL";
    size_t iters = 30;
    double soft_const_p = 1e5;
    double exp_factor = 5.0;

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
    while (iters-- != 0) {
        timer.start();
        slim_solve(sData, 1); // 1 iter
        static int iter = 1;
        std::cout << "iter = " << iter++ << " time = " << timer.getElapsedTime() << std::endl;
        WriteVtk(sData, TETRAHEDRA);
        WriteHexVtk(sData, HEXAHEDRA);
//        if (iters == 0)
//            WriteHexVtk(sData, HEXAHEDRA, result.c_str());
        GetVertices(sData.V_o, mesh.V);
        if (GetHausdorffError() < m_hausdorffError)
            return;
    }
}

AutoMeshOpt::AutoMeshOpt(const Mesh& mesh, const Mesh& origMesh)
: LocalMeshOpt(mesh)
, origMesh(origMesh)
, m_hausdorffError(1e-2)
{
    // TODO Auto-generated constructor stub

}

AutoMeshOpt::~AutoMeshOpt()
{
    // TODO Auto-generated destructor stub
}

void GetSurfacVertices(const std::vector<Vertex>& V, Eigen::MatrixXd& V_o)
{
    size_t numOfSurfaceVertices = 0;
    for (size_t i = 0; i < V.size(); i++) {
        const Vertex& v = V.at(i);
        if (v.isBoundary)
            numOfSurfaceVertices++;
    }
    V_o.resize(numOfSurfaceVertices, 3);
    numOfSurfaceVertices = 0;
    for (size_t i = 0; i < V.size(); i++) {
        const Vertex& v = V.at(i);
        if (v.isBoundary) {
            for (size_t j = 0; j < V_o.cols(); j++)
                V_o(numOfSurfaceVertices, j) = v[j];
            numOfSurfaceVertices++;
        }
    }
}

void GetSurfaceFaces(const std::vector<Face>& F, Eigen::MatrixXi& F_o)
{
    size_t numOfSurfaceFaces = 0;
    for (size_t i = 0; i < F.size(); i++) {
        const Face& f = F.at(i);
        if (f.isBoundary)
            numOfSurfaceFaces++;
    }
    F_o.resize(numOfSurfaceFaces, F[0].Vids.size());
    numOfSurfaceFaces = 0;
    for (size_t i = 0; i < F.size(); i++) {
        const Face& f = F.at(i);
        if (f.isBoundary) {
            for (size_t j = 0; j < F_o.cols(); j++)
                F_o(numOfSurfaceFaces, j) = f.Vids[j];
            numOfSurfaceFaces++;
        }
    }
}

void AutoMeshOpt::Run()
{

    // 1. Smooth and map the surface to origin triangle surface. Obtain target surface.
    // 2. Use a large parameters to make hard constraint on surface to untangle the mesh.
    // 3. If the mesh is not untangled, half decrease the parameters until the mesh is untangled.
    // 4. Use inversion free deformation (SLIM, etc.) to get the target surface
    // 5. Use a large parameters to make hard constraint on surface until the mesh reaches the minimum Scaled Jacobian (plus 0.05 every time)
    // 6. If the mesh does not reach the minimum Scaled Jacobian, half decrease the parameters until the mesh reaches the minimum Scaled Jacobian.
    //    Of course we have to guarantee the target Hausdorff Distance error. If the error > target, map the surface to origin triangle surface.
    //    If the error still > target or quality decreases, stop and output the hexahedral mesh. Otherwise go to step 5.

    SetAlpha(1000);
    SetBeta(100);
    SetGamma(1.0);
    SetStepSize(0.9);
    SetUseAverageTargetLength(false);
    SetAnisotropy(0.2);
    SetMinScaledJacobian(0.0);
    SetRecoverable(true);
    SetAllowBigStep(true);
    //while (true)
    {
        while (true) {
            LocalMeshOpt::Run(10, 25);
            double minSJ = 0;
            std::vector<size_t> badCellIds;
            m_numOfInvertdElements = GetMinScaledJacobian(mesh, minSJ, badCellIds, this->minScaledJacobian);
            if (m_numOfInvertdElements == 0) {
                InversionFreeDeformToTargetMesh();
                break;
            }
            else {
                this->alpha /= 2;
                this->beta /= 2;
            }
        }
    }

    //while (true)
    {
        while (true) {
            LocalMeshOpt::Run(10, 25);
            double minSJ = 0;
            std::vector<size_t> badCellIds;
            m_numOfInvertdElements = GetMinScaledJacobian(mesh, minSJ, badCellIds, this->minScaledJacobian);
            if (m_numOfInvertdElements == 0) {
                if (GetHausdorffError() < m_hausdorffError) {
                    std::string filename = std::string("MSJ=") + std::to_string(minSJ) + ".vtk";
                    MeshFileWriter writer(mesh, filename.c_str());
                    writer.WriteFile();
                    this->minScaledJacobian += 0.05;
                }
                else {
                    mesh.ProjectTo(origMesh);
                    m_numOfInvertdElements = GetMinScaledJacobian(mesh, minSJ, badCellIds, this->minScaledJacobian);
                    if (m_numOfInvertdElements == 0 && minSJ >= this->minScaledJacobian) {
                        if (GetHausdorffError() < m_hausdorffError) {
                            std::string filename = std::string("MSJ=") + std::to_string(minSJ) + ".vtk";
                            MeshFileWriter writer(mesh, filename.c_str());
                            writer.WriteFile();
                            this->minScaledJacobian += 0.05;
                            continue;
                        }
                    }
                    break;
                }
            }
            else {
                this->alpha /= 2;
                this->beta /= 2;
                if (this->beta < 2 || this->alpha < 2)
                    break;
            }
        }
    }
}

double AutoMeshOpt::GetHausdorffError() const
{
    double local_hausdorffErr = 0;

    Eigen::MatrixXd VA;
    Eigen::MatrixXi FA;
    Eigen::MatrixXd VB;
    Eigen::MatrixXi FB;
    GetSurfacVertices(mesh.V, VA);
    GetSurfaceFaces(mesh.F, FA);
    GetSurfacVertices(origMesh.V, VB);
    GetSurfaceFaces(origMesh.F, FB);
    igl::hausdorff(VA, FA, VB, FB, local_hausdorffErr);

    return local_hausdorffErr;
}
void AutoMeshOpt::SetHausdorffError(const double value/* = 1e-2*/)
{
    m_hausdorffError = value;
}
