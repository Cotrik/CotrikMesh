/*
 * LocalMeshOpt.cpp
 *
 *  Created on: Nov 25, 2016
 *      Author: cotrik
 */

#include "LocalMeshOpt.h"
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

//LocalMeshOpt::LocalMeshOpt()
LocalMeshOpt::LocalMeshOpt(const Mesh& mesh)
: MeshOpt(mesh)
{
    // TODO Auto-generated constructor stub

}

LocalMeshOpt::~LocalMeshOpt()
{
    // TODO Auto-generated destructor stub
}


double GetAvgEdgeLength(Mesh& mesh)
{
    double sumEdgeLength = 0.0;
    size_t numOfBoundaryEdges = 0;
    for (size_t i = 0; i < mesh.E.size(); i++)
    {
        if (mesh.E.at(i).isBoundary){
            const Vertex& v1 = mesh.V.at(mesh.E.at(i).Vids.at(0));
            const Vertex& v2 = mesh.V.at(mesh.E.at(i).Vids.at(1));
            mesh.E[i].length = glm::length(v1 - v2);
            sumEdgeLength += mesh.E[i].length;
            numOfBoundaryEdges++;
        }
    }
    return sumEdgeLength / numOfBoundaryEdges;
}

//bool LocalMeshOpt::UntangleLocalMesh(Mesh& localMesh, const size_t localIters/* = 20*/)
//{
//    MeshOpt meshOpt(localMesh);
//    meshOpt.SetAlpha(alpha);
//    meshOpt.SetBeta(beta);
//    meshOpt.SetGamma(gamma);
//    meshOpt.SetStepSize(stepSize);
//    meshOpt.SetUseAverageTargetLength(useAverageTargetLength);
//    meshOpt.SetAnisotropy(anisotropy);
//    meshOpt.SetMinScaledJacobian(minScaledJacobian);
//    meshOpt.SetRecoverable(recoverable);
//    meshOpt.SetAllowBigStep(allowBigStep);
//    meshOpt.SetChangeBoundary(changeBoundary);
//    size_t iter = meshOpt.Run(localIters);
//
//    std::string filename = std::string("MeshOpt.") + std::to_string((int)iter) + ".vtk";
//    double minimumScaledJacobian;
//    size_t numOfInvertedElements = GetScaledJacobianVerdict(localMesh, minimumScaledJacobian, this->minScaledJacobian);
//
//    static int count = 0;
//    count++;
//    if (numOfInvertedElements < m_numOfInvertdElements) {
//        MeshFileReader reader(filename.c_str());
//        Mesh& localOptMesh = (Mesh&)reader.GetMesh();
//        filename = std::string("LocalOpt.") + std::to_string(count) + ".vtk";
//        MeshFileWriter writerT(localOptMesh, filename.c_str());
//        writerT.WriteFile();
//        return true;
//    }
//    return false;
//}

void LocalMeshOpt::Untangle(const std::vector<size_t>& cellIds, const size_t localIters/* = 20*/)
{
    Mesh localMesh(mesh, cellIds);
    localMesh.RemoveUselessVertices();
    localMesh.BuildAllConnectivities();
    localMesh.ExtractBoundary();
    localMesh.ExtractSingularities();
    localMesh.SetCosAngleThreshold(0.984807753);
    localMesh.LabelSurface();
    localMesh.BuildParallelE();
    localMesh.BuildConsecutiveE();
    localMesh.BuildOrthogonalE();
    localMesh.GetNormalOfSurfaceFaces();
    localMesh.GetNormalOfSurfaceVertices();
    localMesh.ExtractTwoRingNeighborSurfaceFaceIdsForEachVertex(2);
    localMesh.GetAvgEdgeLength();

    Mesh bestMesh;
    if (UntangleLocalMesh(localMesh, bestMesh, localIters)) {
        Mesh& localOptMesh = bestMesh;
        ModifyMeshFrom(localOptMesh, cellIds);
    }
}

bool LocalMeshOpt::UntangleLocalMesh(Mesh& localMesh, Mesh& bestMesh, const size_t localIters/* = 20*/)
{
    MeshOpt meshOpt(localMesh);
    meshOpt.SetAlpha(alpha);
    meshOpt.SetBeta(beta);
    meshOpt.SetGamma(gamma);
    meshOpt.SetStepSize(stepSize);
    meshOpt.SetUseAverageTargetLength(useAverageTargetLength);
    meshOpt.SetAnisotropy(anisotropy);
    meshOpt.SetMinScaledJacobian(minScaledJacobian);
    meshOpt.SetRecoverable(recoverable);
    meshOpt.SetAllowBigStep(allowBigStep);
    meshOpt.SetChangeBoundary(changeBoundary);
    meshOpt.SetRefMesh(mesh);
    meshOpt.SetTargetSurfaceMesh(*m_targetSurfaceMesh);
    meshOpt.SetProjectToTargetSurface(projectToTargetSurface);
    size_t iter = meshOpt.Run(localIters);

    double minimumScaledJacobian = 0.0;
    size_t numOfInvertedElements = GetMinScaledJacobianVerdict(meshOpt.bestmesh, minimumScaledJacobian, this->minScaledJacobian);
    if (numOfInvertedElements < m_numOfInvertdElements) {
        bestMesh = meshOpt.bestmesh;
        return true;
    }
    return false;
}

void LocalMeshOpt::ModifyMeshFrom(const Mesh& localMesh, const std::vector<size_t>& badCellIdsT)
{
    std::vector<size_t> v_real_index;
    size_t c_size = 0;
    for (size_t i = 0; i < badCellIdsT.size(); i++) {
        const Cell& cell = mesh.C.at(badCellIdsT.at(i));
        for (size_t j = 0; j < cell.Vids.size(); j++)
            v_real_index.push_back(cell.Vids.at(j));
        c_size++;
    }

    std::sort(v_real_index.begin(), v_real_index.end());
    std::vector<size_t>::iterator iter = std::unique(v_real_index.begin(), v_real_index.end());
    v_real_index.resize(std::distance(v_real_index.begin(), iter));

    for (size_t i = 0; i < v_real_index.size(); i++) {
        Vertex& v = mesh.V.at(v_real_index.at(i));
        const Vertex& newv = localMesh.V.at(i);
        v.x = newv.x;
        v.y = newv.y;
        v.z = newv.z;
    }
}

void LocalMeshOpt::Run(const size_t iters/* = 1*/, const size_t localIters/* = 20*/)
{
    int iter = 0;
    while (iter++ != iters) {
        double minimumScaledJacobian = 0.0;
        std::vector<size_t> badCellIds;
        m_numOfInvertdElements = GetMinScaledJacobianVerdict(mesh, minimumScaledJacobian, badCellIds, this->minScaledJacobian);
        if (m_numOfInvertdElements == 0) break;

        if (!useSmallBlock) {
            std::vector<size_t> cellIds = GetBadCellsAndExtendCells(badCellIds, blockSize);
            Untangle(cellIds, localIters);
        }
        else {
            std::vector<std::vector<size_t> > regions;
            DivideIntoMultipleRegions(badCellIds, regions, blockSize);
            std::cout << "#badCellIds = " << badCellIds.size() << " #regions = " << regions.size() << "\n";
            if (parallel) {
#pragma omp parallel for
                for (int i = 0; i < regions.size(); i++) {
                    std::vector<size_t> cellIds = GetBadCellsAndExtendCells(regions[i]);
                    Untangle(cellIds, localIters);
                }
            } else {
                for (int i = 0; i < regions.size(); i++) {
                    std::vector<size_t> cellIds = GetBadCellsAndExtendCells(regions[i]);
                    Untangle(cellIds, localIters);
                }
            }
        }
    }

    std::string filename = std::string("BestLocalOpt.vtk");
    MeshFileWriter writerT(mesh, filename.c_str());
    writerT.WriteFile();
}

//void LocalMeshOpt::Run(const size_t iters/* = 1*/, const size_t localIters/* = 20*/)
//{
//    int count = iters;
//    int iter = 0;
//    std::string filename;
//    while (count-- != 0)
//    {
//        double minimumScaledJacobian = 0.0;
//        std::vector<size_t> badCellIds;
//        m_numOfInvertdElements = GetMinScaledJacobianVerdict(mesh, minimumScaledJacobian, badCellIds, this->minScaledJacobian);
//        if (m_numOfInvertdElements == 0) break;
//        GetLocalRegions(badCellIds, 2);
//        UntangleLocalRegions(localIters);
//        ++iter;
//    }
//
//    filename = std::string("BestLocalOpt.vtk");
//    MeshFileWriter writerT(mesh, filename.c_str());
//    writerT.WriteFile();
//}

void LocalMeshOpt::GetLocalRegions(const std::vector<size_t>& badCellIds, int extendLayers/* = 2*/)
{
    m_localCellIds = GetBadCellsAndExtendCells(badCellIds, extendLayers);
    m_localVertexIds.clear();
    for (size_t i = 0; i < m_localCellIds.size(); i++) {
        const Cell& cell = mesh.C.at(m_localCellIds.at(i));
        for (size_t j = 0; j < cell.Vids.size(); j++)
            m_localVertexIds.push_back(cell.Vids.at(j));
    }

    std::sort(m_localVertexIds.begin(), m_localVertexIds.end());
    std::vector<size_t>::iterator iter = std::unique(m_localVertexIds.begin(), m_localVertexIds.end());
    m_localVertexIds.resize(std::distance(m_localVertexIds.begin(), iter));
}

void LocalMeshOpt::DivideIntoMultipleRegions(const std::vector<size_t>& badCellIds, std::vector<std::vector<size_t> >& regions, const int N)
{
    //const int N = 50;
    const size_t n = badCellIds.size();
    for (size_t i = 0; i < 1000000; i++) {
        size_t begin = i * N;
        size_t end = (i + 1) * N;
        if (n > begin && n <= end) {
            end = n;
        }
        std::vector<size_t> region(end - begin);
        std::copy(badCellIds.begin() + begin, badCellIds.begin() + end, region.begin());
        regions.push_back(region);
        if (end != (i + 1) * N || end == n) break;
    }
}

void LocalMeshOpt::SetUseSmallBlock(bool value/* = true*/)
{
    useSmallBlock = value;
}

void LocalMeshOpt::SetBlockSize(size_t value/* = 50*/)
{
    blockSize = value;
}

void LocalMeshOpt::SetParallel(bool value/* = true*/)
{
    parallel = value;
}
