/*
 * LocalMeshOptASJ.cpp
 *
 *  Created on: Nov 25, 2016
 *      Author: cotrik
 */

#include "LocalMeshOptASJ.h"
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
LocalMeshASJOpt::LocalMeshASJOpt(const Mesh& mesh)
: MeshASJOpt(mesh)
{
    // TODO Auto-generated constructor stub

}

LocalMeshASJOpt::~LocalMeshASJOpt()
{
    // TODO Auto-generated destructor stub
}

bool LocalMeshASJOpt::UntangleLocalMesh(Mesh& localMesh, const size_t localIters/* = 20*/)
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
    size_t iter = meshOpt.Run(localIters);

    double minimumScaledJacobian = 0.0;
    std::vector<size_t> badCellIds;
    size_t numOfInvertedElements = GetMinScaledJacobianVerdict(meshOpt.bestmesh, minimumScaledJacobian, badCellIds, this->minScaledJacobian);
    static int count = 0;
    count++;
    if (numOfInvertedElements < m_numOfInvertdElements) {
        bestmesh = meshOpt.bestmesh;
        return true;
    }
    return false;
}

bool LocalMeshASJOpt::UntangleLocalMesh(Mesh& localMesh, Mesh& bestMesh, const size_t localIters/* = 20*/)
{
    MeshASJOpt meshOpt(localMesh);
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
    size_t iter = meshOpt.Run(localIters);

    double minimumScaledJacobian = 0.0;
    size_t numOfInvertedElements = GetMinScaledJacobianVerdict(meshOpt.bestmesh, minimumScaledJacobian, this->minScaledJacobian);
    static int count = 0;
    count++;
    if (numOfInvertedElements < m_numOfInvertdElements) {
        bestMesh = meshOpt.bestmesh;
        return true;
    }
    return false;
}

void LocalMeshASJOpt::ModifyMeshFrom(const Mesh& localMesh, const std::vector<size_t>& badCellIdsT)
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
struct IntBool {
    IntBool(): id(0), exist(true) {

    }
    size_t id;
    bool exist;
    bool operator == (const IntBool& rhs) {
        return (id == rhs.id);
    }
};

void LocalMeshASJOpt::DivideIntoMultipleRegions(const std::vector<size_t>& badCellIds, std::vector<std::vector<size_t> >& regions)
{
    const int N = 50;
    const size_t n = badCellIds.size();
    for (size_t i = 0; i < 1000000; i++) {
        size_t begin = i*N;
        size_t end = (i+1)*N;
        if (n > begin && n <= end) {
            end = n;
        }
        std::vector<size_t> region(end - begin);
        std::copy(badCellIds.begin() + begin, badCellIds.begin() + end, region.begin());
        regions.push_back(region);
        if (end != (i+1)*N)
            break;
    }
}

void LocalMeshASJOpt::Run(const size_t iters/* = 1*/, const size_t localIters/* = 20*/)
{
    int count = iters;
    int iter = 0;
    std::string filename;
    while (count-- != 0)
    {
        double minimumScaledJacobian = 0.0;  std::vector<size_t> badCellIds;
        m_numOfInvertdElements = GetMinScaledJacobianVerdict(mesh, minimumScaledJacobian, badCellIds, this->minScaledJacobian);
        if (m_numOfInvertdElements == 0)
            break;
        std::vector<std::vector<size_t> > regions;
        DivideIntoMultipleRegions(badCellIds, regions);
        std::cout << "#badCellIds = " << badCellIds.size() << " #regions = " << regions.size() << "\n";
        ++iter;

        //MeshFileReader reader(filename.c_str());
        //Mesh& localMesh = (Mesh&)reader.GetMesh();
//#pragma omp parallel for
        for (int i = 0; i < regions.size(); i++)
        {
            std::vector<size_t> cellIds = GetBadCellsAndExtendCells(regions[i]);
            //cellIds = regions[0];
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

            Mesh bestMesh;
            if (UntangleLocalMesh(localMesh, bestMesh, localIters)){
                Mesh& localOptMesh = bestMesh;
                ModifyMeshFrom(localOptMesh, cellIds);
            }
        }
    }

    filename = std::string("BestLocalOpt.vtk");
    MeshFileWriter writerT(mesh, filename.c_str());
    writerT.WriteFile();
}
