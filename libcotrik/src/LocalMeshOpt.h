/*
 * MeshOpt.h
 *
 *  Created on: Nov 25, 2016
 *      Author: cotrik
 */

#ifndef LOCAL_MESH_OPT_H_
#define LOCAL_MESH_OPT_H_

#include "Mesh.h"
#include "MeshOpt.h"

class LocalMeshOpt : public MeshOpt
{
public:
    LocalMeshOpt(const Mesh& mesh);
    virtual ~LocalMeshOpt();

protected:
    LocalMeshOpt();
    LocalMeshOpt(const LocalMeshOpt& meshOpt);

public:
    void Run(const size_t iters = 1, const size_t localIters = 20);
//    virtual bool Optimize();
    virtual void Untangle(const std::vector<size_t>& cellIds, const size_t localIters = 20);
    //virtual bool UntangleLocalMesh(Mesh& localMesh, const size_t localIters = 20);
    virtual bool UntangleLocalMesh(Mesh& localMesh, Mesh& bestMesh, const size_t localIters = 20);
    virtual void ModifyMeshFrom(const Mesh& localMesh, const std::vector<size_t>& badCellIdsT);
    void DivideIntoMultipleRegions(const std::vector<size_t>& badCellIds, std::vector<std::vector<size_t> >& regions, const int N = 50);
    void GetLocalRegions(const std::vector<size_t>& badCellIds, int extendLayers = 2);
    void SetUseSmallBlock(bool value = true);
    void SetBlockSize(size_t value = 50);
    void SetParallel(bool value = true);
protected:
//    void ComputeAvgEdgeLength();
//    void ComputeMeshTargetLength();
//    virtual bool Optimize();
//    virtual size_t OptimizeSurfaceVertices(std::vector<Trip>& coefficients, std::vector<float>& B, size_t& m_id);  // * alpha
//    virtual void OptimizeEdgeOrthogonality(std::vector<Trip>& coefficients, std::vector<float>& B, size_t& m_id);  // * gamma
//    virtual void OptimizeEdgeStraightness(std::vector<Trip>& coefficients, std::vector<float>& B, size_t& m_id);   // * gamma
//    virtual void OptimizeSingularity(std::vector<Trip>& coefficients, std::vector<float>& B, size_t& m_id);        // * gamma

    std::vector<size_t> m_localVertexIds;
    std::vector<size_t> m_localCellIds;
    bool useSmallBlock = false;
    bool parallel = false;
    size_t blockSize = 50;
};

#endif /* LOCAL_MESH_OPT_H_ */
