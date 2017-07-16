/*
 * MeshOpt.h
 *
 *  Created on: Nov 25, 2016
 *      Author: cotrik
 */

#ifndef LOCAL_MESH_ASJ_OPT_H_
#define LOCAL_MESH_ASJ_OPT_H_

#include "Mesh.h"
#include "MeshASJOpt.h"

class LocalMeshASJOpt : public MeshASJOpt
{
public:
    LocalMeshASJOpt(const Mesh& mesh);
    virtual ~LocalMeshASJOpt();

protected:
    LocalMeshASJOpt();
    LocalMeshASJOpt(const LocalMeshASJOpt& meshOpt);

public:
    void Run(const size_t iters = 1, const size_t localIters = 20);
//    virtual bool Optimize();
    virtual bool UntangleLocalMesh(Mesh& localMesh, const size_t localIters = 20);
    virtual bool UntangleLocalMesh(Mesh& localMesh, Mesh& bestMesh, const size_t localIters = 20);
    virtual void ModifyMeshFrom(const Mesh& localMesh, const std::vector<size_t>& badCellIdsT);
    void DivideIntoMultipleRegions(const std::vector<size_t>& badCellIds, std::vector<std::vector<size_t> >& regions);

//private:
//    size_t OptimizeSurfaceVertices(std::vector<Trip>& coefficients, std::vector<float>& B, size_t& m_id);  // * alpha
//    void OptimizeEdgeOrthogonality(std::vector<Trip>& coefficients, std::vector<float>& B, size_t& m_id);  // * gamma
//    void OptimizeEdgeStraightness(std::vector<Trip>& coefficients, std::vector<float>& B, size_t& m_id);   // * gamma
//    void OptimizeSingularity(std::vector<Trip>& coefficients, std::vector<float>& B, size_t& m_id);        // * gamma
};

#endif /* LOCAL_MESH_ASJ_OPT_H_ */
