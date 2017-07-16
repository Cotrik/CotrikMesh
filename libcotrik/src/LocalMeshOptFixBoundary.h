/*
 * LocalMeshOptBoundary.h
 *
 *  Created on: Nov 25, 2016
 *      Author: cotrik
 */

#ifndef LOCAL_MESH_OPT_FIX_BOUNDARY_H_
#define LOCAL_MESH_OPT_FIX_BOUNDARY_H_

#include "MeshOptFixBoundary.h"

class LocalMeshOptFixBoundary : public MeshOptFixBoundary
{
public:
    LocalMeshOptFixBoundary(const Mesh& mesh);
    virtual ~LocalMeshOptFixBoundary();

protected:
    LocalMeshOptFixBoundary();
    LocalMeshOptFixBoundary(const LocalMeshOptFixBoundary& meshOpt);

public:
    void Run(const size_t iters = 1, const size_t localIters = 20);
    virtual bool UntangleLocalMesh(Mesh& localMesh, Mesh& bestMesh, const size_t localIters = 20);
    virtual void ModifyMeshFrom(const Mesh& localMesh, const std::vector<size_t>& badCellIdsT);
    void DivideIntoMultipleRegions(const std::vector<size_t>& badCellIds, std::vector<std::vector<size_t> >& regions);
};

#endif /* LOCAL_MESH_OPT_FIX_BOUNDARY_H_ */
