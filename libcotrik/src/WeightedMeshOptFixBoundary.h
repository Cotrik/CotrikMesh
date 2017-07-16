/*
 * WeightedMeshOptFixBoundary.h
 *
 *  Created on: Nov 25, 2016
 *      Author: cotrik
 */

#ifndef WEIGHTED_MESH_OPT_FIX_BOUNDARY_H_
#define WEIGHTED_MESH_OPT_FIX_BOUNDARY_H_

#include "MeshOptFixBoundary.h"

class WeightedMeshOptFixBoundary: public MeshOptFixBoundary
{
public:
    WeightedMeshOptFixBoundary(const Mesh& mesh);
    virtual ~WeightedMeshOptFixBoundary();

private:
    WeightedMeshOptFixBoundary();
    WeightedMeshOptFixBoundary(const WeightedMeshOptFixBoundary& meshOpt);

public:
    virtual size_t Run(const size_t iters = 1);
    virtual bool Optimize();
    virtual void OptimizeEdgeOrthogonality(std::vector<Trip>& coefficients, std::vector<float>& B, size_t& m_id);  // * gamma
    virtual void OptimizeEdgeStraightness(std::vector<Trip>& coefficients, std::vector<float>& B, size_t& m_id);   // * gamma
    virtual void OptimizeSingularity(std::vector<Trip>& coefficients, std::vector<float>& B, size_t& m_id);        // * gamma
    virtual void ComputeMeshTargetLength();

private:

};

#endif /* WEIGHTED_MESH_OPT_FIX_BOUNDARY_H_ */
