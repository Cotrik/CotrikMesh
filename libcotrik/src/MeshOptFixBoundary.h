/*
 * MeshOptFixBoundary.h
 *
 *  Created on: Nov 25, 2016
 *      Author: cotrik
 */

#ifndef MESH_OPT_FIX_BOUNDARY_H_
#define MESH_OPT_FIX_BOUNDARY_H_

#include "MeshOpt.h"

class MeshOptFixBoundary: public MeshOpt
{
public:
    MeshOptFixBoundary(const Mesh& mesh);
    virtual ~MeshOptFixBoundary();

protected:
    MeshOptFixBoundary();
    MeshOptFixBoundary(const MeshOptFixBoundary& meshOpt);

public:
    virtual size_t Run(const size_t iters = 1);
    virtual bool Optimize();
    virtual void OptimizeEdgeOrthogonality(std::vector<Trip>& coefficients, std::vector<float>& B, size_t& m_id);  // * gamma
    virtual void OptimizeEdgeStraightness(std::vector<Trip>& coefficients, std::vector<float>& B, size_t& m_id);   // * gamma
    virtual void OptimizeSingularity(std::vector<Trip>& coefficients, std::vector<float>& B, size_t& m_id);        // * gamma

private:

};

#endif /* MESH_OPT_FIX_BOUNDARY_H_ */
