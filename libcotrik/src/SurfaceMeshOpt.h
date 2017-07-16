/*
 * SurfaceMeshOpt.h
 *
 *  Created on: Apr 8, 2017
 *      Author: cotrik
 */

#ifndef SRC_SURFACEMESHOPT_H_
#define SRC_SURFACEMESHOPT_H_

#include "MeshOpt.h"
class SurfaceMeshOpt : public MeshOpt
{
public:
    SurfaceMeshOpt(Mesh& mesh);
    virtual ~SurfaceMeshOpt();
private:
    SurfaceMeshOpt();
    SurfaceMeshOpt(const SurfaceMeshOpt& meshOpt);
public:
    virtual size_t Run(const size_t iters = 1);
    virtual bool Optimize();
protected:
    virtual size_t OptimizeSurfaceVertices(std::vector<Trip>& A_Entries, std::vector<float>& b, size_t& row);
    virtual void OptimizeEdgeOrthogonality(std::vector<Trip>& A_Entries, std::vector<float>& b, size_t& row);
    virtual void OptimizeEdgeStraightness(std::vector<Trip>& A_Entries, std::vector<float>& b, size_t& row);
    virtual void OptimizeSingularity(std::vector<Trip>& A_Entries, std::vector<float>& b, size_t& row);
private:
    //std::vector<size_t> surfaceVids;
};

#endif /* SRC_SURFACEMESHOPT_H_ */
