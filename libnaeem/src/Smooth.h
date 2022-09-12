/*
* Smooth.h
*
*  Created on: November 15, 2021
*      Author: https://github.com/naeem014
*/

#ifndef SMOOTH_H_
#define SMOOTH_H_

#include "Mesh.h"
#include "QuadSurfaceMapper.h"
#include "MeshUtil.h"
// #include "ParallelFor.h"

class Smoother {
    public:
        Smoother();
        Smoother(Mesh& mesh_, MeshUtil& mu_);
        Smoother(const Smoother& r);
        ~Smoother();

        void SetMesh(Mesh& mesh_);

        void Smooth(Mesh& mesh_, std::vector<size_t>& V);
        void SetIters(int iters_) {iters = iters_;}
    private:
        Mesh& mesh = Mesh();
        MeshUtil& mu = MeshUtil();
        SurfaceMapper sm;

        int iters = 0;
        
        void CheckValidity();
        void GetVerticesToSmooth(int iter, Mesh& mesh_, std::vector<size_t>& V);
        void GetOptimizedPositions(int iter, Mesh& mesh_, std::vector<size_t>& V, std::vector<glm::dvec3>& centers);
        void SetOptimizedPositions(int iter, Mesh& mesh_, std::vector<size_t>& V);
        void SetPosition(Mesh& mesh_, Vertex& v);
        void SetPositionBoundary(Mesh& mesh_, Vertex& v);
        void SetCoords(Mesh& mesh_, size_t vid, glm::dvec3& c);
};

#endif