/*
* Smooth.h
*
*  Created on: November 15, 2021
*      Author: https://github.com/naeem014
*/

#ifndef SMOOTH_H_
#define SMOOTH_H_

#include "Mesh.h"
#include "MeshUtil.h"
// #include "ParallelFor.h"

class Smoother {
    public:
        Smoother();
        Smoother(Mesh& mesh_, MeshUtil& mu_);
        Smoother(const Smoother& r);
        ~Smoother();

        void SetMembers(Mesh& mesh_, MeshUtil& mu_);

        void Smooth(std::vector<size_t> V);
        void SetIters(int iters_) {iters = iters_;}
    private:
        Mesh* mesh;
        MeshUtil* mu;

        int iters = 0;
        
        void CheckValidity();
        void GetVerticesToSmooth(int iter, std::vector<size_t>& V);
        void GetOptimizedPositions(int iter, std::vector<size_t>& V, std::vector<glm::dvec3>& centers);
        void SetOptimizedPositions(int iter, std::vector<size_t>& V);
        void SetPosition(Vertex& v);
        void SetPositionBoundary(Vertex& v);
        void SetCoords(size_t vid, glm::dvec3& c);
};

#endif