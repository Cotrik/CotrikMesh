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
// #include "ParallelFor.h"

class Smoother {
    public:
        Smoother();
        Smoother(Mesh& mesh_);
        Smoother(const Smoother& r);
        ~Smoother();

        void SetMesh(Mesh& mesh_);

        void Smooth(Mesh& mesh_, std::vector<size_t>& V);
    private:
        Mesh& mesh = Mesh();
        SurfaceMapper sm;

        void CheckValidity();
        void GetVerticesToSmooth(int iter, Mesh& mesh_, std::vector<size_t>& V);
        void GetOptimizedPositions(int iter, Mesh& mesh_, std::vector<size_t>& V, std::vector<glm::dvec3>& centers);
        void SetCoords(Mesh& mesh_, size_t vid, glm::dvec3& c);
};

#endif