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

class Smoother {
    public:
        Smoother();
        Smoother(Mesh& mesh_);
        ~Smoother();

        void SetMesh(Mesh& mesh_);

        void Smooth(std::vector<size_t>& V);
    private:
        Mesh& mesh = Mesh();
        SurfaceMapper sm;

        void CheckValidity();
};

#endif