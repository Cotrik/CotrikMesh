/*
* SemiGlobalSimplifier.h
*
*  Created on: October 25, 2021
*      Author: https://github.com/naeem014
*/

#ifndef SEMI_GLOBAL_SIMPLIFIER_H_
#define SEMI_GLOBAL_SIMPLIFIER_H_

#include <glm/glm.hpp>

#include "Mesh.h"
// #include "SimplificationOperation.h"
#include "DiagonalCollapse.h"
#include "DirectSeparatrixCollapse.h"
#include "MeshUtil.h"
#include "Smooth.h"

class SemiGlobalSimplifier {
    public:
        // Constructors and Destructor
        SemiGlobalSimplifier();
        SemiGlobalSimplifier(Mesh& mesh_);
        ~SemiGlobalSimplifier();

        // MeshUtil setters and getters
        void SetMesh(Mesh& mesh_);

        // Simplification Operations
        void SetSimplificationOperations();
        void SetDiagonalCollapseOperations();
        void SetDirectSeparatrixOperations();

        MeshUtil mu;
        
    private:
        Mesh& mesh = Mesh();
        Smoother smoother;

        void CheckValidity();
        std::vector<std::unique_ptr<SimplificationOperation>> Ops;
};

#endif