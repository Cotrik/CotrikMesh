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
#include "BaseComplexQuad.h"
// #include "SimplificationOperation.h"
#include "DiagonalCollapse.h"
#include "DirectSeparatrixCollapse.h"
#include "MeshUtil.h"
#include "Smooth.h"
#include "PQueue.h"
#include "PQueue.cpp"

class SemiGlobalSimplifier {
    public:
        // Constructors and Destructor
        SemiGlobalSimplifier();
        SemiGlobalSimplifier(Mesh& mesh_, MeshUtil& mu_, Smoother& smoother_);
        ~SemiGlobalSimplifier();

        // MeshUtil setters and getters
        void SetMesh(Mesh& mesh_);

        // Simplification Operations
        void SetSimplificationOperations();
        void SetDiagonalCollapseOperations();
        void SetDirectSeparatrixOperations();
        void SetSeparatrixOperations();
        void TraceSingularityLinks(Vertex& v, BaseComplexQuad& bc);
        void PerformGlobalOperations();

        MeshUtil& mu = MeshUtil();
        
    private:
        Mesh& mesh = Mesh();
        Smoother& smoother = Smoother();

        void CheckValidity();
        size_t GetFaceId(size_t vid, size_t exclude_vid);
        size_t GetDiagonalV(size_t vid, size_t fid);
        std::vector<std::shared_ptr<SimplificationOperation>> Ops;
        PQueue<std::shared_ptr<SimplificationOperation>> Op_Q;
};

#endif