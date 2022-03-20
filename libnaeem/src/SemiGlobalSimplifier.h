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
#include "BaseComplexSheetQuad.h"
#include "SheetSimplifier.h"
#include "SingleSheetSimplifier.h"
#include "ChordExtractor.h"
#include "ChordCollapse.h"
// #include "SimplificationOperation.h"
#include "DiagonalCollapse.h"
#include "DirectSeparatrixCollapse.h"
#include "SeparatrixCollapse.h"
#include "EdgeRotation.h"
#include "VertexRotation.h"
#include "EdgeCollapse.h"
#include "VertexSplit.h"
#include "MeshUtil.h"
#include "Smooth.h"
#include "PQueue.h"
#include "PQueue.cpp"

struct SingularityLink {
    std::vector<size_t> linkVids;
    std::vector<size_t> linkEids;
};

class SemiGlobalSimplifier {
    public:
        // Constructors and Destructor
        SemiGlobalSimplifier();
        SemiGlobalSimplifier(Mesh& mesh_, MeshUtil& mu_, Smoother& smoother_);
        ~SemiGlobalSimplifier();

        // MeshUtil setters and getters
        void SetMesh(Mesh& mesh_);
        void SetIters(int iters_);

        // Simplification Operations
        void FixBoundary();
        void SetSimplificationOperations();
        void SetDiagonalCollapseOperations();
        void SetDirectSeparatrixOperations();
        void SetSeparatrixOperations();
        void SetHalfSeparatrixOperations();
        void SetChordCollapseOperations();
        void SetEdgeRotationOperations();
        void SetVertexRotationOperations();
        void SetEdgeCollapseOperations();
        void SetVertexSplitOperations();
        std::vector<SingularityLink> TraceSingularityLinks(Vertex& v, BaseComplexQuad& bc);
        void PerformGlobalOperations();
        void Smooth();

        MeshUtil& mu = MeshUtil();
        int iters = 0;
        
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