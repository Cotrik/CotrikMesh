/*
* SemiGlobalSimplifier.h
*
*  Created on: October 25, 2021
*      Author: https://github.com/naeem014
*/

#ifndef SEMI_GLOBAL_SIMPLIFIER_H_
#define SEMI_GLOBAL_SIMPLIFIER_H_

#include <glm/glm.hpp>
#include <queue>

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
#include "QuadSplit.h"
#include "ThreeFivePair.h"
#include "MeshUtil.h"
#include "Smooth.h"
#include "PQueue.h"
#include "PQueue.cpp"

struct SingularityLink {
    std::vector<size_t> linkVids;
    std::vector<size_t> linkEids;
    int frontId;
    int backId;
    int id;
    int a = 0;
    int b = 0;
};

struct SingularityGroup {
    SingularityLink l1;
    SingularityLink l2;
    int rank;
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
        void SetBoundaryDirectSeparatrixOperations(bool looseCollapse);
        void SetDirectSeparatrixOperations(bool looseCollapse);
        void SetSeparatrixOperations();
        void SetBoundarySeparatrixOperations();
        void SetHalfSeparatrixOperations();
        void SetChordCollapseOperations();
        void SetEdgeRotationOperations();
        void SetVertexRotationOperations();
        void SetEdgeCollapseOperations();
        void SetVertexSplitOperations();
        void SetQuadSplitOperations();
        void GetSingularityPairs();
        void ResolveSingularityPairs();
        void CheckAndResolveThreeFivePair(size_t vid);
        void CheckAndResolveFiveThreePair(size_t vid);
        std::vector<SingularityLink> TraceSingularityLinks(Vertex& v, BaseComplexQuad& bc);
        void PerformGlobalOperations();
        void Smooth();

        bool doesCrossBoundary(std::vector<size_t> in, bool isVertex);

        MeshUtil& mu = MeshUtil();
        int iters = 0;
        
    private:
        Mesh& mesh = Mesh();
        Smoother& smoother = Smoother();

        void CheckValidity();
        bool cmpLink(SingularityLink left, SingularityLink right) {return left.a + left.b < right.a + right.b;}
        size_t GetFaceId(size_t vid, size_t exclude_vid);
        size_t GetDiagonalV(size_t vid, size_t fid);
        std::vector<std::shared_ptr<SimplificationOperation>> Ops;
        PQueue<std::shared_ptr<SimplificationOperation>> Op_Q;
};

#endif