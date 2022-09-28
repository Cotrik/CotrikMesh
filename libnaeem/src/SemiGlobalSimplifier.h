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
#include <memory>

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
    double rank = 1.0;
    int rot = 0;
};

struct LinkComparator {
    public:

    bool operator()(SingularityLink& l, SingularityLink& r) {
        if (l.rank == r.rank) {
            return (l.a+l.b) > (r.a+r.b);
        }
        return l.rank > r.rank;
        // return (l.a+l.b) >= (r.a+r.b) && l.rank >= r.rank;
    }
};

struct SingularityGroup {
    SingularityLink l1;
    SingularityLink l2;
    double rank = 0.0;
};

struct GroupComparator {
    public:

    bool operator()(SingularityGroup& l, SingularityGroup& r) {
        return l.rank > r.rank;
    }
};

class SemiGlobalSimplifier {
    public:
        // Constructors and Destructor
        SemiGlobalSimplifier();
        SemiGlobalSimplifier(Mesh& mesh_, MeshUtil& mu_, Smoother& smoother_);
        ~SemiGlobalSimplifier();

        // MeshUtil setters and getters
        void SetMembers(Mesh& mesh_, MeshUtil& mu_, Smoother& Smoother_);
        void SetIters(int iters_);

        // Simplification Operations
        bool FixBoundary();
        bool FixValences();
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
        void SetSingularityLinks(std::vector<SingularityLink>& SingularityLinks, std::vector<std::vector<size_t>>& SingularityMap, BaseComplexQuad& bc);
        void SelectSingularityGroups(std::vector<SingularityGroup>& Groups, std::vector<SingularityLink>& SingularityLinks, std::vector<std::vector<size_t>>& SingularityMap);
        void ResolveSingularityGroups(std::vector<SingularityGroup>& Groups, BaseComplexQuad& bc);
        void ResolveSingularityPairs();
        void CheckAndResolveThreeFivePair(size_t vid);
        void CheckAndResolveFiveThreePair(size_t vid);
        std::vector<SingularityLink> TraceSingularityLinks(Vertex& v, BaseComplexQuad& bc);
        void PerformGlobalOperations();
        void Smooth();

        bool ResolveHighValences();
        void AlignSingularities();
        void AlignSingularities(size_t vid, std::queue<size_t>& Singularities, std::vector<bool>& isAvailable, BaseComplexQuad& bc, bool checkValence = true);
        void AlignAndResolveSingularities(bool checkValence = true);
        bool IsPair(size_t vid);
        std::vector<size_t> GetPair(size_t vid);
        int MovePair(std::vector<size_t> threeFiveIds, std::vector<size_t>& secondaryPath, bool checkValence = false);
        bool ResolveIsolatedSingularities(BaseComplexQuad& bc);
        void GenerateSingularityPair(SingularityLink& l1, SingularityLink& l2);
        void ResolveSingularities();
        void GetSingularityGroups(std::vector<size_t> Singularities, BaseComplexQuad& bc);
        int PullSingularity(SingularityLink& l1, SingularityLink& l2);
        std::vector<int> GetTraverseInfo(SingularityLink& l1, SingularityLink& l2);
        std::vector<size_t> TraversePath(size_t prev, size_t current, std::vector<int> Rots);
        std::vector<size_t> GetSecondaryPath(int offset, std::vector<size_t>& mainPath, BaseComplexQuad& bc);
        std::vector<SingularityLink> GetLinks(size_t sid, BaseComplexQuad& bc, bool checkValence = true);
        std::vector<SingularityLink> SelectLinks(std::vector<SingularityLink> links, int valence, bool checkValence = true);
        std::vector<SingularityLink> GetCrossLinks(SingularityLink& l, BaseComplexQuad& bc, bool checkValence = true);
        bool ValidateLink(SingularityLink& l);
        bool ValidatePath(std::vector<size_t> p);
        int ResolveSingularity(size_t sid, BaseComplexQuad& bc);
        int MoveSingularity(SingularityGroup& sg);
        int MoveSingularity(int offset, std::vector<size_t>& mainPath, std::vector<size_t>& secondaryPath);
        bool IsExclusive(SingularityLink& l1, SingularityLink& l2);
        bool IsExclusive(size_t vid, std::vector<size_t> a, std::vector<size_t> b);

        bool doesCrossBoundary(std::vector<size_t> in, bool isVertex);

        void PrototypeBoundary(bool checkValence = true);
        void PrototypeA(size_t vid, BaseComplexQuad& bc, bool checkValence = true);
        void PrototypeB();
        void PrototypeC();
        void PrototypeD();
        void PrototypeE();
        int PrototypeGetRotations(size_t vid, size_t start, size_t end);
        void PrototypeSaveMesh(SingularityLink& l1, SingularityLink& l2, std::string in);
        void PrototypeResolveGroup(SingularityLink& l1, SingularityLink& l2);
        int PrototypeCancelThreeFivePair(SingularityLink& l1, SingularityLink& l2);
        int PrototypeCancelSingularity(size_t vid, BaseComplexQuad& bc);
        bool PrototypeCheckBoundarySingularity(size_t vid);
        int PrototypeCancelSingularityPair(SingularityLink& l, BaseComplexQuad& bc);
        SingularityLink PrototypeGetLink(size_t vid, BaseComplexQuad& bc, size_t vertexToSkip = 0, std::vector<size_t> edgesToCheck = {}, bool checkValence = true, bool boundary = true);

        MeshUtil* mu;
        int iters = 0;
        
    private:
        Mesh* mesh;
        Smoother* smoother;

        void CheckValidity();
        size_t GetFaceId(size_t vid, size_t exclude_vid);
        size_t GetDiagonalV(size_t vid, size_t fid);
        std::vector<size_t> GetThreeFivePairIds(size_t vid, size_t mainId, size_t secondaryId);
        void MoveSingularities(size_t& toMoveId, size_t& sourceId, size_t& secondaryId, size_t& sourceDir, size_t& secondaryDir, std::vector<size_t>& secondaryPath);
        void SetSecondaryPath(size_t& secondaryId, size_t& toMoveId, size_t& sourceId, std::vector<size_t>& mainPath, std::vector<size_t>& secondaryPath, BaseComplexQuad& bc);
        std::vector<std::shared_ptr<SimplificationOperation>> Ops;
        PQueue<std::shared_ptr<SimplificationOperation>> Op_Q;
};

#endif