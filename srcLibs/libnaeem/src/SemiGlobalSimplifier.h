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
#include <mutex>
#include <map>

#include "Mesh.h"
#include "BaseComplexQuad.h"
#include "ChordExtractor.h"
#include "ChordCollapse.h"
// #include "SimplificationOperation.h"
#include "DiagonalCollapse.h"
#include "DiagonalThreeFivePair.h"
#include "DirectSeparatrixCollapse.h"
#include "SeparatrixCollapse.h"
#include "EdgeRotation.h"
#include "VertexRotation.h"
#include "EdgeCollapse.h"
#include "VertexSplit.h"
#include "QuadSplit.h"
#include "ThreeFivePair.h"
#include "MeshUtil.h"
#include "SingularityPair.h"
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
    int rots = 0;
    int volt = 0;
    bool diagonal = false;
    double* delta;
};

struct Separatrix {
    size_t threeId;
    size_t fiveId;
    size_t frontId;
    size_t backId;
    int b1 = 0;
    int b2 = 0;
    std::vector<size_t> vidsA;
    std::vector<size_t> vidsB;
    bool empty = true;
    int minSize = -1;
    int rots = 0;
    std::map<std::string, int> dirMap;
    /*std::vector<int> frontdirMap;
    std::vector<int> backdirMap;
    std::vector<int> score;
    std::vector<int> frontScore;
    std::vector<int> backScore;*/
};

struct SeparatrixComparator {
    public:

    bool operator()(Separatrix& l, Separatrix& r) {
        return l.b1+l.b2 > r.b1+r.b2;
    }
};

struct LinkComparator {
    public:

    bool operator()(SingularityLink& l, SingularityLink& r) {
        // if (l.rank == r.rank) {
        //     return (l.a+l.b) > (r.a+r.b);
        // }
        // return fabs(*l.delta + l.rank) > fabs(*l.delta + r.rank);
        return fabs(l.rank) > fabs(r.rank);
        // return (l.a+l.b) >= (r.a+r.b) && l.rank >= r.rank;
    }
};

struct SingularityGroup {
    SingularityLink l1;
    SingularityLink l2;
    double rank = 0.0;
    int l1l2Rots = 0;
    int l2l1Rots = 0;
    int l1l2Volt = 0;
    int l2l1Volt = 0;
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
        bool RemoveDoublets();
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
        void TraceSingularityLinks(size_t vid, std::vector<SingularityLink>& links, bool traceDiagonals = false);
        void TraceLink(const Vertex& v, const Edge& edge, std::vector<SingularityLink>& links, std::vector<size_t> vids = {}, std::vector<size_t> eids = {});
        void TraceLink(const Vertex& v, const Face& face, std::vector<SingularityLink>& links, std::vector<size_t> vids = {});
        void TraceAlongEdge(const Vertex& start_vertex, const Edge& start_edge, std::vector<size_t>& vids_link, bool checkBoundary = true);
        void TraceAlongEdge(const Vertex& start_vertex, const Edge& start_edge, std::vector<size_t>& vids_link, std::vector<size_t>& eids_link);
        void TraceAlongDiagonal(const Vertex& start_vertex, const Face& start_face, std::vector<size_t>& vids_link);
        void GetDiagonalPath(SingularityLink& l);
        int GetEdgeRots(size_t eid1, size_t eid2);
        int GetFaceRots(size_t fid1, size_t fid2, size_t vid);
        void SelectLinks(std::vector<SingularityLink>& links);
        void SelectDirectPairLink(std::vector<size_t> threeFiveIds, std::vector<SingularityLink>& links, std::vector<size_t> verticesToAvoid = {});
        void SelectDiagonalPairLink(std::vector<size_t> threeFiveIds, std::vector<SingularityLink>& links, std::vector<size_t> verticesToAvoid = {});
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
        bool ValidatePath(std::vector<size_t> p, bool checkValences = true);
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
        void PrototypeF(int idxOffset = 2);
        void PrototypeG(int vid, BaseComplexQuad& bc);
        bool PrototypeH(int idxOffset = 2);
        bool PrototypeI();
        void PrototypeJ();
        void PrototypeMoveSingularity(SingularityLink& l);

        SingularityLink PrototypeGetLink(size_t vid);
        SingularityLink PrototypePairLink(std::vector<size_t> threeFiveIds, std::vector<size_t>& verticesToAvoid, bool isDiagonal = false);
        bool PrototypePairIds(std::vector<size_t>& threeFiveIds, std::vector<size_t>& vids, size_t toMoveIdx, size_t destIdx);
        void PrototypeResolvePairIds(std::vector<size_t>& threeFiveIds, std::vector<size_t>& vids, size_t toMoveIdx, size_t destIdx);
        void PrototypeMovePair(std::vector<size_t>& threeFiveIds, SingularityLink& l, bool isDiagonal = false);
        bool PrototypeIsLinkValid(SingularityLink& l, std::vector<size_t> verticesToAvoid = {});

        int PrototypeGetRotations(size_t vid, size_t start, size_t end);
        void PrototypeSaveMesh(SingularityLink& l1, SingularityLink& l2, std::string in);
        void PrototypeSaveMesh(std::vector<SingularityLink>& links, std::string in);
        void PrototypeResolveGroup(SingularityLink& l1, SingularityLink& l2);
        int PrototypeCancelThreeFivePair(SingularityLink& l1, SingularityLink& l2);
        int PrototypeCancelSingularity(size_t vid, BaseComplexQuad& bc);
        bool PrototypeCheckBoundarySingularity(size_t vid);
        int PrototypeCancelSingularityPair(SingularityLink& l, BaseComplexQuad& bc);
        SingularityLink PrototypeGetLink(size_t vid, BaseComplexQuad& bc, size_t vertexToSkip = 0, std::vector<size_t> edgesToCheck = {}, bool checkValence = true, bool boundary = true);
        bool PrototypeIsLinkValid(SingularityGroup& s);
        int PrototypeGetElementPrediction(SingularityGroup& s);


        void PrototypeK();
        void PrototypeTraceSeparatrix(size_t vid, size_t eid, Separatrix& s, bool checkValence = false);
        void PrototypeTransportSeparatrix(Separatrix& s, bool includeIters, std::vector<std::vector<size_t>>& vec);
        void PrototypeTransportDiagonalPair(size_t threeId, size_t fiveId, bool includeIters, std::vector<size_t> verticesToAvoid = {});
        void PrototypeTraceDiagSeparatrix(size_t vid, size_t eid, Separatrix& s, bool checkBoundary, std::vector<size_t> verticesToAvoid);
        void PrototypeTransportDirectPair(size_t threeId, size_t fiveId, bool includeIters, std::vector<size_t> verticesToAvoid = {});
        void PrototypeTraceDirectSeparatrix(size_t vid, size_t eid, Separatrix& s, std::vector<size_t> verticesToAvoid, int rotOffset, bool isRotEdge);
        void PrototypeTransportDirectLink(Separatrix& primary, std::vector<std::vector<size_t>>& vec);
        void PrototypeTraceDirectLinkSeparatrix(size_t vid, size_t eid, Separatrix& s, size_t vertexToAvoid, int rotOffset, bool isRotEdge);
        void PrototypeMoveLinkPair(Separatrix& p, Separatrix& s, size_t singularity, size_t moveDir, size_t sourceDir);
        void PrototypeTransportLink(Separatrix& primary);
        void PrototypeMoveLinkPair(size_t singularityId, size_t moveDir, std::vector<size_t> verticesToAvoid);
        size_t GetEdgeId(size_t vid, size_t vid2);
        size_t GetNextCCedge(size_t vid, size_t startEid);
        size_t GetCCedgeAt(size_t vid, size_t eid, int counter);
        void TraceAlongEdge(const Vertex& start_vertex, const Edge& start_edge, std::vector<size_t>& vids_link, int b1, int b2, int rotOffset, bool checkBoundary = true);
        void PrototypeSaveSeparatrices(std::vector<std::vector<size_t>> vec, std::string in = "test");

        void SaveEdgeMesh();


        void PrototypeExecute();
        std::vector<bool> PrototypeSetSingularities();
        void PrototypeExtractSingularityLinks(size_t vid, std::vector<Separatrix>& separatrices);
        void PrototypeTraceLinks(size_t vid, std::vector<Separatrix>& seps);
        bool PrototypeGetPairIds(size_t singularityId, size_t moveDir, std::vector<size_t>& threeFiveIds, int secVid = -1);
        void PrototypeSetDirMap(Separatrix& s);
        void PrototypeSet33DirMap(Separatrix& s);
        void PrototypeSet55DirMap(Separatrix& s);
        void PrototypeSet35DirMap(Separatrix& s);

        size_t GettCCFaceAt(size_t vid, size_t eid, int counter);


        bool CheckMeshValidity();

        MeshUtil* mu;
        int iters = 0;
        
    private:
        Mesh* mesh;
        Smoother* smoother;

        void CheckValidity();
        size_t GetFaceID(size_t vid, size_t exclude_vid);
        size_t GetFaceId(size_t vid1, size_t vid2);
        size_t GetDiagonalV(size_t vid, size_t fid);
        size_t GetFaceV(size_t vid, size_t fid, int offset);
        bool Contains(std::vector<size_t> v, size_t val);
        bool Contains(std::vector<size_t> v, std::vector<size_t> v2);
        std::vector<size_t> GetThreeFivePairIds(size_t vid, size_t mainId, size_t secondaryId);
        void MoveSingularities(size_t& toMoveId, size_t& sourceId, size_t& secondaryId, size_t& sourceDir, size_t& secondaryDir, std::vector<size_t> secondaryPath);
        void SetSecondaryPath(size_t& secondaryId, size_t& toMoveId, size_t& sourceId, std::vector<size_t>& mainPath, std::vector<size_t>& secondaryPath, BaseComplexQuad& bc);
        std::vector<std::shared_ptr<SimplificationOperation>> Ops;
        PQueue<std::shared_ptr<SimplificationOperation>> Op_Q;
        std::recursive_mutex mtx;
        double delta = 0.0;
};

#endif