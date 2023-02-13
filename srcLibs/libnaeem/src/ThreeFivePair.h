#ifndef THREE_FIVE_PAIR_H
#define THREE_FIVE_PAIR_H

#include <glm/glm.hpp>
#include "SingularityPair.h"

// #include "SimplificationOperation.h"
// #include "DiagonalCollapse.h"
// #include "EdgeRotation.h"
// #include "QuadSplit.h"
// #include "VertexSplit.h"

// #include "Mesh.h"
// #include "MeshUtil.h"
// #include "Smooth.h"

class ThreeFivePair : public SingularityPair {
    public:
        ThreeFivePair();
        ThreeFivePair(Mesh& mesh_, MeshUtil& mu_, Smoother& smoother_, size_t threeId_, size_t fiveId_);
        ~ThreeFivePair();

        void SetMembers(Mesh& mesh_, MeshUtil& mu_, Smoother& smoother_, size_t threeId_, size_t fiveId_);

        // void MoveUpperLeft();
        // void MoveUpperRight();
        // void MoveLeft();
        // void MoveRight();
        // void MoveLowerLeft();
        // void MoveLowerRight();
        void SplitSixSingularity();
        int Move(size_t dest, double& delta, bool skipCheck = true);
        void SetResolvedSingularity(size_t dest, int ref);
        int GetResolvedSingularity();
        std::vector<size_t> GetPairIds();
        
    private:
        // void CheckValidity();
        // bool IsOperationValid();
        void SetEdgeId();
        // void AddContents(std::vector<size_t>& a, std::vector<size_t> b);
        // void UpdateContents(std::vector<size_t>& a, std::vector<size_t> b);
        // std::vector<size_t> GetDifference(std::vector<size_t> a, std::vector<size_t> b);
        // std::vector<size_t> GetUnion(std::vector<size_t> a, std::vector<size_t> b);
        // std::vector<size_t> GetIntersection(std::vector<size_t> a, std::vector<size_t> b);
    
        // Mesh* mesh;
        // MeshUtil* mu;
        // size_t threeId, fiveId, 
        size_t edgeId;
        std::vector<size_t> pairEdges;
        int resolvedSingularity = -1;
        std::map<std::string, int> refMap = {{"Upper3", 1}, {"Upper5", 2}, {"Middle3", 3}, {"Middle5", 4}, {"Lower3", 5}, {"Lower5", 6}};
        // Smoother* smoother;
};

#endif