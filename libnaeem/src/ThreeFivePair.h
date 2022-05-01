#ifndef THREE_FIVE_PAIR_H
#define THREE_FIVE_PAIR_H

#include <glm/glm.hpp>

#include "SimplificationOperation.h"
#include "DiagonalCollapse.h"
#include "EdgeRotation.h"
#include "QuadSplit.h"

#include "Mesh.h"
#include "MeshUtil.h"

class ThreeFivePair {
    public:
        ThreeFivePair();
        ThreeFivePair(Mesh& mesh_, MeshUtil& mu_, size_t threeId_, size_t fiveId_);
        ~ThreeFivePair();

        void MoveUpperLeft();
        void MoveUpperRight();
        void MoveLeft();
        void MoveRight();
        void MoveLowerLeft();
        void MoveLowerRight();
        void SplitSixSingularity();
        void Move(size_t dest);
        
    private:
        void CheckValidity();
        bool IsOperationValid();
        void AddContents(std::vector<size_t>& a, std::vector<size_t>& b);
        void UpdateContents(std::vector<size_t>& a, std::vector<size_t>& b);
        std::vector<size_t> GetDifference(std::vector<size_t>& a, std::vector<size_t>& b);
        std::vector<size_t> GetUnion(std::vector<size_t>& a, std::vector<size_t>& b);
        std::vector<size_t> GetIntersection(std::vector<size_t>& a, std::vector<size_t>& b);
    
        Mesh& mesh = Mesh();
        MeshUtil& mu = MeshUtil();
        size_t threeId, fiveId, edgeId;
        std::vector<size_t> pairEdges;
};

#endif