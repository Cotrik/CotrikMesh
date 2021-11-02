/*
* DiagonalCollapse.h
*
*  Created on: October 25, 2021
*      Author: https://github.com/naeem014
*/

#ifndef DIAGONAL_COLLAPSE_H_
#define DIAGONAL_COLLAPSE_H_

#include <glm/glm.hpp>

#include "SimplificationOperation.h"

class DiagonalCollapse : public SimplificationOperation {
    public:
        // Constructors and Destructor
        DiagonalCollapse();
        DiagonalCollapse(Mesh& mesh_, size_t f, size_t d_idx1_, size_t d_idx2_);
        ~DiagonalCollapse();

        void SetRanking(MeshUtil& mu);
        bool IsOperationValid();
        void PerformOperation();

    private:
        void UpdateNeighborInfo(Vertex& target, Vertex& source);
        void UpdateContents(std::vector<size_t>& a, std::vector<size_t>& b);
        std::vector<size_t> GetDifference(std::vector<size_t>& a, std::vector<size_t>& b);
        size_t fId, d_idx1, d_idx2;
};

#endif