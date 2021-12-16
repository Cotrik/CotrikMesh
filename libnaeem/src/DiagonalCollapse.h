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
#include "MeshUtil.h"

class DiagonalCollapse : public SimplificationOperation {
    public:
        // Constructors and Destructor
        DiagonalCollapse();
        DiagonalCollapse(Mesh& mesh_, MeshUtil& mu_, size_t f, size_t d_idx1_, size_t d_idx2_);
        ~DiagonalCollapse();

        void SetRanking(glm::dvec3 d = glm::dvec3(0, 0, 0));
        bool IsOperationValid();
        void PerformOperation();
        glm::dvec3 GetLocation() {return glm::dvec3(0, 0, 0);}
        size_t GetCenterId() {return fId;}

    private:
        void UpdateNeighborInfo(Vertex& target, Vertex& source);
        size_t fId, d_idx1, d_idx2;
};

#endif