/*
* DirectSeparatrixCollapse.h
*
*  Created on: November 10, 2021
*      Author: https://github.com/naeem014
*/

#ifndef DIRECT_SEPARATRIX_COLLAPSE_H_
#define DIRECT_SEPARATRIX_COLLAPSE_H_

#include <glm/glm.hpp>

#include "SimplificationOperation.h"
#include "MeshUtil.h"

class DirectSeparatrixCollapse : public SimplificationOperation {
    public:
        // Constructors and Destructor
        DirectSeparatrixCollapse();
        DirectSeparatrixCollapse(Mesh& mesh_, MeshUtil& mu_, Smoother& smoother_, size_t cid_, std::vector<size_t> s1_, std::vector<size_t> s2_, bool looseCollapse_);
        ~DirectSeparatrixCollapse();

        void SetRanking(glm::dvec3 d = glm::dvec3(0, 0, 0));
        bool IsOperationValid();
        void PerformOperation();
        glm::dvec3 GetLocation();
        size_t GetCenterId() {return cid;}
    private:
        void SetUpdateElements(size_t vid);
        double GetDistance(glm::dvec3 a);

        size_t cid;
        std::vector<size_t> s1, s2;
        bool looseCollapse = false;
};

#endif