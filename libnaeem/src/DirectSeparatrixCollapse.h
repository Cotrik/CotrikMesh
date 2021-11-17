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
        DirectSeparatrixCollapse(Mesh& mesh_, MeshUtil& mu_, size_t cid_, std::vector<size_t> s1_, std::vector<size_t> s2_);
        ~DirectSeparatrixCollapse();

        void SetRanking();
        bool IsOperationValid();
        void PerformOperation();

    private:
        size_t cid;
        std::vector<size_t> s1, s2;
};

#endif