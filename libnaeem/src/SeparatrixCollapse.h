/*
* SeparatrixCollapse.h
*
*  Created on: December 9, 2021
*      Author: https://github.com/naeem014
*/

#ifndef SEPARATRIX_COLLAPSE_H_
#define SEPARATRIX_COLLAPSE_H_

#include <glm/glm.hpp>

#include "SimplificationOperation.h"
#include "MeshUtil.h"

class SeparatrixCollapse : public SimplificationOperation {
    public:
        // Constructors and Destructor
        SeparatrixCollapse();
        SeparatrixCollapse(Mesh& mesh_, MeshUtil& mu_, Smoother& smoother_, std::vector<size_t> linkV, std::vector<size_t> linkE, bool half_ = false);
        ~SeparatrixCollapse();

        void SetRanking(glm::dvec3 d = glm::dvec3(0, 0, 0));
        bool IsOperationValid();
        void PerformOperation();
        glm::dvec3 GetLocation() {return glm::dvec3(0, 0, 0);}
        size_t GetCenterId() {return target.front();}

    private:
        std::vector<size_t> target;
        std::vector<std::vector<size_t>> collapse;

        void BuildSeparatrix(std::vector<size_t> linkV, std::vector<size_t> linkE);
        void Collapse(size_t targetId, size_t collapseId);
        std::vector<size_t> GetCollapseVids(size_t vid, size_t eid);
        void SetUpdateElements();
        double GetRanking();
        bool half = false;
};

#endif