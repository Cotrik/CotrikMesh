/*
* EdgeCollapse.h
*
*  Created on: February 15, 2022
*      Author: https://github.com/naeem014
*/

#ifndef EDGE_COLLAPSE_H_
#define EGE_COLLAPSE_H_

#include "SimplificationOperation.h"

class EdgeCollapse : public SimplificationOperation {
    public:
        EdgeCollapse();
        EdgeCollapse(Mesh& mesh_, MeshUtil& mu_, Smoother& smoother_, std::shared_ptr<SimplificationOperation> vr_, std::shared_ptr<SimplificationOperation> dc_);
        ~EdgeCollapse();

        void SetRanking(glm::dvec3 d = glm::dvec3(0, 0, 0));
        bool IsOperationValid();
        void PerformOperation();
        glm::dvec3 GetLocation() {return glm::dvec3(0.0, 0.0, 0.0);}
        size_t GetCenterId() {return -1;}

    private:
        std::shared_ptr<SimplificationOperation> vr;
        std::shared_ptr<SimplificationOperation> dc;
};

#endif