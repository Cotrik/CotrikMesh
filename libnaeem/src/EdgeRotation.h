/*
* EdgeRotation.h
*
*  Created on: January 31, 2022
*      Author: https://github.com/naeem014
*/

#ifndef EDGE_ROTATION_H_
#define EDGE_ROTATION_H_

#include <glm/glm.hpp>

#include "SimplificationOperation.h"
#include "MeshUtil.h"

class EdgeRotation : public SimplificationOperation {
    public:
        EdgeRotation();
        EdgeRotation(Mesh& mesh_, MeshUtil& mu_, Smoother& smoother_, size_t eid_, bool clockwise_);
        ~EdgeRotation();
        
        void SetRanking(glm::dvec3 d = glm::dvec3(0, 0, 0));
        bool IsOperationValid();
        void PerformOperation();
        glm::dvec3 GetLocation() {return glm::dvec3(0.0, 0.0, 0.0);}
        size_t GetCenterId() {return -1;}

    private:
        std::vector<size_t> GetVertices(Edge& e, size_t v1, size_t v2);

        size_t eid;
        bool clockwise;
};

#endif