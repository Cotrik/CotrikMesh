/*
* VertexRotation.h
*
*  Created on: February 7, 2022
*      Author: https://github.com/naeem014
*/

#ifndef VERTEX_ROTATION_H_
#define VERTEX_ROTATION_H_

#include "SimplificationOperation.h"

class VertexRotation : public SimplificationOperation {
    public:
        VertexRotation();
        VertexRotation(Mesh& mesh_, MeshUtil& mu_, size_t vid_);
        ~VertexRotation();

        
        void SetRanking(glm::dvec3 d = glm::dvec3(0, 0, 0));
        bool IsOperationValid();
        void PerformOperation();
        glm::dvec3 GetLocation() {return glm::dvec3(0.0, 0.0, 0.0);}
        size_t GetCenterId() {return -1;}
        
    private:
        size_t vid;

};

#endif