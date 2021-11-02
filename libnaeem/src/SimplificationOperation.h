/*
* SimplificationOperation.h
*
*  Created on: October 25, 2021
*      Author: https://github.com/naeem014
*/

#ifndef SIMPLIFICATION_OPERATION_H_
#define SIMPLIFICATION_OPERATION_H_

#include <glm/glm.hpp>

#include "Mesh.h"
#include "MeshUtil.h"

class SimplificationOperation {
    public:
        // Constructors and Destructor
        SimplificationOperation();
        SimplificationOperation(Mesh& mesh_);
        ~SimplificationOperation();

        // MeshUtil setters and getters
        void SetMesh(Mesh& mesh_);

        virtual void SetRanking(MeshUtil& mu) = 0;
        virtual bool IsOperationValid() = 0;
        virtual void PerformOperation() = 0;

        double ranking = -1.0;

    protected:
        void CheckValidity();
        Mesh& mesh = Mesh();
};

#endif