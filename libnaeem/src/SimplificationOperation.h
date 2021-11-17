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
        SimplificationOperation(Mesh& mesh_, MeshUtil& mu_);
        ~SimplificationOperation();

        // setters and getters
        void SetMesh(Mesh& mesh_);

        virtual void SetRanking() = 0;
        virtual bool IsOperationValid() = 0;
        virtual void PerformOperation() = 0;

        double ranking = -1.0;
        std::vector<size_t> smoothV;

    protected:
        void CheckValidity();
        void AddContents(std::vector<size_t>& a, std::vector<size_t>& b);
        void UpdateContents(std::vector<size_t>& a, std::vector<size_t>& b);
        std::vector<size_t> GetDifference(std::vector<size_t>& a, std::vector<size_t>& b);
        
        Mesh& mesh = Mesh();
        MeshUtil& mu = MeshUtil();
};

#endif