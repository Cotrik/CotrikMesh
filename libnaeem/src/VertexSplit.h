/*
* VertexSplit.h
*
*  Created on: March 7, 2022
*      Author: https://github.com/naeem014
*/

#ifndef VERTEX_SPLIT_H_
#define VERTEX_SPLIT_H_

#include "SimplificationOperation.h"

class VertexSplit : public SimplificationOperation {
    public:
        VertexSplit();
        VertexSplit(Mesh& mesh_, MeshUtil& mu_, Smoother& smoother_, size_t vid_, std::vector<size_t> edgesToSplit_ = {});
        ~VertexSplit();

        
        void SetRanking(glm::dvec3 d = glm::dvec3(0, 0, 0));
        bool IsOperationValid();
        void PerformOperation();
        glm::dvec3 GetLocation() {return glm::dvec3(0.0, 0.0, 0.0);}
        size_t GetCenterId() {return -1;}
        double CalculateRanking() {return 0.0;}
        
    private:
        size_t vid;
        std::vector<size_t> splitEdges;
        std::vector<size_t> SelectEdgesToSplit();

};

#endif