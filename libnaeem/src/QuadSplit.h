#ifndef QUAD_SPLIT_H
#define QUAD_SPLIT_H

#include <glm/glm.hpp>

#include "SimplificationOperation.h"
#include "MeshUtil.h"

class QuadSplit : public SimplificationOperation {
    public:
        QuadSplit();
        QuadSplit(Mesh& mesh_, MeshUtil& mu_, size_t vid_, std::vector<size_t> verticesToSplit_, std::vector<size_t> verticesToChange_);
        ~QuadSplit();

        void SetRanking(glm::dvec3 d = glm::dvec3(0, 0, 0));
        bool IsOperationValid();
        void PerformOperation();
        glm::dvec3 GetLocation() {return glm::dvec3(0, 0, 0);}
        size_t GetCenterId() {return -1;}

    private:
        size_t vid;
        std::vector<size_t> verticesToSplit;
        std::vector<size_t> verticesToChange;

};

#endif