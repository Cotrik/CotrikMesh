/*
* ChordCollapse.h
*
*  Created on: January 19, 2022
*      Author: https://github.com/naeem014
*/

#ifndef CHORD_COLLAPSE_H_
#define CHORD_COLLAPSE_H_

#include <glm/glm.hpp>

#include "SimplificationOperation.h"
#include "MeshUtil.h"
#include "ChordExtractor.h"

class ChordCollapse : public SimplificationOperation {
    public:
        ChordCollapse();
        ChordCollapse(Mesh& mesh_, MeshUtil& mu_, Smoother& smoother_, ChordExtractor& ce_, size_t chordId_);
        ~ChordCollapse();

        void SetRanking(glm::dvec3 d = glm::dvec3(0, 0, 0));
        bool IsOperationValid();
        void PerformOperation();
        glm::dvec3 GetLocation() {return glm::dvec3(0.0, 0.0, 0.0);}
        size_t GetCenterId() {return -1;}
        double CalculateRanking();
        glm::dvec4 CalculateQEM(size_t v1_id, size_t v2_id);
        double CalculateAreaDistance(size_t v1_id, size_t v2_id);
        double CalculateValenceTerm(size_t v1_id, size_t v2_id);

    private:
        void CollapseEdge(Edge& e, Edge& edgeToRemove);
        void UpdateEdges(size_t fid, Edge& edgeToKeep, Edge& edgeToRemove);

        ChordExtractor& ce = ChordExtractor();
        Chord& chord = Chord();
        std::vector<size_t> edgesIds;
};

#endif