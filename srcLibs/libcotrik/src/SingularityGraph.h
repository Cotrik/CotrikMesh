/*
 * SingularityGraph.h
 *
 *  Created on: Jul 25, 2017
 *      Author: cotrik
 */

#ifndef LIBCOTRIK_SRC_SINGULARITYGRAPH_H_
#define LIBCOTRIK_SRC_SINGULARITYGRAPH_H_

#include "BaseComplex.h"

struct SingularVertex{
    size_t id;
    size_t id_mesh;
    std::vector<size_t> neighborSingularVertexIds;
    std::vector<size_t> neighborSingularEdgeIds;
};

struct SingularEdge{
    size_t id = MAXID;
    size_t valence = MAXID;
    std::vector<size_t> vids_link;
    std::vector<size_t> eids_link;
    std::vector<size_t> singularVertexIds;
    std::vector<size_t> directlyLinkedSingularEdgeIds;
    std::vector<size_t> linkedByOneComponentEdgeSingularEdgeIds;
};

class SingularityGraph
{
public:
    SingularityGraph(const BaseComplex& baseComplex);
    SingularityGraph(const SingularityGraph& rhs);
    virtual ~SingularityGraph();
private:
    SingularityGraph();
public:
    const BaseComplex& GetBaseComplex() const;
    void Build();
    void WriteMatrixVTK(const char* filename, const std::vector<std::vector<size_t>>& m);
    void WriteMatrixMat(const char* filename, const std::vector<std::vector<size_t>>& m);
private:
    void BuildV();
    void BuildE();
    void BuildV_V();
    void BuildV_E();
    void BuildE_directlyLinkedSingularEdgeIds();
    void BuildE_linkedByOneComponentEdgeSingularEdgeIds();
    //void BuildE_linkedByOneComponentEdgeSingularEdgeIds();
    void BuildE_notDirectlyLinked_And_NotLinkedByOneComponentEdge_But_ParallelOnTheSameFacesPatchSingularEdgeIds();
    void BuildE_parallelDirectionSingularEdgeIds();
    void BuildE_orthogonalDirectionSingularEdgeIds();
    void BuildE_onTheSameFacesPatchSingularEdgeIds();
    void Build_regularSingularEdgeIds();
    void Build_irregularSingularEdgeIds();
public:
    std::vector<SingularVertex> V;
    std::vector<SingularEdge> E;
    std::vector<std::vector<size_t>> directlyLinkedSingularEdgeIds;
    std::vector<std::vector<size_t>> linkedByOneComponentEdgeSingularEdgeIds;
    std::vector<std::vector<size_t>> notDirectlyLinked_And_NotLinkedByOneComponentEdge_But_ParallelOnTheSameFacesPatchSingularEdgeIds;
    //std::vector<std::vector<size_t>> theSameDirectionSingularEdgeIds;
    std::vector<std::vector<size_t>> onTheSameFacesPatchSingularEdgeIds;
    std::vector<std::vector<size_t>> parallelDirectionSingularEdgeIds;
    std::vector<std::vector<size_t>> orthogonalDirectionSingularEdgeIds;
    std::vector<size_t> regularSingularEdgeIds;
    std::vector<size_t> irregularSingularEdgeIds;

    const BaseComplex& m_baseComplex;
};

#endif /* LIBCOTRIK_SRC_SINGULARITYGRAPH_H_ */
