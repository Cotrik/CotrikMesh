/*
 * SingularityGraph.cpp
 *
 *  Created on: Jul 25, 2017
 *      Author: cotrik
 */

#include "SingularityGraph.h"
#include <iostream>
#include <fstream>
#include <algorithm>

SingularityGraph::SingularityGraph(const BaseComplex& baseComplex)
: m_baseComplex(baseComplex)
{

}

SingularityGraph::SingularityGraph(const SingularityGraph& rhs)
: m_baseComplex(rhs.GetBaseComplex())
, V(rhs.V)
, E(rhs.E)
, directlyLinkedSingularEdgeIds(rhs.directlyLinkedSingularEdgeIds)
, linkedByOneComponentEdgeSingularEdgeIds(rhs.linkedByOneComponentEdgeSingularEdgeIds)
, parallelDirectionSingularEdgeIds(rhs.parallelDirectionSingularEdgeIds)
, orthogonalDirectionSingularEdgeIds(rhs.orthogonalDirectionSingularEdgeIds)
, onTheSameFacesPatchSingularEdgeIds(rhs.onTheSameFacesPatchSingularEdgeIds)
, regularSingularEdgeIds(rhs.regularSingularEdgeIds)
, irregularSingularEdgeIds(rhs.irregularSingularEdgeIds)
{

}

SingularityGraph::~SingularityGraph()
{
    // TODO Auto-generated destructor stub
}

const BaseComplex& SingularityGraph::GetBaseComplex() const
{
    return m_baseComplex;
}

void SingularityGraph::Build()
{
    BuildV();
    BuildE();
    BuildV_V();
    BuildV_E();
    BuildE_directlyLinkedSingularEdgeIds();
    BuildE_linkedByOneComponentEdgeSingularEdgeIds();
    BuildE_notDirectlyLinked_And_NotLinkedByOneComponentEdge_But_ParallelOnTheSameFacesPatchSingularEdgeIds();
    BuildE_parallelDirectionSingularEdgeIds();
    BuildE_orthogonalDirectionSingularEdgeIds();
    BuildE_onTheSameFacesPatchSingularEdgeIds();
    Build_regularSingularEdgeIds();
    Build_irregularSingularEdgeIds();
}

void SingularityGraph::BuildV()
{
    int id = 0;
    for (const auto& singularV : m_baseComplex.SingularityI.V) {
        SingularVertex v;
        v.id = id++;
        v.id_mesh = singularV.id_mesh;
        V.push_back(v);
    }
}

void SingularityGraph::BuildE()
{
    int id = 0;
    for (const auto& singularE : m_baseComplex.SingularityI.E) {
        SingularEdge e;
        e.id = id++;
        e.vids_link = singularE.vs_link;
        e.eids_link = singularE.es_link;
        E.push_back(e);
    }
}

void SingularityGraph::BuildV_V()
{
    for (const auto& singularE : m_baseComplex.SingularityI.E) {
        const size_t vid1 = singularE.vs_link.front();
        const size_t vid2 = singularE.vs_link.back();
        size_t singularVid1 = MAXID;
        size_t singularVid2 = MAXID;
        if (vid1 != vid2) {
            for (const auto& singularV : V)
                if (singularV.id_mesh == vid1) {
                    singularVid1 = singularV.id;
                    break;
                }
            for (const auto& singularV : V)
                if (singularV.id_mesh == vid2) {
                    singularVid2 = singularV.id;
                    break;
                }
            if (singularVid1 == MAXID || singularVid2 == MAXID) {
                std::cerr << "singular Edge Id = " << singularE.id  << " vid1 = " << vid1 << " vid2 = " << vid2 << "\n";

                std::cerr << "Error in void SingularityGraph::BuildV_V()\n";
                continue;
            }
            V[singularVid1].neighborSingularVertexIds.push_back(singularVid2);
            V[singularVid2].neighborSingularVertexIds.push_back(singularVid1);

//            const size_t vid1 = singularE.startend_Vid[0];
//            const size_t vid2 = singularE.startend_Vid[1];
//            if (vid1 == MAXID || vid2 == MAXID) {
//                std::cerr << "singular Edge Id = " << singularE.id  << " vid1 = " << vid1 << " vid2 = " << vid2 << "\n";
//                std::cerr << "Error in void SingularityGraph::BuildV_V()\n";
//                continue;
//            }
//            V[singularVid1].neighborSingularVertexIds.push_back(vid1);
//            V[singularVid2].neighborSingularVertexIds.push_back(vid2);
        }
        else // cicular
        {
            ;
        }
    }
}

void SingularityGraph::BuildV_E()
{
    for (const auto& singularE : m_baseComplex.SingularityI.E) {
        const size_t vid1 = singularE.vs_link.front();
        const size_t vid2 = singularE.vs_link.back();
        size_t singularVid1 = MAXID;
        size_t singularVid2 = MAXID;
        if (vid1 != vid2) {
            for (const auto& singularV : V)
                if (singularV.id_mesh == vid1) {
                    singularVid1 = singularV.id;
                    break;
                }
            for (const auto& singularV : V)
                if (singularV.id_mesh == vid2) {
                    singularVid2 = singularV.id;
                    break;
                }
            if (singularVid1 == MAXID || singularVid2 == MAXID) {
                std::cerr << "Error in void SingularityGraph::BuildV_E()\n";
                continue;
            }
            V[singularVid1].neighborSingularEdgeIds.push_back(singularE.id);
            V[singularVid2].neighborSingularEdgeIds.push_back(singularE.id);

//            const size_t vid1 = singularE.startend_Vid[0];
//            const size_t vid2 = singularE.startend_Vid[1];
//            V[vid1].neighborSingularEdgeIds.push_back(singularE.id);
//            V[vid2].neighborSingularEdgeIds.push_back(singularE.id);
        } else // cicular
        {
            ;
        }
    }
}
void combine(std::vector<size_t>& com, std::vector<std::vector<size_t> > &res, int n, int k, int start) {
    if (k == com.size()) {
        res.push_back(com);
        return;
    }
    for (int i = start; i < n; ++i) {
        com.push_back(i);
        combine(com, res, n, k, i + 1);
        com.pop_back();
    }
}

std::vector<std::vector<size_t>> combine(int n, int k) {
    std::vector<std::vector<size_t>> res;
    std::vector<size_t> com;
    combine(com, res, n, k, 0);
    return res;
}

void SingularityGraph::BuildE_directlyLinkedSingularEdgeIds()
{
    directlyLinkedSingularEdgeIds.clear();
    directlyLinkedSingularEdgeIds.resize(E.size(), std::vector<size_t>(E.size(), 0));
    for (const auto& v : V) {
        //auto neighborSingularEdgeIds = v.neighborSingularEdgeIds;
        std::vector<std::vector<size_t>> com = combine(v.neighborSingularEdgeIds.size(), 2);
        for (const auto& c : com) {
//            directlyLinkedSingularEdgeIds[v.neighborSingularEdgeIds[c[0]]].push_back(v.neighborSingularEdgeIds[c[1]]);
//            directlyLinkedSingularEdgeIds[v.neighborSingularEdgeIds[c[1]]].push_back(v.neighborSingularEdgeIds[c[0]]);
            directlyLinkedSingularEdgeIds[v.neighborSingularEdgeIds[c[0]]][v.neighborSingularEdgeIds[c[1]]] = 1;
            directlyLinkedSingularEdgeIds[v.neighborSingularEdgeIds[c[1]]][v.neighborSingularEdgeIds[c[0]]] = 1;
        }
    }
}

void SingularityGraph::BuildE_linkedByOneComponentEdgeSingularEdgeIds()
{
    linkedByOneComponentEdgeSingularEdgeIds.clear();
    linkedByOneComponentEdgeSingularEdgeIds.resize(E.size(), std::vector<size_t>(E.size(), 0));
    for (auto& componentEdge : m_baseComplex.componentE) {
        const auto& meshV = m_baseComplex.mesh.V;
        if (meshV.at(componentEdge.vids_link.front()).isSingularity && meshV.at(componentEdge.vids_link.back()).isSingularity) {
            std::vector<size_t> frontSingularEdgeIds;
            std::vector<size_t> backSingularEdgeIds;
            for (const auto& singularE : m_baseComplex.SingularityI.E)
                for (auto vid: singularE.vs_link)
                    if (vid == componentEdge.vids_link.front()) frontSingularEdgeIds.push_back(singularE.id);
                    else if (vid == componentEdge.vids_link.back()) backSingularEdgeIds.push_back(singularE.id);

            for (auto frontSingularEdgeId : frontSingularEdgeIds)
                for (auto backSingularEdgeId : backSingularEdgeIds) {
//                    linkedByOneComponentEdgeSingularEdgeIds.at(frontSingularEdgeId).push_back(backSingularEdgeId);
//                    linkedByOneComponentEdgeSingularEdgeIds.at(backSingularEdgeId).push_back(frontSingularEdgeId);
                    if (directlyLinkedSingularEdgeIds[frontSingularEdgeId][backSingularEdgeId] == 0) {
                    linkedByOneComponentEdgeSingularEdgeIds[frontSingularEdgeId][backSingularEdgeId] = 1;
                    linkedByOneComponentEdgeSingularEdgeIds[backSingularEdgeId][frontSingularEdgeId] = 1;
                    }
                }
        }
    }
    for (int i = 0; i < E.size(); ++i)
        linkedByOneComponentEdgeSingularEdgeIds[i][i] = 0;

    // two parallel edges in a component is considered to be linkedByOneComponentEdge;
    const auto& mesh = m_baseComplex.mesh;
    for (const auto& componentCell : m_baseComplex.componentC) {
        std::vector<size_t> singularEdgeIds;
        for (const auto componentEid : componentCell.Eids) {
            if (mesh.E.at(m_baseComplex.componentE.at(componentEid).eids_link.front()).isSingularity) {
                singularEdgeIds.push_back(mesh.E.at(m_baseComplex.componentE.at(componentEid).eids_link.front()).singularEid);
            }
        }
        if (singularEdgeIds.size() < 2) continue;
        std::vector<std::vector<size_t>> com = combine(singularEdgeIds.size(), 2);
        for (const auto& c : com) {
            const size_t seid1 = singularEdgeIds[c[0]];
            const size_t seid2 = singularEdgeIds[c[1]];
            if (directlyLinkedSingularEdgeIds[seid1][seid2] == 0) {
                linkedByOneComponentEdgeSingularEdgeIds[seid1][seid2] = 1;
                linkedByOneComponentEdgeSingularEdgeIds[seid2][seid1] = 1;
            }
        }
    }
}

void SingularityGraph::BuildE_notDirectlyLinked_And_NotLinkedByOneComponentEdge_But_ParallelOnTheSameFacesPatchSingularEdgeIds()
{

}

size_t GetTovisitedId(const std::vector<std::string>& visited)
{
    for (size_t i = 0; i != visited.size(); ++i)
        if (visited[i].empty()) return i;
    return visited.size();
}
void SingularityGraph::BuildE_parallelDirectionSingularEdgeIds()
{
    parallelDirectionSingularEdgeIds.clear();
    parallelDirectionSingularEdgeIds.resize(E.size(), std::vector<size_t>(E.size(), 0));
    orthogonalDirectionSingularEdgeIds.clear();
    orthogonalDirectionSingularEdgeIds.resize(E.size(), std::vector<size_t>(E.size(), 0));
    const auto& mesh = m_baseComplex.mesh;
    for (const auto & singularEdge : m_baseComplex.SingularityI.E) {
        for (size_t column = 0; column < directlyLinkedSingularEdgeIds.at(singularEdge.id).size(); ++column)
            if (directlyLinkedSingularEdgeIds[singularEdge.id][column] != 0) {
                orthogonalDirectionSingularEdgeIds[singularEdge.id][column] = 1;
                orthogonalDirectionSingularEdgeIds[column][singularEdge.id] = 1;
            }
    }

//    for (const auto & singularEdge : m_baseComplex.SingularityI.E) {
//        for (int i = 0; i < singularEdge.separatedFacePatchIds.size(); ++i) {
//            const std::vector<size_t>&  separatedComponentEids = singularEdge.separatedComponentEids.at(i);  // already sorted in Base Complex extraction
//            std::vector<size_t> otherSingularEdgeIdsInThePatch;
//            for (auto componentEid : separatedComponentEids)
//                if (mesh.E.at(m_baseComplex.componentE.at(componentEid).eids_link.front()).isSingularity &&
//                    mesh.E.at(m_baseComplex.componentE.at(componentEid).eids_link.front()).singularEid != singularEdge.id)
//                    otherSingularEdgeIdsInThePatch.push_back(mesh.E.at(m_baseComplex.componentE.at(componentEid).eids_link.front()).singularEid);
//            if (otherSingularEdgeIdsInThePatch.empty()) continue;
//            std::sort(separatedComponentEids.begin(), separatedComponentEids.end());
//            const auto iter = std::unique(separatedComponentEids.begin(), separatedComponentEids.end());
//            separatedComponentEids.resize(std::distance(separatedComponentEids.begin(), iter));
//            std::vector<std::string> relationForOtherSingularEdgeIdsInThePatch(otherSingularEdgeIdsInThePatch.size());
//            for (int j = 0; j < otherSingularEdgeIdsInThePatch.size(); ++j) {
//                auto otherSingularEdgeId = otherSingularEdgeIdsInThePatch.at(j);
//                if (orthogonalDirectionSingularEdgeIds[singularEdge.id][otherSingularEdgeId] != 0)
//                    relationForOtherSingularEdgeIdsInThePatch[j] = "orthogonal";
//            }
////            int next;
////            while
//        }
//    }
}

void SingularityGraph::BuildE_orthogonalDirectionSingularEdgeIds()
{
    orthogonalDirectionSingularEdgeIds.clear();
    orthogonalDirectionSingularEdgeIds.resize(E.size(), std::vector<size_t>(E.size(), 0));
    const Mesh& mesh = m_baseComplex.GetMesh();
    for (const auto& v : V) {
        //auto neighborSingularEdgeIds = v.neighborSingularEdgeIds;
        std::vector<std::vector<size_t>> com = combine(v.neighborSingularEdgeIds.size(), 2);
        for (const auto& c : com) {
            const size_t seid1 = v.neighborSingularEdgeIds[c[0]];
            const size_t seid2 = v.neighborSingularEdgeIds[c[1]];
            const SingularEdge& se1 = E.at(seid1);
            const SingularEdge& se2 = E.at(seid2);
            const bool isCircular = se1.vids_link.front() == se1.vids_link.back() || se2.vids_link.front() == se2.vids_link.back();
            if (isCircular) {
                directlyLinkedSingularEdgeIds[seid1][seid2] = 1;
                directlyLinkedSingularEdgeIds[seid1][seid1] = 1;
            } else {
                size_t eid1 = se1.eids_link.front();
                size_t eid2 = se2.eids_link.front();
                if (mesh.E.at(eid1).Vids[0] != v.id_mesh && mesh.E.at(eid1).Vids[1] != v.id_mesh)
                    eid1 = se1.eids_link.back();
                if (mesh.E.at(eid2).Vids[0] != v.id_mesh && mesh.E.at(eid2).Vids[1] != v.id_mesh)
                    eid2 = se2.eids_link.back();
                if (IsTwoEdgesOnOneFace(mesh, mesh.E.at(eid1), mesh.E.at(eid2), mesh.V.at(v.id_mesh))) {
                    orthogonalDirectionSingularEdgeIds[seid1][seid2] = 1;
                    orthogonalDirectionSingularEdgeIds[seid2][seid1] = 1;
                }
            }
        }
    }
}

void SingularityGraph::BuildE_onTheSameFacesPatchSingularEdgeIds()
{
    onTheSameFacesPatchSingularEdgeIds.clear();
    onTheSameFacesPatchSingularEdgeIds.resize(E.size(), std::vector<size_t>(E.size(), 0));
    const Mesh& mesh = m_baseComplex.GetMesh();
    int seid = 0;
    for (const auto& singularityEdge : m_baseComplex.SingularityI.E) {
        for (const auto& componentEids : singularityEdge.separatedComponentEids) {
            std::vector<size_t> singularEdgeIds;
            for (const auto componentEid : componentEids) {
                if (mesh.E.at(m_baseComplex.componentE.at(componentEid).eids_link.front()).isSingularity) {
                    singularEdgeIds.push_back(mesh.E.at(m_baseComplex.componentE.at(componentEid).eids_link.front()).singularEid);
                }
            }
            std::sort(singularEdgeIds.begin(), singularEdgeIds.end());
            singularEdgeIds.resize(std::distance(singularEdgeIds.begin(), std::unique(singularEdgeIds.begin(), singularEdgeIds.end())));
            //std::vector<std::vector<size_t>> combinations = combine(v.neighborSingularEdgeIds.size(), 2);
            for (const auto singularEdgeId : singularEdgeIds) {
                onTheSameFacesPatchSingularEdgeIds[seid][singularEdgeId] = 1;
                onTheSameFacesPatchSingularEdgeIds[singularEdgeId][seid] = 1;
            }
        }
        ++seid;
    }
}

void SingularityGraph::Build_regularSingularEdgeIds()
{

}

void SingularityGraph::Build_irregularSingularEdgeIds()
{

}

static void GetQuadVertexIds(size_t i, size_t j, int n, size_t quadVids[]){
    quadVids[0] = i * (n + 1) + j;
    quadVids[1] = i * (n + 1) + j + 1;
    quadVids[2] = (i + 1) * (n + 1) + j + 1;
    quadVids[3] = (i + 1) * (n + 1) + j;
}

void SingularityGraph::WriteMatrixVTK(const char* filename, const std::vector<std::vector<size_t>>& m)
{
    const size_t n = m.size();
    const size_t numOfVertices = (n + 1) * (n + 1);
    const size_t numOfCells = n * n;
    std::vector<glm::vec2> V(numOfVertices);
    for (size_t i = 0; i <= n; ++i)
        for (size_t j = 0; j <= n; ++j)
            V[i * (n + 1) + j] = glm::vec2(i, j);

    std::ofstream ofs(filename);
    ofs << "# vtk DataFile Version 3.0\n"
        << filename << "\n"
        << "ASCII\n\n"
        << "DATASET POLYDATA\n";

    ofs << "POINTS " << numOfVertices << " double\n";
    for (size_t i = 0; i < V.size(); i++)
        ofs << V.at(i).x << " " << V.at(i).y << " 0\n";

    ofs << "POLYGONS " << numOfCells << " " << 5 * numOfCells << "\n";
    size_t quadVids[4] = {0};
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j) {
            GetQuadVertexIds(i, j, n, quadVids);
            ofs << "4 " << quadVids[0] << " " << quadVids[1] << " " << quadVids[2] << " " << quadVids[3] << "\n";
        }

    ofs << "CELL_DATA " << numOfCells << "\n"
        << "SCALARS " << "label" << " int 1\n"
        << "LOOKUP_TABLE default\n";
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            ofs << m[i][j] << "\n";
}

void SingularityGraph::WriteMatrixMat(const char* filename, const std::vector<std::vector<size_t>>& m)
{
    std::ofstream ofs(filename);
    const size_t n = m.size();
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j)
            ofs << m[i][j] << "\t";
        ofs << "\n";
    }
}
