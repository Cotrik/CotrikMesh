/*
 * EdgeLines.cpp
 *
 *  Created on: Nov 25, 2016
 *      Author: cotrik
 */

#include "Mesh.h"
#include "EdgeLines.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>

EdgeLine::EdgeLine()
: id (MAXID)
{
    // TODO Auto-generated constructor stub

}

EdgeLine::~EdgeLine()
{
    // TODO Auto-generated destructor stub
}

//EdgeLines::EdgeLines()
//{
//    // TODO Auto-generated constructor stub
//
//}
EdgeLines::EdgeLines(const Mesh& mesh)
: mesh(mesh)
{
    // TODO Auto-generated constructor stub

}

EdgeLines::~EdgeLines()
{
    // TODO Auto-generated destructor stub
}

size_t EdgeLine::BuildFrom(const Mesh& mesh, size_t startVId, size_t startEdgeId)
{
    size_t prevEdgeId = startEdgeId;
    while (true)
    {
        size_t nextVId = mesh.E.at(startEdgeId).Vids.at(1);
        if (mesh.E.at(startEdgeId).Vids.at(1) == startVId)
            nextVId = mesh.E.at(startEdgeId).Vids.at(0);

        Vids.push_back(startVId);
        Eids.push_back(startEdgeId);

        const bool isBoundaryNextV = mesh.V.at(nextVId).isBoundary;
        const bool isBoundaryStartV = mesh.V.at(startVId).isBoundary;
        const bool isSingularityNextV = mesh.V.at(nextVId).isSingularity;
        const bool isBoundaryBeginV = mesh.V.at(Vids.at(0)).isBoundary;
        const bool isSingularityBeginV = mesh.V.at(Vids.at(0)).isSingularity;
        if (nextVId == Vids.at(0) || (isBoundaryNextV && !isBoundaryStartV) || (isSingularityNextV && isBoundaryStartV)
                || (isSingularityNextV && (isBoundaryBeginV || isSingularityBeginV)) )
        {
            Vids.push_back(nextVId);
            break;
        }

        size_t nextEdgeId = mesh.E.at(startEdgeId).consecutiveEids.at(0);
        if (nextEdgeId == prevEdgeId)
            nextEdgeId = mesh.E.at(startEdgeId).consecutiveEids.at(1);

        prevEdgeId = startEdgeId;
        startVId = nextVId;
        startEdgeId = nextEdgeId;
    }

    return Eids.at(Eids.size() - 1);
}

size_t EdgeLine::BuildFrom(const Mesh& mesh, size_t startVId, size_t startEdgeId, size_t startFaceId)
{
    size_t prevEdgeId = startEdgeId;
    while (true)
    {
        const Face& startFace = mesh.F.at(startFaceId);
        size_t cellId = startFace.N_Cids.at(0);
        const size_t nextFaceId = GetoppositeFaceId(mesh, cellId, startFaceId);
        const Face& nextFace = mesh.F.at(nextFaceId);

        size_t nextVId = mesh.E.at(startEdgeId).Vids.at(1);
        if (mesh.E.at(startEdgeId).Vids.at(1) == startVId)
            nextVId = mesh.E.at(startEdgeId).Vids.at(0);

        Vids.push_back(startVId);
        Eids.push_back(startEdgeId);

        const bool isBoundaryNextV = mesh.V.at(nextVId).isBoundary;
        const bool isSingularityNextV = mesh.V.at(nextVId).isSingularity;
        if (nextVId == Vids.at(0) || isBoundaryNextV || isSingularityNextV)
        {
            Vids.push_back(nextVId);
            break;
        }

//        size_t layerj = 0;
//        for (size_t j = 0; j < nextFace.N_Ortho_4Eids.size(); j++)
//            for (size_t k = 0; k < nextFace.N_Ortho_4Eids.at(j).size(); k++)
//                if (nextFace.N_Ortho_4Eids.at(j).at(k) == startEdgeId) {
//                    layerj = j == 0 ? 1 : 0;
//                    break;
//                }

        const Edge& startEdge = mesh.E.at(startEdgeId);
        size_t nextEdgeId = startEdge.consecutiveEids.at(0);
//
//        for (size_t k = 0; k < nextFace.N_Ortho_4Eids.at(layerj).size(); k++)
//            for (size_t i = 0; i < mesh.E.at(startEdgeId).consecutiveEids.size(); i++)
//                if (mesh.E.at(startEdgeId).consecutiveEids.at(i) == nextFace.N_Ortho_4Eids.at(layerj).at(k)) {
//                    nextEdgeId = mesh.E.at(startEdgeId).consecutiveEids.at(i);
//                    break;
//                }

        if (nextEdgeId == prevEdgeId)
            nextEdgeId = mesh.E.at(startEdgeId).consecutiveEids.at(1);

        prevEdgeId = startEdgeId;
        startVId = nextVId;
        startEdgeId = nextEdgeId;
    }

    return Eids.at(Eids.size() - 1);
}

//void EdgeLines::Build()
//{
//    size_t edgelineId = 0;
//    std::vector<size_t> visitedEdgeIds;
//    std::vector<size_t> edgesLabels(mesh.E.size(), MAXID);
//    for (size_t i = 0; i < mesh.V.size(); i++) {
//        const size_t startVId = mesh.V.at(i).id;
//        const Vertex& v = mesh.V.at(startVId);
//        if (v.isBoundary) {
//            for (size_t j = 0; j < v.N_Eids.size(); j++) {
//                const size_t startEdgeId = v.N_Eids.at(j);
//                if (!Find(visitedEdgeIds, startEdgeId)) {
//                    visitedEdgeIds.push_back(startEdgeId);
//                    EdgeLine edgeline;
//                    edgeline.id = edgelineId++;
//                    const size_t endEdgeId = edgeline.BuildFrom(mesh, startVId, startEdgeId);
//                    for (size_t k = 0; k < edgeline.Eids.size(); k++)
//                        edgesLabels.at(edgeline.Eids.at(k)) = edgeline.id;
//                    edgeLines.push_back(edgeline);
//                    if (!Find(visitedEdgeIds, endEdgeId)) {
//                        visitedEdgeIds.push_back(endEdgeId);
//                    }
//                }
//            }
//        }
//    }
//
//    for (size_t i = 0; i < mesh.E.size(); i++) {
//        const size_t startEdgeId = mesh.E.at(i).id;
//        const Edge& startEdge = mesh.E.at(i);
//        if (edgesLabels.at(startEdgeId) == MAXID) {
//            const size_t v1Id = startEdge.Vids[0];
//            const size_t v2Id = startEdge.Vids[1];
//            const Vertex& v1 = mesh.V.at(v1Id);
//            const size_t startVId = v1.isSingularity ? v1Id : v2Id;
//            const Vertex& v = mesh.V.at(startVId);
//            if (v.isSingularity)
//            if (!Find(visitedEdgeIds, startEdgeId)) {
//                visitedEdgeIds.push_back(startEdgeId);
//                EdgeLine edgeline;
//                edgeline.id = edgelineId++;
//                const size_t endEdgeId = edgeline.BuildFrom(mesh, startVId, startEdgeId);
//                for (size_t k = 0; k < edgeline.Eids.size(); k++)
//                    edgesLabels.at(edgeline.Eids.at(k)) = edgeline.id;
//                edgeLines.push_back(edgeline);
//                if (!Find(visitedEdgeIds, endEdgeId)) {
//                    visitedEdgeIds.push_back(endEdgeId);
//                }
//            }
//        }
//    }
//}

void EdgeLines::Build()
{
    size_t edgelineId = 0;
    std::vector<size_t> visitedEdgeIds;
    std::vector<size_t> edgesLabels(mesh.E.size(), MAXID);
    for (size_t i = 0; i < mesh.F.size(); i++) {
        const Face& face = mesh.F.at(i);
        if (!face.isBoundary)
            continue;
        for (size_t j = 0; j < face.N_Ortho_4Eids.size(); j++) {
            for (size_t k = 0; k < face.N_Ortho_4Eids.at(j).size(); k++) {
                const size_t startEdgeId = face.N_Ortho_4Eids.at(j).at(k);
                const Edge& startEdge = mesh.E.at(startEdgeId);
                const size_t startVId = mesh.V.at(startEdge.Vids.at(0)).isBoundary ? startEdge.Vids.at(0) : startEdge.Vids.at(1);
                if (!Find(visitedEdgeIds, startEdgeId)) {
                    visitedEdgeIds.push_back(startEdgeId);
                    EdgeLine edgeline;
                    edgeline.id = edgelineId++;
                    const size_t endEdgeId = edgeline.BuildFrom(mesh, startVId, startEdgeId, face.id);
                    for (size_t k = 0; k < edgeline.Eids.size(); k++)
                        edgesLabels.at(edgeline.Eids.at(k)) = edgeline.id;
                    edgeLines.push_back(edgeline);
                    if (!Find(visitedEdgeIds, endEdgeId)) {
                        visitedEdgeIds.push_back(endEdgeId);
                    }
                }
            }
        }
    }

    for (size_t i = 0; i < mesh.E.size(); i++) {
        const size_t startEdgeId = mesh.E.at(i).id;
        const Edge& startEdge = mesh.E.at(i);
        if (edgesLabels.at(startEdgeId) == MAXID && !startEdge.isBoundary) {
            const size_t v1Id = startEdge.Vids[0];
            const size_t v2Id = startEdge.Vids[1];
            const Vertex& v1 = mesh.V.at(v1Id);
            const size_t startVId = v1.isSingularity ? v1Id : v2Id;
            const Vertex& v = mesh.V.at(startVId);
            if (v.isSingularity)
            if (!Find(visitedEdgeIds, startEdgeId)) {
                visitedEdgeIds.push_back(startEdgeId);
                EdgeLine edgeline;
                edgeline.id = edgelineId++;
                const size_t endEdgeId = edgeline.BuildFrom(mesh, startVId, startEdgeId);
                for (size_t k = 0; k < edgeline.Eids.size(); k++)
                    edgesLabels.at(edgeline.Eids.at(k)) = edgeline.id;
                edgeLines.push_back(edgeline);
                if (!Find(visitedEdgeIds, endEdgeId)) {
                    visitedEdgeIds.push_back(endEdgeId);
                }
            }
        }
    }
}

void EdgeLines::WriteFile(const char* filename)
{
    std::ofstream ofs(filename);
    ofs << "# vtk DataFile Version 2.0" << std::endl
        << filename << std::endl
        << "ASCII" << std::endl << std::endl
        << "DATASET POLYDATA" << std::endl;

    ofs << "POINTS " << mesh.V.size() << " double" << std::endl;
    ofs << std::fixed << std::setprecision(7);
    for (size_t i = 0; i < mesh.V.size(); i++)
        ofs << mesh.V.at(i).x << " " << mesh.V.at(i).y << " " << mesh.V.at(i).z << std::endl;

    size_t numOfBoundaryV;
//    ofs << "VERTICES " << mesh.V.size() << " " << 2 * mesh.V.size() << std::endl;
//    for (size_t i = 0; i < mesh.V.size(); i++){
//        if (mesh.V[i].isBoundary){
//            ofs << 1 << " " << i << std::endl;
//            numOfBoundaryV++;
//        }
//    }

    size_t sum = 0;
    for (size_t i = 0; i < edgeLines.size(); i++)
        sum += edgeLines.at(i).Vids.size() + 1;

    ofs << "LINES " << edgeLines.size() << " " << sum << std::endl;
    for (size_t i = 0; i < edgeLines.size(); i++){
        const EdgeLine& edgeline = edgeLines.at(i);
        ofs << edgeline.Vids.size();
        for (size_t j = 0; j < edgeline.Vids.size(); j++)
            ofs << " " << edgeline.Vids.at(j);
        ofs << std::endl;
    }

    ofs << "CELL_DATA " << edgeLines.size() + numOfBoundaryV << std::endl
            << "SCALARS " << " Label" << " int 1" << std::endl
            << "LOOKUP_TABLE default" << std::endl;
//    for (size_t i = 0; i < mesh.V.size(); i++)
//        if (mesh.V[i].isBoundary)
//            ofs << 0 << std::endl;
    for (size_t i = 0; i < edgeLines.size(); i++)
        ofs << i << std::endl;
}
