/*
 * PolyLine.cpp
 *
 *  Created on: Nov 11, 2016
 *      Author: cotrik
 */

#include "PolyLine.h"
#include "Mesh.h"
#include "FrameField.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>

PolyLine::PolyLine()
{
    // TODO Auto-generated constructor stub

}

PolyLine::~PolyLine()
{
    // TODO Auto-generated destructor stub
}

//PolyLines::PolyLines()
//{
//    // TODO Auto-generated constructor stub
//
//}
PolyLines::PolyLines(const Mesh& mesh, const FrameField& framefield)
: mesh(mesh)
, framefield(framefield)
{
    // TODO Auto-generated constructor stub

}
PolyLines::~PolyLines()
{
    // TODO Auto-generated destructor stub
}

size_t PolyLine::BuildFrom(const Mesh& mesh, const FrameField& framefield, size_t startFrameNodeId, size_t startFrameEdgeId)
{
    while (true)
    {
        size_t nextFrameNodeId = framefield.frameEdges.at(startFrameEdgeId).Vids.at(1);
        if (framefield.frameEdges.at(startFrameEdgeId).Vids.at(1) == startFrameNodeId)
            nextFrameNodeId = framefield.frameEdges.at(startFrameEdgeId).Vids.at(0);

        //const bool isBoundaryStartFramNode = framefield.frameNodes.at(startFrameNodeId).isBoundary;
        Vids.push_back(startFrameNodeId);
        Eids.push_back(startFrameEdgeId);

        const bool isBoundaryNextFramNode = framefield.frameNodes.at(nextFrameNodeId).isBoundary;
        if (nextFrameNodeId == Vids.at(0) || isBoundaryNextFramNode)
        {
            Vids.push_back(nextFrameNodeId);
            break;
        }

        const size_t cellId = nextFrameNodeId;
        const size_t faceId = framefield.frameEdges.at(startFrameEdgeId).id;
        const size_t nextFrameEdgeId = Util::GetOppositeFaceId(mesh, cellId, faceId);

        startFrameNodeId = nextFrameNodeId;
        startFrameEdgeId = nextFrameEdgeId;
    }

    return Eids.at(Eids.size() - 1);
}

void PolyLines::BuildFrom(const Mesh& mesh)
{

}

//void PolyLines::BuildFrom(const Mesh& mesh, const FrameField& framefield)
void PolyLines::Build()
{
    size_t polylineId = 0;
    std::vector<size_t> visitedFrameEdgeIds;
    std::vector<size_t> edgesLabels(framefield.frameEdges.size(), MAXID);
    for (size_t i = 0; i < framefield.frameEdges.size(); i++){
        const size_t startFrameEdgeId = framefield.frameEdges.at(i).id;
        const Frame& frameNode1 = framefield.frameNodes.at(framefield.frameEdges.at(i).Vids.at(0));
        const Frame& frameNode2 = framefield.frameNodes.at(framefield.frameEdges.at(i).Vids.at(1));
        if (frameNode1.isBoundary || frameNode2.isBoundary){
            const Frame& frameNode = frameNode1.isBoundary ? frameNode1 : frameNode2;
            const size_t startFrameNodeId = frameNode1.isBoundary ? frameNode1.id : frameNode2.id;
            for (size_t j = 0; j < frameNode.N_Eids.size(); j++){
                const size_t startFrameEdgeId = frameNode.N_Eids.at(j);
                if (!Util::Find(visitedFrameEdgeIds, startFrameEdgeId)){
                    visitedFrameEdgeIds.push_back(startFrameEdgeId);
                    PolyLine polyline;
                    polyline.id = polylineId++;
                    const size_t endFrameEdgeId = polyline.BuildFrom(mesh, framefield, startFrameNodeId, startFrameEdgeId);
//                    std::vector<FrameEdge>& FE = (std::vector<FrameEdge>&)framefield.frameEdges;
//                    for (size_t k = 0; k < polyline.Eids.size(); k++)
//                        FE.at(polyline.Eids.at(k)).polylineId = polyline.id;
                    std::vector<FrameEdge>& FE = (std::vector<FrameEdge>&) framefield.frameEdges;
                    for (size_t k = 0; k < polyline.Eids.size(); k++)
                        edgesLabels.at(polyline.Eids.at(k)) = polyline.id;
                    polyLines.push_back(polyline);
                    if (!Util::Find(visitedFrameEdgeIds, endFrameEdgeId)){
                        visitedFrameEdgeIds.push_back(endFrameEdgeId);
                    }
                }
            }
        }
    }

    for (size_t i = 0; i < framefield.frameEdges.size(); i++){
        const size_t startFrameEdgeId = framefield.frameEdges.at(i).id;
        const Frame& frameNode1 = framefield.frameNodes.at(framefield.frameEdges.at(i).Vids.at(0));
        const Frame& frameNode2 = framefield.frameNodes.at(framefield.frameEdges.at(i).Vids.at(1));
        if (edgesLabels.at(startFrameEdgeId) == MAXID){
            const Frame& frameNode = frameNode1;
            const size_t startFrameNodeId = frameNode1.id;
            //for (size_t j = 0; j < frameNode.N_Eids.size(); j++)
            {
                //const size_t startFrameEdgeId = frameNode.N_Eids.at(j);
                if (!Util::Find(visitedFrameEdgeIds, startFrameEdgeId)){
                    visitedFrameEdgeIds.push_back(startFrameEdgeId);
                    PolyLine polyline;
                    polyline.id = polylineId++;
                    const size_t endFrameEdgeId = polyline.BuildFrom(mesh, framefield, startFrameNodeId, startFrameEdgeId);
                    std::vector<FrameEdge>& FE = (std::vector<FrameEdge>&) framefield.frameEdges;
                    for (size_t k = 0; k < polyline.Eids.size(); k++)
                        edgesLabels.at(polyline.Eids.at(k)) = polyline.id;
                    polyLines.push_back(polyline);
                    if (!Util::Find(visitedFrameEdgeIds, endFrameEdgeId)){
                        visitedFrameEdgeIds.push_back(endFrameEdgeId);
                    }
                }
            }
        }
    }
}

void PolyLines::WriteFile(const char* filename)
{
    std::ofstream ofs(filename);
    ofs << "# vtk DataFile Version 2.0" << std::endl
        << filename << std::endl
        << "ASCII" << std::endl << std::endl
        << "DATASET POLYDATA" << std::endl;

    ofs << "POINTS " << framefield.frameNodes.size() << " double" << std::endl;
    ofs << std::fixed << std::setprecision(7);
    for (size_t i = 0; i < framefield.frameNodes.size(); i++)
        ofs << framefield.frameNodes.at(i).x << " " << framefield.frameNodes.at(i).y << " " << framefield.frameNodes.at(i).z << std::endl;

    ofs << "VERTICES " << framefield.frameNodes.size() << " " << 2 * framefield.frameNodes.size() << std::endl;
    for (size_t i = 0; i < framefield.frameNodes.size(); i++){
        ofs << 1 << " " << i << std::endl;
    }

    size_t sum = 0;
    for (size_t i = 0; i < polyLines.size(); i++)
        sum += polyLines.at(i).Vids.size() + 1;

    ofs << "LINES " << polyLines.size() << " " << sum << std::endl;
    for (size_t i = 0; i < polyLines.size(); i++){
        const PolyLine& polyline = polyLines.at(i);
        ofs << polyline.Vids.size();
        for (size_t j = 0; j < polyline.Vids.size(); j++)
            ofs << " " << polyline.Vids.at(j);
        ofs << std::endl;
    }

    ofs << "CELL_DATA " << polyLines.size() + framefield.frameNodes.size() << std::endl
            << "SCALARS " << " Label" << " int 1" << std::endl
            << "LOOKUP_TABLE default" << std::endl;
    for (size_t i = 0; i < framefield.frameNodes.size(); i++)
        ofs << 0 << std::endl;

    for (size_t i = 0; i < polyLines.size(); i++)
        ofs << i << std::endl;

    ofs << "SCALARS " << " Color" << " int 1" << std::endl
        << "LOOKUP_TABLE default" << std::endl;
    for (size_t i = 0; i < framefield.frameNodes.size(); i++)
        ofs << 0 << std::endl;
    for (size_t i = 0; i < polyLines.size(); i++)
        ofs << i % 6 << std::endl;
}
