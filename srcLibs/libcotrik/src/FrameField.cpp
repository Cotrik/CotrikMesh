/*
 * FrameField.cpp
 *
 *  Created on: Nov 10, 2016
 *      Author: cotrik
 */

#include "FrameField.h"
#include "Frame.h"
#include "PolyLine.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>

FrameField::FrameField()
{
    // TODO Auto-generated constructor stub

}

FrameField::~FrameField()
{
    // TODO Auto-generated destructor stub
}


//void FrameField::BuildFrom(const Mesh& mesh)
//{
//    frameNodes.resize(mesh.C.size());
//    for (size_t i = 0; i < mesh.C.size(); i++)
//        frameNodes.at(i).BuildFrom(mesh, i);
//
//    size_t innerFaceSize = 0;
//    for (size_t i = 0; i < mesh.F.size(); i++)
//        if (!mesh.F.at(i).isBoundary)
//            innerFaceSize++;
//
//    frameEdges.resize(innerFaceSize);
//    size_t id = 0;
//    for (size_t i = 0; i < mesh.F.size(); i++){
//        if (!mesh.F.at(i).isBoundary){
//            frameEdges.at(id).id = id;
//            frameEdges.at(id).Vids = mesh.F.at(i).N_Cids;
//            id++;
//        }
//    }
//}

void FrameField::BuildFrom(const Mesh& mesh)
{
    size_t surfaceSize = 0;
    size_t innerFaceSize = 0;
    for (size_t i = 0; i < mesh.F.size(); i++)
        if (mesh.F.at(i).isBoundary)
            surfaceSize++;
        else
            innerFaceSize++;

    frameNodes.resize(surfaceSize + mesh.C.size());

    for (size_t i = 0; i < mesh.C.size(); i++)
        frameNodes.at(i).BuildFrom(mesh, i, false);
//    for (size_t i = 0; i < surfaceSize; i++)
//        frameNodes.at(i + mesh.C.size()).BuildFrom(mesh, i, true);

    size_t cellId = mesh.C.size();
    for (size_t i = 0; i < mesh.F.size(); i++)
        if (mesh.F.at(i).isBoundary){
            frameNodes.at(cellId).BuildFrom(mesh, cellId, true, i);
            cellId++;
        }

    frameEdges.resize(mesh.F.size());
    size_t innerFaceId = 0;
    size_t surfaceId = 0;
//    for (size_t i = 0; i < mesh.F.size(); i++){
//        if (!mesh.F.at(i).isBoundary){
//            frameEdges.at(innerFaceId).id = innerFaceId;
//            frameEdges.at(innerFaceId).Vids = mesh.F.at(i).N_Cids;
//            innerFaceId++;
//        }
//        else{
//            frameEdges.at(innerFaceSize + surfaceId).id = innerFaceSize + surfaceId;
//            frameEdges.at(innerFaceSize + surfaceId).Vids.push_back(mesh.F.at(i).N_Cids.at(0));
//            frameEdges.at(innerFaceSize + surfaceId).Vids.push_back(mesh.C.size() + surfaceId);
//            surfaceId++;
//        }
//    }
    for (size_t i = 0; i < mesh.F.size(); i++){
        frameEdges.at(i).id = mesh.F.at(i).id;
        if (!mesh.F.at(i).isBoundary)
            frameEdges.at(i).Vids = mesh.F.at(i).N_Cids;
        else{
            frameEdges.at(i).Vids.push_back(mesh.C.size() + surfaceId++);
            frameEdges.at(i).Vids.push_back(mesh.F.at(i).N_Cids.at(0));
        }
    }
}

void FrameField::SetPolylineIds(const PolyLines& polyLines)
{
    for (size_t i = 0; i < polyLines.polyLines.size(); i++) {
        const PolyLine& polyline = polyLines.polyLines.at(i);
        for (size_t j = 0; j < polyline.Eids.size(); j++) {
            FrameEdge& frameEdge = frameEdges.at(polyline.Eids.at(j));
            frameEdge.SetPolylineId(polyline.id);
        }
    }
}

void FrameField::BuildOrthogonal4EidsForEachFrame()
{
    std::vector<size_t> ortho4Eids(4);
    std::vector<size_t> ortho4Vids(4);
    for (size_t i = 0; i < frameNodes.size(); i++) {
        Frame& frame = frameNodes.at(i);
        if (!frame.isBoundary)
        for (size_t j = 0; j < frame.N_Eids.size(); j++) {
            const FrameEdge& frameEdge = frameEdges.at(frame.N_Eids.at(j));
            int count = 0;
            for (size_t k = 0; k < frame.N_Eids.size(); k++) {
                const size_t frameEdgeId = frame.N_Eids.at(k);
                const size_t frameNodeId = frame.N_Vids.at(k);
                const FrameEdge& N_frameEdge = frameEdges.at(frameEdgeId);
                if (N_frameEdge.polylineId != frameEdge.polylineId){
                    ortho4Eids[count] = frameEdgeId;
                    ortho4Vids[count] = frameNodeId;
                    count++;
                }
            }
            frame.ortho4Eids.push_back(ortho4Eids);
            frame.ortho4Vids.push_back(ortho4Vids);
        }
    }
}

void FrameField::WriteFile(const char* filename)
{
    std::ofstream ofs(filename);
    ofs << "# vtk DataFile Version 2.0" << std::endl
        << filename << std::endl
        << "ASCII" << std::endl << std::endl
        << "DATASET POLYDATA" << std::endl;

    ofs << "POINTS " << frameNodes.size() << " double" << std::endl;
    ofs << std::fixed << std::setprecision(7);
    for (size_t i = 0; i < frameNodes.size(); i++)
        ofs << frameNodes.at(i).x << " " << frameNodes.at(i).y << " " << frameNodes.at(i).z << std::endl;

    ofs << "LINES " << frameEdges.size() << " " << 3 * frameEdges.size() << std::endl;
    for (size_t i = 0; i < frameEdges.size(); i++){
        ofs << "2 " << frameEdges.at(i).Vids[0] << " " << frameEdges.at(i).Vids[1] << std::endl;
    }

    std::string filename2 = std::string(filename).substr(0, std::string(filename).size() - 15) + ".ff";
    std::ofstream ofs2(filename2);
    ofs2 << frameEdges.size() << " 4" << std::endl;
    ofs2 << std::fixed << std::setprecision(7);
    for (size_t i = 0; i < frameEdges.size(); i++){
        ofs2 << frameNodes.at(frameEdges.at(i).Vids[0]).x << " " << frameNodes.at(frameEdges.at(i).Vids[0]).y << " " << frameNodes.at(frameEdges.at(i).Vids[0]).z << std::endl;
        ofs2 << frameNodes.at(frameEdges.at(i).Vids[1]).x << " " << frameNodes.at(frameEdges.at(i).Vids[1]).y << " " << frameNodes.at(frameEdges.at(i).Vids[1]).z << std::endl;
        // ofs << "2 " << frameEdges.at(i).Vids[0] << " " << frameEdges.at(i).Vids[1] << std::endl;
    }
}

void FrameField::WriteFile(const char* filename, const std::vector<size_t>& frameIds)
{
    std::vector<size_t> frameEdgeIds;
    for (size_t i = 0; i < frameIds.size(); i++) {
        const Frame& frame = frameNodes.at(frameIds.at(i));
        std::copy(frame.N_Eids.begin(), frame.N_Eids.end(), back_inserter(frameEdgeIds));
    }
    std::sort(frameEdgeIds.begin(), frameEdgeIds.end());
    std::vector<size_t>::iterator iter = std::unique(frameEdgeIds.begin(), frameEdgeIds.end());
    frameEdgeIds.resize(std::distance(frameEdgeIds.begin(), iter));

    size_t numberOfBoundaryVertices = 0;
    for (size_t i = 0; i < frameEdgeIds.size(); i++) {
        const FrameEdge& frameEdge = frameEdges.at(frameEdgeIds.at(i));
        for (size_t j = 0; j < frameEdge.Vids.size(); j++) {
            const Frame& frame = frameNodes.at(frameEdge.Vids.at(j));
            if (frame.isBoundary)
                numberOfBoundaryVertices++;
        }
    }

    std::ofstream ofs(filename);
    ofs << "# vtk DataFile Version 2.0" << std::endl
        << filename << std::endl
        << "ASCII" << std::endl << std::endl
        << "DATASET POLYDATA" << std::endl;

    ofs << "POINTS " << frameNodes.size() << " double" << std::endl;
    ofs << std::fixed << std::setprecision(7);
    for (size_t i = 0; i < frameNodes.size(); i++)
        ofs << frameNodes.at(i).x << " " << frameNodes.at(i).y << " " << frameNodes.at(i).z << std::endl;

    ofs << "VERTICES " << numberOfBoundaryVertices << " " << 2 * numberOfBoundaryVertices << std::endl;
    for (size_t i = 0; i < frameEdgeIds.size(); i++) {
        const FrameEdge& frameEdge = frameEdges.at(frameEdgeIds.at(i));
        for (size_t j = 0; j < frameEdge.Vids.size(); j++) {
            const Frame& frame = frameNodes.at(frameEdge.Vids.at(j));
            if (frame.isBoundary)
                ofs << 1 << " " << frameEdge.Vids.at(j) << std::endl;
        }
    }

    ofs << "LINES " << frameEdgeIds.size() << " " << 3 * frameEdgeIds.size() << std::endl;
    for (size_t i = 0; i < frameEdgeIds.size(); i++){
        const FrameEdge& frameEdge = frameEdges.at(frameEdgeIds.at(i));
        ofs << "2 " << frameEdge.Vids[0] << " " << frameEdge.Vids[1] << std::endl;
    }

    ofs << "CELL_DATA " << numberOfBoundaryVertices + frameEdgeIds.size() << std::endl
        << "SCALARS " << " Label" << " int 1" << std::endl
        << "LOOKUP_TABLE default" << std::endl;
    for (size_t i = 0; i < frameEdgeIds.size(); i++) {
        const FrameEdge& frameEdge = frameEdges.at(frameEdgeIds.at(i));
        for (size_t j = 0; j < frameEdge.Vids.size(); j++) {
            const Frame& frame = frameNodes.at(frameEdge.Vids.at(j));
            if (frame.isBoundary)
                ofs << 0 << std::endl;
        }
    }
    for (size_t i = 0; i < frameEdgeIds.size(); i++)
        ofs << frameEdges.at(frameEdgeIds.at(i)).polylineId << std::endl;
}

void FrameField::WriteColorFile(const char* filename)
{
    std::ofstream ofs(filename);
    ofs << "# vtk DataFile Version 2.0" << std::endl
        << filename << std::endl
        << "ASCII" << std::endl << std::endl
        << "DATASET POLYDATA" << std::endl;

    ofs << "POINTS " << frameNodes.size() << " double" << std::endl;
    ofs << std::fixed << std::setprecision(7);
    for (size_t i = 0; i < frameNodes.size(); i++)
        ofs << frameNodes.at(i).x << " " << frameNodes.at(i).y << " " << frameNodes.at(i).z << std::endl;

    ofs << "LINES " << frameEdges.size() << " " << 3 * frameEdges.size() << std::endl;
    for (size_t i = 0; i < frameEdges.size(); i++){
        ofs << "2 " << frameEdges.at(i).Vids[0] << " " << frameEdges.at(i).Vids[1] << std::endl;
    }

    ofs << "CELL_DATA " << frameEdges.size() << std::endl
            << "SCALARS " << " Label" << " int 1" << std::endl
            << "LOOKUP_TABLE default" << std::endl;
    for (size_t i = 0; i < frameEdges.size(); i++)
        ofs << frameEdges.at(i).polylineId << std::endl;
}
