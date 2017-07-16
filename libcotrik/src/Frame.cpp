/*
 * Frame.cpp
 *
 *  Created on: Nov 10, 2016
 *      Author: cotrik
 */

#include "Frame.h"

Frame::Frame()
{
    // TODO Auto-generated constructor stub

}

Frame::~Frame()
{
    // TODO Auto-generated destructor stub
}

//Frame::Frame(size_t vnum)
//: Cell(vnum)
//{
//
//}
//Frame::Frame(size_t vnum, size_t eNum, size_t fNum)
//: Cell (vnum, eNum, fNum)
//{
//
//}
//Frame::Frame(const std::vector<size_t>& Vids)
//: Cell(Vids)
//{
//}
//
//Frame::Frame(const std::vector<size_t>& Vids, const std::vector<size_t> Eids, const std::vector<size_t> Fids)
//: Cell(Vids, Eids, Fids)
//{
//}

//void Frame::BuildFrom(const Mesh& mesh, const size_t cellId, bool isBoundary)
//{
//    this->isBoundary = isBoundary;
//    if (!isBoundary){
//        id = cellId;
//        const Cell& cell = mesh.C.at(cellId);
//        const size_t facesSize = cell.N_Fids.size();
//
//        // Compute center x,y,z coordinates
//        for (size_t i = 0; i < facesSize; i++)
//        {
//            const Face& face = mesh.F.at(cell.N_Fids.at(i));
//            glm::vec3 faceCenter(0.0, 0.0, 0.0);
//            const size_t faceVidsSize = face.Vids.size();
//            for (size_t j = 0; j < faceVidsSize; j++)
//            {
//                const Vertex& v = mesh.V.at(face.Vids.at(j));
//                faceCenter.x += v.x;
//                faceCenter.y += v.y;
//                faceCenter.z += v.z;
//            }
//            x += faceCenter.x / faceVidsSize;
//            y += faceCenter.y / faceVidsSize;
//            z += faceCenter.z / faceVidsSize;
//        }
//        x /= facesSize;
//        y /= facesSize;
//        z /= facesSize;
//        // Compute N_Eids and N_Vids
//        for (size_t i = 0; i < facesSize; i++)
//        {
//            const Face& face = mesh.F.at(cell.N_Fids.at(i));
//            const size_t faceVidsSize = face.Vids.size();
//
//            for (size_t j = 0; j < face.N_Cids.size(); j++)
//                if (face.N_Cids.at(j) != cellId){
//                    N_Vids.push_back(face.N_Cids.at(j));
//                    break;
//                }
//            N_Eids.push_back(face.id);
//        }
//    }
//    else{
//        // cellId is faceId
//        size_t surfaceId = 0;
//        size_t i = 0;
//        for (; i < mesh.F.size(); i++){
//            if (mesh.F.at(i).isBoundary){
//                if (surfaceId == cellId){
//                    const Face& face = mesh.F.at(i);
//                    glm::vec3 faceCenter(0.0, 0.0, 0.0);
//                    const size_t faceVidsSize = face.Vids.size();
//                    for (size_t j = 0; j < faceVidsSize; j++)
//                    {
//                        const Vertex& v = mesh.V.at(face.Vids.at(j));
//                        faceCenter.x += v.x;
//                        faceCenter.y += v.y;
//                        faceCenter.z += v.z;
//                    }
//                    x += faceCenter.x / faceVidsSize;
//                    y += faceCenter.y / faceVidsSize;
//                    z += faceCenter.z / faceVidsSize;
//                    id = mesh.C.size() + cellId;
//                    N_Vids.push_back(face.N_Cids.at(0));
//                    N_Eids.push_back(mesh.C.size() + surfaceId);
//                    break;
//                }
//                surfaceId++;
//            }
//        }
//    }
//}

void Frame::BuildFrom(const Mesh& mesh, const size_t cellId, bool isBoundary, size_t surfaceId)
{
    this->id = cellId;
    this->isBoundary = isBoundary;
    if (!isBoundary){
        const Cell& cell = mesh.C.at(cellId);
        const size_t facesSize = cell.N_Fids.size();
        // Compute center x,y,z coordinates
        for (size_t i = 0; i < facesSize; i++)
        {
            const Face& face = mesh.F.at(cell.N_Fids.at(i));
            glm::vec3 faceCenter(0.0, 0.0, 0.0);
            const size_t faceVidsSize = face.Vids.size();
            for (size_t j = 0; j < faceVidsSize; j++)
            {
                const Vertex& v = mesh.V.at(face.Vids.at(j));
                faceCenter.x += v.x;
                faceCenter.y += v.y;
                faceCenter.z += v.z;
            }
            x += faceCenter.x / faceVidsSize;
            y += faceCenter.y / faceVidsSize;
            z += faceCenter.z / faceVidsSize;
        }
        x /= facesSize;
        y /= facesSize;
        z /= facesSize;
        // Compute N_Eids and N_Vids
        for (size_t i = 0; i < facesSize; i++)
        {
            const Face& face = mesh.F.at(cell.N_Fids.at(i));
            if (!face.isBoundary){
                for (size_t j = 0; j < face.N_Cids.size(); j++)
                    if (face.N_Cids.at(j) != cellId){
                        N_Vids.push_back(face.N_Cids.at(j));
                        break;
                    }
            }
            else {
                size_t surfaceNodeId = mesh.C.size();
                for (size_t j = 0; j < mesh.F.size(); j++)
                    if (mesh.F.at(j).isBoundary){
                        if (face.id == mesh.F.at(j).id)
                        {
                            N_Vids.push_back(surfaceNodeId);
                            break;
                        }
                        surfaceNodeId++;
                    }
            }
            N_Eids.push_back(face.id);
        }
    }
    else{
        const Face& face = mesh.F.at(surfaceId);
        glm::vec3 faceCenter(0.0, 0.0, 0.0);
        const size_t faceVidsSize = face.Vids.size();
        for (size_t j = 0; j < faceVidsSize; j++)
        {
            const Vertex& v = mesh.V.at(face.Vids.at(j));
            faceCenter.x += v.x;
            faceCenter.y += v.y;
            faceCenter.z += v.z;
        }
        x += faceCenter.x / faceVidsSize;
        y += faceCenter.y / faceVidsSize;
        z += faceCenter.z / faceVidsSize;

        N_Vids.push_back(face.N_Cids.at(0));
        N_Eids.push_back(surfaceId);
    }
}

