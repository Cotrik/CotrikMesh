/*
 * Frame.h
 *
 *  Created on: Nov 10, 2016
 *      Author: cotrik
 */

#ifndef FRAME_H_
#define FRAME_H_

#include "Mesh.h"

typedef glm::dvec3 FrameNode;
class FrameEdge : public GeoInfo
{
public:
    FrameEdge()
    : polylineId(MAXID)
    , length(0.0)
    {}

    FrameEdge(size_t vnum)
    : polylineId(MAXID)
    , length(0.0)
    {
        Vids.resize(vnum);
    }
    FrameEdge(const std::vector<size_t>& Vids)
    : Vids(Vids)
    , polylineId(MAXID)
    , length(0.0)
    {
    }
    virtual ~FrameEdge()
    {}
    void SetPolylineId(const size_t value)
    {
        polylineId = value;
    }

public:
    std::vector<size_t> Vids;    // ids of two cells
    size_t polylineId;           // indicates to which polyline it belong;
    double length;
};

class Frame: public FrameNode, public GeoInfo
{
public:
    Frame();
//    Frame(size_t vnum);
//    Frame(size_t vnum, size_t eNum, size_t fNum);
//    Frame(const std::vector<size_t>& Vids);
//    Frame(const std::vector<size_t>& Vids, const std::vector<size_t> Eids, const std::vector<size_t> Fids);
    virtual ~Frame();

public:
    void BuildFrom(const Mesh& mesh, const size_t cellId, bool isBoundary = false, size_t surfaceId = MAXID);

public:
    std::vector<size_t> N_Vids;          // Neighboring Cids(Cell ids);
    std::vector<size_t> N_Eids;          // Neighboring Fids(Face ids);
    std::vector<std::vector<size_t> > ortho4Vids; // 4Vids that are orthogonal to one of N_Vids
    std::vector<std::vector<size_t> > ortho4Eids; // 4Eids that are orthogonal to one of N_Eids
};


#endif /* FRAME_H_ */
