/*
 * EdgeLines.h
 *
 *  Created on: Nov 25, 2016
 *      Author: cotrik
 */

#ifndef EDGE_LINES_H
#define EDGE_LINES_H

#include <vector>
#include "Mesh.h"

class EdgeLine
{
public:
    EdgeLine();
    virtual ~EdgeLine();

public:
    // startV is a Vertex on the boundary; startEdgeId tell the direction, return endEdgeId
    size_t BuildFrom(const Mesh& mesh, size_t startFrameNodeId, size_t startFrameEdgeId);
    size_t BuildFrom(const Mesh& mesh, size_t startFrameNodeId, size_t startFrameEdgeId, size_t startFaceId);

public:
    size_t id;
    std::vector<size_t> Vids; // Vertex ids
    std::vector<size_t> Eids; // Edge ids
};

class EdgeLines
{
public:
    EdgeLines(const Mesh& mesh);
    virtual ~EdgeLines();
private:
    EdgeLines();
    EdgeLines(const EdgeLines& edgeLines);
public:
    void Build();
    void WriteFile(const char* filename);

public:
    std::vector<EdgeLine> edgeLines;

private:
    const Mesh& mesh;
};

#endif /* EDGE_LINES_H */
