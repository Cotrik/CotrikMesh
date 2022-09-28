/*
 * PolyLine.h
 *
 *  Created on: Nov 11, 2016
 *      Author: cotrik
 */

#ifndef POLY_LINE_H_
#define POLY_LINE_H_

#include <vector>
#include "FrameField.h"

class PolyLine
{
public:
    PolyLine();
    virtual ~PolyLine();

public:
    // startFrameNode is a Node on the boundary; startFrameEdgeId tell the direction, return endFrameEdgeId
    size_t BuildFrom(const Mesh& mesh, const FrameField& framefield, size_t startFrameNodeId, size_t startFrameEdgeId);

public:
    size_t id;
    std::vector<size_t> Vids; // Vertex ids
    std::vector<size_t> Eids; // Edge ids
};

class PolyLines
{
public:
    PolyLines(const Mesh& mesh, const FrameField& framefield);
    virtual ~PolyLines();
private:
    PolyLines();
    PolyLines(const PolyLines& polyLines);
public:
    void BuildFrom(const Mesh& mesh);
    void BuildFrom(const FrameField& framefield);
    void BuildFrom(const Mesh& mesh, const FrameField& framefield);
    void Build();
    void WriteFile(const char* filename);

public:
    std::vector<PolyLine> polyLines;

private:
    const Mesh& mesh;
    const FrameField& framefield;
};

#endif /* POLY_LINE_H_ */
