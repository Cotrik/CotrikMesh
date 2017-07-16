/*
 * FrameField.h
 *
 *  Created on: Nov 10, 2016
 *      Author: cotrik
 */

#ifndef FRAMEFIELD_H_
#define FRAMEFIELD_H_

#include <vector>
#include "Frame.h"

class Mesh;
class FrameEdge;
class PolyLines;

class FrameField
{
public:
    FrameField();
    virtual ~FrameField();

public:
    void BuildFrom(const Mesh& mesh);
    void SetPolylineIds(const PolyLines& polyLines);
    void BuildOrthogonal4EidsForEachFrame();                    // Run after SetPolylineIds;
    void BuildOrthogonal4VidsForEachFrame();                    // Run after BuildOrthogonal4EidsForEachFrame;
    void WriteFile(const char* filename);
    void WriteFile(const char* filename, const std::vector<size_t>& frameIds);
    void WriteColorFile(const char* filename);

public:
    std::vector<Frame> frameNodes;
    std::vector<FrameEdge> frameEdges;
};

#endif /* FRAMEFIELD_H_ */
