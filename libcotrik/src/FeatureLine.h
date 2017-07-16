/*
 * FeatureLine.h
 *
 *  Created on: Jan 3, 2017
 *      Author: cotrik
 */

#ifndef SRC_FEATURELINE_H_
#define SRC_FEATURELINE_H_

#include "Mesh.h"

class FeatureLine
{
public:
    FeatureLine(const Mesh& mesh);
    FeatureLine(const FeatureLine& r);
    virtual ~FeatureLine();
private:
    FeatureLine();
    FeatureLine& operator = (const FeatureLine& r);
public:
    void Extract(const size_t label = 0);
    size_t FindNextVid(const size_t vid, std::vector<bool>& visited, size_t& nextEid);
public:
    std::vector<size_t> Vids;  // consecutive vertex ID set
    std::vector<size_t> Eids;  // consecutive edge   ID set
    const Mesh& mesh;
};

#endif /* SRC_FEATURELINE_H_ */
