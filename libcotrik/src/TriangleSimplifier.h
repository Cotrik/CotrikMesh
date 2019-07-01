/*
 * TriangleSimplifier.h
 *
 *  Created on: Jan 20, 2019
 *      Author: cotrik
 */

#ifndef LIBCOTRIK_SRC_TRIANGLE_SIMPLIFIER_H_
#define LIBCOTRIK_SRC_TRIANGLE_SIMPLIFIER_H_

#include "Simplifier.h"
class TriangleSimplifier : public Simplifier {
public:
    TriangleSimplifier(Mesh& mesh);
    virtual ~TriangleSimplifier();
private:
    TriangleSimplifier();
    TriangleSimplifier(const TriangleSimplifier&);
    TriangleSimplifier& operator = (const TriangleSimplifier&);
public:
    void Run(std::set<size_t>& canceledFids);
    void Collapse(const Vertex& v);
    int GetCollapseValence(const Vertex& v);
private:
    bool AllNeighboringVerticesRegular(const Vertex& v);
    bool AllNeighboringVerticesInterior(const Vertex& v);
    bool AllDiagonalVerticesSingular(const Vertex& v);
};

#endif /* LIBCOTRIK_SRC_TRIANGLE_SIMPLIFIER_H_ */
