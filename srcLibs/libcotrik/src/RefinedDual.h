/*
 * RefineDual.h
 *
 *  Created on: Nov 21, 2017
 *      Author: cotrik
 */

#ifndef LIBCOTRIK_SRC_REFINEDDUAL_H_
#define LIBCOTRIK_SRC_REFINEDDUAL_H_

#include "Dual.h"

class RefinedDual {
public:
    RefinedDual(Mesh& mesh);
    virtual ~RefinedDual();
    void WriteFaces() const;
    void WriteEdges() const;
private:
    RefinedDual(Dual&);
    RefinedDual();
public:
    void Build();
private:
    void BuildV();
    void BuildE();
    void BuildF();
    void BuildC();
public:
    std::vector<DualVertex> V;
    std::vector<DualEdge> E;
    std::vector<DualFace> F;
    std::vector<DualCell> C;
private:
    Mesh& mesh;
    //Mesh refinedMesh;
};

Mesh GetRefineMesh2(const Mesh& hex_mesh, int clockwise = 0);
Mesh GetRefineMesh3(const Mesh& hex_mesh, int clockwise = 0);
#endif /* LIBCOTRIK_SRC_REFINEDDUAL_H_ */
