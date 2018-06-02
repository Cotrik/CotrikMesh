/*
 * RefineDualQuad.h
 *
 *  Created on: Jan 12, 2018
 *      Author: cotrik
 */

#ifndef LIBCOTRIK_SRC_REFINEDDUALQUAD_H_
#define LIBCOTRIK_SRC_REFINEDDUALQUAD_H_

#include "Dual.h"

class RefinedDualQuad {
public:
    RefinedDualQuad(Mesh& mesh);
    virtual ~RefinedDualQuad();
    void WriteFaces() const;
    void WriteEdges() const;
private:
    RefinedDualQuad(Dual&);
    RefinedDualQuad();
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

Mesh GetRefineQuadMesh(const Mesh& quad_mesh, int clockwise = 0);
#endif /* LIBCOTRIK_SRC_REFINEDDUALQUAD_H_ */
