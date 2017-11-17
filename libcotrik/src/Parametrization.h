/*
 * Parametrization.h
 *
 *  Created on: Nov 10, 2017
 *      Author: cotrik
 */

#ifndef LIBCOTRIK_SRC_PARAMETRIZATION_H_
#define LIBCOTRIK_SRC_PARAMETRIZATION_H_

#include "BaseComplex.h"

class Parametrization {
public:
    Parametrization(BaseComplex& baseComplex);
    virtual ~Parametrization();
public:
    void Build();
    void WriteMeshFile(const char* filename);
private:
    Parametrization();
    Parametrization(const Parametrization&);
private:
    void BuildV();
    void BuildE();
    void BuildF();
    void BuildC();
    // Tao Ju. <<Mean Value Coordinates for Closed Triangular Meshes>>
    void Parametrize3DInnerPoint(const Vertex& innerPoint,
            const std::vector<Vertex>& boundaryVertices,
            const std::vector<Cell>& boundaryTriangles,
            std::vector<glm::vec3>& w,
            const double eps = 1e-8);
    BaseComplex& baseComplex;
    Mesh newMesh;
};

#endif /* LIBCOTRIK_SRC_PARAMETRIZATION_H_ */
