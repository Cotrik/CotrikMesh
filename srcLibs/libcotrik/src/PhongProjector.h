/*
 * PhongProjector.h
 *
 *  Created on: May 5, 2017
 *      Author: cotrik
 */

#ifndef LIBCOTRIK_SRC_PHONGPROJECTOR_H_
#define LIBCOTRIK_SRC_PHONGPROJECTOR_H_

#include "Mesh.h"

class PhongProjector
{
public:
    PhongProjector(const Mesh& mesh);
    virtual ~PhongProjector();

public:
    void ProjectPointOnQuad(const Face& quad, const glm::dvec3& p, glm::dvec2& uv, glm::dvec3& interpolP, glm::dvec3& interpolN);
    void ProjectPointOnTriangle(const Face& tri, const glm::dvec3& p, glm::dvec2& uv, glm::dvec3& interpolP, glm::dvec3& interpolN);
private:
    const Mesh& mesh;
};

glm::dvec3 bilinear(const glm::dvec3& v1, const glm::dvec3& v2, const glm::dvec3& v3, const glm::dvec3& v4, const glm::dvec2& uv);
glm::dvec3 barycentric(const glm::dvec3& v1, const glm::dvec3& v2, const glm::dvec3& v3, const glm::dvec2& uv);

#endif /* LIBCOTRIK_SRC_PHONGPROJECTOR_H_ */
