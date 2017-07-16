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
    void ProjectPointOnQuad(const Face& quad, const glm::vec3& p, glm::vec2& uv, glm::vec3& interpolP, glm::vec3& interpolN);
    void ProjectPointOnTriangle(const Face& tri, const glm::vec3& p, glm::vec2& uv, glm::vec3& interpolP, glm::vec3& interpolN);
private:
    const Mesh& mesh;
};

glm::vec3 bilinear(const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& v3, const glm::vec3& v4, const glm::vec2& uv);
glm::vec3 barycentric(const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& v3, const glm::vec2& uv);

#endif /* LIBCOTRIK_SRC_PHONGPROJECTOR_H_ */
