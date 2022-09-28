/*
 * PhongProjector.cpp
 *
 *  Created on: May 5, 2017
 *      Author: cotrik
 */

#include "PhongProjector.h"

#include <Eigen/Core>
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Cholesky>
#include <Eigen/LU>
#include <Eigen/SVD>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#include <Eigen/SparseCholesky>

#define GLM_FORCE_RADIANS
#include <glm/mat3x2.hpp>
#include <glm/gtc/matrix_transform.hpp>

PhongProjector::PhongProjector(const Mesh& mesh)
: mesh(mesh)
{
    // TODO Auto-generated constructor stub

}

PhongProjector::~PhongProjector()
{
    // TODO Auto-generated destructor stub
}

void PhongProjector::ProjectPointOnQuad(const Face& quad, const glm::dvec3& p, glm::dvec2& uv, glm::dvec3& interpolP, glm::dvec3& interpolN)
{
    const Vertex& v1 = mesh.V[quad.Vids[0]];
    const Vertex& v2 = mesh.V[quad.Vids[1]];
    const Vertex& v3 = mesh.V[quad.Vids[2]];
    const Vertex& v4 = mesh.V[quad.Vids[3]];

    uv.x = 0.5;
    uv.y = 0.5;

    //Newton iteration for projection on quad
    //objective function:
    // F(u, v) = (p - interpolP) x interpolN = 0

    for (int i = 0; i < 4; ++i) {
        interpolP = bilinear(v1.xyz(), v2.xyz(), v3.xyz(), v4.xyz(), uv);
        interpolN = bilinear(v1.normal, v2.normal, v3.normal, v4.normal, uv);

        glm::dvec3 dPdu = (1 - uv.y) * v2.xyz() + uv.y * v3.xyz() - ((1 - uv.y) * v1.xyz() + uv.y * v4.xyz());
        glm::dvec3 dPdv = (1 - uv.x) * v4.xyz() + uv.x * v3.xyz() - ((1 - uv.x) * v1.xyz() + uv.x * v2.xyz());
        glm::dvec3 dNdu = (1 - uv.y) * v2.normal + uv.y * v3.normal - ((1 - uv.y) * v1.normal + uv.y * v4.normal);
        glm::dvec3 dNdv = (1 - uv.x) * v4.normal + uv.x * v3.normal - ((1 - uv.x) * v1.normal + uv.x * v2.normal);

        glm::dvec3 f = glm::cross(p - interpolP, interpolN);
        glm::dvec3 dFdu = glm::cross(-dPdu, interpolN) + glm::cross(p - interpolP, dNdu);
        glm::dvec3 dFdv = glm::cross(-dPdv,interpolN) + glm::cross(p - interpolP, dNdv);

        Eigen::Matrix<double, 3, 2> jacobian;
        jacobian(0,0) = dFdu.x; jacobian(1,0) = dFdu.y; jacobian(2,0) = dFdu.z;
        jacobian(0,1) = dFdv.x; jacobian(1,1) = dFdv.y; jacobian(2,1) = dFdv.z;

        //std::cout << uv.transpose() << " => " << F.transpose() << std::endl;
        Eigen::Vector3d F(f.x, f.y, f.z);
        Eigen::Vector2d rhs = -jacobian.transpose() * F;
        auto lhs = jacobian.transpose() * jacobian;
        double norm = 1.0f / (lhs(0, 0) * lhs(1, 1) - lhs(0, 1) * lhs(1, 0));

        uv += glm::dvec2(lhs(1, 1) * rhs.x() - lhs(0, 1) * rhs.y(), -lhs(1, 0) * rhs.x() + lhs(0, 0) * rhs.y()) * norm;
    }

    interpolP = bilinear(v1.xyz(), v2.xyz(), v3.xyz(), v4.xyz(), uv);
    interpolN = bilinear(v1.normal, v2.normal, v3.normal, v4.normal, uv);
}

void PhongProjector::ProjectPointOnTriangle(const Face& tri, const glm::dvec3& p, glm::dvec2& uv, glm::dvec3& interpolP, glm::dvec3& interpolN)
{
    const Vertex& v1 = mesh.V[tri.Vids[0]];
    const Vertex& v2 = mesh.V[tri.Vids[1]];
    const Vertex& v3 = mesh.V[tri.Vids[2]];

    uv.x = 0.333f;
    uv.y = 0.333f;

    //Newton iteration for projection on triangle

    //objective function:
    // F(u, v) = (p - interpolP) x interpolN = 0
    Eigen::Vector3d F;

    glm::dvec3 dPdu = v1.xyz() - v3.xyz();
    glm::dvec3 dPdv = v2.xyz() - v3.xyz();
    glm::dvec3 dNdu = v1.normal - v3.normal;
    glm::dvec3 dNdv = v2.normal - v3.normal;

    for (int i = 0; i < 4; ++i) {
        interpolP = barycentric(v1.xyz(), v2.xyz(), v3.xyz(), uv);
        interpolN = barycentric(v1.normal, v2.normal, v3.normal, uv);

        glm::dvec3 f = glm::cross(p - interpolP, interpolN);
        glm::dvec3 dFdu = glm::cross(-dPdu, interpolN) + glm::cross(p - interpolP, dNdu);
        glm::dvec3 dFdv = glm::cross(-dPdv,interpolN) + glm::cross(p - interpolP, dNdv);

        Eigen::Matrix<double, 3, 2> jacobian;
        jacobian(0,0) = dFdu.x; jacobian(1,0) = dFdu.y; jacobian(2,0) = dFdu.z;
        jacobian(0,1) = dFdv.x; jacobian(1,1) = dFdv.y; jacobian(2,1) = dFdv.z;

        //std::cout << uv.transpose() << " => " << F.transpose() << std::endl;
        Eigen::Vector3d F(f.x, f.y, f.z);
        Eigen::Vector2d rhs = -jacobian.transpose() * F;
        auto lhs = jacobian.transpose() * jacobian;
        double norm = 1.0f / (lhs(0, 0) * lhs(1, 1) - lhs(0, 1) * lhs(1, 0));

        uv += glm::dvec2(lhs(1, 1) * rhs.x() - lhs(0, 1) * rhs.y(), -lhs(1, 0) * rhs.x() + lhs(0, 0) * rhs.y()) * norm;
    }

    interpolP = barycentric(v1.xyz(), v2.xyz(), v3.xyz(), uv);
    interpolN = barycentric(v1.normal, v2.normal, v3.normal, uv);
}

glm::dvec3 bilinear(const glm::dvec3& v1, const glm::dvec3& v2, const glm::dvec3& v3, const glm::dvec3& v4, const glm::dvec2& uv)
{
    return (1 - uv.x) * ((1 - uv.y) * v1 + uv.y * v4) + uv.x * ((1 - uv.y) * v2 + uv.y * v3);
}
glm::dvec3 barycentric(const glm::dvec3& v1, const glm::dvec3& v2, const glm::dvec3& v3, const glm::dvec2& uv)
{
    return uv.x * v1 + uv.y * v2 + (1 - uv.x - uv.y) * v3;
}
