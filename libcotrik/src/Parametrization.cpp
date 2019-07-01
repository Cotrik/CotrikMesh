/*
 * Parametrization.cpp
 *
 *  Created on: Nov 10, 2017
 *      Author: cotrik
 */

#include "Parametrization.h"

Parametrization::Parametrization(BaseComplex& basecomplex): baseComplex(basecomplex) {
    // TODO Auto-generated constructor stub

}

Parametrization::~Parametrization() {
    // TODO Auto-generated destructor stub
}

void Parametrization::Build() {
    BuildV();
    BuildE();
    BuildF();
    BuildC();
}


void Parametrization::BuildV() {
//    size_t numOfComponentV = baseComplex.componentV.size();
//    size_t numOfComponentE = baseComplex.componentE.size();
//    size_t numOfComponentF = baseComplex.componentF.size();
//    size_t numOfComponentC = baseComplex.componentC.size();
//    size_t totalV = numOfComponentV + 9 * numOfComponentE + 9 * 9 * numOfComponentF + 9 * 9 * 9 * numOfComponentC;
//    newMesh.V.resize(totalV);
//    size_t id = 0;
//    for (const auto& v : baseComplex.componentV) {
//        newMesh.V.at(id) = v.xyz();
//        newMesh.V.at(id).id = id++;
//    }
//
//    for (const auto& e : baseComplex.componentE) {
//        const auto& v1 = baseComplex.componentV.at(e.eids_link.front());
//        const auto& v2 = baseComplex.componentV.at(e.eids_link.back());
//        for (float t = 0.1f; t < 0.95; t += 0.1f) {
//            const auto v = v1 + t * (v2.xyz() - v1.xyz());
//            newMesh.V.at(id) = v;
//            newMesh.V.at(id).id = id++;
//        }
//    }
//
//    for (const auto& f : baseComplex.componentF) {
//        const auto& v1 = baseComplex.componentF.at(e.eids_link.front());
//        const auto& v2 = baseComplex.componentV.at(e.eids_link.back());
//        for (float t = 0.1f; t < 0.95; t += 0.1f) {
//            const auto v = v1 + t * (v2.xyz() - v1.xyz());
//            newMesh.V.at(id) = v;
//            newMesh.V.at(id).id = id++;
//        }
//    }

    size_t numOfComponentV = baseComplex.componentV.size();
    size_t numOfComponentE = baseComplex.componentE.size();
    size_t numOfComponentF = baseComplex.componentF.size();
    size_t numOfComponentC = baseComplex.componentC.size();
    size_t totalV = 9 * 9 * 9 * numOfComponentC;
    newMesh.V.resize(totalV);
    size_t id = 0;
    /*
                      3__________________2
                      /|                 /|
                     / |                / |
                    /  |               /  |
                0  /___|_____________1/   |
                   |   |              |   |
                   |   |              |   |
                   |   |              |   |
                   |   |______________|___|
                   |   / 7            |  /6
                   |  /               | /
                   | /                |/
                   |/_________________/
                 4                    5
    */

    std::vector<Vertex> cubeV(8);
    cubeV[0].x = 0.0f, cubeV[0].y = 1.0f, cubeV[0].z = 1.0f;
    cubeV[1].x = 1.0f, cubeV[1].y = 1.0f, cubeV[1].z = 1.0f;
    cubeV[2].x = 1.0f, cubeV[2].y = 1.0f, cubeV[2].z = 0.0f;
    cubeV[3].x = 0.0f, cubeV[3].y = 1.0f, cubeV[3].z = 0.0f;
    cubeV[4].x = 0.0f, cubeV[4].y = 0.0f, cubeV[4].z = 1.0f;
    cubeV[5].x = 1.0f, cubeV[5].y = 0.0f, cubeV[5].z = 1.0f;
    cubeV[6].x = 1.0f, cubeV[6].y = 0.0f, cubeV[6].z = 0.0f;
    cubeV[7].x = 0.0f, cubeV[7].y = 0.0f, cubeV[7].z = 0.0f;
    std::vector<Vertex> innerV(9 * 9 * 9);
    for (float x = 0.1f; x < 0.95f; x += 0.1f)
        for (float y = 0.1f; y < 0.95f; y += 0.1f)
            for (float z = 0.1f; z < 0.95f; z += 0.1f)
                innerV[id++] = glm::dvec3(x, y, z);
    id = 0;
    std::vector<std::vector<glm::dvec3>> w(9 * 9 * 9);
    for (auto i = 0; i < innerV.size(); ++i)
        Parametrize3DInnerPoint(innerV.at(i), cubeV, innerV.C, w, 1e-6);

    for (const auto& c : baseComplex.componentC) {
        std::vector<Vertex> deformedTetV();
        for (int i = 0; i < origTetMesh.V.size(); i++)
        {
            std::vector<glm::dvec3> w;
            parametrizer.Parametrize3DInnerPoint(origTetMesh.V.at(i), origTriMesh.V, origTriMesh.C, w, 1e-6);
            double total_w = 0;
            glm::dvec3 total_p(0.0, 0.0, 0.0);
            for (int j = 0; j < origTriMesh.C.size(); j++)
            {
                const Cell& c = origTriMesh.C.at(j);
                const double wi[3] = {w[j].x, w[j].y, w[j].z};
                for (int k = 0; k < 3; k++)
                {
                    total_w += wi[k];
                    total_p.x += wi[k] * deformedMesh.V.at(c.at(k)).x;
                    total_p.y += wi[k] * deformedMesh.V.at(c.at(k)).y;
                    total_p.z += wi[k] * deformedMesh.V.at(c.at(k)).z;
                }
            }
            const Vertex v(total_p.x/total_w, total_p.y/total_w, total_p.z/total_w);
            deformedTetV.at(i).x = v.x;
            deformedTetV.at(i).y = v.y;
            deformedTetV.at(i).z = v.z;
        }
    }
}

void MapBoundaryPointsToUnitSphere(const Vertex& center, const std::vector<Vertex>& boundaryVertices, std::vector<Vertex>& V) {
    const std::vector<Vertex>& VB = boundaryVertices;
    for (size_t i = 0; i < VB.size(); i++) {
        const Vertex& v = VB.at(i);
        const glm::dvec3 d(v.x - center.x, v.y - center.y, v.z - center.z);
        const double l = glm::length(d);
        const double r = 1.0 / l;
        const glm::dvec3 newV(center.x + r * d.x, center.y + r * d.y, center.z + r * d.z);
        V.at(i).x = newV.x;
        V.at(i).y = newV.y;
        V.at(i).z = newV.z;
    }
}

int GetBoundaryVertexIndexThatIsCloseToInnerPoint(const Vertex& innerPoint, const std::vector<Vertex>& boundaryVertices, const double eps = 1e-8) {
    const std::vector<Vertex>& VB = boundaryVertices;
    for (size_t i = 0; i < VB.size(); i++) {
        const Vertex& v = VB.at(i);
        const glm::dvec3 d(v.x - innerPoint.x, v.y - innerPoint.y, v.z - innerPoint.z);
        const double l = glm::length(d);
        if (l < eps) {
//          std::cout << "Boundary Vertex index <" << i << ">is Very Close To InnerPoint" << std::endl;
            return i;
        }
    }

    return -1;
}
const double PI = 3.1415926536;
void GetWeightsForAlltriangles(const Vertex& center, const std::vector<Vertex>& VB, const std::vector<Vertex>& V, const std::vector<Cell>& C, std::vector<glm::dvec3>& W,
        const double eps = 1e-8) {
    glm::dvec3 x(center.x, center.y, center.z);
    glm::dvec3 cc = x;
    for (size_t j = 0; j < C.size(); j++) {
        const Cell& cell = C.at(j);

        glm::dvec3 u[3], p[3];
        double d[3];

        for (size_t i = 0; i < 3; i++) {
            const Vertex& vbi = VB.at(cell.at(i));
            p[i].x = vbi.x;
            p[i].y = vbi.y;
            p[i].z = vbi.z;
            d[i] = glm::length(p[i] - x);
            u[i] = p[i] - x;
            u[i].x /= d[i];
            u[i].y /= d[i];
            u[i].z /= d[i];
        }

        double theta[3];
        for (size_t i = 0; i < 3; i++) {
            const double l_i = glm::length(u[(i + 1) % 3] - u[(3 + i - 1) % 3]);
            theta[i] = 2.0 * asin(0.5 * l_i);
        }

        const double h = (theta[0] + theta[1] + theta[2]) * 0.5;
        if (fabs(PI - h) < eps) {
            // x lies on t, use 2D barycentric coordinates
            double w[3] = { 0, 0, 0 };
            for (size_t i = 0; i < 3; i++) {
                w[i] = sin(theta[i]) * d[(3 + i - 1) % 3] * d[(i + 1) % 3];
            }
            const glm::dvec3 wi(w[0], w[1], w[2]);
            W.clear();
            W.resize(C.size(), glm::dvec3(0, 0, 0));
            W.at(j) = wi;
//          std::cout << "point lies on t, use 2D barycentric coordinates" << std::endl;
//          std::cout << "(cell_index = " << j << " w1 = " << w[0] << " w2 = " << w[1] << " w3 = " << w[2] << ")\n";
            return;
        }

        const glm::mat3x3 m(u[0], u[1], u[2]);
        const double det = glm::determinant(m);
        double sign_det = 1;
        if (det < 0) sign_det = -1;
        double c[3], s[3];
        for (size_t i = 0; i < 3; i++) {
            const double num = 2.0 * sin(h) * sin(h - theta[i]);
            const double denum = sin(theta[(i + 1) % 3]) * sin(theta[(3 + i - 1) % 3]);
            c[i] = num / denum - 1.0;
            if (c[i] > 1.0) c[i] = 1.0;
            if (c[i] < -1.0) c[i] = -1.0;
            //c[i] = (2.0 * sin(h) * sin(h - theta[i]))/(sin(theta[(i + 4)%3]) * sin(theta[(i+2)%3])) - 1.0;
            s[i] = sign_det * sqrt(1.0 - c[i] * c[i]);
        }
        bool isOutside = false;
        for (size_t i = 0; i < 3; i++) {
            if (fabs(s[i]) < eps) isOutside = true;
        }
        if (isOutside) continue;
        double w[3] = { 0, 0, 0 };
        for (size_t i = 0; i < 3; i++) {
            const double num = (theta[i] - c[(i + 1) % 3] * theta[(3 + i - 1) % 3] - c[(3 + i - 1) % 3] * theta[(i + 1) % 3]);
            const double denum = (d[i] * sin(theta[(i + 1) % 3]) * s[(3 + i - 1) % 3]);
            w[i] = num / denum;
//          if (i == 2)
            //std::cout << "num = " << num << " denum = " << denum << "\n";
            //w[i] = (theta[i] - c[(i+1)%3]*theta[(3+i-1)%3] - c[(3+i-1)%3]*theta[(i+1)%3])/(di*sin(theta[i+1])*s[(3+i-1)%3]);
        }
        glm::dvec3 wi(w[0], w[1], w[2]);
//      if (wi.x < 0 || wi.y < 0 || wi.z < 0)
//      {
//          wi.x = 0;
//          wi.y = 0;
//          wi.z = 0;
//      }
        //std::cout << "(--cell_index = " << j << " w1 = " << w[0] << " w2 = " << w[1] << " w3 = " << w[2] << ")\n";
        W.at(j) = wi;
    }
}

void Parametrization::Parametrize3DInnerPoint(const Vertex& innerPoint, const std::vector<Vertex>& boundaryVertices, const std::vector<Cell>& boundaryTriangles,
        std::vector<glm::dvec3>& w, const double eps) {
    w.clear();
    w.resize(boundaryTriangles.size(), glm::dvec3(0, 0, 0));
    std::vector < Vertex > V = boundaryVertices;
    const std::vector<Cell>& C = boundaryTriangles;
    int index = GetBoundaryVertexIndexThatIsCloseToInnerPoint(innerPoint, boundaryVertices, eps);
    if (index != -1) {
        double wt[3] = { 0, 0, 0 };
        size_t j;
        for (j = 0; j < C.size(); j++) {
            bool foundIndex = false;
            size_t i;
            for (i = 0; i < 3; i++) {
                if (C.at(j).at(i) == index) {
                    wt[i] = 1.0;
                    foundIndex = true;
                    break;
                }
            }
            if (foundIndex) {
                break;
            }
        }
        const glm::dvec3 wi(wt[0], wt[1], wt[2]);
        w.at(j) = wi;
//      std::cout << "cell_index = " << j << " w1 = " << w[j].x << " w2 = " << w[j].y << " w3 = " << w[j].z << ")\n";
        return;
    }
    // Step 1, map boundary to unit sphere, the center of sphere is at inner
    MapBoundaryPointsToUnitSphere(innerPoint, boundaryVertices, V);

    // Step 2. Get wi for each triangle;
    GetWeightsForAlltriangles(innerPoint, boundaryVertices, V, C, w, eps);
}

void Parametrization::BuildE() {

}

void Parametrization::BuildF() {

}

void Parametrization::BuildC() {

}
