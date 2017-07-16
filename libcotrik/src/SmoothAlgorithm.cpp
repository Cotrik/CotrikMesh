/*
 * SmoothAlgorithm.cpp
 *
 *  Created on: Mar 29, 2017
 *      Author: cotrik
 */

#include "SmoothAlgorithm.h"
#include "MeshQuality.h"
#include <iostream>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Cholesky>
#include <eigen3/Eigen/LU>
#include <eigen3/Eigen/SVD>
#include <eigen3/Eigen/SparseCore>
#include <eigen3/Eigen/SparseLU>
#include <eigen3/Eigen/SparseQR>
#include <eigen3/Eigen/SparseCholesky>

SmoothAlgorithm::SmoothAlgorithm(Mesh& mesh, const Smooth_Algorithm smoothAlgorithm/* = LAPLACIAN*/)
: mesh(mesh)
, m_smoothAlgorithm(smoothAlgorithm)
{
    // TODO Auto-generated constructor stub

}

SmoothAlgorithm::~SmoothAlgorithm()
{
    // TODO Auto-generated destructor stub
}

int SmoothAlgorithm::Run(const size_t iters/* = 1e2*/, const double eps/* = 1e-4*/)
{
    int error_code = 0;
    bool converged = false;
    double last_energy = 1e10;
    int iter = 0;
    std::cout << "===== Smooth max iters = " << iters << " eps = " << eps << "\n";
    while (!converged && iter++ != iters) {
        double current = SmoothVolume();
        converged = fabs(current - last_energy)/last_energy < eps;
        last_energy = current;
        std::cout << "Smooth iter = " << iter << "\n";
    }
    return error_code;
}

glm::vec3 SmoothAlgorithm::LapLace(const Vertex& v)
{
    {
        glm::vec3 sum(0.0, 0.0, 0.0);
        int count = 0;
        for (size_t j = 0; j < v.N_Vids.size(); j++) {
            const Vertex& n_v = mesh.V.at(v.N_Vids[j]);
            sum += n_v.xyz();
            count++;
        }
        return glm::vec3(sum.x/count, sum.y/count, sum.z/count);
    }
}
const int SJP[8][3] = {
    { 1, 3, 4 },
    { 2, 0, 5 },
    { 3, 1, 6 },
    { 0, 2, 7 },
    { 7, 5, 0 },
    { 4, 6, 1 },
    { 5, 7, 2 },
    { 6, 4, 3 }
};

glm::vec3 SmoothAlgorithm::ScaledJacobian(const Vertex& v)
{
    {
        const std::vector<Vertex>& V = mesh.V;
        const std::vector<Cell>& C = mesh.C;
//        glm::vec3 b;
//        glm::mat3x3 m;
        Eigen::MatrixXd m(3, 3);
        Eigen::VectorXd b(3);
        glm::vec3 sum(0.0, 0.0, 0.0);
        glm::vec3 x;
        int count = 0;
        for (size_t i = 0; i < v.N_Cids.size(); i++) {
            const Cell& c = C[v.N_Cids.at(i)];
            int vid = 0;
            for (int j = 0; j < c.Vids.size(); j++)
                if (c.Vids.at(j) == v.id) {
                    vid = j;
                    break;
                }
            const glm::vec3& v1 = V.at(c.Vids.at(SJP[vid][0]));
            const glm::vec3& v2 = V.at(c.Vids.at(SJP[vid][1]));
            const glm::vec3& v3 = V.at(c.Vids.at(SJP[vid][2]));
//            m[0][0] = v2.x - v3.x; m[0][1] = v2.y - v3.y; m[0][2] = v2.z - v3.z;
//            m[1][0] = v1.x - v3.x; m[1][1] = v1.y - v3.y; m[1][2] = v1.z - v3.z;
//            m[2][0] = v1.x - v2.x; m[2][1] = v1.y - v2.y; m[2][2] = v1.z - v2.z;
//            b.x = v1.x * (v2.x - v3.x) + v1.y * (v2.y - v3.y) + v1.z * (v2.z - v3.z);
//            b.y = v2.x * (v1.x - v3.x) + v2.y * (v1.y - v3.y) + v2.z * (v1.z - v3.z);
//            b.z = v3.x * (v1.x - v2.x) + v3.y * (v1.y - v2.y) + v3.z * (v1.z - v2.z);
//            const glm::vec3 x = glm::inverse(m)*b;

            m(0,0) = v2.x - v3.x; m(0,1) = v2.y - v3.y; m(0,2) = v2.z - v3.z;
            m(1,0) = v1.x - v3.x; m(1,1) = v1.y - v3.y; m(1,2) = v1.z - v3.z;
            m(2,0) = v1.x - v2.x; m(2,1) = v1.y - v2.y; m(2,2) = v1.z - v2.z;
            b(0) = v1.x * (v2.x - v3.x) + v1.y * (v2.y - v3.y) + v1.z * (v2.z - v3.z);
            b(1) = v2.x * (v1.x - v3.x) + v2.y * (v1.y - v3.y) + v2.z * (v1.z - v3.z);
            b(2) = v3.x * (v1.x - v2.x) + v3.y * (v1.y - v2.y) + v3.z * (v1.z - v2.z);
            double b0 = b(0);
            double b1 = b(1);
            double b2 = b(2);

            Eigen::VectorXd X = m.llt().solve(b);
            x.x = X(0);
            x.y = X(1);
            x.z = X(2);
            sum += x;
            count++;
        }
        return glm::vec3(sum.x/count, sum.y/count, sum.z/count);
    }
}

double SmoothAlgorithm::SmoothVolume(const Smooth_Algorithm smoothMethod/* = LAPLACE_EDGE*/)
{
    std::vector<Vertex>& V = mesh.V;
    std::vector<glm::vec3> oldV(mesh.V.size());
    std::vector<glm::vec3> newV(mesh.V.size());
    for (size_t i = 0; i < V.size(); i++) {
        oldV[i] = V[i];
    }
    //while (iters-- != 0)
    {
//#pragma omp parallel for
        for (size_t i = 0; i < V.size(); i++) {
            const Vertex& v = V.at(i);
            if (v.isBoundary)
                continue;
            if (m_smoothAlgorithm == LAPLACIAN) {
                newV.at(i) = LapLace(v);
            }else if (m_smoothAlgorithm == SCALED_JACOBIAN)
                newV.at(i) = ScaledJacobian(v);
        }

        for (size_t i = 0; i < V.size(); i++) {
            Vertex& v = V.at(i);
            if (v.isBoundary)
                continue;
            v = newV[i];
        }
    }

    double energy = 0;
    for (size_t i = 0; i < V.size(); i++) {
        const Vertex& v = V.at(i);
        if (v.isBoundary)
            continue;
        const glm::vec3& oldv = oldV.at(i);
        const double distance = glm::length(v.xyz() - oldv);
        energy += distance * distance;
    }

    std::cout << "Volume Energy = " << energy << std::endl;
    return energy;
}
