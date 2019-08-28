/*
 * SurfaceMeshOpt.cpp
 *
 *  Created on: Apr 8, 2017
 *      Author: cotrik
 */

#include "SurfaceMeshOpt.h"
#include "MeshQuality.h"
#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <queue>
#include <iostream>
#include <math.h>

SurfaceMeshOpt::SurfaceMeshOpt(Mesh& mesh)
: MeshOpt(mesh)
{
    // TODO Auto-generated constructor stub

}

SurfaceMeshOpt::~SurfaceMeshOpt()
{
    // TODO Auto-generated destructor stub
}


size_t SurfaceMeshOpt::Run(const size_t iters/* = 1*/)
{
    ComputeAvgEdgeLength();
    //----------------------------------------------
//    double minimumScaledJacobian = 0.0;
//    m_numOfInvertdElements = GetMinScaledJacobianVerdict(mesh, minimumScaledJacobian, this->minScaledJacobian);
//    std::cout << "iter = " << 0 << " #inverted = " << m_numOfInvertdElements << " MSJ = " << minimumScaledJacobian << std::endl;
    //----------------------------------------------
    int iter = 0;
    bool converged = false;
    double initStepSize = stepSize;
    while (!converged && iter++ < iters) {
        mesh.ClearLabelOfSurface();
        mesh.LabelSurface();
        converged = Optimize();
        stepSize *= initStepSize;
        std::string filename = std::string("MeshSurfaceOpt.") + std::to_string(iter) + ".vtk";
        //std::cout << "------ iter " << "\033[1;32m" << iter << "\033[0m" << std::endl;
        std::cout << "------ iter " <<  iter  << std::endl;
//        MeshFileWriter optwriter(mesh, filename.c_str());
//        optwriter.WriteFile();
    }

    for (size_t i = 0; i < mesh.V.size(); i++) {
        m_refMesh->V[mesh.m_refIds[i]] = mesh.V[i].xyz();
    }
    MeshFileWriter optwriter(*m_refMesh, "SurfaceMeshOpt.vtk");
    optwriter.WriteFile();
    return iter - 1;
}

bool SurfaceMeshOpt::Optimize()
{
    std::vector<Trip> A_Entries;
    std::vector<float> b;
    size_t row = 0;
    OptimizeSingularity(A_Entries, b, row);
    OptimizeEdgeOrthogonality(A_Entries, b, row);
    OptimizeEdgeStraightness(A_Entries, b, row);
    size_t a_id = OptimizeSurfaceVertices(A_Entries, b, row);

    const size_t col = 3 * mesh.V.size() + a_id;

    Eigen::SparseMatrix<float> A(row, col);
    A.setFromTriplets(A_Entries.begin(), A_Entries.end());
    Eigen::SparseMatrix<float> AT = A.transpose();
    Eigen::SparseMatrix<float> ATA = AT * A;

    //  MatrixXf ;
    VectorXf B(row);
    for (size_t i = 0; i < row; i++)
        B(i) = b[i];

    VectorXf ATB = AT * B;
    Eigen::SimplicialLDLT<SpMat> chol(ATA);
    //Eigen::SparseLU<SpMat> chol(ATA);
    VectorXf X(col);
    X = chol.solve(ATB);

    bool converged = true;
    if (recoverable) {
        //Mesh targetMesh(mesh);
//        std::vector<glm::dvec3> oldV(mesh.V.size());
//        for (size_t i = 0; i < mesh.V.size(); i++) {
//            oldV[i].x = mesh.V[i].x;            oldV[i].y = mesh.V[i].y;            oldV[i].z = mesh.V[i].z;
//        }

        const double maxEdgeLength = 2.0 * avgMeshEdgeLength;
        for (size_t i = 0; i < mesh.V.size(); i++) {
            Vertex& v = mesh.V.at(i);
            //if (!changeBoundary && m_refMesh->V[mesh.m_refIds[i]].isBoundary) continue;
            const double& x = X[3 * i + 0]; const double& y = X[3 * i + 1]; const double& z = X[3 * i + 2];
            if (fabs(double(v.x - x)) > 1e-6 || fabs(double(v.y - y)) > 1e-6 || fabs(double(v.z - z)) > 1e-6) converged = false;
            if (!allowBigStep) if (fabs(double(v.x - x)) > maxEdgeLength || fabs(double(v.y - y)) > maxEdgeLength || fabs(double(v.z - z)) > maxEdgeLength) continue;
            v.x = stepSize * x + (1.0 - stepSize) * v.x;
            v.y = stepSize * y + (1.0 - stepSize) * v.y;
            v.z = stepSize * z + (1.0 - stepSize) * v.z;
        }
        //mesh.ProjectTo(targetMesh);

//        double minimumScaledJacobian = 0.0;
//        size_t InvertedElements = GetMinScaledJacobianVerdict(mesh, minimumScaledJacobian, this->minScaledJacobian);
//
//        if (InvertedElements > m_numOfInvertdElements) {
//            std::cout << "Recover previous mesh\n";
//            for (size_t i = 0; i < mesh.V.size(); i++) {
//                mesh.V[i].x = oldV[i].x;
//                mesh.V[i].y = oldV[i].y;
//                mesh.V[i].z = oldV[i].z;
//            }
//        }
    }
    else {
        for (size_t i = 0; i < mesh.V.size(); i++) {
            if (fabs(double(mesh.V[i].x - X[3 * i + 0])) > 1e-6 || fabs(double(mesh.V[i].y - X[3 * i + 1])) > 1e-6 || fabs(double(mesh.V[i].z - X[3 * i + 2])) > 1e-6) {
                converged = false;
            }
            if (!allowBigStep) if (fabs(double(mesh.V[i].x - X[3 * i + 0])) > 2.0 * avgMeshEdgeLength || fabs(double(mesh.V[i].y - X[3 * i + 1])) > 2.0 * avgMeshEdgeLength
                    || fabs(double(mesh.V[i].z - X[3 * i + 2])) > 2.0 * avgMeshEdgeLength) {
                continue;
            }
            mesh.V[i].x = stepSize * X[3 * i + 0] + (1.0 - stepSize) * mesh.V[i].x;
            mesh.V[i].y = stepSize * X[3 * i + 1] + (1.0 - stepSize) * mesh.V[i].y;
            mesh.V[i].z = stepSize * X[3 * i + 2] + (1.0 - stepSize) * mesh.V[i].z;
        }
    }
    return converged;
    return converged;
}

size_t SurfaceMeshOpt::OptimizeSurfaceVertices(std::vector<Trip>& A_Entries, std::vector<float>& b, size_t& row)
{
    size_t a_id = 0;
    for (size_t i = 0; i < mesh.V.size(); i++) {
        const Vertex& v = mesh.V.at(i);
        if (!v.isBoundary)
            continue;
        if (v.type == REGULAR) {
            //beta(nv + d)
            const glm::dvec3& implicit_n = v.normal;
            float implicit_d = glm::dot(implicit_n, v.xyz());
            for (size_t k = 0; k < 3; k++) {
                A_Entries.push_back(Trip(row, 3 * i + k, beta * implicit_n[k]));
            }
            b.push_back(beta * implicit_d);
            row++;
        } else if (v.type == FEATURE) {
            int a_index = a_id + 3 * mesh.V.size();
            for (size_t k = 0; k < 3; k++) {
                A_Entries.push_back(Trip(row, 3 * i + k, alpha));
                A_Entries.push_back(Trip(row, a_index, -alpha * v.tangent[k]));
                b.push_back(alpha * v[k]);
                row++;
            }
        } else if (v.type == CORNER) {
            for (size_t k = 0; k < 3; k++) {
                A_Entries.push_back(Trip(row, 3 * i + k, alpha));
                b.push_back(alpha * v[k]);
                row++;
            }
        }
    }
    a_id = 0;
    for (size_t i = 0; i < mesh.V.size(); i++) {
        const Vertex& v = mesh.V.at(i);
        if (v.isBoundary && v.type == FEATURE) {
            int a_index = a_id++ + 3 * mesh.V.size();
            A_Entries.push_back(Trip(row, a_index, alpha));
            b.push_back(0);
            row++;
        }
    }

    return a_id;
}

void SurfaceMeshOpt::OptimizeEdgeOrthogonality(std::vector<Trip>& A_Entries, std::vector<float>& b, size_t& row)
{
    for (size_t k = 0; k < mesh.E.size(); k++){
        Edge& edge = mesh.E.at(k);
        edge.energyOrthogonality = 0;
    }
    double energy = 0.0;
    for (size_t k = 0; k < mesh.E.size(); k++){
        /*const */Edge& edge1 = mesh.E.at(k);
        for (size_t j = 0; j < edge1.orthogonalEids.size(); j++) {
            const size_t edgeId2 = edge1.orthogonalEids.at(j);
            /*const */ Edge& edge2 = mesh.E.at(edgeId2);
            const size_t vId1_1 = edge1.Vids.at(0);
            const size_t vId1_2 = edge1.Vids.at(1);
            const size_t vId2_1 = edge2.Vids.at(0);
            const size_t vId2_2 = edge2.Vids.at(1);
            size_t shareVId = vId1_1;
            if (vId1_2 == vId2_1 || vId1_2 == vId2_2)
                shareVId = vId1_2;
            const Vertex& vetex1 = shareVId == vId1_1 ? mesh.V.at(vId1_2) : mesh.V.at(vId1_1);
            const Vertex& vetexC = mesh.V.at(shareVId);
            const Vertex& vetex2 = shareVId == vId2_1 ? mesh.V.at(vId2_2) : mesh.V.at(vId2_1);
            const glm::dvec3 v2 = glm::dvec3(vetexC.x - vetex2.x, vetexC.y - vetex2.y, vetexC.z - vetex2.z);
            const glm::dvec3 v2n = glm::normalize(v2);
            double length_v1 = edge1.length;
            if (!useAverageTargetLength) {
                length_v1 = glm::length(glm::dvec3(vetexC.x - vetex1.x, vetexC.y - vetex1.y, vetexC.z - vetex1.z));
                // length_v1 = length_v1 > 1e-4 ? length_v1 : 1e-4;
                length_v1 = length_v1 > avgMeshEdgeLength*anisotropy ? length_v1 : avgMeshEdgeLength*anisotropy;
                length_v1 = std::isnan(length_v1) ? avgMeshEdgeLength*anisotropy : length_v1;
            }
            ///////////////////////////////////////////////
            const double length_v1_ = 1.0/length_v1;
            const glm::dvec3 v1n = glm::dvec3((vetexC.x - vetex1.x)*length_v1_, (vetexC.y - vetex1.y)*length_v1_, (vetexC.z - vetex1.z)*length_v1_);
            const double energy_orthognality = (glm::dot(v1n, v2n) - 0.0) * (glm::dot(v1n, v2n) - 0.0) * gamma;
            energy += energy_orthognality;
            edge1.energyOrthogonality += energy_orthognality;
            ///////////////////////////////////////////////
            // v1 is variable, and use v2 as constant
            // <v1/length_v1, v2n> = 0.0
            // double weight = energy_orthognality > gamma ? energy_orthognality : gamma;
            // if (stepSize > 0.9)
            double weight = gamma;
            for (size_t n = 0; n < 3; n++) {
                double a = weight * v2n[n] / length_v1;
                a = std::isnan(a) ? 1.0 : a;
                A_Entries.push_back(Trip(row, 3 * vetexC.id + n, a));
                A_Entries.push_back(Trip(row, 3 * vetex1.id + n, -a));
            }
            b.push_back(0.0 * weight);
            row++;
        }
    }
    EOrthogonality.push_back(energy);
}

void SurfaceMeshOpt::OptimizeEdgeStraightness(std::vector<Trip>& A_Entries, std::vector<float>& b, size_t& row)
{
    for (size_t k = 0; k < mesh.E.size(); k++){
        Edge& edge = mesh.E.at(k);
        edge.energyStraightness = 0;
    }
    double initgamma = gamma;
    gamma = 1.0;
    double energy = 0.0;
    for (size_t k = 0; k < mesh.E.size(); k++) {
        /*const */ Edge& edge1 = mesh.E.at(k);
        for (size_t j = 0; j < edge1.consecutiveEids.size(); j++) {
            const size_t edgeId2 = edge1.consecutiveEids.at(j);
            /*const */ Edge& edge2 = mesh.E.at(edgeId2);
            const size_t vId1_1 = edge1.Vids.at(0);
            const size_t vId1_2 = edge1.Vids.at(1);
            const size_t vId2_1 = edge2.Vids.at(0);
            const size_t vId2_2 = edge2.Vids.at(1);
            size_t shareVId = vId1_1;
            if (vId1_2 == vId2_1 || vId1_2 == vId2_2)
                shareVId = vId1_2;
            const Vertex& vetex1 = shareVId == vId1_1 ? mesh.V.at(vId1_2) : mesh.V.at(vId1_1);
            const Vertex& vetexC = mesh.V.at(shareVId);
            const Vertex& vetex2 = shareVId == vId2_1 ? mesh.V.at(vId2_2) : mesh.V.at(vId2_1);
            const glm::dvec3 v2 = glm::dvec3(vetexC.x - vetex2.x, vetexC.y - vetex2.y, vetexC.z - vetex2.z);
            const glm::dvec3 v2n = glm::normalize(v2);
            double length_v1 = edge1.length;
            if (!useAverageTargetLength) {
                length_v1 = glm::length(glm::dvec3(vetexC.x - vetex1.x, vetexC.y - vetex1.y, vetexC.z - vetex1.z));
                // length_v1 = length_v1 > 1e-4 ? length_v1 : 1e-4;
                length_v1 = length_v1 > avgMeshEdgeLength*anisotropy ? length_v1 : avgMeshEdgeLength*anisotropy;
                length_v1 = std::isnan(length_v1) ? avgMeshEdgeLength*anisotropy : length_v1;
            }
            ///////////////////////////////////////////////
            const double length_v1_ = 1.0/length_v1;
            const glm::dvec3 v1n = glm::dvec3((vetexC.x - vetex1.x)*length_v1_, (vetexC.y - vetex1.y)*length_v1_, (vetexC.z - vetex1.z)*length_v1_);
            const double energy_straightness = (glm::dot(v1n, v2n) + 1.0) * (glm::dot(v1n, v2n) + 1.0) * gamma;
            energy += energy_straightness;
            edge1.energyStraightness += energy_straightness;
            ///////////////////////////////////////////////
            // v1 is variable, and use v2 as constant
            // <v1/length_v1, v2n> = -1
            // double weight = energy_straightness * 10 > gamma ? energy_straightness  * 10: gamma;
            // if (stepSize > 0.9)
            double weight = gamma;
            for (size_t n = 0; n < 3; n++) {
                double a = weight * v2n[n] / length_v1;
                a = std::isnan(a) ? 1.0 : a;
                A_Entries.push_back(Trip(row, 3 * vetexC.id + n, a));
                A_Entries.push_back(Trip(row, 3 * vetex1.id + n, -a));
            }
            b.push_back(-1.0 * weight);
            row++;
        }
    }
    EStraightness.push_back(energy);
    gamma = initgamma;
}

void SurfaceMeshOpt::OptimizeSingularity(std::vector<Trip>& A_Entries, std::vector<float>& b, size_t& row)
{
    for (size_t k = 0; k < mesh.E.size(); k++){
        Edge& edge = mesh.E.at(k);
        edge.energySingularity = 0;
    }
    double energy = 0.0;
    for (size_t k = 0; k < mesh.E.size(); k++){
        const Edge& edge = mesh.E.at(k);
        if (!edge.isSingularity)
            continue;
        if (edge.N_Cids.size() != 5 && edge.N_Cids.size() != 3)
            continue;
        for (size_t j = 0; j < edge.orthogonalEids.size(); j++) {
            for (size_t i = 0; i < edge.orthogonalEids.size(); i++) {
                if (j == i)
                    continue;
                const size_t edgeId1 = edge.orthogonalEids.at(j);
                const size_t edgeId2 = edge.orthogonalEids.at(i);
                /*const */ Edge& edge1 = mesh.E.at(edgeId1);
                /*const */ Edge& edge2 = mesh.E.at(edgeId2);
                const size_t vId1_1 = edge1.Vids.at(0);
                const size_t vId1_2 = edge1.Vids.at(1);
                const size_t vId2_1 = edge2.Vids.at(0);
                const size_t vId2_2 = edge2.Vids.at(1);
                if (vId1_1 != vId2_1 &&
                    vId1_1 != vId2_2 &&
                    vId1_2 != vId2_1 &&
                    vId1_2 != vId2_2)
                    continue;

                bool inTheSameCell = false;
                for (size_t n = 0; n < edge.N_Cids.size(); n++)
                    if (Util::IsEdgeInCell(mesh, edge.N_Cids.at(n), edgeId1) && Util::IsEdgeInCell(mesh, edge.N_Cids.at(n), edgeId2)) {
                        inTheSameCell = true;
                        break;
                    }

                size_t shareVId = vId1_1;
                if (vId1_2 == vId2_1 || vId1_2 == vId2_2)
                    shareVId = vId1_2;
                const Vertex& vetex1 = shareVId == vId1_1 ? mesh.V.at(vId1_2) : mesh.V.at(vId1_1);
                const Vertex& vetexC = mesh.V.at(shareVId);
                const Vertex& vetex2 = shareVId == vId2_1 ? mesh.V.at(vId2_2) : mesh.V.at(vId2_1);
                const glm::dvec3 v2 = glm::dvec3(vetexC.x - vetex2.x, vetexC.y - vetex2.y, vetexC.z - vetex2.z);
                const glm::dvec3 v2n = glm::normalize(v2);
                double length_v1 = edge1.length;
                if (!useAverageTargetLength) {
                    length_v1 = glm::length(glm::dvec3(vetexC.x - vetex1.x, vetexC.y - vetex1.y, vetexC.z - vetex1.z));
                    // length_v1 = length_v1 > 1e-4 ? length_v1 : 1e-4;
                    length_v1 = length_v1 > avgMeshEdgeLength*anisotropy ? length_v1 : avgMeshEdgeLength*anisotropy;
                    length_v1 = std::isnan(length_v1) ? avgMeshEdgeLength*anisotropy : length_v1;
                }
                double b_term = -0.5;
                if (edge.N_Cids.size() == 3)
                    ; // b.push_back(-0.5 * gamma);
                else if (edge.N_Cids.size() == 5) {
                    if (inTheSameCell)
                        b_term = 0.309;
                    else
                        b_term = -0.809;
                }
                else
                    std::cout << "Error! edge.N_Cids.size() = " << edge.N_Cids.size() << std::endl;
                ///////////////////////////////////////////////
                const double length_v1_ = 1.0/length_v1;
                const glm::dvec3 v1n = glm::dvec3((vetexC.x - vetex1.x)*length_v1_, (vetexC.y - vetex1.y)*length_v1_, (vetexC.z - vetex1.z)*length_v1_);
                const double energy_sigularity = (glm::dot(v1n, v2n) - b_term) * (glm::dot(v1n, v2n) - b_term) * gamma;
                energy += energy_sigularity;
                edge1.energySingularity += energy_sigularity;
                edge2.energySingularity += energy_sigularity;
                ///////////////////////////////////////////////
                // v1 is variable, and use v2 as constant
                // <v1/length_v1, v2n> = -0.5
                // double weight = energy_sigularity > gamma ? energy_sigularity : gamma;
                // if (stepSize > 0.9)
                double weight = gamma;
                if (stepSize < 0.5)
                    weight = gamma * 2;
                if (m_numOfInvertdElements == 0)
                    weight = gamma * 10;
                for (size_t n = 0; n < 3; n++) {
                    double a = weight * v2n[n] / length_v1;
                    a = std::isnan(a) ? 1.0 : a;
                    A_Entries.push_back(Trip(row, 3 * vetexC.id + n, a));
                    A_Entries.push_back(Trip(row, 3 * vetex1.id + n, -a));
                }
                b.push_back(b_term * weight);
                row++;
            }
        }
    }
    ESingularity.push_back(energy);
}

