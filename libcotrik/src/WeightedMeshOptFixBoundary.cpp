/*
 * WeightedMeshOptFixBoundary.cpp
 *
 *  Created on: Nov 25, 2016
 *      Author: cotrik
 */

#include "WeightedMeshOptFixBoundary.h"
#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "MeshQuality.h"

#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <queue>
#include <iostream>
#include <cmath>
#include <math.h>
#include "verdict.h"

//MeshOpt::MeshOpt()
WeightedMeshOptFixBoundary::WeightedMeshOptFixBoundary(const Mesh& mesh)
: MeshOptFixBoundary(mesh)
{
    // TODO Auto-generated constructor stub

}

WeightedMeshOptFixBoundary::~WeightedMeshOptFixBoundary()
{
    // TODO Auto-generated destructor stub
}

size_t WeightedMeshOptFixBoundary::Run(const size_t iters/* = 1*/)
{
    //----------------------------------------------
    double minimumScaledJacobian = 0.0;
    m_numOfInvertdElements = GetMinScaledJacobianVerdict(mesh, minimumScaledJacobian, this->minScaledJacobian);
    std::cout << "iter = " << 0 << " #inverted = " << m_numOfInvertdElements << " MSJ = " << minimumScaledJacobian << std::endl;
    //----------------------------------------------
    int iter = 0;
    double prevMinimumScaledJacobian = -1.0;
    bool converged = false;
    double initStepSize = stepSize;
    bool initUseAverageTargetLength = useAverageTargetLength;
    bool untangled = false;
    MeshOpt::ComputeMeshTargetLength();
    static int global_count = 0;
    global_count++;
    double last_total_energy = 100000000000000000;
    double eps = 1e-2;
    while (!converged && iter++ < iters) {
        if (useAverageTargetLength) MeshOpt::ComputeMeshTargetLength();
        MeshOptFixBoundary::Optimize();
        double total_energy = EOrthogonality.at(iter - 1) + EStraightness.at(iter - 1) + ESingularity.at(iter - 1);
        double rate = fabs(total_energy - last_total_energy) / last_total_energy;
        //if (rate < eps && iter > 1) converged = true;
        last_total_energy = total_energy;
        stepSize *= initStepSize;

        m_numOfInvertdElements = GetMinScaledJacobianVerdict(mesh, minimumScaledJacobian, this->minScaledJacobian);

        if (iter > 1) std::cout << "iter = " << iter << " #inverted = " << m_numOfInvertdElements << " MSJ = " << minimumScaledJacobian << " converge = " << rate << std::endl;
        else
            std::cout << "iter = " << iter << " #inverted = " << m_numOfInvertdElements << " MSJ = " << minimumScaledJacobian << std::endl;
        if (prevMinimumScaledJacobian > this->minScaledJacobian && minimumScaledJacobian <= prevMinimumScaledJacobian) {
            std::cout << "*************************" << std::endl;
            std::cout << "\033[1;32mBest Mesh is " << "opt.vtk" << "\033[0m" << std::endl;
            std::cout << "*************************" << std::endl;
            MeshFileWriter optwriter(bestmesh, "opt.vtk");
            optwriter.WriteFile();
            untangled = true;
            break;
        }
        else if (converged) {
            std::cout << "*************************" << std::endl;
            std::cout << "Converged at iter " << iter << std::endl;
            std::cout << "*************************" << std::endl;
            MeshFileWriter optwriter(bestmesh, "converged.vtk");
            optwriter.WriteFile();
        }
        prevMinimumScaledJacobian = minimumScaledJacobian;
        bestmesh = mesh;
    }

    mesh = bestmesh;
    EOrthogonality.clear();
    EStraightness.clear();
    ESingularity.clear();
    converged = false;
    iter = 0;
    initStepSize = 0.98;
    eps = 1e-6;
    prevMinimumScaledJacobian = -1.0;
    while (!converged && iter++ < iters) {
        ComputeMeshTargetLength();
        Optimize();
        double total_energy = EOrthogonality.at(iter - 1) + EStraightness.at(iter - 1) + ESingularity.at(iter - 1);
        double rate = 0;
        rate = fabs(total_energy - last_total_energy) / last_total_energy;
        if (rate < eps) converged = true;
        last_total_energy = total_energy;
        stepSize *= initStepSize;

        m_numOfInvertdElements = GetMinScaledJacobianVerdict(mesh, minimumScaledJacobian, this->minScaledJacobian);

        std::cout << "iter = " << iter << " #inverted = " << m_numOfInvertdElements << " MSJ = " << minimumScaledJacobian << " converge = " << rate << std::endl;
//        if (prevMinimumScaledJacobian > this->minScaledJacobian && minimumScaledJacobian <= prevMinimumScaledJacobian) {
//            std::cout << "*************************" << std::endl;
//            std::cout << "\033[1;32mBest Mesh is " << "opt_further.vtk" << "\033[0m" << std::endl;
//            std::cout << "*************************" << std::endl;
//            MeshFileWriter optwriter(bestmesh, "opt_further.vtk");
//            optwriter.WriteFile();
//            untangled = true;
//            break;
//        }
//        else if (converged) {
//            std::cout << "*************************" << std::endl;
//            std::cout << "Converged at iter " << iter << std::endl;
//            std::cout << "*************************" << std::endl;
//        }
//        prevMinimumScaledJacobian = minimumScaledJacobian;
//        bestmesh = mesh;
    }

    MeshFileWriter optwriter(mesh, "last.vtk");
    optwriter.WriteFile();
    return iter - 1;
}

bool WeightedMeshOptFixBoundary::Optimize()
{
    //std::cout << "Constraint-------\t#Entries\t #row\n";
    //std::cout << "-----------------\t---------------------------\t\n";
    std::vector<Trip> A_Entries;
    std::vector<float> b;
    size_t row = 0;
    size_t currentRow = row;
    size_t currentNumOfA_Entries = A_Entries.size();

    currentRow = row;
    currentNumOfA_Entries = A_Entries.size();
    OptimizeSingularity(A_Entries, b, row);
    //std::cout << "EdgeSingularity--\t" << A_Entries.size() - currentNumOfA_Entries << "\t -------- " << row - currentRow << std::endl;

    currentRow = row;
    currentNumOfA_Entries = A_Entries.size();
    OptimizeEdgeOrthogonality(A_Entries, b, row);
    //std::cout << "EdgeOrthogonality\t" << A_Entries.size() - currentNumOfA_Entries << "\t -------- " << row - currentRow << std::endl;

    currentRow = row;
    currentNumOfA_Entries = A_Entries.size();
    OptimizeEdgeStraightness(A_Entries, b, row);
    //std::cout << "EdgeStraightness-\t" << A_Entries.size() - currentNumOfA_Entries << "\t -------- " << row - currentRow << std::endl;

    currentRow = row;
    currentNumOfA_Entries = A_Entries.size();
    //size_t a_id = OptimizeSurfaceVertices(A_Entries, b, row);
    size_t a_id = 0;OptimizeBoundaryVertices(A_Entries, b, row);
    //std::cout << "SurfaceVertices--\t" << A_Entries.size() - currentNumOfA_Entries << "\t -------- " << row - currentRow << std::endl;

    const size_t col = 3 * mesh.V.size() + a_id;

    //std::cout << "----------- A Info ------------\n";
    //std::cout << "Entries = " << A_Entries.size() << " row = " << row << " col = " << col << std::endl;

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
        std::vector<glm::dvec3> oldV(mesh.V.size());
        for (size_t i = 0; i < mesh.V.size(); i++) {
            oldV[i].x = mesh.V[i].x;
            oldV[i].y = mesh.V[i].y;
            oldV[i].z = mesh.V[i].z;
        }

        const double maxEdgeLength = 2.0 * avgMeshEdgeLength;
        for (size_t i = 0; i < mesh.V.size(); i++) {
            Vertex& v = mesh.V.at(i);
            //if (!changeBoundary && m_refMesh->V[mesh.m_refIds[i]].isBoundary)
            if (v.isBoundary)
                continue;
            const double& x = X[3 * i + 0];
            const double& y = X[3 * i + 1];
            const double& z = X[3 * i + 2];
            if (fabs(double(v.x - x)) > 1e-6 || fabs(double(v.y - y)) > 1e-6 || fabs(double(v.z - z)) > 1e-6)
                converged = false;
            if (!allowBigStep)
                if (fabs(double(v.x - x)) > maxEdgeLength || fabs(double(v.y - y)) > maxEdgeLength || fabs(double(v.z - z)) > maxEdgeLength)
                    continue;
            v.x = stepSize * x + (1.0 - stepSize) * v.x;
            v.y = stepSize * y + (1.0 - stepSize) * v.y;
            v.z = stepSize * z + (1.0 - stepSize) * v.z;
        }
        //mesh.ProjectTo(targetMesh);

        double minimumScaledJacobian = 0.0;
        size_t InvertedElements = GetMinScaledJacobianVerdict(mesh, minimumScaledJacobian, this->minScaledJacobian);

        if (InvertedElements > m_numOfInvertdElements) {
            std::cout << "Recover previous mesh\n";
            for (size_t i = 0; i < mesh.V.size(); i++) {
                mesh.V[i].x = oldV[i].x;
                mesh.V[i].y = oldV[i].y;
                mesh.V[i].z = oldV[i].z;
            }
        }
    }
    else {
        const double maxEdgeLength = 2.0 * avgMeshEdgeLength;
        for (size_t i = 0; i < mesh.V.size(); i++) {
            Vertex& v = mesh.V.at(i);
            if (v.isBoundary) continue;
            const double& x = X[3 * i + 0];
            const double& y = X[3 * i + 1];
            const double& z = X[3 * i + 2];
            if (fabs(double(v.x - x)) > 1e-6 || fabs(double(v.y - y)) > 1e-6 || fabs(double(v.z - z)) > 1e-6) converged = false;
            if (!allowBigStep) if (fabs(double(v.x - x)) > maxEdgeLength || fabs(double(v.y - y)) > maxEdgeLength || fabs(double(v.z - z)) > maxEdgeLength) continue;
            v.x = stepSize * x + (1.0 - stepSize) * v.x;
            v.y = stepSize * y + (1.0 - stepSize) * v.y;
            v.z = stepSize * z + (1.0 - stepSize) * v.z;
        }
    }
    return converged;
}

void WeightedMeshOptFixBoundary::OptimizeSingularity(std::vector<Trip>& A_Entries, std::vector<float>& b, size_t& row)
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
        if (edge.isBoundary)
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
                    if (IsEdgeInCell(mesh, edge.N_Cids.at(n), edgeId1) && IsEdgeInCell(mesh, edge.N_Cids.at(n), edgeId2)) {
                        inTheSameCell = true;
                        break;
                    }

                size_t shareVId = vId1_1;
                if (vId1_2 == vId2_1 || vId1_2 == vId2_2)
                    shareVId = vId1_2;
                const Vertex& vetex1 = shareVId == vId1_1 ? mesh.V.at(vId1_2) : mesh.V.at(vId1_1);
                const Vertex& vetexC = mesh.V.at(shareVId);
                const Vertex& vetex2 = shareVId == vId2_1 ? mesh.V.at(vId2_2) : mesh.V.at(vId2_1);
                if (vetexC.isBoundary)
                    continue;
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

void WeightedMeshOptFixBoundary::OptimizeEdgeOrthogonality(std::vector<Trip>& A_Entries, std::vector<float>& b, size_t& row)
{
    for (size_t k = 0; k < mesh.E.size(); k++){
        Edge& edge = mesh.E.at(k);
        edge.energyOrthogonality = 0;
    }
    double energy = 0.0;
    for (size_t k = 0; k < mesh.E.size(); k++){
        /*const */Edge& edge1 = mesh.E.at(k);
//        if (edge1.isBoundary)
//            continue;
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
//            if (vetexC.isBoundary)
//                continue;
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

void WeightedMeshOptFixBoundary::OptimizeEdgeStraightness(std::vector<Trip>& A_Entries, std::vector<float>& b, size_t& row)
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
//        if (edge1.isBoundary)
//            continue;
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

void WeightedMeshOptFixBoundary::ComputeMeshTargetLength()
{
    std::cout << "=============================\n";
    std::cout << "   Computing TargetLength    \n";

    std::vector<Trip> coefficients;
    std::vector<float> b;
    size_t row  = 0;
    size_t numOfparallelEids = 0;
    size_t numOfconsecutiveEids = 0;
    for (size_t i = 0; i < mesh.E.size(); i++) {
        coefficients.push_back(Trip(row, i, 1));
        if (mesh.E[i].isBoundary) b.push_back(mesh.E[i].length);
        else
            b.push_back(avgMeshEdgeLength);
        row++;

        for (size_t j = 0; j < mesh.E[i].parallelEids.size(); j++) {
            size_t N_Eid = mesh.E[i].parallelEids[j];
            coefficients.push_back(Trip(row, i, 1));
            coefficients.push_back(Trip(row, N_Eid, -1));
            b.push_back(0);
            numOfparallelEids++;
            row++;
        }
        for (size_t j = 0; j < mesh.E[i].consecutiveEids.size(); j++) {
            size_t N_Eid = mesh.E[i].consecutiveEids[j];
            coefficients.push_back(Trip(row, i, 1));
            coefficients.push_back(Trip(row, N_Eid, -1));
            b.push_back(0);
            numOfconsecutiveEids++;
            row++;
        }
    }

    for (size_t k = 0; k < mesh.E.size(); k++) {
        const Edge& edge = mesh.E.at(k);
        if (!edge.isSingularity) continue;
        if (edge.N_Cids.size() != 5 && edge.N_Cids.size() != 3 && edge.N_Cids.size() != 6) continue;
        for (size_t j = 0; j < edge.orthogonalEids.size(); j++) {
            for (size_t i = 0; i < edge.orthogonalEids.size(); i++) {
                if (j == i) continue;
                const size_t edgeId1 = edge.orthogonalEids.at(j);
                const size_t edgeId2 = edge.orthogonalEids.at(i);
                /*const */Edge& edge1 = mesh.E.at(edgeId1);
                /*const */Edge& edge2 = mesh.E.at(edgeId2);
                const size_t vId1_1 = edge1.Vids.at(0);
                const size_t vId1_2 = edge1.Vids.at(1);
                const size_t vId2_1 = edge2.Vids.at(0);
                const size_t vId2_2 = edge2.Vids.at(1);
                if (vId1_1 != vId2_1 && vId1_1 != vId2_2 && vId1_2 != vId2_1 && vId1_2 != vId2_2) continue;

                coefficients.push_back(Trip(row, edge1.id, 1));
                coefficients.push_back(Trip(row, edge2.id, -1));
                b.push_back(0);
                row++;
            }
        }
    }

    size_t col = mesh.E.size();

    Eigen::SparseMatrix<float> A(row, col);
    A.setFromTriplets(coefficients.begin(), coefficients.end());
    Eigen::SparseMatrix<float> A_T = A.transpose();
    Eigen::SparseMatrix<float> ATA = A_T * A;
    //  MatrixXf ;
    VectorXf X(col), B(row);
    for (size_t i = 0; i < row; i++)
        B(i) = b[i];

    VectorXf ATB = A_T * B;
    Eigen::SimplicialLDLT<SpMat> chol(ATA);
    X = chol.solve(ATB);
//    for (size_t i = 0; i < mesh.E.size(); i++)
//        mesh.E[i].length = X[i];

    const double exp = 2.7182818284590452353602875;
    const double sigma = 1.0/0.15;
    std::vector<double> w(mesh.E.size(), 1);
    double coordinates[8][3];
    for (size_t i = 0; i < mesh.E.size(); i++) {
        const Edge& e = mesh.E.at(i);
        if (e.isBoundary) continue;
        double minSJ = 1.0;
        for (size_t k = 0; k < e.N_Cids.size(); k++) {
            const Cell& c = mesh.C.at(e.N_Cids.at(k));
            for (size_t j = 0; j < 8; j++) {
                const Vertex& v = mesh.V[c.Vids[j]];
                coordinates[j][0] = v.x;
                coordinates[j][1] = v.y;
                coordinates[j][2] = v.z;
            }
            double scaledJacobian = v_hex_scaled_jacobian(8, coordinates);
            if (scaledJacobian < minSJ)
                minSJ = scaledJacobian;
        }
        w.at(i) = pow(exp, -0.5* pow(minSJ * sigma, 2));
        //w.at(i) = sqrt(1 - minSJ);
        //std::cout << w.at(i) << " ";//std::endl;
    }

    for (size_t i = 0; i < mesh.E.size(); i++) {
        Edge& e = mesh.E[i];
        e.length = (1 - w[i]) * X[i] + w[i] * e.length;
    }
    std::cout << "=============================\n";
}
