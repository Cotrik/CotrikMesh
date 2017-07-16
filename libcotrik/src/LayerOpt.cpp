/*
 * LayerOpt.cpp
 *
 *  Created on: Nov 25, 2016
 *      Author: cotrik
 */

#include "MeshOpt.h"
#include "Mesh.h"
#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "MeshQuality.h"

#include <algorithm>
#include <vector>
#include <iostream>
#include <math.h>

//MeshOpt::MeshOpt()
LayerOpt::LayerOpt(const Mesh& mesh,
        const double alpha/* = 100.0*/, const double beta/* = 1.0*/, const double gamma/* = 1.0*/)
: mesh((Mesh&)mesh)
, alpha(alpha)
, beta(beta)
, gamma(gamma)
, anisotropy(0.05)
, avgMeshEdgeLength(1.0)
, stepSize(1.0)
, useAverageTargetLength(false)
, recoverable(true)
, m_numOfInvertdElements(1000000000)
{
    // TODO Auto-generated constructor stub

}

LayerOpt::~LayerOpt()
{
    // TODO Auto-generated destructor stub
}

void LayerOpt::Run(const size_t iters/* = 1*/)
{
    if (stepSize < 0.1){
        std::cout << "ERROR! stepSize < 0.1\n";
        return;
    }
    double sumEdgeLength = 0.0;
    size_t numOfBoundaryEdges = 0;
    for (size_t i = 0; i < mesh.E.size(); i++)
    {
        if (mesh.E.at(i).isBoundary){
            const Vertex& v1 = mesh.V.at(mesh.E.at(i).Vids.at(0));
            const Vertex& v2 = mesh.V.at(mesh.E.at(i).Vids.at(1));
            mesh.E[i].length = glm::length(v1 - v2);
            sumEdgeLength += mesh.E[i].length;
            //sumEdgeLength += glm::length(v1 - v2);
            numOfBoundaryEdges++;
        }
    }
    avgMeshEdgeLength = sumEdgeLength / numOfBoundaryEdges;
    std::cout << "Average Surface Edge Length = " << avgMeshEdgeLength << std::endl;

    int iter = 0;
    double prevMinimumScaledJacobian = -1.0;
    bool converged = false;
    double initStepSize = stepSize;
    bool initUseAverageTargetLength = useAverageTargetLength;
    while (!converged && iter++ < iters)
    {
        std::string filename = std::string("EdgesEnergy.") + std::to_string(iter) + ".vtk";
        MeshFileWriter edgeWriter(mesh, filename.c_str());
        edgeWriter.WriteEdgesVtk();

        if (!initUseAverageTargetLength  && iter == 1)
            useAverageTargetLength = true;
        std::cout << "stepSize = " << stepSize << std::endl;
        if (useAverageTargetLength)
            ComputeMeshTargetLength();
        converged = Optimize();
        if (!initUseAverageTargetLength && iter == 1)
            useAverageTargetLength = false;
//        if (iter % 30 == 0)
//            stepSize = initStepSize;
        //stepSize = initStepSize - (initStepSize - 0.1) / iters;
        stepSize *= initStepSize;

        filename = std::string("MeshOpt.") + std::to_string(iter) + ".vtk";
        MeshFileWriter writer(mesh, filename.c_str());
        writer.WriteFile();

        double minimumScaledJacobian = 0.0;
        double averageScaledJacobian = 0.0;
        double maximumScaledJacobian = 0.0;
        std::vector<size_t> badCellIds;
        m_numOfInvertdElements = GetQuality(filename.c_str(), minimumScaledJacobian, averageScaledJacobian, maximumScaledJacobian, badCellIds);
        filename = std::string("BadCells.") + std::to_string(iter) + ".vtk";
        OutputBadCells(badCellIds, filename.c_str());

        if (prevMinimumScaledJacobian > 0 && minimumScaledJacobian <= prevMinimumScaledJacobian)
        {
            filename = std::string("MeshOpt.") + std::to_string(iter - 1) + ".vtk";
            std::cout << "*************************" << std::endl;
            std::cout << "Best Mesh is " << filename << std::endl;
            std::cout << "*************************" << std::endl;
            MeshFileReader reader(filename.c_str());
            const Mesh& bestmesh = reader.GetMesh();
            MeshFileWriter optwriter(bestmesh, "opt.vtk");
            optwriter.WriteFile();
            break;
        }
        else if (converged)
        {
            std::cout << "*************************" << std::endl;
            std::cout << "Converged at " << filename << std::endl;
            std::cout << "*************************" << std::endl;
        }
        prevMinimumScaledJacobian = minimumScaledJacobian;
    }

    std::vector<double> E_total(ESingularity.size());
    for (size_t i = 0; i < ESingularity.size(); i++)
        E_total.at(i) = ESingularity.at(i) + EOrthogonality.at(i) + EStraightness.at(i);
    std::cout << "*************************" << std::endl;
    std::cout << "E_total";
    for (size_t i = 0; i < E_total.size(); i++)
        std::cout << "\t" << E_total.at(i);
    std::cout << "\nESingularity";
    for (size_t i = 0; i < ESingularity.size(); i++)
        std::cout << "\t" << ESingularity.at(i);
    std::cout << "\nEOrthogonality";
    for (size_t i = 0; i < EOrthogonality.size(); i++)
        std::cout << "\t" << EOrthogonality.at(i);
    std::cout << "\nEStraightness";
    for (size_t i = 0; i < EStraightness.size(); i++)
        std::cout << "\t" << EStraightness.at(i);
    std::cout << "\n*************************\n" << std::endl;
}

bool LayerOpt::Optimize()
{
    std::vector<Trip> A_Entries;
    std::vector<float> b;
    size_t row = 0;

    OptimizeBoundaryVertices(A_Entries, b, row);
    std::cout << "OptimizeBoundaryVertices" << std::endl;
    //OptimizeInnerVertices(A_Entries, b, row);
    //std::cout << "OptimizeFrameInnerVertices" << std::endl;
    //std::cout << "Boundary row = " << row << std::endl;
    //std::cout << "Boundary A_Entries = " << A_Entries.size() << std::endl;
//    if (m_numOfInvertdElements < 10)
//    OptimizeEdgeComformalty(A_Entries, b, row);
//    std::cout << "OptimizeEdgeComformalty" << std::endl;

    OptimizeSingularity(A_Entries, b, row);
    std::cout << "OptimizeSingularity" << std::endl;
    //std::cout << "Singularity row = " << row << std::endl;
    //std::cout << "Singularity A_Entries = " << A_Entries.size() << std::endl;
    OptimizeEdgeOrthogonality(A_Entries, b, row);
    std::cout << "OptimizeEdgeOrthogonality" << std::endl;
    //std::cout << "Orthogonality row = " << row << std::endl;
    //std::cout << "Orthogonality A_Entries = " << A_Entries.size() << std::endl;
    OptimizeEdgeStraightness(A_Entries, b, row);
    std::cout << "OptimizeEdgeStraightness" << std::endl;
    //std::cout << "Straightness row = " << row << std::endl;
    //std::cout << "Straightness A_Entries = " << A_Entries.size() << std::endl;
    const size_t col = 3 * mesh.V.size();
    //std::cout << " col = " << col << std::endl;

    Eigen::SparseMatrix<float> A(row, col);
    A.setFromTriplets(A_Entries.begin(), A_Entries.end());
//    std::ofstream ofs("myA.txt");
//    ofs << A;
    Eigen::SparseMatrix<float> AT = A.transpose();
    Eigen::SparseMatrix<float> ATA = AT * A;

    //  MatrixXf ;
    VectorXf B(row);
    for (size_t i = 0; i < row; i++)
        B(i) = b[i];

    VectorXf ATB = AT * B;

    Eigen::SparseLU<SpMat> chol(ATA);
    VectorXf X(col);
    X = chol.solve(ATB);

    bool converged = true;
    if (recoverable) {
        std::vector<glm::vec3> oldV(mesh.V.size());
        for (size_t i = 0; i < mesh.V.size(); i++)
        {
            oldV[i].x = mesh.V[i].x;
            oldV[i].y = mesh.V[i].y;
            oldV[i].z = mesh.V[i].z;
        }

        for (size_t i = 0; i < mesh.V.size(); i++)
        {
            if (fabs(double(mesh.V[i].x - X[3 * i + 0])) > 1e-6
                || fabs(double(mesh.V[i].y - X[3 * i + 1])) > 1e-6
                || fabs(double(mesh.V[i].z - X[3 * i + 2])) > 1e-6)
            {
                converged = false;
            }
            if (fabs(double(mesh.V[i].x - X[3 * i + 0])) > 2.0 * avgMeshEdgeLength
                || fabs(double(mesh.V[i].y - X[3 * i + 1])) > 2.0 * avgMeshEdgeLength
                || fabs(double(mesh.V[i].z - X[3 * i + 2])) > 2.0 * avgMeshEdgeLength)
            {
                continue;
            }
            mesh.V[i].x = stepSize * X[3 * i + 0] + (1.0 - stepSize) * mesh.V[i].x;
            mesh.V[i].y = stepSize * X[3 * i + 1] + (1.0 - stepSize) * mesh.V[i].y;
            mesh.V[i].z = stepSize * X[3 * i + 2] + (1.0 - stepSize) * mesh.V[i].z;
        }

        MeshFileWriter writer(mesh, "temp.vtk");
        writer.WriteFile();
        double minimumScaledJacobian = 0.0;
        double averageScaledJacobian = 0.0;
        double maximumScaledJacobian = 0.0;
        std::vector<size_t> badCellIds1;
        size_t InvertedElements = GetQuality("temp.vtk", minimumScaledJacobian, averageScaledJacobian, maximumScaledJacobian, badCellIds1);

        if (InvertedElements > m_numOfInvertdElements){
            std::cout << "Recover previous mesh\n";
            for (size_t i = 0; i < mesh.V.size(); i++)
            {
                mesh.V[i].x = oldV[i].x;
                mesh.V[i].y = oldV[i].y;
                mesh.V[i].z = oldV[i].z;
            }
        }
    }
    else {
        for (size_t i = 0; i < mesh.V.size(); i++)
        {
            if (fabs(double(mesh.V[i].x - X[3 * i + 0])) > 1e-6
                || fabs(double(mesh.V[i].y - X[3 * i + 1])) > 1e-6
                || fabs(double(mesh.V[i].z - X[3 * i + 2])) > 1e-6)
            {
                converged = false;
            }
            if (fabs(double(mesh.V[i].x - X[3 * i + 0])) > 2.0 * avgMeshEdgeLength
                || fabs(double(mesh.V[i].y - X[3 * i + 1])) > 2.0 * avgMeshEdgeLength
                || fabs(double(mesh.V[i].z - X[3 * i + 2])) > 2.0 * avgMeshEdgeLength)
            {
                continue;
            }
            mesh.V[i].x = stepSize * X[3 * i + 0] + (1.0 - stepSize) * mesh.V[i].x;
            mesh.V[i].y = stepSize * X[3 * i + 1] + (1.0 - stepSize) * mesh.V[i].y;
            mesh.V[i].z = stepSize * X[3 * i + 2] + (1.0 - stepSize) * mesh.V[i].z;
        }
    }
    return converged;
}

bool IsInFace(const Mesh& mesh, const Edge& edge, const size_t vid1, const size_t vid2)
{
    bool bRet = false;
    for (size_t i = 0; i < edge.N_Cids.size(); i++) {
        const Cell& cell = mesh.C.at(edge.N_Cids.at(i));
        for (size_t j = 0; j < cell.Fids.size(); j++) {
            const Face& face = mesh.F.at(cell.N_Fids.at(j));
            bool isHaveVid1 = false;
            bool isHaveVid2 = false;
            for (size_t k = 0; k < face.Vids.size(); k++) {
                const size_t vid = face.Vids.at(k);
                if (vid == vid1)
                    isHaveVid1 = true;
                if (vid == vid2)
                    isHaveVid2 = true;
            }
            if (isHaveVid1 && isHaveVid2)
                return true;
        }
    }
    return bRet;
}
void LayerOpt::OptimizeEdgeComformalty(std::vector<Trip>& A_Entries, std::vector<float>& b, size_t& row)
{
//    for (size_t k = 0; k < mesh.E.size(); k++){
//        Edge& edge = mesh.E.at(k);
//        edge.energyStraightness = 0;
//    }

    double energy = 0.0;
    for (size_t k = 0; k < mesh.E.size(); k++) {
        /*const */ Edge& edge2 = mesh.E.at(k);
        if (!edge2.isBoundary)
            continue;
        for (size_t j = 0; j < edge2.parallelEids.size(); j++) {
            const size_t edgeId1 = edge2.parallelEids.at(j);
            /*const */ Edge& edge1 = mesh.E.at(edgeId1);
            size_t vId1_1 = edge1.Vids.at(0);
            size_t vId1_2 = edge1.Vids.at(1);
            size_t vId2_1 = edge2.Vids.at(0);
            size_t vId2_2 = edge2.Vids.at(1);
            if (!IsInFace(mesh, edge1, vId1_1, vId2_1))
                std::swap(vId1_1, vId1_2);

            const Vertex& vetex1_1 = mesh.V.at(vId1_1);
            const Vertex& vetex1_2 = mesh.V.at(vId1_2);
            const Vertex& vetex2_1 = mesh.V.at(vId2_1);
            const Vertex& vetex2_2 = mesh.V.at(vId2_2);

            const glm::vec3 v2 = glm::vec3(vetex2_1.x - vetex2_2.x, vetex2_1.y - vetex2_2.y, vetex2_1.z - vetex2_2.z);
            const glm::vec3 v2n = glm::normalize(v2);
            double length_v1 = edge1.length;
            if (!useAverageTargetLength) {
                length_v1 = glm::length(glm::vec3(vetex1_1.x - vetex1_2.x, vetex1_1.y - vetex1_2.y, vetex1_1.z - vetex1_2.z));
                // length_v1 = length_v1 > 1e-4 ? length_v1 : 1e-4;
                length_v1 = length_v1 > avgMeshEdgeLength*anisotropy ? length_v1 : avgMeshEdgeLength*anisotropy;
            }
            ///////////////////////////////////////////////
//            const double length_v1_ = 1.0/length_v1;
//            const glm::vec3 v1n = glm::vec3((vetex1_1.x - vetex1_2.x)*length_v1_, (vetex1_1.y - vetex1_2.y)*length_v1_, (vetex1_1.z - vetex1_2.z)*length_v1_);
//            const double energy_straightness = (glm::dot(v1n, v2n) + 1.0) * (glm::dot(v1n, v2n) + 1.0) * gamma;
//            energy += energy_straightness;
//            edge1.energyStraightness += energy_straightness;
            ///////////////////////////////////////////////
            // v1 is variable, and use v2 as constant
            // <v1/length_v1, v2n> = 1
            // double weight = energy_straightness * 10 > gamma ? energy_straightness  * 10: gamma;
            // if (stepSize > 0.9)
            double weight = 2*gamma;
            for (size_t n = 0; n < 3; n++) {
                const double a = weight * v2n[n] / length_v1;
                A_Entries.push_back(Trip(row, 3 * vetex1_1.id + n, a));
                A_Entries.push_back(Trip(row, 3 * vetex1_2.id + n, -a));
            }
            b.push_back(1.0 * weight);
            row++;
        }
    }
//    EStraightness.push_back(energy);
}

void LayerOpt::OptimizeSingularity(std::vector<Trip>& A_Entries, std::vector<float>& b, size_t& row)
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
                const glm::vec3 v2 = glm::vec3(vetexC.x - vetex2.x, vetexC.y - vetex2.y, vetexC.z - vetex2.z);
                const glm::vec3 v2n = glm::normalize(v2);
                double length_v1 = edge1.length;
                if (!useAverageTargetLength) {
                    length_v1 = glm::length(glm::vec3(vetexC.x - vetex1.x, vetexC.y - vetex1.y, vetexC.z - vetex1.z));
                    // length_v1 = length_v1 > 1e-4 ? length_v1 : 1e-4;
                    length_v1 = length_v1 > avgMeshEdgeLength*anisotropy ? length_v1 : avgMeshEdgeLength*anisotropy;
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
                const glm::vec3 v1n = glm::vec3((vetexC.x - vetex1.x)*length_v1_, (vetexC.y - vetex1.y)*length_v1_, (vetexC.z - vetex1.z)*length_v1_);
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
                    const double a = weight * v2n[n] / length_v1;
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

void LayerOpt::OptimizeEdgeOrthogonality(std::vector<Trip>& A_Entries, std::vector<float>& b, size_t& row)
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
            const glm::vec3 v2 = glm::vec3(vetexC.x - vetex2.x, vetexC.y - vetex2.y, vetexC.z - vetex2.z);
            const glm::vec3 v2n = glm::normalize(v2);
            double length_v1 = edge1.length;
            if (!useAverageTargetLength) {
                length_v1 = glm::length(glm::vec3(vetexC.x - vetex1.x, vetexC.y - vetex1.y, vetexC.z - vetex1.z));
                // length_v1 = length_v1 > 1e-4 ? length_v1 : 1e-4;
                length_v1 = length_v1 > avgMeshEdgeLength*anisotropy ? length_v1 : avgMeshEdgeLength*anisotropy;
            }
            ///////////////////////////////////////////////
            const double length_v1_ = 1.0/length_v1;
            const glm::vec3 v1n = glm::vec3((vetexC.x - vetex1.x)*length_v1_, (vetexC.y - vetex1.y)*length_v1_, (vetexC.z - vetex1.z)*length_v1_);
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
                const double a = weight * v2n[n] / length_v1;
                A_Entries.push_back(Trip(row, 3 * vetexC.id + n, a));
                A_Entries.push_back(Trip(row, 3 * vetex1.id + n, -a));
            }
            b.push_back(0.0 * weight);
            row++;
        }
    }
    EOrthogonality.push_back(energy);
}

void LayerOpt::OptimizeEdgeStraightness(std::vector<Trip>& A_Entries, std::vector<float>& b, size_t& row)
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
            const glm::vec3 v2 = glm::vec3(vetexC.x - vetex2.x, vetexC.y - vetex2.y, vetexC.z - vetex2.z);
            const glm::vec3 v2n = glm::normalize(v2);
            double length_v1 = edge1.length;
            if (!useAverageTargetLength) {
                length_v1 = glm::length(glm::vec3(vetexC.x - vetex1.x, vetexC.y - vetex1.y, vetexC.z - vetex1.z));
                // length_v1 = length_v1 > 1e-4 ? length_v1 : 1e-4;
                length_v1 = length_v1 > avgMeshEdgeLength*anisotropy ? length_v1 : avgMeshEdgeLength*anisotropy;
            }
            ///////////////////////////////////////////////
            const double length_v1_ = 1.0/length_v1;
            const glm::vec3 v1n = glm::vec3((vetexC.x - vetex1.x)*length_v1_, (vetexC.y - vetex1.y)*length_v1_, (vetexC.z - vetex1.z)*length_v1_);
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
                const double a = weight * v2n[n] / length_v1;
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

void LayerOpt::OptimizeBoundaryVertices(std::vector<Trip>& A_Entries, std::vector<float>& b, size_t& row)
{
    for (size_t k = 0; k < mesh.V.size(); k++){
        const Vertex& v = mesh.V.at(k);
        if (v.isBoundary){
            for (size_t n = 0; n < 3; n++)
            {
                A_Entries.push_back(Trip(row, 3 * v.id + n, alpha));
                b.push_back(alpha * v[n]);
                row++;
            }
        }
    }
}

void LayerOpt::OptimizeInnerVertices(std::vector<Trip>& A_Entries, std::vector<float>& b, size_t& row)
{
    /*
     *      3
     *     /|
     *    / |
     *   /| |
     * 0/ | |
     * |  |/|
     * | f*-|-------->* v(innerV)
     * |/ | |
     * |  |/ 2
     * |  /
     * | /
     * |/
     * 1
     * */
    MeshFileReader reader("hex.vtk");
    const Mesh& hex = reader.GetMesh();

    //std::vector<size_t> frameIds;
    for (size_t k = 0; k < mesh.V.size(); k++) {
        const Vertex& v = mesh.V.at(k);
        if (v.isBoundary) {
            for (size_t j = 0; j < v.N_Vids.size(); j++) {
                Vertex& innerV = mesh.V.at(v.N_Vids.at(j));
                if (innerV.isBoundary) {
                    //frameIds.push_back(innerFrame.id);
                    for (size_t n = 0; n < 3; n++)
                    {
                        A_Entries.push_back(Trip(row, 3 * innerV.id + n, alpha));
                        b.push_back(alpha * v[n]);
                        row++;
                    }
                }
            }
        }
    }
    //framefield.WriteFile("innerframe.vtk", frameIds);
}

void LayerOpt::ComputeMeshTargetLength()
{
    std::cout << "=============================\n";
    std::cout << "   Computing TargetLength    \n";
//    double sumEdgeLength = 0.0;
//    size_t numOfBoundaryEdges = 0;
//    for (size_t i = 0; i < mesh.E.size(); i++)
//    {
//        if (mesh.E.at(i).isBoundary){
//            const Vertex& v1 = mesh.V.at(mesh.E.at(i).Vids.at(0));
//            const Vertex& v2 = mesh.V.at(mesh.E.at(i).Vids.at(1));
//            mesh.E[i].length = glm::length(v1 - v2);
//            sumEdgeLength += mesh.E[i].length;
//            //sumEdgeLength += glm::length(v1 - v2);
//            numOfBoundaryEdges++;
//        }
//    }
//    double avgMeshEdgeLength = sumEdgeLength / numOfBoundaryEdges;
//    std::cout << "Average Surface Edge Length = " << avgMeshEdgeLength << std::endl;

    std::vector<Trip> coefficients;
    std::vector<float> b;
    size_t row  = 0;
    size_t numOfparallelEids = 0;
    size_t numOfconsecutiveEids = 0;
    for (size_t i = 0; i < mesh.E.size(); i++)
    {
        coefficients.push_back(Trip(row, i, 1));
        if (mesh.E[i].isBoundary)
            b.push_back(mesh.E[i].length);
//        else if (mesh.E[i].isSingularity)
//            b.push_back(0.5 * avgMeshEdgeLength);
        else
            b.push_back(avgMeshEdgeLength);
        row++;

        for (size_t j = 0; j < mesh.E[i].parallelEids.size(); j++)
        {
            size_t N_Eid = mesh.E[i].parallelEids[j];
            coefficients.push_back(Trip(row, i, 1));
            coefficients.push_back(Trip(row, N_Eid, -1));
            b.push_back(0);
            numOfparallelEids++;
            row++;
        }
        for (size_t j = 0; j < mesh.E[i].consecutiveEids.size(); j++)
        {
            size_t N_Eid = mesh.E[i].consecutiveEids[j];
            coefficients.push_back(Trip(row, i, 1));
            coefficients.push_back(Trip(row, N_Eid, -1));
            b.push_back(0);
            numOfconsecutiveEids++;
            row++;
        }
    }

    size_t col = mesh.E.size();
//    std::cout << "Edge Length Row = " << row << std::endl;
//    std::cout << "Edge Length Col = " << col << std::endl;
//    std::cout << "numOfparallelEids = " << numOfparallelEids << std::endl;
//    std::cout << "numOfconsecutiveEids = " << numOfconsecutiveEids << std::endl;

    Eigen::SparseMatrix<float> A(row, col);
    A.setFromTriplets(coefficients.begin(), coefficients.end());
    Eigen::SparseMatrix<float> A_T = A.transpose();
    Eigen::SparseMatrix<float> ATA = A_T * A;
    //  MatrixXf ;
    VectorXf X(col), B(row);
    for (size_t i = 0; i < row; i++)
        B(i) = b[i];

    VectorXf ATB = A_T * B;
    Eigen::SparseLU<SpMat> chol(ATA);
    X = chol.solve(ATB);
    for (size_t i = 0; i < mesh.E.size(); i++)
        mesh.E[i].length = X[i];

    std::cout << "=============================\n";
}

void LayerOpt::SetAlpha(const double value/* = 100.0*/)
{
    alpha = value;
}

void LayerOpt::SetBeta(const double value/* = 100.0*/)
{
    beta = value;
}

void LayerOpt::SetGamma(const double value/* = 100.0*/)
{
    gamma = value;
}

void LayerOpt::SetStepSize(const double value/* = 100.0*/)
{
    stepSize = value;
}

void MeshOpt::SetAnisotropy(const double value/* = 100.0*/)
{
    anisotropy = value;
}

void LayerOpt::SetUseAverageTargetLength(bool value)
{
    useAverageTargetLength = value;
}

void LayerOpt::SetRecoverable(bool value)
{
    recoverable = value;
}

void LayerOpt::OutputBadCells(const std::vector<size_t>& badCellIds, const char* filename)
{
    std::vector<Cell> cells(badCellIds.size());
    for (size_t i = 0; i < badCellIds.size(); i++)
        cells.at(i) = mesh.C.at(badCellIds.at(i));
    MeshFileWriter writer(mesh.V, cells, filename);
    writer.WriteFile();
}
