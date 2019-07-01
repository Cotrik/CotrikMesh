/*
 * MeshASJOpt.cpp
 *
 *  Created on: March 23, 2017
 *      Author: cotrik
 */

#include "MeshASJOpt.h"
#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "MeshQuality.h"

#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <queue>
#include <iostream>
#include <math.h>

MeshASJOpt::MeshASJOpt(const Mesh& mesh)
: MeshOpt(mesh)
{
    // TODO Auto-generated constructor stub

}

MeshASJOpt::~MeshASJOpt()
{
    // TODO Auto-generated destructor stub
}

size_t MeshASJOpt::Run(const size_t iters/* = 1*/)
{
    //ConvertSurfaceToTriangleMesh();
    Eigen::initParallel();
    if (stepSize < 0.1){
        std::cout << "ERROR! stepSize < 0.1\n";
        return 0;
    }
    double sumEdgeLength = 0.0;
    size_t numOfBoundaryEdges = 0;
    for (size_t i = 0; i < mesh.E.size(); i++) {
        if (mesh.E.at(i).isBoundary) {
            const Vertex& v1 = mesh.V.at(mesh.E.at(i).Vids.at(0));
            const Vertex& v2 = mesh.V.at(mesh.E.at(i).Vids.at(1));
            mesh.E[i].length = glm::length(v1 - v2);
            sumEdgeLength += mesh.E[i].length;
            numOfBoundaryEdges++;
        }
    }
    avgMeshEdgeLength = sumEdgeLength / numOfBoundaryEdges;
    std::cout << "Average Surface Edge Length = " << avgMeshEdgeLength << std::endl;
    //----------------------------------------------
    this->minScaledJacobian = 0.0;
    m_numOfInvertdElements = GetScaledJacobianVerdict(mesh, this->currentMSJ, this->currentASJ, this->minScaledJacobian);
    this->minScaledJacobian = this->currentMSJ;
    std::cout <<  "initial MSJ = " << this->currentMSJ << "initial ASJ = " << this->currentASJ << std::endl;
    //----------------------------------------------
    int iter = 0;
    bool converged = false;
    double initStepSize = stepSize;
    double lastEnergy = 10000000000;
    ComputeMeshTargetLength();

    while (!converged && iter++ < iters)
    {
        double energy = 0;
        Optimize(energy);
        double rate = fabs(lastEnergy - energy)/lastEnergy;
        converged = rate < 1e-6;
        lastEnergy = energy;
        stepSize *= initStepSize;
        double minimumScaledJacobian = 0.0;
        double averageScaledJacobian = 0.0;
        m_numOfInvertdElements = GetScaledJacobianVerdict(mesh, minimumScaledJacobian, averageScaledJacobian, this->minScaledJacobian);
        this->minScaledJacobian = minimumScaledJacobian;
        std::cout << "iter = " << iter << " MSJ = " << minimumScaledJacobian << " ASJ = " << averageScaledJacobian << " converged = " << rate <<std::endl;
        if (converged) {
            std::cout << "*************************" << std::endl;
            std::cout << "Converged at iter " << iter << std::endl;
            std::cout << "*************************" << std::endl;
            break;
        }
    }
//    MeshFileWriter writer(mesh, "opt.vtk");
//    writer.WriteFile();
    return iter;
}

bool MeshASJOpt::Optimize(double& energy)
{
    std::vector<Trip> A_Entries;
    std::vector<float> b;
    size_t row = 0;
    size_t currentRow = row;
    size_t currentNumOfA_Entries = A_Entries.size();

    currentRow = row;
    currentNumOfA_Entries = A_Entries.size();
    OptimizeSingularity(A_Entries, b, row);

    currentRow = row;
    currentNumOfA_Entries = A_Entries.size();
    OptimizeEdgeOrthogonality(A_Entries, b, row);

    currentRow = row;
    currentNumOfA_Entries = A_Entries.size();
    OptimizeEdgeStraightness(A_Entries, b, row);

    currentRow = row;
    currentNumOfA_Entries = A_Entries.size();
    size_t a_id = 0;
    OptimizeBoundaryVertices(A_Entries, b, row);

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

    energy = 0;
    if (recoverable) {
//#pragma omp parallel for
        for (size_t i = 0; i < mesh.V.size(); i++) {
            Vertex& v = mesh.V.at(i);
            if (v.isBoundary) continue;
            const double& x = X[3 * i + 0];
            const double& y = X[3 * i + 1];
            const double& z = X[3 * i + 2];
            const glm::dvec3 oldv = v.xyz();
            double oldMinSJ = 0;
            double oldAvgSJ = 0;
            mesh.GetQuality(v, oldMinSJ, oldAvgSJ);
            double newMinSJ = 0;
            double newAvgSJ = 0;
            const glm::dvec3 newv(x,y,z);
            v = newv;
            mesh.GetQuality(v, newMinSJ, newAvgSJ);
            //if (newMinSJ < this->minScaledJacobian || newAvgSJ < oldAvgSJ){
            if (newMinSJ < oldMinSJ){
                v = oldv;
            } else {
                double l = glm::length(newv - oldv);
                energy += l*l;
            }
        }
    }
    return energy;
}
