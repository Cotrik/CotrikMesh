/*
 * MeshOpt.cpp
 *
 *  Created on: Nov 25, 2016
 *      Author: cotrik
 */

#include "MeshOpt.h"
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

//MeshOpt::MeshOpt()
MeshOpt::MeshOpt(const Mesh& mesh)
: mesh((Mesh&)mesh)
//, pTriMesh(NULL)
, alpha(1000000.0)
, beta(1000000.0)
, gamma(1.0)
, anisotropy(0.05)
, minScaledJacobian(0.0)
, avgMeshEdgeLength(1.0)
, stepSize(1.0)
, useAverageTargetLength(false)
, recoverable(true)
, allowBigStep(false)
, changeBoundary(false)
, m_numOfInvertdElements(MAXID)
{
    // TODO Auto-generated constructor stub

}

MeshOpt::~MeshOpt()
{
    // TODO Auto-generated destructor stub
}

void MeshOpt::ComputeAvgEdgeLength()
{
    //ConvertSurfaceToTriangleMesh();
    Eigen::initParallel();
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
}

size_t MeshOpt::Run(const size_t iters/* = 1*/)
{
    ComputeAvgEdgeLength();
    //----------------------------------------------
//    MeshFileWriter writerT(mesh, "temp.vtk");
//    writerT.WriteFile();
//    double minimumScaledJacobianT = 0.0;
//    double averageScaledJacobianT = 0.0;
//    double maximumScaledJacobianT = 0.0;
//    std::vector<size_t> badCellIdsT;
//    m_numOfInvertdElements = GetQuality("temp.vtk", minimumScaledJacobianT, averageScaledJacobianT, maximumScaledJacobianT,
//                                        badCellIdsT, true, this->minScaledJacobian);
//    m_numOfInvertdElements = GetMinScaledJacobian(mesh, minimumScaledJacobianT, badCellIdsT, this->minScaledJacobian);
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
    //if (useAverageTargetLength)
        ComputeMeshTargetLength();
        static int global_count = 0;
        global_count++;
    while (!converged && iter++ < iters)
    {
//        std::string filename = std::string("EdgesEnergy.") + std::to_string(iter) + ".vtk";
//        MeshFileWriter edgeWriter(mesh, filename.c_str());
//        edgeWriter.WriteEdgesVtk();

        if (!initUseAverageTargetLength  && iter == 1)
            useAverageTargetLength = true;
//        std::cout << "stepSize = " << stepSize << std::endl;
//        std::cout << "beta = " << beta << std::endl;
//        if (useAverageTargetLength)
//            ComputeMeshTargetLength();
        /////////////////////////////////////////////
        mesh.ClearLabelOfSurface();
        mesh.LabelSurface();
        std::string filename = std::string("Faces.") + std::to_string(iter) + ".vtk";
//        MeshFileWriter facesFileWriter(mesh, filename.c_str());
//        facesFileWriter.WriteFacesVtk();
        /////////////////////////////////////////////
        converged = Optimize();
        if (!initUseAverageTargetLength && iter == 1)
            useAverageTargetLength = false;
        stepSize *= initStepSize;
        //beta *= initStepSize;

//        filename = std::string("MeshOpt.") + std::to_string(iter) + ".vtk";
//        MeshFileWriter writer(mesh, filename.c_str());
//        writer.WriteFile();

//        double minimumScaledJacobian = 0.0;
//        double averageScaledJacobian = 0.0;
//        double maximumScaledJacobian = 0.0;
//        std::vector<size_t> badCellIds;
//        m_numOfInvertdElements = GetQuality(filename.c_str(), minimumScaledJacobian, averageScaledJacobian, maximumScaledJacobian,
//                                            badCellIds, true, this->minScaledJacobian);
//        filename = std::to_string(global_count) + std::string("BadCells.") + std::to_string(iter) + ".vtk";
//        OutputBadCells(badCellIds, filename.c_str());
//        filename = std::to_string(global_count) + std::string("BadCells.ex.") + std::to_string(iter) + ".vtk";
//        OutputBadCellsAndExtendCells(badCellIds, filename.c_str());
//        filename = std::string("EdgeOfBadCells.") + std::to_string(iter) + ".vtk";
//        OutputEdgesOfBadCells(badCellIds, filename.c_str());
//        double minimumScaledJacobian = 0.0; std::vector<size_t> badCellIds;
//        m_numOfInvertdElements = GetMinScaledJacobian(mesh, minimumScaledJacobian, badCellIds, this->minScaledJacobian);

        m_numOfInvertdElements = GetMinScaledJacobianVerdict(mesh, minimumScaledJacobian, this->minScaledJacobian);
        std::cout << "iter = " << iter << " #inverted = " << m_numOfInvertdElements << " MSJ = " << minimumScaledJacobian << std::endl;
        if (prevMinimumScaledJacobian > this->minScaledJacobian && minimumScaledJacobian <= prevMinimumScaledJacobian)
        {
            filename = std::string("MeshOpt.") + std::to_string(iter - 1) + ".vtk";
            std::cout << "*************************" << std::endl;
            std::cout << "\033[1;32mBest Mesh is " << filename << "\033[0m" << std::endl;
            std::cout << "*************************" << std::endl;
//            MeshFileReader reader(filename.c_str());
//            const Mesh& bestmesh = reader.GetMesh();
//            MeshFileWriter optwriter(bestmesh, "opt.vtk");
//            optwriter.WriteFile();
            untangled = true;
            break;
        }
        else if (converged)
        {
            std::cout << "*************************" << std::endl;
            std::cout << "Converged at " << filename << std::endl;
            std::cout << "*************************" << std::endl;
        }
        prevMinimumScaledJacobian = minimumScaledJacobian;
        bestmesh = mesh;
    }

//    std::vector<double> E_total(ESingularity.size());
//    for (size_t i = 0; i < ESingularity.size(); i++)
//        E_total.at(i) = ESingularity.at(i) + EOrthogonality.at(i) + EStraightness.at(i);
//    std::cout << "*************************" << std::endl;
//    std::cout << "E_total";
//    for (size_t i = 0; i < E_total.size(); i++)
//        std::cout << "\t" << E_total.at(i);
//    std::cout << "\nESingularity";
//    for (size_t i = 0; i < ESingularity.size(); i++)
//        std::cout << "\t" << ESingularity.at(i);
//    std::cout << "\nEOrthogonality";
//    for (size_t i = 0; i < EOrthogonality.size(); i++)
//        std::cout << "\t" << EOrthogonality.at(i);
//    std::cout << "\nEStraightness";
//    for (size_t i = 0; i < EStraightness.size(); i++)
//        std::cout << "\t" << EStraightness.at(i);
//    std::cout << "\n*************************\n" << std::endl;

//    if (untangled)
        return iter - 1;
//    return iter;
}

bool MeshOpt::Optimize()
{
    //std::cout << "Constraint-------\t#Entries\t #row\n";
    //std::cout << "-----------------\t---------------------------\t\n";
    std::vector<Trip> A_Entries;
    std::vector<float> b;
    size_t row = 0;
    size_t currentRow = row;
    size_t currentNumOfA_Entries = A_Entries.size();
    //OptimizeBoundaryVertices(A_Entries, b, row);
    //std::cout << "BoundaryVertices-\t" << A_Entries.size() - currentNumOfA_Entries << "\t -------- " << row - currentRow << std::endl;

    //currentRow = row;
    //currentNumOfA_Entries = A_Entries.size();
    //OptimizeFacePlanarity(A_Entries, b, row);
    //std::cout << "FacePlanarity----\t" << A_Entries.size() - currentNumOfA_Entries << "\t -------- " << row - currentRow << std::endl;

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

    //currentRow = row;
    //currentNumOfA_Entries = A_Entries.size();
    //OptimizeSmoothness(A_Entries, b, row);
    //std::cout << "Smoothness------\t" << A_Entries.size() - currentNumOfA_Entries << "\t -------- " << row - currentRow << std::endl;

    currentRow = row;
    currentNumOfA_Entries = A_Entries.size();
    size_t a_id = OptimizeSurfaceVertices(A_Entries, b, row);
    //std::cout << "SurfaceVertices--\t" << A_Entries.size() - currentNumOfA_Entries << "\t -------- " << row - currentRow << std::endl;

    const size_t col = 3 * mesh.V.size() + a_id;

//    std::cout << "----------- A Info ------------\n";
//    std::cout << "Entries = " << A_Entries.size() << " row = " << row << " col = " << col << std::endl;

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

    // ------- Output File -----------
//    std::ofstream ofs2("B.txt");
//    ofs2 << B;
//    std::ofstream ofs1("ATB.txt");
//    ofs1 << ATB;
//    std::ofstream ofs("jX.txt");
//    ofs << X;
    // ------- ----------- -----------

    bool converged = true;
    if (recoverable) {
        Mesh targetMesh(mesh);
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
            if (!allowBigStep)
            if (fabs(double(mesh.V[i].x - X[3 * i + 0])) > 2.0 * avgMeshEdgeLength
                || fabs(double(mesh.V[i].y - X[3 * i + 1])) > 2.0 * avgMeshEdgeLength
                || fabs(double(mesh.V[i].z - X[3 * i + 2])) > 2.0 * avgMeshEdgeLength)
            {
                continue;
            }
//            if (mesh.V[i].isCorner)
//                continue;
            mesh.V[i].x = stepSize * X[3 * i + 0] + (1.0 - stepSize) * mesh.V[i].x;
            mesh.V[i].y = stepSize * X[3 * i + 1] + (1.0 - stepSize) * mesh.V[i].y;
            mesh.V[i].z = stepSize * X[3 * i + 2] + (1.0 - stepSize) * mesh.V[i].z;
//            if (std::isnan(mesh.V[i].x) || std::isnan(mesh.V[i].y) || std::isnan(mesh.V[i].z))
//            {
//                std::cout << "Recover previous mesh\n";
//                for (size_t i = 0; i < mesh.V.size(); i++)
//                {
//                    mesh.V[i].x = oldV[i].x;
//                    mesh.V[i].y = oldV[i].y;
//                    mesh.V[i].z = oldV[i].z;
//                }
//                break;
//            }
        }
        if (useProjection) {
        if (!projectToTargetSurface)
            mesh.ProjectTo(targetMesh);
        else
            mesh.FastProjectTo(targetMesh);
        }
            //mesh.ProjectToTargetSurface(*m_targetSurfaceMesh, *m_targetSurfaceMesh);

//        MeshFileWriter writer(mesh, "temp.vtk");
//        writer.WriteFile();
//        double minimumScaledJacobian = 0.0;
//        double averageScaledJacobian = 0.0;
//        double maximumScaledJacobian = 0.0;
//        std::vector<size_t> badCellIds1;
//        size_t InvertedElements = GetQuality("temp.vtk", minimumScaledJacobian, averageScaledJacobian, maximumScaledJacobian,
//                                             badCellIds1, true, this->minScaledJacobian);

        double minimumScaledJacobian = 0.0;
        size_t InvertedElements = GetMinScaledJacobianVerdict(mesh, minimumScaledJacobian, this->minScaledJacobian);

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
            if (!allowBigStep)
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

//bool MeshOpt::Optimize()
//{
//    std::cout << "Constraint-------\t#Entries\t #row\n";
//    std::cout << "-----------------\t---------------------------\t\n";
//    std::vector<Trip> A_Entries;
//    std::vector<float> b;
//    size_t row = 0;
//    size_t currentRow = row;
//    size_t currentNumOfA_Entries = A_Entries.size();
//    OptimizeBoundaryVertices(A_Entries, b, row);
//    std::cout << "BoundaryVertices-\t" << A_Entries.size() - currentNumOfA_Entries << "\t -------- " << row - currentRow << std::endl;
//
//    currentRow = row;
//    currentNumOfA_Entries = A_Entries.size();
//    OptimizeScaledJacobian(A_Entries, b, row);
//    std::cout << "ScaledJacobian---\t" << A_Entries.size() - currentNumOfA_Entries << "\t -------- " << row - currentRow << std::endl;
//
////    currentRow = row;
////    currentNumOfA_Entries = A_Entries.size();
////    OptimizeFacePlanarity(A_Entries, b, row);
////    std::cout << "FacePlanarity----\t" << A_Entries.size() - currentNumOfA_Entries << "\t -------- " << row - currentRow << std::endl;
////
////    currentRow = row;
////    currentNumOfA_Entries = A_Entries.size();
////    OptimizeSingularity(A_Entries, b, row);
////    std::cout << "EdgeSingularity--\t" << A_Entries.size() - currentNumOfA_Entries << "\t -------- " << row - currentRow << std::endl;
////
////    currentRow = row;
////    currentNumOfA_Entries = A_Entries.size();
////    OptimizeEdgeOrthogonality(A_Entries, b, row);
////    std::cout << "EdgeOrthogonality\t" << A_Entries.size() - currentNumOfA_Entries << "\t -------- " << row - currentRow << std::endl;
////
////    currentRow = row;
////    currentNumOfA_Entries = A_Entries.size();
////    OptimizeEdgeStraightness(A_Entries, b, row);
////    std::cout << "EdgeStraightness\t" << A_Entries.size() - currentNumOfA_Entries << "\t -------- " << row - currentRow << std::endl;
//    const size_t col = 3 * mesh.V.size();
//
//    std::cout << "----------- A Info ------------\n";
//    std::cout << "Entries = " << A_Entries.size() << " row = " << row << " col = " << col << std::endl;
//
//    Eigen::SparseMatrix<float> A(row, col);
//    A.setFromTriplets(A_Entries.begin(), A_Entries.end());
//    Eigen::SparseMatrix<float> AT = A.transpose();
//    Eigen::SparseMatrix<float> ATA = AT * A;
//
//    //  MatrixXf ;
//    VectorXf B(row);
//    for (size_t i = 0; i < row; i++)
//        B(i) = b[i];
//
//    VectorXf ATB = AT * B;
//    Eigen::SimplicialLDLT<SpMat> chol(ATA);
//    VectorXf X(col);
//    X = chol.solve(ATB);
//
//    bool converged = true;
//    if (recoverable) {
//        std::vector<glm::vec3> oldV(mesh.V.size());
//        for (size_t i = 0; i < mesh.V.size(); i++)
//        {
//            oldV[i].x = mesh.V[i].x;
//            oldV[i].y = mesh.V[i].y;
//            oldV[i].z = mesh.V[i].z;
//        }
//
//        for (size_t i = 0; i < mesh.V.size(); i++)
//        {
//            if (fabs(double(mesh.V[i].x - X[3 * i + 0])) > 1e-6
//                || fabs(double(mesh.V[i].y - X[3 * i + 1])) > 1e-6
//                || fabs(double(mesh.V[i].z - X[3 * i + 2])) > 1e-6)
//            {
//                converged = false;
//            }
//            if (!allowBigStep)
//            if (fabs(double(mesh.V[i].x - X[3 * i + 0])) > 2.0 * avgMeshEdgeLength
//                || fabs(double(mesh.V[i].y - X[3 * i + 1])) > 2.0 * avgMeshEdgeLength
//                || fabs(double(mesh.V[i].z - X[3 * i + 2])) > 2.0 * avgMeshEdgeLength)
//            {
//                continue;
//            }
//            mesh.V[i].x = stepSize * X[3 * i + 0] + (1.0 - stepSize) * mesh.V[i].x;
//            mesh.V[i].y = stepSize * X[3 * i + 1] + (1.0 - stepSize) * mesh.V[i].y;
//            mesh.V[i].z = stepSize * X[3 * i + 2] + (1.0 - stepSize) * mesh.V[i].z;
//        }
//
//        MeshFileWriter writer(mesh, "temp.vtk");
//        writer.WriteFile();
//        double minimumScaledJacobian = 0.0;
//        double averageScaledJacobian = 0.0;
//        double maximumScaledJacobian = 0.0;
//        std::vector<size_t> badCellIds1;
//        size_t InvertedElements = GetQuality("temp.vtk", minimumScaledJacobian, averageScaledJacobian, maximumScaledJacobian, badCellIds1);
//
//        if (InvertedElements > m_numOfInvertdElements){
//            std::cout << "Recover previous mesh\n";
//            for (size_t i = 0; i < mesh.V.size(); i++)
//            {
//                mesh.V[i].x = oldV[i].x;
//                mesh.V[i].y = oldV[i].y;
//                mesh.V[i].z = oldV[i].z;
//            }
//        }
//    }
//    else {
//        for (size_t i = 0; i < mesh.V.size(); i++)
//        {
//            if (fabs(double(mesh.V[i].x - X[3 * i + 0])) > 1e-6
//                || fabs(double(mesh.V[i].y - X[3 * i + 1])) > 1e-6
//                || fabs(double(mesh.V[i].z - X[3 * i + 2])) > 1e-6)
//            {
//                converged = false;
//            }
//            if (!allowBigStep)
//            if (fabs(double(mesh.V[i].x - X[3 * i + 0])) > 2.0 * avgMeshEdgeLength
//                || fabs(double(mesh.V[i].y - X[3 * i + 1])) > 2.0 * avgMeshEdgeLength
//                || fabs(double(mesh.V[i].z - X[3 * i + 2])) > 2.0 * avgMeshEdgeLength)
//            {
//                continue;
//            }
//            mesh.V[i].x = stepSize * X[3 * i + 0] + (1.0 - stepSize) * mesh.V[i].x;
//            mesh.V[i].y = stepSize * X[3 * i + 1] + (1.0 - stepSize) * mesh.V[i].y;
//            mesh.V[i].z = stepSize * X[3 * i + 2] + (1.0 - stepSize) * mesh.V[i].z;
//        }
//    }
//    return converged;
//}

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
void MeshOpt::OptimizeEdgeComformalty(std::vector<Trip>& A_Entries, std::vector<float>& b, size_t& row)
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

void MeshOpt::OptimizeSingularity(std::vector<Trip>& A_Entries, std::vector<float>& b, size_t& row)
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

void MeshOpt::OptimizeFacePlanarity(std::vector<Trip>& A_Entries, std::vector<float>& B, size_t& row)
{
//    for (size_t k = 0; k < mesh.E.size(); k++){
//        Edge& edge = mesh.E.at(k);
//        edge.energyOrthogonality = 0;
//    }
//    double energy = 0.0;
    double weight = beta;
    for (size_t k = 0; k < mesh.F.size(); k++) {
        const Face& face = mesh.F.at(k);
        if (!changeBoundary && face.isBoundary)
            continue;
        for (size_t j = 0; j < face.Vids.size(); j++) {
            const Vertex& v1 = mesh.V.at(face.Vids.at((j + 0)%4));  // (x1, y1, z1)
            const Vertex& v2 = mesh.V.at(face.Vids.at((j + 1)%4));  // (x2, y2, z2)
            const Vertex& v3 = mesh.V.at(face.Vids.at((j + 2)%4));  // (x3, y3, z3)
            const Vertex& v4 = mesh.V.at(face.Vids.at((j + 3)%4));  // (x4, y4, z4)
            const Vertex& v = v4;                        // (x , y , z )
            //if (v1.isBoundary || v2.isBoundary || v3.isBoundary || v4.isBoundary)
            if (!changeBoundary && v.isBoundary)
                continue;
            /*
             * | x -x1   y -y1   z -z1 |
             * | x2-x1   y2-y1   z2-z1 | = 0
             * | x3-x1   y3-y1   z3-z1 |
             * ------------------------------------------------------
             * (x - x1)((y2 - y1)(z3 - z1) - (y3 - y1)(z2 - z1)) +
             * (y - y1)((z2 - z1)(x3 - x1) - (z3 - z1)(x2 - x1)) +
             * (z - z1)((x2 - x1)(y3 - y1) - (x3 - x1)(y2 - y1)) = 0
             * ------------------------------------------------------
             * a = (y2 - y1)(z3 - z1) - (y3 - y1)(z2 - z1)
             * b = (z2 - z1)(x3 - x1) - (z3 - z1)(x2 - x1)
             * c = (x2 - x1)(y3 - y1) - (x3 - x1)(y2 - y1)
             * ------------------------------------------------------
             * a(x - x1) + b(y - y1) + c(z - z1) = 0
             * ax + by + cz = ax1 + by1 + cz1
             * */
            const double a = (v2.y - v1.y)*(v3.z - v1.z) - (v3.y - v1.z)*(v2.z - v1.z);
            const double b = (v2.z - v1.z)*(v3.x - v1.x) - (v3.z - v1.x)*(v2.z - v1.x);
            const double c = (v2.x - v1.x)*(v3.y - v1.y) - (v3.x - v1.y)*(v2.z - v1.y);
            A_Entries.push_back(Trip(row, 3 * v.id + 0, weight * a)); // ax
            A_Entries.push_back(Trip(row, 3 * v.id + 1, weight * b)); // by
            A_Entries.push_back(Trip(row, 3 * v.id + 2, weight * c)); // cz
            B.push_back((a*v1.x + b*v1.y + c*v1.z) * weight);
            row++;
        }
    }
//    EOrthogonality.push_back(energy);
}

void MeshOpt::OptimizeScaledJacobian(std::vector<Trip>& A_Entries, std::vector<float>& B, size_t& row)
{
    double weight = beta;
    for (size_t k = 0; k < mesh.C.size(); k++) {
        const Cell& cell = mesh.C.at(k);
        for (size_t j = 0; j < cell.Vids.size(); j++) {
            const Vertex& v1 = mesh.V.at(HexPoint_Points[j][0]);    // (x1, y1, z1)
            const Vertex& v2 = mesh.V.at(HexPoint_Points[j][1]);    // (x2, y2, z2)
            const Vertex& v3 = mesh.V.at(HexPoint_Points[j][2]);    // (x3, y3, z3)
            const Vertex& v4 = mesh.V.at(cell.Vids.at(j));          // (x4, y4, z4)
            const Vertex& v = v4;                                   // (x , y , z )

            const Edge& e1 = mesh.E.at(cell.Eids[HexPoint_Edges[j][0]]);
            const Edge& e2 = mesh.E.at(cell.Eids[HexPoint_Edges[j][1]]);
            const Edge& e3 = mesh.E.at(cell.Eids[HexPoint_Edges[j][2]]);

            if (!changeBoundary && v.isBoundary)
                continue;
            /*
             * | x1-x   y1-y   z1-z |         1        1        1
             * | x2-x   y2-y   z2-z | = 1 * ------ * ------ * ------
             * | x3-x   y3-y   z3-z |       ||e1||   ||e2||   ||e3||
             * ------------------------------------------------------
             * (x1 - x)((y2 - y)(z3 - z) - (y3 - y)(z2 - z)) +
             * (y1 - y)((z2 - z)(x3 - x) - (z3 - z)(x2 - x)) +
             * (z1 - z)((x2 - x)(y3 - y) - (x3 - x)(y2 - y))
             * ------------------------------------------------------
             * a = (−y2*z3+y1*z3+y3*z2−y1*z2−y3*z1+y2*z1) =  (y1-y2)z3 + (y3-y1)z2 + (y2-y3)z1
             * b = (+x2*z3−x1*z3−x3*z2+x1*z2+x3*z1−x2*z1) = -(x1-x2)z3 - (x3-x1)z2 - (x2-x3)z1
             * c = (−x2*y3+x1*y3+x3*y2−x1*y2−x3*y1+x2*y1) =  (x1-x2)y3 + (x3-x1)y2 + (x2-x3)y1
             * d = x1*y2*z3−x2*y1*z3−x1*y3*z2+x3*y1*z2+x2*y3*z1−x3*y2*z1 = (x1y2-x2y1)z3 - (x1y3-x3y1)z2 + (x2y3-x3y2)z1
             * ------------------------------------------------------
             *                          1        1        1
             * ax + by + cz = d + 1 * ------ * ------ * ------
             *                        ||e1||   ||e2||   ||e3||
             * */
            const double a =  (v1.y - v2.y) * v3.z + (v3.y - v1.y) * v2.z + (v2.y - v3.y) * v1.z;
            const double b = -(v1.x - v2.x) * v3.z - (v3.x - v1.x) * v2.z - (v2.x - v3.x) * v1.z;
            const double c =  (v1.x - v2.x) * v3.y + (v3.x - v1.x) * v2.y + (v2.x - v3.x) * v1.y;
//            const double d = v1.x * v2.y * v3.z − v2.x * v1.y * v3.z
//                           − v1.x * v3.y * v2.z + v3.x * v1.y * v2.z
//                           + v2.x * v3.y * v1.z − v1.x * v2.y * v1.z;
            const double d = (v1.x*v2.y-v2.x*v1.y)*v3.z - (v1.x*v3.y-v3.x*v1.y)*v2.z + (v2.x*v3.y-v3.x*v2.y)*v1.z;
            double length_v1 = e1.length;
            if (!useAverageTargetLength) {
                length_v1 = glm::length(glm::vec3(v1.x - v.x, v1.y - v.y, v1.z - v.z));
                length_v1 = length_v1 > avgMeshEdgeLength*anisotropy ? length_v1 : avgMeshEdgeLength*anisotropy;
                length_v1 = std::isnan(length_v1) ? avgMeshEdgeLength*anisotropy : length_v1;
            }

            double length_v2 = e2.length;
            if (!useAverageTargetLength) {
                length_v2 = glm::length(glm::vec3(v2.x - v.x, v2.y - v.y, v2.z - v.z));
                length_v2 = length_v2 > avgMeshEdgeLength*anisotropy ? length_v2 : avgMeshEdgeLength*anisotropy;
                length_v2 = std::isnan(length_v2) ? avgMeshEdgeLength*anisotropy : length_v2;
            }

            double length_v3 = e3.length;
            if (!useAverageTargetLength) {
                length_v3 = glm::length(glm::vec3(v3.x - v.x, v3.y - v.y, v3.z - v.z));
                length_v3 = length_v3 > avgMeshEdgeLength*anisotropy ? length_v3 : avgMeshEdgeLength*anisotropy;
                length_v3 = std::isnan(length_v3) ? avgMeshEdgeLength*anisotropy : length_v3;
            }

            A_Entries.push_back(Trip(row, 3 * v.id + 0, weight * a)); // ax
            A_Entries.push_back(Trip(row, 3 * v.id + 1, weight * b)); // by
            A_Entries.push_back(Trip(row, 3 * v.id + 2, weight * c)); // cz
            //B.push_back((d + 1.0/length_v1 + 1.0/length_v2 + 1.0/length_v3) * weight);
            B.push_back((d + 1.0/e1.length + 1.0/e2.length + 1.0/e3.length) * weight);
            row++;
        }
    }
}

void MeshOpt::OptimizeSmoothness(std::vector<Trip>& A_Entries, std::vector<float>& B, size_t& row)
{
    double weight = beta;
    for (size_t k = 0; k < mesh.V.size(); k++) {
        const Vertex& v = mesh.V.at(k);
        glm::vec3 avg(0.0, 0.0, 0.0);
        for (size_t j = 0; j < v.N_Vids.size(); j++) {
            const Vertex& n_v = mesh.V.at(v.N_Vids.at(j));
            avg.x += n_v.x;
            avg.y += n_v.y;
            avg.z += n_v.z;
        }
        avg.x /= v.N_Vids.size();
        avg.y /= v.N_Vids.size();
        avg.z /= v.N_Vids.size();
        for (size_t n = 0; n < 3; n++) {
            A_Entries.push_back(Trip(row, 3 * v.id + n, weight));
            B.push_back(weight * avg[n]);
            row++;
        }
    }
}

void MeshOpt::OptimizeEdgeOrthogonality(std::vector<Trip>& A_Entries, std::vector<float>& b, size_t& row)
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
                length_v1 = std::isnan(length_v1) ? avgMeshEdgeLength*anisotropy : length_v1;
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

void MeshOpt::OptimizeEdgeStraightness(std::vector<Trip>& A_Entries, std::vector<float>& b, size_t& row)
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
                length_v1 = std::isnan(length_v1) ? avgMeshEdgeLength*anisotropy : length_v1;
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

void MeshOpt::OptimizeBoundaryVertices(std::vector<Trip>& A_Entries, std::vector<float>& b, size_t& row)
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

float uctet(const glm::vec3& a, glm::vec3& b, const glm::vec3& c, glm::vec3& d)
{
    glm::vec3 x = b - a;
    glm::vec3 y = c - a;
    glm::vec3 z = d - a;
    float res = -((x[0] * y[1] * z[2] + x[1] * y[2] * z[0] + x[2] * y[0] * z[1]) - (x[2] * y[1] * z[0] + x[1] * y[0] * z[2] + x[0] * y[2] * z[1]));
    return res;
}

void orient_triangle_mesh_index(std::vector<Vertex>& Vs, std::vector<Cell>& Ts)
{
    std::vector<glm::vec3 > p(Vs.size());
    std::vector<std::vector<size_t> > f(Ts.size());
    for (int i = 0; i < Vs.size(); i++) {
        const Vertex& v = Vs.at(i);
        p[i].x = v.x;
        p[i].y = v.y;
        p[i].z = v.z;
    }
    for (int i = 0; i < Ts.size(); i++)
        f[i] = Ts[i].Vids;

    std::map<std::set<size_t>, std::vector<size_t> > edge_2_neb_tri;
    std::set<std::vector<size_t> > direct_edges;
    std::set<std::set<size_t> > no_direct_edges;
    std::vector<bool> whe_tri_in(f.size(), false);
    std::vector<std::vector<size_t> > nf = f;

    for (int i = 0; i < f.size(); i++) {
        std::vector<size_t> ct;
        ct.push_back(i);
        for (size_t j = 0; j < 3; j++) {
            std::set<size_t> x;
            x.insert(f[i][TriEdge[j][0]]);
            x.insert(f[i][TriEdge[j][1]]);
            if (edge_2_neb_tri.find(x) == edge_2_neb_tri.end())
                edge_2_neb_tri.insert(std::pair<std::set<size_t>, std::vector<size_t> >(x, ct));
            else
                edge_2_neb_tri[x].push_back(i);
        }
    }

    std::set<size_t> xa, xb, xc;
    std::vector<size_t> ya, yb, yc;

    xa.insert(f[0][0]);
    xa.insert(f[0][1]);

    xb.insert(f[0][1]);
    xb.insert(f[0][2]);

    xc.insert(f[0][2]);
    xc.insert(f[0][0]);

    ya.push_back(f[0][0]);
    ya.push_back(f[0][1]);

    yb.push_back(f[0][1]);
    yb.push_back(f[0][2]);

    yc.push_back(f[0][2]);
    yc.push_back(f[0][0]);

    no_direct_edges.insert(xa);
    no_direct_edges.insert(xb);
    no_direct_edges.insert(xc);

    direct_edges.insert(ya);
    direct_edges.insert(yb);
    direct_edges.insert(yc);

    whe_tri_in[0] = true;
    std::queue<size_t> queue_loop;

    for (int i = 0; i < 2; i++) {
        if (!whe_tri_in[edge_2_neb_tri[xa][i]]) {
            queue_loop.push(edge_2_neb_tri[xa][i]);
            whe_tri_in[edge_2_neb_tri[xa][i]] = true;
        }

        if (!whe_tri_in[edge_2_neb_tri[xb][i]]) {
            queue_loop.push(edge_2_neb_tri[xb][i]);

            whe_tri_in[edge_2_neb_tri[xb][i]] = true;
        }

        if (!whe_tri_in[edge_2_neb_tri[xc][i]]) {
            queue_loop.push(edge_2_neb_tri[xc][i]);

            whe_tri_in[edge_2_neb_tri[xc][i]] = true;
        }
    }

    while (!queue_loop.empty())
    {
        xa.clear();
        xb.clear();
        xc.clear();

        ya.clear();
        yb.clear();
        yc.clear();

        int c;

        c = queue_loop.front();

        xa.insert(f[c][0]);
        xa.insert(f[c][1]);

        xb.insert(f[c][1]);
        xb.insert(f[c][2]);

        xc.insert(f[c][2]);
        xc.insert(f[c][0]);

        ya.push_back(f[c][0]);
        ya.push_back(f[c][1]);

        yb.push_back(f[c][1]);
        yb.push_back(f[c][2]);

        yc.push_back(f[c][2]);
        yc.push_back(f[c][0]);

        int cnt, ct;

        cnt = 0;
        ct = 0;

        if (no_direct_edges.find(xa) != no_direct_edges.end())
            cnt++;

        if (no_direct_edges.find(xb) != no_direct_edges.end())
            cnt++;

        if (no_direct_edges.find(xc) != no_direct_edges.end())
            cnt++;

        if (direct_edges.find(ya) != direct_edges.end())
            ct++;

        if (direct_edges.find(yb) != direct_edges.end())
            ct++;

        if (direct_edges.find(yc) != direct_edges.end())
            ct++;

        if (cnt == 0) {
            std::cout << "Error in triangle direction solving!" << std::endl;
            exit(0);
        }

        if (ct != 0) {
            ya.clear();
            yb.clear();
            yc.clear();

            ya.push_back(f[c][1]);
            ya.push_back(f[c][0]);

            yb.push_back(f[c][2]);
            yb.push_back(f[c][1]);

            yc.push_back(f[c][0]);
            yc.push_back(f[c][2]);

            nf[c][0] = f[c][1];
            nf[c][1] = f[c][0];
        }

        no_direct_edges.insert(xa);
        no_direct_edges.insert(xb);
        no_direct_edges.insert(xc);

        direct_edges.insert(ya);
        direct_edges.insert(yb);
        direct_edges.insert(yc);

        for (size_t i = 0; i < 2; i++)
        {
            if (!whe_tri_in[edge_2_neb_tri[xa][i]]) {
                queue_loop.push(edge_2_neb_tri[xa][i]);
                whe_tri_in[edge_2_neb_tri[xa][i]] = true;
            }

            if (!whe_tri_in[edge_2_neb_tri[xb][i]]) {
                queue_loop.push(edge_2_neb_tri[xb][i]);
                whe_tri_in[edge_2_neb_tri[xb][i]] = true;
            }

            if (!whe_tri_in[edge_2_neb_tri[xc][i]]) {
                queue_loop.push(edge_2_neb_tri[xc][i]);
                whe_tri_in[edge_2_neb_tri[xc][i]] = true;
            }
        }

        queue_loop.pop();
    }

    float res = 0;
    glm::vec3 ori(0.0, 0.0, 0.0);
    for (size_t i = 0; i < nf.size(); i++)
        res += uctet(ori, p[nf[i][0]], p[nf[i][1]], p[nf[i][2]]);

    if (res > 0) {
        int tmi;
        for (int i = 0; i < nf.size(); i++) {
            tmi = nf[i][0];
            nf[i][0] = nf[i][1];
            nf[i][1] = tmi;
        }
    }

    for (int i = 0; i < Ts.size(); i++)
        Ts[i].Vids = nf[i];
}
void MeshOpt::ConvertSurfaceToTriangleMesh()
{
    std::vector<bool> tags(mesh.V.size(), false);
    std::vector<int> tag_indices(mesh.V.size(), -1);

    std::vector<Vertex> V;
    std::vector<Cell> T;   // triangles
    for (size_t i = 0; i < mesh.F.size(); i++)
    {
        const Face& face = mesh.F.at(i);
        if (face.isBoundary)
        {
            for (size_t j = 0; j < 4; j++)
                tags[face.Vids[j]] = true;

            Cell tri(3);
            tri.id = T.size();
            tri.Vids[0] = face.Vids[0];
            tri.Vids[1] = face.Vids[1];
            tri.Vids[2] = face.Vids[2];
            T.push_back(tri);

            tri.id = T.size();
            tri.Vids[0] = face.Vids[2];
            tri.Vids[1] = face.Vids[3];
            tri.Vids[2] = face.Vids[0];
            T.push_back(tri);
        }
    }
    for (size_t i = 0; i < tags.size(); i++)
    {
        if (tags[i])
        {
            Vertex v;
            v.id = V.size();
            v.x = mesh.V[i].x;
            v.y = mesh.V[i].y;
            v.z = mesh.V[i].z;

            v.hvid = i;
            //v.type=0;
            v.type = 2;

            mesh.V[i].triVid = v.id;   // Need Fixed
            V.push_back(v);
            tag_indices[i] = v.id;
        }
    }
    for (size_t i = 0; i < T.size(); i++)
        for (size_t j = 0; j < 3; j++)
            T[i].Vids[j] = tag_indices[T[i].Vids[j]];

    orient_triangle_mesh_index(V, T);
    //build_tri_mesh_data_structure(tmi_sur);
    //Mesh triMesh(V, T, TRIANGLE);
    triMesh.V = V;
    triMesh.C = T;
    triMesh.m_cellType = TRIANGLE;
    triMesh.BuildAllConnectivities();
    triMesh.ExtractBoundary();
    triMesh.GetNormalOfSurfaceFaces();
    triMesh.GetNormalOfSurfaceVertices();
    triMesh.ClassifyVertexTypes();

    MeshFileWriter writer(triMesh, "tri.vtk");
    writer.WriteFile();
    std::vector<int> vertexTypes(triMesh.V.size());
    for (size_t i = 0; i < triMesh.V.size(); i++)
        vertexTypes[i] = triMesh.V[i].type;
    writer.WritePointData(vertexTypes, "tri.vtk");
}

//size_t MeshOpt::OptimizeSurfaceVertices(std::vector<Trip>& A_Entries, std::vector<float>& b, size_t& row)
//{
////    static int iter = 1;
////    std::string name = std::string("typeAb.") + std::to_string(iter++) + ".txt";
////    std::ofstream ofs(name.c_str());
//    //float beta_root = sqrt(beta), alpha_root = sqrt(alpha);
//    size_t a_id = 0;
//    for (size_t i = 0; i < mesh.V.size(); i++) {
//        const Vertex& v = mesh.V.at(i);
//        if (v.isBoundary) {
//            int tvid = mesh.V[i].triVid;
//            const Vertex& triV = triMesh.V[tvid];
//            if (triV.type == REGULAR) {
//                //beta(nv + d)
//                const glm::vec3& implicit_n = triV.normal;
//                float implicit_d = glm::dot(implicit_n, triV.xyz());
//                for (size_t k = 0; k < 3; k++) {
//                    A_Entries.push_back(Trip(row, 3 * i + k, beta * implicit_n[k]));
//                }
//                //std::cout << "implicit_d = " << implicit_d << "\n";
//                b.push_back(beta * implicit_d);
//                row++;
//
////                for (size_t k = 0; k < 3; k++) {
////                    A_Entries.push_back(Trip(row, 3 * i + k, beta));
////                    b.push_back(beta * v[k]);
////                    row++;
////                }
//                //ofs << "triVid = " << tvid << " implicit_n = (" << implicit_n.x << "," << implicit_n.y << "," << implicit_n.z << ")\n";
//            }
//            else if (triV.type == FEATURE) {
//                int a_index = a_id + 3 * mesh.V.size();
//                for (size_t k = 0; k < 3; k++) {
//                    A_Entries.push_back(Trip(row, 3 * i + k, alpha));
//                    A_Entries.push_back(Trip(row, a_index, -alpha * triV.tangent[k]));
//                    b.push_back(alpha * triV[k]);
//                    row++;
//                }
//                //ofs << "triVid = " << tvid << " tangent = (" << triV.tangent.x << "," << triV.tangent.y << "," << triV.tangent.z << ")\n";
//            }
//            else if (triV.type == CORNER) {
//                for (size_t k = 0; k < 3; k++) {
//                    A_Entries.push_back(Trip(row, 3 * i + k, alpha));
//                    b.push_back(alpha * v[k]);
//                    row++;
//                }
//                //ofs << "triVid = " << tvid << " corner = (" << triV.z << "," << triV.y << "," << triV.z << ")\n";
//            }
//        }
//    }
//    a_id = 0;
//    for (size_t i = 0; i < mesh.V.size(); i++) {
//        const Vertex& v = mesh.V.at(i);
//        size_t tvid = mesh.V[i].triVid;
//        const Vertex& triV = triMesh.V[tvid];
//        if (v.isBoundary && triV.type == FEATURE) {
//            int a_index = a_id++ + 3 * mesh.V.size();
//            A_Entries.push_back(Trip(row, a_index, alpha));
//            b.push_back(0);
//            row++;
//        }
//    }
//
//    return a_id;
//}
size_t MeshOpt::OptimizeSurfaceVertices(std::vector<Trip>& A_Entries, std::vector<float>& b, size_t& row)
{
    size_t a_id = 0;
    for (size_t i = 0; i < mesh.V.size(); i++) {
        const Vertex& v = mesh.V.at(i);
        if (!v.isBoundary)
            continue;
        if (v.type == REGULAR) {
            //beta(nv + d)
            const glm::vec3& implicit_n = v.normal;
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

void MeshOpt::AdjustOutlayer()
{
    for (size_t i = 0; i < mesh.V.size(); i++) {
        Vertex& v = mesh.V.at(i);
        if (v.isBoundary)
            continue;
        for (size_t j = 0; j < v.N_Vids.size(); j++) {
            const Vertex& bv = mesh.V.at(v.N_Vids.at(j));
            if (bv.isBoundary) {
                double targetLen = 0.0;
                int count = 0;
                for (size_t k = 0; k < bv.N_Vids.size(); k++) {
                    const Vertex& bv1 = mesh.V.at(bv.N_Vids.at(k));
                    if (bv.isBoundary) {
                        targetLen += glm::length(glm::vec3(bv.x - bv1.x, bv.y - bv1.y, bv.z - bv1.z));
                        count++;
                    }
                }
                targetLen /= 2*count;

                glm::vec3 start(bv.x, bv.y, bv.z);
                v.x = start.x - bv.normal.x*targetLen;
                v.y = start.y - bv.normal.y*targetLen;
                v.z = start.z - bv.normal.z*targetLen;
            }
        }
    }
}

void MeshOpt::OptimizeInnerVertices(std::vector<Trip>& A_Entries, std::vector<float>& b, size_t& row)
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

    for (size_t i = 0; i < mesh.V.size(); i++) {
        Vertex& v = mesh.V.at(i);
        if (v.isBoundary)
            continue;
        for (size_t j = 0; j < v.N_Vids.size(); j++) {
            const Vertex& bv = mesh.V.at(v.N_Vids.at(j));
            if (bv.isBoundary) {
                double targetLen = 0.0;
                int count = 0;
                for (size_t k = 0; k < bv.N_Vids.size(); k++) {
                    const Vertex& bv1 = mesh.V.at(bv.N_Vids.at(k));
                    if (bv.isBoundary) {
                        targetLen += glm::length(glm::vec3(bv.x - bv1.x, bv.y - bv1.y, bv.z - bv1.z));
                        count++;
                    }
                }
                targetLen /= 2*count;

                glm::vec3 start(bv.x, bv.y, bv.z);
//                v.x = start.x - bv.normal.x*targetLen;
//                v.y = start.y - bv.normal.y*targetLen;
//                v.z = start.z - bv.normal.z*targetLen;
                const glm::vec3 innerV(start.x - bv.normal.x*targetLen, start.y - bv.normal.y*targetLen, start.z - bv.normal.z*targetLen);
                for (size_t n = 0; n < 3; n++)
                {
                    A_Entries.push_back(Trip(row, 3 * v.id + n, alpha));
                    b.push_back(alpha * innerV[n]);
                    row++;
                }
            }
        }
    }
//    MeshFileReader reader("hex.vtk");
//    const Mesh& hex = reader.GetMesh();
//
//    //std::vector<size_t> frameIds;
//    for (size_t k = 0; k < mesh.V.size(); k++) {
//        const Vertex& v = mesh.V.at(k);
//        if (v.isBoundary) {
//            for (size_t j = 0; j < v.N_Vids.size(); j++) {
//                Vertex& innerV = mesh.V.at(v.N_Vids.at(j));
//                if (innerV.isBoundary) {
//                    //frameIds.push_back(innerFrame.id);
//                    for (size_t n = 0; n < 3; n++)
//                    {
//                        A_Entries.push_back(Trip(row, 3 * innerV.id + n, alpha));
//                        b.push_back(alpha * v[n]);
//                        row++;
//                    }
//                }
//            }
//        }
//    }
    //framefield.WriteFile("innerframe.vtk", frameIds);
}

void MeshOpt::ComputeMeshTargetLength()
{
    std::cout << "=============================\n";
    std::cout << "   Computing TargetLength    \n";

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

    for (size_t k = 0; k < mesh.E.size(); k++){
        const Edge& edge = mesh.E.at(k);
        if (!edge.isSingularity)
            continue;
        if (edge.N_Cids.size() != 5 && edge.N_Cids.size() != 3 && edge.N_Cids.size() != 6)
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

                coefficients.push_back(Trip(row, edge1.id, 1));
                coefficients.push_back(Trip(row, edge2.id, -1));
                b.push_back(0);
                row++;
            }
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
    Eigen::SimplicialLDLT<SpMat> chol(ATA);
    X = chol.solve(ATB);

    // ------- Output File -----------
//    std::ofstream ofs("iX.txt");
//    //ofs << A;
//    ofs << X;
    // ------- ----------- -----------

    for (size_t i = 0; i < mesh.E.size(); i++)
        mesh.E[i].length = X[i];

    std::cout << "=============================\n";
}

void MeshOpt::SetTriMesh(Mesh& triMesh)
{
    this->triMesh = triMesh;
}

void MeshOpt::SetRefMesh(Mesh& refMesh)
{
    this->m_refMesh = &refMesh;
}

void MeshOpt::SetTargetSurfaceMesh(Mesh& refMesh)
{
    this->m_targetSurfaceMesh = &refMesh;
}

void MeshOpt::SetAlpha(const double value/* = 100.0*/)
{
    alpha = value;
}

void MeshOpt::SetBeta(const double value/* = 100.0*/)
{
    beta = value;
}

void MeshOpt::SetGamma(const double value/* = 100.0*/)
{
    gamma = value;
}

void MeshOpt::SetStepSize(const double value/* = 100.0*/)
{
    stepSize = value;
}

void MeshOpt::SetAnisotropy(const double value/* = 100.0*/)
{
    anisotropy = value;
}

void MeshOpt::SetMinScaledJacobian(const double value/* = 100.0*/)
{
    minScaledJacobian = value;
}

void MeshOpt::SetUseAverageTargetLength(bool value)
{
    useAverageTargetLength = value;
}

void MeshOpt::SetRecoverable(bool value)
{
    recoverable = value;
}

void MeshOpt::SetAllowBigStep(bool value)
{
    allowBigStep = value;
}

void MeshOpt::SetProjectToTargetSurface(bool value)
{
    projectToTargetSurface = value;
}

void MeshOpt::SetUseProjection(bool value /*value = false*/)
{
    useProjection = value;
}

void MeshOpt::SetChangeBoundary(bool value)
{
    changeBoundary = value;
}

void MeshOpt::OutputBadCells(const std::vector<size_t>& badCellIds, const char* filename)
{
    std::vector<Cell> cells(badCellIds.size());
    for (size_t i = 0; i < badCellIds.size(); i++)
        cells.at(i) = mesh.C.at(badCellIds.at(i));
    MeshFileWriter writer(mesh.V, cells, filename);
    writer.WriteFile();
}

std::vector<size_t> MeshOpt::OutputBadCellsAndExtendCells(const std::vector<size_t>& badCellIds, const char* filename)
{
    std::vector<size_t> BadCellIds = badCellIds;
    std::vector<size_t> cellIds;
    size_t count = 2;
    while (count-- != 0)
    {
        std::vector<size_t> badCellVIds;
        for (size_t i = 0; i < BadCellIds.size(); i++)
            std::copy(mesh.C.at(BadCellIds.at(i)).Vids.begin(), mesh.C.at(BadCellIds.at(i)).Vids.end(), back_inserter(badCellVIds));
        std::sort(badCellVIds.begin(), badCellVIds.end());
        std::vector<size_t>::iterator iter = std::unique(badCellVIds.begin(), badCellVIds.end());
        badCellVIds.resize(std::distance(badCellVIds.begin(), iter));

        //std::vector<size_t> cellIds;
        cellIds.clear();
        for (size_t i = 0; i < badCellVIds.size(); i++) {
            const Vertex& v = mesh.V.at(badCellVIds.at(i));
            std::copy(v.N_Cids.begin(), v.N_Cids.end(), back_inserter(cellIds));
        }
        std::sort(cellIds.begin(), cellIds.end());
        iter = std::unique(cellIds.begin(), cellIds.end());
        cellIds.resize(std::distance(cellIds.begin(), iter));

        BadCellIds = cellIds;
    }
    std::vector<Cell> cells(cellIds.size());
    for (size_t i = 0; i < cellIds.size(); i++)
        cells.at(i) = mesh.C.at(cellIds.at(i));

    MeshFileWriter writer(mesh.V, cells, filename);
    writer.WriteFile();

    return cellIds;
}

std::vector<size_t> MeshOpt::GetBadCellsAndExtendCells(const std::vector<size_t>& badCellIds, int extendLayers/* = 2*/)
{
    std::vector<size_t> BadCellIds = badCellIds;
    std::vector<size_t> cellIds;
    size_t count = extendLayers;
    while (count-- != 0)
    {
        std::vector<size_t> badCellVIds;
        for (size_t i = 0; i < BadCellIds.size(); i++)
            std::copy(mesh.C.at(BadCellIds.at(i)).Vids.begin(), mesh.C.at(BadCellIds.at(i)).Vids.end(), back_inserter(badCellVIds));
        std::sort(badCellVIds.begin(), badCellVIds.end());
        std::vector<size_t>::iterator iter = std::unique(badCellVIds.begin(), badCellVIds.end());
        badCellVIds.resize(std::distance(badCellVIds.begin(), iter));

        cellIds.clear();
        for (size_t i = 0; i < badCellVIds.size(); i++) {
            const Vertex& v = mesh.V.at(badCellVIds.at(i));
            std::copy(v.N_Cids.begin(), v.N_Cids.end(), back_inserter(cellIds));
        }
        std::sort(cellIds.begin(), cellIds.end());
        iter = std::unique(cellIds.begin(), cellIds.end());
        cellIds.resize(std::distance(cellIds.begin(), iter));

        BadCellIds = cellIds;
    }

    return cellIds;
}

void MeshOpt::OutputEdgesOfBadCells(const std::vector<size_t>& badCellIds, const char* filename)
{
    std::vector<size_t> Vids;
    for (size_t i = 0; i < badCellIds.size(); i++) {
        const Cell& cell = mesh.C.at(badCellIds.at(i));
        std::copy(cell.Vids.begin(), cell.Vids.end(), back_inserter(Vids));
    }
    std::sort(Vids.begin(), Vids.end());
    std::vector<size_t>::iterator iter = std::unique(Vids.begin(), Vids.end());
    Vids.resize(std::distance(Vids.begin(), iter));

    std::vector<size_t> badEdgeIds;
    for (size_t i = 0; i < Vids.size(); i++) {
        const size_t vid = Vids.at(i);
        for (size_t k = 0; k < mesh.E.size(); k++) {
            const Edge& e = mesh.E.at(k);
            if (  e.Vids[0] == vid || e.Vids[1] == vid)
                badEdgeIds.push_back(e.id);
        }
    }

    MeshFileWriter writer(mesh, filename);
    writer.WriteEdgesVtk(badEdgeIds);
}
