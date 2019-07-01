/*
 * FrameOpt.cpp
 *
 *  Created on: Nov 15, 2016
 *      Author: cotrik
 */

#include "FrameOpt.h"
#include "Mesh.h"
#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "MeshQuality.h"

#include <algorithm>
#include <vector>
#include <iostream>
#include <math.h>

//FrameOpt::FrameOpt()
FrameOpt::FrameOpt(const Mesh& mesh, const FrameField& framefield, const PolyLines& polyLines,
        const double alpha/* = 1000000.0*/, const double beta/* = 1.0*/, const double gamma/* = 1.0*/)
: mesh((Mesh&)mesh)
, framefield((FrameField&)framefield)
, polyLines((PolyLines&)polyLines)
, alpha(alpha)
, beta(beta)
, gamma(gamma)
, anisotropy(0.05)
, stepSize(1.0)
, avgMeshEdgeLength(1.0)
, avgFrameEdgeLength(1.0)
, useAverageTargetLength(false)
, recoverable(true)
, m_numOfInvertdElements(MAXID)
{
    // TODO Auto-generated constructor stub

}

FrameOpt::~FrameOpt()
{
    // TODO Auto-generated destructor stub
}

void FrameOpt::Run(const size_t iters/* = 1*/)
{
    double sumEdgeLength = 0.0;
    size_t numOfBoundaryEdges = 0;
    for (size_t i = 0; i < mesh.E.size(); i++)
    {
        if (mesh.E.at(i).isBoundary){
            const Vertex& v1 = mesh.V.at(mesh.E.at(i).Vids.at(0));
            const Vertex& v2 = mesh.V.at(mesh.E.at(i).Vids.at(1));
            mesh.E[i].length = glm::length(v1 - v2);
            sumEdgeLength += mesh.E[i].length;
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
    //if (useAverageTargetLength)
        ComputeMeshTargetLength();
    while (!converged && iter++ < iters)
    {
        if (!initUseAverageTargetLength  && iter == 1)
            useAverageTargetLength = true;
        std::cout << "stepSize = " << stepSize << std::endl;
        converged = OptimizeFrame();
        UpdateMeshFromFrameField();
        //UpdateFrameFieldFromMesh();

        if (!initUseAverageTargetLength && iter == 1)
            useAverageTargetLength = false;
        stepSize *= initStepSize;

        std::string filename = std::string("FrameOpt.") + std::to_string(iter) + ".vtk";
        framefield.WriteColorFile(filename.c_str());

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
        filename = std::string("FrameOfBadCells.") + std::to_string(iter) + ".vtk";
        OutputFramesOfBadCells(badCellIds, filename.c_str());
        if (prevMinimumScaledJacobian > 0 && minimumScaledJacobian < prevMinimumScaledJacobian)
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
}

bool FrameOpt::Optimize()
{
    std::vector<Trip> A_Entries;
    std::vector<float> b;
    size_t row = 0;

    OptimizeBoundary(A_Entries, b, row);
    OptimizeOrthogonality(A_Entries, b, row);
    std::cout << "Orthogonality row = " << row << std::endl;
    std::cout << "Orthogonality A_Entries = " << A_Entries.size() << std::endl;
    OptimizeStraightness(A_Entries, b, row);
    std::cout << "Straightness row = " << row << std::endl;
    std::cout << "Straightness A_Entries = " << A_Entries.size() << std::endl;
    const size_t col = 3 * mesh.V.size();
    std::cout << " col = " << col << std::endl;

    Eigen::SparseMatrix<float> A(row, col);
    A.setFromTriplets(A_Entries.begin(), A_Entries.end());
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
    std::vector<float> changes;
    for (size_t i = 0; i < mesh.V.size(); i++)
    {
        if (fabs(mesh.V[i].x - X[3 * i + 0]) > 1e-6
            || fabs(mesh.V[i].y - X[3 * i + 1]) > 1e-6
            || fabs(mesh.V[i].z - X[3 * i + 2]) > 1e-6)
        {
            converged = false;
        }
        mesh.V[i].x = stepSize * X[3 * i + 0] + (1.0 - stepSize) * mesh.V[i].x;
        mesh.V[i].y = stepSize * X[3 * i + 1] + (1.0 - stepSize) * mesh.V[i].y;
        mesh.V[i].z = stepSize * X[3 * i + 2] + (1.0 - stepSize) * mesh.V[i].z;
    }

    return converged;
}

bool FrameOpt::OptimizeFrame()
{
    std::vector<Trip> A_Entries;
    std::vector<float> b;
    size_t row = 0;

    OptimizeFrameBoundary(A_Entries, b, row);
    std::cout << "OptimizeFrameBoundary" << std::endl;
    //OptimizeFrameInnerVertices(A_Entries, b, row);
    //std::cout << "OptimizeFrameInnerVertices" << std::endl;
    //std::cout << "Boundary row = " << row << std::endl;
    //std::cout << "Boundary A_Entries = " << A_Entries.size() << std::endl;
    OptimizeFrameOrthogonality(A_Entries, b, row);
    std::cout << "OptimizeFrameOrthogonality" << std::endl;
    //std::cout << "Orthogonality row = " << row << std::endl;
    //std::cout << "Orthogonality A_Entries = " << A_Entries.size() << std::endl;
    OptimizeFrameStraightness(A_Entries, b, row);
    std::cout << "OptimizeFrameStraightness" << std::endl;
    //std::cout << "Straightness row = " << row << std::endl;
    //std::cout << "Straightness A_Entries = " << A_Entries.size() << std::endl;
    const size_t col = 3 * framefield.frameNodes.size();
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

    std::vector<glm::dvec3> oldV(mesh.V.size());
    for (size_t i = 0; i < mesh.V.size(); i++)
    {
        oldV[i].x = mesh.V[i].x;
        oldV[i].y = mesh.V[i].y;
        oldV[i].z = mesh.V[i].z;
    }

    bool converged = true;

    if (recoverable) {
        for (size_t i = 0; i < framefield.frameNodes.size(); i++)
        {
            if (fabs(double(framefield.frameNodes[i].x - X[3 * i + 0])) > 1e-6
                || fabs(double(framefield.frameNodes[i].y - X[3 * i + 1])) > 1e-6
                || fabs(double(framefield.frameNodes[i].z - X[3 * i + 2])) > 1e-6)
            {
                converged = false;
            }
            if (!allowBigStep)
            if (fabs(double(mesh.V[i].x - X[3 * i + 0])) > 2.0 * avgFrameEdgeLength
                || fabs(double(mesh.V[i].y - X[3 * i + 1])) > 2.0 * avgFrameEdgeLength
                || fabs(double(mesh.V[i].z - X[3 * i + 2])) > 2.0 * avgFrameEdgeLength)
            {
                continue;
            }
            framefield.frameNodes[i].x = stepSize * X[3 * i + 0] + (1.0 - stepSize) * framefield.frameNodes[i].x;
            framefield.frameNodes[i].y = stepSize * X[3 * i + 1] + (1.0 - stepSize) * framefield.frameNodes[i].y;
            framefield.frameNodes[i].z = stepSize * X[3 * i + 2] + (1.0 - stepSize) * framefield.frameNodes[i].z;
        }

        MeshFileWriter writer(mesh, "temp.vtk");
        writer.WriteFile();
        double minimumScaledJacobian = 0.0;
        double averageScaledJacobian = 0.0;
        double maximumScaledJacobian = 0.0;
        std::vector<size_t> badCellIds1;
        size_t InvertedElements = GetQuality("temp.vtk", minimumScaledJacobian, averageScaledJacobian, maximumScaledJacobian, badCellIds1);

        if (InvertedElements > m_numOfInvertdElements){
            for (size_t i = 0; i < mesh.V.size(); i++)
            {
                mesh.V[i].x = oldV[i].x;
                mesh.V[i].y = oldV[i].y;
                mesh.V[i].z = oldV[i].z;
            }
            UpdateFrameFieldFromMesh();
        }
    }
    else {
        for (size_t i = 0; i < framefield.frameNodes.size(); i++)
        {
            if (fabs(double(framefield.frameNodes[i].x - X[3 * i + 0])) > 1e-6
                || fabs(double(framefield.frameNodes[i].y - X[3 * i + 1])) > 1e-6
                || fabs(double(framefield.frameNodes[i].z - X[3 * i + 2])) > 1e-6)
            {
                converged = false;
            }
            if (!allowBigStep)
            if (fabs(double(mesh.V[i].x - X[3 * i + 0])) > 2.0 * avgFrameEdgeLength
                || fabs(double(mesh.V[i].y - X[3 * i + 1])) > 2.0 * avgFrameEdgeLength
                || fabs(double(mesh.V[i].z - X[3 * i + 2])) > 2.0 * avgFrameEdgeLength)
            {
                continue;
            }
            framefield.frameNodes[i].x = stepSize * X[3 * i + 0] + (1.0 - stepSize) * framefield.frameNodes[i].x;
            framefield.frameNodes[i].y = stepSize * X[3 * i + 1] + (1.0 - stepSize) * framefield.frameNodes[i].y;
            framefield.frameNodes[i].z = stepSize * X[3 * i + 2] + (1.0 - stepSize) * framefield.frameNodes[i].z;
        }
    }
    return converged;
}

void FrameOpt::OptimizeOrthogonality(std::vector<Trip>& A_Entries, std::vector<float>& b, size_t& row)
{
    for (size_t k = 0; k < framefield.frameNodes.size(); k++){
        const Frame& frame = framefield.frameNodes.at(k);
        if (frame.isBoundary)
            continue;
        for (size_t i = 0; i < frame.N_Eids.size(); i++) {
            const size_t iId = frame.N_Vids.at(i);
            const Frame& framei = framefield.frameNodes.at(iId);
            const double length_vi = framefield.frameEdges.at(frame.N_Eids.at(i)).length;
            // vi is variable, and use vj as constant
            for (size_t j = 0; j < frame.ortho4Eids.at(i).size(); j++) {
                const size_t jId = frame.ortho4Vids.at(i).at(j);
                const Frame& framej = framefield.frameNodes.at(jId);
                const glm::dvec3 vj = glm::normalize(glm::dvec3(frame.x - framej.x, frame.y - framej.y, frame.z - framej.z));
                // const double length_vi = glm::length(glm::dvec3(frame.x - framei.x, frame.y - framei.y, frame.z - framei.z));
                // <vi/length_vi, vj> = 0
                for (size_t n = 0; n < 3; n++){
                    //if (!frame.isBoundary)
                    {
                        const Cell& cell = mesh.C.at(frame.id);
                        const double v = vj[n] / length_vi * 0.125;
                        for (size_t m = 0; m < 8; m++)
                            A_Entries.push_back(Trip(row, 3 * cell.Vids[m] + n, v));
                    }
                    if (!framei.isBoundary) {
                        const Cell& celli = mesh.C.at(framei.id);
                        const double v = -vj[n] / length_vi * 0.125;
                        for (size_t m = 0; m < 8; m++)
                            A_Entries.push_back(Trip(row, 3 * celli.Vids[m] + n, v));
                    }
                    else {
                        const size_t eiId = frame.N_Eids.at(i);
                        const FrameEdge& frameEdgei = framefield.frameEdges.at(eiId);
                        const Face& facei = mesh.F.at(frameEdgei.id);
                        const double v = -vj[n] / length_vi * 0.25;
                        for (size_t m = 0; m < 4; m++)
                            A_Entries.push_back(Trip(row, 3 * facei.Vids[m] + n, v));
                    }
                }
                b.push_back(0);
                row++;
            }
        }
    }
}
void FrameOpt::OptimizeBoundary(std::vector<Trip>& A_Entries, std::vector<float>& b, size_t& row)
{
    for (size_t k = 0; k < framefield.frameNodes.size(); k++){
        const Frame& frame = framefield.frameNodes.at(k);
        if (frame.isBoundary){
            for (size_t n = 0; n < 3; n++)
            {
                const size_t eiId = frame.N_Eids.at(0);
                const FrameEdge& frameEdgei = framefield.frameEdges.at(eiId);
                const Face& facei = mesh.F.at(frameEdgei.id);
                for (size_t m = 0; m < facei.Vids.size(); m++)
                {
                    A_Entries.push_back(Trip(row, 3 * facei.Vids[m] + n, alpha));
                    b.push_back(alpha * mesh.V.at(facei.Vids[m])[n]);
                    row++;
                }
            }
        }
    }
}
void FrameOpt::OptimizeStraightness(std::vector<Trip>& A_Entries, std::vector<float>& b, size_t& row)
{
    for (size_t k = 0; k < polyLines.polyLines.size(); k++) {
        const PolyLine& polyline = polyLines.polyLines.at(k);
        for (size_t j = 0; j < polyline.Eids.size() - 1; j++) {
            const size_t frameEdgeId1 = polyline.Eids.at(j);
            const size_t frameEdgeId2 = polyline.Eids.at(j + 1);
            const FrameEdge& frameEdge1 = framefield.frameEdges.at(frameEdgeId1);
            const FrameEdge& frameEdge2 = framefield.frameEdges.at(frameEdgeId2);
            const Frame& frame1 = framefield.frameNodes.at(polyline.Vids.at(j));
            const Frame& frameC = framefield.frameNodes.at(polyline.Vids.at(j + 1));
            const Frame& frame2 = framefield.frameNodes.at(polyline.Vids.at(j + 2));
            const glm::dvec3 v2 = glm::dvec3(frameC.x - frame2.x, frameC.y - frame2.y, frameC.z - frame2.z);
            const glm::dvec3 v2n = glm::normalize(v2);
            const double length_v1 = frameEdge1.length;
            // const double length_v1 = glm::length(glm::dvec3(frameC.x - frame1.x, frameC.y - frame1.y, frameC.z - frame1.z));
            // v1 is variable, and use v2 as constant
            // <v1/length_v1, v2n> = -1
            for (size_t n = 0; n < 3; n++){
            //if (!frameC.isBoundary)
                {
                    const Cell& cellC = mesh.C.at(frameC.id);
                    const double v = v2n[n] / length_v1 * 0.125;
                    for (size_t m = 0; m < 8; m++)
                        A_Entries.push_back(Trip(row, 3 * cellC.Vids[m] + n, v));
                }

                if (!frame1.isBoundary) {
                    const Cell& cell1 = mesh.C.at(frame1.id);
                    const double v = -v2n[n] / length_v1 * 0.125;
                    for (size_t m = 0; m < 8; m++)
                        A_Entries.push_back(Trip(row, 3 * cell1.Vids[m] + n, v));
                }
                else {
                    const Face& face1 = mesh.F.at(frameEdge1.id);
                    const double v = -v2n[n] / length_v1 * 0.25;
                    for (size_t m = 0; m < 4; m++)
                        A_Entries.push_back(Trip(row, 3 * face1.Vids[m] + n, v));
                }
            }
            b.push_back(-1);
            row++;
        }
    }
}


void FrameOpt::OptimizeFrameOrthogonality(std::vector<Trip>& A_Entries, std::vector<float>& b, size_t& row)
{
    const size_t prevA_EntriesSize = A_Entries.size();
    A_Entries.resize(prevA_EntriesSize + mesh.C.size() * 6 * 4 * 3 * 2);
    size_t id = prevA_EntriesSize;
    for (size_t k = 0; k < framefield.frameNodes.size(); k++){
        const Frame& frame = framefield.frameNodes.at(k);
        if (frame.isBoundary)
            continue;
        const Cell& cell = mesh.C.at(frame.id);
        for (size_t i = 0; i < frame.N_Eids.size(); i++) {                                                                    /*6*/
            const size_t iId = frame.N_Vids.at(i);
            const Frame& framei = framefield.frameNodes.at(iId);
            double length_vi = framefield.frameEdges.at(frame.N_Eids.at(i)).length;
            if (!useAverageTargetLength) {
                length_vi = glm::length(glm::dvec3(frame.x - framei.x, frame.y - framei.y, frame.z - framei.z));
                // length_vi = length_vi > 1e-4 ? length_vi : 1e-4;
                length_vi = length_vi > avgFrameEdgeLength*anisotropy ? length_vi : avgFrameEdgeLength*anisotropy;
            }
            // vi is variable, and use vj as constant
            for (size_t j = 0; j < frame.ortho4Eids.at(i).size(); j++) {                                                      /*4*/
                const size_t jId = frame.ortho4Vids.at(i).at(j);
                const Frame& framej = framefield.frameNodes.at(jId);
                const glm::dvec3 vj = glm::dvec3(frame.x - framej.x, frame.y - framej.y, frame.z - framej.z);
                const glm::dvec3 vjn = glm::normalize(vj);
                //double length_vi = glm::length(glm::dvec3(frame.x - framei.x, frame.y - framei.y, frame.z - framei.z));
                //length_vi = length_vi > 1e-4 ? length_vi : 1e-4;
                //double length_vi = framefield.frameEdges.at(frame.ortho4Eids.at(i).at(j)).length;
                // <vi/length_vi, vjn> = 0
                for (size_t n = 0; n < 3; n++) {                                                                              /*3*/
                    double v = gamma * vjn[n] / length_vi;
                    if (std::isnan(v)){
                        v = 0.0;
                    }
//                    A_Entries.push_back(Trip(row, 3 * frame.id + n, v));
//                    A_Entries.push_back(Trip(row, 3 * framei.id + n, -v));
                    A_Entries.at(id++) = Trip(row, 3 * frame.id + n, v);                                                      /*2*/
                    A_Entries.at(id++) = Trip(row, 3 * framei.id + n, -v);
                }
                b.push_back(0);
                row++;
            }
        }
    }
}

void FrameOpt::OptimizeFrameStraightness(std::vector<Trip>& A_Entries, std::vector<float>& b, size_t& row)
{
    for (size_t k = 0; k < polyLines.polyLines.size(); k++) {
        const PolyLine& polyline = polyLines.polyLines.at(k);
        for (size_t j = 0; j < polyline.Eids.size() - 1; j++) {
            const size_t frameEdgeId1 = polyline.Eids.at(j);
            const size_t frameEdgeId2 = polyline.Eids.at(j + 1);
            const FrameEdge& frameEdge1 = framefield.frameEdges.at(frameEdgeId1);
            const FrameEdge& frameEdge2 = framefield.frameEdges.at(frameEdgeId2);
            const Frame& frame1 = framefield.frameNodes.at(polyline.Vids.at(j));
            const Frame& frameC = framefield.frameNodes.at(polyline.Vids.at(j + 1));
            const Frame& frame2 = framefield.frameNodes.at(polyline.Vids.at(j + 2));
            const glm::dvec3 v2 = glm::dvec3(frameC.x - frame2.x, frameC.y - frame2.y, frameC.z - frame2.z);
            const glm::dvec3 v2n = glm::normalize(v2);
            double length_v1 = frameEdge1.length;
            if (!useAverageTargetLength) {
                length_v1 = glm::length(glm::dvec3(frameC.x - frame1.x, frameC.y - frame1.y, frameC.z - frame1.z));
                //length_v1 = length_v1 > 1e-4 ? length_v1 : 1e-4;
                length_v1 = length_v1 > avgFrameEdgeLength*anisotropy ? length_v1 : avgFrameEdgeLength*anisotropy;
            }
            // v1 is variable, and use v2 as constant
            // <v1/length_v1, v2n> = -1
            for (size_t n = 0; n < 3; n++) {
                const double v = gamma * v2n[n] / length_v1;
                A_Entries.push_back(Trip(row, 3 * frameC.id + n, v));
                A_Entries.push_back(Trip(row, 3 * frame1.id + n, -v));
            }
            b.push_back(-1 * gamma);
            row++;
        }
    }
}

void FrameOpt::OptimizeFrameBoundary(std::vector<Trip>& A_Entries, std::vector<float>& b, size_t& row)
{
    for (size_t k = 0; k < framefield.frameNodes.size(); k++){
        const Frame& frame = framefield.frameNodes.at(k);
        if (frame.isBoundary){
            for (size_t n = 0; n < 3; n++)
            {
                A_Entries.push_back(Trip(row, 3 * frame.id + n, alpha));
                b.push_back(alpha * frame[n]);
                row++;
            }
        }
    }
}

//void FrameOpt::OptimizeFrameInnerVertices(std::vector<Trip>& A_Entries, std::vector<float>& b, size_t& row)
//{
//    /*
//     *      3
//     *     /|
//     *    / |
//     *   /  |
//     * 0/   |
//     * |    |
//     * | f*-|-------->* v(innerFrame)
//     * |    |
//     * |   / 2
//     * |  /
//     * | /
//     * |/
//     * 1
//     * */
//    std::vector<size_t> frameIds;
//    for (size_t k = 0; k < framefield.frameNodes.size(); k++) {
//        const Frame& frame = framefield.frameNodes.at(k);
//        if (frame.isBoundary) {
//            const Face& face = mesh.F.at(frame.id);
//            const Vertex& v0 = mesh.V.at(face.Vids[0]);
//            const Vertex& v1 = mesh.V.at(face.Vids[1]);
//            const Vertex& v2 = mesh.V.at(face.Vids[2]);
//            const Vertex& v3 = mesh.V.at(face.Vids[3]);
//            const float l0 = glm::length(v1 - v0);
//            const float l1 = glm::length(v2 - v1);
//            const float l2 = glm::length(v3 - v2);
//            const float l3 = glm::length(v0 - v3);
//            const double avgl = 0.25 * (l0 + l1 + l2 + l3);
//            const glm::dvec3 dir = glm::cross(v2 - v1, v0 - v1); // + glm::cross(v0 - v3, v2 - v3);
//            const glm::dvec3 ndir = glm::normalize(dir);
//            const glm::dvec3 f(frame.x, frame.y, frame.z);
//            const double fvlen = 0.5 * avgl;
//            const glm::dvec3 v(f + glm::dot(glm::dvec3(fvlen, fvlen, fvlen), ndir));  // innerFrame
//            Frame& innerFrame = framefield.frameNodes.at(frame.N_Vids.at(0)); // only 1 neighbor for the boundary frame;
//            frameIds.push_back(innerFrame.id);
//            for (size_t n = 0; n < 3; n++)
//            {
//                innerFrame[n] = v[n];
//                A_Entries.push_back(Trip(row, 3 * innerFrame.id + n, alpha));
//                b.push_back(alpha * v[n]);
//                row++;
//            }
//        }
//    }
//    framefield.WriteFile("innerframe.vtk", frameIds);
//}

void FrameOpt::OptimizeFrameInnerVertices(std::vector<Trip>& A_Entries, std::vector<float>& b, size_t& row)
{
    /*
     *      3
     *     /|
     *    / |
     *   /  |
     * 0/   |
     * |    |
     * | f*-|-------->* v(innerFrame)
     * |    |
     * |   / 2
     * |  /
     * | /
     * |/
     * 1
     * */
    MeshFileReader reader("hex.vtk");
    const Mesh& hex = reader.GetMesh();

    //std::vector<size_t> frameIds;
    for (size_t k = 0; k < framefield.frameNodes.size(); k++) {
        const Frame& frame = framefield.frameNodes.at(k);
        if (frame.isBoundary) {
            Frame& innerFrame = framefield.frameNodes.at(frame.N_Vids.at(0)); // only 1 neighbor for the boundary frame;
            //frameIds.push_back(innerFrame.id);
            const Cell& c = hex.C.at(innerFrame.id);
            glm::dvec3 v(0.0, 0.0, 0.0);
            for (size_t i = 0; i < c.Vids.size(); i++)
                v = v + hex.V.at(c.Vids[i]);
            v.x /= 8;
            v.y /= 8;
            v.z /= 8;
            for (size_t n = 0; n < 3; n++)
            {
                A_Entries.push_back(Trip(row, 3 * innerFrame.id + n, alpha));
                b.push_back(alpha * v[n]);
                row++;
            }
        }
    }
    //framefield.WriteFile("innerframe.vtk", frameIds);
}


void FrameOpt::UpdateMeshFromFrameField()
{
    for (size_t i = 0; i < mesh.V.size(); i++)
    {
        if (mesh.V[i].isBoundary)
            continue;
        for (size_t j = 0; j < 3; j++)
            mesh.V[i][j] = 0;
        for (size_t j = 0; j < mesh.V[i].N_Cids.size(); j++)
        {
            int hid = mesh.V[i].N_Cids[j];
            for (size_t k = 0; k < 3; k++)
                mesh.V[i][k] += framefield.frameNodes[hid][k];
        }
        for (size_t j = 0; j < 3; j++)
            mesh.V[i][j] /= mesh.V[i].N_Cids.size();
    }
}

void FrameOpt::UpdateFrameFieldFromMesh()
{
    for (size_t i = 0; i < framefield.frameNodes.size(); i++)
    {
        Frame& frame = framefield.frameNodes.at(i);
        if (frame.isBoundary) {
            const Face& face = mesh.F.at(frame.id);
            frame.x = 0.0;
            frame.y = 0.0;
            frame.z = 0.0;
            for (size_t j = 0; j < face.Vids.size(); j++) {
                const Vertex& v = mesh.V.at(face.Vids.at(j));
                frame.x += v.x;
                frame.y += v.y;
                frame.z += v.z;
            }
            frame.x /= face.Vids.size();
            frame.y /= face.Vids.size();
            frame.z /= face.Vids.size();
        }
        else {
            const Cell& cell = mesh.C.at(frame.id);
            frame.x = 0.0;
            frame.y = 0.0;
            frame.z = 0.0;
            for (size_t j = 0; j < cell.Vids.size(); j++) {
                const Vertex& v = mesh.V.at(cell.Vids.at(j));
                frame.x += v.x;
                frame.y += v.y;
                frame.z += v.z;
            }
            frame.x /= cell.Vids.size();
            frame.y /= cell.Vids.size();
            frame.z /= cell.Vids.size();
        }

    }
}

void FrameOpt::ComputeMeshTargetLength()
{
    std::cout << "=============================\n";
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

//    std::vector<float> sortedX(mesh.E.size());
//    for (size_t i = 0; i < mesh.E.size(); i++)
//        sortedX[i] = X[i];
//    std::sort(sortedX.begin(), sortedX.end());
//    std::ofstream ofs("XuE_Length.txt");
//    //ofs << X;
//    for (size_t i = 0; i < mesh.E.size(); i++)
//        ofs << sortedX[i] << std::endl;

    double sum_len = 0.0;
    std::vector<float> length(framefield.frameEdges.size());
    for (size_t i = 0; i < framefield.frameEdges.size(); i++)
    {
        FrameEdge& frameEdge = framefield.frameEdges.at(i);

        const Frame& frame1 = framefield.frameNodes.at(frameEdge.Vids[0]);
        const Frame& frame2 = framefield.frameNodes.at(frameEdge.Vids[1]);

        std::vector<size_t> hs;
        size_t shared_fid = frameEdge.id;
        if (!frame1.isBoundary)
            hs.push_back(frame1.id);
        if (!frame2.isBoundary)
            hs.push_back(frame2.id);

        std::vector<size_t> es_f;
        std::vector<size_t> vs_f;
        for (size_t j = 0; j < 4; j++)
        {
            es_f.push_back(mesh.F[shared_fid].Eids[j]);
            vs_f.push_back(mesh.F[shared_fid].Vids[j]);
        }

        for (size_t j = 0; j < hs.size(); j++)
        {
            double ave_len = 0;
            std::vector<size_t> es;
            for (size_t k = 0; k < 6; k++)
                for (size_t m = 0; m < 4; m++)
                    es.push_back(mesh.F[mesh.C[hs[j]].Fids[k]].Eids[m]);
            set_redundent_clearn(es);
            std::vector<size_t> es_left;
            set_exclusion(es, es_f, es_left);
            es.clear();
            for (size_t k = 0; k < es_left.size(); k++)
            {
                if (set_contain(vs_f, mesh.E[es_left[k]].Vids[0]) || set_contain(vs_f, mesh.E[es_left[k]].Vids[1]))
                    es.push_back(es_left[k]);
            }
            for (size_t k = 0; k < es.size(); k++)
                ave_len += mesh.E[es[k]].length;
            ave_len /= es.size() * 2;
            frameEdge.length += ave_len;       //*0.5;
            sum_len += ave_len;
        }
//        length[i] = frameEdge.length;
    }
//    static int iter = 1;
//    std::sort(length.begin(), length.end());
//    std::string filename = std::string("XuFE_Length") + std::to_string(iter++) + ".txt";
//    std::ofstream ofsf(filename.c_str());
//    for (size_t i = 0; i < length.size(); i++)
//        ofsf << length[i] << std::endl;
    avgFrameEdgeLength = sum_len / framefield.frameEdges.size();
    std::cout << "average frame edge length " << avgFrameEdgeLength << std::endl;
    std::cout << "=============================\n";
}

void FrameOpt::SetAlpha(const double value/* = 100.0*/)
{
    alpha = value;
}

void FrameOpt::SetBeta(const double value/* = 100.0*/)
{
    beta = value;
}

void FrameOpt::SetGamma(const double value/* = 100.0*/)
{
    gamma = value;
}

void FrameOpt::SetStepSize(const double value/* = 1.0*/)
{
    stepSize = value;
}

void FrameOpt::SetAnisotropy(const double value/* = 100.0*/)
{
    anisotropy = value;
}

void FrameOpt::SetUseAverageTargetLength(bool value)
{
    useAverageTargetLength = value;
}

void FrameOpt::SetRecoverable(bool value)
{
    recoverable = value;
}

void FrameOpt::SetAllowBigStep(bool value)
{
    allowBigStep = value;
}

void FrameOpt::OutputBadCells(const std::vector<size_t>& badCellIds, const char* filename)
{
    std::vector<Cell> cells(badCellIds.size());
    for (size_t i = 0; i < badCellIds.size(); i++)
        cells.at(i) = mesh.C.at(badCellIds.at(i));
    MeshFileWriter writer(mesh.V, cells, filename);
    writer.WriteFile();
}

void FrameOpt::OutputFramesOfBadCells(const std::vector<size_t>& badCellIds, const char* filename)
{
    std::vector<Cell> cells(badCellIds.size());
    std::vector<size_t> frameIds;
    for (size_t i = 0; i < badCellIds.size(); i++) {
        const Cell& cell = mesh.C.at(badCellIds.at(i));
        for (size_t j = 0; j < cell.Vids.size(); j++) {
            const Vertex& v = mesh.V.at(cell.Vids.at(j));
//            if (v.isBoundary)
//                continue;
            std::copy(v.N_Cids.begin(), v.N_Cids.end(), back_inserter(frameIds));
        }
    }

    std::sort(frameIds.begin(), frameIds.end());
    std::vector<size_t>::iterator iter = std::unique(frameIds.begin(), frameIds.end());
    frameIds.resize(std::distance(frameIds.begin(), iter));

    framefield.WriteFile(filename, frameIds);
}
