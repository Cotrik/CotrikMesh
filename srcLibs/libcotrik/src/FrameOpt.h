/*
 * FrameOpt.h
 *
 *  Created on: Nov 15, 2016
 *      Author: cotrik
 */

#ifndef FRAME_OPT_H_
#define FRAME_OPT_H_

#include "PolyLine.h"

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Cholesky>
#include <eigen3/Eigen/LU>
#include <eigen3/Eigen/SVD>
#include <eigen3/Eigen/SparseCore>
using namespace Eigen;
typedef SparseMatrix<float> SpMat;
typedef Triplet<float> Trip;

class FrameOpt
{
public:
    FrameOpt(const Mesh& mesh, const FrameField& framefield, const PolyLines& polyLines,
            const double alpha = 1000000.0, const double beta = 1.0, const double gamma = 1.0);
    virtual ~FrameOpt();

private:
    FrameOpt();
    FrameOpt(const FrameOpt& frameOpt);

public:
    void Run(const size_t iters = 1);
    bool Optimize();
    bool OptimizeFrame();
    void SetAlpha(const double value = 100.0);
    void SetBeta(const double value = 100.0);
    void SetGamma(const double value = 1.0);
    void SetUseAverageTargetLength(bool value = true);
    void SetRecoverable(bool value = true);
    void SetAllowBigStep(bool value = true);
    void SetAnisotropy(const double value = 0.05);
    void SetStepSize(const double value = 1.0);
    void ComputeMeshTargetLength();
    void UpdateMeshFromFrameField();
    void UpdateFrameFieldFromMesh();

private:
    void OptimizeFrameOrthogonality(std::vector<Trip>& coefficients, std::vector<float>& B, size_t& m_id);  // * beta
    void OptimizeFrameStraightness(std::vector<Trip>& coefficients, std::vector<float>& B, size_t& m_id);   // * gamma
    void OptimizeFrameBoundary(std::vector<Trip>& coefficients, std::vector<float>& B, size_t& m_id);       // * alpha
    void OptimizeFrameInnerVertices(std::vector<Trip>& coefficients, std::vector<float>& B, size_t& m_id);  // * alpha
    void OptimizeOrthogonality(std::vector<Trip>& coefficients, std::vector<float>& B, size_t& m_id);
    void OptimizeStraightness(std::vector<Trip>& coefficients, std::vector<float>& B, size_t& m_id);
    void OptimizeBoundary(std::vector<Trip>& coefficients, std::vector<float>& B, size_t& m_id);

    void OutputBadCells(const std::vector<size_t>& badCellIds, const char* filename);
    void OutputFramesOfBadCells(const std::vector<size_t>& badCellIds, const char* filename);

private:
    Mesh& mesh;
    FrameField& framefield;
    PolyLines& polyLines;

    double alpha;
    double beta;
    double gamma;
    double anisotropy;
    double stepSize;
    double avgMeshEdgeLength;
    double avgFrameEdgeLength;
    bool useAverageTargetLength;
    bool recoverable;
    bool allowBigStep;
    size_t m_numOfInvertdElements;
};

#endif /* FRAME_OPT_H_ */
