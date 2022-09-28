/*
 * MeshOpt.h
 *
 *  Created on: Nov 25, 2016
 *      Author: cotrik
 */

#ifndef LAYER_OPT_H_
#define LAYER_OPT_H_

#include "Mesh.h"

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

class LayerOpt
{
public:
    LayerOpt(const Mesh& mesh,
            const double alpha = 100.0, const double beta = 1.0, const double gamma = 1.0);
    virtual ~LayerOpt();

private:
    LayerOpt();
    LayerOpt(const LayerOpt& meshOpt);

public:
    void Run(const size_t iters = 1);
    bool Optimize();
    void SetAlpha(const double value = 100.0);
    void SetBeta(const double value = 100.0);
    void SetGamma(const double value = 1.0);
    void SetUseAverageTargetLength(bool value = true);
    void SetRecoverable(bool value = true);
    void SetStepSize(const double value = 1.0);
    void SetAnisotropy(const double value = 0.05);
    void ComputeMeshTargetLength();

private:
    void OptimizeBoundaryVertices(std::vector<Trip>& coefficients, std::vector<float>& B, size_t& m_id);   // * alpha
    void OptimizeEdgeLength(std::vector<Trip>& coefficients, std::vector<float>& B, size_t& m_id);         // * gamma
    void OptimizeInnerVertices(std::vector<Trip>& coefficients, std::vector<float>& B, size_t& m_id);      // * alpha
    void OptimizeEdgeOrthogonality(std::vector<Trip>& coefficients, std::vector<float>& B, size_t& m_id);  // * beta
    void OptimizeEdgeStraightness(std::vector<Trip>& coefficients, std::vector<float>& B, size_t& m_id);   // * gamma
    void OptimizeSingularity(std::vector<Trip>& coefficients, std::vector<float>& B, size_t& m_id);        // * gamma
    void OptimizeEdgeComformalty(std::vector<Trip>& coefficients, std::vector<float>& B, size_t& m_id);    // * gamma

    void OutputBadCells(const std::vector<size_t>& badCellIds, const char* filename);
    void OutputFramesOfBadCells(const std::vector<size_t>& badCellIds, const char* filename);

private:
    Mesh& mesh;

    double alpha;
    double beta;
    double gamma;
    double anisotropy;

    double stepSize;
    double avgMeshEdgeLength;
    bool useAverageTargetLength;
    bool recoverable;
    size_t m_numOfInvertdElements;

    std::vector<double> ESingularity;
    std::vector<double> EOrthogonality;
    std::vector<double> EStraightness;
};

#endif /* LAYER_OPT_H_ */
