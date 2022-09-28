/*
 * MeshOpt.h
 *
 *  Created on: Nov 25, 2016
 *      Author: cotrik
 */

#ifndef MESH_OPT_H_
#define MESH_OPT_H_

#include "Mesh.h"

#include <Eigen/Core>
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Cholesky>
#include <Eigen/LU>
#include <Eigen/SVD>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#include <Eigen/SparseCholesky>
using namespace Eigen;
typedef Eigen::SparseMatrix<float> SpMat;
typedef Eigen::Triplet<float> Trip;

class MeshOpt
{
public:
    MeshOpt(const Mesh& mesh);
    virtual ~MeshOpt();

protected:
    MeshOpt();
    MeshOpt(const MeshOpt& meshOpt);

public:
    virtual size_t Run(const size_t iters = 1);
    virtual bool Optimize();
    void SetTriMesh(Mesh& triMesh);
    void SetRefMesh(Mesh& refMesh);
    void SetTargetSurfaceMesh(Mesh& refMesh);
    void SetAlpha(const double value = 100.0);
    void SetBeta(const double value = 100.0);
    void SetGamma(const double value = 1.0);
    void SetUseAverageTargetLength(bool value = true);
    void SetRecoverable(bool value = true);
    void SetAllowBigStep(bool value = true);
    void SetProjectToTargetSurface(bool value = true);
    void SetUseProjection(bool value = false);
    void SetChangeBoundary(bool value = true);
    void SetStepSize(const double value = 1.0);
    void SetAnisotropy(const double value = 0.05);
    void SetMinScaledJacobian(const double value = 0.00);
    void ComputeMeshTargetLength();
    void AdjustOutlayer();
    void ConvertSurfaceToTriangleMesh();
    void ComputeAvgEdgeLength();

protected:
    void OutputBadCells(const std::vector<size_t>& badCellIds, const char* filename);
    std::vector<size_t> OutputBadCellsAndExtendCells(const std::vector<size_t>& badCellIds, const char* filename);
    std::vector<size_t> GetBadCellsAndExtendCells(const std::vector<size_t>& badCellIds, int extendLayers = 2);
    void OutputEdgesOfBadCells(const std::vector<size_t>& badCellIds, const char* filename);
    void OutputFramesOfBadCells(const std::vector<size_t>& badCellIds, const char* filename);

    void OptimizeBoundaryVertices(std::vector<Trip>& coefficients, std::vector<float>& B, size_t& m_id);   // * alpha
    size_t OptimizeSurfaceVertices(std::vector<Trip>& coefficients, std::vector<float>& B, size_t& m_id);    // * alpha
    void OptimizeEdgeLength(std::vector<Trip>& coefficients, std::vector<float>& B, size_t& m_id);         // * gamma
    void OptimizeInnerVertices(std::vector<Trip>& coefficients, std::vector<float>& B, size_t& m_id);      // * alpha
    void OptimizeEdgeOrthogonality(std::vector<Trip>& coefficients, std::vector<float>& B, size_t& m_id);  // * gamma
    void OptimizeEdgeStraightness(std::vector<Trip>& coefficients, std::vector<float>& B, size_t& m_id);   // * gamma
    void OptimizeSingularity(std::vector<Trip>& coefficients, std::vector<float>& B, size_t& m_id);        // * gamma
    void OptimizeEdgeComformalty(std::vector<Trip>& coefficients, std::vector<float>& B, size_t& m_id);    // * gamma
    void OptimizeFacePlanarity(std::vector<Trip>& coefficients, std::vector<float>& B, size_t& m_id);      // * beta
    void OptimizeScaledJacobian(std::vector<Trip>& coefficients, std::vector<float>& B, size_t& m_id);     // * beta
    void OptimizeSmoothness(std::vector<Trip>& coefficients, std::vector<float>& B, size_t& m_id);         // * gamma
public:
    Mesh bestmesh;

protected:
    Mesh& mesh;
    Mesh* m_refMesh;
    Mesh triMesh;
    Mesh* m_targetSurfaceMesh = NULL;

    double alpha;
    double beta;
    double gamma;
    double anisotropy;
    double minScaledJacobian;

    double stepSize;
    double avgMeshEdgeLength;
    bool useAverageTargetLength;
    bool recoverable;
    bool allowBigStep;
    bool changeBoundary;
    bool projectToTargetSurface = false;
    bool useProjection = true;
    size_t m_numOfInvertdElements;

    std::vector<double> ESingularity;
    std::vector<double> EOrthogonality;
    std::vector<double> EStraightness;
};

bool IsInFace(const Mesh& mesh, const Edge& edge, const size_t vid1, const size_t vid2);
#endif /* MESH_OPT_H_ */
