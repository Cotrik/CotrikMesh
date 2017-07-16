#ifndef MESH_QUALITY_H
#define MESH_QUALITY_H
#include "Mesh.h"
size_t GetQuality(const char* filename,
        double& minValue, double& maxValue, double& avgValue, std::vector<size_t>& badCellIds,
        const bool output = true, const double minSJ = 0.0);
size_t GetQuality(const char* filename, double& minValue, double& avgValue, const double minSJ = 0.0);
size_t GetScaledJacobianVerdict(const Mesh& mesh, std::vector<double>& scaledJacobian, const double minSJ = 0.0);
size_t GetScaledJacobianVerdict(const Mesh& mesh, double& MinScaledJacobian, double& AvgScaledJacobian, const double minSJ = 0.0);
size_t GetMinScaledJacobian(const Mesh& mesh, double& MinScaledJacobian, std::vector<size_t>& badCellIds, const double minSJ = 0.0);
size_t GetMinScaledJacobianVerdict(const Mesh& mesh, double& MinScaledJacobian, const double minSJ = 0.0);
size_t GetMinScaledJacobianVerdict(const Mesh& mesh, double& MinScaledJacobian, std::vector<size_t>& badCellIds, const double minSJ = 0.0);
size_t GetMinScaledJacobianVerdict(const Mesh& mesh, double& MinScaledJacobian, double& AvgScaledJacobian, std::vector<size_t>& badCellIds, const double minSJ = 0.0);
//double GetMinScaledJacobian(const Mesh& mesh);
//const float GetScaledJacobian(const Mesh& mesh, const Cell& c);
size_t GetQuality(const char* filename, std::vector<size_t>& badCellIds, std::vector<size_t>& warningCellIds,
        std::vector<size_t>& goodCellIds, std::vector<size_t>& highCellIds, std::vector<size_t>& excellentCellIds);

size_t GetMinScaledJacobianQuad(const Mesh& mesh, double& MinScaledJacobian, const double minSJ = 0.0);

#endif // MESH_QUALITY_H
