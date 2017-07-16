#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkDataSet.h>
#include <vtkMeshQuality.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkTriangleFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkProperty.h>
#include <vtkScalarBarActor.h>
#include <vtkLookupTable.h>
#include <vtkColorTransferFunction.h>
#include <vtkThreshold.h>
#include <vtkExtractEdges.h>

#include "MeshQuality.h"
void GetMetrics(vtkMeshQuality* qualityFilter, double& minValue, double& maxValue, double& avgValue)
{
    vtkSmartPointer<vtkDoubleArray> qualityArray = vtkDoubleArray::SafeDownCast(qualityFilter->GetOutput()->GetCellData()->GetArray("Quality"));
    minValue = qualityArray->GetValue(0);
    maxValue = qualityArray->GetValue(0);
    avgValue = 0;
    for (vtkIdType i = 0; i < qualityArray->GetNumberOfTuples(); i++)
    {
        double val = qualityArray->GetValue(i);
        if (minValue > val) minValue = val;
        if (maxValue < val) maxValue = val;
        avgValue += val;
        //std::cout << "value " << i << " : " << val << std::endl;
    }
    avgValue /= qualityArray->GetNumberOfTuples();
}

size_t GetQuality(const char* filename,
        double& minValue, double& maxValue, double& avgValue, std::vector<size_t>& badCellIds,
        const bool output, const double minSJ)
{
    vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName(filename);
    reader->Update();
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGridMesh = reader->GetOutput();
    vtkSmartPointer<vtkMeshQuality> qualityFilter = vtkSmartPointer<vtkMeshQuality>::New();
#if VTK_MAJOR_VERSION <= 5
    qualityFilter->SetInputConnection(mesh->GetProducerPort());
#else
    ///qualityFilter->SetInputData(mesh);
    qualityFilter->SetInputConnection(reader->GetOutputPort());
#endif
    minValue = 0;
    maxValue = 0;
    avgValue = 0;
    std::vector<double> metrics;
    //qualityFilter->SetTriangleQualityMeasureToArea();
    //qualityFilter->SetHexQualityMeasureToEdgeRatio()
    //qualityFilter->SetHexQualityMeasureToMedAspectFrobenius()
    //qualityFilter->SetHexQualityMeasureToMaxAspectFrobenius()
    //qualityFilter->SetHexQualityMeasureToMaxEdgeRatios()
    //qualityFilter->SetHexQualityMeasureToSkew()
    //qualityFilter->SetHexQualityMeasureToTaper()
    //qualityFilter->SetHexQualityMeasureToVolume()
    //qualityFilter->SetHexQualityMeasureToStretch()
    //qualityFilter->SetHexQualityMeasureToDiagonal()
    //qualityFilter->SetHexQualityMeasureToDimension()
    //qualityFilter->SetHexQualityMeasureToOddy()
    //qualityFilter->SetHexQualityMeasureToCondition()
    //qualityFilter->SetHexQualityMeasureToJacobian()
    qualityFilter->SetHexQualityMeasureToScaledJacobian();
    //qualityFilter->SetHexQualityMeasureToShear()
    //qualityFilter->SetHexQualityMeasureToShape()
    //qualityFilter->SetHexQualityMeasureToRelativeSizeSquared()
    //qualityFilter->SetHexQualityMeasureToShapeAndSize()
    //qualityFilter->SetHexQualityMeasureToShearAndSize()
    //qualityFilter->SetHexQualityMeasureToDistortion()
    qualityFilter->Update();

    vtkSmartPointer<vtkDoubleArray> qualityArray = vtkDoubleArray::SafeDownCast(qualityFilter->GetOutput()->GetCellData()->GetArray("Quality"));
    minValue = qualityArray->GetValue(0);
    maxValue = qualityArray->GetValue(0);
    avgValue = 0;
    size_t numOfInvertedElements = 0;
    for (vtkIdType i = 0; i < qualityArray->GetNumberOfTuples(); i++)
    {
        double val = qualityArray->GetValue(i);
        if (minValue > val) minValue = val;
        if (maxValue < val) maxValue = val;
        avgValue += val;

        if (val < minSJ) {
            numOfInvertedElements++;
            badCellIds.push_back(i);
        }
        else if (val > 1) numOfInvertedElements++;
        //std::cout << "value " << i << " : " << val << std::endl;
    }
    avgValue /= qualityArray->GetNumberOfTuples();
    if (output) {
        std::cout << "------------------------------------ " << std::endl;
        std::cout << "\033[1;36mmesh = " << filename << "\033[0m"<< std::endl;
        std::cout << "\033[1;31mmin scaled jacobian\033[0m = " << minValue << std::endl;
        //std::cout << "max scaled jacobian = " << maxValue << std::endl;
        //std::cout << "avg scaled jacobian = " << avgValue << std::endl;
        std::cout << "\033[1;31m#InvertedElements\033[0m = " << numOfInvertedElements << std::endl;
        std::cout << "------------------------------------ " << std::endl;
    }
    return numOfInvertedElements;
}

size_t GetQuality(const char* filename, double& minValue, double& avgValue, const double minSJ/* = 0.0*/)
{
    vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName(filename);
    reader->Update();
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGridMesh = reader->GetOutput();
    vtkSmartPointer<vtkMeshQuality> qualityFilter = vtkSmartPointer<vtkMeshQuality>::New();
#if VTK_MAJOR_VERSION <= 5
    qualityFilter->SetInputConnection(mesh->GetProducerPort());
#else
    ///qualityFilter->SetInputData(mesh);
    qualityFilter->SetInputConnection(reader->GetOutputPort());
#endif
    minValue = 0;
    double maxValue = 0;
    avgValue = 0;
    std::vector<double> metrics;
    //qualityFilter->SetTriangleQualityMeasureToArea();
    //qualityFilter->SetHexQualityMeasureToEdgeRatio()
    //qualityFilter->SetHexQualityMeasureToMedAspectFrobenius()
    //qualityFilter->SetHexQualityMeasureToMaxAspectFrobenius()
    //qualityFilter->SetHexQualityMeasureToMaxEdgeRatios()
    //qualityFilter->SetHexQualityMeasureToSkew()
    //qualityFilter->SetHexQualityMeasureToTaper()
    //qualityFilter->SetHexQualityMeasureToVolume()
    //qualityFilter->SetHexQualityMeasureToStretch()
    //qualityFilter->SetHexQualityMeasureToDiagonal()
    //qualityFilter->SetHexQualityMeasureToDimension()
    //qualityFilter->SetHexQualityMeasureToOddy()
    //qualityFilter->SetHexQualityMeasureToCondition()
    //qualityFilter->SetHexQualityMeasureToJacobian()
    qualityFilter->SetHexQualityMeasureToScaledJacobian();
    //qualityFilter->SetHexQualityMeasureToShear()
    //qualityFilter->SetHexQualityMeasureToShape()
    //qualityFilter->SetHexQualityMeasureToRelativeSizeSquared()
    //qualityFilter->SetHexQualityMeasureToShapeAndSize()
    //qualityFilter->SetHexQualityMeasureToShearAndSize()
    //qualityFilter->SetHexQualityMeasureToDistortion()
    qualityFilter->Update();

    vtkSmartPointer<vtkDoubleArray> qualityArray = vtkDoubleArray::SafeDownCast(qualityFilter->GetOutput()->GetCellData()->GetArray("Quality"));
    minValue = qualityArray->GetValue(0);
    maxValue = qualityArray->GetValue(0);
    avgValue = 0;
    size_t numOfInvertedElements = 0;
    for (vtkIdType i = 0; i < qualityArray->GetNumberOfTuples(); i++)
    {
        double val = qualityArray->GetValue(i);
        if (minValue > val) minValue = val;
        if (maxValue < val) maxValue = val;
        avgValue += val;

        if (val < minSJ) {
            numOfInvertedElements++;
            //badCellIds.push_back(i);
        }
        else if (val > 1) numOfInvertedElements++;
        //std::cout << "value " << i << " : " << val << std::endl;
    }
    avgValue /= qualityArray->GetNumberOfTuples();
    return numOfInvertedElements;
}
size_t GetQuality(const char* filename, std::vector<size_t>& badCellIds, std::vector<size_t>& warningCellIds,
        std::vector<size_t>& goodCellIds, std::vector<size_t>& highCellIds, std::vector<size_t>& excellentCellIds)
{
    vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName(filename);
    reader->Update();
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGridMesh = reader->GetOutput();
    vtkSmartPointer<vtkMeshQuality> qualityFilter = vtkSmartPointer<vtkMeshQuality>::New();
#if VTK_MAJOR_VERSION <= 5
    qualityFilter->SetInputConnection(mesh->GetProducerPort());
#else
    ///qualityFilter->SetInputData(mesh);
    qualityFilter->SetInputConnection(reader->GetOutputPort());
#endif
    double minValue = 0;
    double maxValue = 0;
    double avgValue = 0;
    std::vector<double> metrics;
    //qualityFilter->SetTriangleQualityMeasureToArea();
    //qualityFilter->SetHexQualityMeasureToEdgeRatio()
    //qualityFilter->SetHexQualityMeasureToMedAspectFrobenius()
    //qualityFilter->SetHexQualityMeasureToMaxAspectFrobenius()
    //qualityFilter->SetHexQualityMeasureToMaxEdgeRatios()
    //qualityFilter->SetHexQualityMeasureToSkew()
    //qualityFilter->SetHexQualityMeasureToTaper()
    //qualityFilter->SetHexQualityMeasureToVolume()
    //qualityFilter->SetHexQualityMeasureToStretch()
    //qualityFilter->SetHexQualityMeasureToDiagonal()
    //qualityFilter->SetHexQualityMeasureToDimension()
    //qualityFilter->SetHexQualityMeasureToOddy()
    //qualityFilter->SetHexQualityMeasureToCondition()
    //qualityFilter->SetHexQualityMeasureToJacobian()
    qualityFilter->SetHexQualityMeasureToScaledJacobian();
    //qualityFilter->SetHexQualityMeasureToShear()
    //qualityFilter->SetHexQualityMeasureToShape()
    //qualityFilter->SetHexQualityMeasureToRelativeSizeSquared()
    //qualityFilter->SetHexQualityMeasureToShapeAndSize()
    //qualityFilter->SetHexQualityMeasureToShearAndSize()
    //qualityFilter->SetHexQualityMeasureToDistortion()
    qualityFilter->Update();

    vtkSmartPointer<vtkDoubleArray> qualityArray = vtkDoubleArray::SafeDownCast(qualityFilter->GetOutput()->GetCellData()->GetArray("Quality"));
    minValue = qualityArray->GetValue(0);
    maxValue = qualityArray->GetValue(0);
    avgValue = 0;
    size_t numOfInvertedElements = 0;
    for (vtkIdType i = 0; i < qualityArray->GetNumberOfTuples(); i++)
    {
        double val = qualityArray->GetValue(i);
        if (minValue > val) minValue = val;
        if (maxValue < val) maxValue = val;
        avgValue += val;

        if (val <= 0) {
            numOfInvertedElements++;
            badCellIds.push_back(i);
        }
        else if (val > 1 || val < -1)  {
            numOfInvertedElements++;
            badCellIds.push_back(i);
        }
        else if (val > 0 && val < 0.2) {
            warningCellIds.push_back(i);
        }
        else if (val >= 0.2 && val < 0.4) {
            goodCellIds.push_back(i);
        }
        else if (val >= 0.4 && val < 0.6) {
            highCellIds.push_back(i);
        }
        else if (val >= 0.6 && val <= 1) {
            excellentCellIds.push_back(i);
        }
    }
//    avgValue /= qualityArray->GetNumberOfTuples();
//    if (output) {
//        std::cout << "------------------------------------ " << std::endl;
//        std::cout << "\033[1;36mmesh = " << filename << "\033[0m"<< std::endl;
//        std::cout << "\033[1;31mmin scaled jacobian\033[0m = " << minValue << std::endl;
//        //std::cout << "max scaled jacobian = " << maxValue << std::endl;
//        //std::cout << "avg scaled jacobian = " << avgValue << std::endl;
//        std::cout << "\033[1;31m#InvertedElements\033[0m = " << numOfInvertedElements << std::endl;
//        std::cout << "------------------------------------ " << std::endl;
//    }
    return numOfInvertedElements;
}

const int V_T[8][4] =
{
    {0, 3, 4, 1},
    {1, 0, 5, 2},
    {2, 1, 6, 3},
    {3, 2, 7, 0},
    {4, 7, 5, 0},
    {5, 4, 6, 1},
    {6, 5, 7, 2},
    {7, 6, 4, 3}
};

static float cal_volume_Tet_real(float v0[3],float v1[3],float v2[3],float v3[3])
{
    float v1v0[3], v2v0[3], v3v0[3];
    for (int i = 0; i < 3; i++) {
        v1v0[i] = v1[i] - v0[i];
        v2v0[i] = v2[i] - v0[i];
        v3v0[i] = v3[i] - v0[i];
    }

    float norm1 = sqrt(v1v0[0] * v1v0[0] + v1v0[1] * v1v0[1] + v1v0[2] * v1v0[2]);
    float norm2 = sqrt(v2v0[0] * v2v0[0] + v2v0[1] * v2v0[1] + v2v0[2] * v2v0[2]);
    float norm3 = sqrt(v3v0[0] * v3v0[0] + v3v0[1] * v3v0[1] + v3v0[2] * v3v0[2]);

    float volume = v1v0[0] * (v2v0[1] * v3v0[2] - v2v0[2] * v3v0[1]) - v1v0[1] * (v2v0[0] * v3v0[2] - v2v0[2] * v3v0[0]) + v1v0[2] * (v2v0[0] * v3v0[1] - v2v0[1] * v3v0[0]);
    return volume;
}

static bool JudgeDirection(const Mesh& mesh, const Cell& c)
{
    float v[8][3];
    for (int i = 0; i < c.Vids.size(); i++) {
        v[i][0] = mesh.V.at(c.Vids.at(i)).x;
        v[i][1] = mesh.V.at(c.Vids.at(i)).y;
        v[i][2] = mesh.V.at(c.Vids.at(i)).z;
    }

    float VL[8];
    for (int i = 0; i < 8; i++) {
        const int* p = V_T[i];
        VL[i] = cal_volume_Tet_real(v[p[0]], v[p[1]], v[p[2]], v[p[3]]);
    }

    if (VL[0] + VL[1] + VL[2] + VL[3] + VL[4] + VL[5] + VL[6] + VL[7] < 0) return false;
    else
        return true;
}

static float _GetScaledJacobian(const glm::vec3& i, const glm::vec3& j, const glm::vec3& k)
{
    const glm::mat3x3 m(i, j, k);
    return glm::determinant(m);
}

static const float GetScaledJacobian(const Mesh& mesh, const Cell& c)
{
    Cell c1(c);
    if (!JudgeDirection(mesh, c)) {
        std::swap(c1.Vids[0], c1.Vids[3]);
        std::swap(c1.Vids[1], c1.Vids[2]);
        std::swap(c1.Vids[4], c1.Vids[7]);
        std::swap(c1.Vids[5], c1.Vids[6]);
    }
    float minScaledJacobian = 1;
    for (int n = 0; n < 7; n++) {
        const Vertex& o = mesh.V.at(c1.Vids.at(n));

        const Vertex& i = mesh.V.at(c1.Vids.at(HexPoint_Points[n][0]));
        const Vertex& j = mesh.V.at(c1.Vids.at(HexPoint_Points[n][1]));
        const Vertex& k = mesh.V.at(c1.Vids.at(HexPoint_Points[n][2]));

        const glm::vec3 ei(i.x - o.x, i.y - o.y, i.z - o.z);
        const glm::vec3 ej(j.x - o.x, j.y - o.y, j.z - o.z);
        const glm::vec3 ek(k.x - o.x, k.y - o.y, k.z - o.z);

        const float length_i = glm::length(ei);
        const float length_j = glm::length(ej);
        const float length_k = glm::length(ek);

        const glm::vec3 ni(ei.x / length_i, ei.y / length_i, ei.z / length_i);
        const glm::vec3 nj(ej.x / length_j, ej.y / length_j, ej.z / length_j);
        const glm::vec3 nk(ek.x / length_k, ek.y / length_k, ek.z / length_k);

        float scaledJacobian = _GetScaledJacobian(ni, nj, nk);
        minScaledJacobian = minScaledJacobian < scaledJacobian ? minScaledJacobian : scaledJacobian;
    }

    return minScaledJacobian;
}

size_t GetMinScaledJacobian(const Mesh& mesh, double& MinScaledJacobian, std::vector<size_t>& badCellIds, const double minSJ/* = 0.0*/)
{
    size_t numOfInvertedElements = 0;
    badCellIds.clear();
    float minScaledJacobian = 1;
    for (int i = 0; i < mesh.C.size(); i++) {
        float scaledJacobian = GetScaledJacobian(mesh, mesh.C.at(i));
        minScaledJacobian = minScaledJacobian < scaledJacobian ? minScaledJacobian : scaledJacobian;
        if (scaledJacobian < minSJ) {
            numOfInvertedElements++;
            badCellIds.push_back(i);
        }
    }
    return numOfInvertedElements;
}

#include "verdict.h"
size_t GetMinScaledJacobianVerdict(const Mesh& mesh, double& MinScaledJacobian, const double minSJ/* = 0.0*/)
{
    size_t numOfInvertedElements = 0;
    double minScaledJacobian = 1;
    const std::vector<Vertex>& V = mesh.V;
    const std::vector<Cell>& C = mesh.C;
    double coordinates[8][3];
    for (size_t i = 0; i < C.size(); i++) {
        const Cell& c = C[i];
        for (size_t j = 0; j < 8; j++) {
            const Vertex& v = V[c.Vids[j]];
            coordinates[j][0] = v.x;
            coordinates[j][1] = v.y;
            coordinates[j][2] = v.z;
        }
        double scaledJacobian = v_hex_scaled_jacobian(8, coordinates);
        minScaledJacobian = minScaledJacobian < scaledJacobian ? minScaledJacobian : scaledJacobian;
        if (scaledJacobian < minSJ) {
            numOfInvertedElements++;
        }
    }
    MinScaledJacobian = minScaledJacobian;

    return numOfInvertedElements;
}

size_t GetMinScaledJacobianQuad(const Mesh& mesh, double& MinScaledJacobian, const double minSJ/* = 0.0*/)
{
    size_t numOfInvertedElements = 0;
    double minScaledJacobian = 1;
    const std::vector<Vertex>& V = mesh.V;
    const std::vector<Face>& F = mesh.F;
    double coordinates[4][3];
    for (size_t i = 0; i < F.size(); i++) {
        const Face& f = F[i];
        for (size_t j = 0; j < 4; j++) {
            const Vertex& v = V[f.Vids[j]];
            coordinates[j][0] = v.x;
            coordinates[j][1] = v.y;
            coordinates[j][2] = v.z;
        }
        double scaledJacobian = v_quad_scaled_jacobian(4, coordinates);
        minScaledJacobian = minScaledJacobian < scaledJacobian ? minScaledJacobian : scaledJacobian;
        if (scaledJacobian < minSJ) {
            numOfInvertedElements++;
        }
    }
    MinScaledJacobian = minScaledJacobian;

    return numOfInvertedElements;
}
size_t GetScaledJacobianVerdict(const Mesh& mesh, std::vector<double>& scaledJacobian, const double minSJ/* = 0.0*/)
{
    size_t numOfInvertedElements = 0;
    //AvgScaledJacobian = 0;
    //double minScaledJacobian = 1;
    scaledJacobian.resize(mesh.C.size());
    const std::vector<Vertex>& V = mesh.V;
    const std::vector<Cell>& C = mesh.C;
    double coordinates[8][3];
    for (size_t i = 0; i < C.size(); i++) {
        const Cell& c = C[i];
        for (size_t j = 0; j < 8; j++) {
            const Vertex& v = V[c.Vids[j]];
            coordinates[j][0] = v.x;
            coordinates[j][1] = v.y;
            coordinates[j][2] = v.z;
        }
        scaledJacobian[i] = v_hex_scaled_jacobian(8, coordinates);;
        //double scaledJacobian = v_hex_scaled_jacobian(8, coordinates);
        //minScaledJacobian = minScaledJacobian < scaledJacobian ? minScaledJacobian : scaledJacobian;
        if (scaledJacobian[i] < minSJ) {
            numOfInvertedElements++;
        }
        //AvgScaledJacobian +=  scaledJacobian;
    }
    //MinScaledJacobian = minScaledJacobian;
    //AvgScaledJacobian /= C.size();

    return numOfInvertedElements;
}
size_t GetScaledJacobianVerdict(const Mesh& mesh, double& MinScaledJacobian, double& AvgScaledJacobian, const double minSJ/* = 0.0*/)
{
    size_t numOfInvertedElements = 0;
    AvgScaledJacobian = 0;
    double minScaledJacobian = 1;
    const std::vector<Vertex>& V = mesh.V;
    const std::vector<Cell>& C = mesh.C;
    double coordinates[8][3];
    for (size_t i = 0; i < C.size(); i++) {
        const Cell& c = C[i];
        for (size_t j = 0; j < 8; j++) {
            const Vertex& v = V[c.Vids[j]];
            coordinates[j][0] = v.x;
            coordinates[j][1] = v.y;
            coordinates[j][2] = v.z;
        }
        double scaledJacobian = v_hex_scaled_jacobian(8, coordinates);
        minScaledJacobian = minScaledJacobian < scaledJacobian ? minScaledJacobian : scaledJacobian;
        if (scaledJacobian < minSJ) {
            numOfInvertedElements++;
        }
        AvgScaledJacobian +=  scaledJacobian;
    }
    MinScaledJacobian = minScaledJacobian;
    AvgScaledJacobian /= C.size();

    return numOfInvertedElements;
}

size_t GetMinScaledJacobianVerdict(const Mesh& mesh, double& MinScaledJacobian, std::vector<size_t>& badCellIds, const double minSJ/* = 0.0*/)
{
    size_t numOfInvertedElements = 0;
    badCellIds.clear();
    double minScaledJacobian = 1;
    const std::vector<Vertex>& V = mesh.V;
    const std::vector<Cell>& C = mesh.C;
    double coordinates[8][3];
    for (size_t i = 0; i < C.size(); i++) {
        const Cell& c = C[i];
        for (size_t j = 0; j < 8; j++) {
            const Vertex& v = V[c.Vids[j]];
            coordinates[j][0] = v.x;
            coordinates[j][1] = v.y;
            coordinates[j][2] = v.z;
        }
        double scaledJacobian = v_hex_scaled_jacobian(8, coordinates);
        minScaledJacobian = minScaledJacobian < scaledJacobian ? minScaledJacobian : scaledJacobian;
        if (scaledJacobian < minSJ) {
            numOfInvertedElements++;
            badCellIds.push_back(i);
        }
    }
    MinScaledJacobian = minScaledJacobian;

    return numOfInvertedElements;
}

size_t GetMinScaledJacobianVerdict(const Mesh& mesh, double& MinScaledJacobian, double& AvgScaledJacobian, std::vector<size_t>& badCellIds, const double minSJ/* = 0.0*/)
{
    size_t numOfInvertedElements = 0;
    AvgScaledJacobian = 0;
    badCellIds.clear();
    double minScaledJacobian = 1;
    const std::vector<Vertex>& V = mesh.V;
    const std::vector<Cell>& C = mesh.C;
    double coordinates[8][3];
    for (size_t i = 0; i < C.size(); i++) {
        const Cell& c = C[i];
        for (size_t j = 0; j < 8; j++) {
            const Vertex& v = V[c.Vids[j]];
            coordinates[j][0] = v.x;
            coordinates[j][1] = v.y;
            coordinates[j][2] = v.z;
        }
        double scaledJacobian = v_hex_scaled_jacobian(8, coordinates);
        minScaledJacobian = minScaledJacobian < scaledJacobian ? minScaledJacobian : scaledJacobian;
        if (scaledJacobian < minSJ) {
            numOfInvertedElements++;
            badCellIds.push_back(i);
        }
        AvgScaledJacobian +=  scaledJacobian;
    }
    MinScaledJacobian = minScaledJacobian;
    AvgScaledJacobian /= C.size();

    return numOfInvertedElements;
}

static double GetMinScaledJacobian(const Mesh& mesh)
{
    double minScaledJacobian = 1;
    for (int i = 0; i < mesh.C.size(); i++) {
        float scaledJacobian = GetScaledJacobian(mesh, mesh.C.at(i));
        minScaledJacobian = minScaledJacobian < scaledJacobian ? minScaledJacobian : scaledJacobian;
    }
    return minScaledJacobian;
}
