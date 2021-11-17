#include <vtkNew.h>
#include <vtkSmartPointer.h>
#include <vtkPolygon.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkGenericCell.h>
#include <vtkCellArray.h>
#include <vtkMatrix4x4.h>
#include <vtkPolyDataNormals.h>
#include <vtkImplicitPolyDataDistance.h>
#include <vtkDistancePolyDataFilter.h>
#include <vtkTriangleFilter.h>

#include "QuadSurfaceMapper.h"

SurfaceMapper::SurfaceMapper() {}

SurfaceMapper::SurfaceMapper(Mesh& source_, Mesh& target_) : source(source_), target(target_) {
    source_polyData = GetPolyDataFromMesh(source);
    target_polyData = GetPolyDataFromMesh(target);
}

SurfaceMapper::SurfaceMapper(Mesh& target_) : target(target_) {
    target_polyData = GetPolyDataFromMesh(target);
    SetCellLocator();
}

SurfaceMapper::~SurfaceMapper() {}

void SurfaceMapper::SetSource(Mesh& mesh) {
    source = mesh;
    source_polyData = GetPolyDataFromMesh(mesh);
}

void SurfaceMapper::SetTarget(Mesh& mesh) {
    target = mesh;
    target_polyData = GetPolyDataFromMesh(mesh);
    SetCellLocator();
}

void SurfaceMapper::SetCellLocator() {
    auto tf = vtkSmartPointer<vtkTriangleFilter>::New();
    tf->PassVertsOff();
    tf->PassLinesOff();
    tf->SetInputData(target_polyData);
    tf->Update();

    target_cellPolyData = tf->GetOutput();
    target_cellPolyData->BuildLinks();


    int N = target_cellPolyData->GetNumberOfCells() * 0.1;

    cell_locator->SetDataSet(target_cellPolyData);
    cell_locator->SetTolerance(1e-12);
    cell_locator->SetNumberOfCellsPerBucket(N);
    cell_locator->CacheCellBoundsOn();
    cell_locator->AutomaticOn();
    cell_locator->BuildLocator();

}

vtkSmartPointer<vtkPolyData> SurfaceMapper::GetPolyDataFromMesh(Mesh& mesh) {
    mu.SetMesh(mesh);
    return mu.GetPolyData();
}

glm::dvec3 SurfaceMapper::GetClosestPoint(glm::dvec3 p) {
    double point[] = {p.x, p.y, p.z};
    double closestPoint[3];

    vtkIdType cellId;
    int subId;
    double dist2;
    cell_locator->FindClosestPoint(point, closestPoint, cellId, subId, dist2);
    
    return glm::dvec3(closestPoint[0], closestPoint[1], closestPoint[2]);
}
