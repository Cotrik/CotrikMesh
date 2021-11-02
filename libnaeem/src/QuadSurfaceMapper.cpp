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
    vtkNew<vtkPoints> points;
    vtkNew<vtkCellArray> cells;
    auto polyData = vtkSmartPointer<vtkPolyData>::New();

    for (auto& v: mesh.V) {
        points->InsertNextPoint(v.x, v.y, v.z);
    }
    for (auto& c: mesh.C) {
        vtkNew<vtkPolygon> p;
        p->GetPointIds()->SetNumberOfIds(c.Vids.size());
        for (int i = 0; i < c.Vids.size(); i++) {
            p->GetPointIds()->SetId(i, c.Vids.at(i));
        }
        cells->InsertNextCell(p);
    }
    
    polyData->SetPoints(points);
    polyData->SetPolys(cells);

    return polyData;
}

void SurfaceMapper::Map() {
    // auto p1 = GetPolyDataFromMesh(source);
    // auto p2 = GetPolyDataFromMesh(target);

    auto nf = vtkSmartPointer<vtkPolyDataNormals>::New();
    nf->SetInputData(source_polyData);
    nf->ComputePointNormalsOn();
    nf->ComputeCellNormalsOff();
    nf->SplittingOff();
    nf->FlipNormalsOff();
    nf->ConsistencyOn();
    nf->AutoOrientNormalsOn();
    nf->Update();

    auto df = vtkSmartPointer<vtkDistancePolyDataFilter>::New();
    df->SetInputDataObject(0, source_polyData);
    df->SetInputDataObject(1, target_polyData);
    df->Update();

    auto points = source_polyData->GetPointData();
    auto normalsData = nf->GetOutput();
    auto normals = normalsData->GetPointData()->GetNormals();
    auto distanceData = df->GetOutput();
    auto scalars = distanceData->GetPointData()->GetScalars();

    vtkIdType numPoints = source_polyData->GetNumberOfPoints();
    for (vtkIdType i = 0; i < numPoints; i++) {
        double p[3];
        source_polyData->GetPoint(i, p);
        double n[3];
        normals->GetTuple(i, n);
        double d = scalars->GetTuple1(i);

        glm::dvec3 point(p[0], p[1], p[2]);
        glm::dvec3 normal = glm::normalize(glm::dvec3(n[0], n[1], n[2]));

        glm::dvec3 updatedPoint = point + (normal * -d);
        source_polyData->GetPoints()->SetPoint(i, updatedPoint.x, updatedPoint.y, updatedPoint.z);
        source.V.at(i) = updatedPoint;
        // std::cout << source.V.at(i).x << " " << source.V.at(i).y << " " << source.V.at(i).z << " " << std::endl;
    }
}

glm::dvec3 SurfaceMapper::MapPoint(glm::dvec3 p) {
    double closestPoint[3];
    double point[] = {p.x, p.y, p.z};
    auto df = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
    df->SetInput(target_polyData);
    df->EvaluateFunctionAndGetClosestPoint(point, closestPoint);

    return glm::dvec3(closestPoint[0], closestPoint[1], closestPoint[2]);
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
