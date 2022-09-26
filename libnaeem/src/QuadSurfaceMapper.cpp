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


#include "QuadSurfaceMapper.h"

SurfaceMapper::SurfaceMapper() {}

SurfaceMapper::SurfaceMapper(Mesh& source_, Mesh& target_) : source(source_), target(target_) {
    source_polyData = GetPolyDataFromMesh(source);
    target_polyData = GetPolyDataFromMesh(target);
}

SurfaceMapper::SurfaceMapper(Mesh& target_) : target(target_) {
    target_polyData = GetPolyDataFromMesh(target);
    // point_finder->SetInput(target_polyData);
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
    // point_finder->SetInput(target_polyData);
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


    int N = target_cellPolyData->GetNumberOfCells() > 100000 ? target_cellPolyData->GetNumberOfCells() * 0.001 : 100;

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

void SurfaceMapper::ExecutePolyDistanceFilter(Mesh& mesh) {
    vtkSmartPointer<vtkPolyData> src = GetPolyDataFromMesh(mesh);
    polyDistanceFilter->SetInputData(0, src);
    polyDistanceFilter->SetInputData(1, target_polyData);
    polyDistanceFilter->Update();
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

// double SurfaceMapper::ExecuteHaursdorffDsitanceFilter(Mesh& m1, Mesh& m2) {
//     vtkSmartPointer<vtkPolyData> p1 = GetPolyDataFromMesh(m1);
//     vtkSmartPointer<vtkPolyData> p2 = GetPolyDataFromMesh(m2);

//     haursdorffDistanceFilter->SetInput(0, p1);
//     haursdorffDistanceFilter->SetInput(1, p2);
//     haursdorffDistanceFilter->Update();

//     return haursdorffDistanceFilter->GetHausdorffDistance();
// }

void SurfaceMapper::RemapVertex(Mesh& mesh_, size_t vid, glm::dvec3 c) {
    // vtkSmartPointer<vtkPolyData> poly = vtkSmartPointer<vtkPolyData>::New();
    poly->Reset();
    points->Reset();
    cells->Reset();
    vertices.clear();
    faces = mesh_.V.at(vid).N_Fids;
    double radius = -1.0;
    for (auto fid: faces) {
        auto& f = mesh_.F.at(fid);
        if (f.N_Fids.empty() || f.Vids.empty()) continue;
        mu.AddContents(vertices, std::vector<size_t>(f.Vids.begin(), f.Vids.end()));
        p->GetPointIds()->SetNumberOfIds(f.Vids.size());
        for (int i = 0; i < f.Vids.size(); i++) {
            p->GetPointIds()->SetId(i, f.Vids.at(i));
        }
        cells->InsertNextCell(p);
    }
    for (auto vid: vertices) {
        auto& v = mesh_.V.at(vid);
        points->InsertPoint(vid, v.x, v.y, v.z);
        double r = glm::distance(glm::dvec3(0.0, 0.0, 0.0), v.xyz());
        if (radius < 0 || r < radius) {
            radius = r;
        }
    }
    if (points->GetNumberOfPoints() < 1 || cells->GetNumberOfCells() < 1) return;
    poly->SetPoints(points);
    poly->SetPolys(cells);
    
    // auto tf = vtkSmartPointer<vtkTriangleFilter>::New();
    // tf->PassVertsOff();
    // tf->PassLinesOff();
    // tf->SetInputData(poly);
    // tf->Update();

    // tp = tf->GetOutput();
    // return;
    // tp->BuildLinks();
    // if (poly->GetNumberOfVerts() < 1 || poly->GetNumberOfPolys() < 1) return;

    // vtkSmartPointer<vtkCellLocator> cl = vtkSmartPointer<vtkCellLocator>::New();
    cl->SetDataSet(poly);
    // cl->SetDataSet(poly);
    cl->SetTolerance(1e-12);
    cl->CacheCellBoundsOn();
    cl->AutomaticOn();
    cl->BuildLocator();

    double point[] = {c.x, c.y, c.z};
    double closestPoint[3];

    vtkIdType cellId;
    int subId;
    double dist2;
    auto& v = mesh_.V.at(vid);
    if (cl->FindClosestPointWithinRadius(point, radius, closestPoint, cellId, subId, dist2) == 1) {
        v.x = closestPoint[0];
        v.y = closestPoint[1];
        v.z = closestPoint[2];
    } else {
        v.x = c.x;
        v.y = c.y;
        v.z = c.z;
    }
}

void SurfaceMapper::SetLocator(Mesh& mesh_, std::vector<size_t> V) {
    vtkSmartPointer<vtkPolyData> poly = mu.GetPolyData(mesh_, V);

    auto tf = vtkSmartPointer<vtkTriangleFilter>::New();
    tf->PassVertsOff();
    tf->PassLinesOff();
    tf->SetInputData(poly);
    tf->Update();

    target_cellPolyData = tf->GetOutput();
    target_cellPolyData->BuildLinks();

    cell_locator->SetDataSet(target_cellPolyData);
    // cell_locator->SetDataSet(poly);
    cell_locator->SetTolerance(1e-12);
    cell_locator->CacheCellBoundsOn();
    cell_locator->AutomaticOn();
    cell_locator->BuildLocator();
}

