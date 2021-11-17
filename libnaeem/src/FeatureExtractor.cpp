#include "FeatureExtractor.h"
#include <vtkFeatureEdges.h>
#include <vtkPointData.h>
#include <vtkPointLocator.h>

FeatureExtractor::FeatureExtractor() {}
FeatureExtractor::FeatureExtractor(Mesh& mesh_, double angle_threshold_) : mesh(mesh_) {
    angle_threshold = angle_threshold_;
    mesh_polyData = GetPolyDataFromMesh();
}
FeatureExtractor::~FeatureExtractor() {}

void FeatureExtractor::SetMesh(Mesh& mesh_) {
    mesh = mesh_;
    mesh_polyData = GetPolyDataFromMesh();
}

vtkSmartPointer<vtkPolyData> FeatureExtractor::GetPolyDataFromMesh() {
    mu.SetMesh(mesh);
    return mu.GetPolyData();
}

void FeatureExtractor::Extract() {
    vtkSmartPointer<vtkFeatureEdges> fe = vtkSmartPointer<vtkFeatureEdges>::New();
    fe->BoundaryEdgesOn();
    fe->FeatureEdgesOn();
    fe->SetFeatureAngle(angle_threshold);
    fe->SetInputData(mesh_polyData);
    fe->Update();

    int numBucketPoints = mesh_polyData->GetNumberOfPoints() * 0.01;
    vtkSmartPointer<vtkPointLocator> pl = vtkSmartPointer<vtkPointLocator>::New();
    pl->AutomaticOn();
    pl->SetNumberOfPointsPerBucket(numBucketPoints);
    pl->SetDataSet(mesh_polyData);
    pl->BuildLocator();

    auto res = fe->GetOutput()->GetPoints();
    auto n = res->GetNumberOfPoints();

    for (vtkIdType i = 0; i < n; i++) {
        double point[3];
        res->GetPoint(i, point);
        vtkIdType id = pl->FindClosestPoint(point);
        mesh.V.at(id).type = FEATURE;
    }
}
