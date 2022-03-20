#include "FeatureExtractor.h"
#include "ParallelFor.h"

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
    fe->SetFeatureAngle(angle_threshold);
    fe->SetInputData(mesh_polyData);

    vtkSmartPointer<vtkPointLocator> pl = vtkSmartPointer<vtkPointLocator>::New();
    pl->AutomaticOn();
    pl->SetDataSet(mesh_polyData);
    pl->BuildLocator();

    
    fe->FeatureEdgesOn();
    fe->BoundaryEdgesOff();
    fe->Update();
    
    vtkSmartPointer<vtkPoints> res = fe->GetOutput()->GetPoints();
    int n = res->GetNumberOfPoints();
    if (n > 0) {
        std::cout << "Feature points: " << n << std::endl;
        PARALLEL_FOR_BEGIN(n) {
            SetFeatures(i, res, pl, false);
        } PARALLEL_FOR_END();
    }

    fe->BoundaryEdgesOn();
    fe->FeatureEdgesOff();
    fe->Update();

    res = fe->GetOutput()->GetPoints();
    n = res->GetNumberOfPoints();
    if (n > 0) {
        std::cout << "Boundary Points: " << n << std::endl;
        PARALLEL_FOR_BEGIN(n) {
            SetFeatures(i, res, pl, true);
        } PARALLEL_FOR_END();
    }
}

void FeatureExtractor::SetFeatures(int i, vtkSmartPointer<vtkPoints> res, vtkSmartPointer<vtkPointLocator> pl, bool boundary) {
    double point[3];
    res->GetPoint(i, point);
    vtkIdType id = pl->FindClosestPoint(point);
    if (boundary) {
        mesh.V.at(id).isBoundary = true;
    } else {
        mesh.V.at(id).type = FEATURE;
    }
}
