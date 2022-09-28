#include "FeatureExtractor.h"
#include "ParallelFor.h"

#include <math.h> 
#define PI 3.14159265

FeatureExtractor::FeatureExtractor() {}
FeatureExtractor::FeatureExtractor(Mesh& mesh_, double angle_threshold_, MeshUtil& mu_) : mesh(&mesh_), mu(&mu_) {
    angle_threshold = angle_threshold_;
    mesh_polyData = GetPolyDataFromMesh();
}
FeatureExtractor::~FeatureExtractor() {}

void FeatureExtractor::SetMembers(Mesh& mesh_, MeshUtil& mu_) {
    mesh = &mesh_;
    mu = &mu_;
    mesh_polyData = GetPolyDataFromMesh();
}

void FeatureExtractor::CheckValidity() {
    if (mesh == NULL) {
        std::cout << "No mesh to use for Feature Extractor." << std::endl;
        exit(0);
    }
    if (mu == NULL) {
        std::cout << "MeshUtil is not initialized for Feature Extractor." << std::endl;
        exit(0);
    }
    if (mesh->V.size() == 0 || mesh->F.size() == 0 || mesh->C.size() == 0) {
        std::cout << "No mesh to use for Feature Extractor." << std::endl;
        exit(0);
    }
}

vtkSmartPointer<vtkPolyData> FeatureExtractor::GetPolyDataFromMesh() {
    return mu->GetPolyData();
}

void FeatureExtractor::Extract() {
    vtkSmartPointer<vtkFeatureEdges> fe = vtkSmartPointer<vtkFeatureEdges>::New();
    fe->SetFeatureAngle(angle_threshold);
    fe->SetInputData(mesh_polyData);

    vtkSmartPointer<vtkPointLocator> pl = vtkSmartPointer<vtkPointLocator>::New();
    pl->AutomaticOn();
    pl->SetDataSet(mesh_polyData);
    pl->BuildLocator();

    
    // fe->FeatureEdgesOn();
    // fe->BoundaryEdgesOff();
    // fe->Update();
    
    // vtkSmartPointer<vtkPoints> res = fe->GetOutput()->GetPoints();
    // int n = res->GetNumberOfPoints();
    // bool setPlanar = n > 0 ? false : true;
    // if (n > 0) {
    //     std::cout << "Feature points: " << n << std::endl;
    //     PARALLEL_FOR_BEGIN(n) {
    //     // for (int i = 0; i < n; i++) {
    //         SetFeatures(i, res, pl, false);
    //     // }
    //     } PARALLEL_FOR_END();
    // }

    fe->BoundaryEdgesOn();
    fe->FeatureEdgesOff();
    fe->Update();

    vtkSmartPointer<vtkPoints> res = fe->GetOutput()->GetPoints();
    int n = res->GetNumberOfPoints();
    // mesh->isPlanar = setPlanar && n > 0 ? true : false;
    if (n > 0) {
        std::cout << "Boundary Points: " << n << std::endl;
        PARALLEL_FOR_BEGIN(n) {
        // for (int i = 0; i < n; i++) {
            SetFeatures(i, res, pl, true);
        // }
        } PARALLEL_FOR_END();
    }

    PARALLEL_FOR_BEGIN(mesh->V.size()) {
        if (mesh->V.at(i).type == FEATURE || mesh->V.at(i).isBoundary) {
            mesh->SetIdealValence(mesh->V.at(i).id);
        }
    } PARALLEL_FOR_END();
    
    // for (int i = 0; i < mesh->V.size(); i++) {
    //     if (mesh->V.at(i).N_Vids.size() != 4) continue;
    //     if (mesh->V.at(i).type == FEATURE || mesh->V.at(i).isBoundary) {
    //         SetIdealValence(mesh->V.at(i).id);
    //         break;
    //     }
    // }
}

void FeatureExtractor::SetFeatures(int i, vtkSmartPointer<vtkPoints> res, vtkSmartPointer<vtkPointLocator> pl, bool boundary) {
    double point[3];
    res->GetPoint(i, point);
    vtkIdType id = pl->FindClosestPoint(point);
    if (boundary) {
        mesh->V.at(id).isBoundary = true;
    } else {
        mesh->V.at(id).type = FEATURE;
    }
}
