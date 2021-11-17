#include "SimplificationOperation.h"

SimplificationOperation::SimplificationOperation() {}

SimplificationOperation::SimplificationOperation(Mesh& mesh_, MeshUtil& mu_) : mesh(mesh_), mu(mu_) {}

SimplificationOperation::~SimplificationOperation() {}

void SimplificationOperation::CheckValidity() {
    if (mesh.V.size() == 0 || mesh.F.size() == 0 || mesh.C.size() == 0) {
        std::cout << "No mesh to use for Simplification Operation." << std::endl;
        exit(0);
    }
}

void SimplificationOperation::SetMesh(Mesh& mesh_) {
    mesh = mesh_;
}

std::vector<size_t> SimplificationOperation::GetDifference(std::vector<size_t>& a, std::vector<size_t>& b) {
    return mu.GetDifference(a, b);
}

void SimplificationOperation::AddContents(std::vector<size_t>& a, std::vector<size_t>& b) {
    mu.AddContents(a, b);
}

void SimplificationOperation::UpdateContents(std::vector<size_t>& a, std::vector<size_t>& b) {
    mu.UpdateContents(a, b);
}


