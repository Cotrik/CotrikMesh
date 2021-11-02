#include "SimplificationOperation.h"

SimplificationOperation::SimplificationOperation() {}

SimplificationOperation::SimplificationOperation(Mesh& mesh_) : mesh(mesh_) {}

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

