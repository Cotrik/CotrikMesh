#include <algorithm>
#include "SemiGlobalSimplifier.h"

SemiGlobalSimplifier::SemiGlobalSimplifier() {}

SemiGlobalSimplifier::SemiGlobalSimplifier(Mesh& mesh_) : mesh(mesh_) {
    mu.SetMesh(mesh);
}

SemiGlobalSimplifier::~SemiGlobalSimplifier() {}

void SemiGlobalSimplifier::CheckValidity() {
    if (mesh.V.size() == 0 || mesh.F.size() == 0 || mesh.C.size() == 0) {
        std::cout << "No mesh to use for Semi Global Simplifier." << std::endl;
        exit(0);
    }
}

void SemiGlobalSimplifier::SetMesh(Mesh& mesh_) {
    mesh = mesh_;
}

void SemiGlobalSimplifier::SetSimplificationOperations() {
    CheckValidity();
    
    SetDiagonalCollapseOperations();
}

void SemiGlobalSimplifier::SetDiagonalCollapseOperations() {
    CheckValidity();

    for (auto& f: mesh.F) {
        std::unique_ptr<SimplificationOperation> dc1 = std::make_unique<DiagonalCollapse>(mesh, f.id, 0, 2);
        dc1->SetRanking(mu);
        Ops.push_back(std::move(dc1));

        std::unique_ptr<SimplificationOperation> dc2 = std::make_unique<DiagonalCollapse>(mesh, f.id, 1, 3);
        dc2->SetRanking(mu);
        Ops.push_back(std::move(dc2));
    }

    int i = 0;
    for (auto& op: Ops) {
        op->PerformOperation();
        std::cout << op->ranking << std::endl;
        // if (i == 10) {
            break;
        // }
        // i += 1;
    }
}
