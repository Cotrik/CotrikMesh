#include <algorithm>
#include "SemiGlobalSimplifier.h"

SemiGlobalSimplifier::SemiGlobalSimplifier() {}

SemiGlobalSimplifier::SemiGlobalSimplifier(Mesh& mesh_) : mesh(mesh_) {
    mu.SetMesh(mesh);
    smoother.SetMesh(mesh);
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
    mu.SetMesh(mesh);
    smoother.SetMesh(mesh);
}

void SemiGlobalSimplifier::SetSimplificationOperations() {
    CheckValidity();
    
    // SetDiagonalCollapseOperations();
    SetDirectSeparatrixOperations();
}

void SemiGlobalSimplifier::SetDiagonalCollapseOperations() {
    CheckValidity();

    for (auto& f: mesh.F) {
        std::unique_ptr<SimplificationOperation> dc1 = std::make_unique<DiagonalCollapse>(mesh, mu, f.id, 0, 2);
        dc1->SetRanking();
        Ops.push_back(std::move(dc1));

        std::unique_ptr<SimplificationOperation> dc2 = std::make_unique<DiagonalCollapse>(mesh, mu, f.id, 1, 3);
        dc2->SetRanking();
        Ops.push_back(std::move(dc2));
    }

    int i = 0;
    for (auto& op: Ops) {
        // op->PerformOperation();
        std::cout << op->ranking << std::endl;
        // if (i == 10) {
        //     break;
        // }
        // i += 1;
    }
    std::cout << Ops.size() << std::endl;

}

void SemiGlobalSimplifier::SetDirectSeparatrixOperations() {
    CheckValidity();

    Smoother meshSmoother(mesh);
    for (auto& v: mesh.V) {
        if (v.N_Vids.size() != 4 || v.type == FEATURE) continue;
        std::vector<size_t> c1, c2;
        for (auto vid: v.N_Vids) mesh.V.at(vid).N_Vids.size() == 3 ? c1.push_back(vid) : c2.push_back(vid);
        if (c1.size() < 2) continue;
        auto& s3_v1 = mesh.V.at(c1.at(0));
        auto& s3_v2 = mesh.V.at(c1.at(1));

        if (mu.GetDifference(s3_v1.N_Fids, s3_v2.N_Fids).size() != s3_v1.N_Fids.size()) continue;
        std::unique_ptr<SimplificationOperation> ds = std::make_unique<DirectSeparatrixCollapse>(mesh, mu, v.id, c1, c2);
        ds->SetRanking();
        Ops.push_back(std::move(ds));
    }

    std::cout << Ops.size() << std::endl;
    while (!Ops.empty()) {
        auto& op = Ops.back();
        std::cout << op->ranking << std::endl;
        op->PerformOperation();
        // meshSmoother.Smooth(op->smoothV);
        Ops.pop_back();
    }

    std::cout << "Smoothing mesh" << std::endl;
    std::vector<size_t> smoothv;
    for (auto& v: mesh.V) smoothv.push_back(v.id);
    meshSmoother.Smooth(smoothv);

}