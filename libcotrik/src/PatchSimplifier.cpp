#include "PatchSimplifier.h"
#include "EdgeRotateSimplifier.h"
#include "DoubletSimplifier.h"
#include "SheetSimplifier.h"
#include "TriangleSimplifier.h"
#include "GlobalSheetSimplifier.h"
#include "SingleSheetSimplifier.h"
#include "SheetSplitSimplifier.h"
#include "DiagnalCollapseSimplifier.h"


PatchSimplifier::PatchSimplifier(Mesh& mesh) : Simplifier(mesh) {

}

PatchSimplifier::~PatchSimplifier() {

}

void PatchSimplifier::Run() {
    auto maxValence_copy = Simplifier::maxValence;
    if (maxValence_copy > 5) Simplifier::maxValence = 5;
    int iter = 0;
    while (Simplifier::maxValence <= maxValence_copy) {
        while (iters-- > 0) {
            if (!Simplify(iter)) break;
        }
        ++Simplifier::maxValence;
    }
}

bool PatchSimplifier::Simplify(int& iter) {
    std::set<size_t> canceledFids;
    init();
    if (iter == 0 && featurePreserved) get_feature();
    if (iter == 0)
    {
        auto eids = get_rotate_eids();
        MeshFileWriter writer(mesh, "rotate_eids.vtk");
        writer.WriteEdgesVtk(eids);
    }
    if (checkCorner && !CheckCorners()) {
        std::cout << "----------------------- Failed in CheckCorners! -----------------------\n" << std::endl;
        update(canceledFids);
        return false;
    }

    // Step 1 -- doublet removal
    if (canceledFids.empty() && Simplifier::REMOVE_DOUBLET) {
        DoubletSimplifier doubletSimplifier(mesh);
        doubletSimplifier.Run(canceledFids);
        if (!canceledFids.empty()) std::cout << "remove_doublet" << std::endl;
    }
    // Step 2 -- doublet splitting
    if (canceledFids.empty() && Simplifier::SHEET_SPLIT) {
        SheetSplitSimplifier sheetSplitSimplifier(mesh);
        sheetSplitSimplifier.Run(canceledFids);
        if (!canceledFids.empty()) std::cout << "remove_doublet from sheetSplitSimplifier" << std::endl;
    }
    // Step -- triplet splitting (optional)
    if (canceledFids.empty() && Simplifier::TRIP) {
        TriangleSimplifier triangleSimplifier(mesh);
        triangleSimplifier.Run(canceledFids);
        if (!canceledFids.empty()) std::cout << "collapse faces from TriangleSimplifier" << std::endl;
    }
    // Step 3 -- edge rotation
    if (canceledFids.empty() && Simplifier::ROTATE) {
        EdgeRotateSimplifier edgeRotateSimplifier(mesh);
        edgeRotateSimplifier.Run(canceledFids);
        if (!canceledFids.empty()) std::cout << "rotate_edge" << std::endl;
    }
    static bool aligned = false;
    if (canceledFids.empty() && !aligned) {
        aligned = true;
        std::cout << "writing rotate.vtk " << std::endl;
        MeshFileWriter writer(mesh, "rotate.vtk");
        writer.WriteFile();
    }
    // Step 4 -- singlet collapsing
    if (canceledFids.empty()) {
        DiagnalCollapseSimplifier diagnalCollapseSimplifier(mesh);
        diagnalCollapseSimplifier.Run3(canceledFids);
        if (!canceledFids.empty()) std::cout << "singlet collapsing" << std::endl;
    }
    // Step 5 -- <separatrix splitting> and <separatrix splitting (optional)>
    if (canceledFids.empty() && (Simplifier::COLLAPSE || Simplifier::SPLIT)) {
        BaseComplexQuad baseComplex(mesh);
        baseComplex.ExtractSingularVandE();
        baseComplex.BuildE();
        strict_simplify(baseComplex, canceledFids);
        if (canceledFids.empty()) {
            loose_simplify(baseComplex, canceledFids);
            // loose_simplify_random(baseComplex, canceledFids);
            if (!canceledFids.empty()) std::cout << "loose_simplify\n";
//            else if (canceledFids.empty() && Simplifier::HALF) {
//                half_simplify(baseComplex, canceledFids);
//                if (!canceledFids.empty()) std::cout << "half_simplify\n";
//            }
        } else std::cout << "strict_simplify\n";
    }
    // Step 6 -- chord collapsing
    if (canceledFids.empty() && Simplifier::GLOBAL) {
        update(canceledFids);
        init();
        SingleSheetSimplifier sheetSimplifier(mesh);
        sheetSimplifier.Run(canceledFids);
        if (!canceledFids.empty()) std::cout << "chord collapsing" << std::endl;
    }
    // Step 7 -- half separatrix collapsing
    if (canceledFids.empty() && Simplifier::HALF) {
        update(canceledFids);
        init();
        BaseComplexQuad baseComplex(mesh);
        baseComplex.ExtractSingularVandE();
        baseComplex.BuildE();
        half_simplify(baseComplex, canceledFids);
        if (!canceledFids.empty()) std::cout << "half_simplify\n";
    }
    // Step 8 -- diagonal collapsing
    if (canceledFids.empty() && Simplifier::COLLAPSE_DIAGNAL) {
        DiagnalCollapseSimplifier diagnalCollapseSimplifier(mesh);
        diagnalCollapseSimplifier.Run(canceledFids);
        if (!canceledFids.empty()) std::cout << "collapse_diagnal" << std::endl;
    }

//    if (canceledFids.empty() && Simplifier::TRIP) {
//        TriangleSimplifier triangleSimplifier(mesh);
//        triangleSimplifier.Run(canceledFids);
//        if (!canceledFids.empty()) std::cout << "collapse faces from TriangleSimplifier" << std::endl;
//    }
//    if (canceledFids.empty() && Simplifier::GLOBAL) {
//        update(canceledFids);
//        init();
//        SheetSimplifier sheetSimplifier(mesh);
//        sheetSimplifier.Run(canceledFids);
//    }
//    if (canceledFids.empty() && Simplifier::HALF) {
//        update(canceledFids);
//        init();
//        BaseComplexQuad baseComplex(mesh);
//        baseComplex.ExtractSingularVandE();
//        baseComplex.BuildE();
//        half_simplify(baseComplex, canceledFids);
//        if (!canceledFids.empty()) std::cout << "half_simplify\n";
//    }
    if (canceledFids.empty()) {
        update(canceledFids);
        return false;
    }
    update(canceledFids);
    if (Simplifier::writeFile)
    {
        auto num = std::to_string(iter);
        while (num.size() < 3) num.insert(num.begin(), '0');
        {
            std::string fname = std::string("iter") + num + ".vtk";
            std::cout << "writing " << fname << std::endl;
            MeshFileWriter writer(mesh, fname.c_str());
            writer.WriteFile();
        }
//        {
//            std::string fname = std::string("VertexFeature") + num + ".vtk";
//            std::cout << "writing " << fname << std::endl;
//            MeshFileWriter writer(mesh, fname.c_str());
//            writer.WriteVertexFeatureVtk();
//        }
    }
    std::cout << "iter = " << iter++ << std::endl;
    std::cout << "---------------------------------------------------\n";
    return true;
}

bool PatchSimplifier::CheckCorners() {
    for (auto& v : mesh.V)
        if (v.isCorner && v.N_Fids.size() < v.idealValence)
            return false;
    return true;
}
