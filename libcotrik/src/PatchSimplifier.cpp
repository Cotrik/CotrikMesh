#include "PatchSimplifier.h"
#include "EdgeRotateSimplifier.h"
#include "DoubletSimplifier.h"
#include "SheetSimplifier.h"
#include "TriangleSimplifier.h"
#include "GlobalSheetSimplifier.h"
#include "SingleSheetSimplifier.h"
#include "SheetSplitSimplifier.h"
#include "DiagnalCollapseSimplifier.h"
#include "AngleBasedSmoothQuadMesh.h"


PatchSimplifier::PatchSimplifier(Mesh& mesh) : Simplifier(mesh) {
    // smoothing_algorithm = new SmoothAlgorithm(mesh, mesh, 200, 1, true, true);
}

PatchSimplifier::~PatchSimplifier() {

}

void PatchSimplifier::Run() {
    originalFaces = mesh.F.size();
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

/*void PatchSimplifier::Run() {
    init();
    if (featurePreserved) {
        get_feature();
        auto eids = get_rotate_eids();
        MeshFileWriter writer(mesh, "rotate_eids.vtk");
        writer.WriteEdgesVtk(eids);
    }
    auto maxValence_copy = Simplifier::maxValence;
    while (iters >= 0) {
        int iter = 0;
        if (maxValence_copy > 5) Simplifier::maxValence = 5;
        while (Simplifier::maxValence <= maxValence_copy) {
            std::cout << "Edge Rotation Loop" << std::endl;
            while (iters-- > 0) {
                if (!EdgeRotateSimplify(iter)) break;
                // if (!Simplify(iter)) break;
            }
            ++Simplifier::maxValence;
        }
        // Simplifier::maxValence = maxValence_copy;
        // if (maxValence_copy > 5) Simplifier::maxValence = 5;
        // iter = 0;
        // while (Simplifier::maxValence <= maxValence_copy) {
        //     std::cout << "Chord Collapse Loop" << std::endl;
        //     while (iters-- > 0) {
        //         if (!ChordCollapseSimplify(iter)) break;
        //         // if (!Simplify(iter)) break;
        //     }
        //     ++Simplifier::maxValence;
        //     // std::cout << "maxValence copy: " << maxValence_copy << " maxValence: " << Simplifier::maxValence << " iter: " << iter << std::endl;
        // }
        // // std::cout << "end chord collapsing, iter = " << iter << std::endl;
        // Simplifier::maxValence = maxValence_copy;
        // if (maxValence_copy > 5) Simplifier::maxValence = 5;
        // iter = 0;
        // while (Simplifier::maxValence <= maxValence_copy) {
        //     std::cout << "Separatrix Collapse Loop" << std::endl;
        //     while (iters-- > 0) {
        //         if (!SeparatrixCollapseSimplify(iter)) break;
        //         // if (!Simplify(iter)) break;
        //     }
        //     ++Simplifier::maxValence;
        // }
        // Simplifier::maxValence = maxValence_copy;
        // if (maxValence_copy > 5) Simplifier::maxValence = 5;
        // iter = 0;
        // while (Simplifier::maxValence <= maxValence_copy) {
        //     std::cout << "Half Separatrix Collapse Loop" << std::endl;
        //     while (iters-- > 0) {
        //         if (!HalfSeparatrixCollapseSimplify(iter)) break;
        //         // if (!Simplify(iter)) break;
        //     }
        //     ++Simplifier::maxValence;
        // }
        if (iter == 0) {
            break;
        }
    }
}*/

bool PatchSimplifier::ChordCollapseSimplify(int& iter) {
    std::set<size_t> canceledFids;
    init();
    // if (iter == 0 && featurePreserved) get_feature();
    // if (iter == 0)
    // {
    //     auto eids = get_rotate_eids();
    //     MeshFileWriter writer(mesh, "rotate_eids.vtk");
    //     writer.WriteEdgesVtk(eids);
    // }

    if (checkCorner && !CheckCorners()) {
        std::cout << "----------------------- Failed in CheckCorners! -----------------------\n" << std::endl;
        update(canceledFids);
        return false;
    }

    if (canceledFids.empty() && Simplifier::REMOVE_DOUBLET) {
        DoubletSimplifier doubletSimplifier(mesh);
        doubletSimplifier.RunCollective(canceledFids);
        if (!canceledFids.empty()) std::cout << "remove_doublet" << std::endl;
    }

    if (canceledFids.empty() && Simplifier::ROTATE) {
        EdgeRotateSimplifier edgeRotateSimplifier(mesh);
        edgeRotateSimplifier.RunCollective(canceledFids);
        if (!canceledFids.empty()) std::cout << "rotate_edge" << std::endl;
    }    

    // Step 6 -- chord collapsing
    if (canceledFids.empty() && Simplifier::GLOBAL) {
        update(canceledFids);
        init();
        SingleSheetSimplifier sheetSimplifier(mesh);
        sheetSimplifier.ExtractAndCollapse(canceledFids);
        if (!canceledFids.empty()) std::cout << "chord collapsing" << std::endl;
    }

    // if (canceledFids.empty() && Simplifier::REMOVE_DOUBLET) {
    //     DoubletSimplifier doubletSimplifier(mesh);
    //     doubletSimplifier.RunCollective(canceledFids);
    //     if (!canceledFids.empty()) std::cout << "remove_doublet" << std::endl;
    // }

    // if (canceledFids.empty() && Simplifier::ROTATE) {
    //     EdgeRotateSimplifier edgeRotateSimplifier(mesh);
    //     edgeRotateSimplifier.RunCollective(canceledFids);
    //     if (!canceledFids.empty()) std::cout << "rotate_edge" << std::endl;
    // }


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
    }
    
    std::cout << "iter = " << iter++ << std::endl;
    std::cout << "---------------------------------------------------\n";
    return true;
}

bool PatchSimplifier::SeparatrixSplitSimplify(int& iter) {
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

    if (canceledFids.empty() && Simplifier::REMOVE_DOUBLET) {
        DoubletSimplifier doubletSimplifier(mesh);
        doubletSimplifier.RunCollective(canceledFids);
        if (!canceledFids.empty()) std::cout << "remove_doublet" << std::endl;
    }

    if (canceledFids.empty() && Simplifier::ROTATE) {
        EdgeRotateSimplifier edgeRotateSimplifier(mesh);
        edgeRotateSimplifier.RunCollective(canceledFids);
        if (!canceledFids.empty()) std::cout << "rotate_edge" << std::endl;
    }    

    if (canceledFids.empty() && Simplifier::SPLIT) {
        BaseComplexQuad baseComplex(mesh);
        baseComplex.ExtractSingularVandE();
        baseComplex.BuildE();

        five_connections_split(baseComplex, canceledFids, false);
        if (canceledFids.empty()) {
            five_connections_split(baseComplex, canceledFids, true);
            if (!canceledFids.empty()) std::cout << "loose_simplify\n";
        } else std::cout << "strict_simplify\n";
    }

    // if (canceledFids.empty() && Simplifier::REMOVE_DOUBLET) {
    //     DoubletSimplifier doubletSimplifier(mesh);
    //     doubletSimplifier.RunCollective(canceledFids);
    //     if (!canceledFids.empty()) std::cout << "remove_doublet" << std::endl;
    // }

    // if (canceledFids.empty() && Simplifier::ROTATE) {
    //     EdgeRotateSimplifier edgeRotateSimplifier(mesh);
    //     edgeRotateSimplifier.RunCollective(canceledFids);
    //     if (!canceledFids.empty()) std::cout << "rotate_edge" << std::endl;
    // } 

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
    }

    std::cout << "iter = " << iter++ << std::endl;
    std::cout << "---------------------------------------------------\n";
    return true;
}

bool PatchSimplifier::SeparatrixCollapseSimplify(int& iter) {
    std::set<size_t> canceledFids;
    init();
    // if (iter == 0 && featurePreserved) get_feature();
    // if (iter == 0)
    // {
    //     auto eids = get_rotate_eids();
    //     MeshFileWriter writer(mesh, "rotate_eids.vtk");
    //     writer.WriteEdgesVtk(eids);
    // }

    if (checkCorner && !CheckCorners()) {
        std::cout << "----------------------- Failed in CheckCorners! -----------------------\n" << std::endl;
        update(canceledFids);
        return false;
    }
    
    if (canceledFids.empty() && Simplifier::REMOVE_DOUBLET) {
        DoubletSimplifier doubletSimplifier(mesh);
        doubletSimplifier.RunCollective(canceledFids);
        if (!canceledFids.empty()) std::cout << "remove_doublet" << std::endl;
    }

    if (canceledFids.empty() && Simplifier::ROTATE) {
        EdgeRotateSimplifier edgeRotateSimplifier(mesh);
        edgeRotateSimplifier.RunCollective(canceledFids);
        if (!canceledFids.empty()) std::cout << "rotate_edge" << std::endl;
    }    

    if (canceledFids.empty() && Simplifier::COLLAPSE) {
        BaseComplexQuad baseComplex(mesh);
        baseComplex.ExtractSingularVandE();
        baseComplex.BuildE();

        three_connections_collapse(baseComplex, canceledFids, false);
        if (canceledFids.empty()) {
            three_connections_collapse(baseComplex, canceledFids, true);
            if (!canceledFids.empty()) std::cout << "loose_simplify\n";
        } else std::cout << "strict_simplify\n";
    }

    // if (canceledFids.empty() && Simplifier::REMOVE_DOUBLET) {
    //     DoubletSimplifier doubletSimplifier(mesh);
    //     doubletSimplifier.RunCollective(canceledFids);
    //     if (!canceledFids.empty()) std::cout << "remove_doublet" << std::endl;
    // }

    // if (canceledFids.empty() && Simplifier::ROTATE) {
    //     EdgeRotateSimplifier edgeRotateSimplifier(mesh);
    //     edgeRotateSimplifier.RunCollective(canceledFids);
    //     if (!canceledFids.empty()) std::cout << "rotate_edge" << std::endl;
    // }

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
    }

    std::cout << "iter = " << iter++ << std::endl;
    std::cout << "---------------------------------------------------\n";
    return true;
}

bool PatchSimplifier::HalfSeparatrixCollapseSimplify(int& iter) {
    std::set<size_t> canceledFids;
    init();
    // if (iter == 0 && featurePreserved) get_feature();
    // if (iter == 0)
    // {
    //     auto eids = get_rotate_eids();
    //     MeshFileWriter writer(mesh, "rotate_eids.vtk");
    //     writer.WriteEdgesVtk(eids);
    // }

    if (checkCorner && !CheckCorners()) {
        std::cout << "----------------------- Failed in CheckCorners! -----------------------\n" << std::endl;
        update(canceledFids);
        return false;
    }

    if (canceledFids.empty() && Simplifier::REMOVE_DOUBLET) {
        DoubletSimplifier doubletSimplifier(mesh);
        doubletSimplifier.RunCollective(canceledFids);
        if (!canceledFids.empty()) std::cout << "remove_doublet" << std::endl;
    }

    if (canceledFids.empty() && Simplifier::ROTATE) {
        EdgeRotateSimplifier edgeRotateSimplifier(mesh);
        edgeRotateSimplifier.RunCollective(canceledFids);
        if (!canceledFids.empty()) std::cout << "rotate_edge" << std::endl;
    }    

    if (canceledFids.empty() && Simplifier::HALF) {
        update(canceledFids);
        init();
        BaseComplexQuad baseComplex(mesh);
        baseComplex.ExtractSingularVandE();
        baseComplex.BuildE();
        half_separatrix_collapse(baseComplex, canceledFids);
        if (!canceledFids.empty()) std::cout << "half_simplify\n";
    }

    // if (canceledFids.empty() && Simplifier::REMOVE_DOUBLET) {
    //     DoubletSimplifier doubletSimplifier(mesh);
    //     doubletSimplifier.RunCollective(canceledFids);
    //     if (!canceledFids.empty()) std::cout << "remove_doublet" << std::endl;
    // }

    // if (canceledFids.empty() && Simplifier::ROTATE) {
    //     EdgeRotateSimplifier edgeRotateSimplifier(mesh);
    //     edgeRotateSimplifier.RunCollective(canceledFids);
    //     if (!canceledFids.empty()) std::cout << "rotate_edge" << std::endl;
    // }


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
    }

    std::cout << "iter = " << iter++ << std::endl;
    std::cout << "---------------------------------------------------\n";
    return true;
}

bool PatchSimplifier::EdgeRotateSimplify(int& iter) {
    std::set<size_t> canceledFids;
    init();
    // if (iter == 0 && featurePreserved) get_feature();
    // if (iter == 0)
    // {
    //     auto eids = get_rotate_eids();
    //     MeshFileWriter writer(mesh, "rotate_eids.vtk");
    //     writer.WriteEdgesVtk(eids);
    // }

    if (checkCorner && !CheckCorners()) {
        std::cout << "----------------------- Failed in CheckCorners! -----------------------\n" << std::endl;
        update(canceledFids);
        return false;
    }

    // if (canceledFids.empty() && Simplifier::REMOVE_DOUBLET) {
    //     DoubletSimplifier doubletSimplifier(mesh);
    //     doubletSimplifier.RunCollective(canceledFids);
    //     if (!canceledFids.empty()) std::cout << "remove_doublet" << std::endl;
    // }

    if (canceledFids.empty() && Simplifier::ROTATE) {
        EdgeRotateSimplifier edgeRotateSimplifier(mesh);
        edgeRotateSimplifier.RunCollective(canceledFids);
        if (!canceledFids.empty()) std::cout << "rotate_edge" << std::endl;
    }

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
    }

    std::cout << "iter = " << iter++ << std::endl;
    std::cout << "---------------------------------------------------\n";
    return true;
}

bool PatchSimplifier::Simplify(int& iter) {
    std::set<size_t> canceledFids;
    init();
    // if (findCrossQuads()) {
    //     std::cout << "After cross quads check" << std::endl;
    //     mesh.smoothGlobal = false;
    //     smooth_project(2);
    //     mesh.smoothGlobal = true;
    //     // SmoothAlgorithm algorithm(mesh, mesh, 1000, 1, true, true);
    //     // algorithm.smoothMesh();
    //     // update(canceledFids);
    //     // return false;
    // }
    if (iter == 0 && featurePreserved) get_feature();
    if (iter == 0)
    {
        auto eids = get_rotate_eids();
        MeshFileWriter writer(mesh, "rotate_eids.vtk");
        writer.WriteEdgesVtk(eids);
    }
    // if (checkCorner && !CheckCorners()) {
    //     std::cout << "----------------------- Failed in CheckCorners! -----------------------\n" << std::endl;
    //     update(canceledFids);
    //     return false;
    // }
    // if (mesh.F.size() <= originalFaces * 0.5) {
        // std::cout << "new Faces: " << mesh.F.size() << ", original Faces: " << originalFaces << std::endl;
        // smooth_project();
        // {
            // RefineMesh();        
        // }
    // }
    // smooth_project();
    // {
        // SmoothAlgorithm SmoothAlgo(mesh, mesh, 200, 1, true, true);
        // SmoothAlgo.smoothMesh();
    // }
    
    // Step 1 -- doublet removal
    if (canceledFids.empty() && Simplifier::REMOVE_DOUBLET) {
        // update(canceledFids);
        // init();
        DoubletSimplifier doubletSimplifier(mesh);
        doubletSimplifier.Run(canceledFids);
        // doubletSimplifier.RunCollective(canceledFids);
        if (!canceledFids.empty()) std::cout << "remove_doublet" << std::endl;
    }
    // Step 2 -- doublet splitting
    if (canceledFids.empty() && Simplifier::SHEET_SPLIT) {
        // update(canceledFids);
        // init();
        SheetSplitSimplifier sheetSplitSimplifier(mesh);
        sheetSplitSimplifier.Run(canceledFids);
        if (!canceledFids.empty()) std::cout << "remove_doublet from sheetSplitSimplifier" << std::endl;
    }
    // Step -- triplet splitting (optional)
    // if (canceledFids.empty() && Simplifier::TRIP) {
    //     update(canceledFids);
    //     init();
    //     TriangleSimplifier triangleSimplifier(mesh);
    //     triangleSimplifier.Run(canceledFids);
    //     if (!canceledFids.empty()) std::cout << "collapse faces from TriangleSimplifier" << std::endl;
    // }
    // Step 3 -- edge rotation
    if (canceledFids.empty() && Simplifier::ROTATE) {
        // update(canceledFids);
        // init();
        EdgeRotateSimplifier edgeRotateSimplifier(mesh);
        edgeRotateSimplifier.Run(canceledFids);
        // edgeRotateSimplifier.RunCollective(canceledFids);
        if (!canceledFids.empty()) std::cout << "rotate_edge" << std::endl;
    }
    
    // static bool aligned = false;
    // if (canceledFids.empty() && !aligned) {
    //     aligned = true;
    //     std::cout << "writing rotate.vtk " << std::endl;
    //     MeshFileWriter writer(mesh, "rotate.vtk");
    //     writer.WriteFile();
    // }
    // Step 4 -- singlet collapsing
    if (canceledFids.empty()) {
        // update(canceledFids);
        // init();
        DiagnalCollapseSimplifier diagnalCollapseSimplifier(mesh);
        diagnalCollapseSimplifier.Run3(canceledFids);
        if (!canceledFids.empty()) std::cout << "singlet collapsing" << std::endl;
    }
    // Step 5 -- <separatrix splitting> and <separatrix splitting (optional)>
    if (canceledFids.empty() && (Simplifier::COLLAPSE || Simplifier::SPLIT)) {
        // update(canceledFids);
        // init();
        BaseComplexQuad baseComplex(mesh);
        baseComplex.ExtractSingularVandE();
        baseComplex.BuildE();
        // strict_simplify(baseComplex, canceledFids);
        // three_connections_strict_collapse(baseComplex, canceledFids);
        three_connections_collapse(baseComplex, canceledFids, false);
        if (canceledFids.empty()) {
            // loose_simplify(baseComplex, canceledFids);
            // three_connections_loose_collapse(baseComplex, canceledFids);
            three_connections_collapse(baseComplex, canceledFids, true);
            // loose_simplify_random(baseComplex, canceledFids);
            if (!canceledFids.empty()) std::cout << "loose_simplify\n";
//            else if (canceledFids.empty() && Simplifier::HALF) {
//                half_simplify(baseComplex, canceledFids);
//                if (!canceledFids.empty()) std::cout << "half_simplify\n";
//            }
        } else std::cout << "strict_simplify\n";
    }
    // Step 6 -- chord collapsing
    // if (canceledFids.empty() && Simplifier::GLOBAL) {
    //     update(canceledFids);
    //     init();
    //     SingleSheetSimplifier sheetSimplifier(mesh);
    //     sheetSimplifier.Run(canceledFids);
    //     if (!canceledFids.empty()) std::cout << "chord collapsing" << std::endl;
    // }
    // Step 7 -- half separatrix collapsing
    // if (canceledFids.empty() && Simplifier::HALF) {
    //     update(canceledFids);
    //     init();
    //     BaseComplexQuad baseComplex(mesh);
    //     baseComplex.ExtractSingularVandE();
    //     baseComplex.BuildE();
    //     half_simplify(baseComplex, canceledFids);
    //     if (!canceledFids.empty()) std::cout << "half_simplify\n";
    // }
    // Step 8 -- diagonal collapsing
    if (canceledFids.empty() && Simplifier::COLLAPSE_DIAGNAL) {
        // update(canceledFids);
        // init();
        DiagnalCollapseSimplifier diagnalCollapseSimplifier(mesh);
        diagnalCollapseSimplifier.Run(canceledFids);
        if (!canceledFids.empty()) std::cout << "collapse_diagnal" << std::endl;
    }
    

   if (canceledFids.empty() && Simplifier::TRIP) {
    //    update(canceledFids);
    //    init();
       TriangleSimplifier triangleSimplifier(mesh);
       triangleSimplifier.Run(canceledFids);
       if (!canceledFids.empty()) std::cout << "collapse faces from TriangleSimplifier" << std::endl;
   }
   if (canceledFids.empty() && Simplifier::GLOBAL) {
    //    update(canceledFids);
    //    init();
       SheetSimplifier sheetSimplifier(mesh);
    //    sheetSimplifier.Run(canceledFids);
       sheetSimplifier.ExtractAndCollapse(canceledFids);
   }
   if (canceledFids.empty() && Simplifier::HALF) {
    //    update(canceledFids);
    //    init();
       BaseComplexQuad baseComplex(mesh);
       baseComplex.ExtractSingularVandE();
       baseComplex.BuildE();
    //    half_simplify(baseComplex, canceledFids);
       half_separatrix_collapse(baseComplex, canceledFids);

       if (!canceledFids.empty()) std::cout << "half_simplify\n";
   }
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

bool PatchSimplifier::SimplifyCollective(int& iter) {
    std::set<size_t> canceledFids;
    init();
    if (iter == 0 && featurePreserved) get_feature();
    if (iter == 0)
    {
        auto eids = get_rotate_eids();
        MeshFileWriter writer(mesh, "rotate_eids.vtk");
        writer.WriteEdgesVtk(eids);
    }

    // if (iter % 50 == 0)
    // {
    //     std::string fname = std::string("simplified_") +  std::to_string(iter) + ".vtk";
    //     MeshFileWriter writer(mesh, fname.c_str());
    //     writer.WriteFile();
    // }
    
    // smoothing_algorithm->resampleBoundaryVertices();
    // smoothing_algorithm->smoothLaplacianScaleBased();
    // if (iter % 50 == 0)
    // {
    //     std::string fname = std::string("smoothed_") +  std::to_string(iter) + ".vtk";
    //     MeshFileWriter writer(mesh, fname.c_str());
    //     writer.WriteFile();
    // }
    if (checkCorner && !CheckCorners()) {
        std::cout << "----------------------- Failed in CheckCorners! -----------------------\n" << std::endl;
        update(canceledFids);
        return false;
    }

    // Step 1 -- doublet removal
    if (canceledFids.empty() && Simplifier::REMOVE_DOUBLET) {
        DoubletSimplifier doubletSimplifier(mesh);
        doubletSimplifier.RunCollective(canceledFids);
        // doubletSimplifier.Run(canceledFids);
        if (!canceledFids.empty()) std::cout << "remove_doublet" << std::endl;
    }
    // Step 2 -- doublet splitting
    /*if (canceledFids.empty() && Simplifier::SHEET_SPLIT) {
        SheetSplitSimplifier sheetSplitSimplifier(mesh);
        sheetSplitSimplifier.Run(canceledFids);
        if (!canceledFids.empty()) std::cout << "remove_doublet from sheetSplitSimplifier" << std::endl;
    }*/
    // Step -- triplet splitting (optional)
    /*if (canceledFids.empty() && Simplifier::TRIP) {
        TriangleSimplifier triangleSimplifier(mesh);
        triangleSimplifier.Run(canceledFids);
        if (!canceledFids.empty()) std::cout << "collapse faces from TriangleSimplifier" << std::endl;
    }*/
    // Step 3 -- edge rotation
    if (canceledFids.empty() && Simplifier::ROTATE) {
        EdgeRotateSimplifier edgeRotateSimplifier(mesh);
        edgeRotateSimplifier.RunCollective(canceledFids);
        // edgeRotateSimplifier.Run(canceledFids);
        if (!canceledFids.empty()) std::cout << "rotate_edge" << std::endl;
    }
    /*static bool aligned = false;
    if (canceledFids.empty() && !aligned) {
        aligned = true;
        std::cout << "writing rotate.vtk " << std::endl;
        MeshFileWriter writer(mesh, "rotate.vtk");
        writer.WriteFile();
    }*/
    // Step 4 -- singlet collapsing
    /*if (canceledFids.empty()) {
        DiagnalCollapseSimplifier diagnalCollapseSimplifier(mesh);
        diagnalCollapseSimplifier.Run3(canceledFids);
        if (!canceledFids.empty()) std::cout << "singlet collapsing" << std::endl;
    }*/
    // Step 5 -- <separatrix splitting> and <separatrix splitting (optional)>
    if (canceledFids.empty() && (Simplifier::COLLAPSE)) {
        BaseComplexQuad baseComplex(mesh);
        baseComplex.ExtractSingularVandE();
        baseComplex.BuildE();
        // strict_simplify(baseComplex, canceledFids);
        three_connections_collapse(baseComplex, canceledFids, false);
        // three_connections_strict_collapse(baseComplex, canceledFids);
        if (canceledFids.empty()) {
            // loose_simplify(baseComplex, canceledFids);
            // three_connections_loose_collapse(baseComplex, canceledFids);
            three_connections_collapse(baseComplex, canceledFids, true);
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
        // sheetSimplifier.Run(canceledFids);
        sheetSimplifier.ExtractAndCollapse(canceledFids);
        if (!canceledFids.empty()) std::cout << "chord collapsing" << std::endl;
    }
    
    // std::cout << iter << " " << canceledFids.size() << std::endl;

    // Step 7 -- half separatrix collapsing
    if (canceledFids.empty() && Simplifier::HALF) {
        update(canceledFids);
        init();
        BaseComplexQuad baseComplex(mesh);
        baseComplex.ExtractSingularVandE();
        baseComplex.BuildE();
        // half_simplify(baseComplex, canceledFids);
        half_separatrix_collapse(baseComplex, canceledFids);
        if (!canceledFids.empty()) std::cout << "half_simplify\n";
    }
    // Step 8 -- diagonal collapsing
    /*if (canceledFids.empty() && Simplifier::COLLAPSE_DIAGNAL) {
        DiagnalCollapseSimplifier diagnalCollapseSimplifier(mesh);
        diagnalCollapseSimplifier.Run(canceledFids);
        if (!canceledFids.empty()) std::cout << "collapse_diagnal" << std::endl;
    }*/

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
    // init();
    // std::cout << mesh.V.size() << " " << mesh.E.size() << " " << mesh.F.size() << " " << mesh.C.size() << std::endl;
    // if (smoothing_algorithm->isMeshNonManifold()) {
    //     std::cout << "iter: " << iter << std::endl;
    //     return false;
    // }
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

    // init();
    // smoothing_algorithm->smoothMesh();
    
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

void PatchSimplifier::smoothMesh(int iters_, bool global) {
    // mesh.SetOneRingNeighborhood();
    int it = 0;
    while (it < iters_) {
        // std::cout << "it: " << it << std::endl;
        std::vector<glm::dvec3> new_coords(mesh.V.size(), glm::dvec3(0.0, 0.0, 0.0));
        int v_pos = 0;
        for (auto& v: mesh.V) {
            if (v.isBoundary) {
                v_pos++;
                continue;
            }
            if (!global && !v.smoothLocal) {
                v_pos++;
                continue;
            }
            for (auto vid: v.N_Vids) {
                Vertex& v_prime = mesh.V.at(vid);
                glm::dvec3 direction(v_prime.x - v.x, v_prime.y - v.y, v_prime.z - v.z);
                // new_coords.at(v_pos) += glm::dvec3(v_prime.x, v_prime.y, v_prime.z);
                new_coords.at(v_pos) += direction;
            }
            new_coords.at(v_pos).x = (new_coords.at(v_pos).x / (v.N_Vids.size()));
            new_coords.at(v_pos).y = (new_coords.at(v_pos).y / (v.N_Vids.size()));
            v_pos++;
        }
        v_pos = 0;
        for (auto& v: mesh.V) {
            if (v.isBoundary) {
                v_pos++;
                continue;
            }
            if (!global && !v.smoothLocal) {
                v_pos++;
                continue;
            }
            v.x += (new_coords.at(v_pos).x);
            v.y += (new_coords.at(v_pos).y);
            v_pos++;
        }
        it++;
    }
}

void PatchSimplifier::RefineMesh() {
    std::set<size_t> canceledFids;
    for (auto& e: mesh.E) {
        Vertex& v1 = mesh.V.at(e.Vids.at(0));
        Vertex& v2 = mesh.V.at(e.Vids.at(1));
        glm::dvec3 a(v1.x, v1.y, v1.z);
        glm::dvec3 b(v2.x - v1.x, v2.y - v1.y, v2.z - v1.z);
        std::cout << "v1: " << v1.x << " " << v1.y << " " << v1.z << std::endl;
        std::cout << "v2: " << v2.x << " " << v2.y << " " << v2.z << std::endl;
        std::cout << "a: " << a.x << " " << a.y << " " << a.z << std::endl;
        std::cout << "b: " << b.x << " " << b.y << " " << b.z << std::endl;

        Vertex newV((a + (glm::normalize(b) * glm::length(b) * 0.5)));
        if (glm::length(b) == 0) {
            newV.x = a.x;
            newV.y = a.y;
            newV.z = a.z;
        }
        std::cout << "new Vertex: " << newV.x << " " << newV.y << " " << newV.z << std::endl;
        std::cout << "---------------------------------------------------" << std::endl;

        newV.id = mesh.V.size();
        mesh.V.push_back(newV);
        for (auto fid: e.N_Fids) {
            Face& f = mesh.F.at(fid);
            for (int i = 0; i < f.Vids.size(); i++) {
                if (std::find(e.Vids.begin(), e.Vids.end(), f.Vids.at(i)) != e.Vids.end() &&
                    std::find(e.Vids.begin(), e.Vids.end(), f.Vids.at((i + 1) % f.Vids.size())) != e.Vids.end()) {
                        f.Vids.insert(f.Vids.begin() + (i + 1) % f.Vids.size(), newV.id);
                        break;
                }
            }
        }
    }
    int n = mesh.F.size();
    for (int i = 0; i < n; i++) {
        Face& f = mesh.F.at(i);

        Vertex& v1 = mesh.V.at(f.Vids.at(0));
        Vertex& v2 = mesh.V.at(f.Vids.at(2));
        Vertex& v3 = mesh.V.at(f.Vids.at(4));
        Vertex& v4 = mesh.V.at(f.Vids.at(6));

        glm::dvec3 a1(v1.x, v1.y, v1.z);
        glm::dvec3 b1(v3.x - v1.x, v3.y - v1.y, v3.z - v1.z);
        glm::dvec3 c1((a1 + (glm::normalize(b1) * glm::length(b1) * 0.5)));

        
        glm::dvec3 a2(v2.x, v2.y, v2.z);
        glm::dvec3 b2(v4.x - v2.x, v4.y - v2.y, v4.z - v2.z);
        glm::dvec3 c2((a2 + (glm::normalize(b2) * glm::length(b2) * 0.5)));

        glm::dvec3 a3(c1.x, c1.y, c1.z);
        glm::dvec3 b3(c2.x - c1.x, c2.y - c1.y, c2.z - c1.z);

        Vertex newV((a3 + (glm::normalize(b3) * glm::length(b3) * 0.5)));
        if (glm::length(b3) == 0) {
            newV.x = a3.x;
            newV.y = a3.y;
            newV.z = a3.z;
        }

        newV.id = mesh.V.size();
        mesh.V.push_back(newV);

        Face f1;
        f1.Vids = {newV.id, f.Vids.at(0), f.Vids.at(1), f.Vids.at(2)};
        
        Face f2;
        f2.Vids = {newV.id, f.Vids.at(2), f.Vids.at(3), f.Vids.at(4)};
        
        Face f3;
        f3.Vids = {newV.id, f.Vids.at(4), f.Vids.at(5), f.Vids.at(6)};
        
        Face f4;
        f4.Vids = {newV.id, f.Vids.at(6), f.Vids.at(7), f.Vids.at(0)};

        f1.id = mesh.F.size();
        mesh.F.push_back(f1);

        f2.id = mesh.F.size();
        mesh.F.push_back(f2);

        f3.id = mesh.F.size();
        mesh.F.push_back(f3);

        f4.id = mesh.F.size();
        mesh.F.push_back(f4);
        canceledFids.insert(i);
    }
    update(canceledFids);
    init();
    // for (auto& v: mesh.V) {
    //     v.type = NULL;
	// 	v.isCorner = NULL;
	// 	v.isConvex = NULL;
	// 	v.label = NULL;
	// 	v.patch_id = NULL;
	// 	v.isSpecial = NULL;
	// 	v.labels.clear();
	// 	v.patch_ids.clear();
	// 	v.idealValence = NULL;
	// 	v.prescribed_length = NULL;
	// 	v.smoothLocal = NULL;
    // }
    // get_feature();
}

int PatchSimplifier::GetBoundaryEdges() {
    int numberBoundaryEdges = 0;
    for (auto& e: mesh.E) {
        if (e.isBoundary) {
            numberBoundaryEdges += 1;
        }
    }
    return numberBoundaryEdges;
}

bool PatchSimplifier::findCrossQuads() {
    std::vector<size_t> faceIds;
    bool foundCrossQuad = false;
    int nCrossQuads = 0;
    for (auto& v: mesh.V) {
        v.smoothLocal = false;
    }
    for (auto& f: mesh.F) {
        int sign = 0;
        for (int i = 0; i < f.Vids.size(); i++) {
            Vertex& a = mesh.V.at(f.Vids.at(i));
            Vertex& b = mesh.V.at(f.Vids.at((i + 1) % f.Vids.size()));
            Vertex& c = mesh.V.at(f.Vids.at((i + 2) % f.Vids.size()));
            double det = ((b.x - a.x) * (c.y - a.y)) - ((c.x - a.x) * (b.y - a.y));
            if (det > 0) {
                sign += 1;
            } else if (det <= 0) {
                sign -= 1;
            }
        }
        if (sign == 0) {
            nCrossQuads += 1;
            foundCrossQuad = true;
            faceIds.push_back(f.id);
            for (auto vid: f.Vids) {
                auto& v = mesh.V.at(vid);
                v.smoothLocal = true;
                for (auto nvid: v.N_Vids) {
                    mesh.V.at(nvid).smoothLocal = true;
                }
            }
        }
    }
    // if (foundCrossQuad) {
    //     std::ofstream ofs("cross_quads.vtk");
    //     ofs << "# vtk DataFile Version 3.0\n"
    //         << "cross_quads.vtk\n"
    //         << "ASCII\n\n"
    //         << "DATASET UNSTRUCTURED_GRID\n";
    //     ofs << "POINTS " << mesh.V.size() << " double\n";

    //     for (auto& v: mesh.V) {
    //         ofs << v.x << " " << v.y << " " << v.z << std::endl;
    //     }

    //     ofs << "CELLS " << faceIds.size() << " " << faceIds.size() * 5 << std::endl;
    //     for (auto fid: faceIds) {
    //         Face& f = mesh.F.at(fid);
    //         ofs << "4 " << f.Vids.at(0) << " " << f.Vids.at(1) << " " << f.Vids.at(2) << " " << f.Vids.at(3) <<  std::endl; 
    //         std::cout << "writing face: " << f.id << ", has " << f.Vids.size() << " vertices" << std::endl;
    //     }
    //     ofs << "CELL_TYPES " << faceIds.size() << "\n";
    //     for (auto fid: faceIds) {
    //         ofs << "9 " << std::endl; 
    //     }
    // }
    std::cout << "# cross quads: " << nCrossQuads << std::endl;
    return foundCrossQuad;
}