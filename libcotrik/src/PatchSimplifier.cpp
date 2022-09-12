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
#include "MeshQuality.h"


PatchSimplifier::PatchSimplifier(Mesh& mesh) : Simplifier(mesh) {
    // smoothing_algorithm = new SmoothAlgorithm(mesh, mesh, 200, 1, true, true);
}

PatchSimplifier::~PatchSimplifier() {

}

bool OpSortAscending(SimplificationOperationStruct op1, SimplificationOperationStruct op2) {
    return op1.profitability < op2.profitability;
}

bool OpSortDescending(SimplificationOperationStruct op1, SimplificationOperationStruct op2) {
    return op1.profitability > op2.profitability;
}

static void refineVertexInFaces(Mesh& mesh, std::vector<Vertex>& refinedV, int resolution = 3) {
    Vertex vertex;
    for (auto& f : mesh.F) {
        for (double u = 0; u < resolution; ++u)
            for (double v = 0; v < resolution; ++v) {
                auto base = 1.0 / (resolution + 1);
                auto v01 = (1.0 + u) * base * mesh.V.at(f.Vids[0]).xyz() + (resolution - u) * base * mesh.V.at(f.Vids[1]).xyz();
                auto v32 = (1.0 + u) * base * mesh.V.at(f.Vids[3]).xyz() + (resolution - u) * base * mesh.V.at(f.Vids[2]).xyz();
                
				vertex = (1.0 + v) * base * v01 + (resolution - v) * base * v32;
                // vertex.id = mesh.V.size();
                vertex.id = refinedV.size();
                vertex.patch_id = f.label;
                // mesh.V.push_back(vertex);
                refinedV.push_back(vertex);
            }
    }
}

static void refineVertexInEdges(Mesh& mesh, std::vector<Vertex>& refinedV, int resolution = 3) {
    auto& origMesh = mesh;
    Vertex vertex;
    for (auto& e : origMesh.E) {
        for (double u = 0; u < resolution; ++u) {
            auto base = 1.0 / (resolution + 1);
            vertex = (1.0 + u) * base * origMesh.V.at(e.Vids[0]).xyz() + (resolution - u) * base * origMesh.V.at(e.Vids[1]).xyz();
            // vertex.id = origMesh.V.size();
            vertex.id = refinedV.size();

            if (!e.isSharpFeature) {
                vertex.patch_id = origMesh.V.at(e.Vids[0]).type == REGULAR ? origMesh.V.at(e.Vids[0]).patch_id : origMesh.V.at(e.Vids[1]).patch_id;
                vertex.label = MAXID;
            } else {
                if (origMesh.V.at(e.Vids[0]).label == origMesh.V.at(e.Vids[1]).label)
                    vertex.label = origMesh.V.at(e.Vids[0]).label;
                else vertex.label = origMesh.V.at(e.Vids[0]).isCorner ? origMesh.V.at(e.Vids[1]).label : origMesh.V.at(e.Vids[0]).label;
                vertex.type = FEATURE;
            }
			vertex.isBoundary = e.isBoundary;
            // origMesh.V.push_back(vertex);
            refinedV.push_back(vertex);
        }
    }
}

void PatchSimplifier::Run() {
    originalFaces = mesh.F.size();
    auto maxValence_copy = Simplifier::maxValence;
    if (maxValence_copy > 5) Simplifier::maxValence = 5;
    int iter = 0;
    while (Simplifier::maxValence <= maxValence_copy) {
        while (iters-- > 0) {
            if (!Simplify(iter)) break;
            // if (!SimplifyMesh(iter)) break;
        }
        ++Simplifier::maxValence;
    }
    
    // init();
    
    // RefineMesh();
    // smoothGlobal = true;
    // SmoothMesh(true);
}

bool PatchSimplifier::SimplifyMesh(int& iter) {
    std::set<size_t> canceledFids;
    init();

    if (iter == 0) {
        featurePreserved ? get_feature() : origMesh = mesh;
        for (auto& v: origMesh.V) {
            if (v.isBoundary && !v.isCorner) {
                origBoundaryVids.push_back(v.id);
            }
        }
        // smoothGlobal = true;
        mesh.totalArea = mesh.GetQuadMeshArea();
        SmoothMesh(true);
    }
    else {
        // smoothGlobal = true;
        SmoothMesh(false);
    }
    mesh.prescribed_length = sqrt(mesh.totalArea / mesh.F.size());    
    // if (mesh.F.size() <= 0.25 * refinementFactor * origMesh.F.size()) {
    //     std::cout << "current Faces: " << mesh.F.size() << " original Faces: " << origMesh.F.size() << std::endl;
        
    //     RefineMesh();
    //     smoothGlobal = true;
    //     SmoothMesh(true);
    // }

    bool(*fn_pt1)(SimplificationOperationStruct, SimplificationOperationStruct) = OpSortDescending;
    bool(*fn_pt2)(SimplificationOperationStruct, SimplificationOperationStruct) = OpSortAscending;
    std::multiset<SimplificationOperationStruct, bool(*)(SimplificationOperationStruct, SimplificationOperationStruct)> SimplificationOps(fn_pt2);


    ///////////////// Cleaning Operations (Local Operations) ////////////////////
    // -- singlet collapsing
    if (canceledFids.empty()) {
        DiagnalCollapseSimplifier diagnalCollapseSimplifier(mesh);
        diagnalCollapseSimplifier.CollapseSinglets(canceledFids);
        if (!canceledFids.empty()) std::cout << "singlet collapsing" << std::endl;
    }
    // -- doublet removal
    if (canceledFids.empty() && Simplifier::REMOVE_DOUBLET) {
        DoubletSimplifier doubletSimplifier(mesh);
        doubletSimplifier.RunCollective(canceledFids);
        if (!canceledFids.empty()) std::cout << "remove_doublet" << std::endl;
    }
    // -- doublet splitting
    if (canceledFids.empty() && Simplifier::SHEET_SPLIT) {
        SheetSplitSimplifier sheetSplitSimplifier(mesh);
        sheetSplitSimplifier.Run(canceledFids);
        if (!canceledFids.empty()) std::cout << "remove_doublet from sheetSplitSimplifier" << std::endl;
    }

    // -- triplet splitting (optional)
    if (canceledFids.empty() && Simplifier::TRIP) {
        TriangleSimplifier triangleSimplifier(mesh);
        triangleSimplifier.Run(canceledFids);
        if (!canceledFids.empty()) std::cout << "collapse faces from TriangleSimplifier" << std::endl;
    }
    ////////////////////////////////////////////////////////////////////////////

    ////////////////// Boundary Conformality (Local Operations) ////////////////
    // -- edge rotation
    if (canceledFids.empty() && Simplifier::ROTATE) {
        EdgeRotateSimplifier edgeRotateSimplifier(mesh);
        edgeRotateSimplifier.RunCollective(canceledFids);
        if (!canceledFids.empty()) std::cout << "rotate_edge" << std::endl;
    }
    
    ////////////////////////////////////////////////////////////////////////////

    /////// Simplification Operations (Local and Semi-global Operations) ///////
    
    
    
    // if (canceledFids.empty()) {
    //     GetOperations(SimplificationOps, "Diagonal_Collapse");
    //     // GetOperations(SimplificationOps, "Edge_Collapse");
    //     GetOperations(SimplificationOps, "Edge_Rotate");
    //     // GetOperations(SimplificationOps, "Vertex_Rotate");
    //     if (!SimplificationOps.empty()) {
    //         PerformOperations(SimplificationOps, canceledFids);
    //     }
    // }

    if (canceledFids.empty()) {
        GetOperations(SimplificationOps, "Strict_Diagonal_Collapse");
        // GetOperations(SimplificationOps, "Separatrix_Collapse");
        // GetOperations(SimplificationOps, "Half_Separatrix_Collapse");
        // GetOperations(SimplificationOps, "Chord_Collapse");
        if (!SimplificationOps.empty()) {
            PerformOperations(SimplificationOps, canceledFids);
        }
    }

    if (canceledFids.empty()) {
        GetOperations(SimplificationOps, "Separatrix_Collapse");
        if (!SimplificationOps.empty()) {
            PerformOperations(SimplificationOps, canceledFids);
        }
    }

    if (canceledFids.empty()) {
        GetOperations(SimplificationOps, "Loose_Separatrix_Collapse");
        if (!SimplificationOps.empty()) {
            PerformOperations(SimplificationOps, canceledFids);
        }
    }

    if (canceledFids.empty()) {
        GetOperations(SimplificationOps, "Chord_Collapse");
        if (!SimplificationOps.empty()) {
            PerformOperations(SimplificationOps, canceledFids);
        }
    }

    // if (canceledFids.empty() && Simplifier::GLOBAL) {
    //    SingleSheetSimplifier sheetSimplifier(mesh);
    // //    sheetSimplifier.Run(canceledFids);
    //    sheetSimplifier.ExtractAndCollapse(canceledFids);
    // }
    
    
    if (canceledFids.empty()) {
        GetOperations(SimplificationOps, "Half_Separatrix_Collapse");
        if (!SimplificationOps.empty()) {
            PerformOperations(SimplificationOps, canceledFids);
        }
    }

    ////////////////////////////////////////////////////////////////////////////

    if (canceledFids.empty()) {
        update(canceledFids);
        return false;
    }

    std::set<size_t> tempVids;
    for (auto fid: canceledFids) {
        auto& f = mesh.F.at(fid);
        tempVids.insert(f.Vids.begin(), f.Vids.end());
        // for (auto nfid: f.N_Fids) {
        //     tempVids.insert(mesh.F.at(nfid).Vids.begin(), mesh.F.at(nfid).Vids.end());
        // } 
    }
    for (auto vid: tempVids) mesh.V.at(vid).smoothLocal = true;

    update(canceledFids);
    if (Simplifier::writeFile) {
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

void PatchSimplifier::GetOperations(std::multiset<SimplificationOperationStruct, bool(*)(SimplificationOperationStruct, SimplificationOperationStruct)>& SimplificationOps, std::string OpType) {
    BaseComplexQuad baseComplex(mesh);

    if (OpType == "Edge_Rotate") {
        EdgeRotate(SimplificationOps);
    }
    
    if (OpType == "Vertex_Rotate") {
        VertexRotate(SimplificationOps);
    }

    if (OpType == "Edge_Collapse") {
        EdgeCollapse(SimplificationOps);
    }

    if (OpType == "Diagonal_Collapse") {
        DiagonalCollapse(SimplificationOps);
    }

    if (OpType == "Strict_Diagonal_Collapse" && Simplifier::COLLAPSE_DIAGNAL) {
        DiagnalCollapseSimplifier diagonalCollapseSimplifier(mesh);
        diagonalCollapseSimplifier.GetDiagonalCollapseOps(SimplificationOps);
    }

    if (OpType == "Separatrix_Collapse" && Simplifier::COLLAPSE) {
        baseComplex.ExtractSingularVandE();
        baseComplex.BuildE();
        GetSeparatrixCollapseOps(baseComplex, false, SimplificationOps);
        // GetSeparatrixCollapseOps(baseComplex, true, SimplificationOps);
    }

    if (OpType == "Loose_Separatrix_Collapse" && Simplifier::COLLAPSE) {
        baseComplex.ExtractSingularVandE();
        baseComplex.BuildE();
        GetSeparatrixCollapseOps(baseComplex, true, SimplificationOps);
    }

    if (OpType == "Half_Separatrix_Collapse" && Simplifier::HALF) {
        baseComplex.ExtractSingularVandE();
        baseComplex.BuildE();
        GetHalfSeparatrixOps(baseComplex, SimplificationOps);
    }

    if (OpType == "Chord_Collapse" && Simplifier::GLOBAL) {
        baseComplex.Build();
        SheetSimplifier sheetSimplifier(mesh);
        sheetSimplifier.GetChordCollapseOps(baseComplex, SimplificationOps);
    }
    
}

void PatchSimplifier::PerformOperations(std::multiset<SimplificationOperationStruct, bool(*)(SimplificationOperationStruct, SimplificationOperationStruct)>& SimplificationOps, std::set<size_t>& canceledFids) {
    std::set<size_t> processedFids;
    for(auto op: SimplificationOps) {
        bool negativeFacePresent = false;
        for (auto& f: op.newFaces) {
            if (GetScaledJacobianQuad(mesh, f) < 0.0) {
                negativeFacePresent = true;
                break;
            }
            for (auto vid: f.Vids) {
                auto& v = mesh.V.at(vid);
                for (auto fid: v.N_Fids) {
                    if (GetScaledJacobianQuad(mesh, mesh.F.at(fid)) < 0.0) {
                        negativeFacePresent = true;
                        break;
                    }
                }
                if (negativeFacePresent) {
                    break;
                }
            }
            if (negativeFacePresent) {
                break;
            }
        }
        if (negativeFacePresent) {
            // std::cout << op.type << std::endl;
            // std::cout << "======================" << std::endl;
            continue;
        }
        bool overlappedOp = false;
        for (auto id: op.canceledFids) {
            if (std::find(processedFids.begin(), processedFids.end(), id) != processedFids.end()) {
                overlappedOp = true;
                break;
            }
        }
        if (overlappedOp) continue;
        std::cout << op.type << std::endl;
        for (auto f: op.newFaces) {
            f.id = mesh.F.size();
            mesh.F.push_back(f);
            for (auto vid: f.Vids) mesh.V.at(vid).smoothLocal = true;
        }
        // for (int i = 0; i < op.updateVertexIds.size(); i++) {
        //     if (mesh.V.at(op.updateVertexIds.at(i)).type > FEATURE) continue;
        //     mesh.V.at(op.updateVertexIds.at(i)) = op.updatedVertexPos.at(i);
        // }
        canceledFids.insert(op.canceledFids.begin(), op.canceledFids.end());
        processedFids.insert(op.canceledFids.begin(), op.canceledFids.end());
        
        for (auto fid: op.canceledFids) {
            Face& f = mesh.F.at(fid);
            processedFids.insert(f.N_Fids.begin(), f.N_Fids.end());
            for (auto n_fid: f.N_Fids) {
                Face& nf = mesh.F.at(n_fid);
                processedFids.insert(nf.N_Fids.begin(), nf.N_Fids.end());
            }
        }
    }
}

bool PatchSimplifier::Simplify(int& iter) {
    std::set<size_t> canceledFids;
    init();
    
    if (iter == 0 && featurePreserved) {
        get_feature();
        std::cout << "After get_feature" << std::endl;
        // refinedV.insert(refinedV.begin(), origMesh.V.begin(), origMesh.V.end());
        // refineVertexInFaces(origMesh, refinedV, 1);
	    // refineVertexInEdges(origMesh, refinedV, 1);
        for (auto& v: origMesh.V) {
            if (v.isBoundary && !v.isCorner) {
                origBoundaryVids.push_back(v.id);
            }
        }
    }
    if (iter == 0)
    {
        auto eids = get_rotate_eids();
        MeshFileWriter writer(mesh, "rotate_eids.vtk");
        writer.WriteEdgesVtk(eids);
        smoothGlobal = true;
        std::cout << "Begin simplification" << std::endl;
        SmoothMesh();
    }
    else {
        smoothGlobal = false;
        // std::cout << "Before smoothing local" << std::endl;
        SmoothMesh();
        // std::cout << "After smoothing local" << std::endl;
        // for (auto& v: mesh.V) {
        //     std::cout << v.x << " " << v.y << " " << v.z << std::endl; 
        // }
        // SurfaceMapper sm(mesh, origMesh);
        // sm.Map();
        // smooth_project();
    }
    /*if (mesh.F.size() <= 0.25 * origMesh.F.size()) {
    //     std::cout << "current Faces: " << mesh.F.size() << " original Faces: " << origMesh.F.size() << std::endl;
        
    //     smoothGlobal = true;
    //     SmoothMesh();
        // RefineMesh();
        // SmoothMesh(true);
        mesh = RefineWithFeaturePreserved(mesh, 0);
        size_t id = 0;
        for (auto& c : mesh.C) {
            c.cellType = VTK_QUAD;
            c.id = id++;
        }
        for (auto& v : mesh.V) {
	        v.N_Vids.clear();
            v.N_Eids.clear();
            v.N_Fids.clear();
            v.N_Cids.clear();
	    }
        init();
        std::cout << "Refined Mesh" << std::endl;
        smoothGlobal = true;
        SmoothMesh();
    }*/
    
   

    // Step 1 -- doublet removal
    if (canceledFids.empty() && Simplifier::REMOVE_DOUBLET) {
        DoubletSimplifier doubletSimplifier(mesh);
        // doubletSimplifier.Run(canceledFids);
        doubletSimplifier.RunCollective(canceledFids);
        if (!canceledFids.empty()) std::cout << "remove_doublet" << std::endl;
    }
    // Step 2 -- doublet splitting
    // if (canceledFids.empty() && Simplifier::SHEET_SPLIT) {
    //     SheetSplitSimplifier sheetSplitSimplifier(mesh);
    //     sheetSplitSimplifier.Run(canceledFids);
    //     if (!canceledFids.empty()) std::cout << "remove_doublet from sheetSplitSimplifier" << std::endl;
    // }
    // Step 4 -- singlet collapsing
    if (canceledFids.empty()) {
        DiagnalCollapseSimplifier diagnalCollapseSimplifier(mesh);
        diagnalCollapseSimplifier.Run3(canceledFids);
        // diagnalCollapseSimplifier.CollapseSinglets(canceledFids);
        if (!canceledFids.empty()) std::cout << "singlet collapsing" << std::endl;
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
        // edgeRotateSimplifier.Run(canceledFids);
        edgeRotateSimplifier.RunCollective(canceledFids);
        if (!canceledFids.empty()) std::cout << "rotate_edge" << std::endl;
    }
    
    // static bool aligned = false;
    // if (canceledFids.empty() && !aligned) {
    //     aligned = true;
    //     std::cout << "writing rotate.vtk " << std::endl;
    //     MeshFileWriter writer(mesh, "rotate.vtk");
    //     writer.WriteFile();
    // }    
    
         // Step 8 -- diagonal collapsing
    if (canceledFids.empty() && Simplifier::COLLAPSE_DIAGNAL) {
        DiagnalCollapseSimplifier diagnalCollapseSimplifier(mesh);
        // diagnalCollapseSimplifier.Run(canceledFids);
        diagnalCollapseSimplifier.RunCollective(canceledFids);
        if (!canceledFids.empty()) std::cout << "collapse_diagnal" << std::endl;
    }


    // Step 5 -- <separatrix splitting> and <separatrix splitting (optional)>
    if (canceledFids.empty() && (Simplifier::COLLAPSE || Simplifier::SPLIT)) {
        BaseComplexQuad baseComplex(mesh);
        baseComplex.ExtractSingularVandE();
        baseComplex.BuildE();
        // strict_simplify(baseComplex, canceledFids);
        three_connections_collapse(baseComplex, canceledFids, false);
        if (canceledFids.empty()) {
            // loose_simplify(baseComplex, canceledFids);
            three_connections_collapse(baseComplex, canceledFids, true);
            // loose_simplify_random(baseComplex, canceledFids);
            if (!canceledFids.empty()) std::cout << "loose_simplify\n";
//            else if (canceledFids.empty() && Simplifier::HALF) {
//                half_simplify(baseComplex, canceledFids);
//                if (!canceledFids.empty()) std::cout << "half_simplify\n";
//            }
        } else std::cout << "strict_simplify\n";
    }

    // if (canceledFids.empty() && Simplifier::GLOBAL) {
    //     // global_simplify(canceledFids);
    //    SheetSimplifier sheetSimplifier(mesh);
    // //    sheetSimplifier.Run(canceledFids);
    //    sheetSimplifier.ExtractAndCollapse(canceledFids);
    //     if (!canceledFids.empty()) std::cout << "chord_collapse" << std::endl;

    // }
    
    if (canceledFids.empty() && Simplifier::HALF) {
       BaseComplexQuad baseComplex(mesh);
       baseComplex.ExtractSingularVandE();
       baseComplex.BuildE();
    //    half_simplify(baseComplex, canceledFids);
        half_separatrix_collapse(baseComplex, canceledFids);

        if (!canceledFids.empty()) std::cout << "half_simplify\n";
    }

    if (canceledFids.empty() && Simplifier::GLOBAL) {
        SingleSheetSimplifier sheetSimplifier(mesh);
    //    sheetSimplifier.Run(canceledFids);
        sheetSimplifier.ExtractAndCollapse(canceledFids);
    }

   

    if (canceledFids.empty()) {
        update(canceledFids);
        return false;
    }
    std::set<size_t> tempVids;
    for (auto fid: canceledFids) {
        auto& f = mesh.F.at(fid);
        tempVids.insert(f.Vids.begin(), f.Vids.end());
        for (auto nfid: f.N_Fids) {
            tempVids.insert(mesh.F.at(nfid).Vids.begin(), mesh.F.at(nfid).Vids.end());
            // for (auto nnfid: mesh.F.at(nfid).N_Fids) tempVids.insert(mesh.F.at(nnfid).Vids.begin(), mesh.F.at(nnfid).Vids.end());
        } 
    }
    for (auto vid: tempVids) mesh.V.at(vid).smoothLocal = true;
    // std::cout << tempVids.size() << std::endl;
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

/*void PatchSimplifier::SmoothMesh(bool SmoothGlobal_) {
    int n = 100;
    for (auto& v: mesh.V) {
        // if (n < 0) break;
        if (v.N_Vids.empty() || v.N_Eids.empty() || v.N_Fids.empty()) continue;
        if (SmoothGlobal_) {
            smoothVids.push_back(v.id);
        } else {
            if (v.smoothLocal) smoothVids.push_back(v.id);
        }
        n -= 1;
    }
    int iters = 1;
    while (iters--) {
        for (auto vid: smoothVids) {
            auto& v = mesh.V.at(vid);
            if (v.type == FEATURE || v.isBoundary) {
                SetPositionBoundary(v);
            } else {
                SetPosition(v);
            }
        }
    }    
    for (auto& v: mesh.V) v.smoothLocal = false;
    smoothVids.clear();
}*/

void PatchSimplifier::SetPosition(Vertex& v) {
    double polyArea = 0.0;
    glm::dvec3 centroid(0.0, 0.0, 0.0);
    size_t startE = v.N_Eids.at(0);
    
    std::cout << "Setting Postion for vertex: " << v.id << " nvids: " << v.N_Vids.size() << " neids: " << v.N_Eids.size() << " nfids: " << v.N_Fids.size() << std::endl;
    for (int i = 0; i < v.N_Eids.size(); i++) {
        auto& edge = mesh.E.at(startE);
        size_t ev = edge.Vids.at(0) == v.id ? edge.Vids.at(1) : edge.Vids.at(0);
        size_t ev_plus1;
        size_t ev_minus1;
        std::cout << "edge nfids: " << edge.N_Fids.size() << std::endl;
        for (auto fid: edge.N_Fids) {
            auto& f = mesh.F.at(fid);
            std::cout << "face vids: " << f.Vids.at(0) << " " << f.Vids.at(1) << " " << f.Vids.at(2) << " " << f.Vids.at(3) << std::endl;
            int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
            if (f.Vids.at((idx+1)%f.Vids.size()) == ev) {
                ev_plus1 = f.Vids.at((idx+3)%f.Vids.size());
                startE = mu.GetDifference(mu.GetIntersection(f.Eids, v.N_Eids), std::vector<size_t>{edge.id}).at(0);
            }
            if (f.Vids.at((idx+3)%f.Vids.size()) == ev) {
                ev_minus1 = f.Vids.at((idx+1)%f.Vids.size());
            }
        }

        std::cout << "vid: " << v.id << " ev: " << ev << " ev_plus1: " << ev_plus1 << " ev_minus1: " << ev_minus1 << std::endl;
        auto& v2 = mesh.V.at(ev);
        auto& v3 = mesh.V.at(ev_plus1);
        auto& v4 = mesh.V.at(ev_minus1);

 
        std::cout << "v: " << v.x << " " << v.y << " " << v.z << std::endl;
        std::cout << "v2: " << v2.x << " " << v2.y << " " << v2.z << std::endl;
        std::cout << "v3: " << v3.x << " " << v3.y << " " << v3.z << std::endl;
        std::cout << "v4: " << v4.x << " " << v4.y << " " << v4.z << std::endl;
        glm::dvec3 AB = v3.xyz() - v2.xyz();
        glm::dvec3 BC = v4.xyz() - v3.xyz();
        glm::dvec3 CA = v2.xyz() - v4.xyz();
        glm::dvec3 AC = v4.xyz() - v2.xyz();

        double a = glm::length(BC);
        double b = glm::length(CA);
        double c = glm::length(AB);
        glm::dvec3 incenter = ((a * v2.xyz()) + (b * v3.xyz()) + (c * v4.xyz())) / (a + b + c);
        
        double area = 0.5 * glm::length(glm::cross(AB, AC));
        centroid += (area * incenter); 
        polyArea += area;
    }
    std::cout << "centroid: " << centroid.x << " " << centroid.y << " " << centroid.z << std::endl;
    std::cout << "polyArea: " << polyArea << std::endl;
    centroid /= polyArea;
    std::cout << "centroid: " << centroid.x << " " << centroid.y << " " << centroid.z << std::endl;
    v.x = centroid.x;
    v.y = centroid.y;
    v.z = centroid.z;
    // std::cout << "**************************" << std::endl;
    /*for (auto fid: v.N_Fids) {
        auto& f = mesh_.F.at(fid);
        int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
        if (idx == -1) continue;
        auto& v2 = mesh_.V.at(f.Vids.at((idx+1)%f.Vids.size()));
        auto& v3 = mesh_.V.at(f.Vids.at((idx+3)%f.Vids.size()));

        glm::dvec3 AB = v2.xyz() - v.xyz();
        glm::dvec3 AC = v3.xyz() - v.xyz();
        glm::dvec3 BC = v3.xyz() - v2.xyz();
        glm::dvec3 CA = v.xyz() - v3.xyz();

        glm::dvec3 T_cross = glm::cross(AB, AC);

        glm::dvec3 normal = glm::normalize(T_cross);
        glm::dvec3 temp = centroid - v.xyz();
        double dist = glm::dot(temp, normal);
        glm::dvec3 projected_point = centroid - (dist * normal);
        glm::dvec3 AP = projected_point - v.xyz();
        glm::dvec3 BP = projected_point - v2.xyz();

        double T_area = 0.5 * glm::length(T_cross);
        if (0.5 * (glm::length(glm::cross(AB, AP)) + glm::length(glm::cross(AC, AP)), glm::length(glm::cross(BP, BC))) > T_area) continue;
        v.x = projected_point.x;
        v.y = projected_point.y;
        v.z = projected_point.z;
        break;
    }*/
}

void PatchSimplifier::SetPositionBoundary(Vertex& v) {
    if (v.isCorner || v.N_Fids.size() < 2 || v.N_Fids.size() > 2) return;
    std::vector<size_t> boundaryVertices;
    for (auto vid: v.N_Vids) {
        if (mesh.V.at(vid).type == FEATURE || mesh.V.at(vid).isBoundary) boundaryVertices.push_back(vid);
    }
    auto& b1 = mesh.V.at(boundaryVertices.at(0));
    auto& b2 = mesh.V.at(boundaryVertices.at(1));
    std::cout << "b1: " << b1.x << " " << b1.y << " " << b1.z << std::endl;
    std::cout << "b2: " << b2.x << " " << b2.y << " " << b2.z << std::endl;
    glm::dvec3 vb1 = glm::normalize(b1.xyz() - v.xyz());
    glm::dvec3 vb2 = glm::normalize(b2.xyz() - v.xyz());

    double angle = acos(glm::dot(vb1, vb2)) * 180.0 / PI;
    glm::dvec3 centroid(0.0, 0.0, 0.0);
    int k = 0;
    for (auto vid: v.N_Vids) {
        auto& nv = mesh.V.at(vid);
        if (vid == b1.id || vid == b2.id) continue;

        k += 1;
        glm::dvec3 A(v.xyz() - nv.xyz());
        glm::dvec3 B(b1.xyz() - nv.xyz());
        glm::dvec3 C(b2.xyz() - nv.xyz());

        double a1 = glm::dot(A, B) / (glm::length(A) * glm::length(B));
        double a2 = glm::dot(A, C) / (glm::length(A) * glm::length(C));
        if (a1 < -1.0) a1 = -1.0;
        if (a1 > 1.0) a1 = 1.0;
        if (a2 < -1.0) a2 = -1.0;
        if (a2 > 1.0) a2 = 1.0;

        double alpha1 = acos(a1);
        double alpha2 = acos(a2);

        double beta = (alpha2 - alpha1) / 2;

        glm::dvec3 r1 = v.xyz();
        double l = 0;
        if (beta > 0) {
            r1 = b2.xyz() - v.xyz();
            l = fabs(beta / alpha2) * glm::length(r1); 
        } else if (beta < 0) {
            r1 = b1.xyz() - v.xyz();
            l = fabs(beta / alpha1) * glm::length(r1);
        }
        r1 = glm::normalize(r1);
        centroid += (v.xyz() + (l*r1));
    }
    std::cout << "centroid: " << centroid.x << " " << centroid.y << " " << centroid.z << std::endl;
    std::cout << "k: " << k << std::endl;
    if (k > 0) {
        centroid /= (double) k;
        std::cout << "centroid: " << centroid.x << " " << centroid.y << " " << centroid.z << std::endl;
        v.x = centroid.x;
        v.y = centroid.y;
        v.z = centroid.z;
    }
}

void PatchSimplifier::SmoothMesh() {
    mesh.SetOneRingNeighborhood();
    for (auto& v: mesh.V) {
        if (smoothGlobal) {
            smoothVids.push_back(v.id);
        } else {
            if (v.smoothLocal) smoothVids.push_back(v.id);
        }
    }
    std::vector<glm::dvec3> current_coords(smoothVids.size(), glm::dvec3(0.0, 0.0, 0.0));
    for (int i = 0; i < smoothVids.size(); i++) {
        auto& v = mesh.V.at(smoothVids.at(i));
        if (!v.isBoundary || v.isCorner) continue;
        current_coords.at(i) = v.xyz();
    }
    RemapBoundaryVertices(current_coords);
    for (int i = 0; i < smoothVids.size(); i++) {
        Vertex& v = mesh.V.at(smoothVids.at(i));
        if (!v.isBoundary || v.isCorner) continue;
        v = current_coords.at(i);
    }

    int iters = 10;
    while (iters--) {
        std::vector<glm::dvec3> delta_coords;
        delta_coords.resize(smoothVids.size(), glm::dvec3(0.0, 0.0, 0.0));
        double currentE = GetMeshEnergy();
        AngleBasedSmoothing(delta_coords);
        for (int i = 0; i < smoothVids.size(); i++) {
            Vertex& v = mesh.V.at(smoothVids.at(i));
            if (v.isBoundary) continue;
            glm::dvec3 temp_coord = v.xyz();
            v = delta_coords.at(i);
            delta_coords.at(i) = temp_coord;
        }
        double newE = GetMeshEnergy();
        // std::cout << "interior iter: " << iters << " oldE: " << currentE << " newE: " << newE << std::endl;
        if (currentE - newE < 1e-4) {
            for (int i = 0; i < smoothVids.size(); i++) {
                Vertex& v = mesh.V.at(smoothVids.at(i));
                if (v.isBoundary) continue;
                v = delta_coords.at(i);
            }
            break;
        }
        SmoothBoundary();
    }
    for (auto& v: mesh.V) v.smoothLocal = false;
    smoothVids.clear();
    // mesh.E.clear();
    // mesh.F.clear();
    if (smoothGlobal) smoothGlobal = false;
    std::set<size_t> bogus;
    update(bogus);
    init();
}

void PatchSimplifier::SmoothBoundary() {
    int iters = 10;
    while (iters--) {
        std::vector<glm::dvec3> delta_coords;
        delta_coords.resize(smoothVids.size(), glm::dvec3(0.0, 0.0, 0.0));
        double currentE = GetMeshEnergy();
        ResampleBoundaryVertices(delta_coords);
        RemapBoundaryVertices(delta_coords);
        for (int i = 0; i < smoothVids.size(); i++) {
            Vertex& v = mesh.V.at(smoothVids.at(i));
            if (!v.isBoundary || v.isCorner) continue;
            glm::dvec3 temp_coord = v.xyz();
            v = delta_coords.at(i);
            delta_coords.at(i) = temp_coord;
        }
        // RemapBoundaryVertices(delta_coords);
        double newE = GetMeshEnergy();
        // std::cout << "boundary iter: " << iters << " oldE: " << currentE << " newE: " << newE << std::endl;
        if (currentE - newE < 1e-4) {
            // RemapBoundaryVertices(delta_coords);
            for (int i = 0; i < smoothVids.size(); i++) {
                Vertex& v = mesh.V.at(smoothVids.at(i));
                if (!v.isBoundary || v.isCorner) continue;
                v = delta_coords.at(i);
            }
            break;
        }
    }
}

void PatchSimplifier::AngleBasedSmoothing(std::vector<glm::dvec3>& delta_coords) {
    std::vector<double> vertex_weights(smoothVids.size(), 0);
    std::vector<glm::dvec3> new_coords(smoothVids.size(), glm::dvec3(0.0, 0.0, 0.0));

    for (int i = 0; i < smoothVids.size(); i++) {
        Vertex& v_i = mesh.V.at(smoothVids.at(i));
        if (v_i.isBoundary) continue;
        std::vector<size_t> one_ring_neighbors = v_i.N_Vids;
        std::vector<glm::dvec3> neighbor_coords(one_ring_neighbors.size(), glm::dvec3(0.0, 0.0, 0.0));
        std::vector<double> neighbor_weights(one_ring_neighbors.size(), 0.0);
        double weight_agg = 0.0;
        int k = one_ring_neighbors.size();
        for (int j = 0; j < one_ring_neighbors.size(); j++) {
            int a = j;
            int b = one_ring_neighbors.size();
            int index = one_ring_neighbors.at((a % b + b) % b);
            Vertex& v_j = mesh.V.at(index);

            a = j - 1;
            index = one_ring_neighbors.at((a % b + b) % b);
            Vertex& v_j_prev = mesh.V.at(index);

            a = j + 1;
            index = one_ring_neighbors.at((a % b + b) % b);
            Vertex& v_j_next = mesh.V.at(index);

            glm::dvec3 V_j = v_i.xyz() - v_j.xyz();
            glm::dvec3 V_j_minus_1 = v_j_prev.xyz() - v_j.xyz();
            glm::dvec3 V_j_plus_1 = v_j_next.xyz() - v_j.xyz();

            double alpha1 = acos(glm::dot(V_j, V_j_plus_1) / (glm::length(V_j) * glm::length(V_j_plus_1)));
            double alpha2 = acos(glm::dot(V_j, V_j_minus_1) / (glm::length(V_j) * glm::length(V_j_minus_1)));
            double beta = (alpha2 - alpha1) / 2;
            
            neighbor_coords.at(j).x = (v_j.x + (V_j.x * cos(beta)) - (V_j.y * sin(beta))); 
            neighbor_coords.at(j).y = (v_j.y + (V_j.x * sin(beta)) + (V_j.y * cos(beta)));
            neighbor_weights.at(j) = fabs(beta);
            weight_agg += neighbor_weights.at(j);
        }
        if (weight_agg == 0) {
            new_coords.at(i) = glm::dvec3(0, 0, 0);
            continue;
        }
        for (int j = 0; j < one_ring_neighbors.size(); j++) {
            glm::dvec3 current_v = neighbor_coords.at(j);
            double weight = (1 - neighbor_weights.at(j) / weight_agg);
            vertex_weights.at(i) += weight;
            new_coords.at(i) += (weight * current_v);
        }

        new_coords.at(i) = ((1 / vertex_weights.at(i)) * new_coords.at(i)) - v_i.xyz();
    }

    for (int i = 0; i < smoothVids.size(); i++) {
        Vertex& v = mesh.V.at(smoothVids.at(i));
        if (v.isBoundary) continue;
        delta_coords.at(i) += v.xyz() + new_coords.at(i);
    }
}

void PatchSimplifier::ResampleBoundaryVertices(std::vector<glm::dvec3>& delta_coords) {
    // bool hasNegativeElementsPresentAlready = false;
    // for (auto& v: mesh.V) {
    //     if (!v.isBoundary || v.isCorner) {
    //         continue;
    //     }
    //     for (auto fid: v.N_Fids) {
    //         if (IsFaceNegative(fid, v.id, glm::dvec3(v.x, v.y, 0.0))) {
    //             hasNegativeElementsPresentAlready = true;
    //             break;
    //         }
    //     }
    // }
    for (int i = 0; i <  smoothVids.size(); i++) {
        auto& v = mesh.V.at(smoothVids.at(i));
        if (!v.isBoundary || v.isCorner) continue;
        int nb_n_index = -1;
        int k = 0;
        std::vector<size_t> neighbors = v.N_Vids;
        if (v.N_Fids.size() < 2) {
            neighbors = v.oneRingNeighborVertices;
        }
        if (v.N_Fids.size() > 2) {
            delta_coords.at(i) = v.xyz();
            continue;
        }
        for (int j = 0; j < neighbors.size(); j++) {
            if (!mesh.V.at(neighbors.at(j)).isBoundary) {
                nb_n_index = j;
                int a_ = nb_n_index;
                int b_ = neighbors.size();
                int index = (a_ % b_ + b_) % b_;
                Vertex& b0 = mesh.V.at(neighbors.at(index));
                
                k += 1;
                a_ = nb_n_index - 1;
                index = (a_ % b_ + b_) % b_;
                Vertex& b1 = mesh.V.at(neighbors.at(index));

                a_ = nb_n_index + 1;
                index = (a_ % b_ + b_) % b_;
                Vertex& b2 = mesh.V.at(neighbors.at(index));

                glm::dvec3 V_j = v.xyz() - b0.xyz();
                glm::dvec3 V_j_minus_1 = b1.xyz() - b0.xyz();
                glm::dvec3 V_j_plus_1 = b2.xyz() - b0.xyz();

                double a1 = glm::dot(V_j, V_j_plus_1) / (glm::length(V_j) * glm::length(V_j_plus_1));
                double a2 = glm::dot(V_j, V_j_minus_1) / (glm::length(V_j) * glm::length(V_j_minus_1));
                if (a1 < -1.0) a1 = -1.0;
                if (a1 > 1.0) a1 = 1.0;
                if (a2 < -1.0) a2 = -1.0;
                if (a2 > 1.0) a2 = 1.0;
                
                double alpha1 = acos(a1);
                double alpha2 = acos(a2);

                double beta = (alpha2 - alpha1) / 2;

                glm::dvec3 r1 = v.xyz();
                double l = 0;
                if (beta > 0) {
                    r1 = b1.xyz() - v.xyz();
                    l = fabs(beta / alpha2) * glm::length(r1); 
                } else if (beta < 0) {
                    r1 = b2.xyz() - v.xyz();
                    l = fabs(beta / alpha1) * glm::length(r1);
                }
                r1 = glm::normalize(r1);
                glm::dvec3 new_coords = v.xyz() + (l * r1);
                
                // bool isNegativeElementPresent = false;
                // for (auto fid: v.N_Fids) {
                //     if (IsFaceNegative(fid, v.id, new_coords)) {
                //         isNegativeElementPresent = true;
                //     }
                // }
                // if (!hasNegativeElementsPresentAlready && isNegativeElementPresent) {
                //     new_coords = v.xyz();
                // }
                delta_coords.at(i) += new_coords;
            }           
        }
        if (k > 0) {
            delta_coords.at(i) = ((double) 1/k) * delta_coords.at(i);
        } else {
            delta_coords.at(i) = v.xyz();
        }
        // double currentE = GetVertexEnergy(v.id);
        // glm::dvec3 temp_coord = v.xyz();
        // v = delta_coords.at(i);
        // double newE = GetVertexEnergy(v.id);
        // v = temp_coord;
        // if (newE >= currentE) {
        //     delta_coords.at(i) = v.xyz();
        // }     
    }
}

void PatchSimplifier::RemapBoundaryVertices(std::vector<glm::dvec3>& delta_coords) {
    for (int i = 0; i < smoothVids.size(); i++) {
        if (!mesh.V.at(smoothVids.at(i)).isBoundary || mesh.V.at(smoothVids.at(i)).isCorner) continue;
        
        glm::dvec3 new_v = delta_coords.at(i);
        double min_length = std::numeric_limits<double>::max();
        int min_index = -1;
        for (int j = 0; j < origBoundaryVids.size(); j++) {
            double length = glm::length(new_v - origMesh.V.at(origBoundaryVids.at(j)).xyz());
            if (length < min_length) {
                min_length = length;
                min_index = origBoundaryVids.at(j);
            }
        }
        if (min_index == -1) {
            continue;
        }
        Vertex& v = origMesh.V.at(min_index);
        std::vector<size_t> boundary_neighbors;
        for (int j = 0; j < v.N_Vids.size(); j++) {
            if (origMesh.V.at(v.N_Vids.at(j)).isBoundary) {
                boundary_neighbors.push_back(v.N_Vids.at(j));
            }                    
        }
        Vertex& b1 = origMesh.V.at(boundary_neighbors.at(0));
        Vertex& b2 = origMesh.V.at(boundary_neighbors.at(1));
        
        glm::dvec3 a = new_v - v.xyz();
        glm::dvec3 b = b1.xyz() - v.xyz();
        glm::dvec3 c = b2.xyz() - v.xyz();

        double length_a = glm::dot(a, a);
        double length_b = glm::dot(b, b);
        double length_c = glm::dot(c, c);

        double dot_a_b = glm::dot(a, b) / glm::length(b);
        double dot_a_c = glm::dot(a, c) / glm::length(c);
        
        glm::dvec3 new_coords = new_v;
        if (dot_a_b >= 0) {
            b = glm::normalize(b);
            new_coords = v.xyz() + (dot_a_b * b);
        } else if (dot_a_c >= 0){
            c = glm::normalize(c);
            new_coords = v.xyz() + (dot_a_c * c);
        }
        
        // bool isNegativeElementPresent = false;
        // for (auto fid: mesh.V.at(i).N_Fids) {
        //     if (IsFaceNegative(fid, mesh.V.at(i).id, new_coords)) {
        //         isNegativeElementPresent = true;
        //     }
        // }
        // if (isNegativeElementPresent) {
        //     delta_coords.at(i) = mesh.V.at(i).xyz();
        //     continue;
        // }
        delta_coords.at(i) = new_coords;
        // double currentE = GetVertexEnergy(mesh.V.at(i).id);
        // glm::dvec3 temp_coord = mesh.V.at(i);
        // mesh.V.at(i) = new_coords;
        // double newE = GetVertexEnergy(mesh.V.at(i).id);
        // mesh.V.at(i) = temp_coord;
        // if (newE >= currentE) {
        //     delta_coords.at(i) = mesh.V.at(i).xyz();
        // }
    }
}

double PatchSimplifier::GetMeshEnergy() {
    double E = 0.0;
    for (int i = 0; i < smoothVids.size(); i++) {
        auto& v = mesh.V.at(smoothVids.at(i));
        E += GetVertexEnergy(v.id);
    }
    return E;
}

double PatchSimplifier::GetVertexEnergy(int vid) {
    Vertex& v = mesh.V.at(vid);
    double E = 0.0;
    std::vector<size_t> neighbors = v.N_Vids;
    for (int j = 0; j < neighbors.size(); j++) {
        int a = j;
        int b = neighbors.size();
        int index = neighbors.at((a % b + b) % b);
        Vertex& v_j = mesh.V.at(index);
        if (v.isBoundary && v_j.isBoundary) {
            continue;
        }

        a = j - 1;
        index = neighbors.at((a % b + b) % b);
        Vertex& v_j_prev = mesh.V.at(index);

        a = j + 1;
        index = neighbors.at((a % b + b) % b);
        Vertex& v_j_next = mesh.V.at(index);

        glm::dvec3 V_j = glm::normalize(v.xyz() - v_j.xyz());
        glm::dvec3 V_j_minus_1 = glm::normalize(v_j_prev.xyz() - v_j.xyz());
        glm::dvec3 V_j_plus_1 = glm::normalize(v_j_next.xyz() - v_j.xyz());

        double a1 = glm::dot(V_j, V_j_plus_1) / (glm::length(V_j) * glm::length(V_j_plus_1));
        double a2 = glm::dot(V_j, V_j_minus_1) / (glm::length(V_j) * glm::length(V_j_minus_1));
        if (a1 < -1.0) a1 = -1.0;
        if (a1 > 1.0) a1 = 1.0;
        if (a2 < -1.0) a2 = -1.0;
        if (a2 > 1.0) a2 = 1.0;

        double alpha1 = acos(a1);
        double alpha2 = acos(a2);
        E += ((alpha1 * alpha1) / 2) + ((alpha2 * alpha2) / 2);
    }
    return E;
}

void PatchSimplifier::RefineMesh() {
    std::vector<Vertex> newV(mesh.V.size());
	std::vector<Cell> newC;
	for (size_t i = 0; i < mesh.V.size(); ++i) {
	    auto& v = mesh.V.at(i);
	    auto& newv = newV.at(i);
		newv.id = i;
		newv = v.xyz();
		newv.type = v.type;
		newv.isCorner = v.isCorner;
		newv.isConvex = v.isConvex;
		newv.label = v.label;
		newv.patch_id = v.patch_id;
		newv.isSpecial = v.isSpecial;
		newv.labels = v.labels;
		newv.patch_ids = v.patch_ids;
		newv.idealValence = v.idealValence;
		newv.prescribed_length = v.prescribed_length;
		newv.smoothLocal = v.smoothLocal;
	}
    for (auto& e: mesh.E) {
        Vertex new_v;
        new_v = 0.5 * (mesh.V.at(e.Vids.at(0)).xyz() + mesh.V.at(e.Vids.at(1)).xyz());

        new_v.id = newV.size();
        for (auto fid: e.N_Fids) {
            Face& f = mesh.F.at(fid);
            for (int i = 0; i < f.Vids.size(); i++) {
                if (std::find(e.Vids.begin(), e.Vids.end(), f.Vids.at(i)) != e.Vids.end() &&
                    std::find(e.Vids.begin(), e.Vids.end(), f.Vids.at((i + 1) % f.Vids.size())) != e.Vids.end()) {
                        f.Vids.insert(f.Vids.begin() + (i + 1) % f.Vids.size(), new_v.id);
                        break;
                }
            }
        }
        newV.push_back(new_v);
    }
    for (int i = 0; i < mesh.F.size(); i++) {
        Face& f = mesh.F.at(i);

        Vertex new_v;
        new_v = 0.25 * (mesh.V.at(f.Vids.at(1)).xyz() + mesh.V.at(f.Vids.at(3)).xyz() + mesh.V.at(f.Vids.at(5)).xyz() + mesh.V.at(f.Vids.at(7)).xyz());
        new_v.id = newV.size();
        newV.push_back(new_v);

        Cell f1;
        f1.Vids = {new_v.id, f.Vids.at(0), f.Vids.at(1), f.Vids.at(2)};
        f1.cellType = VTK_QUAD;
        
        Cell f2;
        f2.Vids = {new_v.id, f.Vids.at(2), f.Vids.at(3), f.Vids.at(4)};
        f2.cellType = VTK_QUAD;

        Cell f3;
        f3.Vids = {new_v.id, f.Vids.at(4), f.Vids.at(5), f.Vids.at(6)};
        f3.cellType = VTK_QUAD;

        Cell f4;
        f4.Vids = {new_v.id, f.Vids.at(6), f.Vids.at(7), f.Vids.at(0)};
        f4.cellType = VTK_QUAD;

        f1.id = newC.size();
        newC.push_back(f1);

        f2.id = newC.size();
        newC.push_back(f2);

        f3.id = newC.size();
        newC.push_back(f3);

        f4.id = newC.size();
        newC.push_back(f4);
    }
    mesh.V.clear();
	mesh.E.clear();
	mesh.F.clear();
	mesh.C.clear();

	mesh.V = newV;
	mesh.C = newC;
    init();
}

void PatchSimplifier::SetOriginalRefinedMesh() {
    centerVertices.clear();
	centerVertices.insert(centerVertices.begin(), origMesh.V.begin(), origMesh.V.end());

    for (auto& f : origMesh.F) {
		glm::dvec3 center(0, 0, 0);
		for (auto nvid : origMesh.F.at(f.id).Vids)
			center += origMesh.V.at(nvid).xyz();
		center *= 0.25;
		Vertex vertex = center;
		vertex.id = centerVertices.size();
		vertex.patch_id = f.label;
		centerVertices.push_back(vertex);
	}

    for (auto& e : origMesh.E) {
		glm::dvec3 center(0, 0, 0);
		for (auto nvid : e.Vids)
			center += origMesh.V.at(nvid).xyz();
		center *= 0.5;
		Vertex vertex = center;
		vertex.id = centerVertices.size();

		if (!e.isSharpFeature) {
			vertex.patch_id = origMesh.V.at(e.Vids[0]).type == REGULAR ? origMesh.V.at(e.Vids[0]).patch_id : origMesh.V.at(e.Vids[1]).patch_id;
		} else {
			continue;
			vertex.label = origMesh.V.at(e.Vids[0]).label == MAXID ? origMesh.V.at(e.Vids[1]).label : origMesh.V.at(e.Vids[0]).label;
			vertex.type = FEATURE;
		}
		vertex.isBoundary = e.isBoundary;
		centerVertices.push_back(vertex);
	}

    for (auto& v : centerVertices)
		if (v.label != MAXID) origLabel_vids[v.label].insert(v.id);
	for (auto& v : centerVertices) {
		if (v.label != MAXID) {
			origLabel_vids[v.label].insert(v.id);
			for (auto& nvid : v.N_Vids) {
				auto& nv = centerVertices.at(nvid);
				auto& lineVids = origLabel_vids[v.label];
				if (nv.isCorner || lineVids.find(nvid) != lineVids.end()) origSharpEdgeVid_NVids[v.id].insert(nvid);
			}
		} else if (v.label == MAXID && !v.isCorner) {
			origPatch_vids[v.patch_id].insert(v.id);
		}
	}

    for (auto& v : mesh.V)
		if (v.label != MAXID) label_vids[v.label].insert(v.id);
	for (auto& v : mesh.V) {
		if (v.label != MAXID) {
			label_vids[v.label].insert(v.id);
			for (auto& nvid : v.N_Vids) {
				auto& nv = mesh.V.at(nvid);
				auto& lineVids = label_vids[v.label];
				if (nv.isCorner || lineVids.find(nvid) != lineVids.end()) sharpEdgeVid_NVids[v.id].insert(nvid);
			}
		}
	}
}

void PatchSimplifier::SmoothMesh(bool smoothGlobal_) {
    // smooth and project
    for (auto& v: mesh.V) {
        if (smoothGlobal_) {
            smoothVids.push_back(v.id);
        } else {
            if (v.smoothLocal) smoothVids.push_back(v.id);
        }
    }
	int iters = 10;

    /*int it1 = 0;
    while (it1++ < iters) {
        double currentE = GetMeshEnergy();
        std::vector<glm::dvec3> centers(mesh.V.size(), glm::dvec3(0.0, 0.0, 0.0));
        std::vector<glm::dvec3> temp(mesh.V.size(), glm::dvec3(0.0, 0.0, 0.0));
        for (auto vid: smoothVids) {
            auto& v = mesh.V.at(vid);
            if (v.type < FEATURE) {
                centers.at(v.id) = v.xyz();
                double n = 0.0;
                glm::dvec3 center(0.0, 0.0, 0.0);
                for (int i = 0; i < v.N_Eids.size(); i++) {
                    auto& e = mesh.E.at(v.N_Eids.at(i));
                    size_t nvid = e.Vids[0] != v.id ? e.Vids[0] : e.Vids[1];
                    std::vector<size_t> nvids;
                    for (auto fid: e.N_Fids) {
                        auto& f = mesh.F.at(fid);
                        for (auto fvid: f.Vids) {
                            if (fvid != v.id && fvid != nvid && std::find(v.N_Vids.begin(), v.N_Vids.end(), fvid) != v.N_Vids.end()) {
                                nvids.push_back(fvid);
                                break;
                            }
                        }
                    }
                    auto& v_a = mesh.V.at(nvid);
                    auto& v_b = mesh.V.at(nvids.at(0));
                    auto& v_c = mesh.V.at(nvids.at(1));
                    glm::dvec3 V_j = v.xyz() - v_a.xyz();
                    glm::dvec3 V_j_minus_1 = v_b.xyz() - v_a.xyz(); 
                    glm::dvec3 V_j_plus_1 = v_c.xyz() - v_a.xyz();

                    double r = glm::length(V_j);

                    glm::dvec3 p1 = v_a.xyz() + (r * glm::normalize(V_j_minus_1));
                    glm::dvec3 p2 = v_a.xyz() + (r * glm::normalize(V_j_plus_1));
                    if (glm::length((0.5 * (p1 + p2)) - v_a.xyz()) == 0) continue;
                    glm::dvec3 direction = glm::normalize((0.5 * (p1 + p2)) - v_a.xyz());
                    glm::dvec3 newPoint_p = v_a.xyz() + (r * direction);
                    glm::dvec3 newPoint_n = v_a.xyz() + (r * (-direction));
                    
                    center += glm::distance(v.xyz(), newPoint_p) < glm::distance(v.xyz(), newPoint_n) ? newPoint_p : newPoint_n;
                    n += 1;
                }
                centers.at(vid) = (center / n);
            }
        }
        for (auto vid: smoothVids) {
            auto& v = mesh.V.at(vid);
            if (v.type < FEATURE) {
                temp.at(v.id) = v.xyz();
                v = centers.at(v.id);
                bool negativeFacePresent = false;
                for (auto nfid: v.N_Fids) {
                    if (GetScaledJacobianQuad(mesh, mesh.F.at(nfid)) < 0) {
                        negativeFacePresent = true;
                        break;
                    }
                }
                if (negativeFacePresent) v = temp.at(v.id);
            }
        }
        double newE = GetMeshEnergy();
        if (currentE - newE < 1e-4) {
            for (auto vid: smoothVids) {
                auto& v = mesh.V.at(vid);
                if (v.type < FEATURE) {
                    v = temp.at(v.id);
                }
            }
            break;                
        }
        int it2 = 0;
        while (it2++ < iters) {
            currentE = GetMeshEnergy();
            for (auto vid: smoothVids) {
                auto& v = mesh.V.at(vid);
                if (v.type == FEATURE) {
                    centers.at(v.id) = v.xyz();
                    std::vector<size_t> neighbors = v.N_Vids;
                    std::vector<size_t> boundary_neighbors;
                    for (auto nvid: v.N_Vids) {
                        if (mesh.V.at(nvid).isBoundary) boundary_neighbors.push_back(nvid);
                    }
                    auto& v_b = mesh.V.at(boundary_neighbors[0]);
                    auto& v_c = mesh.V.at(boundary_neighbors[1]);
                    glm::dvec3 center(0, 0, 0);
                    double n = 0;
                    for (auto nvid: v.N_Vids) {
                        if (mesh.V.at(nvid).isBoundary) continue;
                        auto& v_a = mesh.V.at(nvid);
                        glm::dvec3 V_j = v.xyz() - v_a.xyz();
                        glm::dvec3 V_j_minus_1 = v_b.xyz() - v_a.xyz();
                        glm::dvec3 V_j_plus_1 = v_c.xyz() - v_a.xyz();

                        double r = glm::length(V_j);

                        glm::dvec3 p1 = v_a.xyz() + (r * glm::normalize(V_j_minus_1));
                        glm::dvec3 p2 = v_a.xyz() + (r * glm::normalize(V_j_plus_1));
                        if (glm::length((0.5 * (p1 + p2)) - v_a.xyz()) == 0) continue;
                        glm::dvec3 direction = glm::normalize((0.5 * (p1 + p2)) - v_a.xyz());
                        glm::dvec3 newPoint_p = v_a.xyz() + (r * direction);
                        glm::dvec3 newPoint_n = v_a.xyz() + (r * (-direction));
                        
                        center += glm::distance(v.xyz(), newPoint_p) < glm::distance(v.xyz(), newPoint_n) ? newPoint_p : newPoint_n;
                        n += 1;
                    }
                    if (n > 0) {
                        center /= n;
                        centers.at(v.id) = center;
                    }

                    double min_length = std::numeric_limits<double>::max();
                    int min_index = -1;
                    for (int i = 0; i < origBoundaryVids.size(); i++) {
                        double length = glm::length(centers.at(v.id) - origMesh.V.at(origBoundaryVids.at(i)).xyz());
                        if (length < min_length) {
                            min_length = length;
                            min_index = origBoundaryVids.at(i);
                        }
                    }
                    if (min_index == -1) {
                        continue;
                    }

                    Vertex& origv = origMesh.V.at(min_index);
                    boundary_neighbors.clear();
                    for (int i = 0; i < origv.N_Vids.size(); i++) {
                        if (origMesh.V.at(origv.N_Vids.at(i)).isBoundary) {
                            boundary_neighbors.push_back(origv.N_Vids.at(i));
                        }                    
                    }
                    Vertex& b1 = origMesh.V.at(boundary_neighbors.at(0));
                    Vertex& b2 = origMesh.V.at(boundary_neighbors.at(1));
                    
                    glm::dvec3 a = centers.at(v.id) - origv.xyz();
                    glm::dvec3 b = b1.xyz() - origv.xyz();
                    glm::dvec3 c = b2.xyz() - origv.xyz();

                    double length_a = glm::dot(a, a);
                    double length_b = glm::dot(b, b);
                    double length_c = glm::dot(c, c);

                    double dot_a_b = glm::dot(a, b) / glm::length(b);
                    double dot_a_c = glm::dot(a, c) / glm::length(c);
                    
                    glm::dvec3 new_coords = centers.at(v.id);
                    if (dot_a_b >= 0) {
                        b = glm::normalize(b);
                        new_coords = origv.xyz() + (dot_a_b * b);
                    } else if (dot_a_c >= 0){
                        c = glm::normalize(c);
                        new_coords = origv.xyz() + (dot_a_c * c);
                    }
                    centers.at(v.id) = new_coords;
                }
            }
            for (auto vid: smoothVids) {
                auto& v = mesh.V.at(vid);
                if (v.type == FEATURE) {
                    temp.at(v.id) = v.xyz();
                    v = centers.at(v.id);
                    bool negativeFacePresent = false;
                    for (auto nfid: v.N_Fids) {
                        if (GetScaledJacobianQuad(mesh, mesh.F.at(nfid)) < 0) {
                            negativeFacePresent = true;
                            break;
                        }
                    }
                    if (negativeFacePresent) v = temp.at(vid);
                }
            }
            newE = GetMeshEnergy();
            if (currentE - newE < 1e-4) {
                for (auto vid: smoothVids) {
                    auto& v = mesh.V.at(vid);
                    if (v.type == FEATURE) {
                        v = temp.at(v.id);
                    }
                }
                break;                
            }
        }
    }*/

    // SurfaceMapper sm(mesh, origMesh);
    // sm.SetTarget(origMesh);
	while (iters--) {
        std::vector<glm::dvec3> centers(mesh.V.size(), glm::dvec3(0.0, 0.0, 0.0));
        for (auto vid: smoothVids) {
            auto& v = mesh.V.at(vid);
            if (v.type < FEATURE) {
                centers.at(v.id) = v.xyz();
                double n = 0.0;
                glm::dvec3 center(0.0, 0.0, 0.0);
                for (int i = 0; i < v.N_Eids.size(); i++) {
                    auto& e = mesh.E.at(v.N_Eids.at(i));
                    size_t nvid = e.Vids[0] != v.id ? e.Vids[0] : e.Vids[1];
                    std::vector<size_t> nvids;
                    for (auto fid: e.N_Fids) {
                        auto& f = mesh.F.at(fid);
                        for (auto fvid: f.Vids) {
                            if (fvid != v.id && fvid != nvid && std::find(v.N_Vids.begin(), v.N_Vids.end(), fvid) != v.N_Vids.end()) {
                                nvids.push_back(fvid);
                                break;
                            }
                        }
                    }
                    auto& v_a = mesh.V.at(nvid);
                    auto& v_b = mesh.V.at(nvids.at(0));
                    auto& v_c = mesh.V.at(nvids.at(1));
                    glm::dvec3 V_j = v.xyz() - v_a.xyz();
                    glm::dvec3 V_j_minus_1 = v_b.xyz() - v_a.xyz(); 
                    glm::dvec3 V_j_plus_1 = v_c.xyz() - v_a.xyz();

                    double r = glm::length(V_j);
                    glm::dvec3 p1 = v_a.xyz() + (r * glm::normalize(V_j_minus_1));
                    glm::dvec3 p2 = v_a.xyz() + (r * glm::normalize(V_j_plus_1));
                    if (glm::length((0.5 * (p1 + p2)) - v_a.xyz()) == 0) continue;
                    glm::dvec3 direction = glm::normalize((0.5 * (p1 + p2)) - v_a.xyz());
                    glm::dvec3 newPoint_p = v_a.xyz() + (r * direction);
                    glm::dvec3 newPoint_n = v_a.xyz() + (r * (-direction));
                    
                    center += glm::distance(v.xyz(), newPoint_p) < glm::distance(v.xyz(), newPoint_n) ? newPoint_p : newPoint_n;
                    n += 1;
                }
                centers.at(vid) = (center / n);
            }
        }
        // std::vector<glm::dvec3> temp(mesh.V.size(), glm::dvec3(0.0, 0.0, 0.0));
        for (auto vid: smoothVids) {
            auto& v = mesh.V.at(vid);
            if (v.type < FEATURE) {
                // temp.at(v.id) = v.xyz();
                v = centers.at(v.id);
                // v = sm.MapPoint(centers.at(v.id));
                // v = sm.GetClosestPoint(centers.at(v.id));
                // bool negativeFacePresent = false;
                // for (auto nfid: v.N_Fids) {
                //     if (GetScaledJacobianQuad(mesh, mesh.F.at(nfid)) < 0) {
                //         negativeFacePresent = true;
                //         break;
                //     }
                // }
                // if (negativeFacePresent) v = temp.at(v.id);
            }
        }
        // sm.Map();
        /*for (auto vid: smoothVids) {
            auto& v = mesh.V.at(vid);
			if (v.type == FEATURE) {
                centers.at(v.id) = v.xyz();
                std::vector<size_t> neighbors = v.N_Vids;
                std::vector<size_t> boundary_neighbors;
                for (auto nvid: v.N_Vids) {
                    if (mesh.V.at(nvid).isBoundary) boundary_neighbors.push_back(nvid);
                }
                auto& v_b = mesh.V.at(boundary_neighbors[0]);
                auto& v_c = mesh.V.at(boundary_neighbors[1]);
                glm::dvec3 center(0, 0, 0);
                double n = 0;
                for (auto nvid: v.N_Vids) {
                    if (mesh.V.at(nvid).isBoundary) continue;
                    auto& v_a = mesh.V.at(nvid);
                    glm::dvec3 V_j = v.xyz() - v_a.xyz();
                    glm::dvec3 V_j_minus_1 = v_b.xyz() - v_a.xyz();
                    glm::dvec3 V_j_plus_1 = v_c.xyz() - v_a.xyz();

                    double r = glm::length(V_j);

                    glm::dvec3 p1 = v_a.xyz() + (r * glm::normalize(V_j_minus_1));
                    glm::dvec3 p2 = v_a.xyz() + (r * glm::normalize(V_j_plus_1));
                    if (glm::length((0.5 * (p1 + p2)) - v_a.xyz()) == 0) continue;
                    glm::dvec3 direction = glm::normalize((0.5 * (p1 + p2)) - v_a.xyz());
                    glm::dvec3 newPoint_p = v_a.xyz() + (r * direction);
                    glm::dvec3 newPoint_n = v_a.xyz() + (r * (-direction));
                    
                    center += glm::distance(v.xyz(), newPoint_p) < glm::distance(v.xyz(), newPoint_n) ? newPoint_p : newPoint_n;
                    n += 1;
                }
                if (n > 0) {
                    center /= n;
                    centers.at(v.id) = center;
                }
                double min_length = std::numeric_limits<double>::max();
                int min_index = -1;
                for (int i = 0; i < origBoundaryVids.size(); i++) {
                    double length = glm::length(centers.at(v.id) - origMesh.V.at(origBoundaryVids.at(i)).xyz());
                    if (length < min_length) {
                        min_length = length;
                        min_index = origBoundaryVids.at(i);
                    }
                }
                if (min_index == -1) {
                    continue;
                }
                Vertex& origv = origMesh.V.at(min_index);
                boundary_neighbors.clear();
                for (int i = 0; i < origv.N_Vids.size(); i++) {
                    if (origMesh.V.at(origv.N_Vids.at(i)).isBoundary) {
                        boundary_neighbors.push_back(origv.N_Vids.at(i));
                    }                    
                }
                Vertex& b1 = origMesh.V.at(boundary_neighbors.at(0));
                Vertex& b2 = origMesh.V.at(boundary_neighbors.at(1));
                
                glm::dvec3 a = centers.at(v.id) - origv.xyz();
                glm::dvec3 b = b1.xyz() - origv.xyz();
                glm::dvec3 c = b2.xyz() - origv.xyz();

                double length_a = glm::dot(a, a);
                double length_b = glm::dot(b, b);
                double length_c = glm::dot(c, c);

                double dot_a_b = glm::dot(a, b) / glm::length(b);
                double dot_a_c = glm::dot(a, c) / glm::length(c);
                
                glm::dvec3 new_coords = centers.at(v.id);
                if (dot_a_b >= 0) {
                    b = glm::normalize(b);
                    new_coords = origv.xyz() + (dot_a_b * b);
                } else if (dot_a_c >= 0){
                    c = glm::normalize(c);
                    new_coords = origv.xyz() + (dot_a_c * c);
                }
                centers.at(v.id) = new_coords;
            }
        }

        for (auto vid: smoothVids) {
            auto& v = mesh.V.at(vid);
            if (v.type == FEATURE) {
                glm::dvec3 temp = v.xyz();
                v = centers.at(v.id);
                // bool negativeFacePresent = false;
                // for (auto nfid: v.N_Fids) {
                //     if (GetScaledJacobianQuad(mesh, mesh.F.at(nfid)) < 0) {
                //         negativeFacePresent = true;
                //         break;
                //     }
                // }
                // if (negativeFacePresent) v = temp;
            }
        }*/
	}
    for (auto& v: mesh.V) v.smoothLocal = false;
    smoothVids.clear();
}

void PatchSimplifier::VertexRotate(std::multiset<SimplificationOperationStruct, bool(*)(SimplificationOperationStruct, SimplificationOperationStruct)>& SimplificationOps) {
    for (auto& v: mesh.V) {
        if (v.isBoundary) continue;
        
        SimplificationOperationStruct op;
        op.type = "Vertex_Rotate";
        double sumEdges = 0;
        double sumDiagonals = 0;
        for (auto fid: v.N_Fids) {
            Face& f = mesh.F.at(fid);
            for (auto vid: f.Vids) {
                if (vid == v.id) continue;
                if (std::find(v.N_Vids.begin(), v.N_Vids.end(), vid) == v.N_Vids.end()) {
                    if (mesh.prescribed_length > glm::length(mesh.V.at(vid).xyz() - v.xyz())) {
                        op.profitability += mesh.prescribed_length - glm::length(mesh.V.at(vid).xyz() - v.xyz());
                        op.n += 1;
                    }
                }
            }
        }
        op.n == v.N_Fids.size() ? op.profitability /= op.n : op.profitability = 0;
        for (auto eid: v.N_Eids) {
            std::vector<size_t> newVids;
            Edge& e = mesh.E.at(eid);
            newVids.push_back(e.Vids.at(0));
            Face& f1 = mesh.F.at(e.N_Fids.at(0));
            op.canceledFids.insert(f1.id);
            for (int i = 0; i < f1.Vids.size(); i++) {
                if (f1.Vids.at(i) == v.id) continue;
                if (std::find(v.N_Vids.begin(), v.N_Vids.end(), f1.Vids.at(i)) == v.N_Vids.end()) {
                    newVids.push_back(f1.Vids.at(i));
                    break;
                }
            }
            newVids.push_back(e.Vids.at(1));
            Face& f2 = mesh.F.at(e.N_Fids.at(1));
            op.canceledFids.insert(f2.id);
            for (int i = 0; i < f2.Vids.size(); i++) {
                if (f2.Vids.at(i) == v.id) continue;
                if (std::find(v.N_Vids.begin(), v.N_Vids.end(), f2.Vids.at(i)) == v.N_Vids.end()) {
                    newVids.push_back(f2.Vids.at(i));
                    break;
                }
            }
            Face newF;
            newF.Vids = newVids;
        }
        if (op.profitability > 0) {
            SimplificationOps.insert(op);
        }
    }
}

void PatchSimplifier::EdgeRotate(std::multiset<SimplificationOperationStruct, bool(*)(SimplificationOperationStruct, SimplificationOperationStruct)>& SimplificationOps) {
    for (auto& e: mesh.E) {
        bool isBoundary = false;
        if (e.isBoundary || mesh.V.at(e.Vids.at(0)).isBoundary || mesh.V.at(e.Vids.at(1)).isBoundary) isBoundary = true;
        for (auto vid: mesh.V.at(e.Vids[0]).N_Vids) {
            if (mesh.V.at(vid).isBoundary) {
                isBoundary = true;
                break;
            }
        }
        for (auto vid: mesh.V.at(e.Vids[1]).N_Vids) {
            if (mesh.V.at(vid).isBoundary) {
                isBoundary = true;
                break;
            }
        }
        if (isBoundary) continue;
        SimplificationOperationStruct op;
        op.type = "Edge_Rotate";

        Vertex& v1 = mesh.V.at(e.Vids.at(0));
        Vertex& v2 = mesh.V.at(e.Vids.at(1));
        Face& f1 = mesh.F.at(e.N_Fids.at(0));
        Face& f2 = mesh.F.at(e.N_Fids.at(1));

        std::vector<size_t> v1_nVids;
        std::vector<size_t> v2_nVids;
        for (auto vid: f1.Vids) {
            if (vid == v1.id || vid == v2.id) continue;
            if (std::find(v1.N_Vids.begin(), v1.N_Vids.end(), vid) != v1.N_Vids.end()) {
                v1_nVids.push_back(vid);
            }
            if (std::find(v2.N_Vids.begin(), v2.N_Vids.end(), vid) != v2.N_Vids.end()) {
                v2_nVids.push_back(vid);
            }
        }
        for (auto vid: f2.Vids) {
            if (vid == v1.id || vid == v2.id) continue;
            if (std::find(v1.N_Vids.begin(), v1.N_Vids.end(), vid) != v1.N_Vids.end()) {
                v1_nVids.push_back(vid);
            }
            if (std::find(v2.N_Vids.begin(), v2.N_Vids.end(), vid) != v2.N_Vids.end()) {
                v2_nVids.push_back(vid);
            }
        }

        double current_edge_length = glm::length(v1.xyz() - v2.xyz());

        Face new_f1;
        Face new_f2;
        Vertex& new_v1 = mesh.V.at(v1_nVids.at(0));
        new_f1.Vids.push_back(new_v1.id);
        new_f2.Vids.push_back(new_v1.id);
        double new_edge_length1 = 0;
        double diag1 = 0;
        double diag2 = 0;
        double new_diag1 = glm::length(mesh.V.at(v1_nVids.at(0)).xyz() - mesh.V.at(v1_nVids.at(1)).xyz());
        double new_diag2 = glm::length(mesh.V.at(v2_nVids.at(0)).xyz() - mesh.V.at(v2_nVids.at(1)).xyz());
        if (std::find(new_v1.N_Vids.begin(), new_v1.N_Vids.end(), v2_nVids.at(0)) == new_v1.N_Vids.end()) {
            new_f1.Vids.push_back(v2_nVids.at(1));
            new_f1.Vids.push_back(v2.id);
            new_f1.Vids.push_back(v2_nVids.at(0));

            new_f2.Vids.push_back(v2_nVids.at(0));
            new_f2.Vids.push_back(v1_nVids.at(1));
            new_f2.Vids.push_back(v1.id);

            new_edge_length1 = glm::length(mesh.V.at(v1_nVids.at(0)).xyz() - mesh.V.at(v2_nVids.at(1)).xyz());
            diag1 = glm::length(mesh.V.at(v1_nVids.at(1)).xyz() - v2.xyz());
            diag2 = glm::length(mesh.V.at(v2_nVids.at(0)).xyz() - v1.xyz());
        } else {
            new_f1.Vids.push_back(v2_nVids.at(0));
            new_f1.Vids.push_back(v2.id);
            new_f1.Vids.push_back(v2_nVids.at(1));
            
            new_f2.Vids.push_back(v2_nVids.at(1));
            new_f2.Vids.push_back(v1_nVids.at(1));
            new_f2.Vids.push_back(v1.id);
            
            new_edge_length1 = glm::length(mesh.V.at(v1_nVids.at(0)).xyz() - mesh.V.at(v2_nVids.at(0)).xyz());
            diag1 = glm::length(mesh.V.at(v1_nVids.at(1)).xyz() - v2.xyz());
            diag2 = glm::length(mesh.V.at(v2_nVids.at(1)).xyz() - v1.xyz());
        }
        double rotation1_profitability = -1;
        // if (mesh.prescribed_length > new_edge_length1 && mesh.prescribed_length > new_diag1 && mesh.prescribed_length > new_diag2) {
        //     rotation1_profitability = (mesh.prescribed_length - new_edge_length1) + (mesh.prescribed_length - new_diag1) + (mesh.prescribed_length - new_diag2);
        //     rotation1_profitability /= 3;
        // }
        if (current_edge_length > new_edge_length1 && diag1 > new_diag1 && diag2 > new_diag2) {
            rotation1_profitability = (current_edge_length - new_edge_length1) + (diag1 - new_diag1) + (diag2 - new_diag2);
        }
        Face new_f3;
        Face new_f4;
        Vertex& new_v2 = mesh.V.at(v1_nVids.at(1));
        new_f3.Vids.push_back(new_v2.id);
        new_f4.Vids.push_back(new_v2.id);
        double new_edge_length2 = 0;
        double diag3 = 0;
        double diag4 = 0;
        if (std::find(new_v2.N_Vids.begin(), new_v2.N_Vids.end(), v2_nVids.at(0)) == new_v2.N_Vids.end()) {
            new_f3.Vids.push_back(v2_nVids.at(1));
            new_f3.Vids.push_back(v2.id);
            new_f3.Vids.push_back(v2_nVids.at(0));

            new_f4.Vids.push_back(v2_nVids.at(0));
            new_f4.Vids.push_back(v1_nVids.at(0));
            new_f4.Vids.push_back(v1.id);
                    
            new_edge_length2 = glm::length(mesh.V.at(v1_nVids.at(1)).xyz() - mesh.V.at(v2_nVids.at(1)).xyz());
            diag3 = glm::length(mesh.V.at(v1_nVids.at(0)).xyz() - v2.xyz());
            diag4 = glm::length(mesh.V.at(v2_nVids.at(0)).xyz() - v1.xyz());
        } else {
            new_f3.Vids.push_back(v2_nVids.at(0));
            new_f3.Vids.push_back(v2.id);
            new_f3.Vids.push_back(v2_nVids.at(1));
            
            new_f4.Vids.push_back(v2_nVids.at(1));
            new_f4.Vids.push_back(v1_nVids.at(0));
            new_f4.Vids.push_back(v1.id);
                    
            new_edge_length2 = glm::length(mesh.V.at(v1_nVids.at(1)).xyz() - mesh.V.at(v2_nVids.at(0)).xyz());
            diag3 = glm::length(mesh.V.at(v1_nVids.at(0)).xyz() - v2.xyz());
            diag4 = glm::length(mesh.V.at(v2_nVids.at(1)).xyz() - v1.xyz());
        }
        double rotation2_profitability = -1;
        // if (mesh.prescribed_length > new_edge_length2 && mesh.prescribed_length > new_diag1 && mesh.prescribed_length > new_diag2) {
        //     rotation2_profitability = (mesh.prescribed_length - new_edge_length2) + (mesh.prescribed_length - new_diag1) + (mesh.prescribed_length - new_diag2);
        //     rotation2_profitability /= 3;
        // }
        if (current_edge_length > new_edge_length2 && diag3 > new_diag1 && diag4 > new_diag2) {
            rotation2_profitability = (current_edge_length - new_edge_length2) + (diag3 - new_diag1) + (diag4 - new_diag2);
        }

        if (rotation1_profitability > rotation2_profitability) {
            op.newFaces.push_back(new_f1);
            op.newFaces.push_back(new_f2);
            op.profitability = rotation1_profitability;
        } else {
            op.newFaces.push_back(new_f3);
            op.newFaces.push_back(new_f4);
            op.profitability = rotation2_profitability;
        }
        op.canceledFids.insert(e.N_Fids.begin(), e.N_Fids.end());
        if (op.profitability > 0) { 
            SimplificationOps.insert(op);
        }
    }
}

void PatchSimplifier::EdgeCollapse(std::multiset<SimplificationOperationStruct, bool(*)(SimplificationOperationStruct, SimplificationOperationStruct)>& SimplificationOps) {
    for (auto& e: mesh.E) {
        bool isBoundary = false;
        if (e.isBoundary || mesh.V.at(e.Vids.at(0)).isBoundary || mesh.V.at(e.Vids.at(1)).isBoundary) isBoundary = true;
        for (auto vid: e.N_Vids) {
            if (mesh.V.at(vid).isBoundary) {
                isBoundary = true;
                break;
            }
        }
        if (isBoundary) continue;
        SimplificationOperationStruct op;
        op.type = "Edge_Collapse";
        op.profitability = mesh.prescribed_length - glm::length(mesh.V.at(e.Vids[0]).xyz() - mesh.V.at(e.Vids[1]).xyz());
        size_t source = e.Vids.at(0);
        size_t target = e.Vids.at(1);
        if (mesh.V.at(target).isBoundary) {
            source = e.Vids.at(1);
            target = e.Vids.at(0);
        }
        Vertex& v1 = mesh.V.at(source);
        Vertex& v2 = mesh.V.at(target);
        op.canceledFids.insert(v1.N_Fids.begin(), v1.N_Fids.end());
        for (auto eid: v1.N_Eids) {
            if (eid == e.id) continue;
            std::vector<size_t> newVids;
            Edge& e_n = mesh.E.at(eid);
            newVids.push_back(e_n.Vids.at(0));
            Face& f1 = mesh.F.at(e_n.N_Fids.at(0));
            for (int i = 0; i < f1.Vids.size(); i++) {
                if (f1.Vids.at(i) == v1.id) continue;
                if (std::find(v1.N_Vids.begin(), v1.N_Vids.end(), f1.Vids.at(i)) == v1.N_Vids.end()) {
                    newVids.push_back(f1.Vids.at(i));
                    break;
                }
            }
            newVids.push_back(e_n.Vids.at(1));
            Face& f2 = mesh.F.at(e_n.N_Fids.at(1));
            for (int i = 0; i < f2.Vids.size(); i++) {
                if (f2.Vids.at(i) == v1.id) continue;
                if (std::find(v1.N_Vids.begin(), v1.N_Vids.end(), f2.Vids.at(i)) == v1.N_Vids.end()) {
                    newVids.push_back(f2.Vids.at(i));
                    break;
                }
            }
            Face newF;
            newF.Vids = newVids;
            op.newFaces.push_back(newF);
        }
        for (auto& f: op.newFaces) {
            for (int i = 0; i < f.Vids.size(); i++) {
                if (f.Vids.at(i) == v1.id) {
                    f.Vids.at(i) = v2.id;
                }
            }
        }
        if (op.profitability > 0) {
            SimplificationOps.insert(op);
        }
    }
}

void PatchSimplifier::DiagonalCollapse(std::multiset<SimplificationOperationStruct, bool(*)(SimplificationOperationStruct, SimplificationOperationStruct)>& SimplificationOps) {
    for (auto& f: mesh.F) {
        bool isBoundary = false;
        for (auto vid: f.Vids) {
            if (mesh.V.at(vid).isBoundary) {
                isBoundary = true;
                break;
            }
        }
        if (isBoundary) continue;
        SimplificationOperationStruct op;
        op.type = "Diagonal_Collapse";
        double l1 = glm::length(mesh.V.at(f.Vids[0]).xyz() - mesh.V.at(f.Vids[2]).xyz());
        double l2 = glm::length(mesh.V.at(f.Vids[1]).xyz() - mesh.V.at(f.Vids[3]).xyz());
        size_t target_id = 0;
        size_t source_id = 0;
        if (l1 < l2) {
            op.profitability = l1 / (sqrt(2) * mesh.prescribed_length);
            target_id = f.Vids[0];
            source_id = f.Vids[2];
        } else {
            op.profitability = l2 / (sqrt(2) * mesh.prescribed_length);
            
            target_id = f.Vids[1];
            source_id = f.Vids[3];
        }
        Vertex& target = mesh.V.at(target_id);
        Vertex& source = mesh.V.at(source_id);
        // op.updateVertexIds.push_back(target_id);
        // op.updatedVertexPos.push_back(0.5 * (source.xyz() + target.xyz()));
        for (auto fid: source.N_Fids) {
            if (std::find(target.N_Fids.begin(), target.N_Fids.end(), fid) == target.N_Fids.end()) {
                Face& n_f = mesh.F.at(fid);
                Face newF;
                for (int i = 0; i < n_f.Vids.size(); i++) {
                    if (n_f.Vids.at(i) == source.id) {
                        newF.Vids.push_back(target.id);
                    } else {
                        newF.Vids.push_back(n_f.Vids.at(i));
                    }
                }
                op.newFaces.push_back(newF);
                op.canceledFids.insert(n_f.id);
            } else {
                op.canceledFids.insert(fid);
            }
        }
        if (op.profitability > 0 && op.profitability < 1) {
            SimplificationOps.insert(op);
        }
    }
}

// TODO:
// Writing section 7
// Edit current images in results section
// Generating results

