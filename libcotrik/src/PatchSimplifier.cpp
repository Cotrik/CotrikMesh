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

bool OpSortAscending(SimplificationOperation op1, SimplificationOperation op2) {
    return op1.profitability < op2.profitability;
}

bool OpSortDescending(SimplificationOperation op1, SimplificationOperation op2) {
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
            // if (!Simplify(iter)) break;
            if (!SimplifyMesh(iter)) break;
        }
        ++Simplifier::maxValence;
    }
}

bool PatchSimplifier::SimplifyMesh(int& iter) {
    std::set<size_t> canceledFids;
    init();

    if (iter == 0) {
        mesh.GetQuadMeshArea();
        featurePreserved ? get_feature() : origMesh = mesh;
        for (auto& v: origMesh.V) {
            if (v.isBoundary && !v.isCorner) {
                origBoundaryVids.push_back(v.id);
            }
        }
    }

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
    
    if (canceledFids.empty()) {
        bool(*fn_pt1)(SimplificationOperation, SimplificationOperation) = OpSortDescending;
        bool(*fn_pt2)(SimplificationOperation, SimplificationOperation) = OpSortAscending;
        std::multiset<SimplificationOperation, bool(*)(SimplificationOperation, SimplificationOperation)> SimplificationOps(fn_pt1);
        GetOperations(SimplificationOps);

        if (!SimplificationOps.empty()) {
            std::cout << "# of Simplification Operations: " << SimplificationOps.size() << std::endl;
            for (auto op: SimplificationOps) {
                std::cout << op.type << ": " << op.profitability << std::endl;
            }
        }
    }
    
    ////////////////////////////////////////////////////////////////////////////

    if (canceledFids.empty()) {
        update(canceledFids);
        return false;
    }
    
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

void PatchSimplifier::GetOperations(std::multiset<SimplificationOperation, bool(*)(SimplificationOperation, SimplificationOperation)>& SimplificationOps) {
    BaseComplexQuad baseComplex(mesh);
    baseComplex.ExtractSingularVandE();
    baseComplex.BuildE();
    
    
    if (Simplifier::COLLAPSE_DIAGNAL) {
        DiagnalCollapseSimplifier diagonalCollapseSimplifier(mesh);
        diagonalCollapseSimplifier.GetDiagonalCollapseOps(SimplificationOps);
    }

    if (Simplifier::COLLAPSE) {
        GetSeparatrixCollapseOps(baseComplex, false, SimplificationOps);
        GetSeparatrixCollapseOps(baseComplex, true, SimplificationOps);
    }

    if (Simplifier::HALF) {
        GetHalfSeparatrixOps(baseComplex, SimplificationOps);
    }

    
}

bool PatchSimplifier::Simplify(int& iter) {
    std::set<size_t> canceledFids;
    init();
    
    if (iter == 0 && featurePreserved) {
        get_feature();
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
        SmoothMesh();
    }
    else {
        smoothGlobal = true;
        SmoothMesh();  
    }
    // if (mesh.F.size() <= 0.25 * refinementFactor * origMesh.F.size()) {
    //     std::cout << "current Faces: " << mesh.F.size() << " original Faces: " << origMesh.F.size() << std::endl;
        
    //     smoothGlobal = true;
    //     SmoothMesh();
    //     RefineMesh();
    //     smoothGlobal = true;
    //     SmoothMesh();
    // }
    
    // Step 1 -- doublet removal
    if (canceledFids.empty() && Simplifier::REMOVE_DOUBLET) {
        DoubletSimplifier doubletSimplifier(mesh);
        // doubletSimplifier.Run(canceledFids);
        doubletSimplifier.RunCollective(canceledFids);
        if (!canceledFids.empty()) std::cout << "remove_doublet" << std::endl;
    }
    // Step 2 -- doublet splitting
    if (canceledFids.empty() && Simplifier::SHEET_SPLIT) {
        SheetSplitSimplifier sheetSplitSimplifier(mesh);
        sheetSplitSimplifier.Run(canceledFids);
        if (!canceledFids.empty()) std::cout << "remove_doublet from sheetSplitSimplifier" << std::endl;
    }
    // Step 4 -- singlet collapsing
    if (canceledFids.empty()) {
        DiagnalCollapseSimplifier diagnalCollapseSimplifier(mesh);
        // diagnalCollapseSimplifier.Run3(canceledFids);
        diagnalCollapseSimplifier.CollapseSinglets(canceledFids);
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

    // Step 8 -- diagonal collapsing
    if (canceledFids.empty() && Simplifier::COLLAPSE_DIAGNAL) {
        DiagnalCollapseSimplifier diagnalCollapseSimplifier(mesh);
        // diagnalCollapseSimplifier.Run(canceledFids);
        diagnalCollapseSimplifier.RunCollective(canceledFids);
        if (!canceledFids.empty()) std::cout << "collapse_diagnal" << std::endl;
    }
    
    // static bool aligned = false;
    // if (canceledFids.empty() && !aligned) {
    //     aligned = true;
    //     std::cout << "writing rotate.vtk " << std::endl;
    //     MeshFileWriter writer(mesh, "rotate.vtk");
    //     writer.WriteFile();
    // }
    
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

    if (canceledFids.empty() && Simplifier::HALF) {
       BaseComplexQuad baseComplex(mesh);
       baseComplex.ExtractSingularVandE();
       baseComplex.BuildE();
    //    half_simplify(baseComplex, canceledFids);
        half_separatrix_collapse(baseComplex, canceledFids);

        if (!canceledFids.empty()) std::cout << "half_simplify\n";
   }

   if (canceledFids.empty() && Simplifier::GLOBAL) {
       SheetSimplifier sheetSimplifier(mesh);
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

bool PatchSimplifier::IsFaceNegative(int fid, int vid, glm::dvec3 false_coord) {
    Face& f = mesh.F.at(fid);
    int sign = 0;
    glm::dvec3 temp_coord(0.0, 0.0, 0.0);
    for (auto id: f.Vids) {
        if (id == vid) {
            temp_coord = mesh.V.at(id).xyz();
            mesh.V.at(id) = false_coord;
        }
    }    
    for (int i = 0; i < f.Vids.size(); i++) {
        Vertex& a = mesh.V.at(f.Vids.at(i));
        Vertex& b = mesh.V.at(f.Vids.at((i + 1) % f.Vids.size()));
        Vertex& c = mesh.V.at(f.Vids.at((i + 2) % f.Vids.size()));
        double det = ((b.x - a.x) * (c.y - a.y)) - ((c.x - a.x) * (b.y - a.y));
        if (det > 0) {
            sign += 1;
        } else if (det < 0) {
            sign -= 1;
        }
    }
    for (auto id: f.Vids) {
        if (id == vid) {
            mesh.V.at(id) = temp_coord;
        }
    }
    if (abs(sign) == 4 || sign == 0) {
        return false;
    }
    return true;
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
