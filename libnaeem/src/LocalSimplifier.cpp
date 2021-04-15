/*
* Simplifier.h
*
*  Created on: October 7, 2020
*      Author: https://github.com/naeem014
*/

// #include <bits/stdc++.h>

#include "LocalSimplifier.h"
#include "DoubletSimplifier.h"
#include "DiagnalCollapseSimplifier.h"
#include "AngleBasedSmoothQuadMesh.h"

#include <ctime>

bool OpSortAscending(LocalOperation op1, LocalOperation op2) {
    return op1.profitability < op2.profitability;
}

bool OpSortDescending(LocalOperation op1, LocalOperation op2) {
    return op1.profitability > op2.profitability;
}

LocalSimplifier::LocalSimplifier(Mesh& mesh) : Simplifier(mesh) {}

LocalSimplifier::~LocalSimplifier() {}

void LocalSimplifier::Simplify() {
    SmoothAlgorithm smooth_algo(mesh, mesh, 1000, 1, true, false);
    init();
    // get_feature();
    std::vector<int> valence_histogram(10, 0);
    for (auto& v: mesh.V) {
        if (v.isBoundary) continue;
        valence_histogram.at(v.N_Vids.size()) += 1;
    }
    for (int i = 0; i < valence_histogram.size(); i++) {
        if (valence_histogram.at(i) > 0) {
            std::cout << i << " " << valence_histogram.at(i) << std::endl;
        }
    }
    return;
    // RemoveBoundaryDegenerates();
    // return;
    // std::ofstream ofs("feature_vertices.vtk");
	// ofs << "# vtk DataFile Version 3.0\n"
	// 	<< "feature_vertices.vtk\n"
	// 	<< "ASCII\n\n"
	// 	<< "DATASET UNSTRUCTURED_GRID\n";
	// ofs << "POINTS " << mesh.V.size() << " double\n";

    // std::vector<size_t> feature_edges;
    // for (auto& e: mesh.E) {
    //     if (e.isBoundary) {
    //         feature_edges.push_back(e.id);
    //     } 
    // }

    // for (auto& v: mesh.V) {
    //     ofs << v.x << " " << v.y << " " << v.z << std::endl;
    // }

	// ofs << "CELLS " << feature_edges.size() << " " << feature_edges.size() * 3 << std::endl;
    // for (auto eid: feature_edges) {
    //     Edge& e = mesh.E.at(eid);
    //     ofs << "2 " << e.Vids.at(0) << " " << e.Vids.at(1) << std::endl; 
    // }
	// ofs << "CELL_TYPES " << feature_edges.size() << "\n";
    // for (auto eid: feature_edges) {
    //     ofs << "3 " << std::endl; 
    // }
    
    // return;
    int it = 0;
    mesh_prescribed_length = getPrescribedLength();
    getLocalPrescribedLength();
    std::cout << "Prescribed length: " << mesh_prescribed_length << std::endl;
    // std::set<size_t> Vids;
    // for (auto& v: mesh.V) {
    //     // if (v.isBoundary) continue;
    //     Vids.insert(v.id);
    // }
    smoothMesh(1000, true);
    // Vids.clear();
    // smoothMeshGlobal();
    // return;
    std::clock_t start;
    double duration;
    start = std::clock();
    while (it < iters) {
        std::cout << "it: " << it << std::endl;
        std::cout << "V: " << mesh.V.size() << std::endl;
        std::cout << "E: " << mesh.E.size() << std::endl;
        std::cout << "F: " << mesh.F.size() << std::endl;
        std::set<size_t> canceledFids;
        // mesh_prescribed_length = getPrescribedLength();
        // getLocalPrescribedLength();
        // double current_energy = getMeshEnergy();
        // bool doublets_removed = false;
        // while (!doublets_removed) {
        //     RemoveDoublets(canceledFids);
        //     if (!canceledFids.empty()) {
        //         update(canceledFids);
        //         init();
        //     } else {
        //         doublets_removed = true;
        //     }
        // }
        // break;
        bool(*fn_pt1)(LocalOperation, LocalOperation) = OpSortDescending;
        bool(*fn_pt2)(LocalOperation, LocalOperation) = OpSortAscending;
        std::multiset<LocalOperation, bool(*)(LocalOperation, LocalOperation)> OptimizationOps(fn_pt1);
        std::multiset<LocalOperation, bool(*)(LocalOperation, LocalOperation)> CollapseOps(fn_pt2);
        getOperationsPriorities(OptimizationOps, CollapseOps);
        // std::cout << "Optimization Ops: " << CollapseOps.size() << std::endl;
        // for (auto op: CollapseOps) {
        //     std::cout << "(" << op.type << " " << op.profitability << ") ";
        // }
        // std::cout << std::endl;
        // return;
        std::set<size_t> processedFids;
        for(auto op: OptimizationOps) {
            if (op.profitability < 0) continue;
            bool overlappedOp = false;
            for (auto id: op.canceledFids) {
                if (std::find(processedFids.begin(), processedFids.end(), id) != processedFids.end()) {
                    overlappedOp = true;
                    break;
                }
            }
            if (overlappedOp) continue;
            std::cout << "it: " << it << " Optimization Op: " << op.type << " profitability: " << op.profitability << std::endl;
            for (auto f: op.newFaces) {
                f.id = mesh.F.size();
                // Vids.insert(f.Vids.begin(), f.Vids.end());
                // for (auto id: f.Vids) {
                //     Vids.insert(mesh.V.at(id).N_Vids.begin(), mesh.V.at(id).N_Vids.end());
                // }
                for (auto vid: f.Vids) {
                    Vertex& f_v = mesh.V.at(vid);
                    f_v.smoothLocal = true;
                    for (auto n_vid: f_v.N_Vids) {
                        mesh.V.at(n_vid).smoothLocal = true;
                    }
                }
                mesh.F.push_back(f);
            }
            canceledFids.insert(op.canceledFids.begin(), op.canceledFids.end());
            processedFids.insert(op.canceledFids.begin(), op.canceledFids.end());
            for (auto fid: op.canceledFids) {
                Face& f = mesh.F.at(fid);
                processedFids.insert(f.N_Fids.begin(), f.N_Fids.end());
            }
        }
        // if (!canceledFids.empty()) {
        //     update(canceledFids);
        //     init();
        // }
        // bool doublets_removed = false;
        // while (!doublets_removed) {
        //     RemoveDoublets(canceledFids);
        //     if (!canceledFids.empty()) {
        //         update(canceledFids);
        //         init();
        //     } else {
        //         doublets_removed = true;
        //     }
        // }
        // OptimizationOps.clear();
        // CollapseOps.clear();
        // processedFids.clear();
        // if (it % 10 == 0) {
            // mesh_prescribed_length = getPrescribedLength();
            // smoothMesh(Vids, 1000);
            // Vids.clear();
        // }
        // smoothMesh(Vids, 1000);
        // Vids.clear();
        // std::cout << Vids.size() << std::endl;
        // it++;
        // continue;
        // getOperationsPriorities(OptimizationOps, CollapseOps);
        for (auto op: CollapseOps) {
            // if (op.type.compare("Diagonal Collapse") == 0) continue;
            bool overlappedOp = false;
            for (auto id: op.canceledFids) {
                if (std::find(processedFids.begin(), processedFids.end(), id) != processedFids.end()) {
                    overlappedOp = true;
                    break;
                }
            }
            if (overlappedOp) continue;
            std::cout << "it: " << it << " Collapse Op: " << op.type << " profitability: " << op.profitability <<  std::endl;
            for (auto f: op.newFaces) {
                f.id = mesh.F.size();
                // Vids.insert(f.Vids.begin(), f.Vids.end());
                // for (auto id: f.Vids) {
                //     Vids.insert(mesh.V.at(id).N_Vids.begin(), mesh.V.at(id).N_Vids.end());
                // }
                for (auto vid: f.Vids) {
                    Vertex& f_v = mesh.V.at(vid);
                    f_v.smoothLocal = true;
                    for (auto n_vid: f_v.N_Vids) {
                        mesh.V.at(n_vid).smoothLocal = true;
                    }
                }
                mesh.F.push_back(f);
            }
            canceledFids.insert(op.canceledFids.begin(), op.canceledFids.end());
            processedFids.insert(op.canceledFids.begin(), op.canceledFids.end());
            for (auto fid: op.canceledFids) {
                Face& f = mesh.F.at(fid);
                processedFids.insert(f.N_Fids.begin(), f.N_Fids.end());
            }
            // break;
        }
        if (!canceledFids.empty()) {
            update(canceledFids);
            init();
        } else {
            break;
        }
        bool doublets_removed = false;
        while (!doublets_removed) {
            RemoveDoublets(canceledFids);
            if (!canceledFids.empty()) {
                update(canceledFids);
                init();
            } else {
                doublets_removed = true;
            }
        }
        // if (Vids.empty()) {
        //     break;
        // }
        // if (it % 10 == 0) {
            // mesh_prescribed_length = getPrescribedLength();
            
            smoothMesh(1000, false);
            // Vids.clear();
        // }
        // smoothMesh(Vids, 20);
        // mesh_prescribed_length = getPrescribedLength();
        // getLocalPrescribedLength();
        // double new_energy = getMeshEnergy();
        // if (current_energy - new_energy < 0) {
        //     break;
        // }
        // smoothMeshGlobal();
        it++;
    }
    duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
    std::cout << "Simplification time: " << duration << " seconds" << std::endl;
    // Vids.clear();
    // RemoveBoundaryDegenerates();
    // for (auto& v: mesh.V) {
    //     Vids.insert(v.id);
    // }
    smoothMesh(1000, true);
    // smoothMeshGlobal();
}

double LocalSimplifier::getPrescribedLength() {
    double p_length = 0;
    double total_area = 0;
    int n_faces = mesh.F.size();
    // double min_x = mesh.V.at(0).x;
    // double min_y = mesh.V.at(0).y;
    // double max_x = mesh.V.at(0).x;
    // double max_y = mesh.V.at(0).y;
    // for (auto& v: mesh.V) {
    //     if (v.x < min_x) min_x = v.x;
    //     if (v.y < min_y) min_y = v.y;
    //     if (v.x > max_x) max_x = v.x;
    //     if (v.y > max_y) max_y = v.y;
    // }
    // double x1 = (max_x - min_x);
    // double x2 = (max_y - min_y);
    // double X = x1 > x2 ? x1 : x2;
    // total_area = x1 * x2;

    // std::cout << "min_x: " << min_x << " min_y: " << min_y << std::endl;
    // std::cout << "max_x: " << max_x << " max_y: " << max_y << std::endl;
    // std::cout << "x1: " << x1 << " x2: " << x2 << std::endl;
    // std::cout << "X: " << X << std::endl;
    // std::cout << "total_area: " << total_area << std::endl;
    // std::cout << "# faces: " << n_faces << std::endl;
    // std::ofstream ofs("ideal_square.vtk");
	// ofs << "# vtk DataFile Version 3.0\n"
	// 	<< "output.vtk\n"
	// 	<< "ASCII\n\n"
	// 	<< "DATASET UNSTRUCTURED_GRID\n";
	// ofs << "POINTS " << 4 << " double\n";
	
    // ofs << min_x << " " << min_y  << " " <<  0 << "\n";
    // ofs << max_x << " " <<  min_y << " " <<  0 << "\n";
    // ofs << max_x << " " << min_y + (max_x - min_x)  << " " <<  0 << "\n";
    // ofs << min_x << " " << min_y + (max_x - min_x)  << " " <<  0 << "\n";
	// ofs << "CELLS " << 1 << " " << 5 << std::endl;
    // ofs << "4 0 1 2 3" << std::endl;
	// ofs << "CELL_TYPES " << 1 << "\n";
    // ofs << "9" << std::endl;
    for (auto& f: mesh.F) {
        total_area += getArea(f.Vids);
    }
    p_length = sqrt(total_area / n_faces);

    // p_agg = 0;
    // lambda_p.clear();
    // lambda_p.resize(mesh.V.size(), 0);
    // for (int i = 0; i < mesh.V.size(); i++) {
    //     Vertex& v = mesh.V.at(i);
    //     double local_area = 0;
    //     for (auto fid: v.N_Fids) {
    //         local_area += getArea(mesh.F.at(fid).Vids);
    //     }
    //     double local_prescribed_length = sqrt(local_area / v.N_Vids.size());
    //     // double avg_length = 0;
    //     // for (auto vid: v.N_Vids) {
    //     //     double e_length = getLength(v.id, vid);
    //     //     avg_length += e_length;
    //     // }
    //     // avg_length = avg_length / v.N_Vids.size();
    //     // double p_ratio = avg_length / p_length;
    //     // lambda_p.at(i) = p_ratio;
    //     lambda_p.at(i) = local_prescribed_length;
    //     // lambda_p.at(i) = 1;
    // }
    return p_length;
}

void LocalSimplifier::getLocalPrescribedLength() {
    p_agg = 0;
    lambda_p.clear();
    lambda_p.resize(mesh.V.size(), 0);
    double factor = 1;
    for (int i = 0; i < mesh.V.size(); i++) {
        Vertex& v = mesh.V.at(i);
        double local_area = 0;
        for (auto fid: v.N_Fids) {
            local_area += getArea(mesh.F.at(fid).Vids);
        }
        double local_prescribed_length = sqrt(local_area / (v.N_Vids.size() * factor));
        // double local_prescribed_length = sqrt(local_area / 4);
        // double local_prescribed_length = sqrt(local_area) / 4;
        // double avg_length = 0;
        // for (auto vid: v.N_Vids) {
        //     double e_length = getLength(v.id, vid);
        //     avg_length += e_length;
        // }
        // avg_length = avg_length / v.N_Vids.size();
        // double p_ratio = avg_length / p_length;
        // lambda_p.at(i) = p_ratio;
        lambda_p.at(i) = local_prescribed_length;
        v.prescribed_length = local_prescribed_length;
        // lambda_p.at(i) = 1;
    }
}

double LocalSimplifier::getMeshEnergy() {
    double energy = 0;

    for (auto& e: mesh.E) {
        double e1 = getLength(e.Vids.at(0), e.Vids.at(1)) - mesh_prescribed_length;
        energy += (e1 * e1);
    }

    for (auto& f: mesh.F) {
        double e1 = getLength(f.Vids.at(0), f.Vids.at(2)) - (sqrt(2) * mesh_prescribed_length);
        double e2 = getLength(f.Vids.at(1), f.Vids.at(3)) - (sqrt(2) * mesh_prescribed_length);
        energy += (e1 * e1);
        energy += (e2 * e2);
    }

    return energy;
}

void LocalSimplifier::getOperationsPriorities(std::multiset<LocalOperation, bool(*)(LocalOperation, LocalOperation)>& OptimizationOps, std::multiset<LocalOperation, bool(*)(LocalOperation, LocalOperation)>& CollapseOps) {
    for (auto& v: mesh.V) {
        v.smoothLocal = false;
        // if (v.isBoundary) continue;
        // bool hasBoundaryNeighbor = false;
        // for (auto vid: v.N_Vids) {
        //     if (mesh.V.at(vid).isBoundary) {
        //         hasBoundaryNeighbor = true;
        //     }
        // }
        // if (hasBoundaryNeighbor) continue;
        // if (v.N_Vids.size() == 4) continue;
        VertexRotate(v, OptimizationOps);
        for (auto eid: v.N_Eids) {
            Edge& e = mesh.E.at(eid);
            // bool hasBoundaryVertexNeighbor = false;
            // for (auto vid: e.Vids) {
            //     if (mesh.V.at(vid).isBoundary) {
            //         hasBoundaryVertexNeighbor = true;
            //         break;
            //     }
            // }
            // for (auto vid: e.N_Vids) {
            //     if (mesh.V.at(vid).isBoundary) {
            //         hasBoundaryVertexNeighbor = true;
            //         break;
            //     }
            // }
            // if (hasBoundaryVertexNeighbor) continue;
            EdgeCollapse(v, e, CollapseOps);
            EdgeRotate(v, e, OptimizationOps);
        }

        for (auto fid: v.N_Fids) {
            Face& f = mesh.F.at(fid);
            // bool hasBoundaryVertexNeighbor = false;
            // bool hasRegularVertex = true;
            std::vector<size_t> diag;
            for (auto vid: f.Vids) {
                // if (mesh.V.at(vid).isBoundary) {
                //     hasBoundaryVertexNeighbor = true;
                // }
                // if (mesh.V.at(vid).N_Vids.size() != 4) {
                //     hasRegularVertex = false;
                // }
                if (vid != v.id && std::find(v.N_Vids.begin(), v.N_Vids.end(), vid) == v.N_Vids.end()) {
                    diag.push_back(v.id);
                    diag.push_back(vid);
                }
            }
            // for (auto fid: f.N_Fids) {
            //     for (auto vid: mesh.F.at(fid).Vids)  {       
            //         if (mesh.V.at(vid).isBoundary) {
            //             hasBoundaryVertexNeighbor = true;
            //         }
            //     }
            // }
            // if (hasBoundaryVertexNeighbor) continue;
            DiagonalCollapse(v, diag, CollapseOps);
        }
    }
    
    // Edge Rotate and Edge Collapse
    // for (auto& e: mesh.E) {
    //     bool hasBoundaryVertexNeighbor = false;
    //     for (auto vid: e.Vids) {
    //         if (mesh.V.at(vid).isBoundary) {
    //             hasBoundaryVertexNeighbor = true;
    //             break;
    //         }
    //     }
    //     for (auto vid: e.N_Vids) {
    //         if (mesh.V.at(vid).isBoundary) {
    //             hasBoundaryVertexNeighbor = true;
    //             break;
    //         }
    //     }
    //     if (hasBoundaryVertexNeighbor) continue;
    //     Vertex& v0 = mesh.V.at(e.Vids.at(0));
    //     Vertex& v1 = mesh.V.at(e.Vids.at(1));
    //     if (v0.N_Vids.size() != 4 || v1.N_Vids.size() != 4) {
    //         EdgeCollapse(e, CollapseOps);
    //         EdgeRotate(e, OptimizationOps);
    //     }
        // Face& f1 = mesh.F.at(e.N_Fids.at(0));
        // Face& f2 = mesh.F.at(e.N_Fids.at(1));
        // double e_length = getLength(e.Vids.at(0), e.Vids.at(1));
        // std::vector<double> d_lengths = {
        //     getLength(f1.Vids.at(0), f1.Vids.at(2)) / sqrt(2),
        //     getLength(f1.Vids.at(1), f1.Vids.at(3)) / sqrt(2),
        //     getLength(f2.Vids.at(0), f2.Vids.at(2)) / sqrt(2),
        //     getLength(f2.Vids.at(1), f2.Vids.at(3)) / sqrt(2),
        // };
        // bool isSmallerThanPrescribed = false;
        // for (auto l: d_lengths) {
        //     if (e_length < l) {
        //         isSmallerThanPrescribed = true;
        //         break;
        //     }
        // }
    // }
    // Vertex Rotate
    // for (auto& v: mesh.V) {
    //     if (v.isBoundary) continue;
    //     if (v.N_Vids.size() == 4) continue;
    //     VertexRotate(v, OptimizationOps);
    // }

    // Diagonal Collapse
    // for (auto& f: mesh.F) {
    //     bool hasBoundaryVertexNeighbor = false;
    //     bool hasRegularVertex = true;
    //     for (auto vid: f.Vids) {
    //         if (mesh.V.at(vid).isBoundary) {
    //             hasBoundaryVertexNeighbor = true;
    //             break;
    //         }
    //         if (mesh.V.at(vid).N_Vids.size() != 4) {
    //             hasRegularVertex = false;
    //             break;
    //         }
    //     }
    //     for (auto fid: f.N_Fids) {
    //         for (auto vid: mesh.F.at(fid).Vids)  {       
    //             if (mesh.V.at(vid).isBoundary) {
    //                 hasBoundaryVertexNeighbor = true;
    //                 break;
    //             }
    //         }
    //     }
    //     if (hasBoundaryVertexNeighbor) continue;
    //     if (hasRegularVertex) continue;
    //     DiagonalCollapse(f, CollapseOps);
    // }
}

void LocalSimplifier::sortPriorities(std::multiset<LocalOperation>& OptimizationOps, std::multiset<LocalOperation, bool(*)(LocalOperation, LocalOperation)>& CollapseOps) {
    // Insert Sort
    // int n = OptimizationOps.size();
    // int i, j;
    // LocalOperation key;
    // for (i = 1; i < n; i++) {
    //     key = OptimizationOps.at(i);
    //     j = i - 1;
    //     while (j >= 0 && OptimizationOps.at(j).profitability > key.profitability) {
    //         OptimizationOps.at(j+1) = OptimizationOps.at(j);
    //         j = j - 1;
    //     }
    //     OptimizationOps.at(j+1) = key;
    // }

    // n = CollapseOps.size();
    // for (i = 1; i < n; i++) {
    //     key = CollapseOps.at(i);
    //     j = i - 1;
    //     while (j >= 0 && CollapseOps.at(j).profitability > key.profitability) {
    //         CollapseOps.at(j+1) = CollapseOps.at(j);
    //         j = j - 1;
    //     }
    //     CollapseOps.at(j+1) = key;
    // }
}

void LocalSimplifier::EdgeRotate(Vertex& v, Edge& e, std::multiset<LocalOperation, bool(*)(LocalOperation, LocalOperation)>& Ops) {
    if (v.isBoundary) return;
    bool isRegular = false;
    if (v.N_Vids.size() == 4) {
        isRegular = true;
        for (auto fid: v.N_Fids) {
            auto& f = mesh.F.at(fid);
            for (auto vid: f.Vids) {
                if (vid == v.id) continue;
                auto& n_v = mesh.V.at(vid);
                if (!n_v.isBoundary) {
                    if (n_v.N_Vids.size() != 4) {
                        isRegular = false;
                    }
                }
            }
        }
        // for (auto vid: v.N_Vids) {
        //     auto& n_v = mesh.V.at(vid);
        //     if (!n_v.isBoundary) {
        //         if (n_v.N_Vids.size() != 4) {
        //             isRegular = false;
        //         }
        //     } else {
        //         if (!n_v.isCorner && n_v.N_Vids.size() != 3) {
        //             isRegular = false;
        //         }
        //     }
        // }
    }
    if (isRegular) return;
    if (v.isCorner) return;
    if (e.isBoundary) return;
    // for (auto vid: e.Vids) {
    //     if (mesh.V.at(vid).isBoundary) return;
    // }

    LocalOperation op;
    op.type = "Edge Rotate";

    Vertex& v1 = mesh.V.at(e.Vids.at(0));
    Vertex& v2 = mesh.V.at(e.Vids.at(1));
    Face& f1 = mesh.F.at(e.N_Fids.at(0));
    Face& f2 = mesh.F.at(e.N_Fids.at(1));

    double area_before = getArea(f1.Vids) + getArea(f2.Vids);
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

    double current_edge_length = getLength(e.Vids.at(0), e.Vids.at(1));
    if (current_edge_length / v.prescribed_length < 1.2) {
        // return;
    }
    // if (current_edge_length < mesh_prescribed_length) {
    //     return;
    // }

    Face new_f1;
    Face new_f2;
    Vertex& new_v1 = mesh.V.at(v1_nVids.at(0));
    new_f1.Vids.push_back(new_v1.id);
    new_f2.Vids.push_back(new_v1.id);
    double new_edge_length1 = 0;
    double diag1 = 0;
    double diag2 = 0;
    double new_diag1 = getLength(v1_nVids.at(0), v1_nVids.at(1));
    double new_diag2 = getLength(v2_nVids.at(0), v2_nVids.at(1));
    double value1 = 0;
    value1 += getLength(v1_nVids.at(0), v1_nVids.at(1));
    value1 += getLength(v2_nVids.at(0), v2_nVids.at(1));
    value1 += getLength(new_v1.id, v2.id);
    if (std::find(new_v1.N_Vids.begin(), new_v1.N_Vids.end(), v2_nVids.at(0)) == new_v1.N_Vids.end()) {
        new_f1.Vids.push_back(v2_nVids.at(1));
        new_f1.Vids.push_back(v2.id);
        new_f1.Vids.push_back(v2_nVids.at(0));

        new_f2.Vids.push_back(v2_nVids.at(0));
        new_f2.Vids.push_back(v1_nVids.at(1));
        new_f2.Vids.push_back(v1.id);

        new_edge_length1 = getLength(v1_nVids.at(0), v2_nVids.at(1));
        diag1 = getLength(v1_nVids.at(1), v2.id);
        diag2 = getLength(v2_nVids.at(0), v1.id);
        value1 += getLength(new_v1.id, v2_nVids.at(0));
        value1 += getLength(v1.id, v2_nVids.at(0));
    } else {
        new_f1.Vids.push_back(v2_nVids.at(0));
        new_f1.Vids.push_back(v2.id);
        new_f1.Vids.push_back(v2_nVids.at(1));
        
        new_f2.Vids.push_back(v2_nVids.at(1));
        new_f2.Vids.push_back(v1_nVids.at(1));
        new_f2.Vids.push_back(v1.id);
        
        new_edge_length1 = getLength(v1_nVids.at(0), v2_nVids.at(0));
        diag1 = getLength(v1_nVids.at(1), v2.id);
        diag2 = getLength(v2_nVids.at(1), v1.id);
        value1 += getLength(new_v1.id, v2_nVids.at(1));
        value1 += getLength(v1.id, v2_nVids.at(1));
    }
    double area_after1 = getArea(new_f1.Vids) + getArea(new_f2.Vids);
    double rotation1_profitability = -1;
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
    double value2 = 0;
    value2 += getLength(v1_nVids.at(0), v1_nVids.at(1));
    value2 += getLength(v2_nVids.at(0), v2_nVids.at(1));
    value2 += getLength(new_v2.id, v2.id);
    if (std::find(new_v2.N_Vids.begin(), new_v2.N_Vids.end(), v2_nVids.at(0)) == new_v2.N_Vids.end()) {
        new_f3.Vids.push_back(v2_nVids.at(1));
        new_f3.Vids.push_back(v2.id);
        new_f3.Vids.push_back(v2_nVids.at(0));

        new_f4.Vids.push_back(v2_nVids.at(0));
        new_f4.Vids.push_back(v1_nVids.at(0));
        new_f4.Vids.push_back(v1.id);
                
        new_edge_length2 = getLength(v1_nVids.at(1), v2_nVids.at(1));
        diag3 = getLength(v1_nVids.at(0), v2.id);
        diag4 = getLength(v2_nVids.at(0), v1.id);
        value2 += getLength(new_v2.id, v2_nVids.at(0));
        value2 += getLength(v1.id, v2_nVids.at(0));
    } else {
        new_f3.Vids.push_back(v2_nVids.at(0));
        new_f3.Vids.push_back(v2.id);
        new_f3.Vids.push_back(v2_nVids.at(1));
        
        new_f4.Vids.push_back(v2_nVids.at(1));
        new_f4.Vids.push_back(v1_nVids.at(0));
        new_f4.Vids.push_back(v1.id);
                
        new_edge_length2 = getLength(v1_nVids.at(1), v2_nVids.at(0));
        diag3 = getLength(v1_nVids.at(0), v2.id);
        diag4 = getLength(v2_nVids.at(1), v1.id);
        value2 += getLength(new_v2.id, v2_nVids.at(1));
        value2 += getLength(v1.id, v2_nVids.at(1));
    }
    double area_after2 = getArea(new_f3.Vids) + getArea(new_f4.Vids);
    double rotation2_profitability = -1;
    if (current_edge_length > new_edge_length2 && diag3 > new_diag1 && diag4 > new_diag2) {
        rotation2_profitability = (current_edge_length - new_edge_length2) + (diag3 - new_diag1) + (diag4 - new_diag2);
    }

    double currentValue = 0;
    currentValue += getLength(v1.id, v2.id);
    currentValue += getLength(v1_nVids.at(0), v2.id);
    currentValue += getLength(v1_nVids.at(1), v2.id);
    currentValue += getLength(v2_nVids.at(0), v1.id);
    currentValue += getLength(v2_nVids.at(1), v1.id);

    // if (currentValue - value1 > currentValue - value2) {
    //     op.newFaces.push_back(new_f1);
    //     op.newFaces.push_back(new_f2);
    //     op.profitability = currentValue - value1;
    // } else {
    //     op.newFaces.push_back(new_f3);
    //     op.newFaces.push_back(new_f4);
    //     op.profitability = currentValue - value2;
    // }
    // if (area_before - area_after1 > area_before - area_after2) {
    //     op.newFaces.push_back(new_f1);
    //     op.newFaces.push_back(new_f2);
    //     op.profitability = area_before - area_after1;
    //     std::cout << "area before: " << area_before << " area after: " << area_after1 << std::endl;
    // } else {
    //     op.newFaces.push_back(new_f3);
    //     op.newFaces.push_back(new_f4);
    //     op.profitability = area_before - area_after2;
    //     std::cout << "area before: " << area_before << " area after: " << area_after2 << std::endl;
    // }
    // if (current_edge_length - new_edge_length1 > current_edge_length - new_edge_length2) {
    //     op.newFaces.push_back(new_f1);
    //     op.newFaces.push_back(new_f2);
    //     op.profitability = current_edge_length - new_edge_length1;
    // } else {
    //     op.newFaces.push_back(new_f3);
    //     op.newFaces.push_back(new_f4);
    //     op.profitability = current_edge_length - new_edge_length2;
    // }
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
        Ops.insert(op);
    }
}

void LocalSimplifier::VertexRotate(Vertex& v, std::multiset<LocalOperation, bool(*)(LocalOperation, LocalOperation)>& Ops) {
    if (v.isBoundary) return;
    bool isRegular = false;
    if (v.N_Vids.size() == 4) {
        isRegular = true;
        for (auto fid: v.N_Fids) {
            auto& f = mesh.F.at(fid);
            for (auto vid: f.Vids) {
                if (vid == v.id) continue;
                auto& n_v = mesh.V.at(vid);
                if (!n_v.isBoundary) {
                    if (n_v.N_Vids.size() != 4) {
                        isRegular = false;
                    }
                }
            }
        }
        // for (auto vid: v.N_Vids) {
        //     auto& n_v = mesh.V.at(vid);
        //     if (!n_v.isBoundary) {
        //         if (n_v.N_Vids.size() != 4) {
        //             isRegular = false;
        //         }
        //     } else {
        //         if (!n_v.isCorner && n_v.N_Vids.size() != 3) {
        //             isRegular = false;
        //         }
        //     }
        // }
    }
    if (isRegular) return;
    // for (auto vid: v.N_Vids) {
    //     if (mesh.V.at(vid).isBoundary) return;
    // }
    
    LocalOperation op;
    op.type = "Vertex Rotate";
    bool isSmallerThanPrescribedLength = false;
    for (auto eid: v.N_Eids) {
        Edge& e = mesh.E.at(eid);
        if (getLength(e.Vids.at(0), e.Vids.at(1)) < mesh_prescribed_length) {
            isSmallerThanPrescribedLength = true;
            break;
        }
    }
    // if (isSmallerThanPrescribedLength) {
    //     return;
    // }
    std::set<size_t> star_edges;
    double sumEdges = 0;
    double sumDiagonals = 0;
    for (auto fid: v.N_Fids) {
        Face& f = mesh.F.at(fid);
        for (auto vid: f.Vids) {
            if (vid == v.id) continue;
            if (std::find(v.N_Vids.begin(), v.N_Vids.end(), vid) == v.N_Vids.end()) {
                sumDiagonals += getLength(vid, v.id);
            }
        }
        star_edges.insert(f.Eids.begin(), f.Eids.end());
        // sumDiagonals += getLength(f.Vids.at(0), f.Vids.at(2));
        // sumDiagonals += getLength(f.Vids.at(1), f.Vids.at(3));
    }
    for (auto vid: v.N_Vids) {
        // Edge& e = mesh.E.at(eid);
        // sumEdges += getLength(e.Vids.at(0), e.Vids.at(1));
        sumEdges += getLength(vid, v.id);
    }
    if ((sumEdges / 4) / v.prescribed_length < 1.2) {
        // return;
    }
    op.profitability = sumEdges - sumDiagonals;

    double area_before = 0;
    double area_after = 0;
    for (auto fid: v.N_Fids) {
        Face& f = mesh.F.at(fid);
        area_before += getArea(f.Vids);
    }
    // if (op.profitability < 0) {
    //     return;
    // }
    // double currentValue = 0;
    // double newValue = 0;
    // bool isProfitable = false;
    for (auto eid: v.N_Eids) {
        std::vector<size_t> newVids;
        Edge& e = mesh.E.at(eid);
        double edge_length = getLength(e.Vids.at(0), e.Vids.at(1));
        // currentValue += edge_length;
        // if (edge_length > mesh_prescribed_length) {
        //     isProfitable = true;
        // }
        newVids.push_back(e.Vids.at(0));
        Face& f1 = mesh.F.at(e.N_Fids.at(0));
        op.canceledFids.insert(f1.id);
        for (int i = 0; i < f1.Vids.size(); i++) {
            if (f1.Vids.at(i) == v.id) continue;
            if (std::find(v.N_Vids.begin(), v.N_Vids.end(), f1.Vids.at(i)) == v.N_Vids.end()) {
                newVids.push_back(f1.Vids.at(i));
                // newValue += getLength(v.id, f1.Vids.at(i));
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
                // newValue += getLength(v.id, f2.Vids.at(i));
                break;
            }
        }
        Face newF;
        newF.Vids = newVids;
        area_after += getArea(newF.Vids);
        op.newFaces.push_back(newF);
    }

    // op.profitability = currentValue - (newValue / 2);
    // op.profitability = area_before - area_after;
    // std::cout << "area before: " << area_before << " area after: " << area_after << std::endl;
    if (op.profitability > 0) {
        Ops.insert(op);
    }
}

void LocalSimplifier::EdgeCollapse(Vertex&v, Edge& e, std::multiset<LocalOperation, bool(*)(LocalOperation, LocalOperation)>& Ops) {
    if (v.isBoundary) return;
    bool isRegular = false;
    if (v.N_Vids.size() == 4) {
        isRegular = true;
        for (auto fid: v.N_Fids) {
            auto& f = mesh.F.at(fid);
            for (auto vid: f.Vids) {
                if (vid == v.id) continue;
                auto& n_v = mesh.V.at(vid);
                if (!n_v.isBoundary) {
                    if (n_v.N_Vids.size() != 4) {
                        isRegular = false;
                    }
                }
            }
        }
        // for (auto vid: v.N_Vids) {
        //     auto& n_v = mesh.V.at(vid);
        //     if (!n_v.isBoundary) {
        //         if (n_v.N_Vids.size() != 4) {
        //             isRegular = false;
        //         }
        //     } else {
        //         if (!n_v.isCorner && n_v.N_Vids.size() != 3) {
        //             isRegular = false;
        //         }
        //     }
        // }
    }
    if (isRegular) return;
    if (v.isCorner) return;
    if (e.isBoundary) return;
    // for (auto vid: e.Vids) {
    //     if (mesh.V.at(vid).isBoundary) return;
    // }

    LocalOperation op;
    op.type = "Edge Collapse";
    // double p_ratio = getLength(e.Vids.at(0), e.Vids.at(1)) / (mesh_prescribed_length);
    // double p_ratio = getLength(e.Vids.at(0), e.Vids.at(1)) / (lambda_p.at(v.id));
    double p_ratio = getLength(e.Vids.at(0), e.Vids.at(1)) / (v.prescribed_length);
    // std::cout << "length: " << getLength(e.Vids.at(0), e.Vids.at(1)) << " p_l: " << v.prescribed_length << " ";
    op.profitability = p_ratio;
    // double p_labmda = p_ratio / p_agg;
    // std::cout << "p_ratio: " << p_ratio << std::endl;
    // std::cout << "p_lambda: " << p_labmda << std::endl;
    // op.profitability = getLength(e.Vids.at(0), e.Vids.at(1)) / (p_labmda * mesh_prescribed_length);
    if (op.profitability > 1) {
        return;
    }
    // double value = getLength(e.Vids.at(0), e.Vids.at(1));
    // bool isSmallerThanEdges = false;
    // for (auto fid: e.N_Fids) {
    //     Face& f = mesh.F.at(fid);
    //     if (getLength(f.Vids.at(0), f.Vids.at(2)) / sqrt(2) < value || getLength(f.Vids.at(1), f.Vids.at(3)) / sqrt(2) < value) {
    //         isSmallerThanEdges = true;
    //         break;
    //     }
    // }

    // if (!isSmallerThanEdges) {
    //     return;
    // }
    // op.profitability = value;
    // if (op.profitability / mesh_prescribed_length < 0.8) { 
        size_t source = v.id;
        size_t target = e.Vids.at(0) != v.id ? e.Vids.at(0) : e.Vids.at(1);
        if (v.isBoundary) {
            source = e.Vids.at(0) != v.id ? e.Vids.at(0) : e.Vids.at(1);
            target = v.id;
        }
        Vertex& v1 = mesh.V.at(source);
        Vertex& v2 = mesh.V.at(target);
        // std::cout << "v1: " << v1.id << " v2: " << v2.id << std::endl;
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
        // op.canceledFids.insert(v2.N_Fids.begin(), v2.N_Fids.end());
        // for (auto fid: v2.N_Fids) {
        //     if (std::find(v1.N_Fids.begin(), v1.N_Fids.end(), fid) != v1.N_Fids.end()) continue;
        //     std::vector<size_t> newVids;
        //     for (auto new_vid: mesh.F.at(fid).Vids) {
        //         if (new_vid == v2.id) {
        //             newVids.push_back(v1.id);
        //         } else {
        //             newVids.push_back(new_vid);
        //         }
        //     }
        //     Face newF;
        //     newF.Vids = newVids;
        //     op.newFaces.push_back(newF);
        // }
        for (auto& f: op.newFaces) {
            for (int i = 0; i < f.Vids.size(); i++) {
                if (f.Vids.at(i) == v1.id) {
                    f.Vids.at(i) = v2.id;
                }
            }
        }
        // if (Ops.size() > 0) {
        //     std::multiset<LocalOperation>::iterator iter = Ops.begin();
        //     LocalOperation op1 = *iter;
        //     if (op.profitability < op1.profitability) {
        //         Ops.insert(op);
        //     }
        // } else {
            Ops.insert(op);
        // }
    // }
}

void LocalSimplifier::DiagonalCollapse(Vertex& v, std::vector<size_t> diag, std::multiset<LocalOperation, bool(*)(LocalOperation, LocalOperation)>& Ops) {
    if (v.isBoundary) return;
    bool isRegular = false;
    if (v.N_Vids.size() == 4) {
        isRegular = true;
        for (auto fid: v.N_Fids) {
            auto& f = mesh.F.at(fid);
            for (auto vid: f.Vids) {
                if (vid == v.id) continue;
                auto& n_v = mesh.V.at(vid);
                if (!n_v.isBoundary) {
                    if (n_v.N_Vids.size() != 4) {
                        isRegular = false;
                    }
                }
            }
        }
        // for (auto vid: v.N_Vids) {
        //     auto& n_v = mesh.V.at(vid);
        //     if (!n_v.isBoundary) {
        //         if (n_v.N_Vids.size() != 4) {
        //             isRegular = false;
        //         }
        //     } else {
        //         if (!n_v.isCorner && n_v.N_Vids.size() != 3) {
        //             isRegular = false;
        //         }
        //     }
        // }
    }
    if (isRegular) return;
    if (v.isCorner) return;
    // for (auto vid: v.N_Vids) {
    //     if (mesh.V.at(vid).isBoundary) return;
    // }
    
    LocalOperation op;
    op.type = "Diagonal Collapse";
    // double value1 = getLength(diag.at(0), diag.at(1)) / ((sqrt(2) * mesh_prescribed_length));
    // double value1 = getLength(diag.at(0), diag.at(1)) / ((sqrt(2) * lambda_p.at(v.id)));
    double value1 = getLength(diag.at(0), diag.at(1)) / ((sqrt(2) * v.prescribed_length));
    op.profitability = value1;
    size_t target_id = diag.at(0);
    size_t source_id = diag.at(1);
    if (mesh.V.at(source_id).isBoundary) {
        target_id = diag.at(1);
        source_id = diag.at(0);
    }
    Vertex& target = mesh.V.at(target_id);
    Vertex& source = mesh.V.at(source_id);
    // double p_ratio1 = value1 / (sqrt(2) * mesh_prescribed_length);
    // double p_labmda1 = (p_ratio1 * p_ratio1) / p_agg;
    // op.profitability = value1 / (p_labmda1 * sqrt(2) * mesh_prescribed_length);
    // double value = value1;
    if (op.profitability > 1) {
        return;
    }
    // bool isSmallerThanEdges = false;
    // for (auto eid: f.Eids) { 
    //     Edge& e = mesh.E.at(eid);
    //     double e_length = getLength(e.Vids.at(0), e.Vids.at(1));
    //     if (value < e_length) {
    //         isSmallerThanEdges = true;
    //         break;
    //     }
    // }
    // if (!isSmallerThanEdges) {
    //     return;
    // }
    // if (op.profitability / mesh_prescribed_length < 0.8) {   
        for (auto fid: source.N_Fids) {
            if (std::find(target.N_Fids.begin(), target.N_Fids.end(), fid) == target.N_Fids.end()) {
                Face& n_f = mesh.F.at(fid);
                Face newF;
                // newF.Vids = n_f.Vids;
                for (int i = 0; i < n_f.Vids.size(); i++) {
                    if (n_f.Vids.at(i) == source.id) {
                        // n_f.Vids.at(i) = target;
                        newF.Vids.push_back(target.id);
                        // break;
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
        
        // if (Ops.size() > 0) {
            // std::multiset<LocalOperation>::iterator iter = Ops.begin();
            // LocalOperation op1 = *iter;
            // if (op.profitability < op1.profitability) {
                // Ops.insert(op);
            // }
        // } else {
            Ops.insert(op);
        // }
    // }
}

void LocalSimplifier::RemoveDoublets(std::set<size_t>& canceledFids) {
    for (auto& v: mesh.V) {
        if (v.isBoundary) continue;
        if (v.N_Fids.size() == 2) {
            Face& f1 = mesh.F.at(v.N_Fids.at(0));
            Face& f2 = mesh.F.at(v.N_Fids.at(1));
            canceledFids.insert(f1.id);
            size_t diag_vid = 0;
            for (auto vid: f1.Vids) {
                if (vid == v.id) continue;
                if (std::find(v.N_Vids.begin(), v.N_Vids.end(), vid) == v.N_Vids.end()) {
                    diag_vid = vid;
                    break;
                }
            }
            std::vector<size_t> newVids;
            for (auto vid: f2.Vids) {
                if (vid == v.id) {
                    newVids.push_back(diag_vid);
                } else {
                    newVids.push_back(vid);
                }
            }
            f2.Vids = newVids;
        }
    }
    for (auto& f: mesh.F) {
        int doublet_v_count = 0;
        for (auto vid: f.Vids) {
            if (!mesh.V.at(vid).isBoundary && mesh.V.at(vid).N_Fids.size() == 2) {
                doublet_v_count += 1;
            }
        }
        if (doublet_v_count > 1) {
            canceledFids.insert(f.id);
        }
    }
}

double LocalSimplifier::getLength(size_t vid1, size_t vid2) {
    Vertex& v1 = mesh.V.at(vid1);
    Vertex& v2 = mesh.V.at(vid2);
    double length = glm::length(glm::dvec3(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z));
    return length;
}

double LocalSimplifier::getArea(std::vector<size_t>& Vids) {
    Vertex& a = mesh.V.at(Vids.at(0));
    Vertex& b = mesh.V.at(Vids.at(1));
    Vertex& c = mesh.V.at(Vids.at(2));
    Vertex& d = mesh.V.at(Vids.at(3));
    glm::dvec3 d1(c.x - a.x, c.y - a.y, c.z - a.z);
    glm::dvec3 d2(d.x - b.x, d.y - b.y, d.z - b.z);
    return glm::length(glm::cross(d1, d2)) / 2;
}

void LocalSimplifier::RemoveBoundaryDegenerates() {
    std::set<size_t> canceledFids;
    std::vector<Face> newFaces;
    for (int i = 0; i < mesh.F.size(); i++) {
        // std::cout << i << " " << mesh.F.size() << std::endl;
        // mesh.SetOneRingNeighborhood();
        newFaces.clear();
        canceledFids.clear();
        Face& f = mesh.F.at(i);
        if (!f.isBoundary) continue;
        int boundaryVerticesCount = 0;
        for (auto vid: f.Vids) {
            if (mesh.V.at(vid).isBoundary) {
                boundaryVerticesCount += 1;
            }
        }
        if (boundaryVerticesCount == 3) {
            std::cout << "a three boundary vertex face detected" << std::endl;
            size_t index = 0;
            for (int i = 0; i < f.Vids.size(); i++) {
                if (!mesh.V.at(f.Vids.at(i)).isBoundary) {
                    index = (i + 1) % f.Vids.size();
                    break;
                }
            }
            Vertex& a = mesh.V.at(f.Vids.at(index));
            Vertex& b = mesh.V.at(f.Vids.at((index + 1) % f.Vids.size()));
            Vertex& c = mesh.V.at(f.Vids.at((index + 2) % f.Vids.size()));
            if (b.N_Fids.size() != 1) continue;

            if (a.N_Fids.size() > 2) {
                std::cout << "rotating around a" << std::endl;
                size_t edgeId = 0;
                for (auto eid: a.N_Eids) {
                    Edge& e = mesh.E.at(eid);
                    if (e.isBoundary) continue;
                    if (e.N_Fids.at(0) == f.id || e.N_Fids.at(1) == f.id) {
                        std::cout << "e nfids: " << e.N_Fids.size() << std::endl;
                        Face f1;
                        Face f2;
                        Face& ef1 = mesh.F.at(e.N_Fids.at(0));
                        Face& ef2 = mesh.F.at(e.N_Fids.at(1));
                        std::vector<size_t> newVids;
                        for (int i = 0; i < ef1.Vids.size(); i++) {
                            if (std::find(e.Vids.begin(), e.Vids.end(), ef1.Vids.at(i)) != e.Vids.end() &&
                                std::find(e.Vids.begin(), e.Vids.end(), ef1.Vids.at((i + 1) % ef1.Vids.size())) != e.Vids.end()) {
                                newVids.push_back(ef1.Vids.at((i + 1) % ef1.Vids.size()));
                                newVids.push_back(ef1.Vids.at((i + 2) % ef1.Vids.size()));
                                newVids.push_back(ef1.Vids.at((i + 3) % ef1.Vids.size()));
                            }
                        }
                        for (int i = 0; i < ef2.Vids.size(); i++) {
                            if (std::find(e.Vids.begin(), e.Vids.end(), ef2.Vids.at(i)) != e.Vids.end() &&
                                std::find(e.Vids.begin(), e.Vids.end(), ef2.Vids.at((i + 1) % ef2.Vids.size())) != e.Vids.end()) {
                                newVids.push_back(ef2.Vids.at((i + 1) % ef2.Vids.size()));
                                newVids.push_back(ef2.Vids.at((i + 2) % ef2.Vids.size()));
                                newVids.push_back(ef2.Vids.at((i + 3) % ef2.Vids.size()));
                            }
                        }
                        size_t start = 0;
                        for (int i = 0; i < newVids.size(); i++) {
                            if (newVids.at(i) == b.id) {
                                start = i;
                                break;
                            }
                        }

                        f1.Vids.push_back(newVids.at(start % newVids.size()));
                        f1.Vids.push_back(newVids.at((start + 1) % newVids.size()));
                        f1.Vids.push_back(newVids.at((start + 2) % newVids.size()));
                        f1.Vids.push_back(newVids.at((start + 3) % newVids.size()));
                        f2.Vids.push_back(newVids.at((start + 3) % newVids.size()));
                        f2.Vids.push_back(newVids.at((start + 4) % newVids.size()));
                        f2.Vids.push_back(newVids.at((start + 5) % newVids.size()));
                        f2.Vids.push_back(newVids.at(start % newVids.size()));
                        newFaces.push_back(f1);
                        newFaces.push_back(f2);
                        canceledFids.insert(ef1.id);
                        canceledFids.insert(ef2.id);
                        break;
                    }
                }
            } else if (c.N_Fids.size() > 2) {
                std::cout << "rotating around c" << std::endl;
                size_t edgeId = 0;
                for (auto eid: c.N_Eids) {
                    Edge& e = mesh.E.at(eid);
                    if (e.isBoundary) continue;
                    if (e.N_Fids.at(0) == f.id || e.N_Fids.at(1) == f.id) {
                        std::cout << "e nfids: " << e.N_Fids.size() << std::endl;
                        Face f1;
                        Face f2;
                        Face& ef1 = mesh.F.at(e.N_Fids.at(0));
                        Face& ef2 = mesh.F.at(e.N_Fids.at(1));
                        std::vector<size_t> newVids;
                        for (int i = 0; i < ef1.Vids.size(); i++) {
                            if (std::find(e.Vids.begin(), e.Vids.end(), ef1.Vids.at(i)) != e.Vids.end() &&
                                std::find(e.Vids.begin(), e.Vids.end(), ef1.Vids.at((i + 1) % ef1.Vids.size())) != e.Vids.end()) {
                                newVids.push_back(ef1.Vids.at((i + 1) % ef1.Vids.size()));
                                newVids.push_back(ef1.Vids.at((i + 2) % ef1.Vids.size()));
                                newVids.push_back(ef1.Vids.at((i + 3) % ef1.Vids.size()));
                            }
                        }
                        for (int i = 0; i < ef2.Vids.size(); i++) {
                            if (std::find(e.Vids.begin(), e.Vids.end(), ef2.Vids.at(i)) != e.Vids.end() &&
                                std::find(e.Vids.begin(), e.Vids.end(), ef2.Vids.at((i + 1) % ef2.Vids.size())) != e.Vids.end()) {
                                newVids.push_back(ef2.Vids.at((i + 1) % ef2.Vids.size()));
                                newVids.push_back(ef2.Vids.at((i + 2) % ef2.Vids.size()));
                                newVids.push_back(ef2.Vids.at((i + 3) % ef2.Vids.size()));
                            }
                        }
                        size_t start = 0;
                        for (int i = 0; i < newVids.size(); i++) {
                            if (newVids.at(i) == b.id) {
                                start = i;
                                break;
                            }
                        }

                        f1.Vids.push_back(newVids.at(start % newVids.size()));
                        f1.Vids.push_back(newVids.at((start + 1) % newVids.size()));
                        f1.Vids.push_back(newVids.at((start + 2) % newVids.size()));
                        f1.Vids.push_back(newVids.at((start + 3) % newVids.size()));
                        f2.Vids.push_back(newVids.at((start + 3) % newVids.size()));
                        f2.Vids.push_back(newVids.at((start + 4) % newVids.size()));
                        f2.Vids.push_back(newVids.at((start + 5) % newVids.size()));
                        f2.Vids.push_back(newVids.at(start % newVids.size()));
                        newFaces.push_back(f1);
                        newFaces.push_back(f2);
                        canceledFids.insert(ef1.id);
                        canceledFids.insert(ef2.id);
                        break;
                    }
                }
            } else {
                size_t source_id1 = a.id;
                size_t target_id = b.id;
                size_t source_id2 = c.id;
                if (a.isCorner) {
                    source_id1 = b.id;
                    source_id2 = c.id;
                    target_id = a.id;
                }
                if (c.isCorner) {
                    source_id1 = a.id;
                    source_id2 = b.id;
                    target_id = c.id;
                }
                Vertex& source1 = mesh.V.at(source_id1);
                Vertex& source2 = mesh.V.at(source_id2);
                Vertex& target = mesh.V.at(target_id);

                for (auto nfid: source1.N_Fids) {
                    if (nfid == f.id) continue;
                    Face& nf = mesh.F.at(nfid);
                    for (int i = 0; i < nf.Vids.size(); i++) {
                        if (nf.Vids.at(i) == source1.id) {
                            nf.Vids.at(i) = target.id;
                        }
                    }
                }

                for (auto nfid: source2.N_Fids) {
                    if (nfid == f.id) continue;
                    Face& nf = mesh.F.at(nfid);
                    for (int i = 0; i < nf.Vids.size(); i++) {
                        if (nf.Vids.at(i) == source2.id) {
                            nf.Vids.at(i) = target.id;
                        }
                    }
                }
                canceledFids.insert(f.id);
            }
        } else if (boundaryVerticesCount == 4) {
            // bool isRegular = true;
            // for (auto vid: f.Vids) {
            //     if (mesh.V.at(vid).N_Fids.size() == 1) isRegular = false;
            // }
            // if (isRegular) continue;
            // size_t target_id = 0;
            // size_t source_id = 0;
            // for (auto vid: f.Vids) {
            //     if (mesh.V.at(vid).N_Fids.size() == 1 && mesh.V.at(vid).isCorner) {
            //         target_id = vid;
            //     }
            //     if (mesh.V.at(vid).N_Fids.size() >= 2) {
            //         source_id = vid;
            //     }
            // }
            // Vertex& source = mesh.V.at(source_id);
            // for (auto fid: source.N_Fids) {
            //     if (fid == f.id) continue;
            //     Face& n_f = mesh.F.at(fid);
            //     for (int j = 0; j < n_f.Vids.size(); j++) {
            //         if (n_f.Vids.at(j) == source_id) {
            //             n_f.Vids.at(j) = target_id;
            //         }
            //     }
            // }
            // canceledFids.insert(f.id);
        }
        if (!canceledFids.empty()) {
            std::cout << canceledFids.size() << std::endl;
            std::cout << newFaces.size() << std::endl;
            for (auto nf: newFaces) {
                nf.id = mesh.F.size();
                mesh.F.push_back(nf);
                std::cout << "f: ";
                for (auto fvid: nf.Vids) {
                    std::cout << fvid << " ";
                }
                std::cout << std::endl;
            }
            std::cout << "*************************************************" << std::endl;
            update(canceledFids);
            init();
            bool doublets_removed = false;
            while (!doublets_removed) {
                RemoveDoublets(canceledFids);
                if (!canceledFids.empty()) {
                    update(canceledFids);
                    init();
                } else {
                    doublets_removed = true;
                }
            }
            mesh.SetOneRingNeighborhood();
            i = 0;
            // n = mesh.F.size();
            // break;
        }
    }


    // std::ofstream ofs("degenerateFaces.vtk");
	// ofs << "# vtk DataFile Version 3.0\n"
	// 	<< "degenerateFaces.vtk\n"
	// 	<< "ASCII\n\n"
	// 	<< "DATASET UNSTRUCTURED_GRID\n";
	// ofs << "POINTS " << mesh.V.size() << " double\n";

    // for (auto& v: mesh.V) {
    //     ofs << v.x << " " << v.y << " " << v.z << std::endl;
    // }

	// ofs << "CELLS " << degenerateFaces.size() << " " << degenerateFaces.size() * 5 << std::endl;
    // for (auto fid: degenerateFaces) {
    //     Face& f = mesh.F.at(fid);
    //     ofs << "4 " << f.Vids.at(0) << " " << f.Vids.at(1) << " " << f.Vids.at(2) << " " << f.Vids.at(3) << std::endl; 
    // }
	// ofs << "CELL_TYPES " << degenerateFaces.size() << "\n";
    // for (auto eid: degenerateFaces) {
    //     ofs << "9 " << std::endl; 
    // }
}

/*void LocalSimplifier::smoothMesh(std::set<size_t>& Vertices) {
    mesh_prescribed_length = getPrescribedLength();
    mesh.SetOneRingNeighborhood();
    std::set<size_t> currentVertices;
    if (Vertices.size() < mesh.V.size()) {
        for (auto vid: Vertices) {
            if (vid >= mesh.V.size()) continue;
            currentVertices.insert(vid);
            currentVertices.insert(mesh.V.at(vid).N_Vids.begin(), mesh.V.at(vid).N_Vids.end());
        }
    } else {
        for (auto vid: Vertices) {
            currentVertices.insert(vid);
        }
    }
    int it = 0;
    while (it < 1000) {
        std::vector<double> vertex_weights(currentVertices.size(), 0);
        std::vector<glm::dvec3> new_coords(currentVertices.size(), glm::dvec3(0.0, 0.0, 0.0));

        int v_pos = 0;
        for (auto id: currentVertices) {
            Vertex& v_i = mesh.V.at(id);
            if (v_i.isBoundary) {
                v_pos++;
                continue;
            }
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

                glm::dvec3 V_j(v_i.x - v_j.x, v_i.y - v_j.y, v_i.z - v_j.z);
                glm::dvec3 V_j_minus_1(v_j_prev.x - v_j.x, v_j_prev.y - v_j.y, v_j_prev.z - v_j.z);
                glm::dvec3 V_j_plus_1(v_j_next.x - v_j.x, v_j_next.y - v_j.y, v_j_next.z - v_j.z);

                double alpha1 = acos(glm::dot(V_j, V_j_plus_1) / (glm::length(V_j) * glm::length(V_j_plus_1)));
                double alpha2 = acos(glm::dot(V_j, V_j_minus_1) / (glm::length(V_j) * glm::length(V_j_minus_1)));
                double beta = (alpha2 - alpha1) / 2;
                neighbor_coords.at(j).x = (v_j.x + (V_j.x * cos(beta)) - (V_j.y * sin(beta))); 
                neighbor_coords.at(j).y = (v_j.y + (V_j.x * sin(beta)) + (V_j.y * cos(beta)));
                neighbor_weights.at(j) = fabs(beta);

                weight_agg += neighbor_weights.at(j);
            }
            if (weight_agg == 0) {
                new_coords.at(v_pos).x = 0;
                new_coords.at(v_pos).y = 0;
                v_pos++;
                continue;
            }
            for (int j = 0; j < one_ring_neighbors.size(); j++) {
                glm::dvec3 current_v = neighbor_coords.at(j);
                double weight = (1 - neighbor_weights.at(j) / weight_agg);
                vertex_weights.at(v_pos) += weight;
                new_coords.at(v_pos).x += (weight * current_v.x);
                new_coords.at(v_pos).y += (weight * current_v.y);
            }
            new_coords.at(v_pos).x = (new_coords.at(v_pos).x / vertex_weights.at(v_pos)) - v_i.x;
            new_coords.at(v_pos).y = (new_coords.at(v_pos).y / vertex_weights.at(v_pos)) - v_i.y;
            v_pos++;
        }
        v_pos = 0;
        for (auto id: currentVertices) {
            Vertex& v = mesh.V.at(id);
            if (v.isBoundary) {
                v_pos++;
                continue;
            }
            v.x += new_coords.at(v_pos).x;
            v.y += new_coords.at(v_pos).y;
            v_pos++;
        }
        it++;
    }
}*/

void LocalSimplifier::smoothMeshGlobal() {
    // mesh_prescribed_length = getPrescribedLength();
    double tau = 1e-6;
    // mesh.SetOneRingNeighborhood();
    int it = 0;
    while (it < 10000) {
        // std::cout << "it: " << it << std::endl;
        std::vector<double> vertex_weights(mesh.V.size(), 0);
        std::vector<glm::dvec3> new_coords(mesh.V.size(), glm::dvec3(0.0, 0.0, 0.0));

        int v_pos = 0;
        for (auto& v: mesh.V) {
            if (v.isBoundary) {
                v_pos++;
                continue;
            }
            double n = v.N_Vids.size();
            double weight = 1 / n;
            double weight_agg = 0;
            double x_agg = 0;
            double y_agg = 0;
            for (int j = 0; j < n; j++) {
                Vertex& v_prime = mesh.V.at(v.N_Vids.at(j));
                x_agg += (weight * v_prime.x);
                y_agg += (weight * v_prime.y);
                weight_agg += weight;
            }
            new_coords.at(v_pos).x = (x_agg / weight_agg) - v.x;
            new_coords.at(v_pos).y = (y_agg / weight_agg) - v.y;
            v_pos++;
        }
        v_pos = 0;
        double avg_delta = 0;
        for (auto& v: mesh.V) {
            glm::dvec3 prev(v.x, v.y, v.z);
            if (v.isBoundary) {
                v_pos++;
                continue;
            }
            v.x += new_coords.at(v_pos).x;
            v.y += new_coords.at(v_pos).y;
            glm::dvec3 current(v.x, v.y, v.z);
            avg_delta += glm::length(prev - current);
            v_pos++;
        }
        avg_delta = avg_delta / mesh.V.size();
        if (avg_delta < tau) {
            std::cout << "smoothing stopped at " << it << "th iteration" << std::endl;
            break;
        }
        it++;
    }
}

void LocalSimplifier::smoothMesh(int iters_, bool global) {
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
                double factor = glm::length(direction) - v.prescribed_length;
                // if (std::find(v.N_Vids.begin(), v.N_Vids.end(), vid) == v.N_Vids.end()) {
                //     // glm::dvec3 new_coord(v_prime.x, v_prime.y, v_prime.z);
                //     // double factor = glm::length(direction) - (sqrt(2) * lambda_p.at(v.id));
                //     factor = glm::length(direction) - (sqrt(2) * v.prescribed_length);
                // }
                glm::dvec3 direction_normalized = glm::normalize(direction);
                new_coords.at(v_pos) += (direction_normalized * factor);
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
    std::cout << "Smoothed mesh" << std::endl;
    // init();
    std::cout << "initialized mesh" << std::endl;
}
// global prescribed length for overall operations
// local prescribed length for updating positions or local gaussian curvature
// singularity restriction
// area of quad: 1/2 d1xd2

// send email for discount codes
// shorten algorithm
// register for conferences
// continue working on project

// original Kaoji framework, 2nd framework with integrated operations, 3rd framework with local smoothing.
// continue working on local smoothing

// 1: edge tolerance: 25%, without adaptive simplification, constant ideal edge length, local smoothing iterations: 1000
// terminates after 60 iterations, output faces: 373, local smoothing with changed vertices only


// square:
// time taken:  151.39 seconds
// original: #v = 537, #f = 495, #singularities: 269
// simplified: #v: 292, #f: 250, #singularities: 102

// holes1:
// time taken: 98.996 seconds
// original: #v: 470, #f: 420, #singularities: 230
// out: #v: 248, #f: 198, #singulatities: 68

// holes2_square2: 
// time taken: 2988.82 seconds
// original: #v: 3671, #f: 3465, #singularities: 1825
// out: #v: 2155, #f: 1949, #singulatities: 927

// 5101509533
