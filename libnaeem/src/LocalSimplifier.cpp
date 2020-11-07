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

bool OpSortAscending(LocalOperation op1, LocalOperation op2) {
    return op1.profitability > op2.profitability;
}

bool OpSortDescending(LocalOperation op1, LocalOperation op2) {
    return op1.profitability < op2.profitability;
}

LocalSimplifier::LocalSimplifier(Mesh& mesh) : Simplifier(mesh) {}

LocalSimplifier::~LocalSimplifier() {}

void LocalSimplifier::Simplify() {
    // init();
    // std::vector<LocalOperation> OptimizationOps;
    // std::vector<LocalOperation> CollapseOps;
    // getOperationsPriorities(OptimizationOps, CollapseOps);
    // sortPriorities(OptimizationOps, CollapseOps);
    // for (auto op: OptimizationOps) {
    //     std::cout << op.profitability << std::endl;
    // }
    // for (auto op: CollapseOps) {
    //     std::cout << op.profitability << std::endl;
    // }
    get_feature();
    int it = 0;
    init();
    while (it < iters) {
        std::cout << "it: " << it << std::endl;
        std::cout << "V: " << mesh.V.size() << std::endl;
        std::cout << "E: " << mesh.E.size() << std::endl;
        std::cout << "F: " << mesh.F.size() << std::endl;
        std::set<size_t> canceledFids;
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
        // break;
        bool(*fn_pt1)(LocalOperation, LocalOperation) = OpSortAscending;
        bool(*fn_pt2)(LocalOperation, LocalOperation) = OpSortDescending;
        std::multiset<LocalOperation, bool(*)(LocalOperation, LocalOperation)> OptimizationOps(fn_pt1);
        std::multiset<LocalOperation, bool(*)(LocalOperation, LocalOperation)> CollapseOps(fn_pt2);
        // std::vector<LocalOperation> OptimizationOps;
        // std::vector<LocalOperation> CollapseOps;
        getOperationsPriorities(OptimizationOps, CollapseOps);
        std::cout << "priorities set" << std::endl;
        std::cout << "# optimization ops: " << OptimizationOps.size() << std::endl;
        std::cout << "# collapse ops: " << CollapseOps.size() << std::endl;
        // std::sort(OptimizationOps.begin(), OptimizationOps.end(), OpSort);
        // std::sort(CollapseOps.begin(), CollapseOps.end(), OpSort);
        std::cout << "Sorted" << std::endl;
        // std::multiset<LocalOperation>::iterator iter;
        // for (iter = OptimizationOps.begin(); iter != OptimizationOps.end(); iter++) {
        //     LocalOperation op = *iter;
        //     std::cout << op.profitability << " ";
        // }
        // std::cout << std::endl;
        // for (iter = CollapseOps.begin(); iter != CollapseOps.end(); iter++) {
        //     LocalOperation op = *iter;
        //     std::cout << op.profitability << " ";
        // }
        // std::cout << std::endl;
        // return;
        // sortPriorities(OptimizationOps, CollapseOps);
        std::set<size_t> processedFids;
        for(auto op: OptimizationOps) {
            if (op.profitability < 0) continue;
            // if (op.type.compare("Edge Rotate") == 0) continue;
            // if (op.type.compare("Edge Collapse") == 0 || op.type.compare("Diagonal Collapse") == 0) continue;
            /*if (op.type.compare("Edge Collapse") == 0 || op.type.compare("Diagonal Collapse") == 0) {
                bool overlappedOp = false;
                for (auto id: op.canceledFids) {
                    if (std::find(processedFids.begin(), processedFids.end(), id) != processedFids.end()) {
                        overlappedOp = true;
                        break;
                    }
                }
                if (overlappedOp) continue;
                std::cout << "it: " << it << " Collapse Op: " << op.type <<  std::endl;
                for (auto f: op.newFaces) {
                    f.id = mesh.F.size();
                    std::cout << "f vids: " << f.Vids.size() << ", "; 
                    for (auto id: f.Vids) {
                        std::cout << id << " ";
                    }
                    std::cout << std::endl;
                    mesh.F.push_back(f);
                }
                canceledFids.insert(op.canceledFids.begin(), op.canceledFids.end());
                processedFids.insert(op.canceledFids.begin(), op.canceledFids.end());
                for (auto fid: op.canceledFids) {
                    Face& f = mesh.F.at(fid);
                    processedFids.insert(f.N_Fids.begin(), f.N_Fids.end());
                }
                continue;
            }*/
            // if (op.profitability > 0) {
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
                    std::cout << "f vids: " << f.Vids.size() << ", "; 
                    for (auto id: f.Vids) {
                        std::cout << id << " ";
                    }
                    std::cout << std::endl;
                    for (int i = 0; i < f.Vids.size() - 1; i++) {
                        for (int j = i + 1; j < f.Vids.size(); j++) {
                            if (f.Vids.at(i) == f.Vids.at(j)) {
                                std::cout << "ERROR: VERTICES REPEATED" << std::endl;
                            }
                        }
                    }
                    mesh.F.push_back(f);
                }
                canceledFids.insert(op.canceledFids.begin(), op.canceledFids.end());
                // std::cout << "fids: "; 
                // for (auto id: op.canceledFids) {
                //     std::cout << id << " ";
                // }
                // std::cout << std::endl;
                // std::cout << "canceled Faces: " << canceledFids.size() << std::endl;
                processedFids.insert(op.canceledFids.begin(), op.canceledFids.end());
                for (auto fid: op.canceledFids) {
                    Face& f = mesh.F.at(fid);
                    processedFids.insert(f.N_Fids.begin(), f.N_Fids.end());
                }
                // break;
            // }
        }
        std::cout << "canceled Faces: " << canceledFids.size() << " ";
        for (auto id: canceledFids) {
            std::cout << id << " ";
        }
        std::cout << std::endl;
        // if (!canceledFids.empty()) {
        //     update(canceledFids);
        //     init();
        // }
        // std::cout << "profitability: ";
        // for (auto op: CollapseOps) {
        //     std::cout << op.profitability << " ";
        // }
        // std::cout << std::endl;
        for (auto op: CollapseOps) {
            // if (op.type.compare("Edge Rotate") == 0 || op.type.compare("Vertex Rotate") == 0) continue;
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
                std::cout << "f vids: " << f.Vids.size() << ", "; 
                for (auto id: f.Vids) {
                    std::cout << id << " ";
                }
                std::cout << std::endl;
                for (int i = 0; i < f.Vids.size() - 1; i++) {
                    for (int j = i + 1; j < f.Vids.size(); j++) {
                        if (f.Vids.at(i) == f.Vids.at(j)) {
                            std::cout << "ERROR: VERTICES REPEATED" << std::endl;
                        }
                    }
                }
                mesh.F.push_back(f);
            }
            canceledFids.insert(op.canceledFids.begin(), op.canceledFids.end());
            // processedFids.insert(op.canceledFids.begin(), op.canceledFids.end());
            // for (auto fid: op.canceledFids) {
            //     Face& f = mesh.F.at(fid);
            //     processedFids.insert(f.N_Fids.begin(), f.N_Fids.end());
            // }
            break;
            // std::cout << "canceled Faces: " << canceledFids.size() << std::endl;
            // foundCollapseOp = true;
        }
        std::cout << "canceled Faces: " << canceledFids.size() << std::endl;
        for (auto id: canceledFids) {
            std::cout << id << " ";
        }
        std::cout << std::endl;
        if (!canceledFids.empty()) {
            std::cout << "updating" << std::endl;
            update(canceledFids);
            init();
            std::cout << "new faces: " << mesh.F.size() << std::endl;
        }
        // DoubletSimplifier doubletSimplifier2(mesh);
        // doubletSimplifier2.RunCollective(canceledFids);
        // RemoveDoublets(canceledFids);
        // if (!canceledFids.empty()) {
        //     update(canceledFids);
        //     init();
        // }
        doublets_removed = false;
        while (!doublets_removed) {
            RemoveDoublets(canceledFids);
            if (!canceledFids.empty()) {
                update(canceledFids);
                init();
            } else {
                doublets_removed = true;
            }
        }
        std::cout << "üpdated after doublets removal" << std::endl;
        it++;
    }
}

void LocalSimplifier::getOperationsPriorities(std::multiset<LocalOperation, bool(*)(LocalOperation, LocalOperation)>& OptimizationOps, std::multiset<LocalOperation, bool(*)(LocalOperation, LocalOperation)>& CollapseOps) {
    // Edge Rotate
    // int i = 0;
    for (auto& e: mesh.E) {
        bool hasBoundaryVertexNeighbor = false;
        for (auto vid: e.Vids) {
            if (mesh.V.at(vid).isBoundary) {
                hasBoundaryVertexNeighbor = true;
                break;
            }
        }
        for (auto vid: e.N_Vids) {
            if (mesh.V.at(vid).isBoundary) {
                hasBoundaryVertexNeighbor = true;
                break;
            }
        }
        if (hasBoundaryVertexNeighbor) continue;
        EdgeRotate(e, OptimizationOps);
        EdgeCollapse(e, CollapseOps);
        // i+=1;
        // if (i > 5) break;
    }
    std::cout << "Edge rotate set" << std::endl;
    // Vertex Rotate
    for (auto& v: mesh.V) {
        if (v.isBoundary) continue;
        VertexRotate(v, OptimizationOps);
    }
    std::cout << "Vertex rotate set" << std::endl;

    // Edge Collapse
    /*for (auto& e: mesh.E) {
        bool hasBoundaryVertexNeighbor = false;
        for (auto vid: e.Vids) {
            if (mesh.V.at(vid).isBoundary) {
                hasBoundaryVertexNeighbor = true;
                break;
            }
        }
        for (auto vid: e.N_Vids) {
            if (mesh.V.at(vid).isBoundary) {
                hasBoundaryVertexNeighbor = true;
                break;
            }
        }
        if (hasBoundaryVertexNeighbor) continue;
    }*/
    std::cout << "Edge collapse set" << std::endl;

    // Diagonal Collapse
    for (auto& f: mesh.F) {
        bool hasBoundaryVertexNeighbor = false;
        for (auto vid: f.Vids) {
            if (mesh.V.at(vid).isBoundary) {
                hasBoundaryVertexNeighbor = true;
                break;
            }
        }
        for (auto vid: f.N_Vids) {
            if (mesh.V.at(vid).isBoundary) {
                hasBoundaryVertexNeighbor = true;
                break;
            }
        }
        if (hasBoundaryVertexNeighbor) continue;
        DiagonalCollapse(f, CollapseOps);
    }
    std::cout << "Diagonal collapse set" << std::endl;
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


/*void LocalSimplifier::EdgeRotate(Edge& e, std::vector<LocalOperation>& Ops) {
    LocalOperation op;
    op.type = "Edge Rotate";
    // std::set<size_t> canceledFids;
    std::vector<size_t> Vids;
    Face& f1 = mesh.F.at(e.N_Fids.at(0));
    double currentValue = 0;
    currentValue += getLength(e.Vids.at(0), e.Vids.at(1));
    for (int i = 0; i < f1.Vids.size(); i++) {
        size_t id1 = f1.Vids.at(i);
        size_t id2 = f1.Vids.at((i + 1) % f1.Vids.size());
        if (std::find(e.Vids.begin(), e.Vids.end(), id1) != e.Vids.end() && std::find(e.Vids.begin(), e.Vids.end(), id2) != e.Vids.end()) {
            currentValue += getLength(f1.Vids.at(i), f1.Vids.at((i + 2) % f1.Vids.size()));
            currentValue += getLength(f1.Vids.at((i + 1) % f1.Vids.size()), f1.Vids.at((i + 3) % f1.Vids.size()));
            Vids.push_back(f1.Vids.at((i + 1) % f1.Vids.size()));
            Vids.push_back(f1.Vids.at((i + 2) % f1.Vids.size()));
            Vids.push_back(f1.Vids.at((i + 3) % f1.Vids.size()));
            break;
        }
    }
    Face& f2 = mesh.F.at(e.N_Fids.at(1));
    for (int i = 0; i < f2.Vids.size(); i++) {
        size_t id1 = f2.Vids.at(i);
        size_t id2 = f2.Vids.at((i + 1) % f2.Vids.size());
        if (std::find(e.Vids.begin(), e.Vids.end(), id1) != e.Vids.end() && std::find(e.Vids.begin(), e.Vids.end(), id2) != e.Vids.end()) {
            currentValue += getLength(f2.Vids.at(i), f2.Vids.at((i + 2) % f2.Vids.size()));
            currentValue += getLength(f2.Vids.at((i + 1) % f2.Vids.size()), f2.Vids.at((i + 3) % f2.Vids.size()));
            Vids.push_back(f2.Vids.at((i + 1) % f2.Vids.size()));
            Vids.push_back(f2.Vids.at((i + 2) % f2.Vids.size()));
            Vids.push_back(f2.Vids.at((i + 3) % f2.Vids.size()));
            break;
        }
    }
    double value1 = 0;
    value1 += getLength(Vids.at(1), Vids.at(4));
    value1 += getLength(Vids.at(1), Vids.at(3));
    value1 += getLength(Vids.at(2), Vids.at(4));
    value1 += getLength(Vids.at(4), Vids.at(0));
    value1 += getLength(Vids.at(5), Vids.at(1));

    double value2 = 0;
    value2 += getLength(Vids.at(2), Vids.at(5));
    value2 += getLength(Vids.at(2), Vids.at(4));
    value2 += getLength(Vids.at(3), Vids.at(5));
    value2 += getLength(Vids.at(5), Vids.at(1));
    value2 += getLength(Vids.at(0), Vids.at(2));

    int start = 0;
    if (currentValue - value1 > currentValue - value2) {
        start = 1;
        op.profitability = currentValue - value1;
    } else {
        start = 2;
        op.profitability = currentValue - value2;
    }

    Face new_f1;
    // new_f1.id = mesh.F.size();
    for (int i = start; i < start + 4; i++) {
        new_f1.Vids.push_back(Vids.at(i % Vids.size()));
    }
    op.newFaces.push_back(new_f1);
    // mesh.F.push_back(new_f1);
    Face new_f2;
    // new_f2.id = mesh.F.size();
    start += 3;
    for (int i = start; i < start + 4; i++) {
        new_f2.Vids.push_back(Vids.at(i % Vids.size()));
    }
    op.newFaces.push_back(new_f2);
    // mesh.F.push_back(new_f2);
    op.canceledFids.insert(e.N_Fids.begin(), e.N_Fids.end());
    // canceledFids.insert(e.N_Fids.begin(), e.N_Fids.end());
    Ops.push_back(op);
}*/

void LocalSimplifier::EdgeRotate(Edge& e, std::multiset<LocalOperation, bool(*)(LocalOperation, LocalOperation)>& Ops) {
    LocalOperation op;
    op.type = "Edge Rotate";
    // std::set<size_t> canceledFids;

    Vertex& v1 = mesh.V.at(e.Vids.at(0));
    Vertex& v2 = mesh.V.at(e.Vids.at(1));
    Face& f1 = mesh.F.at(e.N_Fids.at(0));
    Face& f2 = mesh.F.at(e.N_Fids.at(1));

    std::vector<size_t> v1_nVids;
    std::vector<size_t> v2_nVids;
    // std::cout << v1.id << " " << v2.id << std::endl;
    // std::cout << "f1 vids: ";
    for (auto vid: f1.Vids) {
        // std::cout << vid << " ";
        if (vid == v1.id || vid == v2.id) continue;
        if (std::find(v1.N_Vids.begin(), v1.N_Vids.end(), vid) != v1.N_Vids.end()) {
            v1_nVids.push_back(vid);
        }
        if (std::find(v2.N_Vids.begin(), v2.N_Vids.end(), vid) != v2.N_Vids.end()) {
            v2_nVids.push_back(vid);
        }
    }
    // std::cout << std::endl;
    // std::cout << "f2 vids: ";
    for (auto vid: f2.Vids) {
        // std::cout << vid << " ";
        if (vid == v1.id || vid == v2.id) continue;
        if (std::find(v1.N_Vids.begin(), v1.N_Vids.end(), vid) != v1.N_Vids.end()) {
            v1_nVids.push_back(vid);
        }
        if (std::find(v2.N_Vids.begin(), v2.N_Vids.end(), vid) != v2.N_Vids.end()) {
            v2_nVids.push_back(vid);
        }
    }
    // std::cout << std::endl;
    // std::cout << v1_nVids.size() << " " << v2_nVids.size() << std::endl;
    // std::cout << "In Edge Roate" << std::endl;
    Face new_f1;
    Face new_f2;
    Vertex& new_v1 = mesh.V.at(v1_nVids.at(0));
    new_f1.Vids.push_back(new_v1.id);
    new_f2.Vids.push_back(new_v1.id);
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

        value1 += getLength(new_v1.id, v2_nVids.at(0));
        value1 += getLength(v1.id, v2_nVids.at(0));
    } else {
        new_f1.Vids.push_back(v2_nVids.at(0));
        new_f1.Vids.push_back(v2.id);
        new_f1.Vids.push_back(v2_nVids.at(1));
        
        new_f2.Vids.push_back(v2_nVids.at(1));
        new_f2.Vids.push_back(v1_nVids.at(1));
        new_f2.Vids.push_back(v1.id);
        
        value1 += getLength(new_v1.id, v2_nVids.at(1));
        value1 += getLength(v1.id, v2_nVids.at(1));
    }
    // std::cout << "new f1 vids: " << new_f1.Vids.size() << std::endl;
    // for (auto vid: new_f1.Vids) {
    //     std::cout << " " << vid;
    // }
    // std::cout << std::endl;
    // std::cout << "new f2 vids: " << new_f2.Vids.size() << std::endl;
    // for (auto vid: new_f2.Vids) {
    //     std::cout << " " << vid;
    // }
    // std::cout << std::endl;
    Face new_f3;
    Face new_f4;
    Vertex& new_v2 = mesh.V.at(v1_nVids.at(1));
    new_f3.Vids.push_back(new_v2.id);
    new_f4.Vids.push_back(new_v2.id);
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
                
        value2 += getLength(new_v2.id, v2_nVids.at(0));
        value2 += getLength(v1.id, v2_nVids.at(0));
    } else {
        new_f3.Vids.push_back(v2_nVids.at(0));
        new_f3.Vids.push_back(v2.id);
        new_f3.Vids.push_back(v2_nVids.at(1));
        
        new_f4.Vids.push_back(v2_nVids.at(1));
        new_f4.Vids.push_back(v1_nVids.at(0));
        new_f4.Vids.push_back(v1.id);
                
        value2 += getLength(new_v2.id, v2_nVids.at(1));
        value2 += getLength(v1.id, v2_nVids.at(1));
    }
    // std::cout << "new f3 vids: " << new_f3.Vids.size() << std::endl;
    // for (auto vid: new_f3.Vids) {
    //     std::cout << " " << vid;
    // }
    // std::cout << std::endl;
    // std::cout << "new f4 vids: " << new_f4.Vids.size() << std::endl;
    // for (auto vid: new_f4.Vids) {
    //     std::cout << " " << vid;
    // }
    // std::cout << std::endl;

    double currentValue = 0;
    currentValue += getLength(v1.id, v2.id);
    currentValue += getLength(v1_nVids.at(0), v2.id);
    currentValue += getLength(v1_nVids.at(1), v2.id);
    currentValue += getLength(v2_nVids.at(0), v1.id);
    currentValue += getLength(v2_nVids.at(1), v1.id);

    if (currentValue - value1 > currentValue - value2) {
        op.newFaces.push_back(new_f1);
        op.newFaces.push_back(new_f2);
        op.profitability = currentValue - value1;
    } else {
        op.newFaces.push_back(new_f3);
        op.newFaces.push_back(new_f4);
        op.profitability = currentValue - value2;
    }
    // std::cout << "profitability: " << op.profitability << " new faces: " << op.newFaces.size() << std::endl;
    // std::cout << "*********************************************" << std::endl;
    op.canceledFids.insert(e.N_Fids.begin(), e.N_Fids.end());
    if (op.profitability > 0) { 
        Ops.insert(op);
    }
}

void LocalSimplifier::VertexRotate(Vertex& v, std::multiset<LocalOperation, bool(*)(LocalOperation, LocalOperation)>& Ops) {
    LocalOperation op;
    op.type = "Vertex Rotate";
    // std::set<size_t> canceledFids;
    double currentValue = 0;
    double newValue = 0;
    for (auto eid: v.N_Eids) {
        std::vector<size_t> newVids;
        Edge& e = mesh.E.at(eid);
        currentValue += getLength(e.Vids.at(0), e.Vids.at(1));
        // size_t id = e.Vids.at(0);
        // if (id == v.id) {
        //     id = e.Vids.at(1);
        // }
        newVids.push_back(e.Vids.at(0));
        Face& f1 = mesh.F.at(e.N_Fids.at(0));
        op.canceledFids.insert(f1.id);
        for (int i = 0; i < f1.Vids.size(); i++) {
            if (f1.Vids.at(i) == v.id) continue;
            if (std::find(v.N_Vids.begin(), v.N_Vids.end(), f1.Vids.at(i)) == v.N_Vids.end()) {
                newVids.push_back(f1.Vids.at(i));
                // if (std::find(e.Vids.begin(), e.Vids.end(), f1.Vids.at((i + 1) % f1.Vids.size())) != e.Vids.end() &&
                //     std::find(e.Vids.begin(), e.Vids.end(), f1.Vids.at((i + 2) % f1.Vids.size())) != e.Vids.end()) {
                //     newVids.push_back(f1.Vids.at((i + 1) % f1.Vids.size()));
                // } else if (std::find(e.Vids.begin(), e.Vids.end(), f1.Vids.at((i + 2) % f1.Vids.size())) != e.Vids.end() &&
                //            std::find(e.Vids.begin(), e.Vids.end(), f1.Vids.at((i + 3) % f1.Vids.size())) != e.Vids.end()) {
                //     newVids.push_back(f1.Vids.at((i + 2) % f1.Vids.size()));
                // }
                newValue += getLength(v.id, f1.Vids.at(i));
                // if (f1.Vids.at((i + 1) % f1.Vids.size()) == id) {
                //     newVids.insert(newVids.begin(), f1.Vids.at(i));
                // } else if (f1.Vids.at((i + 3) % f1.Vids.size()) == id) {
                //     newVids.insert(newVids.end(), f1.Vids.at(i));
                // }
                // newVids.push_back(f1.Vids.at(i));
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
                // if (std::find(e.Vids.begin(), e.Vids.end(), f2.Vids.at((i + 1) % f2.Vids.size())) != e.Vids.end() &&
                //     std::find(e.Vids.begin(), e.Vids.end(), f2.Vids.at((i + 2) % f2.Vids.size())) != e.Vids.end()) {
                //     newVids.push_back(f2.Vids.at((i + 1) % f2.Vids.size()));
                // } else if (std::find(e.Vids.begin(), e.Vids.end(), f2.Vids.at((i + 2) % f2.Vids.size())) != e.Vids.end() &&
                //            std::find(e.Vids.begin(), e.Vids.end(), f2.Vids.at((i + 3) % f2.Vids.size())) != e.Vids.end()) {
                //     newVids.push_back(f2.Vids.at((i + 2) % f2.Vids.size()));
                // }
                newValue += getLength(v.id, f2.Vids.at(i));
                // if (f2.Vids.at((i + 1) % f2.Vids.size()) == id) {
                //     newVids.insert(newVids.begin(), f2.Vids.at(i));
                // } else if (f2.Vids.at((i + 3) % f2.Vids.size()) == id) {
                //     newVids.insert(newVids.end(), f2.Vids.at(i));
                // }
                // newVids.push_back(f2.Vids.at(i));
                break;
            }
        }
        // newVids.push_back(v.id);
        Face newF;
        newF.Vids = newVids;
        // newF.id = mesh.F.size();
        op.newFaces.push_back(newF);
        // mesh.F.push_back(newF);
        if (newVids.size() < 4) {
            std::cout << "f1 vids " << f1.Vids.size() << ": ";
            for (auto vid: f1.Vids) {
                std::cout << vid << " ";
            }
            std::cout << std::endl;
            std::cout << "f2 vids " << f2.Vids.size() << ": ";
            for (auto vid: f2.Vids) {
                std::cout << vid << " ";
            }
            std::cout << std::endl;
            std::cout << "new vids " << newVids.size() << ": ";
            for (auto vid: newVids) {
                std::cout << vid << " ";
            }
            std::cout << std::endl;
            std::cout << "**********************************************************" << std::endl;
        }
    }
    op.profitability = currentValue - (newValue / 2);
    // op.canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
    // canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
    if (op.profitability > 0) {
        Ops.insert(op);
    }
}

/*void LocalSimplifier::VertexRotate(Vertex& v, std::vector<LocalOperation>& Ops) {
    LocalOperation op;
    op.type = "Vertex Rotate";
    // std::set<size_t> canceledFids;
    double currentValue = 0;
    double newValue = 0;
    std::set<size_t> star_v;
    star_v.insert(v.id);
    for (auto fid: v.N_Fids) {
        Face& f = mesh.F.at(fid);
        star_v.insert(f.Vids.begin(), f.Vids.end());
    }
    
    std::cout << "v: " << v.id << std::endl;
    std::cout << "v neighbors: " << v.N_Vids.size() << std::endl;
    std::cout << "star v: ";
    for (auto id: star_v) {
        std::cout << id << " ";
    }
    std::cout << std::endl;
    std::ofstream ofs("target_vids.vtk");
	ofs << "# vtk DataFile Version 3.0\n"
		<< "output.vtk\n"
		<< "ASCII\n\n"
		<< "DATASET UNSTRUCTURED_GRID\n";
	ofs << "POINTS " << mesh.V.size() << " double\n";
	
	for (size_t i = 0; i < mesh.V.size(); i++) {
		ofs << mesh.V.at(i).x << " " <<  mesh.V.at(i).y << " " <<  mesh.V.at(i).z << "\n";
	}
    std::vector<size_t> target_indices;
    target_indices.push_back(v.id);
    for (auto id: star_v) {
        target_indices.push_back(id);
    }
	ofs << "CELLS " << target_indices.size() << " " << 2 * target_indices.size() << std::endl;
	for (size_t i = 0; i < target_indices.size(); i++) {
		ofs << "1 " << target_indices.at(i) << std::endl;
	}
	ofs << "CELL_TYPES " << target_indices.size() << "\n";
	for (size_t i = 0; i < target_indices.size(); i++) {
		ofs << "1" << std::endl;
	}
    for (auto vid: v.N_Vids) {
        currentValue += getLength(v.id, vid);
        std::vector<size_t> newVids;
        newVids.push_back(v.id);
        Vertex& v_n = mesh.V.at(vid);
        std::vector<size_t> n_vs;
        std::cout << "neighbors: ";
        for (auto id: v_n.N_Vids) {
            std::cout << id << " ";
            if (id == v.id || std::find(star_v.begin(), star_v.end(), id) == star_v.end()) continue;
            n_vs.push_back(id);
            newValue += getLength(v.id, id);
        }
        std::cout << std::endl;
        std::cout << "n_vs size: " << n_vs.size() << std::endl;
        newVids.push_back(n_vs.at(0));
        newVids.push_back(v_n.id);
        newVids.push_back(n_vs.at(1));

        Face newF;
        newF.Vids = newVids;
        op.newFaces.push_back(newF);
    }
    for (auto eid: v.N_Eids) {
        std::vector<size_t> newVids;
        Edge& e = mesh.E.at(eid);
        currentValue += getLength(e.Vids.at(0), e.Vids.at(1));
        size_t id = e.Vids.at(0);
        if (id == v.id) {
            id = e.Vids.at(1);
        }
        newVids.push_back(v.id);
        Face& f1 = mesh.F.at(e.N_Fids.at(0));
        op.canceledFids.insert(f1.id);
        for (int i = 0; i < f1.Vids.size(); i++) {
            if (f1.Vids.at(i) == v.id) continue;
            if (std::find(v.N_Vids.begin(), v.N_Vids.end(), f1.Vids.at(i)) == v.N_Vids.end()) {
                newVids.push_back(f1.Vids.at(i));
                // if (std::find(e.Vids.begin(), e.Vids.end(), f1.Vids.at((i + 1) % f1.Vids.size())) != e.Vids.end() &&
                //     std::find(e.Vids.begin(), e.Vids.end(), f1.Vids.at((i + 2) % f1.Vids.size())) != e.Vids.end()) {
                //     newVids.push_back(f1.Vids.at((i + 1) % f1.Vids.size()));
                // } else if (std::find(e.Vids.begin(), e.Vids.end(), f1.Vids.at((i + 2) % f1.Vids.size())) != e.Vids.end() &&
                //            std::find(e.Vids.begin(), e.Vids.end(), f1.Vids.at((i + 3) % f1.Vids.size())) != e.Vids.end()) {
                //     newVids.push_back(f1.Vids.at((i + 2) % f1.Vids.size()));
                // }
                newValue += getLength(v.id, f1.Vids.at(i));
                // if (f1.Vids.at((i + 1) % f1.Vids.size()) == id) {
                //     newVids.insert(newVids.begin(), f1.Vids.at(i));
                // } else if (f1.Vids.at((i + 3) % f1.Vids.size()) == id) {
                //     newVids.insert(newVids.end(), f1.Vids.at(i));
                // }
                // newVids.push_back(f1.Vids.at(i));
                break;
            }
        }
        newVids.push_back(id);
        Face& f2 = mesh.F.at(e.N_Fids.at(1));
        op.canceledFids.insert(f2.id);
        for (int i = 0; i < f2.Vids.size(); i++) {
            if (f2.Vids.at(i) == v.id) continue;
            if (std::find(v.N_Vids.begin(), v.N_Vids.end(), f2.Vids.at(i)) == v.N_Vids.end()) {
                newVids.push_back(f2.Vids.at(i));
                // if (std::find(e.Vids.begin(), e.Vids.end(), f2.Vids.at((i + 1) % f2.Vids.size())) != e.Vids.end() &&
                //     std::find(e.Vids.begin(), e.Vids.end(), f2.Vids.at((i + 2) % f2.Vids.size())) != e.Vids.end()) {
                //     newVids.push_back(f2.Vids.at((i + 1) % f2.Vids.size()));
                // } else if (std::find(e.Vids.begin(), e.Vids.end(), f2.Vids.at((i + 2) % f2.Vids.size())) != e.Vids.end() &&
                //            std::find(e.Vids.begin(), e.Vids.end(), f2.Vids.at((i + 3) % f2.Vids.size())) != e.Vids.end()) {
                //     newVids.push_back(f2.Vids.at((i + 2) % f2.Vids.size()));
                // }
                newValue += getLength(v.id, f2.Vids.at(i));
                // if (f2.Vids.at((i + 1) % f2.Vids.size()) == id) {
                //     newVids.insert(newVids.begin(), f2.Vids.at(i));
                // } else if (f2.Vids.at((i + 3) % f2.Vids.size()) == id) {
                //     newVids.insert(newVids.end(), f2.Vids.at(i));
                // }
                // newVids.push_back(f2.Vids.at(i));
                break;
            }
        }
        // newVids.push_back(v.id);
        Face newF;
        newF.Vids = newVids;
        // newF.id = mesh.F.size();
        op.newFaces.push_back(newF);
        // mesh.F.push_back(newF);
        if (newVids.size() < 4) {
            std::cout << "f1 vids " << f1.Vids.size() << ": ";
            for (auto vid: f1.Vids) {
                std::cout << vid << " ";
            }
            std::cout << std::endl;
            std::cout << "f2 vids " << f2.Vids.size() << ": ";
            for (auto vid: f2.Vids) {
                std::cout << vid << " ";
            }
            std::cout << std::endl;
            std::cout << "new vids " << newVids.size() << ": ";
            for (auto vid: newVids) {
                std::cout << vid << " ";
            }
            std::cout << std::endl;
            std::cout << "**********************************************************" << std::endl;
        }
    }
    op.profitability = currentValue - (newValue / 2);
    op.canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
    // canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
    Ops.push_back(op);
}*/

void LocalSimplifier::EdgeCollapse(Edge& e, std::multiset<LocalOperation, bool(*)(LocalOperation, LocalOperation)>& Ops) {
    LocalOperation op;
    op.type = "Edge Collapse";
    op.profitability = getLength(e.Vids.at(0), e.Vids.at(1));
    // std::set<size_t> canceledFids;
    Vertex& v = mesh.V.at(e.Vids.at(0));
    op.canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
    /*std::set<size_t> star_v;
    star_v.insert(v.id);
    for (auto fid: v.N_Fids) {
        Face& f = mesh.F.at(fid);
        star_v.insert(f.Vids.begin(), f.Vids.end());
    }
    for (auto vid: v.N_Vids) {
        if (std::find(e.Vids.begin(), e.Vids.end(), v.id) != e.Vids.end() && 
            std::find(e.Vids.begin(), e.Vids.end(), vid) != e.Vids.end()) continue;
        std::vector<size_t> newVids;
        newVids.push_back(v.id);
        Vertex& v_n = mesh.V.at(vid);
        std::vector<size_t> n_vs;
        for (auto id: v_n.N_Vids) {
            if (id == v.id || std::find(star_v.begin(), star_v.end(), id) == star_v.end()) continue;
            n_vs.push_back(id);
        }
        newVids.push_back(n_vs.at(0));
        newVids.push_back(v_n.id);
        newVids.push_back(n_vs.at(1));

        Face newF;
        newF.Vids = newVids;
        op.newFaces.push_back(newF);
    }*/
    for (auto eid: v.N_Eids) {
        if (eid == e.id) continue;
        std::vector<size_t> newVids;
        Edge& e_n = mesh.E.at(eid);
        // size_t id = e_n.Vids.at(0);
        // if (id == v.id) {
        //     id = e_n.Vids.at(1);
        // }
        newVids.push_back(e_n.Vids.at(0));
        Face& f1 = mesh.F.at(e_n.N_Fids.at(0));
        for (int i = 0; i < f1.Vids.size(); i++) {
            if (f1.Vids.at(i) == v.id) continue;
            if (std::find(v.N_Vids.begin(), v.N_Vids.end(), f1.Vids.at(i)) == v.N_Vids.end()) {
                newVids.push_back(f1.Vids.at(i));
                // if (std::find(e_n.Vids.begin(), e_n.Vids.end(), f1.Vids.at((i + 1) % f1.Vids.size())) != e_n.Vids.end() &&
                //     std::find(e_n.Vids.begin(), e_n.Vids.end(), f1.Vids.at((i + 2) % f1.Vids.size())) != e_n.Vids.end()) {
                //     newVids.push_back(f1.Vids.at((i + 1) % f1.Vids.size()));
                // } else if (std::find(e_n.Vids.begin(), e_n.Vids.end(), f1.Vids.at((i + 2) % f1.Vids.size())) != e_n.Vids.end() &&
                //            std::find(e_n.Vids.begin(), e_n.Vids.end(), f1.Vids.at((i + 3) % f1.Vids.size())) != e_n.Vids.end()) {
                //     newVids.push_back(f1.Vids.at((i + 2) % f1.Vids.size()));
                // }
                // if (f1.Vids.at((i + 1) % f1.Vids.size()) == id) {
                //     newVids.insert(newVids.begin(), f1.Vids.at(i));
                // } else if (f1.Vids.at((i + 3) % f1.Vids.size()) == id) {
                //     newVids.insert(newVids.end(), f1.Vids.at(i));
                // }
                // newVids.push_back(f1.Vids.at(i));
                break;
            }
        }
        newVids.push_back(e_n.Vids.at(1));
        Face& f2 = mesh.F.at(e_n.N_Fids.at(1));
        for (int i = 0; i < f2.Vids.size(); i++) {
            if (f2.Vids.at(i) == v.id) continue;
            if (std::find(v.N_Vids.begin(), v.N_Vids.end(), f2.Vids.at(i)) == v.N_Vids.end()) {
                newVids.push_back(f2.Vids.at(i));
                // if (std::find(e_n.Vids.begin(), e_n.Vids.end(), f2.Vids.at((i + 1) % f2.Vids.size())) != e_n.Vids.end() &&
                //     std::find(e_n.Vids.begin(), e_n.Vids.end(), f2.Vids.at((i + 2) % f2.Vids.size())) != e_n.Vids.end()) {
                //     newVids.push_back(f2.Vids.at((i + 1) % f2.Vids.size()));
                // } else if (std::find(e_n.Vids.begin(), e_n.Vids.end(), f2.Vids.at((i + 2) % f2.Vids.size())) != e_n.Vids.end() &&
                //            std::find(e_n.Vids.begin(), e_n.Vids.end(), f2.Vids.at((i + 3) % f2.Vids.size())) != e_n.Vids.end()) {
                //     newVids.push_back(f2.Vids.at((i + 2) % f2.Vids.size()));
                // }
                // if (f2.Vids.at((i + 1) % f2.Vids.size()) == id) {
                //     newVids.insert(newVids.begin(), f2.Vids.at(i));
                // } else if (f2.Vids.at((i + 3) % f2.Vids.size()) == id) {
                //     newVids.insert(newVids.end(), f2.Vids.at(i));
                // }
                // newVids.push_back(f2.Vids.at(i));
                break;
            }
        }
        // newVids.push_back(v.id);
        Face newF;
        newF.Vids = newVids;
        // newF.id = mesh.F.size();
        // mesh.F.push_back(newF);
        op.newFaces.push_back(newF);
    }
    op.canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
    Vertex& v2 = mesh.V.at(e.Vids.at(1));
    op.canceledFids.insert(v2.N_Fids.begin(), v2.N_Fids.end());
    for (auto fid: v2.N_Fids) {
        if (std::find(v.N_Fids.begin(), v.N_Fids.end(), fid) != v.N_Fids.end()) continue;
        std::vector<size_t> newVids;
        for (auto new_vid: mesh.F.at(fid).Vids) {
            if (new_vid == v2.id) {
                newVids.push_back(v.id);
            } else {
                newVids.push_back(new_vid);
            }
        }
        Face newF;
        // newF.id = mesh.F.size();
        newF.Vids = newVids;
        // mesh.F.push_back(newF);
        op.newFaces.push_back(newF);
    }
    if (Ops.size() > 0) {
        std::multiset<LocalOperation>::iterator iter = Ops.begin();
        LocalOperation op1 = *iter;
        if (op.profitability < op1.profitability) {
            Ops.insert(op);
        }
    } else {
        Ops.insert(op);
    }
}

void LocalSimplifier::DiagonalCollapse(Face& f, std::multiset<LocalOperation, bool(*)(LocalOperation, LocalOperation)>& Ops) {
    LocalOperation op;
    op.type = "Diagonal Collapse";
    double value1 = getLength(f.Vids.at(0), f.Vids.at(2));
    double value2 = getLength(f.Vids.at(1), f.Vids.at(3));
    size_t target = f.Vids.at(0);
    size_t source = f.Vids.at(2);
    op.profitability = value1;
    if (value2 < value1) {
        target = f.Vids.at(1);
        source = f.Vids.at(3);
        op.profitability = value2;
    }
    for (auto fid: mesh.V.at(source).N_Fids) {
        if (fid == f.id) continue;
        Face& n_f = mesh.F.at(fid);
        Face newF;
        newF.Vids = n_f.Vids;
        for (int i = 0; i < n_f.Vids.size(); i++) {
            if (n_f.Vids.at(i) == source) {
                newF.Vids.at(i) = target;
                break;
            }
        }
        op.newFaces.push_back(newF);
        op.canceledFids.insert(n_f.id);
    }
    op.canceledFids.insert(f.id);
    
    if (Ops.size() > 0) {
        std::multiset<LocalOperation>::iterator iter = Ops.begin();
        LocalOperation op1 = *iter;
        if (op.profitability < op1.profitability) {
            Ops.insert(op);
        }
    } else {
        Ops.insert(op);
    }
}

void LocalSimplifier::RemoveDoublets(std::set<size_t>& canceledFids) {
    for (auto& v: mesh.V) {
        if (v.isBoundary) continue;
        if (v.N_Fids.size() == 2) {
            // std::cout << "Doublet identified" << std::endl;
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
            std::cout << "old fvids: ";
            for (auto vid: f2.Vids) {
                std::cout << vid << " ";
            }
            std::cout << std::endl;
            // std::cout << "diag vid: " << diag_vid << std::endl;
            std::vector<size_t> newVids;
            for (auto vid: f2.Vids) {
                if (vid == v.id) {
                    newVids.push_back(diag_vid);
                } else {
                    newVids.push_back(vid);
                }
            }
            f2.Vids = newVids;
            std::cout << newVids.size() << " ";
            std::cout << "new fvids: ";
            for (auto vid: f2.Vids) {
                std::cout << vid << " ";
            }
            std::cout << std::endl;
            for (int i = 0; i < f2.Vids.size() - 1; i++) {
                for (int j = i + 1; j < f2.Vids.size(); j++) {
                    if (f2.Vids.at(i) == f2.Vids.at(j)) {
                        std::cout << "ERROR: VERTICES REPEATED" << std::endl;
                    }
                }
            }
        }
    }
    std::cout << "# doublets removed: " << canceledFids.size() << std::endl;
}

double LocalSimplifier::getLength(size_t vid1, size_t vid2) {
    Vertex& v1 = mesh.V.at(vid1);
    Vertex& v2 = mesh.V.at(vid2);
    double length = glm::length(glm::dvec3(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z));
    return length;
}

