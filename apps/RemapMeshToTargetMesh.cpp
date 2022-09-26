#include "Mesh.h"
#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "QuadSurfaceMapper.h"
#include "MeshUtil.h"
#include "SemiGlobalSimplifier.h"
// #include "ParallelFor.h"
#include "FeatureExtractor.h"
#include "Smooth.h"

int main(int argc, char* argv[]) {

    // PQueue<std::string> pq;
    // std::vector<std::string> i = {"am ", "Hello ", "Naeem ", "World ", "I "};
    // std::vector<double> p = {0.6, 0.25, 0.7, 0.3, 0.4};
    // std::vector<int> s = {1, 2, 3, 4, 5};
    // // pq.setMaxQueueOn();
    // for (int j = 0; j < 5; j++) {
    //     pq.insert(p[j], s[j], i[j]);
    // }
    // // pq.update(2.5, 2);
    // // pq.update(1.8, 4);
    // // pq.update(1.5, 5);
    // // pq.update(0.1, 3);
    // // std::cout << pq.size() << std::endl;
    // while (!pq.empty()) {
    //     std::cout << pq.pop() << std::endl;
    // }
    // return 0;
    std::string source_f = argv[1];
    // std::string target_f = argv[2];
    std::string output_f = argv[2];
    int iters = atoi(argv[3]);

    MeshFileReader source_reader(source_f.c_str());
    Mesh& source = (Mesh&) source_reader.GetMesh();
    source.RemoveUselessVertices();
    source.BuildAllConnectivities();
    source.ExtractSingularities();
    // source.ExtractBoundary();
	source.BuildParallelE();

    // MeshFileReader target_reader(output_f.c_str());
    // Mesh& target = (Mesh&) target_reader.GetMesh();
    // target.RemoveUselessVertices();
    // target.BuildAllConnectivities();
    // target.ExtractSingularities();
    // source.ExtractBoundary();
	// target.BuildParallelE();

    // SurfaceMapper sm();
    // std::cout << "Haursdorff distance: " << sm.ExecuteHaursdorffDsitanceFilter(source, target)  << std::endl;
    // return;
    // for (auto& el: source.V) {
    //     std:cout << el.N_Fids.size() << std::endl;
    // }
    // return 0;
    // for (auto& v: source.V) {
    //     std::cout << v.id << ": " << " nvids: " << v.N_Vids.size() << " neids: " << v.N_Eids.size() << " nfids: " << v.N_Fids.size() << std::endl;
    // }
    

    MeshUtil mu(source);
    
    FeatureExtractor fe(source, 20.0, mu);
    fe.Extract();
    std::cout << "Writing output file" << std::endl;
    std::ofstream ofs("FeaturePoints.vtk");
    ofs << "# vtk DataFile Version 3.0\n"
        << "FeaturePoints.vtk.vtk\n"
        << "ASCII\n\n"
        << "DATASET UNSTRUCTURED_GRID\n";
    ofs << "POINTS " << source.V.size() << " double\n";
    std::vector<size_t> c_indices;
    // for (auto& v: source.V) {
        // if (v.type == FEATURE || v.isBoundary) c_indices.push_back(v.id);
        // if ((v.isBoundary) && v.N_Fids.size() != v.idealValence) {
        //     c_indices.push_back(v.id);
        //     std::cout << "vertices: " << v.N_Vids.size() << " edges: " << v.N_Eids.size() << " faces: "  << v.N_Fids.size() << " ideal valence: " << v.idealValence << std::endl;
        // } 
    // }
    // std::vector<size_t> c_indices = {12, 296};
    // std::cout << c_indices.size() << std::endl;
    std::set<size_t> eSet;
    for (auto& e: source.E) {
        auto& v1 = source.V.at(e.Vids.at(0));
        auto& v2 = source.V.at(e.Vids.at(1));
        if ((v1.type == FEATURE || v1.isBoundary) && (v2.type == FEATURE || v2.isBoundary)) eSet.insert(e.id);
    }
    c_indices.insert(c_indices.begin(), eSet.begin(), eSet.end());
    for (size_t i = 0; i < source.V.size(); i++) {
        ofs << std::fixed << std::setprecision(7) <<  source.V.at(i).x << " " <<  source.V.at(i).y << " " <<  source.V.at(i).z << "\n";
    }
    ofs << "CELLS " << c_indices.size() << " " << 3 * c_indices.size() << std::endl;
    for (size_t i = 0; i < c_indices.size(); i++) {
        auto& e = source.E.at(c_indices.at(i));
        ofs << "2 " << e.Vids.at(0) << " " << e.Vids.at(1) << std::endl;
        // ofs << "1 " << c_indices.at(i) << std::endl;
    }
    ofs << "CELL_TYPES " << c_indices.size() << "\n";
    for (size_t i = 0; i < c_indices.size(); i++) {
        ofs << "3" << std::endl;
    }
    return 0;
    // for (auto& v: source.V) {
    //     if ((v.isBoundary || v.type == FEATURE) && v.idealValence != 2) std::cout << v.idealValence << " ";
    // }
    // return 0;
    Smoother sm(source, mu);
    SemiGlobalSimplifier sg(source, mu, sm);
    sg.SetIters(iters);
    // sg.SetVertexSplitOperations();
    // sg.SetQuadSplitOperations();
    // sg.FixBoundary();
    // if (!source.isPlanar) {
        // sg.SetBoundaryDirectSeparatrixOperations(false);
        // sg.SetBoundaryDirectSeparatrixOperations(false);
    // }
    // std::cout << "Input mesh V: " << source.V.size() << std::endl; 
    // std::cout << "Input mesh F: " << source.F.size() << std::endl;
    // return 0;
    bool res = true;
    sg.SetBoundaryDirectSeparatrixOperations(true);
    sg.SetBoundaryDirectSeparatrixOperations(true);
    // res = true;
    while (res) {
        res = sg.FixBoundary();
        // res = sg.FixValences();
    }
    sg.SetDirectSeparatrixOperations(false);
    sg.SetDirectSeparatrixOperations(false);
    sg.SetDirectSeparatrixOperations(true);
    sg.SetDirectSeparatrixOperations(true);
    // sg.SetHalfSeparatrixOperations();
    // sg.SetHalfSeparatrixOperations();
    // sg.SetHalfSeparatrixOperations();
    // sg.SetHalfSeparatrixOperations();
    // sg.SetHalfSeparatrixOperations();
    // sg.SetHalfSeparatrixOperations();
    // sg.SetHalfSeparatrixOperations();
    res = true;
    while (res) {
        // res = sg.FixBoundary();
        res = sg.FixValences();
    }
    // res = true;
    // while (res) {
    // }

    int n = 0;
    std::vector<bool> isVisited(source.V.size(), false);
    for (auto& v: source.V) {
        if (!isVisited.at(v.id) && v.type != FEATURE && !v.isBoundary && (v.N_Vids.size() == 5 || v.N_Vids.size() == 3)) {
            int valenceToCheck = v.N_Vids.size() == 5 ? 3 : 5;
            for (auto fid: v.N_Fids) {
                auto& f = source.F.at(fid);
                int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
                auto& fv = source.V.at(f.Vids.at((idx+2)%f.Vids.size()));
                if (!isVisited.at(fv.id) && fv.type != FEATURE && !fv.isBoundary && fv.N_Vids.size() == valenceToCheck) {
                    n += 1;
                    isVisited.at(v.id) = true;
                    isVisited.at(fv.id) = true;
                    break;
                }
            }
        }
    }
    std::cout << "Diagonal 3-5 pairs: " << n << std::endl;
    sg.Smooth();
    sg.PrototypeE();
    // sg.PrototypeD();
    // sg.PrototypeB();
    // sg.PrototypeB();
    // sg.PrototypeB();
    // sg.PrototypeB();
    // sg.Smooth();

    // sg.PrototypeBoundary();
    // sg.AlignAndResolveSingularities();
    // sg.PrototypeBoundary();
    // sg.AlignAndResolveSingularities();
    // sg.PrototypeBoundary();
    // sg.AlignAndResolveSingularities();
    // sg.PrototypeBoundary();
    // sg.AlignAndResolveSingularities();
    
    // sg.PrototypeBoundary(false);
    // sg.AlignAndResolveSingularities(false);
    // sg.PrototypeBoundary(false);
    // sg.AlignAndResolveSingularities(false);
    // sg.PrototypeBoundary(false);
    // sg.AlignAndResolveSingularities(false);
    // sg.PrototypeBoundary(false);
    // sg.AlignAndResolveSingularities(false);
    
    // sg.ResolveSingularities();
    
    // sg.GetSingularityPairs();
    // sg.AlignSingularities();
    // res = true;
    // while (res) {
    //     res = sg.FixBoundary();
    //     res = sg.ResolveHighValences();
    // }
    // res = true;
    // while (res) {
    //     res = sg.ResolveHighValences();
    //     sg.FixBoundary();
    // }
    // sg.Smooth();
    // for (int i = 0; i < 1; i++) {
    //     sg.ResolveSingularities();
    //     res = true;
    //     while (res) {
    //         res = sg.ResolveHighValences();
    //     }
    //     // sg.Smooth();
    // }
    
    // std::cout << "After GetSingularityPairs" << std::endl;
    // sg.SetQuadSplitOperations();
    // sg.SetDirectSeparatrixOperations();
    // sg.SetBoundarySeparatrixOperations();
    // sg.FixBoundary();
    // sg.SetBoundarySeparatrixOperations();
    // sg.FixBoundary();
    // sg.FixBoundary();
    // for (auto& v: source.V) {
    //     if (v.N_Vids.size() == 4 && v.isSingularity) {
    //         std::cout << "Not a singularity: " << v.id << std::endl;
    //     }
    //     // if (v.N_Vids.size() != v.N_Eids.size() || v.N_Vids.size() != v.N_Fids.size() || v.N_Eids.size() != v.N_Fids.size()) {
    //     //     std::cout << v.id << ": " << v.N_Vids.size() << " " << v.N_Eids.size() << " " << v.N_Fids.size() << std::endl;
    //     // }
    // }
    // sg.PerformGlobalOperations(); 
    // std::cout << "Done with Direct Separatrix Operations" << std::endl; 
    // for (int i = 0; i < 10; i++) {
    //     sg.SetHalfSeparatrixOperations(); 
    // }
    // for (int i = 0; i < 10; i++) {
    //     sg.SetSeparatrixOperations(); 
    // }
    // for (int i = 0; i < 10; i++) {
    //     sg.SetChordCollapseOperations(); 
    // }
    // sg.FixBoundary();
    // sg.SetHalfSeparatrixOperations();
    // sg.SetChordCollapseOperations();
    // std::cout << "Before Smoothing" << std::endl;
    // sg.SetEdgeRotationOperations();
    // sg.SetVertexRotationOperations();
    // sg.SetDiagonalCollapseOperations();
    // sg.SetEdgeCollapseOperations();
    // sg.SetSimplificationOperations();
    // sg.SetDiagonalCollapseOperations();
    // for (auto& v: source.V) {
    //     // vu.GetVertexEnergy(v.id);
    //     std::cout << mu.GetVertexEnergy(v.id) << std::endl;
    // }

    for (auto& f: source.F) {
        if (f.Vids.empty()) continue;
        for (auto eid: f.Eids) {
            auto& e = source.E.at(eid);
            if (e.Vids.size() != 2) {
                std::cout << f.id << " edges: " << f.Eids.size() << std::endl;
                std::cout << e.id << " vertices: " << e.Vids.size() << std::endl;
            }
        }
    }
    for (auto& e: source.E) {
        if (e.Vids.empty()) continue;
        for (auto fid: e.N_Fids) {
            auto& f = source.F.at(fid);
            if (f.Vids.empty()) {
                std::cout << "edge face is empty" << std::endl;
            }
            for (auto feid: f.Eids) {
                auto& fe = source.E.at(feid);
                if (fe.Vids.size() != 2) {
                    std::cout << f.id << " edges: " << f.Eids.size() << std::endl;
                    std::cout << fe.id << " vertices: " << fe.Vids.size() << std::endl;
                }
            }
        }
    }
    for (auto& v: source.V) {
        for (auto eid: v.N_Eids) {
            auto& e = source.E.at(eid);
            // if (e.N_Fids.size() != 2) std::cout << "edge faces are not 2" << std::endl;
            for (auto fid: e.N_Fids) {
                auto& f = source.F.at(fid);
                for (auto feid: f.Eids) {
                    auto& fe = source.E.at(feid);
                    if (fe.Vids.size() != 2) {
                        std::cout << f.id << " edges: " << f.Eids.size() << std::endl;
                        std::cout << fe.id << " vertices: " << fe.Vids.size() << std::endl;
                    }
                }
            }    
        }
    }
    
    // for (auto& v: source.V) {
    //     if ((v.type == FEATURE || v.isBoundary) && v.idealValence != 2) {
    //         std::cout << v.idealValence << " ";
    //     }
    // }
    // std::cout << std::endl;
    

    // MeshFileReader target_reader(target_f.c_str());
    // Mesh& target = (Mesh&) target_reader.GetMesh();
    // target.RemoveUselessVertices();
    // target.BuildAllConnectivities();

    // SurfaceMapper sm(target);
    // sm.Map();

    // source.RemoveUselessVertices();
    std::cout << "# F in input mesh: " << source.C.size() << std::endl;
    // std::cout << source.F.size() << std::endl;
    std::vector<Cell> newC;
    for (auto& f: source.F) {
        if (f.N_Fids.size() == 0) continue;
        Cell c;
        c.id = newC.size();
        c.Vids = f.Vids;
        newC.push_back(c);
    }
    source.E.clear();
    source.F.clear();
    source.C.clear();
    source.C.insert(source.C.begin(), newC.begin(), newC.end());
    // source.BuildAllConnectivities();
    // std::cout << newC.size() << std::endl;
    std::cout <<  "# F in output mesh: " << source.C.size() << std::endl;

    // std::cout << "Writing output file" << std::endl;
    // std::ofstream ofs("FeaturePoints.vtk");
    // ofs << "# vtk DataFile Version 3.0\n"
    //     << "FeaturePoints.vtk.vtk\n"
    //     << "ASCII\n\n"
    //     << "DATASET UNSTRUCTURED_GRID\n";
    // ofs << "POINTS " << source.V.size() << " double\n";
    // std::vector<size_t> c_indices;
    // for (auto& v: source.V) {
    //     // if (v.type == FEATURE) c_indices.push_back(v.id);
    //     if (v.isBoundary) c_indices.push_back(v.id);
    // }
    // // std::vector<size_t> c_indices = {12, 296};
    // // std::cout << c_indices.size() << std::endl;
    // for (size_t i = 0; i < source.V.size(); i++) {
    //     ofs << std::fixed << std::setprecision(7) <<  source.V.at(i).x << " " <<  source.V.at(i).y << " " <<  source.V.at(i).z << "\n";
    // }
    // ofs << "CELLS " << c_indices.size() << " " << 2 * c_indices.size() << std::endl;
    // for (size_t i = 0; i < c_indices.size(); i++) {
    //     ofs << "1 " << c_indices.at(i) << std::endl;
    // }
    // ofs << "CELL_TYPES " << c_indices.size() << "\n";
    // for (size_t i = 0; i < c_indices.size(); i++) {
    //     ofs << "1" << std::endl;
    // }

    MeshFileWriter writer(source, output_f.c_str());
    writer.WriteFile();

    return 0;
}


/* Edge rotation Validation

std::cout << e.Vids.size() << " " << e.N_Fids.size() << std::endl;
for (auto fid: e.N_Fids) {
    auto& f = mesh.F.at(fid);
    for (auto feid: f.Eids) {
        auto& fe = mesh.E.at(feid);
        if (fe.Vids.size() != 2) {
            std::cout << f.id << " edges: " << f.Eids.size() << std::endl;
            std::cout << fe.id << " vertices: " << fe.Vids.size() << std::endl;
            breakLoop = true;
            break;
        }
    }
    if (breakLoop) break;    
}
break;
}
if (breakLoop) break;*/

/* Direct Separatrix Validation
// bool breakLoop = false;
auto& centerV = mesh.V.at(op->GetCenterId());
std::cout << "validating operation" << std::endl;
for (auto fid: centerV.N_Fids) {
    auto& f = mesh.F.at(fid);
    for (auto eid: f.Eids) {
        auto& e = mesh.E.at(eid);
        if (e.Vids.size() != 2) {
            std::cout << f.id << " edges: " << f.Eids.size() << std::endl;
            std::cout << e.id << " vertices: " << e.Vids.size() << std::endl;
            breakLoop = true;
        }
    }
}
for (auto vid: centerV.N_Vids) {
    auto& v = mesh.V.at(vid);
    for (auto fid: v.N_Fids) {
        auto& f = mesh.F.at(fid);
        for (auto eid: f.Eids) {
            auto& e = mesh.E.at(eid);
            if (e.Vids.size() != 2) {
                std::cout << f.id << " edges: " << f.Eids.size() << std::endl;
                std::cout << e.id << " vertices: " << e.Vids.size() << std::endl;
                breakLoop = true;
            }
        }
    }    
}
if (breakLoop) break;*/

// Edge Rotation cases:
// 1. Boundary conformity
// 2. To cennect 3-5 singularities directly
// 3. Edge length update
// 4. Rotate and analyze if it improves base complex configuration

// Diagonal Collapse cases:
// 1. Get out of local minimum or maxima
// 2. Create new singularity to enable new operations

// Vertex Rotation:
// Edge Collapse:

// Proper metric
// Comparison
// Ranking based on configuration

// 3-singularity movements
// Same rotation: element reduction all the way
// Different rotation: element reduction along a and element addition along b
// Same rotation both links: invalid
// Different rotation both links: invalid
// Different rotation one link: a must be greater than other link's a

// 5-singularity movements
