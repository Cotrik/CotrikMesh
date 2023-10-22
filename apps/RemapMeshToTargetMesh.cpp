#include "Mesh.h"
#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "QuadSurfaceMapper.h"
#include "MeshUtil.h"
#include "SemiGlobalSimplifier.h"
// #include "ParallelFor.h"
#include "FeatureExtractor.h"
#include "Smooth.h"
#include "KDTree.h"
#include <ctime>

int Vertex::visitCounter = 0;

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
    mu.Q_A = source.F.size();
    mu.delta = mu.Q_A - source.F.size();
    FeatureExtractor fe(source, 30.0, mu);
    fe.Extract();
    std::cout << "Writing output file" << std::endl;
    std::ofstream ofs("FeaturePoints.vtk");
    ofs << "# vtk DataFile Version 3.0\n"
        << "FeaturePoints.vtk.vtk\n"
        << "ASCII\n\n"
        << "DATASET UNSTRUCTURED_GRID\n";
    ofs << "POINTS " << source.V.size() << " double\n";
    std::vector<size_t> c_indices;
    for (auto& v: source.V) {
        if (v.type == FEATURE) c_indices.push_back(v.id);
        if (v.isBoundary) {
            c_indices.push_back(v.id);
            // std::cout << "v fids: " << v.N_Fids.size() << std::endl;
            // for (auto fid: v.N_Fids) {
            //     auto& f = source.C.at(fid);
            //     for (auto fvid: f.Vids) {
            //         if (fvid == v.id) continue;
            //         auto& fv = source.V.at(fvid);
            //         if (fv.isBoundary) std::cout << fvid << " ";
            //     }
            //     std::cout << std::endl;
            // }
            // std::cout << std::endl;
        }
    }
    std::cout << "Got features" << std::endl;
    // std::vector<size_t> c_indices = {12, 296};
    // std::cout << c_indices.size() << std::endl;
    for (size_t i = 0; i < source.V.size(); i++) {
        ofs << std::fixed << std::setprecision(7) <<  source.V.at(i).x << " " <<  source.V.at(i).y << " " <<  source.V.at(i).z << "\n";
    }
    ofs << "CELLS " << c_indices.size() << " " << 2 * c_indices.size() << std::endl;
    for (size_t i = 0; i < c_indices.size(); i++) {
        ofs << "1 " << c_indices.at(i) << std::endl;
    }
    ofs << "CELL_TYPES " << c_indices.size() << "\n";
    for (size_t i = 0; i < c_indices.size(); i++) {
        ofs << "1" << std::endl;
    }
    Smoother s(source, mu);
    // std::unique_ptr<SurfaceProjector> sp = std::make_unique<SurfaceProjector>(source);
    // SurfaceProjector sp(source);
    KDTree kd(source);
    // Smoother sm(source, mu);
    SemiGlobalSimplifier sg(source, mu, s, kd);
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
    std::clock_t start;
    double duration;
    double startSingularities;
    double nSingularities;
    
    startSingularities = 0.0;
    for (auto& v: source.V) {
        if (v.N_Vids.empty() || v.N_Eids.empty() || v.N_Fids.empty() || v.type == FEATURE || v.isBoundary) continue;
        if (v.N_Vids.size() != 4) startSingularities += 1;
    }

    bool res = true;
    start = std::clock();
    // sg.SetVertexRotationOperations();
    if (iters == 0) {
        while (sg.SetBoundaryDirectSeparatrixOperations(false));
        while (sg.SetBoundaryDirectSeparatrixOperations(true));
        while (sg.SetDirectSeparatrixOperations(false));
        while (sg.SetDirectSeparatrixOperations(true));
        sg.Smooth(nullptr);
    } else {
        while (sg.SetBoundaryDirectSeparatrixOperations(false));
        while (sg.SetBoundaryDirectSeparatrixOperations(true));
        while (sg.SetDirectSeparatrixOperations(false));
        while (sg.SetDirectSeparatrixOperations(true));
        for (int i = 0; i < iters; i++) {
        // for (int i = 0; i < 1; i++) {
            // while(sg.FixValences());
            // while(sg.FixValences());
            // while(sg.FixValences());
            // while(sg.FixValences());
            if (!sg.TestFlips()) break;
            // while(sg.FixValences());
        }
    }
    
    // sg.Smooth();
    duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
    std::cout << "Simplification time: " << duration << " seconds" << std::endl;
    // sg.SetBoundaryDirectSeparatrixOperations(false);
    // sg.SetBoundaryDirectSeparatrixOperations(false);
    // sg.SetBoundaryDirectSeparatrixOperations(true);
    // sg.SetBoundaryDirectSeparatrixOperations(true);
    
    // res = true;
    // while (res) {
    //     res = sg.RemoveDoublets();
        // res = sg.FixBoundary();
        // res = sg.FixValences();
    // }
    // for (int i = 0; i < iters; i++) {
    //     bool res = true;
    //     while (res) {
    //         res = sg.RemoveDoublets();
    //     }
    //     res = true;
    //     while (res) {
    //         res = sg.FixBoundary();
    //     }
    //     // sg.SetSeparatrixOperations();
    //     // sg.SetChordCollapseOperations();
    //     sg.SetDiagonalCollapseOperations();
    // }

    // sg.SetDirectSeparatrixOperations(false);
    // sg.SetDirectSeparatrixOperations(false);
    // sg.SetDirectSeparatrixOperations(true);
    // sg.SetDirectSeparatrixOperations(true);
    
    // sg.SetHalfSeparatrixOperations();
    // sg.SetHalfSeparatrixOperations();
    // sg.SetHalfSeparatrixOperations();
    // sg.SetHalfSeparatrixOperations();
    // sg.SetHalfSeparatrixOperations();
    // sg.SetHalfSeparatrixOperations();
    // sg.SetHalfSeparatrixOperations();
    // while (sg.FixValences());
    // res = true;
    
    // duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
    // std::cout << "Direct Separatrix simplification time: " << duration << " seconds" << std::endl;
    
    nSingularities = 0.0;
    for (auto& v: source.V) {
        if (v.N_Vids.empty() || v.N_Eids.empty() || v.N_Fids.empty() || v.type == FEATURE || v.isBoundary) continue;
        if (v.N_Vids.size() != 4) nSingularities += 1;
    }
    std::cout << "Initial Singularities: " << startSingularities << " Reduced Singularities: " << startSingularities - nSingularities << " ratio: " << ((startSingularities - nSingularities) / startSingularities) * 100.0 << "%" << std::endl;
    // res = true;
    // while (res) {
    // }
    // startSingularities = 0.0;
    // for (auto& v: source.V) {
    //     if (v.N_Vids.empty() || v.N_Eids.empty() || v.N_Fids.empty() || v.type == FEATURE || v.isBoundary) continue;
    //     if (v.N_Vids.size() != 4) startSingularities += 1;
    // }
        // start = std::clock();
        // sg.TestFlips();

    // duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
    // std::cout << "Paths time: " << duration << " seconds" << std::endl;
    // int n = 0;
    // std::vector<bool> isVisited(source.V.size(), false);
    // for (auto& v: source.V) {
    //     if (!isVisited.at(v.id) && v.type != FEATURE && !v.isBoundary && (v.N_Vids.size() == 5 || v.N_Vids.size() == 3)) {
    //         int valenceToCheck = v.N_Vids.size() == 5 ? 3 : 5;
    //         for (auto fid: v.N_Fids) {
    //             auto& f = source.F.at(fid);
    //             int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
    //             auto& fv = source.V.at(f.Vids.at((idx+2)%f.Vids.size()));
    //             if (!isVisited.at(fv.id) && fv.type != FEATURE && !fv.isBoundary && fv.N_Vids.size() == valenceToCheck) {
    //                 n += 1;
    //                 isVisited.at(v.id) = true;
    //                 isVisited.at(fv.id) = true;
    //                 break;
    //             }
    //         }
    //     }
    // }
    // std::cout << "Diagonal 3-5 pairs: " << n << std::endl;
    // sg.Smooth();
    
    // for (int i = 0; i < 10; i++) {
    //     sg.PrototypeK();
    //     sg.FixValences();
    // }
    // for (int i = 0; i < iters; i++) {
    // for (int i = 0; i < 1; i++) {
        // sg.PrototypeExecute();
        
        // while (sg.FixValences());
        // sg.func1();
        // sg.Smooth();
    // }
    // sg.RenderMesh(); 
    // sg.Smooth();
    // sg.PrototypeH(1);
    // sg.PrototypeH();
    // sg.PrototypeH(1);
    // sg.PrototypeH();
    // sg.PrototypeH(1);
    // sg.PrototypeH();
    // sg.PrototypeH(1);
    // sg.PrototypeH();
    // sg.PrototypeH(1);
    // while (sg.PrototypeI());
    // sg.PrototypeJ();
    // sg.PrototypeF();
    // if (!sg.CheckMeshValidity()) return 0;
    // sg.PrototypeF();
    // if (!sg.CheckMeshValidity()) return 0;
    // sg.CheckMeshValidity();
    // sg.PrototypeF();
    // if (!sg.CheckMeshValidity()) return 0;
    // sg.PrototypeF(1);
    // if (!sg.CheckMeshValidity()) return 0;
    // sg.PrototypeF(1);
    // if (!sg.CheckMeshValidity()) return 0;
    // sg.PrototypeF(1);
    // if (!sg.CheckMeshValidity()) return 0;
    // sg.PrototypeF(3);
    // if (!sg.CheckMeshValidity()) return 0;
    // sg.PrototypeF(3);
    // if (!sg.CheckMeshValidity()) return 0;
    // sg.PrototypeF(3);
    // duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
    // std::cout << "Direct pairs simplification time: " << duration << " seconds" << std::endl;
    
    // nSingularities = 0.0;
    // for (auto& v: source.V) {
    //     if (v.N_Vids.empty() || v.N_Eids.empty() || v.N_Fids.empty() || v.type == FEATURE || v.isBoundary) continue;
    //     if (v.N_Vids.size() != 4) nSingularities += 1;
    // }
    // std::cout << "Initial Singularities: " << startSingularities << " Reduced Singularities: " << startSingularities - nSingularities << " ratio: " << ((startSingularities - nSingularities) / startSingularities) * 100.0 << "%" << std::endl;
    // sg.SaveEdgeMesh();
    // sg.Smooth();
    // sg.PrototypeD();
    // sg.PrototypeC();
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
    // sg.Smooth();
    
    // sg.GetSingularityPairs();
    // sg.AlignSingularities();
    // startSingularities = 0.0;
    // for (auto& v: source.V) {
    //     if (v.N_Vids.empty() || v.N_Eids.empty() || v.N_Fids.empty() || v.type == FEATURE || v.isBoundary) continue;
    //     if (v.N_Vids.size() != 4) startSingularities += 1;
    // }
    // start = std::clock();
    
    // res = true;
    // while (res) {
    //     res = sg.FixValences();
    // }
    // // sg.Smooth();
    // for (int i = 0; i < 1; i++) {
    //     sg.ResolveSingularities();
    //     res = true;
    //     while (res) {
    //         res = sg.FixValences();
    //     }
    //     sg.Smooth();
    // }
    // duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
    // std::cout << "singulairty simplification time: " << duration << " seconds" << std::endl;
    
    // nSingularities = 0.0;
    // for (auto& v: source.V) {
    //     if (v.N_Vids.empty() || v.N_Eids.empty() || v.N_Fids.empty() || v.type == FEATURE || v.isBoundary) continue;
    //     if (v.N_Vids.size() != 4) nSingularities += 1;
    // }
    // std::cout << "Initial Singularities: " << startSingularities << " Reduced Singularities: " << startSingularities - nSingularities << " ratio: " << ((startSingularities - nSingularities) / startSingularities) * 100.0 << "%" << std::endl;
    
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

    // for (auto& f: source.F) {
    //     if (f.Vids.empty() || f.N_Fids.empty()) continue;
    //     for (auto eid: f.Eids) {
    //         auto& e = source.E.at(eid);
    //         if (e.Vids.size() != 2) {
    //             std::cout << f.id << " edges: " << f.Eids.size() << std::endl;
    //             std::cout << e.id << " vertices: " << e.Vids.size() << std::endl;
    //         }
    //     }
    // }
    // for (auto& e: source.E) {
    //     if (e.Vids.empty()) continue;
    //     for (auto fid: e.N_Fids) {
    //         auto& f = source.F.at(fid);
    //         if (f.Vids.empty()) {
    //             std::cout << "edge face is empty" << std::endl;
    //         }
    //         for (auto feid: f.Eids) {
    //             auto& fe = source.E.at(feid);
    //             if (fe.Vids.size() != 2) {
    //                 std::cout << f.id << " edges: " << f.Eids.size() << std::endl;
    //                 std::cout << fe.id << " vertices: " << fe.Vids.size() << std::endl;
    //             }
    //         }
    //     }
    // }
    // for (auto& v: source.V) {
    //     for (auto eid: v.N_Eids) {
    //         auto& e = source.E.at(eid);
    //         // if (e.N_Fids.size() != 2) std::cout << "edge faces are not 2" << std::endl;
    //         for (auto fid: e.N_Fids) {
    //             auto& f = source.F.at(fid);
    //             for (auto feid: f.Eids) {
    //                 auto& fe = source.E.at(feid);
    //                 if (fe.Vids.size() != 2) {
    //                     std::cout << f.id << " edges: " << f.Eids.size() << std::endl;
    //                     std::cout << fe.id << " vertices: " << fe.Vids.size() << std::endl;
    //                 }
    //             }
    //         }    
    //     }
    // }
    
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
    // for (int i = 0; i < iters; i++) {
    // }
    std::cout << "# F in input mesh: " << source.C.size() << std::endl;
    // std::cout << source.F.size() << std::endl;
    std::vector<Cell> newC;
    for (auto& f: source.F) {
        if (f.Vids.empty()) continue;
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
    //     if (v.type == FEATURE) c_indices.push_back(v.id);
    //     if (v.isBoundary) {
    //         c_indices.push_back(v.id);
    //         // std::cout << "v fids: " << v.N_Fids.size() << std::endl;
    //         // for (auto fid: v.N_Fids) {
    //         //     auto& f = source.C.at(fid);
    //         //     for (auto fvid: f.Vids) {
    //         //         if (fvid == v.id) continue;
    //         //         auto& fv = source.V.at(fvid);
    //         //         if (fv.isBoundary) std::cout << fvid << " ";
    //         //     }
    //         //     std::cout << std::endl;
    //         // }
    //         // std::cout << std::endl;
    //     }
    // }
    // std::cout << "Got features" << std::endl;
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

// Different rotations both: invalid
// Same rotations both: invalid
// One same one different: valid (if same's a > different's a)

// 5-singularity movements
// Same rotations: invalid


// why quad and hexes mesh over tri and tet mesh
// why start with tri or tet mesh
// include qudriflow reference
// 2021 paper comparison

// list bullet points for comparison in paper
// comparison with different triangle inputs
// comparison with other methods

// rot 5 edge large b1 - small b2 (elements removed) subtract 1 otherwise large b2 - small b1 (elements added)
// nrot 5 edge large b2 - small b1 (elements removed) subtract 1 otherwise large b1 - small b2 (elements added)

// direct three five pair:
// threeEdges rot > 1, rot > 1; rot == 1, rot == 1
// five rot edges rot >, rot > 1; rot == 1, rot == 1
// five n rot edges rot == 2, rot == 1; rot > 2, rot > 1