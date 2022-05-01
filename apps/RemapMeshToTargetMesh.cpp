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
    source.ExtractBoundary();
	source.BuildParallelE();	
    // for (auto& el: source.V) {
    //     std:cout << el.N_Fids.size() << std::endl;
    // }
    // return 0;
    FeatureExtractor fe(source, 20.0);
    fe.Extract();


    MeshUtil mu(source);
    Smoother sm(source);
    SemiGlobalSimplifier sg(source, mu, sm);
    sg.SetIters(iters);
    // sg.SetVertexSplitOperations();
    // sg.SetQuadSplitOperations();
    // sg.FixBoundary();
    sg.SetBoundaryDirectSeparatrixOperations(false);
    sg.SetBoundaryDirectSeparatrixOperations(true);
    sg.FixBoundary();
    sg.SetDirectSeparatrixOperations(false);
    sg.SetDirectSeparatrixOperations(true);
    sg.FixBoundary();
    sg.GetSingularityPairs();
    // sg.SetQuadSplitOperations();
    // sg.SetDirectSeparatrixOperations();
    // sg.SetBoundarySeparatrixOperations();
    // sg.FixBoundary();
    // sg.SetBoundarySeparatrixOperations();
    // sg.FixBoundary();
    // sg.ResolveSingularityPairs();
    // sg.FixBoundary();
    for (auto& v: source.V) {
        if (v.N_Vids.size() == 4 && v.isSingularity) {
            std::cout << "Not a singularity: " << v.id << std::endl;
        }
        // if (v.N_Vids.size() != v.N_Eids.size() || v.N_Vids.size() != v.N_Fids.size() || v.N_Eids.size() != v.N_Fids.size()) {
        //     std::cout << v.id << ": " << v.N_Vids.size() << " " << v.N_Eids.size() << " " << v.N_Fids.size() << std::endl;
        // }
    }
    // sg.PerformGlobalOperations(); 
    // std::cout << "Done with Direct Separatrix Operations" << std::endl; 
    // sg.SetSeparatrixOperations(); 
    // sg.FixBoundary();
    // sg.SetHalfSeparatrixOperations();
    // sg.SetChordCollapseOperations();
    sg.Smooth();
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
            if (e.N_Fids.size() != 2) std::cout << "edge faces are not 2" << std::endl;
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
    //     if (v.type == FEATURE) c_indices.push_back(v.id);
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