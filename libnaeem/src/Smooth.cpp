#include "Smooth.h"
#include "ParallelFor.h"

Smoother::Smoother() {}
Smoother::Smoother(Mesh& mesh_) : mesh(mesh_), sm(mesh) {
    // sm.SetTarget(mesh);
}

Smoother::Smoother(const Smoother& r) {
    mesh = r.mesh;
}

Smoother::~Smoother() {}

void Smoother::CheckValidity() {
    if (mesh.V.size() == 0 || mesh.F.size() == 0 || mesh.C.size() == 0) {
        std::cout << "No mesh to use for MeshUtils." << std::endl;
        exit(0);
    }
}

void Smoother::SetMesh(Mesh& mesh_) {
    mesh = mesh_;
    sm.SetTarget(mesh);
}

void Smoother::Smooth(std::vector<size_t>& V) {
    CheckValidity();

    if (V.empty()) {
        V.resize(mesh.V.size());
        PARALLEL_FOR_BEGIN(mesh.V.size()) {
            GetVerticesToSmooth(i, V);
        } PARALLEL_FOR_END();
        // std::cout << "Setting Vertices to smooth" << std::endl;
        // for (int i = 0; i < mesh.V.size(); i++) {
        //     GetVerticesToSmooth(i, V);
        // }
    }
    // smooth and project
	int iters = 10;

    // SurfaceMapper sm(mesh, origMesh);
    std::vector<glm::dvec3> centers(V.size());
	while (iters--) {
        centers.clear();
        centers.resize(V.size());
        PARALLEL_FOR_BEGIN(V.size()) {
            GetOptimizedPositions(i, V, centers);
        } PARALLEL_FOR_END();
        // std::cout << "Getting optimized positions" << std::endl;
        // for (int i = 0; i < mesh.V.size(); i++) {
        //     GetOptimizedPositions(i, V, centers);
        // }
        PARALLEL_FOR_BEGIN(V.size()) {
            mesh.V.at(V.at(i)) = centers.at(i);
        } PARALLEL_FOR_END();
        // std::cout << "Setting centers" << std::endl;
        // for (int i = 0; i < mesh.V.size(); i++) {
        //     mesh.V.at(V.at(i)) = centers.at(i);
        // }
    }
    for (auto& v: mesh.V) {
        v = sm.GetClosestPoint(v.xyz());
    }
}

void Smoother::GetVerticesToSmooth(int iter, std::vector<size_t>& V) {
    V.at(iter) = mesh.V.at(iter).id;
}

void Smoother::GetOptimizedPositions(int iter, std::vector<size_t>& V, std::vector<glm::dvec3>& centers) {    
    // std::cout << "iter: " << iter << " " << V.at(iter) << std::endl;
    auto& v = mesh.V.at(V.at(iter));
    // std::cout << "v: " << v.id << std::endl;
    centers.at(iter) = v.xyz();
    if (v.N_Vids.empty() || v.N_Eids.empty() || v.N_Fids.empty()) return;
    if (v.type < FEATURE) {
        double n = 0.0;
        glm::dvec3 center(0.0, 0.0, 0.0);
        // std::cout << "vertices neighbor edges access " << v.id << std::endl;
        for (int i = 0; i < v.N_Eids.size(); i++) {
            // std::cout << "getting vertex neighbor number: " << v.N_Eids.size() << std::endl;
            auto& e = mesh.E.at(v.N_Eids.at(i));
            size_t nvid = e.Vids[0] != v.id ? e.Vids[0] : e.Vids[1];
            // std::cout << "e: " << e.Vids[0] << " " << e.Vids[1] << std::endl;
            // std::cout << "nvid: " << nvid << std::endl;
            std::vector<size_t> nvids;
            for (auto fid: e.N_Fids) {
                auto& f = mesh.F.at(fid);
                // std::cout << "f: " << f.Vids.at(0) << " " << f.Vids.at(1) << " " << f.Vids.at(2) << " " << f.Vids.at(3) << std::endl;
                for (auto fvid: f.Vids) {
                    if (fvid != v.id && fvid != nvid && std::find(v.N_Vids.begin(), v.N_Vids.end(), fvid) != v.N_Vids.end()) {
                        nvids.push_back(fvid);
                        break;
                    }
                }
            }
            // std::cout << "neighbor vertices: " << nvids.size() << std::endl;
            auto& v_a = mesh.V.at(nvid);
            auto& v_b = mesh.V.at(nvids.at(0));
            auto& v_c = mesh.V.at(nvids.at(1));
            // std::cout << "got neighbor vertices to smooth" << std::endl;
            glm::dvec3 V_j = v.xyz() - v_a.xyz();
            glm::dvec3 V_j_minus_1 = v_b.xyz() - v_a.xyz(); 
            glm::dvec3 V_j_plus_1 = v_c.xyz() - v_a.xyz();

            double r = glm::length(V_j);
            glm::dvec3 p1 = v_a.xyz() + (r * glm::normalize(V_j_minus_1));
            glm::dvec3 p2 = v_a.xyz() + (r * glm::normalize(V_j_plus_1));
            
            double l = glm::length((0.5 * (p1 + p2)) - v_a.xyz());
            if (l == 0) continue;
            glm::dvec3 direction = glm::normalize((0.5 * (p1 + p2)) - v_a.xyz());
            glm::dvec3 newPoint_p = v_a.xyz() + (r * direction);
            glm::dvec3 newPoint_n = v_a.xyz() + (r * (-direction));
            
            center += glm::distance(v.xyz(), newPoint_p) < glm::distance(v.xyz(), newPoint_n) ? newPoint_p : newPoint_n;
            n += 1;
            // std::cout << "n: " << n << std::endl;
        }
        // std::cout << "centers size: " << centers.size() << std::endl;
        // std::cout << "centers: " << centers.at(iter).x << " " << centers.at(iter).y << " " << centers.at(iter).z << std::endl; 
        // std::cout << "center: " << center.x << " " << center.y << " " << center.z << std::endl; 
        centers.at(iter) = (center / n);
    }
}

