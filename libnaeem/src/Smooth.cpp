#include "Smooth.h"

Smoother::Smoother() {}
Smoother::Smoother(Mesh& mesh_) : mesh(mesh_) {
    sm.SetTarget(mesh);
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


    // smooth and project
	int iters = 10;

    // SurfaceMapper sm(mesh, origMesh);
	while (iters--) {
        // std::vector<glm::dvec3> centers(V.size(), glm::dvec3(0.0, 0.0, 0.0));
        std::vector<glm::dvec3> centers;
        for (auto vid: V) {
            auto& v = mesh.V.at(vid);
            centers.push_back(v.xyz());
            if (v.type < FEATURE) {
                // centers.at(v.id) = v.xyz();
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
                centers.back() = (center / n);
                // centers.at(vid) = (center / n);
            }
        }
        std::reverse(centers.begin(), centers.end());
        for (auto vid: V) {
            auto& v = mesh.V.at(vid);
            if (v.type < FEATURE) {
                v = sm.GetClosestPoint(centers.back());
                // v = centers.back();
            }
            centers.pop_back();
        }
    }
}
