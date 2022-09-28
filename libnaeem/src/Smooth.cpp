#include "Smooth.h"
#include "ParallelFor.h"
#include <math.h>

#define PI 3.14159265

Smoother::Smoother() {}
Smoother::Smoother(Mesh& mesh_, MeshUtil& mu_) : mesh(&mesh_), mu(&mu_) {
}

Smoother::Smoother(const Smoother& r) {
    mesh = r.mesh;
}

Smoother::~Smoother() {}

void Smoother::CheckValidity() {
    if (mesh == NULL) {
        std::cout << "No mesh to use for Smoother." << std::endl;
        exit(0);
    }
    if (mu == NULL) {
        std::cout << "No MeshUtil initialized for Smoother." << std::endl;
        exit(0);
    }
    if (mesh->V.size() == 0 || mesh->F.size() == 0 || mesh->C.size() == 0) {
        std::cout << "No mesh to use for Smoother." << std::endl;
        exit(0);
    }
}

void Smoother::SetMembers(Mesh& mesh_, MeshUtil& mu_) {
    mesh = &mesh_;
    mu = &mu_;
}

void Smoother::Smooth(std::vector<size_t>& V) {
    CheckValidity();
    bool performMapping = V.empty() ? true : false;
    if (V.empty()) {
        V.resize(mesh->V.size());
        if (V.size() >= 2000) {
            PARALLEL_FOR_BEGIN(mesh->V.size()) {
                GetVerticesToSmooth(i, V);
            } PARALLEL_FOR_END();    
        } else {
            for (int i = 0; i < mesh->V.size(); i++) {
                GetVerticesToSmooth(i, V);
            }
        }
    }
    // std::cout << "vertices: " << V.size() << std::endl;
    // sm.SetLocator(V);
    // smooth and project
	int iters_ = 10;
    // std::cout << "Smoothing " << V.size() << " vertices" << std::endl;

    // SurfaceMapper sm(origMesh);
	while (iters_--) {
        // centers.clear();
        // centers.resize(V.size());
        // std::vector<glm::dvec3> centers(V.size());
        // if (V.size() >= 2000) {
            // PARALLEL_FOR_BEGIN(V.size()) {
                // GetOptimizedPositions(i, V, centers);
            // } PARALLEL_FOR_END();
        
            // PARALLEL_FOR_BEGIN(V.size()) {
        //         // mesh->V.at(V.at(i)) = centers.at(i);
                // SetCoords(V.at(i), centers.at(i));
                // sm.RemapVertex(V.at(i), centers.at(i));
            // } PARALLEL_FOR_END();
        
        // } else {
            // for (int i = 0; i < iters; i++) {
            for (int i = 0; i < V.size(); i++) {
                // GetOptimizedPositions(i, V, centers);
                SetOptimizedPositions(i, V);
            }
            // for (int i = 0; i < V.size(); i++) {
                // mesh->V.at(V.at(i)) = centers.at(i);
                // SetCoords(V.at(i), centers.at(i));
                // SetCoords(V.at(i), sm.GetClosestPoint(centers.at(i)));
                // std::cout << "mapping vertex: " << V.at(i) << " from " << centers.at(i).x << " " << centers.at(i).y << " " << centers.at(i).z << std::endl;
                // sm.RemapVertex(V.at(i), centers.at(i));
            // }
        // }
    }
    // if (!performMapping) return;
    // for (int i = 0; i < V.size(); i++) {
    //     // mesh->V.at(V.at(i)) = sm.GetClosestPoint(mesh->V.at(V.at(i)).xyz());
    //     SetCoords(V.at(i), sm.GetClosestPoint(mesh->V.at(V.at(i)).xyz()));
    // }
}

void Smoother::GetVerticesToSmooth(int iter, std::vector<size_t>& V) {
    V.at(iter) = mesh->V.at(iter).id;
}

void Smoother::GetOptimizedPositions(int iter, std::vector<size_t>& V, std::vector<glm::dvec3>& centers) {    
    // std::cout << "iter: " << iter << " " << V.at(iter) << std::endl;
    auto& v = mesh->V.at(V.at(iter));
    // std::cout << "v: " << v.id << std::endl;
    centers.at(iter) = v.xyz();
    if (v.N_Vids.empty() || v.N_Eids.empty() || v.N_Fids.empty()) return;
    if (v.type != FEATURE && !v.isBoundary) {
        double n = 0.0;
        glm::dvec3 center(0.0, 0.0, 0.0);
        // std::cout << "vertices neighbor edges access " << v.id << std::endl;
        for (int i = 0; i < v.N_Eids.size(); i++) {
            // std::cout << "getting vertex neighbor number: " << v.N_Eids.size() << std::endl;
            auto& e = mesh->E.at(v.N_Eids.at(i));
            size_t nvid = e.Vids[0] != v.id ? e.Vids[0] : e.Vids[1];
            // std::cout << "e: " << e.Vids[0] << " " << e.Vids[1] << std::endl;
            // std::cout << "nvid: " << nvid << std::endl;
            std::vector<size_t> nvids;
            for (auto fid: e.N_Fids) {
                auto& f = mesh->F.at(fid);
                // std::cout << "f: " << f.Vids.at(0) << " " << f.Vids.at(1) << " " << f.Vids.at(2) << " " << f.Vids.at(3) << std::endl;
                int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), nvid));
                nvids.push_back(f.Vids.at((idx+2)%f.Vids.size()));
                // for (auto fvid: f.Vids) {
                //     if (fvid != v.id && fvid != nvid && std::find(v.N_Vids.begin(), v.N_Vids.end(), fvid) != v.N_Vids.end()) {
                //         nvids.push_back(fvid);
                //         break;
                //     }
                // }
            }
            // std::cout << "neighbor vertices: " << nvids.size() << std::endl;
            auto& v_a = mesh->V.at(nvid);
            auto& v_b = mesh->V.at(nvids.at(0));
            auto& v_c = mesh->V.at(nvids.at(1));
            // std::cout << "got neighbor vertices to smooth" << std::endl;
            glm::dvec3 V_j = v.xyz() - v_a.xyz();
            // double r = glm::length(V_j);
            glm::dvec3 V_j_minus_1 = v_b.xyz() - v_a.xyz(); 
            // r = glm::length(V_j_minus_1) < r ? glm::length(V_j_minus_1) : r;
            glm::dvec3 V_j_plus_1 = v_c.xyz() - v_a.xyz();
            // r = glm::length(V_j_plus_1) < r ? glm::length(V_j_plus_1) : r;

            double r = glm::length(V_j);

            // double r = glm::length(V_j_minus_1) < glm::length(V_j_plus_1) ? glm::length(V_j_minus_1) : glm::length(V_j_plus_1);
            glm::dvec3 p1 = v_a.xyz() + (r * glm::normalize(V_j_minus_1));
            glm::dvec3 p2 = v_a.xyz() + (r * glm::normalize(V_j_plus_1));
            glm::dvec3 p = (0.5 * (p1 + p2)) - v_a.xyz();
            // double a = glm::length(V_j_minus_1);
            // double b = glm::length(V_j_plus_1);
            // double r = a / (a+b);
            // glm::dvec3 direction = glm::normalize(v_c.xyz()-v_b.xyz());
            // glm::dvec3 p = v_b.xyz() + (r * direction);

            double l = glm::length(p);
            if (l == 0) continue;
            glm::dvec3 direction = glm::normalize(p);
            glm::dvec3 newPoint_p = v_a.xyz() + (r * direction);
            glm::dvec3 newPoint_n = v_a.xyz() + (r * (-direction));
            
            center += glm::distance(v.xyz(), newPoint_p) < glm::distance(v.xyz(), newPoint_n) ? newPoint_p : newPoint_n;
            // center += p;
            n += 1.0;
            // std::cout << "n: " << n << std::endl;
        }
        // std::cout << "centers size: " << centers.size() << std::endl;
        // std::cout << "centers: " << centers.at(iter).x << " " << centers.at(iter).y << " " << centers.at(iter).z << std::endl; 
        // std::cout << "center: " << center.x << " " << center.y << " " << center.z << std::endl;
        if (n > 0) {
            centers.at(iter) = (center / n);
        }
    } else if (v.type == FEATURE || v.isBoundary) {

    }
}

void Smoother::SetOptimizedPositions(int iter, std::vector<size_t>& V) {
    auto& v = mesh->V.at(V.at(iter));
    if (v.N_Vids.empty() || v.N_Eids.empty() || v.N_Fids.empty() || v.N_Vids.size() < 3) return;
    if (v.type == FEATURE || v.isBoundary) {
        // std::cout << "Before setting boundary vertex" << std::endl;
        SetPositionBoundary(v);
        // std::cout << "After setting boundary vertex" << std::endl;
    } else {
        // std::cout << "Before setting non-boundary vertex" << std::endl;
        SetPosition(v);
        // std::cout << "After setting non-boundary vertex" << std::endl;
    }
}

void Smoother::SetPosition(Vertex& v) {
    double polyArea = 0.0;
    glm::dvec3 centroid(0.0, 0.0, 0.0);
    size_t startE = v.N_Eids.at(0);
    
    for (int i = 0; i < v.N_Eids.size(); i++) {
        auto& edge = mesh->E.at(startE);
        size_t ev = edge.Vids.at(0) == v.id ? edge.Vids.at(1) : edge.Vids.at(0);
        size_t ev_plus1;
        size_t ev_minus1;
        for (auto fid: edge.N_Fids) {
            auto& f = mesh->F.at(fid);
            int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
            if (f.Vids.at((idx+1)%f.Vids.size()) == ev) {
                ev_plus1 = f.Vids.at((idx+3)%f.Vids.size());
                startE = mu->GetDifference(mu->GetIntersection(f.Eids, v.N_Eids), std::vector<size_t>{edge.id}).at(0);
            }
            if (f.Vids.at((idx+3)%f.Vids.size()) == ev) {
                ev_minus1 = f.Vids.at((idx+1)%f.Vids.size());
            }
        }

        auto& v2 = mesh->V.at(ev);
        auto& v3 = mesh->V.at(ev_plus1);
        auto& v4 = mesh->V.at(ev_minus1);

 
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
    if (polyArea == 0.0) return;
    centroid /= polyArea;
    for (auto fid: v.N_Fids) {
        auto& f = mesh->F.at(fid);
        int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
        if (idx == -1) continue;
        auto& v2 = mesh->V.at(f.Vids.at((idx+1)%f.Vids.size()));
        auto& v3 = mesh->V.at(f.Vids.at((idx+3)%f.Vids.size()));

        glm::dvec3 AB = v2.xyz() - v.xyz();
        glm::dvec3 AC = v3.xyz() - v.xyz();
        glm::dvec3 BC = v3.xyz() - v2.xyz();
        glm::dvec3 CA = v.xyz() - v3.xyz();

        glm::dvec3 T_cross = glm::cross(AB, AC);
        if (glm::length(T_cross) == 0.0) continue;

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
    }
}

void Smoother::SetPositionBoundary(Vertex& v) {
    std::vector<size_t> boundaryVertices;
    // std::cout << v.N_Vids.size() << std::endl;
    for (auto vid: v.N_Vids) {
        // std::cout << (mesh->V.at(vid).type == FEATURE) << ", " << mesh->V.at(vid).isBoundary << " ";
        if (mesh->V.at(vid).type == FEATURE || mesh->V.at(vid).isBoundary) boundaryVertices.push_back(vid);
    }
    // std::cout << std::endl;
    // std::cout << boundaryVertices.size() << std::endl;
    if (boundaryVertices.size() < 2) {
        // std::cout << "Before setting non-boundary vertex" << std::endl;
        SetPosition(v);
        // std::cout << "After setting non-boundary vertex" << std::endl;
        return;
    }
    if (boundaryVertices.size() > 2) return;
    // std::cout << "Before getting boundary vertices" << std::endl;
    auto& b1 = mesh->V.at(boundaryVertices.at(0));
    auto& b2 = mesh->V.at(boundaryVertices.at(1));
    // std::cout << "Two boundaries present" << std::endl;
    glm::dvec3 vb1 = glm::normalize(b1.xyz() - v.xyz());
    glm::dvec3 vb2 = glm::normalize(b2.xyz() - v.xyz());

    double angle = acos(glm::dot(vb1, vb2)) * 180.0 / PI;
    if (angle < 150) return;
    
    glm::dvec3 centroid(0.0, 0.0, 0.0);
    int k = 0;
    for (auto vid: v.N_Vids) {
        auto& nv = mesh->V.at(vid);
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
    if (k > 0) {
        centroid /= (double) k;
        // std::cout << centroid.x << " " << centroid.y << " " << centroid.z << std::endl;
        v.x = centroid.x;
        v.y = centroid.y;
        v.z = centroid.z;
    }
}

void Smoother::SetCoords(size_t vid, glm::dvec3& c) {
    mesh->V.at(vid).x = c.x;
    mesh->V.at(vid).y = c.y;
    mesh->V.at(vid).z = c.z;
}


