#include <math.h>
#include <algorithm>
#include <glm/glm.hpp>
#include <limits>
#include <numeric>
#include "AngleBasedSmoothQuadMesh.h"

#define PI 3.14159265

SmoothAlgorithm::SmoothAlgorithm(Mesh& mesh, int it, double l_r) : mesh(mesh) {
    iters = it;
    lambda = l_r;
    if (lambda <= 0) {
        lambda = 0.1;
    } else if (lambda > 1) {
        lambda = 1;
    }
    std::cout << "min displacement limit: " << min_displacement_limit << std::endl;
}
SmoothAlgorithm::~SmoothAlgorithm() {}

void SmoothAlgorithm::setOriginalVertices() {
    original_vertices.resize(mesh.V.size(), glm::dvec3(0.0, 0.0, 0.0));
    for (int i = 0; i < mesh.V.size(); i++) {
        if (!mesh.V.at(i).isBoundary) {
            continue;
        }
        original_vertices.at(i).x = mesh.V.at(i).x;
        original_vertices.at(i).y = mesh.V.at(i).y;
        original_vertices.at(i).z = mesh.V.at(i).z;
        original_vertices.at(i).isSingularity = mesh.V.at(i).isSingularity;
        original_vertices.at(i).isBoundary = mesh.V.at(i).isBoundary;
        original_vertices.at(i).isCorner = mesh.V.at(i).isCorner;
        original_vertices.at(i).N_Vids.resize(mesh.V.at(i).N_Vids.size());
        for (int j = 0; j < mesh.V.at(i).N_Vids.size(); j++) {
            original_vertices.at(i).N_Vids.at(j) = mesh.V.at(i).N_Vids.at(j);
        }
    }
}

double SmoothAlgorithm::getMinEdgeLength() {
    double min_edge_length = -1;
    for (auto& e: mesh.E) {
        Vertex& v1 = mesh.V.at(e.Vids[0]);
        Vertex& v2 = mesh.V.at(e.Vids[1]);
        double edge_length = glm::length(glm::dvec3(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z));
        if (min_edge_length == -1) {
            min_edge_length = edge_length;
        } else if (edge_length < min_edge_length) {
            min_edge_length = edge_length;
        }
    }
    return min_edge_length;
}

void SmoothAlgorithm::resampleBoundaryVertices() {
    for (int i = 0; i < mesh.V.size(); i++) {
        Vertex& v = mesh.V.at(i);
        if (v.isBoundary) {
            if (!v.isMovable) {
                delta_coords.at(i).x = v.x;
                delta_coords.at(i).y = v.y;
                continue;
            }
            int nb_n_index = -1;
            for (int j = 0; j < v.N_Vids.size(); j++) {
                if (!mesh.V.at(v.N_Vids.at(j)).isBoundary) {
                    nb_n_index = j;
                    break;
                }                    
            }
            if (nb_n_index < 0) {
                delta_coords.at(i).x = v.x;
                delta_coords.at(i).y = v.y;
                continue;
            }
            int a_ = nb_n_index;
            int b_ = v.N_Vids.size();
            int index = (a_ % b_ + b_) % b_;
            Vertex& b0 = mesh.V.at(v.N_Vids.at(index));

            a_ = nb_n_index - 1;
            index = (a_ % b_ + b_) % b_;
            Vertex& b1 = mesh.V.at(v.N_Vids.at(index));

            a_ = nb_n_index + 1;
            index = (a_ % b_ + b_) % b_;
            Vertex& b2 = mesh.V.at(v.N_Vids.at(index));

            glm::vec3 a(b1.x - v.x, b1.y - v.y, b1.z - v.z);
            glm::vec3 b(b2.x - v.x, b2.y - v.y, b2.z - v.z);

            double angle = acos(glm::dot(a, b) / (glm::length(a) * glm::length(b))) * 180.0 / PI;

            glm::dvec3 V_j(v.x - b0.x, v.y - b0.y, v.z - b0.z);
            glm::dvec3 V_j_minus_1(b1.x - b0.x, b1.y - b0.y, b1.z - b0.z);
            glm::dvec3 V_j_plus_1(b2.x - b0.x, b2.y - b0.y, b2.z - b0.z);

            double alpha1 = acos(glm::dot(V_j, V_j_plus_1) / (glm::length(V_j) * glm::length(V_j_plus_1)));
            double alpha2 = acos(glm::dot(V_j, V_j_minus_1) / (glm::length(V_j) * glm::length(V_j_minus_1)));

            double beta = (alpha2 - alpha1) / 2;
            delta_coords.at(i).x = (b0.x + (V_j.x * cos(beta)) - (V_j.y * sin(beta))); 
            delta_coords.at(i).y = (b0.y + (V_j.x * sin(beta)) + (V_j.y * cos(beta)));
        }
    }
    // std::cout << "Remapping vertices to original boundary" << std::endl;
    for (int i = 0; i < mesh.V.size(); i++) {
        if (mesh.V.at(i).isBoundary) {
            if (delta_coords.at(i).x == mesh.V.at(i).x && delta_coords.at(i).y == mesh.V.at(i).y) {
                continue;
            }
            glm::dvec3 new_v = delta_coords.at(i);
            double min_length = std::numeric_limits<double>::max();
            int min_index = i;
            // std::cout << "original vertices size: " << original_vertices.size() << std::endl;
            // std::cout << "min_index: " << min_index << std::endl;
            int nBoundaryVertices = 0;
            for (int j = 0; j < original_vertices.size(); j++) {
                if (original_vertices.at(j).isBoundary) {
                    nBoundaryVertices += 1;
                    double length = glm::length(glm::vec3(new_v.x - original_vertices.at(j).x, new_v.y - original_vertices.at(j).y, new_v.z - original_vertices.at(j).z));
                    if (length < min_length) {
                        min_length = length;
                        min_index = j;
                    }
                }
            }
            // std::cout << "num bounadry vertices: " << nBoundaryVertices << std::endl;
            // std::cout << "min index" << min_index << std::endl;
            Vertex& v = original_vertices.at(min_index);
            // std::cout << "v neighbors: " << v.N_Vids.size() << std::endl;
            std::vector<size_t> boundary_neighbors;
            for (int j = 0; j < v.N_Vids.size(); j++) {
                if (original_vertices.at(v.N_Vids.at(j)).isBoundary) {
                    boundary_neighbors.push_back(v.N_Vids.at(j));
                }                    
            }
            // std::cout << "boundary neighbor size: " << boundary_neighbors.size() << std::endl;
            Vertex& b1 = original_vertices.at(boundary_neighbors.at(0));
            Vertex& b2 = original_vertices.at(boundary_neighbors.at(1));
            
            glm::vec3 a(new_v.x - v.x, new_v.y - v.y, new_v.z - v.z);
            glm::vec3 b(b1.x - v.x, b1.y - v.y, b1.z - v.z);
            glm::vec3 c(b2.x - v.x, b2.y - v.y, b2.z - v.z);

            double length_a = glm::dot(a, a);
            double length_b = glm::dot(b, b);
            double length_c = glm::dot(c, c);

            double dot_a_b = glm::dot(a, b) / glm::length(b);
            double dot_a_c = glm::dot(a, c) / glm::length(c);

            if (dot_a_b >= 0) {
                b = glm::normalize(b);
                delta_coords.at(i).x = v.x + (dot_a_b * b.x);
                delta_coords.at(i).y = v.y + (dot_a_b * b.y);
            } else if (dot_a_c >= 0){
                c = glm::normalize(c);
                delta_coords.at(i).x = v.x + (dot_a_c * c.x);
                delta_coords.at(i).y = v.y + (dot_a_c * c.y);
            }
        }
    }
}

void SmoothAlgorithm::smoothLaplacianSimple() {
    int it = 0;
    double dt = 0;
    std::vector<Vertex> buffer_V(mesh.V.size());
    glm::dvec3 max_displacement(0, 0, 0);
    while (it < iters) {
        double min_edge_length = getMinEdgeLength();
        // dt = (min_edge_length * min_edge_length) / (2 * lambda);
        dt = 1;
        for (int i = 0; i < mesh.V.size(); i++) {
            Vertex& v = mesh.V.at(i);
            if (v.isBoundary) {
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

            buffer_V.at(i).x = lambda * dt * ((x_agg / weight_agg) - v.x);
            buffer_V.at(i).y = lambda * dt * ((y_agg / weight_agg) - v.y);
            // if (!v.isBoundary && i < mesh.V.size() / 2) {
            //     std::cout << (x_agg / weight_agg) - v.x << std::endl;
            //     std::cout << (y_agg / weight_agg) - v.y << std::endl;
            // }
            // glm::dvec3 disp_(buffer_V.at(i).x, buffer_V.at(i).y, buffer_V.at(i).z);
            // std::cout << glm::length(disp_) << std::endl;
            // if (glm::length(disp_) > glm::length(max_displacement)) {
            //     max_displacement.x = disp_.x;
            //     max_displacement.y = disp_.y;
            //     max_displacement.z = disp_.z;
            // }
        }

        // std::cout << glm::length(max_displacement) << std::endl;
        // if (glm::length(max_displacement) < min_displacement_limit) {
        //     std::cout << "Finished after " << it << " iterations" << std::endl;
        //     break;
        // }

        for (int i = 0; i < mesh.V.size(); i++) {
            mesh.V.at(i).x += buffer_V.at(i).x;
            mesh.V.at(i).y += buffer_V.at(i).y;
            buffer_V.at(i).x = 0;
            buffer_V.at(i).y = 0;
        }
        // std::cout << "-----------------------------------" << std::endl;
        it++;
    }
    calculateMeshAngles();
}

void SmoothAlgorithm::smoothLaplacianScaleBased() {
    calculateMeshAngles();
    int it = 0;
    double dt = 0;
    std::vector<Vertex> buffer_V(mesh.V.size());
    glm::dvec3 max_displacement(0, 0, 0);
    while (it < iters) {
        double min_edge_length = getMinEdgeLength();
        // dt = (min_edge_length * min_edge_length) / (2 * lambda);
        dt = 1;
        for (int i = 0; i < mesh.V.size(); i++) {
            Vertex& v = mesh.V.at(i);
            if (v.isBoundary) {
                continue;
            }
            double n = v.N_Vids.size();
            double weight = 0;
            double weight_agg = 0;
            double x_agg = 0;
            double y_agg = 0;
            double edge_agg = 0;
            for (int j = 0; j < n; j++) {
                Vertex& v_prime = mesh.V.at(v.N_Vids.at(j));
                edge_agg += glm::length(glm::dvec3(v_prime.x - v.x, v_prime.y - v.y, v_prime.z - v.z));
            }
            for (int j = 0; j < n; j++) {
                Vertex& v_prime = mesh.V.at(v.N_Vids.at(j));
                weight = 1 - (glm::length(glm::dvec3(v_prime.x - v.x, v_prime.y - v.y, v_prime.z - v.z)) / edge_agg);
                // weight = (glm::length(glm::dvec3(v_prime.x - v.x, v_prime.y - v.y, v_prime.z - v.z)) / edge_agg);
                // std::cout << weight << std::endl;
                weight_agg += weight;
                x_agg += (weight * v_prime.x);
                y_agg += (weight * v_prime.y);
            }

            buffer_V.at(i).x = lambda * dt * ((x_agg / weight_agg) - v.x);
            buffer_V.at(i).y = lambda * dt * ((y_agg / weight_agg) - v.y);
            // if (!v.isBoundary && i < mesh.V.size() / 2) {
            //     std::cout << (x_agg / weight_agg) - v.x << std::endl;
            //     std::cout << (y_agg / weight_agg) - v.y << std::endl;
            //     std::cout << x_agg << std::endl;
            //     std::cout << y_agg << std::endl;
            //     std::cout << "*****************************" << std::endl;
            // }
            // glm::dvec3 disp_(buffer_V.at(i).x, buffer_V.at(i).y, buffer_V.at(i).z);
            // std::cout << glm::length(disp_) << std::endl;
            // if (glm::length(disp_) > glm::length(max_displacement)) {
            //     max_displacement.x = disp_.x;
            //     max_displacement.y = disp_.y;
            //     max_displacement.z = disp_.z;
            // }
        }

        // std::cout << glm::length(max_displacement) << std::endl;
        // if (glm::length(max_displacement) < min_displacement_limit) {
        //     std::cout << "Finished after " << it << " iterations" << std::endl;
        //     break;
        // }

        for (int i = 0; i < mesh.V.size(); i++) {
            mesh.V.at(i).x += buffer_V.at(i).x;
            mesh.V.at(i).y += buffer_V.at(i).y;
            buffer_V.at(i).x = 0;
            buffer_V.at(i).y = 0;
        }
        // std::cout << "-----------------------------------" << std::endl;
        it++;
    }
    calculateMeshAngles();
}

void SmoothAlgorithm::setBoundaryVerticesMovable() {
    for (auto& v: mesh.V) {
        if (!v.isBoundary) {
            int n_b = 0;
            for (auto id: v.N_Vids) {
                if (mesh.V.at(id).isBoundary) {
                    n_b += 1;
                }
            }
            if (n_b > 1) {
                for (auto id: v.oneRingNeighborVertices) {
                    if (mesh.V.at(id).isBoundary) {
                        mesh.V.at(id).isMovable = false;
                    }
                }
            }
        }
    }
    for (auto& v: mesh.V) {
        if (v.isBoundary) {
            if (v.N_Fids.size() == 2) {
                v.isMovable = true;
            } else {
                v.isMovable = false;
            }
        }
    }
}

void SmoothAlgorithm::findNegativeElements() {
    for (auto& f: mesh.F) {
        int sign = 0;
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
        if (abs(sign) != 4) {
            f.isNegative = true;
        } else {
            f.isNegative = false;
        }
    }
}

void SmoothAlgorithm::fixNegativeElements() {
    std::vector<glm::dvec3> new_coords(mesh.V.size(), glm::dvec3(0.0, 0.0, 0.0));
    for (auto& f: mesh.F) {
        if (f.isNegative) {
            for (auto id: f.Vids) {
                Vertex& v = mesh.V.at(id);
                 if (!v.isBoundary) {
                    continue;
                }
                std::vector<size_t> one_ring_neighbors = v.N_Vids;
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
                    // if (std::find(v.N_Vids.begin(), v.N_Vids.end(), v_j.id) == v.N_Vids.end()) {
                    //     a = j - 1;
                    // }
                    index = one_ring_neighbors.at((a % b + b) % b);
                    Vertex& v_j_prev = mesh.V.at(index);

                    a = j + 1;
                    // if (std::find(v.N_Vids.begin(), v.N_Vids.end(), v_j.id) == v.N_Vids.end()) {
                    //     a = j + 1;
                    // }
                    index = one_ring_neighbors.at((a % b + b) % b);
                    Vertex& v_j_next = mesh.V.at(index);

                    glm::dvec3 V_j(v.x - v_j.x, v.y - v_j.y, v.z - v_j.z);
                    glm::dvec3 V_j_minus_1(v_j_prev.x - v_j.x, v_j_prev.y - v_j.y, v_j_prev.z - v_j.z);
                    glm::dvec3 V_j_plus_1(v_j_next.x - v_j.x, v_j_next.y - v_j.y, v_j_next.z - v_j.z);

                    double alpha1 = acos(glm::dot(V_j, V_j_plus_1) / (glm::length(V_j) * glm::length(V_j_plus_1)));
                    double alpha2 = acos(glm::dot(V_j, V_j_minus_1) / (glm::length(V_j) * glm::length(V_j_minus_1)));
                    
                    double beta = (alpha2 - alpha1) / 2;
                    neighbor_coords.at(j).x = (v_j.x + (V_j.x * cos(beta)) - (V_j.y * sin(beta))); 
                    neighbor_coords.at(j).y = (v_j.y + (V_j.x * sin(beta)) + (V_j.y * cos(beta)));
                    
                    neighbor_weights.at(j) = fabs(beta);
                    if (neighbor_weights.at(j) == 0) {
                        neighbor_weights.at(j) = 1;
                    }
                    weight_agg += neighbor_weights.at(j);
                }
                double vertex_weight = 0.0;
                for (int j = 0; j < one_ring_neighbors.size(); j++) {
                    glm::dvec3 current_v = neighbor_coords.at(j);
                    double weight = (1 - neighbor_weights.at(j) / weight_agg);
                    vertex_weight += weight;
                    new_coords.at(v.id).x += (weight * current_v.x);
                    new_coords.at(v.id).y += (weight * current_v.y);
                }
                new_coords.at(v.id).x = lambda * ((new_coords.at(v.id).x / vertex_weight) - v.x);
                new_coords.at(v.id).y = lambda * ((new_coords.at(v.id).y / vertex_weight) - v.y);
                v.isVisited = true;
            }
        }
    }
    for (int i = 0; i < mesh.V.size(); i++) {
        Vertex& v = mesh.V.at(i);
        if (!v.isVisited) {
            continue;
        }
        v.x = (v.x + new_coords.at(i).x);
        v.y = (v.y + new_coords.at(i).y);
        v.isVisited = false;
        // break;
    }
}

void SmoothAlgorithm::smoothMesh() {
    std::cout << "Started Mesh Smoothing" << std::endl;
    std::cout << mesh.V.size() << std::endl;
    setOriginalVertices();
    setBoundaryVerticesMovable();
    calculateMeshAngles();
    // int it = 0;
    // while (it < iters) {
    //     findNegativeElements();
    //     fixNegativeElements();
    //     it += 1;
    // }

    
    int it = 0;
    while (it < iters) {
        delta_coords.clear();
        delta_coords.resize(mesh.V.size(), glm::dvec3(0.0, 0.0, 0.0));
        // smoothLaplacianCotangentBased();
        SmoothAngleBased();
        // angleBasedSmoothing();        
        for (int i = 0; i < mesh.V.size(); i++) {
            Vertex& v = mesh.V.at(i);
            if (v.isBoundary) {
            //     v.x = delta_coords.at(i).x / 2;
            //     v.y = delta_coords.at(i).y / 2;
                continue;
            } 
            // else {
            v.x = delta_coords.at(i).x;
            v.y = delta_coords.at(i).y;
            // }
            // break;
        }
        it++;
    }
    
    it = 0;
    while (it < iters) {
        delta_coords.clear();
        delta_coords.resize(mesh.V.size(), glm::dvec3(0.0, 0.0, 0.0));
        // smoothLaplacianCotangentBased();
        resampleBoundaryVertices();
        for (int i = 0; i < mesh.V.size(); i++) {
            Vertex& v = mesh.V.at(i);
            if (!v.isBoundary) {
            //     v.x = delta_coords.at(i).x / 2;
            //     v.y = delta_coords.at(i).y / 2;
                continue;
            } 
            // else {
            v.x = delta_coords.at(i).x;
            v.y = delta_coords.at(i).y;
            // }
            // break;
        }
        it++;
    }
    
    std::cout << "Finished Mesh Smoothing" << std::endl;
    calculateMeshAngles();
}

void SmoothAlgorithm::smoothLaplacianCotangentBased() {
    // delta_coords.resize(mesh.V.size(), glm::dvec3(0.0, 0.0, 0.0));
    std::vector<glm::dvec3> new_coords(mesh.V.size(), glm::dvec3(0.0, 0.0, 0.0));
    std::vector<double> vertices_weights(mesh.V.size(), 0.0);
    
    for (int i = 0; i < mesh.V.size(); i++) {
        Vertex& v = mesh.V.at(i);
        if (v.isBoundary) {
            continue;
        }
        
        std::vector<size_t> one_ring_neighbors = v.oneRingNeighborVertices;
        std::vector<double> neighbor_weights(one_ring_neighbors.size(), 0.0);
        double weight_agg = 0.0;
        for (int j = 0; j < one_ring_neighbors.size(); j++) {
            int a = j;
            int b = one_ring_neighbors.size();
            int index = one_ring_neighbors.at((a % b + b) % b);
            Vertex& current_v = mesh.V.at(index);
            if (std::find(v.N_Vids.begin(), v.N_Vids.end(), current_v.id) == v.N_Vids.end()) {
                continue;
            }
            a = j - 2;
            if (std::find(v.N_Vids.begin(), v.N_Vids.end(), current_v.id) == v.N_Vids.end()) {
                a = j - 1;
            }
            index = one_ring_neighbors.at((a % b + b) % b);
            Vertex& prev_v = mesh.V.at(index);

            a = j + 2;
            if (std::find(v.N_Vids.begin(), v.N_Vids.end(), current_v.id) == v.N_Vids.end()) {
                a = j + 1;
            }
            index = one_ring_neighbors.at((a % b + b) % b);
            Vertex& next_v = mesh.V.at(index);

            double e1 = glm::length(glm::dvec3(v.x - current_v.x, v.y - current_v.y, v.z - current_v.z));
            double e2 = glm::length(glm::dvec3(v.x - prev_v.x, v.y - prev_v.y, v.z - prev_v.z));
            double e3 = glm::length(glm::dvec3(prev_v.x - current_v.x, prev_v.y - current_v.y, prev_v.z - current_v.z));

            double cos_alpha = fabs((e3*e3 + e2*e2 - e1*e1) / (2.0*e3*e2));

            double e4 = glm::length(glm::dvec3(v.x - next_v.x, v.y - next_v.y, v.z - next_v.z));
            double e5 = glm::length(glm::dvec3(next_v.x - current_v.x, next_v.y - current_v.y, next_v.z - current_v.z));

            double cos_beta = fabs((e4*e4 + e5*e5 - e1*e1) / (2.0*e4*e5));

            double cot_alpha = cos_alpha / sqrt(1.0 - cos_alpha*cos_alpha);
            double cot_beta = cos_beta / sqrt(1.0 - cos_beta*cos_beta);

            double weight = (cot_alpha + cot_beta) / 2.0;
            neighbor_weights.at(j) = weight;
            weight_agg += weight;

        }
        for (int j = 0; j < one_ring_neighbors.size(); j++) {
            Vertex& current_v = mesh.V.at(one_ring_neighbors.at(j));
            glm::dvec3 current_delta(current_v.x - v.x, current_v.y - v.y, current_v.z - v.z);
            double weight = 1 - (neighbor_weights.at(j) / weight_agg);
            vertices_weights.at(i) += weight;
            new_coords.at(i).x += (weight * current_v.x);
            new_coords.at(i).y += (weight * current_v.y);
        }
        // std::cout << "aggregated weights after normalization: " << vertices_weights.at(i) << std::endl;
        new_coords.at(i).x = lambda * ((new_coords.at(i).x / vertices_weights.at(i)) - v.x);
        new_coords.at(i).y = lambda * ((new_coords.at(i).y / vertices_weights.at(i)) - v.y);
        // break;
    }
    for (int i = 0; i < mesh.V.size(); i++) {
        Vertex& v = mesh.V.at(i);
        if (v.isBoundary) {
            continue;
        }
        delta_coords.at(i).x += (v.x + new_coords.at(i).x);
        delta_coords.at(i).y += (v.y + new_coords.at(i).y);
        // break;
    }
}

void SmoothAlgorithm::SmoothAngleBased() {
    // delta_coords.resize(mesh.V.size(), glm::dvec3(0.0, 0.0, 0.0));
    std::vector<glm::dvec3> new_coords(mesh.V.size(), glm::dvec3(0.0, 0.0, 0.0));
    for (int i = 0; i < mesh.V.size(); i++) {
        Vertex& v_i = mesh.V.at(i);
        if (v_i.isBoundary) {
            continue;
        }
        std::vector<size_t> one_ring_neighbors = v_i.oneRingNeighborVertices;
        int k = one_ring_neighbors.size();
        for (int j = 0; j < one_ring_neighbors.size(); j++) {
            int a = j;
            int b = one_ring_neighbors.size();
            int index = one_ring_neighbors.at((a % b + b) % b);
            Vertex& v_j = mesh.V.at(index);
            
            a = j - 2;
            if (std::find(v_i.N_Vids.begin(), v_i.N_Vids.end(), v_j.id) == v_i.N_Vids.end()) {
                a = j - 1;
            }
            index = one_ring_neighbors.at((a % b + b) % b);
            Vertex& v_j_prev = mesh.V.at(index);

            a = j + 2;
            if (std::find(v_i.N_Vids.begin(), v_i.N_Vids.end(), v_j.id) == v_i.N_Vids.end()) {
                a = j + 1;
            }
            index = one_ring_neighbors.at((a % b + b) % b);
            Vertex& v_j_next = mesh.V.at(index);

            glm::dvec3 V_j(v_i.x - v_j.x, v_i.y - v_j.y, v_i.z - v_j.z);
            glm::dvec3 V_j_minus_1(v_j_prev.x - v_j.x, v_j_prev.y - v_j.y, v_j_prev.z - v_j.z);
            glm::dvec3 V_j_plus_1(v_j_next.x - v_j.x, v_j_next.y - v_j.y, v_j_next.z - v_j.z);

            double alpha1 = acos(glm::dot(V_j, V_j_plus_1) / (glm::length(V_j) * glm::length(V_j_plus_1)));
            double alpha2 = acos(glm::dot(V_j, V_j_minus_1) / (glm::length(V_j) * glm::length(V_j_minus_1)));

            double beta = (alpha2 - alpha1) / 2;
            // delta_coords.at(i).x += (v_j.x + (V_j.x * cos(beta)) - (V_j.y * sin(beta))); 
            // delta_coords.at(i).y += (v_j.y + (V_j.x * sin(beta)) + (V_j.y * cos(beta)));
            new_coords.at(i).x += (v_j.x + (V_j.x * cos(beta)) - (V_j.y * sin(beta))); 
            new_coords.at(i).y += (v_j.y + (V_j.x * sin(beta)) + (V_j.y * cos(beta)));
        }
        // delta_coords.at(i).x = delta_coords.at(i).x / k;
        // delta_coords.at(i).y = delta_coords.at(i).y / k;
        new_coords.at(i).x = new_coords.at(i).x / k;
        new_coords.at(i).y = new_coords.at(i).y / k;
    }
    for (int i = 0; i < mesh.V.size(); i++) {
        Vertex& v = mesh.V.at(i);
        if (v.isBoundary) {
            continue;
        }
        delta_coords.at(i).x += new_coords.at(i).x;
        delta_coords.at(i).y += new_coords.at(i).y;
        // break;
    }
}

void SmoothAlgorithm::angleBasedSmoothing() {
    // delta_coords.resize(mesh.V.size(), glm::dvec3(0.0, 0.0, 0.0));
    std::vector<double> vertex_weights(mesh.V.size(), 0);
    std::vector<glm::dvec3> new_coords(mesh.V.size(), glm::dvec3(0.0, 0.0, 0.0));
    for (int i = 0; i < mesh.V.size(); i++) {
        Vertex& v_i = mesh.V.at(i);
        if (v_i.isBoundary) {
            continue;
        }
        std::vector<size_t> one_ring_neighbors = v_i.oneRingNeighborVertices;
        std::vector<glm::dvec3> neighbor_coords(one_ring_neighbors.size(), glm::dvec3(0.0, 0.0, 0.0));
        std::vector<double> neighbor_weights(one_ring_neighbors.size(), 0.0);
        double weight_agg = 0.0;
        int k = one_ring_neighbors.size();
        for (int j = 0; j < one_ring_neighbors.size(); j++) {
            int a = j;
            int b = one_ring_neighbors.size();
            int index = one_ring_neighbors.at((a % b + b) % b);
            Vertex& v_j = mesh.V.at(index);

            a = j - 2;
            if (std::find(v_i.N_Vids.begin(), v_i.N_Vids.end(), v_j.id) == v_i.N_Vids.end()) {
                a = j - 1;
            }
            index = one_ring_neighbors.at((a % b + b) % b);
            Vertex& v_j_prev = mesh.V.at(index);

            a = j + 2;
            if (std::find(v_i.N_Vids.begin(), v_i.N_Vids.end(), v_j.id) == v_i.N_Vids.end()) {
                a = j + 1;
            }
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
            if (neighbor_weights.at(j) == 0) {
                neighbor_weights.at(j) = 1;
            }
            weight_agg += neighbor_weights.at(j);
        }
        for (int j = 0; j < one_ring_neighbors.size(); j++) {
            glm::dvec3 current_v = neighbor_coords.at(j);
            double weight = (1 - neighbor_weights.at(j) / weight_agg);
            vertex_weights.at(i) += weight;
            new_coords.at(i).x += (weight * current_v.x);
            new_coords.at(i).y += (weight * current_v.y);
        }
        new_coords.at(i).x = lambda * ((new_coords.at(i).x / vertex_weights.at(i)) - v_i.x);
        new_coords.at(i).y = lambda * ((new_coords.at(i).y / vertex_weights.at(i)) - v_i.y);
    }

    for (int i = 0; i < mesh.V.size(); i++) {
        Vertex& v = mesh.V.at(i);
        if (v.isBoundary) {
            continue;
        }
        delta_coords.at(i).x += (v.x + new_coords.at(i).x);
        delta_coords.at(i).y += (v.y + new_coords.at(i).y);
        // break;
    }
}

void SmoothAlgorithm::calculateMeshAngles() {
    std::vector<int> angles_histogram(180, 0);
    double min_angle = -1;
    double max_angle = -1;
    double all_angles = 0;
    double agg_angles = 0;
    for (auto& f: mesh.F) {
        for (int i = 0; i < f.Vids.size(); i++) {
            Vertex& v1 = mesh.V.at(f.Vids.at(i));
            Vertex& v2 = mesh.V.at(f.Vids.at((i-1)%f.Vids.size()));
            Vertex& v3 = mesh.V.at(f.Vids.at((i+1)%f.Vids.size()));

            glm::dvec3 a(v2.x - v1.x, v2.y - v1.y, v2.z - v1.z);
            glm::dvec3 b(v3.x - v1.x, v3.y - v1.y, v3.z - v1.z);
            
            double angle = ceil(acos(glm::dot(a, b) / (glm::length(a) * glm::length(b))) * 180.0 / PI);
            agg_angles += angle;
            all_angles += 1;
            if ((int) angle < angles_histogram.size()) {
                angles_histogram[(int) angle] += 1;
            }
            if (min_angle == -1) {
                min_angle = angle;
            } else if (angle < min_angle) {
                min_angle = angle;
            }

            if (max_angle == -1) {
                max_angle = angle;
            } else if (angle > max_angle) {
                max_angle = angle;
            }
        }
    }
    
    std::cout << "min_angle: " << min_angle << std::endl;
    std::cout << "max_angle: " << max_angle << std::endl;
    std::cout << "average angle: " << agg_angles / all_angles << std::endl;
    // std::cout << "histogram: " << std::endl;
    // for (int i = 0; i < angles_histogram.size(); i++) {
    //     std::cout << angles_histogram[i] << std::endl;
    // }
}

bool SmoothAlgorithm::isMeshNonManifold() {
    std::cout << "checking mesh manifold" << std::endl;
    std::vector<int> vs;
    bool non_mainfold = false;
    for (int i = 0; i < mesh.V.size(); i++) {
        Vertex& v = mesh.V.at(i);
        int nFaces = v.N_Fids.size();
        int nVertices = v.N_Vids.size();
        if (v.isBoundary) {
            nVertices = nVertices - 1;
            // continue;
        }
        if (nFaces != nVertices) {
            vs.push_back(i);
            non_mainfold = true;
            // return true;
        }
    }

    if (non_mainfold) {
        std::cout << "Mesh has non-manifold vertices" << std::endl;
        std::cout << "Writing output file" << std::endl;
        std::ofstream ofs("output_v.vtk");
        ofs << "# vtk DataFile Version 3.0\n"
            << "output.vtk\n"
            << "ASCII\n\n"
            << "DATASET UNSTRUCTURED_GRID\n";
        ofs << "POINTS " << vs.size() << " double\n";
        for (auto id: vs) {
            ofs <<  mesh.V.at(id).x << " " <<  mesh.V.at(id).y << " " <<  mesh.V.at(id).z << "\n";
        } 
        ofs << "CELLS " << vs.size() << " " << vs.size() * 2  << std::endl;
        for (int i = 0; i < vs.size(); i++) {
            ofs << "1 " << i << std::endl;
        }
        ofs << "CELL_TYPES " << vs.size() << "\n";
        for (int i = 0; i < vs.size(); i++) {
            ofs << "1 " << std::endl;
        }
    }

    for (int i = 0; i < mesh.E.size(); i++) {
        Edge& e = mesh.E.at(i);
        if (!e.isBoundary) {
            Face& f1 = mesh.F.at(e.N_Fids.at(0));
            Face& f2 = mesh.F.at(e.N_Fids.at(1));
            int id1 = 0, id2 = 0, id3 = 0, id4 = 0;
            for (int j = 0; j < f1.Vids.size(); j++) {
                int index1 = j;
                int index2 = (j + 1) % f1.Vids.size();
                if (std::find(e.Vids.begin(), e.Vids.end(), f1.Vids.at(index1)) != e.Vids.end() &&
                    std::find(e.Vids.begin(), e.Vids.end(), f1.Vids.at(index2)) != e.Vids.end()) {
                    id1 = f1.Vids.at(j);
                    id2 = f1.Vids.at((j + 1) % f1.Vids.size());
                    break;
                }
            }
            for (int j = 0; j < f2.Vids.size(); j++) {
                int index1 = j;
                int index2 = (j + 1) % f2.Vids.size();
                if (std::find(e.Vids.begin(), e.Vids.end(), f2.Vids.at(index1)) != e.Vids.end() &&
                    std::find(e.Vids.begin(), e.Vids.end(), f2.Vids.at(index2)) != e.Vids.end()) {
                    id3 = f2.Vids.at(j);
                    id4 = f2.Vids.at((j + 1) % f2.Vids.size());
                    break;
                }
            }
            if (id1 != id4 && id2 != id3) {
                std::cout << "Mesh has non-manifold edges" << std::endl;
                return true;
            }
        }
    }

    for (int i = 0; i < mesh.E.size(); i++) {
        Edge& e = mesh.E.at(i);
        if (e.N_Fids.size() > 2) {
            std::cout << "Mesh has non-manifold edges" << std::endl;
            return true;
        }
    }

    if (!mesh.isManifold) {
        return true;
    }
    return false;
}
