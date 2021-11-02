#include <math.h>
#include <algorithm>
#include <glm/glm.hpp>
#include <limits>
#include <numeric>
#include <ctime>
#include "AngleBasedSmoothQuadMesh.h"

#define PI 3.14159265

SmoothAlgorithm::SmoothAlgorithm(Mesh& mesh, Mesh& boundary_mesh, int it, double l_r, bool global_, bool boundary_smoothing) : mesh(mesh), boundary_mesh(boundary_mesh) {
    iters = it;
    lambda = l_r;
    // tau = t;
    if (lambda <= 0) {
        lambda = 0.1;
    } else if (lambda > 1) {
        lambda = 1;
    }
    global = global_;
    boundarySmoothing = boundary_smoothing;
    std::cout << "min displacement limit: " << min_displacement_limit << std::endl;
}
SmoothAlgorithm::~SmoothAlgorithm() {}

void SmoothAlgorithm::setOriginalVertices() {
    original_vertices.resize(boundary_mesh.V.size(), glm::dvec3(0.0, 0.0, 0.0));
    for (int i = 0; i < boundary_mesh.V.size(); i++) {
        if (!boundary_mesh.V.at(i).isBoundary) {
            continue;
        }
        original_vertices.at(i).x = boundary_mesh.V.at(i).x;
        original_vertices.at(i).y = boundary_mesh.V.at(i).y;
        original_vertices.at(i).z = boundary_mesh.V.at(i).z;
        original_vertices.at(i).isSingularity = boundary_mesh.V.at(i).isSingularity;
        original_vertices.at(i).isBoundary = boundary_mesh.V.at(i).isBoundary;
        original_vertices.at(i).isCorner = boundary_mesh.V.at(i).isCorner;
        original_vertices.at(i).N_Vids.resize(boundary_mesh.V.at(i).N_Vids.size());
        for (int j = 0; j < boundary_mesh.V.at(i).N_Vids.size(); j++) {
            original_vertices.at(i).N_Vids.at(j) = boundary_mesh.V.at(i).N_Vids.at(j);
        }
    }
    input_vertices.resize(mesh.V.size(), glm::dvec3(0.0, 0.0, 0.0));
    for (int i = 0; i < mesh.V.size(); i++) {
        if (!mesh.V.at(i).isBoundary) {
            continue;
        }
        input_vertices.at(i).x = mesh.V.at(i).x;
        input_vertices.at(i).y = mesh.V.at(i).y;
        input_vertices.at(i).z = mesh.V.at(i).z;
        input_vertices.at(i).isSingularity = mesh.V.at(i).isSingularity;
        input_vertices.at(i).isBoundary = mesh.V.at(i).isBoundary;
        input_vertices.at(i).isCorner = mesh.V.at(i).isCorner;
        input_vertices.at(i).N_Vids.resize(mesh.V.at(i).N_Vids.size());
        for (int j = 0; j < mesh.V.at(i).N_Vids.size(); j++) {
            input_vertices.at(i).N_Vids.at(j) = mesh.V.at(i).N_Vids.at(j);
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

void SmoothAlgorithm::extractBoundary() {
    for (auto& v: boundary_mesh.V) {
        if (v.isBoundary) {
            v.isVisited = false;
        }
    }
    for (auto& v: boundary_mesh.V) {
        if (v.isBoundary && !v.isVisited) {
            std::vector<int> boundary_patch;
            int current_v_id = v.id;
            while (true) {
                Vertex& current_v = boundary_mesh.V.at(current_v_id);
                current_v.isVisited = true;
                boundary_patch.push_back(current_v_id);
                bool foundVertex = false;
                for (auto id: current_v.N_Vids) {
                    if (boundary_mesh.V.at(id).isBoundary && !boundary_mesh.V.at(id).isVisited) {
                        current_v_id = id;
                        foundVertex = true;
                        break;
                    }
                }
                if (!foundVertex) {
                    original_boundary.push_back(boundary_patch);
                    break;
                }
            }
        }
    }

    for (auto& v: mesh.V) {
        if (v.isBoundary) {
            v.isVisited = false;
        }
    }
    for (auto& v: mesh.V) {
        if (v.isBoundary && !v.isVisited) {
            std::vector<int> boundary_patch;
            int current_v_id = v.id;
            while (true) {
                Vertex& current_v = mesh.V.at(current_v_id);
                current_v.isVisited = true;
                boundary_patch.push_back(current_v_id);
                bool foundVertex = false;
                for (auto id: current_v.N_Vids) {
                    if (mesh.V.at(id).isBoundary && !mesh.V.at(id).isVisited) {
                        current_v_id = id;
                        foundVertex = true;
                        break;
                    }
                }
                if (!foundVertex) {
                    input_boundary.push_back(boundary_patch);
                    break;
                }
            }
        }
    }
}

void SmoothAlgorithm::remapToOriginalBoundary() {
    for (auto patch: input_boundary) {
        delta_coords.clear();
        delta_coords.resize(mesh.V.size(), glm::dvec3(0.0, 0.0, 0.0));
        bool performRemap = true;
        for (auto id: patch) {
            if (mesh.V.at(id).isMovable) {
                glm::dvec3 new_v(mesh.V.at(id).x, mesh.V.at(id).y, 0.0);
                delta_coords.at(id).x = new_v.x;
                delta_coords.at(id).y = new_v.y;
                double min_length = std::numeric_limits<double>::max();
                int min_index = id;
                for (int j = 0; j < original_vertices.size(); j++) {
                    if (original_vertices.at(j).isBoundary) {
                        double length = glm::length(glm::vec3(new_v.x - original_vertices.at(j).x, new_v.y - original_vertices.at(j).y, new_v.z - original_vertices.at(j).z));
                        if (length < min_length) {
                            min_length = length;
                            min_index = j;
                        }
                    }
                }
                if (min_index == id) {
                    continue;
                }
                Vertex& v = original_vertices.at(min_index);
                std::vector<size_t> boundary_neighbors;
                for (int j = 0; j < v.N_Vids.size(); j++) {
                    if (original_vertices.at(v.N_Vids.at(j)).isBoundary) {
                        boundary_neighbors.push_back(v.N_Vids.at(j));
                    }                    
                }
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

                double new_x = 0.0;
                double new_y = 0.0;
                if (dot_a_b > 0) {
                    b = glm::normalize(b);
                    delta_coords.at(id).x = v.x + (dot_a_b * b.x);
                    delta_coords.at(id).y = v.y + (dot_a_b * b.y);
                } else if (dot_a_c > 0){
                    c = glm::normalize(c);
                    delta_coords.at(id).x = v.x + (dot_a_c * c.x);
                    delta_coords.at(id).y = v.y + (dot_a_c * c.y);
                } else {
                    delta_coords.at(id).x = new_v.x;
                    delta_coords.at(id).y = new_v.y;
                }

                for (auto fid: v.N_Fids) {
                    if (isFaceNegative(fid, v.id, glm::dvec3(delta_coords.at(id).x, delta_coords.at(id).y, 0.0))) {
                        performRemap =false;
                        break;
                    }
                }
                if (!performRemap) {
                    break;
                }
            }
        }
        if (!performRemap) {
            continue;
        }
        for (auto id: patch) {
            if (mesh.V.at(id).isMovable) {
                mesh.V.at(id).x = delta_coords.at(id).x;
                mesh.V.at(id).y = delta_coords.at(id).y;
            }
        }
    }
}

void SmoothAlgorithm::resampleBoundaryVertices() {
    bool hasNegativeElementsPresentAlready = false;
    for (auto& v: mesh.V) {
        if (!v.isBoundary || !v.isMovable) {
            continue;
        }
        for (auto fid: v.N_Fids) {
            if (isFaceNegative(fid, v.id, glm::dvec3(v.x, v.y, 0.0))) {
                hasNegativeElementsPresentAlready = true;
                break;
            }
        }
    }
    for (int i = 0; i < mesh.V.size(); i++) {
        Vertex& v = mesh.V.at(i);
        if (v.isBoundary) {
            if (!v.isMovable) {
                delta_coords.at(i).x = v.x;
                delta_coords.at(i).y = v.y;
                continue;
            }
            int nb_n_index = -1;
            int k = 0;
            std::vector<size_t> neighbors = v.N_Vids;
            if (v.N_Fids.size() < 2) {
                neighbors = v.oneRingNeighborVertices;
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

                    glm::dvec3 V_j(v.x - b0.x, v.y - b0.y, v.z - b0.z);
                    glm::dvec3 V_j_minus_1(b1.x - b0.x, b1.y - b0.y, b1.z - b0.z);
                    glm::dvec3 V_j_plus_1(b2.x - b0.x, b2.y - b0.y, b2.z - b0.z);

                    double a1 = glm::dot(V_j, V_j_plus_1) / (glm::length(V_j) * glm::length(V_j_plus_1));
                    double a2 = glm::dot(V_j, V_j_minus_1) / (glm::length(V_j) * glm::length(V_j_minus_1));
                    if (a1 < -1.0) a1 = -1.0;
                    if (a1 > 1.0) a1 = 1.0;
                    if (a2 < -1.0) a2 = -1.0;
                    if (a2 > 1.0) a2 = 1.0;

                    double alpha1 = acos(a1);
                    double alpha2 = acos(a2);

                    // double alpha1 = acos(glm::dot(V_j, V_j_plus_1) / (glm::length(V_j) * glm::length(V_j_plus_1)));
                    // double alpha2 = acos(glm::dot(V_j, V_j_minus_1) / (glm::length(V_j) * glm::length(V_j_minus_1)));
                    double beta = (alpha2 - alpha1) / 2;
                    double new_x = (b0.x + (V_j.x * cos(beta)) - (V_j.y * sin(beta)));
                    double new_y = (b0.y + (V_j.x * sin(beta)) + (V_j.y * cos(beta)));
                    // if (isnan(new_x) || isnan(new_y)) {
                    //     std::cout << "vertex: " << mesh.V.at(i).id << " is faulty" << std::endl;
                    //     std::cout << "beta: " << beta << std::endl;
                    // }
                    
                    // bool isNegativeElementPresent = false;
                    // for (auto fid: v.N_Fids) {
                    //     if (isFaceNegative(fid, v.id, glm::dvec3(new_x, new_y, 0.0))) {
                    //         isNegativeElementPresent = true;
                    //     }
                    // }
                    // if (!hasNegativeElementsPresentAlready && isNegativeElementPresent) {
                    //     new_x = v.x;
                    //     new_y = v.y;
                    // }
                    delta_coords.at(i).x += new_x; 
                    delta_coords.at(i).y += new_y;
                }                    
            }
            if (k > 0) {
                delta_coords.at(i).x = delta_coords.at(i).x / k; 
                delta_coords.at(i).y = delta_coords.at(i).y / k;            
            } else {
                delta_coords.at(i).x = v.x;
                delta_coords.at(i).y = v.y;
            }
            // double currentE = getVertexEnergy(v.id);
            // glm::dvec3 temp_coord(v.x, v.y, 0.0);
            // v.x = delta_coords.at(i).x;
            // v.y = delta_coords.at(i).y;
            // double newE = getVertexEnergy(v.id);
            // v.x = temp_coord.x;
            // v.y = temp_coord.y;
            // if (newE >= currentE) {
            //     delta_coords.at(i).x = v.x;
            //     delta_coords.at(i).y = v.y;
            // }
        }
    }
}

void SmoothAlgorithm::resampleBoundaryVertices1() {
    bool hasNegativeElementsPresentAlready = false;
    for (auto& v: mesh.V) {
        if (!v.isBoundary || !v.isMovable) {
            continue;
        }
        for (auto fid: v.N_Fids) {
            if (isFaceNegative(fid, v.id, glm::dvec3(v.x, v.y, 0.0))) {
                hasNegativeElementsPresentAlready = true;
                break;
            }
        }
    }
    for (int i = 0; i < mesh.V.size(); i++) {
        Vertex& v = mesh.V.at(i);
        // if (!global && !v.smoothLocal) {
        //     continue;
        // }
        if (v.isBoundary && !v.isCorner) {
            // std::vector<size_t> n_bvs;
            // if (v.N_Vids.size() > 2) {
            //     for (auto vid: v.N_Vids) {
            //         if (mesh.V.at(vid).isBoundary) {
            //             n_bvs.push_back(vid);
            //         }
            //     }
            //     Vertex& n_bv1 = mesh.V.at(n_bvs.at(0));
            //     Vertex& n_bv2 = mesh.V.at(n_bvs.at(1));
            //     glm::dvec3 v_bv1(v.x - n_bv1.x, v.y - n_bv1.y, v.z - n_bv1.z);
            //     glm::dvec3 v_bv2(v.x - n_bv2.x, v.y - n_bv2.y, v.z - n_bv2.z);

            //     double angle = (atan2(glm::cross(v_bv1, v_bv2).z, glm::dot(v_bv1, v_bv2))) * 180 / PI;
            //     if (angle < 0) {
            //         angle += 360;
            //     }
            //     if (fabs(angle - 180) > 25) {
            //         continue;
            //     }
            // }
            if (!v.isMovable) {
                delta_coords.at(i).x = v.x;
                delta_coords.at(i).y = v.y;
                continue;
            }
            int nb_n_index = -1;
            int k = 0;
            std::vector<size_t> neighbors = v.N_Vids;
            if (v.N_Fids.size() < 2) {
                neighbors = v.oneRingNeighborVertices;
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

                    glm::dvec3 V_j(v.x - b0.x, v.y - b0.y, v.z - b0.z);
                    glm::dvec3 V_j_minus_1(b1.x - b0.x, b1.y - b0.y, b1.z - b0.z);
                    glm::dvec3 V_j_plus_1(b2.x - b0.x, b2.y - b0.y, b2.z - b0.z);

                    double a1 = glm::dot(V_j, V_j_plus_1) / (glm::length(V_j) * glm::length(V_j_plus_1));
                    double a2 = glm::dot(V_j, V_j_minus_1) / (glm::length(V_j) * glm::length(V_j_minus_1));
                    if (a1 < -1.0) a1 = -1.0;
                    if (a1 > 1.0) a1 = 1.0;
                    if (a2 < -1.0) a2 = -1.0;
                    if (a2 > 1.0) a2 = 1.0;

                    double alpha1 = acos(a1);
                    double alpha2 = acos(a2);

                    // double alpha1 = acos(glm::dot(V_j, V_j_plus_1) / (glm::length(V_j) * glm::length(V_j_plus_1)));
                    // double alpha2 = acos(glm::dot(V_j, V_j_minus_1) / (glm::length(V_j) * glm::length(V_j_minus_1)));
                    double beta = (alpha2 - alpha1) / 2;

                    glm::dvec3 r1(v.x, v.y, v.z);
                    double l = 0;
                    if (beta > 0) {
                        r1 = glm::dvec3(b1.x - v.x, b1.y - v.y, b1.z - v.z);
                        l = fabs(beta / alpha2) * glm::length(r1); 
                    } else if (beta < 0) {
                        r1 = glm::dvec3(b2.x - v.x, b2.y - v.y, b2.z - v.z);
                        l = fabs(beta / alpha1) * glm::length(r1);
                    }
                    r1 = glm::normalize(r1);
                    double new_x = v.x + (l * r1.x);
                    double new_y = v.y + (l * r1.y);

                    // double new_x = (b0.x + (V_j.x * cos(beta)) - (V_j.y * sin(beta)));
                    // double new_y = (b0.y + (V_j.x * sin(beta)) + (V_j.y * cos(beta)));
                    // if (isnan(new_x) || isnan(new_y)) {
                    //     std::cout << "vertex: " << mesh.V.at(i).id << " is faulty" << std::endl;
                    //     std::cout << "beta: " << beta << std::endl;
                    // }
                    
                    bool isNegativeElementPresent = false;
                    for (auto fid: v.N_Fids) {
                        if (isFaceNegative(fid, v.id, glm::dvec3(new_x, new_y, 0.0))) {
                            isNegativeElementPresent = true;
                        }
                    }
                    if (!hasNegativeElementsPresentAlready && isNegativeElementPresent) {
                        new_x = v.x;
                        new_y = v.y;
                    }
                    delta_coords.at(i).x += new_x; 
                    delta_coords.at(i).y += new_y;
                }                    
            }
            if (k > 0) {
                delta_coords.at(i).x = delta_coords.at(i).x / k; 
                delta_coords.at(i).y = delta_coords.at(i).y / k;            
            } else {
                delta_coords.at(i).x = v.x;
                delta_coords.at(i).y = v.y;
            }
            double currentE = getVertexEnergy(v.id);
            glm::dvec3 temp_coord(v.x, v.y, 0.0);
            v.x = delta_coords.at(i).x;
            v.y = delta_coords.at(i).y;
            double newE = getVertexEnergy(v.id);
            v.x = temp_coord.x;
            v.y = temp_coord.y;
            if (newE >= currentE) {
                delta_coords.at(i).x = v.x;
                delta_coords.at(i).y = v.y;
            }
        }
    }
}

void SmoothAlgorithm::remapBoundaryVertices() {
    for (int i = 0; i < mesh.V.size(); i++) {
        if (mesh.V.at(i).isBoundary && mesh.V.at(i).isMovable) {
            glm::dvec3 new_v = delta_coords.at(i);
            // glm::dvec3 new_v(mesh.V.at(i).x, mesh.V.at(i).y, 0.0);
            double min_length = std::numeric_limits<double>::max();
            int min_index = -1;
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
            if (min_index == -1) {
                continue;
            }
            Vertex& v = original_vertices.at(min_index);
            std::vector<size_t> boundary_neighbors;
            for (int j = 0; j < v.N_Vids.size(); j++) {
                if (original_vertices.at(v.N_Vids.at(j)).isBoundary) {
                    boundary_neighbors.push_back(v.N_Vids.at(j));
                }                    
            }
            Vertex& b1 = original_vertices.at(boundary_neighbors.at(0));
            Vertex& b2 = original_vertices.at(boundary_neighbors.at(1));
            
            glm::vec3 a(new_v.x - v.x, new_v.y - v.y, new_v.z - v.z);
            glm::vec3 b(b1.x - v.x, b1.y - v.y, b1.z - v.z);
            glm::vec3 c(b2.x - v.x, b2.y - v.y, b2.z - v.z);

            // if ((a.x * b.y) - (a.y * b.x) <= 1e-4 || (a.x * c.y) - (a.y * c.x) <= 1e-4) {
            //     continue;
            // }

            double length_a = glm::dot(a, a);
            double length_b = glm::dot(b, b);
            double length_c = glm::dot(c, c);

            double dot_a_b = glm::dot(a, b) / glm::length(b);
            double dot_a_c = glm::dot(a, c) / glm::length(c);

            double new_x = new_v.x;
            double new_y = new_v.y;
            if (dot_a_b >= 0) {
                b = glm::normalize(b);
                new_x = v.x + (dot_a_b * b.x);
                new_y = v.y + (dot_a_b * b.y);
            } else if (dot_a_c >= 0){
                c = glm::normalize(c);
                new_x = v.x + (dot_a_c * c.x);
                new_y = v.y + (dot_a_c * c.y);
            }
            
            // delta_coords.at(i).x = v.x + new_x;
            // delta_coords.at(i).y = v.y + new_y;
            // new_x = v.x + new_x;
            // new_y = v.y + new_y;
            // mesh.V.at(i).x = v.x + new_x;
            // mesh.V.at(i).y = v.y + new_y;
            
            bool isNegativeElementPresent = false;
            for (auto fid: mesh.V.at(i).N_Fids) {
                if (isFaceNegative(fid, mesh.V.at(i).id, glm::dvec3(new_x, new_y, 0.0))) {
                    isNegativeElementPresent = true;
                }
            }
            if (isNegativeElementPresent) {
                delta_coords.at(i).x = mesh.V.at(i).x;
                delta_coords.at(i).y = mesh.V.at(i).y;
                continue;
            }

            delta_coords.at(i).x = new_x;
            delta_coords.at(i).y = new_y;
            // mesh.V.at(i).x = new_x;
            // mesh.V.at(i).y = new_y;
            double currentE = getVertexEnergy(mesh.V.at(i).id);
            glm::dvec3 temp_coord(mesh.V.at(i).x, mesh.V.at(i).y, 0.0);
            mesh.V.at(i).x = new_x;
            mesh.V.at(i).y = new_y;
            double newE = getVertexEnergy(mesh.V.at(i).id);
            mesh.V.at(i).x = temp_coord.x;
            mesh.V.at(i).y = temp_coord.y;
            if (newE >= currentE) {
                delta_coords.at(i).x = mesh.V.at(i).x;
                delta_coords.at(i).y = mesh.V.at(i).y;
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
    // for (auto& v: mesh.V) {
        // if (!v.isBoundary) {
        //     int n_b = 0;
        //     for (auto id: v.N_Vids) {
        //         if (mesh.V.at(id).isBoundary) {
        //             n_b += 1;
        //         }
        //     }
        //     if (n_b > 1) {
        //         for (auto id: v.oneRingNeighborVertices) {
        //             if (mesh.V.at(id).isBoundary) {
        //                 mesh.V.at(id).isMovable = false;
        //             }
        //         }
        //     }
        // }
    // }
    for (auto& v: mesh.V) {
        // if (v.isBoundary) {
        //     if (v.N_Fids.size() <= 2) {
        //         v.isMovable = true;
        //     } else {
        //         v.isMovable = false;
        //     }
        // }
        if (v.isBoundary && v.N_Fids.size() != 2) {
            v.isMovable = false;
        }
        if (v.isBoundary) {
            std::vector<int> boundary_neighbors;
            for (auto& id: v.N_Vids) {
                if (mesh.V.at(id).isBoundary) {
                    boundary_neighbors.push_back(id);
                }
            }
            Vertex& v1 = mesh.V.at(boundary_neighbors[0]);
            Vertex& v2 = mesh.V.at(boundary_neighbors[1]);
            glm::dvec3 a(v.x - v1.x, v.y - v1.y, v.z - v1.z);
            glm::dvec3 b(v.x - v2.x, v.y - v2.y, v.z - v2.z);
            if (glm::length(a) > 0) {
                a = glm::normalize(a);
            }
            if (glm::length(b) > 0) {
                b = glm::normalize(b);
            }
            double angle = (atan2(glm::cross(a, b).z, glm::dot(a, b))) * 180 / PI;
            if (angle < 0) {
                angle += 360;
            }
            if (fabs(180 - angle) > 60) {
                v.isMovable = false;
            }
        }
    }
}

void SmoothAlgorithm::findNegativeElements() {
    while (true) {
        std::vector<size_t> faceIds;
        bool foundCrossQuad = false;
        int nCrossQuads = 0;
        for (auto& f: mesh.F) {
            int sign = 0;
            for (int i = 0; i < f.Vids.size(); i++) {
                Vertex& a = mesh.V.at(f.Vids.at(i));
                Vertex& b = mesh.V.at(f.Vids.at((i + 1) % f.Vids.size()));
                Vertex& c = mesh.V.at(f.Vids.at((i + 2) % f.Vids.size()));
                double det = ((b.x - a.x) * (c.y - a.y)) - ((c.x - a.x) * (b.y - a.y));
                if (det > 0) {
                    sign += 1;
                } else if (det <= 0) {
                    sign -= 1;
                }
            }
            if (sign == 0) {
                f.isNegative = true;
                nCrossQuads += 1;
                fixCrossQuad(f.id);
                foundCrossQuad = true;
                faceIds.push_back(f.id);
            } else if (abs(sign) != 4) {
                f.isNegative = true;
            } else {
                f.isNegative = false;
            }
        }
        std::cout << "# cross quads: " << nCrossQuads << std::endl;
        // std::cout << "cross faces: " << faceIds.size() << std::endl;
        // std::ofstream ofs("cross_quads.vtk");
        // ofs << "# vtk DataFile Version 3.0\n"
        // 	<< "cross_quads.vtk\n"
        // 	<< "ASCII\n\n"
        // 	<< "DATASET UNSTRUCTURED_GRID\n";
        // ofs << "POINTS " << mesh.V.size() << " double\n";

        // for (auto& v: mesh.V) {
        //     ofs << v.x << " " << v.y << " " << v.z << std::endl;
        // }

        // ofs << "CELLS " << faceIds.size() << " " << faceIds.size() * 5 << std::endl;
        // for (auto fid: faceIds) {
        //     Face& f = mesh.F.at(fid);
        //     ofs << "4 " << f.Vids.at(0) << " " << f.Vids.at(1) << " " << f.Vids.at(2) << " " << f.Vids.at(3) <<  std::endl; 
        // }
        // ofs << "CELL_TYPES " << faceIds.size() << "\n";
        // for (auto fid: faceIds) {
        //     ofs << "9 " << std::endl; 
        // }
        if (!foundCrossQuad) {
            break;
        }
    }
}

void SmoothAlgorithm::fixCrossQuad(int fid) {
    Face& f = mesh.F.at(fid);
    for (auto id: f.Vids) {
        Vertex& v = mesh.V.at(id);
        glm::dvec3 current_coords(v.x, v.y, v.z);
        glm::dvec3 new_coords(v.x, v.y, v.z);
        int nEdges = 0;
        for (auto n_eid: v.N_Eids) {
            if (std::find(f.Eids.begin(), f.Eids.end(), n_eid) != f.Eids.end()) {
                continue;
            }
            nEdges += 1;
        }
        // std::cout << "vertex: " << v.id << " has # neighbor edges: " << nEdges << std::endl;
        int nIntersections = 0;
        for (auto n_eid: v.N_Eids) {
            if (std::find(f.Eids.begin(), f.Eids.end(), n_eid) != f.Eids.end()) {
                continue;
            }
            Edge& e = mesh.E.at(n_eid);
            for (auto f_eid: f.Eids) {
                Edge& e1 = mesh.E.at(f_eid);
                if (std::find(e1.Vids.begin(), e1.Vids.end(), v.id) != e1.Vids.end()) {
                    continue;
                }

                Vertex& p1 = mesh.V.at(e.Vids.at(0));
                Vertex& p2 = mesh.V.at(e.Vids.at(1));
                Vertex& p3 = mesh.V.at(e1.Vids.at(0));
                Vertex& p4 = mesh.V.at(e1.Vids.at(1));

                double a1 = p2.y - p1.y;
                double b1 = p1.x - p2.x;
                double c1 = (a1 * p1.x) + (b1 * p1.y);

                double a2 = p4.y - p3.y;
                double b2 = p3.x - p4.x;
                double c2 = (a2 * p3.x) + (b2 * p3.y);

                double determinant = (a1 * b2) - (a2 * b1);

                if (determinant == 0) {
                    continue;
                }

                glm::dvec3 intersection_point(((b2 * c1) - (b1 * c2))/ determinant, ((a1 * c2) - (a2 * c1)) / determinant, 0.0);
                double edge_length1 = glm::length(glm::dvec3(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z));
                double lenght_p1_intersection = glm::length(glm::dvec3(p1.x - intersection_point.x, p1.y - intersection_point.y, p1.z - intersection_point.z));
                double lenght_p2_intersection = glm::length(glm::dvec3(p2.x - intersection_point.x, p2.y - intersection_point.y, p2.z - intersection_point.z));
                double edge_length2 = glm::length(glm::dvec3(p4.x - p3.x, p4.y - p3.y, p4.z - p3.z));
                double lenght_p3_intersection = glm::length(glm::dvec3(p3.x - intersection_point.x, p3.y - intersection_point.y, p3.z - intersection_point.z));
                double lenght_p4_intersection = glm::length(glm::dvec3(p4.x - intersection_point.x, p4.y - intersection_point.y, p4.z - intersection_point.z));
                
                if (lenght_p1_intersection / edge_length1 >= 1 || lenght_p2_intersection / edge_length1 >= 1 ||
                    lenght_p3_intersection / edge_length2 >= 1 || lenght_p4_intersection / edge_length2 >= 1) {
                    continue;
                }
                nIntersections += 1;
                glm::dvec3 final_vec = glm::dvec3(p1.x, p1.y, p1.z) - intersection_point;
                if (p1.id == v.id) {
                    final_vec = glm::dvec3(p2.x, p2.y, p2.z) - intersection_point;
                }
                double lenght_final_vec = glm::length(final_vec) / 2;
                final_vec = intersection_point + (lenght_final_vec * glm::normalize(final_vec));
                if (e.isBoundary) {
                    new_coords = final_vec;
                    break;
                }
                if (glm::length(current_coords - new_coords) == 0) {
                    new_coords = final_vec;
                    continue;
                }
                if (glm::length(current_coords - final_vec) < glm::length(current_coords - new_coords)) {
                    new_coords = final_vec;
                }
            }
        }
        // std::cout << "vertex: " << v.id << " has # intersections: " << nIntersections << std::endl;
        if (glm::length(current_coords - new_coords) > 0) {
            v.x = new_coords.x;
            v.y = new_coords.y;
        }
    }
}


void SmoothAlgorithm::smoothMesh() {
    std::cout << "Started Mesh Smoothing" << std::endl;
    std::cout << mesh.V.size() << std::endl;
    setOriginalVertices();
    setBoundaryVerticesMovable();
    // calculateMeshAngles();
    mesh.SetOneRingNeighborhood();
    findNegativeElements();
    mesh.SetOneRingNeighborhood();
    // return;
    // remapBoundaryVertices();
    // setOriginalVertices();
    // return;
    int it = 0;
    // while (it < iters) {
        double E1 = getMeshEnergy(false);
        std::clock_t start;
        double duration;
        start = std::clock();

        while (it < iters) {
            delta_coords.clear();
            delta_coords.resize(mesh.V.size(), glm::dvec3(0.0, 0.0, 0.0));
            double currentE = getMeshEnergy(false);
            angleBasedSmoothing();
            // SmoothAngleBased();
            // smoothLaplacianCotangentBased();
            // resampleBoundaryVertices();
            for (int i = 0; i < mesh.V.size(); i++) {
                Vertex& v = mesh.V.at(i);
                if (v.isBoundary) {
                    continue;
                }
                glm::dvec3 temp_coord(v.x , v.y , 0.0);
                v.x = delta_coords.at(i).x;
                v.y = delta_coords.at(i).y;
                delta_coords.at(i) = temp_coord;
            }
            double newE = getMeshEnergy(false);
            std::cout << "it1: " << it << " oldE: " << currentE << " newE: " << newE << std::endl;
            if (currentE - newE < tau) {
                for (int i = 0; i < mesh.V.size(); i++) {
                    Vertex& v = mesh.V.at(i);
                    if (v.isBoundary) {
                        continue;
                    }
                    v.x = delta_coords.at(i).x;
                    v.y = delta_coords.at(i).y;
                }
                break;                
            }
            // if (boundarySmoothing) {
                int it2 = 0;
                while (it2 < iters) {
                    delta_coords.clear();
                    delta_coords.resize(mesh.V.size(), glm::dvec3(0.0, 0.0, 0.0));
                    currentE = getMeshEnergy(true);
                    resampleBoundaryVertices1();
                    // resampleBoundaryVertices();
                    remapBoundaryVertices();
                    for (int i = 0; i < mesh.V.size(); i++) {
                        Vertex& v = mesh.V.at(i);
                        if (!v.isBoundary || !v.isMovable) {
                            continue;
                        }
                        glm::dvec3 temp_coord(v.x , v.y , 0.0);
                        v.x = delta_coords.at(i).x;
                        v.y = delta_coords.at(i).y;
                        delta_coords.at(i) = temp_coord;
                    }
                    // remapBoundaryVertices();
                    double newE = getMeshEnergy(true);
                    // std::cout << "it2: " << it2 << " oldE: " << currentE << " newE: " << newE << std::endl;
                    if (currentE - newE < tau) {
                        // remapBoundaryVertices();
                        for (int i = 0; i < mesh.V.size(); i++) {
                            Vertex& v = mesh.V.at(i);
                            if (!v.isBoundary || !v.isMovable) {
                                continue;
                            }
                            v.x = delta_coords.at(i).x;
                            v.y = delta_coords.at(i).y;
                        }   
                        break;
                    }
                    it2++;
                }
            // }
            it++;
        }
        std::cout << "-------------------------" << std::endl;
    // }
    duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
    std::cout << "Simplification time: " << duration << " seconds" << std::endl;
    // remapBoundaryVertices();
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
        // int k = one_ring_neighbors.size();
        int k = 0;
        for (int j = 0; j < one_ring_neighbors.size(); j++) {
            int a = j;
            int b = one_ring_neighbors.size();
            int index = one_ring_neighbors.at((a % b + b) % b);
            Vertex& v_j = mesh.V.at(index);
            // if (std::find(v_i.N_Vids.begin(), v_i.N_Vids.end(), v_j.id) != v_i.N_Vids.end()) {
            //     continue;
            // }
            k += 1;
            // a = j - 2;
            // if (std::find(v_i.N_Vids.begin(), v_i.N_Vids.end(), v_j.id) == v_i.N_Vids.end()) {
                a = j - 1;
            // }
            index = one_ring_neighbors.at((a % b + b) % b);
            Vertex& v_j_prev = mesh.V.at(index);

            // a = j + 2;
            // if (std::find(v_i.N_Vids.begin(), v_i.N_Vids.end(), v_j.id) == v_i.N_Vids.end()) {
                a = j + 1;
            // }
            index = one_ring_neighbors.at((a % b + b) % b);
            Vertex& v_j_next = mesh.V.at(index);

            glm::dvec3 V_j(v_i.x - v_j.x, v_i.y - v_j.y, v_i.z - v_j.z);
            glm::dvec3 V_j_minus_1(v_j_prev.x - v_j.x, v_j_prev.y - v_j.y, v_j_prev.z - v_j.z);
            glm::dvec3 V_j_plus_1(v_j_next.x - v_j.x, v_j_next.y - v_j.y, v_j_next.z - v_j.z);

            double alpha1 = acos(glm::dot(V_j, V_j_plus_1) / (glm::length(V_j) * glm::length(V_j_plus_1)));
            double alpha2 = acos(glm::dot(V_j, V_j_minus_1) / (glm::length(V_j) * glm::length(V_j_minus_1)));

            double beta = (alpha2 - alpha1) / 2;
            double newX = (v_j.x + (V_j.x * cos(beta)) - (V_j.y * sin(beta))); 
            double newY = (v_j.y + (V_j.x * sin(beta)) + (V_j.y * cos(beta)));
            // bool isNegativeElementPresent = false;
            // for (auto fid: v_i.N_Fids) {
            //     if (isFaceNegative(fid, v_i.id, glm::dvec3(newX, newY, 0.0))) {
            //         isNegativeElementPresent = true;
            //     }
            // }
            // if (isNegativeElementPresent) {
            //     newX = v_i.x;
            //     newY = v_i.y;
            // }
            new_coords.at(i).x += newX; 
            new_coords.at(i).y += newY;
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
    bool hasNegativeElementsPresentAlready = false;
    for (auto& v: mesh.V) {
        if (v.isBoundary) {
            continue;
        }
        for (auto fid: v.N_Fids) {
            if (isFaceNegative(fid, v.id, glm::dvec3(v.x, v.y, 0.0))) {
                hasNegativeElementsPresentAlready = true;
                break;
            }
        }
    }

    for (int i = 0; i < mesh.V.size(); i++) {
        Vertex& v_i = mesh.V.at(i);
        if (v_i.isBoundary) {
            continue;
        }
        // if (!global && !v_i.smoothLocal) {
        //     continue;
        // }
        bool withOneRing = false;
        std::vector<size_t> one_ring_neighbors = v_i.N_Vids;
        if (withOneRing) {
            one_ring_neighbors = v_i.oneRingNeighborVertices;
        }
        std::vector<glm::dvec3> neighbor_coords(one_ring_neighbors.size(), glm::dvec3(0.0, 0.0, 0.0));
        std::vector<double> neighbor_weights(one_ring_neighbors.size(), 0.0);
        double weight_agg = 0.0;
        int k = one_ring_neighbors.size();
        for (int j = 0; j < one_ring_neighbors.size(); j++) {
            int a = j;
            int b = one_ring_neighbors.size();
            int index = one_ring_neighbors.at((a % b + b) % b);
            Vertex& v_j = mesh.V.at(index);

            if (withOneRing) {
                a = j - 2;
                if (std::find(v_i.N_Vids.begin(), v_i.N_Vids.end(), v_j.id) == v_i.N_Vids.end()) {
                    a = j - 1;
                }
            } else {
                a = j - 1;
            }
            index = one_ring_neighbors.at((a % b + b) % b);
            Vertex& v_j_prev = mesh.V.at(index);

            if (withOneRing) {
                a = j + 2;
                if (std::find(v_i.N_Vids.begin(), v_i.N_Vids.end(), v_j.id) == v_i.N_Vids.end()) {
                    a = j + 1;
                }
            } else {
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
            // std::cout << "neighbor " << j << " id: " << v_j.id << std::endl;
            // std::cout << "v: " << v_j.x << " " << v_j.y << " " << v_j.z << std::endl;
            // std::cout << "neighbor " << j << " weight: " << neighbor_weights.at(j) << std::endl; 
            // std::cout << "neighbor " << j << " coords: " << neighbor_coords.at(j).x << " " << neighbor_coords.at(j).y << " " << neighbor_coords.at(j).z << std::endl; 

            
            // bool isNegativeElementPresent = false;
            // for (auto fid: v_i.N_Fids) {
            //     if (isFaceNegative(fid, v_i.id, glm::dvec3(neighbor_coords.at(j).x, neighbor_coords.at(j).y, 0.0))) {
            //         isNegativeElementPresent = true;
            //         break;
            //     }
            // }
            // if (!hasNegativeElementsPresentAlready && isNegativeElementPresent) {
            //     neighbor_coords.at(j).x = v_i.x;
            //     neighbor_coords.at(j).y = v_i.y;
            //     neighbor_weights.at(j) = 1;
            // }
            weight_agg += neighbor_weights.at(j);
        }
        if (weight_agg == 0) {
            new_coords.at(i).x = 0;
            new_coords.at(i).y = 0;
            continue;
        }
        for (int j = 0; j < one_ring_neighbors.size(); j++) {
            glm::dvec3 current_v = neighbor_coords.at(j);
            double weight = (1 - neighbor_weights.at(j) / weight_agg);
            vertex_weights.at(i) += weight;
            new_coords.at(i).x += (weight * current_v.x);
            new_coords.at(i).y += (weight * current_v.y);
        }
        // std::cout << "neighbours: " << k << std::endl;
        // std::cout << "v id: " << v_i.id << std::endl;
        // std::cout << "v: " << v_i.x << " " << v_i.y << " " << v_i.z << std::endl;
        // std::cout << "new_coords: " << new_coords.at(i).x / vertex_weights.at(i) << " " << new_coords.at(i).y / vertex_weights.at(i) << " " << new_coords.at(i).z / vertex_weights.at(i) << std::endl;
        // std::cout << "vertex weights: " << vertex_weights.at(i) << std::endl;
        // std::cout << "********************************" << std::endl;
        new_coords.at(i).x = lambda * ((new_coords.at(i).x / vertex_weights.at(i)) - v_i.x);
        new_coords.at(i).y = lambda * ((new_coords.at(i).y / vertex_weights.at(i)) - v_i.y);
        // double currentE = getVertexEnergy(v_i.id);
        // glm::dvec3 temp_coord(v_i.x, v_i.y, 0.0);
        // v_i.x += new_coords.at(i).x;
        // v_i.y += new_coords.at(i).y;
        // double newE = getVertexEnergy(v_i.id);
        // v_i.x = temp_coord.x;
        // v_i.y = temp_coord.y;
        // if (newE >= currentE) {
        //      new_coords.at(i).x = 0;
        //      new_coords.at(i).y = 0;
        // }
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
            if (isnan(angle)) {
                std::cout << "f: " << f.id << std::endl;
                std::cout << "a: " << a.x << " " << a.y << " " << a.z << std::endl;
                std::cout << "b: " << b.x << " " << b.y << " " << b.z << std::endl;
                std::cout << "angle: " << angle << std::endl;
                std::cout << "---------------------------------------" << std::endl;
            }
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

double SmoothAlgorithm::getMeshEnergy(bool consider_boundary) {
    double E = 0.0;
    int nV = 1;
    tau = 0.0001;
    for (auto& v: mesh.V) {
        // if (v.isBoundary && !consider_boundary) {
        //     continue;
        // } else if (!v.isBoundary && consider_boundary) {
        //     continue;
        // }
        // if (v.isBoundary && !v.isMovable) {
        //     continue;
        // }
        // nV += 1;
        E += getVertexEnergy(v.id);
    }
    // tau = tau / nV;
    return E / nV;
}

double SmoothAlgorithm::getVertexEnergy(int vid) {
    Vertex& v = mesh.V.at(vid);
    double E = 0.0;
    bool withOneRing = false;
    std::vector<size_t> neighbors = v.N_Vids;
    if (withOneRing && !v.isBoundary) {
        neighbors = v.oneRingNeighborVertices;
    }
    for (int j = 0; j < neighbors.size(); j++) {
        int a = j;
        int b = neighbors.size();
        int index = neighbors.at((a % b + b) % b);
        Vertex& v_j = mesh.V.at(index);
        if (v.isBoundary && v_j.isBoundary) {
            continue;
        }

        if (withOneRing) {
            a = j - 2;
            if (std::find(v.N_Vids.begin(), v.N_Vids.end(), v_j.id) == v.N_Vids.end()) {
                a = j - 1;
            }
        } else {
            a = j - 1;
        }
        index = neighbors.at((a % b + b) % b);
        Vertex& v_j_prev = mesh.V.at(index);

        if (withOneRing) {
            a = j + 2;
            if (std::find(v.N_Vids.begin(), v.N_Vids.end(), v_j.id) == v.N_Vids.end()) {
                a = j + 1;
            }
        } else {
            a = j + 1;
        }
        index = neighbors.at((a % b + b) % b);
        Vertex& v_j_next = mesh.V.at(index);

        glm::dvec3 V_j = glm::normalize(glm::dvec3(v.x - v_j.x, v.y - v_j.y, v.z - v_j.z));
        glm::dvec3 V_j_minus_1 = glm::normalize(glm::dvec3(v_j_prev.x - v_j.x, v_j_prev.y - v_j.y, v_j_prev.z - v_j.z));
        glm::dvec3 V_j_plus_1 = glm::normalize(glm::dvec3(v_j_next.x - v_j.x, v_j_next.y - v_j.y, v_j_next.z - v_j.z));

        double a1 = glm::dot(V_j, V_j_plus_1) / (glm::length(V_j) * glm::length(V_j_plus_1));
        double a2 = glm::dot(V_j, V_j_minus_1) / (glm::length(V_j) * glm::length(V_j_minus_1));
        if (a1 < -1.0) a1 = -1.0;
        if (a1 > 1.0) a1 = 1.0;
        if (a2 < -1.0) a2 = -1.0;
        if (a2 > 1.0) a2 = 1.0;

        double alpha1 = acos(a1);
        double alpha2 = acos(a2);
        // if (isnan(alpha1)) alpha1 = 0;
        // if (isnan(alpha2)) alpha2 = 0;
        E += ((alpha1 * alpha1) / 2) + ((alpha2 * alpha2) / 2);

        // double dp = glm::dot(V_j_plus_1, V_j_minus_1) / (glm::length(V_j_plus_1) * glm::length(V_j_minus_1));
        // double angle = (atan2(glm::cross(V_j_plus_1, V_j_minus_1).z, glm::dot(V_j_minus_1, V_j_plus_1))) * 180 / PI;
        // if (angle < 0) {
        //     angle += 360;
        // }
        // E += ((angle * angle) / 2);
    }
    return E;
}

bool SmoothAlgorithm::isFaceNegative(int fid, int vid, glm::dvec3 false_coord) {
    Face& f = mesh.F.at(fid);
    int sign = 0;
    glm::dvec3 temp_coord(0.0, 0.0, 0.0);
    for (auto id: f.Vids) {
        if (id == vid) {
            temp_coord.x = mesh.V.at(id).x;
            temp_coord.y = mesh.V.at(id).y;
            mesh.V.at(id).x = false_coord.x;
            mesh.V.at(id).y = false_coord.y;
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
            mesh.V.at(id).x = temp_coord.x;
            mesh.V.at(id).y = temp_coord.y;
        }
    }
    if (abs(sign) == 4 || sign == 0) {
        return false;
    }
    return true;
}
