#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "ArgumentManager.h"
#include "AngleBasedSmoothQuadMesh.h"
#include <iostream>
#include <iomanip>
#include <memory>
#include <algorithm>
#include <math.h>
#include <glm/glm.hpp>

#define PI 3.14159265

/*int main(int argc, char* argv[]) {
    // std::cout << atan2(-1, 0) * 180 / PI << std::endl;
    glm::dvec2 source(6.0, 6.0);
    std::cout << "source: " << source.x << " " << source.y << std::endl;
    std::vector<glm::dvec2> one_ring_neighbors {
        glm::dvec2(10.0, 10.0),
        glm::dvec2(8.0, 10.0),
        glm::dvec2(4.0, 8.0),
        glm::dvec2(1.0, 7.0),
        glm::dvec2(1.0, 4.0),
        glm::dvec2(4.0, 1.0),
        glm::dvec2(8.0, 4.0),
        glm::dvec2(10.0, 5.0)
    };
    glm::dvec2 a(one_ring_neighbors.at(0).x - source.x, one_ring_neighbors.at(0).y - source.y);
    std::cout << "a: " << a.x << " " << a.y << std::endl;
    a = glm::normalize(a);
    for (int i = 0; i < one_ring_neighbors.size(); i++) {
        glm::dvec2 b(one_ring_neighbors.at(i).x - source.x, one_ring_neighbors.at(i).y - source.y);
        std::cout << "b: " << b.x << " " << b.y << std::endl;
        b = glm::normalize(b);
        double cross = (a.x * b.y) - (a.y * b.x);
        double dot = (a.x * b.x) + (a.y * b.y);
        std::cout << "cross a, b: " << cross << std::endl;
        std::cout << "dot a, b: " << dot << std::endl;
        double angle = (atan2(cross, dot))  * 180 / PI;
        if (angle < 0) {
            angle += 360;
        }
        std::cout << "neighbor " << i << " : " << angle << std::endl;
        std::cout << "----------------------------------------------" << std::endl;
    }
    /*double x_agg = 0;
    double y_agg = 0;
    for (int i = 0; i < one_ring_neighbors.size(); i += 2) {
        int a = i - 1;
        int b = one_ring_neighbors.size();
        glm::dvec2 v1 = one_ring_neighbors.at((a % b + b) % b);
        a = i;
        glm::dvec2 v2 = one_ring_neighbors.at((a % b + b) % b);
        a = i + 1;
        glm::dvec2 v3 = one_ring_neighbors.at((a % b + b) % b);

        // std::cout << "v1: " << v1.x << " " << v1.y << std::endl;
        // std::cout << "v2: " << v2.x << " " << v2.y << std::endl;
        // std::cout << "v3: " << v3.x << " " << v3.y << std::endl;

        glm::dvec2 v_j = glm::normalize(glm::dvec2(v0.x - v2.x, v0.y - v2.y));
        glm::dvec2 v_j_minus_1 = glm::normalize(glm::dvec2(v1.x - v2.x, v1.y - v2.y));
        glm::dvec2 v_j_plus_1 = glm::normalize(glm::dvec2(v3.x - v2.x, v3.y - v2.y));

        // std::cout << "v_j: " << v_j.x << " " << v_j.y << std::endl;
        // std::cout << "v_j_minus_1: " << v_j_minus_1.x << " " << v_j_minus_1.y << std::endl;
        // std::cout << "v_j_plus_1: " << v_j_plus_1.x << " " << v_j_plus_1.y << std::endl;

        double v_j_dot_v_j_plus_1 = glm::dot(v_j, v_j_plus_1) / glm::length(v_j) * glm::length(v_j_plus_1);
        double v_j_dot_v_j_minus_1 = glm::dot(v_j, v_j_minus_1) / glm::length(v_j) * glm::length(v_j_minus_1);

        // std::cout << "v_j dot v_j_plus_1: " << v_j_dot_v_j_plus_1 << std::endl; 
        // std::cout << "v_j dot v_j_minus_1: " << v_j_dot_v_j_minus_1 << std::endl; 

        double alpha1 = acos(v_j_dot_v_j_plus_1);
        double alpha2 = acos(v_j_dot_v_j_minus_1);

        // double beta = (180 * PI / 180) - (alpha1 + alpha2);
        double beta = (alpha2 - alpha1) / 2;
        std::cout << "beta: " << beta << std::endl;

        x_agg += (v2.x + (v_j.x * cos(beta)) - (v_j.y * sin(beta)));
        y_agg += (v2.y + (v_j.x * sin(beta)) + (v_j.y * cos(beta)));
        std::cout << x_agg << " " << y_agg << std::endl;
        std::cout << "------------------------" << std::endl;
    }

    for (int i = 1; i < one_ring_neighbors.size(); i += 2) {
        int a = i - 1;
        int b = one_ring_neighbors.size();
        glm::dvec2 v1 = one_ring_neighbors.at((a % b + b) % b);
        a = i;
        glm::dvec2 v2 = one_ring_neighbors.at((a % b + b) % b);
        a = i + 1;
        glm::dvec2 v3 = one_ring_neighbors.at((a % b + b) % b);

        glm::dvec2 v_j = glm::normalize(glm::dvec2(v0.x - v2.x, v0.y - v2.y));
        glm::dvec2 v_j_minus_1 = glm::normalize(glm::dvec2(v1.x - v2.x, v1.y - v2.y));
        glm::dvec2 v_j_plus_1 = glm::normalize(glm::dvec2(v3.x - v2.x, v3.y - v2.y));

        double v_j_dot_v_j_plus_1 = glm::dot(v_j, v_j_plus_1) / glm::length(v_j) * glm::length(v_j_plus_1);
        double v_j_dot_v_j_minus_1 = glm::dot(v_j, v_j_minus_1) / glm::length(v_j) * glm::length(v_j_minus_1);

        double alpha1 = acos(v_j_dot_v_j_plus_1);
        double alpha2 = acos(v_j_dot_v_j_minus_1);

        // double beta = (90 * PI / 180) - (alpha1 + alpha2);
        double beta = (alpha2 - alpha1) / 2;
        std::cout << "beta: " << beta << std::endl;

        x_agg += (v2.x + (v_j.x * cos(beta)) - (v_j.y * sin(beta)));
        y_agg += (v2.y + (v_j.x * sin(beta)) + (v_j.y * cos(beta)));
        std::cout << x_agg << " " << y_agg << std::endl;
        std::cout << "------------------------" << std::endl;
    }

    x_agg = x_agg / 10;
    y_agg = y_agg / 10;

    std::cout << x_agg << " " << y_agg << std::endl;
    return 0;
}*/

int main(int argc, char* argv[]) {
    if (argc < 5) {
        std::cout << "Usage: SmoothMesh <input_vtk_file> <output_vtk_file> smoothing_algorithm=<LAPLACIAN,ANGLE_BASED> ";
        std::cout << "iterations=<number_of_iterations> lambda=<value>" << std::endl;
        return -1;
    }

    ArgumentManager argumentManager(argc, argv);
    std::string filename = argv[1];
    std::string output_filename = argv[2];
    std::string smoothing_algorithm = argumentManager.get("smoothing_algorithm");
    int iterations;
    double lambda;
    argumentManager.get("iterations", iterations);
    argumentManager.get("lambda", lambda);

    std::cout << "---------------------------------------" << std::endl;
    std::cout << "filename = " << filename << std::endl;
    std::cout << "output_filename = " << output_filename << std::endl;
    std::cout << "smoothing_algorithm = " << smoothing_algorithm << std::endl;
    std::cout << "iterations = " << iterations << std::endl;
    std::cout << "lambda = " << lambda << std::endl;
    std::cout << "---------------------------------------" << std::endl;

    MeshFileReader reader(filename.c_str());
    Mesh& mesh = (Mesh&)reader.GetMesh();
    mesh.RemoveUselessVertices();
    mesh.BuildAllConnectivities();
    mesh.ExtractBoundary();
    mesh.ExtractSingularities();
    mesh.unifyOrientation();
    mesh.SetOneRingNeighborhood();

    for (auto& v: mesh.V) {
        v.isVisited = false;
    }

    std::vector<size_t> v_;
    if (smoothing_algorithm.compare("LAPLACIAN") == 0) {
        SmoothAlgorithm algorithm(mesh, iterations, lambda);
        // algorithm.setOriginalVertices();
        // algorithm.resampleBoundaryVertices();
        // algorithm.smoothBoundary();
        // algorithm.smoothLaplacianSimple();
        // algorithm.smoothLaplacianScaleBased();
        // algorithm.smoothLaplacianCotangentBased();
        // algorithm.SmoothAngleBased();
        // algorithm.angleBasedSmoothing();
        algorithm.smoothMesh();
        // algorithm.smoothInteriorVertices();
        // algorithm.smoothDiagonalVertices();
        // for (auto& v: mesh.V) {
        //     if (!v.isBoundary && v.N_Vids.size() == 5) {
        //         v_ = v.oneRingNeighborVertices;
        //         // break;
        //     }
        // }
    }

    // std::cout << "Writing output file" << std::endl;
    // std::ofstream ofs("output.vtk");
    // ofs << "# vtk DataFile Version 3.0\n"
    //     << "output.vtk\n"
    //     << "ASCII\n\n"
    //     << "DATASET UNSTRUCTURED_GRID\n";
    // ofs << "POINTS " << mesh.V.size() << " double\n";
    // std::vector<size_t> c_indices = v_;
    // // // std::vector<size_t> c_indices = mesh.F.at(v_.at(iterations)).Vids;
    // for (size_t i = 0; i < mesh.V.size(); i++) {
    //     ofs << std::fixed << std::setprecision(7) <<  mesh.V.at(i).x << " " <<  mesh.V.at(i).y << " " <<  mesh.V.at(i).z << "\n";
    //     // if (mesh.V.at(i).isVisited) {
    //     //     std::cout << mesh.V.at(i).id << " is visited" << std::endl;
    //     //     c_indices.push_back(i);
    //     // }
    // }
    // ofs << "CELLS " << c_indices.size() << " " << 2 * c_indices.size() << std::endl;
    // for (size_t i = 0; i < c_indices.size(); i++) {
    //     ofs << "1 " << c_indices.at(i) << std::endl;
    // }
    // ofs << "CELL_TYPES " << c_indices.size() << "\n";
    // for (size_t i = 0; i < c_indices.size(); i++) {
    //     ofs << "1" << std::endl;
    // }

    // ofs << "CELLS " << 1 << " " << 5  << std::endl;
    // ofs << "4 " << c_indices.at(0) << " " << c_indices.at(1) << " " << c_indices.at(2) << " " << c_indices.at(3) << std::endl;
    // ofs << "CELL_TYPES " << 1 << "\n";
    // ofs << "9" << std::endl;

    // ofs << "CELLS " << 1 << " " << 2  << std::endl;
    // ofs << "1 " << c_indices.at(iterations) << std::endl;
    // ofs << "CELL_TYPES " << 1 << "\n";
    // ofs << "1" << std::endl;

    // std::vector<size_t> c_indices = v_;
    // for (int i = 0; i < v_.size(); i++) {
    //      std::ofstream ofs(("output_" + std::to_string(i) + ".vtk").c_str());
    //      ofs << "# vtk DataFile Version 3.0\n"
    //         << "output.vtk\n"
    //         << "ASCII\n\n"
    //         << "DATASET UNSTRUCTURED_GRID\n";
    //     ofs << "POINTS " << mesh.V.size() << " double\n";
    //     for (size_t i = 0; i < mesh.V.size(); i++) {
    //         ofs << std::fixed << std::setprecision(7) <<  mesh.V.at(i).x << " " <<  mesh.V.at(i).y << " " <<  mesh.V.at(i).z << "\n";
    //     }
    //     ofs << "CELLS " << 1 << " " << 2  << std::endl;
    //     ofs << "1 " << c_indices.at(i) << std::endl;
    //     ofs << "CELL_TYPES " << 1 << "\n";
    //     ofs << "1" << std::endl;
    // }

    MeshFileWriter writer(mesh, output_filename.c_str());
    writer.WriteFile();
    return 0;
}