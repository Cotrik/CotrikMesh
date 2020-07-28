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

void resetMesh(Mesh& mesh) {
    std::vector<size_t> FaceIds;
	FaceIds.reserve(mesh.F.size());
	for (auto& f : mesh.F)
		// if (canceledFids.find(f.id) == canceledFids.end())
			FaceIds.push_back(f.id);
	std::vector<Vertex> newV(mesh.V.size());
//	std::vector<Face> newF(FaceIds.size());
	std::vector<Cell> newC(FaceIds.size());
	for (size_t i = 0; i < mesh.V.size(); ++i) {
	    auto& v = mesh.V.at(i);
	    auto& newv = newV.at(i);
		newv.id = i;
		newv = v.xyz();
		newv.type = v.type;
		newv.isCorner = v.isCorner;
		newv.isConvex = v.isConvex;
		newv.label = v.label;
		newv.patch_id = v.patch_id;
		newv.isSpecial = v.isSpecial;
		newv.labels = v.labels;
		newv.patch_ids = v.patch_ids;
		newv.idealValence = v.idealValence;
	}
//	for (size_t i = 0; i < FaceIds.size(); ++i) {
//		newF.at(i).id = i;
//		newF.at(i).Vids = mesh.F.at(FaceIds[i]).Vids;
//	}
	for (size_t i = 0; i < FaceIds.size(); ++i) {
		newC.at(i).id = i;
		newC.at(i).Vids = mesh.F.at(FaceIds[i]).Vids;
		newC.at(i).cellType = VTK_QUAD;
	}
	mesh.V.clear();
	mesh.E.clear();
	mesh.F.clear();
	mesh.C.clear();

	mesh.V = newV;
//	mesh.F = newF;
	mesh.C = newC;
}

int main(int argc, char* argv[]) {
    if (argc < 6) {
        std::cout << "Usage: SmoothMesh <input_vtk_file> <boundary_file.vtk> <output_vtk_file> smoothing_algorithm=<LAPLACIAN,ANGLE_BASED> ";
        std::cout << "iterations=<number_of_iterations> lambda=<value>" << std::endl;
        return -1;
    }

    ArgumentManager argumentManager(argc, argv);
    std::string filename = argv[1];
    std::string filename2 = argv[2];
    std::string output_filename = argv[3];
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
    mesh.BuildParallelE();
    // mesh.unifyOrientation();
    mesh.SetOneRingNeighborhood();
    // std::cout << "Initialized mesh" << std::endl;
    for (auto& v: mesh.V) {
        v.isVisited = false;
    }

    MeshFileReader reader2(filename2.c_str());
    Mesh& boundary_mesh = (Mesh&)reader2.GetMesh();
    boundary_mesh.RemoveUselessVertices();
    boundary_mesh.BuildAllConnectivities();
    boundary_mesh.ExtractBoundary();
    boundary_mesh.ExtractSingularities();

    std::vector<size_t> v_;
    if (smoothing_algorithm.compare("LAPLACIAN") == 0) {
        SmoothAlgorithm algorithm(mesh, boundary_mesh, iterations, lambda);
        // algorithm.smoothLaplacianSimple();
        // algorithm.smoothLaplacianScaleBased();
        // algorithm.smoothLaplacianCotangentBased();
        std::cout << "Smoothing" << std::endl;
        // algorithm.setOriginalVertices();
        // algorithm.extractBoundary();
        // algorithm.findNegativeElements();
        // mesh.SetOneRingNeighborhood();
        // algorithm.remapToOriginalBoundary();
        algorithm.smoothMesh();
        for (auto& f: mesh.F) {
            if (f.isNegative) {
                v_.push_back(f.id);
            }
        }
        // for (auto& v: mesh.V) {
        //     if (v.isBoundary && v.N_Fids.size() > 2) {
        //         // v_ = v.oneRingNeighborVertices;
        //         v_.push_back(v.id);
        //         // break;
        //     }
        // }
    }

    std::cout << "Writing output file" << std::endl;
    std::ofstream ofs("output.vtk");
    ofs << "# vtk DataFile Version 3.0\n"
        << "output.vtk\n"
        << "ASCII\n\n"
        << "DATASET UNSTRUCTURED_GRID\n";
    ofs << "POINTS " << mesh.V.size() << " double\n";
    // std::vector<size_t> c_indices = v_;
    std::vector<size_t> c_indices = {12, 296};
    std::cout << c_indices.size() << std::endl;
    for (size_t i = 0; i < mesh.V.size(); i++) {
        ofs << std::fixed << std::setprecision(7) <<  mesh.V.at(i).x << " " <<  mesh.V.at(i).y << " " <<  mesh.V.at(i).z << "\n";
    }
    // ofs << "CELLS " << c_indices.size() << " " << 2 * c_indices.size() << std::endl;
    // for (size_t i = 0; i < c_indices.size(); i++) {
    //     ofs << "1 " << c_indices.at(i) << std::endl;
    // }
    // ofs << "CELL_TYPES " << c_indices.size() << "\n";
    // for (size_t i = 0; i < c_indices.size(); i++) {
    //     ofs << "1" << std::endl;
    // }

    ofs << "CELLS " << c_indices.size() << " " << 5 * c_indices.size() << std::endl;
    for (size_t i = 0; i < c_indices.size(); i++) {
        Face& f = mesh.F.at(c_indices.at(i));
        ofs << "4 " << f.Vids.at(0) << " " << f.Vids.at(1) << " " << f.Vids.at(2) << " " << f.Vids.at(3) << std::endl;
    }
    ofs << "CELL_TYPES " << c_indices.size() << "\n";
    for (size_t i = 0; i < c_indices.size(); i++) {
        ofs << "9" << std::endl;
    }

    MeshFileWriter writer(mesh, output_filename.c_str());
    writer.WriteFile();
    return 0;
}