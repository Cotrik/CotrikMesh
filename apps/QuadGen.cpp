/*
 * QuadGen.cpp
 *
 *  Created on: May 29, 2017
 *      Author: cotrik
 */

#include <igl/avg_edge_length.h>
#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>
#include <igl/massmatrix.h>
#include <igl/parula.h>
#include <igl/per_corner_normals.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/principal_curvature.h>
#include <igl/read_triangle_mesh.h>
//#include <igl/viewer/Viewer.h>

#include <igl/barycenter.h>
#include <igl/frame_field_deformer.h>
#include <igl/frame_to_cross_field.h>
#include <igl/jet.h>
#include <igl/local_basis.h>
#include <igl/readDMAT.h>
#include <igl/readOBJ.h>
#include <igl/rotate_vectors.h>
#include <igl/copyleft/comiso/nrosy.h>
#include <igl/copyleft/comiso/miq.h>
#include <igl/copyleft/comiso/frame_field.h>

#include <Eigen/Core>
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Cholesky>
#include <Eigen/LU>
#include <Eigen/SVD>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#include <Eigen/SparseCholesky>

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "ArgumentManager.h"

#include <iostream>

void WriteFrameField(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::MatrixXd& PD1, const Eigen::MatrixXd& PD2, const std::string& output_filename)
{
    const double avg = igl::avg_edge_length(V, F);
    // Make Edge frames
    std::vector<Vertex> frameV(V.rows() * 3, Vertex());
    std::vector<Edge> frameE(V.rows() * 2, Edge());
    for (size_t i = 0; i < V.rows(); ++i) {
        Vertex& v = frameV.at(i);
        v.x = V(i, 0);
        v.y = V(i, 1);
        v.z = V(i, 2);
        v.id = i;
    }
    Eigen::MatrixXd VD1 = V + PD1 * avg;
    Eigen::MatrixXd VD2 = V + PD2 * avg;

    for (size_t i = 0; i < V.rows(); ++i) {
        const size_t vid = i + V.rows();
        Vertex& v = frameV.at(vid);
        v.x = VD1(i, 0);
        v.y = VD1(i, 1);
        v.z = VD1(i, 2);
        v.id = vid;
    }

    for (size_t i = 0; i < V.rows(); ++i) {
        const size_t vid = i + 2 * V.rows();
        Vertex& v = frameV.at(vid);
        v.x = VD2(i, 0);
        v.y = VD2(i, 1);
        v.z = VD2(i, 2);
        v.id = vid;
    }

    for (size_t i = 0; i < V.rows(); ++i) {
        Edge& edge1 = frameE.at(i);
        Edge& edge2 = frameE.at(i + V.rows());
        edge1.Vids = std::vector<size_t> {i, i + V.rows()};
        edge2.Vids = std::vector<size_t> {i, i + V.rows() * 2};
    }
    Mesh newMesh;
    newMesh.V = frameV;
    newMesh.E = frameE;
    MeshFileWriter writer(newMesh, output_filename.c_str());
    writer.WriteFramesVtk();
}
int main(int argc, char *argv[])
{
    if (argc < 3) {
        std::cout << "Usage: QuadGen <input_tri_file> <output_tri_vtk_file>" << std::endl;
        return -1;
    }
    ArgumentManager argumentManager(argc, argv);
    std::string filename = argv[1];
    std::string output_filename = argv[2];

    Eigen::MatrixXd V;  // Vertex Coordinates of Triangle mesh;
    Eigen::MatrixXi F;  // Faces of Triangle mesh

    igl::read_triangle_mesh(filename, V, F);

    // Alternative discrete mean curvature
    Eigen::MatrixXd HN;
    Eigen::SparseMatrix<double> L, M, Minv;
    igl::cotmatrix(V, F, L);
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M);
    igl::invert_diag(M, Minv);
    HN = -Minv * (L * V); // Laplace-Beltrami of position
    Eigen::VectorXd H = HN.rowwise().norm(); // Extract magnitude as mean curvature

    // Compute curvature directions via quadric fitting
    Eigen::MatrixXd PD1, PD2;
    Eigen::VectorXd PV1, PV2;
    igl::principal_curvature(V, F, PD1, PD2, PV1, PV2);
    H = 0.5 * (PV1 + PV2); // mean curvature

    Eigen::MatrixXd C;
    igl::parula(H, true, C); // Compute pseudocolor

    WriteFrameField(V, F, PV1, PV2, output_filename);


    Eigen::MatrixXd B; // Face barycenters

    // Scale for visualizing the fields
    double global_scale;

    // Input frame field constraints
    Eigen::VectorXi b;
    Eigen::MatrixXd bc1 = PD1;
    Eigen::MatrixXd bc2 = PD2;

    // Interpolated frame field
    Eigen::MatrixXd FF1, FF2;

    // Deformed mesh
    Eigen::MatrixXd V_deformed;
    Eigen::MatrixXd B_deformed;

    // Frame field on deformed
    Eigen::MatrixXd FF1_deformed;
    Eigen::MatrixXd FF2_deformed;

    // Cross field on deformed
    Eigen::MatrixXd X1_deformed;
    Eigen::MatrixXd X2_deformed;

    // Global parametrization
    Eigen::MatrixXd V_uv;
    Eigen::MatrixXi F_uv;


    // Interpolate the frame field
    igl::copyleft::comiso::frame_field(V, F, b, bc1, bc2, FF1, FF2);

    // Deform the mesh to transform the frame field in a cross field
    igl::frame_field_deformer(V, F, FF1, FF2, V_deformed, FF1_deformed, FF2_deformed);

    // Compute face barycenters deformed mesh
    igl::barycenter(V_deformed, F, B_deformed);

    // Find the closest crossfield to the deformed frame field
    igl::frame_to_cross_field(V_deformed, F, FF1_deformed, FF2_deformed, X1_deformed);

    // Find a smooth crossfield that interpolates the deformed constraints
    Eigen::MatrixXd bc_x(b.size(), 3);
    for (unsigned i = 0; i < b.size(); ++i)
        bc_x.row(i) = X1_deformed.row(b(i));

    Eigen::VectorXd S;
    igl::copyleft::comiso::nrosy(V, F, b, bc_x, Eigen::VectorXi(), Eigen::VectorXd(), Eigen::MatrixXd(), 4, 0.5, X1_deformed, S);

    // The other representative of the cross field is simply rotated by 90 degrees
    Eigen::MatrixXd B1, B2, B3;
    igl::local_basis(V_deformed, F, B1, B2, B3);
    X2_deformed = igl::rotate_vectors(X1_deformed, Eigen::VectorXd::Constant(1, M_PI / 2), B1, B2);

    // Global seamless parametrization
    igl::copyleft::comiso::miq(V_deformed, F, X1_deformed, X2_deformed, V_uv, F_uv, 60.0, 5.0, false, 2);

    return 0;
}



