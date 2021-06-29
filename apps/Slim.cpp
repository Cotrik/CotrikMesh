/*
 * SLIM.cpp
 *
 *  Created on: Dec 15, 2016
 *      Author: cotrik
 */

#include <iostream>

#include <igl/slim.h>

#include <igl/components.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/Timer.h>

#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/harmonic.h>
#include <igl/serialize.h>
#include <igl/read_triangle_mesh.h>
//#include <igl/viewer/Viewer.h>
#include <igl/flipped_triangles.h>
#include <igl/euler_characteristic.h>
#include <igl/barycenter.h>
#include <igl/adjacency_matrix.h>
#include <igl/is_edge_manifold.h>
#include <igl/doublearea.h>
#include <igl/cat.h>

#include "Mesh.h"
#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include <stdlib.h>

#include <string>
#include <vector>

using namespace std;
using namespace Eigen;

#ifndef TUTORIAL_SHARED_PATH
#define TUTORIAL_SHARED_PATH "../../libigl/tutorial/shared"
#endif

void check_mesh_for_issues(Eigen::MatrixXd& V, Eigen::MatrixXi& F);
void param_2d_demo_iter(/*igl::viewer::Viewer& viewer*/);
void get_soft_constraint_for_circle(Eigen::MatrixXd& V_o, Eigen::MatrixXi& F, Eigen::VectorXi& b, Eigen::MatrixXd& bc);
void soft_const_demo_iter(/*igl::viewer::Viewer& viewer*/);
void deform_3d_demo_iter(/*igl::viewer::Viewer& viewer*/);
void get_cube_corner_constraints(Eigen::MatrixXd& V_o, Eigen::MatrixXi& F, Eigen::VectorXi& b, Eigen::MatrixXd& bc);
void display_3d_mesh(/*igl::viewer::Viewer& viewer*/);
void int_set_to_eigen_vector(const std::set<int>& int_set, Eigen::VectorXi& vec);

Eigen::MatrixXd V;
Eigen::MatrixXi F;
bool first_iter = true;
igl::SLIMData sData;
igl::Timer timer;

double uv_scale_param;

enum DEMO_TYPE
{
    PARAM_2D, SOFT_CONST, DEFORM_3D
};
DEMO_TYPE demo_type;

void GetVertices(const Eigen::MatrixXd& V_o, std::vector<Vertex>& V)
{
    V.resize(V_o.rows());
    for (size_t i = 0; i < V_o.rows(); i++) {
        Vertex& v = V.at(i);
        v.id = i;
        for (size_t j = 0; j < V_o.cols(); j++)
            v[j] = V_o(i, j);
    }
}
void GetCells(const Eigen::MatrixXi& F, std::vector<Cell>& C)
{
    C.resize(F.rows());
    for (size_t i = 0; i < F.rows(); i++) {
        Cell& c = C.at(i);
        c.id = i;
        c.Vids.resize(F.cols());
        for (size_t j = 0; j < F.cols(); j++)
            c.Vids[j] = F(i, j);
    }
}
void WriteVtk(const ElementType cellType = TRIANGLE)
{
    std::vector<Vertex> V;
    GetVertices(sData.V_o, V);
    std::vector<Cell> C;
    GetCells(sData.F, C);
    static int iter = 0;
    std::string filename = std::string("iter.") + std::to_string(iter++) + ".vtk";
    MeshFileWriter writer(V, C, filename.c_str(), cellType);
    writer.WriteFile();
}
const size_t TetToHex[8] = {4, 5, 6, 7, 2, 3, 0, 1};
void GetHexCells(const Eigen::MatrixXi& F, std::vector<Cell>& C)
{
    C.resize(F.rows()/8);
    for (size_t i = 0; i < F.rows(); i+=8) {
        Cell& c = C.at(i/8);
        c.id = i/8;
        c.Vids.resize(8);
        for (size_t j = 0; j < 8; j++)
            c.Vids[TetToHex[j]] = F(i + j, 0);
    }
}
void WriteHexVtk(const ElementType cellType = HEXAHEDRA)
{
    std::vector<Vertex> V;
    GetVertices(sData.V_o, V);
    std::vector<Cell> C;
    GetHexCells(sData.F, C);
    static int iter = 0;
    std::string filename = std::string("Hex.") + std::to_string(iter++) + ".vtk";
    MeshFileWriter writer(V, C, filename.c_str(), cellType);
    writer.WriteFile();
}

void GetVertices(const std::vector<Vertex>& V, Eigen::MatrixXd& V_o)
{
    V_o.resize(V.size(), 3);
    for (size_t i = 0; i < V_o.rows(); i++) {
        const Vertex& v = V.at(i);
        for (size_t j = 0; j < V_o.cols(); j++)
            V_o(i, j) = v[j];
    }
}

void GetCells(const std::vector<Cell>& C, Eigen::MatrixXi& F)
{
    F.resize(C.size(), C[0].Vids.size());
    for (size_t i = 0; i < F.rows(); i++) {
        const Cell& c = C.at(i);
        for (size_t j = 0; j < F.cols(); j++)
            F(i, j) = c.Vids[j];
    }
}
const size_t HexToTet[8][4] =
{
    {4, 6, 7, 3},
    {5, 7, 4, 0},
    {6, 4, 5, 1},
    {7, 5, 6, 2},
    {2, 0, 3, 7},
    {3, 1, 0, 4},
    {0, 2, 1, 5},
    {1, 3, 2, 6}
};

void ConvertHexToTet(const Eigen::MatrixXi& H, Eigen::MatrixXi& T)
{
    T.resize(H.rows() * 8, 4);
    for (size_t i = 0; i < H.rows(); i++)
        for (size_t j = 0; j < 8; j++)
            for (size_t k = 0; k < 4; k++)
                T(i * 8 + j, k) = H(i, HexToTet[j][k]);
}

size_t GetClosestBoundaryVertexId(const Vertex& v, const std::vector<Vertex>& V)
{
    size_t id = 0;
    double distance = MAXID;
    for (size_t i = 0; i < V.size(); i++) {
        const Vertex& tv = V.at(i);
        if (tv.isBoundary) {
            double d = glm::length(v - tv);
            if (d < distance) {
                distance = d;
                id = tv.id;
            }
        }
    }
    return id;
}
void GetPolycubeConstraints(const Mesh& tetOrigMesh, const Mesh& tetPolycubeMesh, const Mesh& hexPolycubeMesh,
                            Eigen::VectorXi& b, Eigen::MatrixXd& bc)
{
    int numOfVerticesOnBondary = 0;
    for (size_t i = 0; i < hexPolycubeMesh.V.size(); i++) {
        const Vertex& v = hexPolycubeMesh.V.at(i);
        if (v.isBoundary)
            numOfVerticesOnBondary++;
    }
    b.resize(numOfVerticesOnBondary);
    bc.resize(numOfVerticesOnBondary, 3);
    numOfVerticesOnBondary = 0;
    for (size_t i = 0; i < hexPolycubeMesh.V.size(); i++) {
        const Vertex& v = hexPolycubeMesh.V.at(i);
        if (v.isBoundary) {
            b(numOfVerticesOnBondary) = v.id;
            const size_t id = GetClosestBoundaryVertexId(v, tetPolycubeMesh.V);
            for (size_t j = 0; j < 3; j++)
                bc(numOfVerticesOnBondary, j) = tetOrigMesh.V[id][j];
            numOfVerticesOnBondary++;
        }
    }
}
//bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifier){
//  if (key == ' ') {
//    switch (demo_type) {
//      case PARAM_2D: {
//        param_2d_demo_iter(viewer);
//        break;
//      }
//      case SOFT_CONST: {
//        soft_const_demo_iter(viewer);
//        break;
//      }
//      case DEFORM_3D: {
//        deform_3d_demo_iter(viewer);
//        break;
//      }
//      default:
//        break;
//    }
//  }
//
//  return false;
//}

void key_down(int iter)
{
    while (iter--) {
        switch (demo_type)
        {
        case PARAM_2D:
        {
            param_2d_demo_iter();
            break;
        }
        case SOFT_CONST:
        {
            soft_const_demo_iter();
            break;
        }
        case DEFORM_3D:
        {
            deform_3d_demo_iter();
            break;
        }
        default:
            break;
        }
    }
}

void param_2d_demo_iter(/*igl::viewer::Viewer& viewer*/)
{
    if (first_iter) {
        timer.start();
        igl::read_triangle_mesh(TUTORIAL_SHARED_PATH "/face.obj", V, F);
        check_mesh_for_issues(V, F);
        cout << "\tMesh is valid!" << endl;

        Eigen::MatrixXd uv_init;
        Eigen::VectorXi bnd;
        Eigen::MatrixXd bnd_uv;
        igl::boundary_loop(F, bnd);
        igl::map_vertices_to_circle(V, bnd, bnd_uv);

        igl::harmonic(V, F, bnd, bnd_uv, 1, uv_init);
        if (igl::flipped_triangles(uv_init, F).size() != 0) {
            igl::harmonic(F, bnd, bnd_uv, 1, uv_init); // use uniform laplacian
        }

        cout << "initialized parametrization" << endl;

        sData.slim_energy = igl::SLIMData::SYMMETRIC_DIRICHLET;
        Eigen::VectorXi b;
        Eigen::MatrixXd bc;
        slim_precompute(V, F, uv_init, sData, igl::SLIMData::SYMMETRIC_DIRICHLET, b, bc, 0);

        uv_scale_param = 15 * (1. / sqrt(sData.mesh_area));
//    viewer.data.set_mesh(V, F);
//    viewer.core.align_camera_center(V,F);
//    viewer.data.set_uv(sData.V_o*uv_scale_param);
//    viewer.data.compute_normals();
//    viewer.core.show_texture = true;
        WriteVtk();
        first_iter = false;
    }
    else {
        slim_solve(sData, 1); // 1 iter
//    viewer.data.set_uv(sData.V_o*uv_scale_param);
        cout << "time = " << timer.getElapsedTime() << endl;
        WriteVtk();
        WriteHexVtk();
    }
}

void soft_const_demo_iter(/*igl::viewer::Viewer& viewer*/)
{
    if (first_iter) {

        igl::read_triangle_mesh(TUTORIAL_SHARED_PATH "/circle.obj", V, F);

        check_mesh_for_issues(V, F);
        cout << "\tMesh is valid!" << endl;
        Eigen::MatrixXd V_0 = V.block(0, 0, V.rows(), 2);

        Eigen::VectorXi b;
        Eigen::MatrixXd bc;
        get_soft_constraint_for_circle(V_0, F, b, bc);
        double soft_const_p = 1e5;
        slim_precompute(V, F, V_0, sData, igl::SLIMData::SYMMETRIC_DIRICHLET, b, bc, soft_const_p);

//    viewer.data.set_mesh(V, F);
//    viewer.core.align_camera_center(V,F);
//    viewer.data.compute_normals();
//    viewer.core.show_lines = true;

        first_iter = false;

    }
    else {
        timer.start();
        slim_solve(sData, 1); // 1 iter
//    viewer.data.set_mesh(sData.V_o, F);
        static int iter = 1;
        cout << "iter = " << iter++ << " time = " << timer.getElapsedTime() << endl;
    }
}

void deform_3d_demo_iter(/*igl::viewer::Viewer& viewer*/)
{
    if (first_iter) {
        //////////////////////////////////////
        MeshFileReader tetOrigReader("airplane1_input_surface.vtk");
        Mesh& tetOrigMesh = (Mesh&)tetOrigReader.GetMesh();
        tetOrigMesh.BuildAllConnectivities();
        tetOrigMesh.ExtractBoundary();

        MeshFileReader tetPolycubeReader("polycube.tet.n.surface.vtk");
        Mesh& tetPolycubeMesh = (Mesh&)tetPolycubeReader.GetMesh();
        tetPolycubeMesh.BuildAllConnectivities();
        tetPolycubeMesh.ExtractBoundary();

        MeshFileReader hexPolycubeReader("polycube.10000.vtk");
        Mesh& hexPolycubeMesh = (Mesh&)hexPolycubeReader.GetMesh();
        hexPolycubeMesh.BuildAllConnectivities();
        hexPolycubeMesh.ExtractBoundary();

        Eigen::MatrixXd hexPolycubeV;
        Eigen::MatrixXi hexPolycubeH;
        GetVertices(hexPolycubeMesh.V, hexPolycubeV);
        GetCells   (hexPolycubeMesh.C, hexPolycubeH);

        Eigen::MatrixXi hexPolycubeT;
        ConvertHexToTet(hexPolycubeH, hexPolycubeT);

        V = hexPolycubeV;
        F = hexPolycubeT;
        Eigen::MatrixXd V_0 = V;
        Eigen::VectorXi b;
        Eigen::MatrixXd bc;
        GetPolycubeConstraints(tetOrigMesh, tetPolycubeMesh, hexPolycubeMesh, b, bc);

        //////////////////////////////////////
//        igl::readOBJ(TUTORIAL_SHARED_PATH "/cube_40k.obj", V, F);
//        Eigen::MatrixXd V_0 = V;
//        Eigen::VectorXi b;
//        Eigen::MatrixXd bc;
//        get_cube_corner_constraints(V_0, F, b, bc);

        double soft_const_p = 1e3;
        sData.exp_factor = 5.0;
        //slim_precompute(V, F, V_0, sData, igl::SLIMData::EXP_CONFORMAL, b, bc, soft_const_p);
        slim_precompute(V, F, V_0, sData, igl::SLIMData::EXP_CONFORMAL, b, bc, soft_const_p);
        cout << "precomputed" << endl;

        first_iter = false;
        WriteVtk(TETRAHEDRA);
        //display_3d_mesh (viewer);
        //display_3d_mesh();
    }
    else {
        timer.start();
        slim_solve(sData, 1); // 1 iter
        static int iter = 1;
        cout << "iter = " << iter++ << " time = " << timer.getElapsedTime() << endl;
        //display_3d_mesh (viewer);
        WriteVtk(TETRAHEDRA);
        WriteHexVtk(HEXAHEDRA);
        //display_3d_mesh();
    }
}

void display_3d_mesh(/*igl::viewer::Viewer& viewer*/)
{
    MatrixXd V_temp;
    MatrixXi F_temp;
    Eigen::MatrixXd Barycenters;

    igl::barycenter(sData.V, sData.F, Barycenters);
    //cout << "Barycenters.rows() = " << Barycenters.rows() << endl;
    //double t = double((key - '1')+1) / 9.0;
    double view_depth = 10.;
    double t = view_depth / 9.;

    VectorXd v = Barycenters.col(2).array() - Barycenters.col(2).minCoeff();
    v /= v.col(0).maxCoeff();

    vector<int> s;

    for (unsigned i = 0; i < v.size(); ++i)
        if (v(i) < t) s.push_back(i);

    V_temp.resize(s.size() * 4, 3);
    F_temp.resize(s.size() * 4, 3);

    for (unsigned i = 0; i < s.size(); ++i) {
        V_temp.row(i * 4 + 0) = sData.V_o.row(sData.F(s[i], 0));
        V_temp.row(i * 4 + 1) = sData.V_o.row(sData.F(s[i], 1));
        V_temp.row(i * 4 + 2) = sData.V_o.row(sData.F(s[i], 2));
        V_temp.row(i * 4 + 3) = sData.V_o.row(sData.F(s[i], 3));
        F_temp.row(i * 4 + 0) << (i * 4) + 0, (i * 4) + 1, (i * 4) + 3;
        F_temp.row(i * 4 + 1) << (i * 4) + 0, (i * 4) + 2, (i * 4) + 1;
        F_temp.row(i * 4 + 2) << (i * 4) + 3, (i * 4) + 2, (i * 4) + 0;
        F_temp.row(i * 4 + 3) << (i * 4) + 1, (i * 4) + 2, (i * 4) + 3;
    }
//  viewer.data.set_mesh(V_temp,F_temp);
//  viewer.core.align_camera_center(V_temp,F_temp);
//  viewer.data.set_face_based(true);
//  viewer.core.show_lines = true;
}

int main(int argc, char *argv[])
{
    cerr << "Press space for running an iteration." << std::endl;
    cerr << "Syntax: " << argv[0] << " demo_number (1 to 3)" << std::endl;
    cerr << "1. 2D unconstrained parametrization" << std::endl;
    cerr << "2. 2D deformation with positional constraints" << std::endl;
    cerr << "3. 3D mesh deformation with positional constraints" << std::endl;

    demo_type = PARAM_2D;

    if (argc == 2) {
        demo_type = DEMO_TYPE(std::atoi(argv[1]) - 1);
        if (demo_type > DEFORM_3D || demo_type < PARAM_2D)
            cerr << "Wrong demo number - Please choose one between 1-3" << std::endl;
    }
    key_down(50);
    return 0;
}

void check_mesh_for_issues(Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{

    Eigen::SparseMatrix<double> A;
    igl::adjacency_matrix(F, A);

    Eigen::MatrixXi C, Ci;
    igl::components(A, C, Ci);

    int connected_components = Ci.rows();
    if (connected_components != 1) {
        cout << "Error! Input has multiple connected components" << endl;
        exit(1);
    }
    int euler_char = igl::euler_characteristic(V, F);
    if (!euler_char) {
        cout << "Error! Input does not have a disk topology, it's euler char is " << euler_char << endl;
        exit(1);
    }
    bool is_edge_manifold = igl::is_edge_manifold(F);
    if (!is_edge_manifold) {
        cout << "Error! Input is not an edge manifold" << endl;
        exit(1);
    }

    Eigen::VectorXd areas;
    igl::doublearea(V, F, areas);
    const double eps = 1e-14;
    for (int i = 0; i < areas.rows(); i++) {
        if (areas(i) < eps) {
            cout << "Error! Input has zero area faces" << endl;
            exit(1);
        }
    }
}

void get_soft_constraint_for_circle(Eigen::MatrixXd& V_o, Eigen::MatrixXi& F, Eigen::VectorXi& b, Eigen::MatrixXd& bc)
{

    Eigen::VectorXi bnd;
    igl::boundary_loop(F, bnd);
    const int B_STEPS = 22; // constraint every B_STEPS vertices of the boundary

    b.resize(bnd.rows() / B_STEPS);
    bc.resize(b.rows(), 2);

    int c_idx = 0;
    for (int i = B_STEPS; i < bnd.rows(); i += B_STEPS) {
        b(c_idx) = bnd(i);
        c_idx++;
    }

    bc.row(0) = V_o.row(b(0)); // keep it there for now
    bc.row(1) = V_o.row(b(2));
    bc.row(2) = V_o.row(b(3));
    bc.row(3) = V_o.row(b(4));
    bc.row(4) = V_o.row(b(5));

    bc.row(0) << V_o(b(0), 0), 0;
    bc.row(4) << V_o(b(4), 0), 0;
    bc.row(2) << V_o(b(2), 0), 0.1;
    bc.row(3) << V_o(b(3), 0), 0.05;
    bc.row(1) << V_o(b(1), 0), -0.15;
    bc.row(5) << V_o(b(5), 0), +0.15;
}

void get_cube_corner_constraints(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::VectorXi& b, Eigen::MatrixXd& bc)
{
    double min_x, max_x, min_y, max_y, min_z, max_z;
    min_x = V.col(0).minCoeff();
    max_x = V.col(0).maxCoeff();
    min_y = V.col(1).minCoeff();
    max_y = V.col(1).maxCoeff();
    min_z = V.col(2).minCoeff();
    max_z = V.col(2).maxCoeff();

    // get all cube corners
    b.resize(8, 1);
    bc.resize(8, 3);
    int x;
    for (int i = 0; i < V.rows(); i++) {
        if (V.row(i) == Eigen::RowVector3d(min_x, min_y, min_z)) b(0) = i;
        if (V.row(i) == Eigen::RowVector3d(min_x, min_y, max_z)) b(1) = i;
        if (V.row(i) == Eigen::RowVector3d(min_x, max_y, min_z)) b(2) = i;
        if (V.row(i) == Eigen::RowVector3d(min_x, max_y, max_z)) b(3) = i;
        if (V.row(i) == Eigen::RowVector3d(max_x, min_y, min_z)) b(4) = i;
        if (V.row(i) == Eigen::RowVector3d(max_x, max_y, min_z)) b(5) = i;
        if (V.row(i) == Eigen::RowVector3d(max_x, min_y, max_z)) b(6) = i;
        if (V.row(i) == Eigen::RowVector3d(max_x, max_y, max_z)) b(7) = i;
    }

    // get all cube edges
    std::set<int> cube_edge1;
    Eigen::VectorXi cube_edge1_vec;
    for (int i = 0; i < V.rows(); i++) {
        if ((V(i, 0) == min_x && V(i, 1) == min_y)) {
            cube_edge1.insert(i);
        }
    }
    Eigen::VectorXi edge1;
    int_set_to_eigen_vector(cube_edge1, edge1);

    std::set<int> cube_edge2;
    Eigen::VectorXi edge2;
    for (int i = 0; i < V.rows(); i++) {
        if ((V(i, 0) == max_x && V(i, 1) == max_y)) {
            cube_edge2.insert(i);
        }
    }
    int_set_to_eigen_vector(cube_edge2, edge2);
    b = igl::cat(1, edge1, edge2);

    std::set<int> cube_edge3;
    Eigen::VectorXi edge3;
    for (int i = 0; i < V.rows(); i++) {
        if ((V(i, 0) == max_x && V(i, 1) == min_y)) {
            cube_edge3.insert(i);
        }
    }
    int_set_to_eigen_vector(cube_edge3, edge3);
    b = igl::cat(1, b, edge3);

    std::set<int> cube_edge4;
    Eigen::VectorXi edge4;
    for (int i = 0; i < V.rows(); i++) {
        if ((V(i, 0) == min_x && V(i, 1) == max_y)) {
            cube_edge4.insert(i);
        }
    }
    int_set_to_eigen_vector(cube_edge4, edge4);
    b = igl::cat(1, b, edge4);

    bc.resize(b.rows(), 3);
    Eigen::Matrix3d m;
    m = Eigen::AngleAxisd(0.3 * M_PI, Eigen::Vector3d(1. / sqrt(2.), 1. / sqrt(2.), 0.)/*Eigen::Vector3d::UnitX()*/);
    int i = 0;
    for (; i < cube_edge1.size(); i++) {
        Eigen::RowVector3d edge_rot_center(min_x, min_y, (min_z + max_z) / 2.);
        bc.row(i) = (V.row(b(i)) - edge_rot_center) * m + edge_rot_center;
    }
    for (; i < cube_edge1.size() + cube_edge2.size(); i++) {
        Eigen::RowVector3d edge_rot_center(max_x, max_y, (min_z + max_z) / 2.);
        bc.row(i) = (V.row(b(i)) - edge_rot_center) * m.transpose() + edge_rot_center;
    }
    for (; i < cube_edge1.size() + cube_edge2.size() + cube_edge3.size(); i++) {
        bc.row(i) = 0.75 * V.row(b(i));
    }
    for (; i < b.rows(); i++) {
        bc.row(i) = 0.75 * V.row(b(i));
    }
}

void int_set_to_eigen_vector(const std::set<int>& int_set, Eigen::VectorXi& vec)
{
    vec.resize(int_set.size());
    int idx = 0;
    for (auto f : int_set) {
        vec(idx) = f;
        idx++;
    }
}

