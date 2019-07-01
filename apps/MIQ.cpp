#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "EdgeLines.h"
#include "BaseComplex.h"
#include "BaseComplexQuad.h"
#include "BaseComplexSheet.h"
#include "BaseComplexSheetQuad.h"
#include "BaseComplexChord.h"
#include "BaseComplexEditor.h"
#include "SingularityGraph.h"
#include "ArgumentManager.h"

#include <igl/avg_edge_length.h>
#include <igl/barycenter.h>
#include <igl/comb_cross_field.h>
#include <igl/comb_frame_field.h>
#include <igl/compute_frame_field_bisectors.h>
#include <igl/cross_field_missmatch.h>
#include <igl/cut_mesh_from_singularities.h>
#include <igl/find_cross_field_singularities.h>
#include <igl/local_basis.h>
#include <igl/readOFF.h>
#include <igl/rotate_vectors.h>
#include <igl/copyleft/comiso/miq.h>
#include <igl/copyleft/comiso/nrosy.h>
#include <igl/writeOBJ.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <igl/serialize.h>

void writeOBJ(const char* filename, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::MatrixXd& UV, const Eigen::MatrixXi& FUV) {
    using namespace std;
    using namespace Eigen;
    assert(V.cols() == 3 && "V should have 3 columns");
    ofstream s(filename);
    s << V.format(IOFormat(FullPrecision, DontAlignCols, " ", "\n", "v ", "", "", "\n"));
    for (int i = 0; i < FUV.rows(); ++i)
        for (int j = 0; j < 3; ++j)
            s << "vt " << UV(FUV(i, j), 0) << " " << UV(FUV(i, j), 1) << "\n";
    for (int i = 0; i < FUV.rows(); ++i) {
        s << "f ";
        for (int j = 0; j < 3; ++j)
            s << F(i, j) + 1 << "/" << 3 * i + j + 1 << " ";
        s << "\n";
    }
}

// Input mesh
Eigen::MatrixXd V;
Eigen::MatrixXi F;

// Face barycenters
Eigen::MatrixXd B;

// Scale for visualizing the fields
double global_scale;
bool extend_arrows = false;

// Cross field
Eigen::MatrixXd X1, X2;

// Bisector field
Eigen::MatrixXd BIS1, BIS2;

// Combed bisector
Eigen::MatrixXd BIS1_combed, BIS2_combed;

// Per-corner, integer mismatches
Eigen::Matrix<int, Eigen::Dynamic, 3> MMatch;

// Field singularities
Eigen::Matrix<int, Eigen::Dynamic, 1> isSingularity, singularityIndex;

// Per corner seams
Eigen::Matrix<int, Eigen::Dynamic, 3> Seams;

// Combed field
Eigen::MatrixXd X1_combed, X2_combed;

// Global parametrization (with seams)
Eigen::MatrixXd UV_seams;
Eigen::MatrixXi FUV_seams;

// Global parametrization
Eigen::MatrixXd UV;
Eigen::MatrixXi FUV;

// Create a texture that hides the integer translation in the parametrization
void line_texture(Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> &texture_R, Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> &texture_G,
        Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> &texture_B) {
    unsigned size = 128;
    unsigned size2 = size / 2;
    unsigned lineWidth = 3;
    texture_R.setConstant(size, size, 255);
    for (unsigned i = 0; i < size; ++i)
        for (unsigned j = size2 - lineWidth; j <= size2 + lineWidth; ++j)
            texture_R(i, j) = 0;
    for (unsigned i = size2 - lineWidth; i <= size2 + lineWidth; ++i)
        for (unsigned j = 0; j < size; ++j)
            texture_R(i, j) = 0;

    texture_G = texture_R;
    texture_B = texture_R;
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        std::cout << "Usage: MIQ inputTri.off outputUV.obj \n";
        return 0;
    }
    // #################
    MeshFileReader reader(argv[1]);
    Mesh& mesh = (Mesh&)reader.GetMesh();
    mesh.BuildAllConnectivities();
    mesh.ExtractBoundary();
    mesh.ExtractSingularities();
    //"Info: cos10 = 0.984807753; cos15 = 0.965925826; cos20 = 0.939692621; cos25 = 0.906307787; cos30 = 0.866025404\n\n";
    const double cosangle = 0.866025404;
    mesh.SetCosAngleThreshold(cosangle);
    mesh.LabelSurface();
    mesh.LabelSharpEdges(true);
    std::vector<int> roundVertexIds;
    std::vector<std::vector<int>> sharpEdgeVertexIds;
    for (auto& v : mesh.V)
        if (v.isCorner) roundVertexIds.push_back(v.id);
    for (auto& e : mesh.E)
        if (e.isBoundary or e.isSharpFeature) sharpEdgeVertexIds.push_back({e.Vids[0], e.Vids[1]});
    {
        std::ofstream ofs("corners.txt");
        for (auto c : roundVertexIds)
            ofs << c << "\n";
    }
    {
        std::ofstream ofs("sharpEdges.txt");
        for (auto& e : sharpEdgeVertexIds)
            ofs << e.front() << " " << e.back() << "\n";
    }
    // #################
    using namespace Eigen;
    igl::readOFF(argv[1], V, F);

    double gradient_size = 200;
    double iter = 0;
    double stiffness = 5.0;
    bool direct_round = 0;

    // Compute face barycenters
    igl::barycenter(V, F, B);

    // Compute scale for visualizing fields
    global_scale = .5 * igl::avg_edge_length(V, F);

    // Contrain one face
    VectorXi b(1);
    b << 0;
    MatrixXd bc(1, 3);
    bc << 1, 0, 0;

    // Create a smooth 4-RoSy field
    VectorXd S;
    igl::copyleft::comiso::nrosy(V, F, b, bc, VectorXi(), VectorXd(), MatrixXd(), 4, 0.5, X1, S);

    // Find the orthogonal vector
    MatrixXd B1, B2, B3;
    igl::local_basis(V, F, B1, B2, B3);
    X2 = igl::rotate_vectors(X1, VectorXd::Constant(1, M_PI / 2), B1, B2);

    // Always work on the bisectors, it is more general
    igl::compute_frame_field_bisectors(V, F, X1, X2, BIS1, BIS2);

    // Comb the field, implicitly defining the seams
    igl::comb_cross_field(V, F, BIS1, BIS2, BIS1_combed, BIS2_combed);

    // Find the integer mismatches
    igl::cross_field_missmatch(V, F, BIS1_combed, BIS2_combed, true, MMatch);

    // Find the singularities
    igl::find_cross_field_singularities(V, F, MMatch, isSingularity, singularityIndex);

    // Cut the mesh, duplicating all vertices on the seams
    igl::cut_mesh_from_singularities(V, F, MMatch, Seams);

    // Comb the frame-field accordingly
    igl::comb_frame_field(V, F, X1, X2, BIS1_combed, BIS2_combed, X1_combed, X2_combed);

//    std::vector<int> roundVertexIds;
//    std::vector<std::vector<int>> sharpEdgeVertexIds;
//    for (auto& v : mesh.V)
//        if (v.isCorner) roundVertexIds.push_back(v.id);
//    for (auto& e : mesh.E)
//        if (e.isBoundary or e.isSharpFeature) sharpEdgeVertexIds.push_back({e.Vids[0], e.Vids[1]});
    // Global parametrization
    igl::copyleft::comiso::miq(V, F, X1_combed, X2_combed, MMatch, isSingularity, Seams, UV, FUV, gradient_size, stiffness, direct_round, iter, 5, true);

// Global parametrization (with seams, only for demonstration)
    igl::copyleft::comiso::miq(V, F, X1_combed, X2_combed, MMatch, isSingularity, Seams, UV_seams, FUV_seams, gradient_size, stiffness, direct_round, iter, 5, false);

    if (argc >= 3) {
        writeOBJ(argv[2], V, F, UV, FUV);
        std::cout << "V info\n";
        std::cout << "rows = " << V.rows() << " cols = " << V.cols() << "\n";
        std::cout << "F info\n";
        std::cout << "rows = " << F.rows() << " cols = " << F.cols() << "\n";
        std::cout << "UV info\n";
        std::cout << "rows = " << UV.rows() << " cols = " << UV.cols() << "\n";
        std::cout << "FUV info\n";
        std::cout << "rows = " << FUV.rows() << " cols = " << FUV.cols() << "\n";
        std::cout << "F(0) = " << F(0, 0) << " " << F(0, 1) << " " << F(0, 2) << "\n";
        std::cout << "FUV(0) = " << FUV(0, 0) << " " << FUV(0, 1) << " " << FUV(0, 2) << "\n";
    }
    if (argc <= 3) return 0;
}
