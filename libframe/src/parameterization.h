#pragma once
//1. dtermine the parametric numbers, three directions for each components.
//2. line, face, and component parameterization,
#include "global_types.h"
#include "global_functions.h"
#include "io.h"
#include <math.h>
#include <time.h>
#include "tetgen.h"
//#include "Eigen/Dense"
#include "Eigen/Sparse"
using namespace Eigen;
typedef SparseMatrix<float> SpMat; // declares a column-major sparse matrix type of double
typedef Triplet<float> TT;
class parameterization
{
	vector<float> Sum_cuboid_hexes;
	vector<int> L_para_Ns, es_Flags;
	vector<int> Us, Vs, Ws;

	bool enlarged;
	float Escalar;
public:
	vector<Hex_V> hvs_parameterized;
	vector<Hex> hhs_parameterized;
public:

	parameterization(void);
	void parameterization_main(char *path);

	void determine_parametric_numbers();
	vector<int> an_edge_area(int e_id, vector<bool> &arrayE_test);

	void parameterization_cuboid_tetgen(int C_Id, Parameterization_Cell &pc);
	void cuboid2triangle_mesh(int C_Id, vector<int> &Vs_arrayIds, vector<Vertex> &vs, vector<Triangle> &fs);
	void translate_to_tetmesh_tetgen(vector<int> &Vs_arrayIds, vector<Vertex> &vs, vector<Triangle> &fs, vector<int> &VTs_arrayIds,
			vector<Vertex> &vts, vector<Triangle> &fts, vector<Tet> &tets);
	bool D2_parameterization_tetgen(int F_Id, Parameterization_Cell &pc, vector<int> &VTs_arrayIds, vector<Vertex> &vts, vector<Triangle> &fts);
	bool D3_parameterization_tetgen(int C_Id, Parameterization_Cell &pc, vector<Vertex> &vts, vector<Edge> &ets, vector<Triangle> &fts,
			vector<Tet> &tets);

//a hex to 24 tets
	void parameterization_cuboid(int C_Id, Parameterization_Cell &pc);
	void reindexing_cuboid_mesh(int C_Id, vector<int> &Vs_arrayIds, vector<Hex_V> &vs, vector<Hex_E> &es, vector<Hex_F> &fs, vector<Hex> &hs);
	void translate_to_tetmesh(vector<int> &Vs_arrayIds, vector<Hex_V> &vs, vector<Hex_E> &es, vector<Hex_F> &fs, vector<Hex> &hs,
			vector<int> &VTs_arrayIds, vector<Vertex> &vts, vector<Edge> &ets, vector<Triangle> &fts, vector<Tet> &tets);
	bool D2_parameterization(int F_Id, Parameterization_Cell &pc, vector<int> &VTs_arrayIds, vector<Vertex> &vts, vector<Edge> &ets,
			vector<Triangle> &fts);
	bool D3_parameterization(int C_Id, Parameterization_Cell &pc, vector<Vertex> &vts, vector<Edge> &ets, vector<Triangle> &fts, vector<Tet> &tets);
//a hex to 24 tets

	void parametric_coords(int C_Id, Parameterization_Cell &pc, vector<Vertex> &vts, vector<Tet> &tets);

	void produce_hex_mesh(vector<Parameterization_Cell> &PCs, char *path);

};

