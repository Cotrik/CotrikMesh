#ifndef __LAPLACIAN_OPTIMIZATION_H__
#define __LAPLACIAN_OPTIMIZATION_H__
//1. base-complex edge chain.
//2. laplacian of singular node
#include "constants.h"
#include "global_types.h"
#include "global_functions.h"
#include "io.h"
#include "parameterization.h"
#include "frame_of_mesh.h"
#include <math.h>
#include <float.h>
#include <time.h>
#include "laplacian_smoothing.h"
#include "Eigen/Sparse"
using namespace Eigen;
typedef SparseMatrix<float> SpMat; // declares a column-major sparse matrix type of double
typedef Triplet<float> TT;
class parametric_optimization
{
	struct edge_chain
	{
		int id;
		int startendIds[2];
		vector<int> bvs;
		vector<int> bes;
		int which_type;
		vector<vector<int>> bvs_1ring;
		vector<vector<int>> bvs_valence_ring;
		vector<vector<int>> bes_1ring;
		vector<int> boundaryes_excepttb;

		vector<vector<para_coor>> ring1_paracoords;

		vector<vector<int>> bfs; //centered valence fs
		vector<vector<int>> bhs;

		vector<int> boundary_fs;
		vector<int> all_fs;
	};
	struct face_center
	{
		int id;
		int f_id;
		vector<int> vs_all; //all nodes
		vector<vector<int>> vs_paralell; //4 in one level

		vector<int> es_all; //all es
		vector<vector<int>> es_paralell; //three or two faces
		vector<int> fs_paralell; //three or two faces
		vector<int> fs_boundary; //except current f

		vector<vector<para_coor>> layer_paracoords;
	};
	//edge
	vector<edge_chain> ECs;
	vector<Parameterization_Chain> ParaCs;
	//face
	vector<face_center> FCs;
	vector<fan_chain> ParaFCs;

	int File_Index;
public:
	parametric_optimization(void);
	void optimization_pipeline(char *fname);
	void optimization_pipelineF(char *fname);

	void optimization_pipeline_edge(int Iter);
	void optimization_pipeline_face(int Iter, bool smooth);

	void optimization_pipeline_edgeface(int Iter, char *fname);
	void find_efcs_cur(vector<int> &ecs_ids, vector<int> &fcs_ids, vector<bool> &ecs_Ind, vector<bool> &fcs_Ind);

//edge
	void edge_chain_extraction();
	void parameter_index_correspondence();
	void parameter_index_correspondence_new();

	void find_ecs_cur(vector<int> &ecs_ids, vector<bool> &ecs_Ind);
	void a_chain_parameterization(int EC_ID);
	void reindexing_cuboid_mesh(int EC_ID, vector<int> &Vs_arrayIds, vector<Hex_V> &vs, vector<Hex_E> &es, vector<Hex_F> &fs, vector<Hex> &hs);
	void parametric_coords(int EC_ID, vector<int> &Vs_arrayIds, Parameterization_Cell &pc, vector<Vertex> &vts, vector<Triangle> &fts,
			vector<Tet> &tets);

//face
	void face_center_extraction();
	void face_parameter_index_correspondence();

	void find_fcs_cur(vector<int> &fcs_ids, vector<bool> &fcs_Ind);
	void face_center_parameterization(int FC_ID);
	void face_reindexing_cuboid_mesh(int FC_ID, vector<int> &Vs_arrayIds, vector<Hex_V> &vs, vector<Hex_E> &es, vector<Hex_F> &fs, vector<Hex> &hs);

	void translate_to_tetmesh(vector<int> &Vs_arrayIds, vector<Hex_V> &vs, vector<Hex_E> &es, vector<Hex_F> &fs, vector<Hex> &hs, vector<Vertex> &vts,
			vector<Edge> &ets, vector<Triangle> &fts, vector<Tet> &tets);
	bool D2_parameterization(int F_Id, Parameterization_Cell &pc, vector<Vertex> &vts, vector<Edge> &ets, vector<Triangle> &fts);
	bool D3_parameterization(Parameterization_Cell &pc, vector<Vertex> &vts, vector<Edge> &ets, vector<Triangle> &fts, vector<Tet> &tets);
	void face_parametric_coords(int EC_ID, Parameterization_Cell &pc, vector<Vertex> &vts, vector<Tet> &tets);
};

#endif // __LAPLACIAN_OPTIMIZATION_H__

