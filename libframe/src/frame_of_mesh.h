#pragma once
#include "global_types.h"
#include "global_functions.h"
#include "io.h"
class frame_of_mesh
{
public:
	frame_of_mesh(void);
	void base_complex_extraction();
	void extract_singular_node_edge();
	void singularity_structure_connectivity();
	void dealwith_circlesingularedge(vector<Singular_E> &circle_ses);
	void extract_base_complex_node_edge();
	void extract_base_complex_face();

	void assign_singular_edge_composedBSEs();
	// assign any hex-mesh components.
	void assign_hex_mesh_component();
	void assign_component();
	void extract_useless_hexes(vector<bool> &hexes_Inds);

	void assign_color_hexfaces_component();

	bool straight_line_test(int h_v1, int h_v2, int h_v3);

	bool is_v_singular(int hvid, int cur_e, int &next_e);
	~frame_of_mesh(void);
};

