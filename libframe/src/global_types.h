#ifndef GLOBAL_TYPE_HEADER
#define GLOBAL_TYPE_HEADER
#include "constants.h"
#include <vector>
using namespace std;

struct Vertex
{
	int index;
	float v[3];
	float normal[3];
	vector<int> neighborv;
	vector<int> neighbort;
	vector<int> neighbor_ot; //opposite t for tet
	vector<int> neighbortet;

	int is_boundary;

	bool on_base_complex;
	int which_F_face;
};
struct Triangle
{
	int index;
	int triangle_v[3];
	float normal[3];
	vector<int> neighbort;
	float area;

	int which_F_face;
};
struct Edge
{
	int index;
	int startend_id[2];
	bool isboundary;
};

struct MapPatch
{
	vector<int> tv_ids;
	vector<int> tt_ids;
};

struct Tet
{
	int index;
	vector<int> vs;
	vector<int> es;
	vector<int> ts;
	bool isboundary;
	int Color_Id;
};
//-------------------------------------------------------------------
//---For Hybrid mesh-------------------------------------------------
struct Hybrid_V
{
	int index;
	vector<float> v;
	vector<int> neighbor_Vs;
	vector<int> neighbor_Es;
	vector<int> neighbor_Fs;
	vector<int> neighbor_Hs;
	bool isboundary;
	int Color_Id;
};
struct Hybrid_E
{
	int index;
	vector<int> startend_Id;
	vector<int> neighbor_Fs;
	vector<int> neighbor_Hs;
	bool isboundary;
	int Color_Id;
};
struct Hybrid_F
{
	int index;
	vector<int> vs;
	vector<int> es;
	vector<int> neighbor_Hs;
	bool isboundary;
	int Color_Id;
};

struct Hybrid
{
	int index;
	vector<int> vs;
	vector<int> es;
	vector<int> fs;
	bool isboundary;
	int Color_Id;
};

//-------------------------------------------------------------------
//-------------------------------------------------------------------

struct Hex_V
{
	int index;
	float v[3];
	float normal[3];

	int where_location; //is inside the surface or not?
	int slice_id;
	int Frame_V_id; //-1 means no. index start from 0
	int is_Extra_ordinary;
	vector<int> neighbor_SEs;

	int on_whichFrame_edge; //-1 means no, index start from 0
	int which_F_face;
	int is_on_base_complex;
	int is_on_base_complex_edge;

	vector<int> neighbor_vs;
	vector<int> neighbor_Es;
	vector<int> neighbor_Fs;
	vector<int> neighbor_Hs;

	bool fixed; //for mesquite

};
struct Hex_E
{
	int index;
	int startend_Id[2];
	float length;
	int is_boundary;
	int checked;

	vector<int> neighbor_Es; //1-opposite ring
	vector<int> neighbor_Fs;
	vector<int> neighbor_Hs;

	int which_F_edge; //which base-complex edge
};
struct Hex_F
{
	int index;
	int is_boundary;
	int cv_Ids[4];
	int ce_Ids[4];

	float normal[3];

	float centroid[3];

	vector<int> neighbor_Es;

	vector<int> neighbor_Fs;

	vector<int> neighbor_Cs;

	int frame_boundary; //for hex-mesh component tagging
	int which_F_face; //which base-complex face
};
struct Hex
{
	int index;
	int at_boundary;
	int F_Ids[6]; //

	int V_Ids[8];
	int E_Ids[12];

	float centroid[3];
	float volume;
	float scalar_jacobian;
	float alpha;

	vector<int> neighbor_FS; //it's 6 faces

	vector<int> neighbor_CS; //neighboring cube

	vector<int> neighbor_ES; //it's 4 edges in sequential

	int frame_component; //for hex-esh component tag

	int Color_ID;
	int Color_E_ID;
	int Color_EF_ID;

};

struct Volume_Piece_Info
{
	int DegreeU;
	int DegreeV;
	int DegreeW;
	int m, n, o;
	vector<Hex_V> Vs;
	vector<float> T_ms;
	vector<float> T_ns;
	vector<float> T_os;
};

struct Frame_V
{
	int index_own;
	int index_hex;

	int what_type; //1 singular, 2, extra-nordinary, 3, regular
	//0 non extra-ordinary node and at surface,1 extra-ordinary node and at surface;
	//2 non extra-ordinary node and inside volume,3 extra-ordinary node and inside volume;
	int velence;

	vector<int> neighbor_vs;
	vector<int> neighbor_Es;
	vector<int> neighbor_Fs;
	vector<int> neighbor_Hs;

	int which_singularity;
	int which_singularity_type;

};
struct Frame_E
{
	int index_own;
	int startend_Id[2];
	int is_boundary;	              //0 interior, 1 boundary, 2 first and last slice boundary
	vector<int> vs_link;	              //v of hex_v

	vector<int> neighbor_Es;
	vector<int> opposite_Es;
	vector<int> neighbor_Fs;
	vector<int> neighbor_Hs;
	int which_type;	              //0,1,2,....
	int singularity_type;	              //0,1,2,....

	int singular_Id;

	int parameterization_U;	              //parameterization
};
struct Frame_F
{
	int index_own;
	int is_boundary;
	int F_location;	              //for rendering

	int fv_Ids[4];
	int fe_Ids[4];
	int U, V;
	vector<vector<int> > hfs_net;
	vector<int> hfs_net_another;

	vector<int> neighbor_Es;
	vector<int> neighbor_Fs;
	vector<int> neighbor_Cs;

	int p_U;	              //parameterization
	int p_V;	              //parameterization

	int Color_ID;
};
struct Frame_H
{
	int index_own;

	int FF_Ids[6];

	int FV_Ids[8];
	int FE_Ids[12];

	int U, V, W;	              //resolutions
	vector<vector<vector<int> > > vs_net;
	vector<int> fs_net;
	vector<int> hs_net;

	vector<int> neighbor_FS;	              //it's 6 faces
	vector<int> neighbor_HS;	              //neighboring cube
	vector<int> neighbor_ES;	              //it's 4 edges in sequential

	int p_U;	              //parameterization
	int p_V;	              //parameterization
	int p_W;	              //parameterization

	int Color_ID;
};

struct Singular_V
{
	int index_own;
	int index_hex;

	int is_boundary;

	vector<int> neighbor_Vs;	              //singular vs
	vector<int> neighbor_Es;	              //singular es

	bool fake;

	int which_singularity;
	int which_singularity_type;
};
struct Singular_E
{
	int index_own;
	int startend_Id[2];	              //singular v
	int is_boundary;	              //-1 interior, 1 boundary,
	vector<int> es_link;	              //hex e
	vector<int> vs_link;	              //hex v

	vector<int> neighbor_Es;	              //singular es
	vector<int> composed_BS_Es;	              //base-complex es, one or more base-complex edges compose one singular edge

	vector<vector<vector<int> > > separatrixSheets;	              //composed of base-complex fs

	int singularity_type;	              //0,1, 2,.... for rendering

	int edge_type;	              //0 for normal, 2 for circle,

	int color_id;
};

struct diag_op
{
	int id;
	vector<int> fs;
	vector<int> hs;
	vector<int> vs1, vs2;

	int is_circular;
};

struct Parameterization_Cell
{
	vector<vector<float> > UVW_coords;
	vector<vector<vector<Hex_V> > > Parametric_coords;
};
struct para_coor
{
	vector<float> coords;
};
struct fan_chain
{
	vector<vector<vector<para_coor> > > UVW_coords;
	vector<vector<vector<int> > > hv_ids;
};
struct Parameterization_Chain
{
	vector<fan_chain> fan_cicle;
};

struct Tri_Mesh_Info
{
	vector<Vertex> Vs;
	vector<Triangle> Ts;
};

struct Quad_Mesh_Info
{

};

struct Tet_Mesh_Info
{
	vector<Vertex> Vs;
	vector<Edge> Es;
	vector<Triangle> Ts;
	vector<Tet> Tets;
};
struct Singularity_Info
{
	vector<Singular_V> SVs;
	vector<Singular_E> SEs;
};

struct Frame_Info
{
	vector<Frame_V> FVs;
	vector<Frame_E> FEs;
	vector<Frame_F> FFs;
	vector<Frame_H> FHs;

	vector<Hex_E> HEs;
	vector<Hex_V> HVs;
};

struct Hybrid_Mesh_Info
{
	vector<Hybrid_V> HVs;
	vector<Hybrid_E> HEs;
	vector<Hybrid_F> HFs;
	vector<Hybrid> HHs;
};
struct Hex_Mesh_Info
{
	vector<Hex_V> HVs;
	vector<Hex_E> HEs;
	vector<Hex_F> HFs;
	vector<Hex> HHs;
	vector<Volume_Piece_Info> VPs;
	vector<Volume_Piece_Info> VPs_spline_fitted;
	vector<Frame_Info> Frs;
	float average_e_len;
};

struct Char_Char
{
	char path[512];
};


struct Component_Chord
{
    vector<int> all_nodes;
    vector<int> nodes[4];

    vector<int> all_edges;
    vector<int> middle_edges[4],edges[4];

    vector<int> all_faces;
    vector<int> middle_faces,faces[4];

    vector<int> components;

    int id;
    int side;
    int which_type;//1 open, 2 close/circle,3 extraordinary (sharing edges, tangling)
    bool self_intersect;//is there one or several components collapse to an edge?
    int which_group;//component sheets in each group have no sharing components
    vector<float> weights;
    float weight;

};
struct Component_Sheet
{
    vector<int> all_nodes;
    vector<int> nodes_left,nodes_right;

    vector<int> all_edges;
    vector<int> middle_edges,edges_left,edges_right;

    vector<int> all_faces;
    vector<int> middle_faces,faces_left,faces_right;

    vector<int> components;

    int id;//which one
    int which_type;//1 open, 2 close
    bool self_intersect;//is there one or several components collapse to an edge?
    int which_group;//component sheets in each group have no sharing components

};
struct Chord_Group
{
    vector<int> component_chords;
    int id;//which one
};

struct Structure_Coverage_chord
{
    vector<int> component_chords;
    vector<int> contained_groups;
    vector<int> shared_components;

    int id;//which one
};
struct Structure_Coverage_Group_chord
{
    vector<int> groups;
    vector<int> shared_components;
    int id;//which one
};

struct Structure_Coverage_Sheet
{
    vector<int> component_sheets;
    vector<int> contained_groups;
    vector<int> shared_components;

    int id;//which one
};
struct Structure_Coverage_Group
{
    vector<int> groups;
    vector<int> shared_components;

    int id;//which one
};


//meshes
extern Tri_Mesh_Info tri_mesh;
extern Hybrid_Mesh_Info hyrbid_mesh;
extern Hex_Mesh_Info hex_mesh;
extern Frame_Info FrameI;
extern Singularity_Info SingularityI;
extern int Para_Total_N;
extern int Para_min_N;

//components
extern vector<Component_Sheet> all_com_sheets;
extern vector<Component_Chord> all_com_chords;
extern int Num_Group;
extern vector<vector<int> > All_groups;
extern vector<Structure_Coverage_Group> All_decompositions;
extern int Num_Group_chord;
extern vector<Chord_Group> All_groups_chord;
extern vector<Structure_Coverage_Group_chord> All_decompositions_chord;
#endif
