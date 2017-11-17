#include "global_types.h"

//meshes
Tri_Mesh_Info tri_mesh;
Hybrid_Mesh_Info hyrbid_mesh;
Hex_Mesh_Info hex_mesh;
Frame_Info FrameI;
Singularity_Info SingularityI;

//components
vector<Component_Sheet> all_com_sheets;
vector<Component_Chord> all_com_chords;
int Num_Group;
vector<vector<int> > All_groups;
vector<Structure_Coverage_Group> All_decompositions;
int Num_Group_chord;
vector<Chord_Group> All_groups_chord;
vector<Structure_Coverage_Group_chord> All_decompositions_chord;

int Para_Total_N;
int Para_min_N;
