#ifndef __COMPONENT_CLASS_H__
#define __COMPONENT_CLASS_H__
#include "global_functions.h"
#include "global_types.h"
#include "io.h"
using namespace std;
class component_class
{
public:
    component_class(void);
//component sheets
    void extract_all_component_sheets();
    vector<int> an_edge_area(int e_id, vector<bool> &arrayE_test);
    void build_component_information(vector<int> &edges_middle, vector<int> &edges_left, vector<int> &edges_right,
            vector<int> &faces_middle, vector<int> &faces_left, vector<int> &faces_right, vector<int> &components,
            vector<int> &nodes_left, vector<int> &nodes_right);

    void classify_component_sheets();
    void find_all_decompositions();
//component chords
    void extract_all_component_chords();
//    Component_Chord an_face_area(int f_id, vector<bool> &arrayF_test);
//    bool Is_component_chord_extraordinary(Component_Chord &cc);
//    void build_chords_information(Component_Chord &cc);

    void classify_component_chords();
    void find_all_decompositions_chords();

    ~component_class(void);
};

#endif //__COMPONENT_CLASS_H__
