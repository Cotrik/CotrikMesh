#ifndef __IO_H__
#define __IO_H__
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include "constants.h"
#include "global_types.h"
#include "global_functions.h"

using namespace std;
class h_io
{
public:
    int counter;
public:
    h_io(void);
    ~h_io(void);

    void read_triangle_mesh_off(vector<Vertex> &Vs, vector<Triangle> &Ts, const char * fname);
    void read_hex_mesh_off(vector<Hex_V> &Vs, vector<Hex> &Hexs, const char * fname);
    void write_hex_mesh_off(vector<Hex_V> &Vs, vector<Hex> &Hexs, const char * fname);

    void write_VTK(const vector<Hex_V> &Vs, const vector<Hex> &Hexs, const char * fname);
    void read_VTK(vector<Hex_V> &Vs, vector<Hex> &Hexs, const char * fname);

    void save_base_complex_VTK(const char *fname);
    void save_BasecomplexNode_VTK(const char *fname);
    void save_BasecomplexEdge_VTK(const char *fname);
    void save_singularG_VTK(const char *fname);
    void save_singularG_VTK_Color_Nodes(const char *fname);
    void save_singularG_VTK_Color_Edges(const char *fname);
    void save_singularNode_VTK(const char *fname);
    void save_component_VTK(const char *fname);
    void save_sheets_VTK(const char *fname);
    void save_sheet_VTK(const char *fname, const unsigned long sheet_index = 0);
    void save_sheet_hex_VTK(const char *fname, const unsigned long sheet_index = 0);
    void save_sheet_geo_VTK(const char *fname, const unsigned long sheet_index = 0);
};

#endif // __IO_H__
