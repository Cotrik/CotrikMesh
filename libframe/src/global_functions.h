#ifndef __GLOBAL_FUNCTION_H__
#define __GLOBAL_FUNCTION_H__
#include "global_types.h"
#include "Eigen/Dense"
#include <omp.h>
#include "io.h"
#include <iostream>
#include <map>
#include <set>
#include <queue>
using namespace Eigen;
using namespace std;

void initialization_parameters();

void construct_Es(vector<Hex_V> &HVs, vector<Hex_E> &HEs, vector<Hex> &HHs);
void construct_Fs(vector<Hex_V> &HVs, vector<Hex_E> &HEs, vector<Hex_F> &HFs, vector<Hex> &HHs);
void determine_boundary_info(vector<Hex_V> &HVs, vector<Hex_E> &HEs, vector<Hex_F> &HFs, vector<Hex> &HHs);
float average_len(vector<Hex_V> &HVs, vector<Hex_E> &HEs);
void calculation_surface_centroid(vector<Hex_V> &HVs, vector<Hex_F> &HFs, vector<Hex> &HHs);
void calculation_volume_centroid(vector<Hex_V> &HVs, vector<Hex> &HHs);
void calculation_centroid(vector<Hex_V> &HVs, vector<Hex_F> &HFs, vector<Hex> &HHs);

void base_complex_e_onering_fs(int eid, vector<int> &fs); //base_complex

void initializeVectorT(vector<int> &Ts, int t, int n);
void initializeVectorT(vector<bool> &Ts, bool t, int n);
void initializeVectorT(vector<float> &Ts, float t, int n);

bool insideVectorT(vector<int> Ts, int t);
//basic functions
void set_exclusion(vector<int> large_set, vector<int> small_set, vector<int> &result_set);
int set_contain(vector<int> large_set, int elm);
bool set_contain(vector<int> large_set, vector<int> small_set);
bool set_contain(vector<int> large_set, vector<int> small_set, int num); //contain num elements of small_set
void set_cross(vector<int> set1, vector<int> set2, vector<int> &result_set);
void set_cross_ids(vector<int> set1, vector<int> set2, vector<int> &result_set);
void set_redundent_clearn(vector<int> &set);
//math
float matrix_determinant(vector<float> &C1, vector<float> &C2, vector<float> &C3, vector<float> &C4);
float sign_of_value(double x);
//vector
void append_vector(vector<int> &v1, vector<int> v2);
void interpolation_vector(vector<float> v1, vector<float> v2, vector<float> &v12, float w);
void reverse_vector(vector<int> &vs);
//geometry
void point_line_projection(vector<float> v1, vector<float> v2, vector<float> v, float &t);
void point_line_projection(float v1[3], float v2[3], float v[3], float &t);
void triangle_coordinates(vector<float> v1, vector<float> v2, vector<float> v3, vector<float> v, vector<float> &ws);
void triangle_coordinates(float v1[3], float v2[3], float v3[3], float v[3], vector<float> &ws);
void triangle_area(vector<float> v1, vector<float> v2, vector<float> v3, float &area);
void triangle_area(float v1[3], float v2[3], float v3[3], float &area);

//topology
int find_opposite_f(int F_f, int F_h);
int find_opposite_f_frame(int F_f, int F_h);

#endif // __GLOBAL_FUNCTION_H__
