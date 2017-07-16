#ifndef __LAPLACIAN_SMOOTHING_H__
#define __LAPLACIAN_SMOOTHING_H__
#include "constants.h"
#include "global_types.h"
#include "global_functions.h"
#include "io.h"
#include <omp.h>
class laplacian_smoothing
{
public:
	int Iteration;
	h_io hi;
	float threshold;
public:
	laplacian_smoothing(void);

	void laplacian_pipeline(vector<Hex_V> &hvs, vector<Hex_F> &hfs, vector<Hex> &hhs);
	void laplacian_global_pipeline(vector<Hex_V> &hvs, vector<Hex_F> &hfs, vector<Hex> &hhs);

	void laplacian_average_vs(int iter, vector<Hex_V> &hvs, vector<Hex_F> &hfs, vector<Hex> &hhs);
	void laplacian_average_fhs(int iter, vector<Hex_V> &hvs, vector<Hex_F> &hfs, vector<Hex> &hhs);

	bool project_a_v(vector<Hex_V> &hvs, int which, bool fixornot);

	void project_surface(vector<Hex_V> &hvs);
	void laplacian_surface(vector<Hex_V> &hvs, vector<Hex_F> &hfs, int which_type);
	void laplacian_inner_volume(vector<Hex_V> &hvs, vector<Hex> &hhs, int which_type);

	void project_surface_global(vector<Hex_V> &hvs);
	void laplacian_surface_global(vector<Hex_V> &hvs, vector<Hex_F> &hfs, int which_type);

	bool projected_v(float A[3], float B[3], float C[3], float v[3], float *pvx, float *pvy, float *pvz, float *dis);
	bool projected_v(vector<float> A, vector<float> B, vector<float> C, vector<float> v, float *pvx, float *pvy, float *pvz, float *dis);

	bool isin_triangle(float A[3], float B[3], float C[3], float P[3]);
	bool isin_triangle(vector<float> A, vector<float> B, vector<float> C, float P[3]);
	~laplacian_smoothing(void);
};

#endif // __LAPLACIAN_SMOOTHING_H__
