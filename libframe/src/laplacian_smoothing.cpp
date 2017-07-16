#include "laplacian_smoothing.h"
//#include <math.h>
#include <cmath>

laplacian_smoothing::laplacian_smoothing(void)
{
	Iteration = 5;
	threshold = 0.2;
}
void laplacian_smoothing::laplacian_pipeline(vector<Hex_V> &hvs, vector<Hex_F> &hfs, vector<Hex> &hhs)
{
	for (int i = 0; i < Iteration; i++)
	{
		printf("laplacian %d\n", i);
		laplacian_average_fhs(i, hvs, hfs, hhs);
	}
}
void laplacian_smoothing::laplacian_global_pipeline(vector<Hex_V> &hvs, vector<Hex_F> &hfs, vector<Hex> &hhs)
{
	Iteration = 5;
	for (int i = 0; i < Iteration; i++)
	{
		printf("laplacian %d\n", i);
		project_surface(hvs);
		laplacian_surface_global(hvs, hfs, 0);
		project_surface(hvs);

		calculation_surface_centroid(hvs, hfs, hhs);
		laplacian_inner_volume(hvs, hhs, 0);
		calculation_volume_centroid(hvs, hhs);
	}
}

void laplacian_smoothing::laplacian_average_vs(int iter, vector<Hex_V> &hvs, vector<Hex_F> &hfs, vector<Hex> &hhs)
{
	laplacian_surface(hvs, hfs, 1);

	project_surface(hvs);
	laplacian_inner_volume(hvs, hhs, 1);
}
void laplacian_smoothing::laplacian_average_fhs(int iter, vector<Hex_V> &hvs, vector<Hex_F> &hfs, vector<Hex> &hhs)
{
	//laplacian_surface(hvs,hfs,0);
	// 		char path2[300];
	// 		sprintf(path2,"%s%d%s",path,i,"before_project.mesh");
	// 		hi.write_hex_mesh_mesh(hvs,hhs,path2);
	//project_surface(hvs);
	// 		sprintf(path2,"%s%d%s",path,i,"after_project.mesh");
	// 		hi.write_hex_mesh_mesh(hvs,hhs,path2);

	//calculation_surface_centroid(hvs,hfs,hhs);
	laplacian_inner_volume(hvs, hhs, 1);

	//calculation_volume_centroid(hvs,hhs);
}

bool laplacian_smoothing::project_a_v(vector<Hex_V> &hvs, int which, bool fixornot)
{
	if (hvs[which].where_location == 1 && !fixornot)
	{
		float Min_dis = 10000;
		int which_one = -1;

		for (int j = 0; j < tri_mesh.Vs.size(); j++)
		{
			float dis;
			DISTANCE(dis, hvs[which].v, tri_mesh.Vs[j].v);
			if (Min_dis > dis)
			{
				Min_dis = dis;
				which_one = j;
			}
		}

		int which_one_inT = -1;
		std::vector<Vertex> projected_vs;
		std::vector<float> Diss;

		std::vector<int> allneighborvs = tri_mesh.Vs[which_one].neighborv;
		allneighborvs.push_back(which_one);
		std::vector<int> allneighborts;
		for (int j = 0; j < allneighborvs.size(); j++)
		{
			for (int k = 0; k < tri_mesh.Vs[allneighborvs[j]].neighbort.size(); k++)
			{
				bool havealready = false;
				for (int m = 0; m < allneighborts.size(); m++)
					if (tri_mesh.Vs[allneighborvs[j]].neighbort[k] == allneighborts[m])
						havealready = true;
				if (!havealready)
					allneighborts.push_back(tri_mesh.Vs[allneighborvs[j]].neighbort[k]);
			}
		}
		for (int j = 0; j < allneighborts.size(); j++)
		{
			int tid = allneighborts[j];

			float x, y, z;
			float dis;
			bool istrue = projected_v(tri_mesh.Vs[tri_mesh.Ts[tid].triangle_v[0]].v, tri_mesh.Vs[tri_mesh.Ts[tid].triangle_v[1]].v,
					tri_mesh.Vs[tri_mesh.Ts[tid].triangle_v[2]].v, hvs[which].v, &x, &y, &z, &dis);
			if (istrue)
			{
				Vertex v;
				v.v[0] = x;
				v.v[1] = y;
				v.v[2] = z;
				projected_vs.push_back(v);
				Diss.push_back(dis);
				if (dis < Min_dis)
				{
					Min_dis = dis;
					which_one_inT = projected_vs.size() - 1;
				}
			}
		}
		if (which_one_inT != -1)
		{
			hvs[which].v[0] = projected_vs[which_one_inT].v[0];
			hvs[which].v[1] = projected_vs[which_one_inT].v[1];
			hvs[which].v[2] = projected_vs[which_one_inT].v[2];
		}
		{
			int which_one_v = -1;
			vector<int> vidnVs = tri_mesh.Vs[which_one].neighborv;
			for (int p = 0; p < vidnVs.size(); p++)
			{
				float kk = -1;
				point_line_projection(tri_mesh.Vs[vidnVs[p]].v, tri_mesh.Vs[which_one].v, hvs[which].v, kk);
				if (kk >= 0 && kk <= 1)
				{
					Vertex v_true;
					for (int q = 0; q < 3; q++)
					{
						v_true.v[q] = tri_mesh.Vs[vidnVs[p]].v[q] + kk * (tri_mesh.Vs[which_one].v[q] - tri_mesh.Vs[vidnVs[p]].v[q]);
					}
					float dis;
					DISTANCE(dis, hvs[which].v, v_true.v);
					projected_vs.push_back(v_true);
					Diss.push_back(dis);
					if (dis < Min_dis)
					{
						Min_dis = dis;
						which_one_v = projected_vs.size() - 1;
					}
				}
			}
			if (which_one_v != -1)
			{
				hex_mesh.HVs[which].v[0] = projected_vs[which_one_v].v[0];
				hex_mesh.HVs[which].v[1] = projected_vs[which_one_v].v[1];
				hex_mesh.HVs[which].v[2] = projected_vs[which_one_v].v[2];
			}
		}
	}

	return true;
}
void laplacian_smoothing::project_surface(vector<Hex_V> &hvs)
{
#pragma omp parallel for //multi-thread	for (int i = 0; i < hvs.size(); i++)	{
		project_a_v(hvs, i, hvs[i].fixed);
	}
}
void laplacian_smoothing::laplacian_surface(vector<Hex_V> &hvs, vector<Hex_F> &hfs, int which_type)
{
	for (int i = 0; i < hvs.size(); i++)
	{
		float v[3];
		v[2] = 0;
		v[1] = 0;
		v[0] = 0;

		if (hvs[i].where_location == 1 && !hvs[i].fixed)
		{
			int num_valid = 0;
			if (which_type == 0)
			{
				for (int j = 0; j < hvs[i].neighbor_Fs.size(); j++)
				{
					if (hfs[hvs[i].neighbor_Fs[j]].is_boundary != 1)
						continue;
					num_valid++;
					v[0] += hfs[hvs[i].neighbor_Fs[j]].centroid[0];
					v[1] += hfs[hvs[i].neighbor_Fs[j]].centroid[1];
					v[2] += hfs[hvs[i].neighbor_Fs[j]].centroid[2];
				}
			}
			else if (which_type == 1)
			{
				for (int j = 0; j < hvs[i].neighbor_vs.size(); j++)
				{
					if (hvs[hvs[i].neighbor_vs[j]].where_location != 1)
						continue;
					num_valid++;
					v[0] += hvs[hvs[i].neighbor_vs[j]].v[0];
					v[1] += hvs[hvs[i].neighbor_vs[j]].v[1];
					v[2] += hvs[hvs[i].neighbor_vs[j]].v[2];
				}
			}
			if (num_valid > 0)
			{
				v[0] /= num_valid;
				v[1] /= num_valid;
				v[2] /= num_valid;
			}

			if (which_type == 1)
			{
				v[0] += (v[0] - hvs[i].v[0]) * threshold;
				v[1] += (v[1] - hvs[i].v[1]) * threshold;
				v[2] += (v[2] - hvs[i].v[2]) * threshold;
			}
			hvs[i].v[0] = v[0];
			hvs[i].v[1] = v[1];
			hvs[i].v[2] = v[2];
		}
	}
}
void laplacian_smoothing::laplacian_inner_volume(vector<Hex_V> &hvs, vector<Hex> &hhs, int which_type)
{
	for (int i = 0; i < hvs.size(); i++)
	{
		float v[3];
		v[2] = 0;
		v[1] = 0;
		v[0] = 0;

		if (hvs[i].where_location != 1)
		{
			if (which_type == 0)
			{
// 				vector<int> neighborhs;
// 				for(int j=0;j<hvs[i].neighbor_vs.size();j++)
// 				{
// 					int nvid=hvs[i].neighbor_vs[j];
// 					for(int k=0;k<hvs[nvid].neighbor_Hs.size();k++)
// 					{
// 						int ncid=hvs[nvid].neighbor_Hs[k];
// 						int whichone=set_contain(neighborhs,ncid);
// 						if(whichone==-1)
// 							neighborhs.push_back(ncid);
// 					}
// 				}
// 				for(int j=0;j<neighborhs.size();j++)
// 				{
// 					v[0]+=hhs[neighborhs[j]].centroid[0];
// 					v[1]+=hhs[neighborhs[j]].centroid[1];
// 					v[2]+=hhs[neighborhs[j]].centroid[2];
// 				}
// 				v[0]/=(neighborhs.size());
// 				v[1]/=(neighborhs.size());
// 				v[2]/=(neighborhs.size());

				for (int j = 0; j < hvs[i].neighbor_Hs.size(); j++)
				{
					v[0] += hhs[hvs[i].neighbor_Hs[j]].centroid[0];
					v[1] += hhs[hvs[i].neighbor_Hs[j]].centroid[1];
					v[2] += hhs[hvs[i].neighbor_Hs[j]].centroid[2];
				}
				if (hvs[i].neighbor_Hs.size())
				{
					v[0] /= (hvs[i].neighbor_Hs.size());
					v[1] /= (hvs[i].neighbor_Hs.size());
					v[2] /= (hvs[i].neighbor_Hs.size());
				}
			}
			else if (which_type == 1)
			{
				for (int j = 0; j < hvs[i].neighbor_vs.size(); j++)
				{
					v[0] += hvs[hvs[i].neighbor_vs[j]].v[0];
					v[1] += hvs[hvs[i].neighbor_vs[j]].v[1];
					v[2] += hvs[hvs[i].neighbor_vs[j]].v[2];
				}
				v[0] /= (hvs[i].neighbor_vs.size());
				v[1] /= (hvs[i].neighbor_vs.size());
				v[2] /= (hvs[i].neighbor_vs.size());

				v[0] = hvs[i].v[0] + (v[0] - hvs[i].v[0]) * threshold;
				v[1] = hvs[i].v[1] + (v[1] - hvs[i].v[1]) * threshold;
				v[2] = hvs[i].v[2] + (v[2] - hvs[i].v[2]) * threshold;
			}

			hvs[i].v[0] = v[0];
			hvs[i].v[1] = v[1];
			hvs[i].v[2] = v[2];
		}
	}
}

void laplacian_smoothing::project_surface_global(vector<Hex_V> &hvs)
{
#pragma omp parallel for //multi-thread	for (int i = 0; i < hvs.size(); i++)	{
// 		if(hvs[i].fixed)
// 			continue;
		project_a_v(hvs, i, false);
	}
}
void laplacian_smoothing::laplacian_surface_global(vector<Hex_V> &hvs, vector<Hex_F> &hfs, int which_type)
{
	for (int i = 0; i < hvs.size(); i++)
	{
		float v[3];
		v[2] = 0;
		v[1] = 0;
		v[0] = 0;
		if (hvs[i].fixed)
			continue;

		if (hvs[i].where_location == 1)
		{
			int num_valid = 0;
			if (which_type == 0)
			{
				for (int j = 0; j < hvs[i].neighbor_Fs.size(); j++)
				{
					if (hfs[hvs[i].neighbor_Fs[j]].is_boundary != 1)
						continue;
					num_valid++;
					v[0] += hfs[hvs[i].neighbor_Fs[j]].centroid[0];
					v[1] += hfs[hvs[i].neighbor_Fs[j]].centroid[1];
					v[2] += hfs[hvs[i].neighbor_Fs[j]].centroid[2];
				}
			}
			else if (which_type == 1)
			{
				for (int j = 0; j < hvs[i].neighbor_vs.size(); j++)
				{
					if (hvs[hvs[i].neighbor_vs[j]].where_location != 1)
						continue;
					num_valid++;
					v[0] += hvs[hvs[i].neighbor_vs[j]].v[0];
					v[1] += hvs[hvs[i].neighbor_vs[j]].v[1];
					v[2] += hvs[hvs[i].neighbor_vs[j]].v[2];
				}
			}
			if (num_valid > 0)
			{
				v[0] /= num_valid;
				v[1] /= num_valid;
				v[2] /= num_valid;
			}

			if (which_type == 1)
			{
				v[0] = hvs[i].v[0] + (v[0] - hvs[i].v[0]) * threshold;
				v[1] = hvs[i].v[1] + (v[1] - hvs[i].v[1]) * threshold;
				v[2] = hvs[i].v[2] + (v[2] - hvs[i].v[2]) * threshold;
			}
			hvs[i].v[0] = v[0];
			hvs[i].v[1] = v[1];
			hvs[i].v[2] = v[2];
		}
	}
}

bool laplacian_smoothing::projected_v(float A[3], float B[3], float C[3], float v[3], float *pvx, float *pvy, float *pvz, float *dis)
{
	//normal of plane
	float vec1[3], vec2[3];
	vec1[0] = B[0] - A[0];
	vec2[0] = C[0] - A[0];
	vec1[1] = B[1] - A[1];
	vec2[1] = C[1] - A[1];
	vec1[2] = B[2] - A[2];
	vec2[2] = C[2] - A[2];

	float normal_f[3];
	CROSSVECTOR3(normal_f, vec1, vec2);
	//a,b,c,d for plane fuc
	float a, b, c, d;
	a = normal_f[0];
	b = normal_f[1];
	c = normal_f[2];
	d = -(normal_f[0] * A[0] + normal_f[1] * A[1] + normal_f[2] * A[2]);

	float t = (a * v[0] + b * v[1] + c * v[2] + d) / (a * a + b * b + c * c);
	if (!isfinite(t))
		return false;
	*pvx = v[0] - a * t;
	*pvy = v[1] - b * t;
	*pvz = v[2] - c * t;
	*dis = sqrt((v[0] - *pvx) * (v[0] - *pvx) + (v[1] - *pvy) * (v[1] - *pvy) + (v[2] - *pvz) * (v[2] - *pvz));

	float P[3];
	P[0] = *pvx;
	P[1] = *pvy;
	P[2] = *pvz;
	return isin_triangle(A, B, C, P);
}
bool laplacian_smoothing::projected_v(vector<float> A, vector<float> B, vector<float> C, vector<float> v, float *pvx, float *pvy, float *pvz,
		float *dis)
{
	//normal of plane
	float vec1[3], vec2[3];
	vec1[0] = B[0] - A[0];
	vec2[0] = C[0] - A[0];
	vec1[1] = B[1] - A[1];
	vec2[1] = C[1] - A[1];
	vec1[2] = B[2] - A[2];
	vec2[2] = C[2] - A[2];

	float normal_f[3];
	CROSSVECTOR3(normal_f, vec1, vec2);
	//a,b,c,d for plane fuc
	float a, b, c, d;
	a = normal_f[0];
	b = normal_f[1];
	c = normal_f[2];
	d = -(normal_f[0] * A[0] + normal_f[1] * A[1] + normal_f[2] * A[2]);

	float t = (a * v[0] + b * v[1] + c * v[2] + d) / (a * a + b * b + c * c);
	if (!isfinite(t))
		return false;
	*pvx = v[0] - a * t;
	*pvy = v[1] - b * t;
	*pvz = v[2] - c * t;
	*dis = sqrt((v[0] - *pvx) * (v[0] - *pvx) + (v[1] - *pvy) * (v[1] - *pvy) + (v[2] - *pvz) * (v[2] - *pvz));

	float P[3];
	P[0] = *pvx;
	P[1] = *pvy;
	P[2] = *pvz;
	return isin_triangle(A, B, C, P);
}
bool laplacian_smoothing::isin_triangle(float A[3], float B[3], float C[3], float P[3])
{
	float v0[3];
	v0[0] = C[0] - A[0];
	v0[1] = C[1] - A[1];
	v0[2] = C[2] - A[2];
	float v1[3];
	v1[0] = B[0] - A[0];
	v1[1] = B[1] - A[1];
	v1[2] = B[2] - A[2];
	float v2[3];
	v2[0] = P[0] - A[0];
	v2[1] = P[1] - A[1];
	v2[2] = P[2] - A[2];

	float dot00 = v0[0] * v0[0] + v0[1] * v0[1] + v0[2] * v0[2];
	float dot01 = v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v1[2];
	float dot02 = v0[0] * v2[0] + v0[1] * v2[1] + v0[2] * v2[2];
	float dot11 = v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2];
	float dot12 = v2[0] * v1[0] + v2[1] * v1[1] + v2[2] * v1[2];

	float inverDeno = 1 / (dot00 * dot11 - dot01 * dot01);

	float u = (dot11 * dot02 - dot01 * dot12) * inverDeno;
	if (u < 0 || u > 1) // if u out of range, return directly
	{
		return false;
	}

	float v = (dot00 * dot12 - dot01 * dot02) * inverDeno;
	if (v < 0 || v > 1) // if v out of range, return directly
	{
		return false;
	}

	return u + v <= 1;
}
bool laplacian_smoothing::isin_triangle(vector<float> A, vector<float> B, vector<float> C, float P[3])
{
	float v0[3];
	v0[0] = C[0] - A[0];
	v0[1] = C[1] - A[1];
	v0[2] = C[2] - A[2];
	float v1[3];
	v1[0] = B[0] - A[0];
	v1[1] = B[1] - A[1];
	v1[2] = B[2] - A[2];
	float v2[3];
	v2[0] = P[0] - A[0];
	v2[1] = P[1] - A[1];
	v2[2] = P[2] - A[2];

	float dot00 = v0[0] * v0[0] + v0[1] * v0[1] + v0[2] * v0[2];
	float dot01 = v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v1[2];
	float dot02 = v0[0] * v2[0] + v0[1] * v2[1] + v0[2] * v2[2];
	float dot11 = v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2];
	float dot12 = v2[0] * v1[0] + v2[1] * v1[1] + v2[2] * v1[2];

	float inverDeno = 1 / (dot00 * dot11 - dot01 * dot01);

	float u = (dot11 * dot02 - dot01 * dot12) * inverDeno;
	if (u < 0 || u > 1) // if u out of range, return directly
	{
		return false;
	}

	float v = (dot00 * dot12 - dot01 * dot02) * inverDeno;
	if (v < 0 || v > 1) // if v out of range, return directly
	{
		return false;
	}

	return u + v <= 1;
}
laplacian_smoothing::~laplacian_smoothing(void)
{
}
