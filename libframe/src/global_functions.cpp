#include "global_functions.h"
#include "global_types.h"

void initialization_parameters()
{
	//meshes
	hex_mesh.HVs.clear();
	hex_mesh.HEs.clear();
	hex_mesh.HFs.clear();
	hex_mesh.HHs.clear();
	FrameI.FVs.clear();
	FrameI.FEs.clear();
	FrameI.FFs.clear();
	FrameI.FHs.clear();
	FrameI.HVs.clear();
	FrameI.HEs.clear();
	SingularityI.SVs.clear();
	SingularityI.SEs.clear();
}
;
void initializeVectorT(vector<int> &Ts, int t, int n)
{
	Ts.clear();
	for (int i = 0; i < n; i++)
		Ts.push_back(t);
}
;
void initializeVectorT(vector<bool> &Ts, bool t, int n)
{
	Ts.clear();
	for (int i = 0; i < n; i++)
		Ts.push_back(t);
}
;
void initializeVectorT(vector<float> &Ts, float t, int n)
{
	Ts.clear();
	for (int i = 0; i < n; i++)
		Ts.push_back(t);
}
;

void construct_Es(vector<Hex_V> &HVs, vector<Hex_E> &HEs, vector<Hex> &HHs)
{
	HEs.clear();
	vector<Hex_E> Temp_HEs;
	for (int i = 0; i < HHs.size(); i++)
	{
		Hex_E e;
		//up
		e.startend_Id[0] = HHs[i].V_Ids[0];
		e.startend_Id[1] = HHs[i].V_Ids[1];
		Temp_HEs.push_back(e);

		e.startend_Id[0] = HHs[i].V_Ids[1];
		e.startend_Id[1] = HHs[i].V_Ids[2];
		Temp_HEs.push_back(e);

		e.startend_Id[0] = HHs[i].V_Ids[2];
		e.startend_Id[1] = HHs[i].V_Ids[3];
		Temp_HEs.push_back(e);

		e.startend_Id[0] = HHs[i].V_Ids[3];
		e.startend_Id[1] = HHs[i].V_Ids[0];
		Temp_HEs.push_back(e);
		//down
		e.startend_Id[0] = HHs[i].V_Ids[4];
		e.startend_Id[1] = HHs[i].V_Ids[5];
		Temp_HEs.push_back(e);

		e.startend_Id[0] = HHs[i].V_Ids[5];
		e.startend_Id[1] = HHs[i].V_Ids[6];
		Temp_HEs.push_back(e);

		e.startend_Id[0] = HHs[i].V_Ids[6];
		e.startend_Id[1] = HHs[i].V_Ids[7];
		Temp_HEs.push_back(e);

		e.startend_Id[0] = HHs[i].V_Ids[7];
		e.startend_Id[1] = HHs[i].V_Ids[4];
		Temp_HEs.push_back(e);
		//connect
		e.startend_Id[0] = HHs[i].V_Ids[0];
		e.startend_Id[1] = HHs[i].V_Ids[4];
		Temp_HEs.push_back(e);

		e.startend_Id[0] = HHs[i].V_Ids[1];
		e.startend_Id[1] = HHs[i].V_Ids[5];
		Temp_HEs.push_back(e);

		e.startend_Id[0] = HHs[i].V_Ids[2];
		e.startend_Id[1] = HHs[i].V_Ids[6];
		Temp_HEs.push_back(e);

		e.startend_Id[0] = HHs[i].V_Ids[3];
		e.startend_Id[1] = HHs[i].V_Ids[7];
		Temp_HEs.push_back(e);
	}

	int E_N = 0;
	for (int i = 0; i < Temp_HEs.size(); i++)
	{
		bool havesame = false;

		int id1 = Temp_HEs[i].startend_Id[0];
		int id2 = Temp_HEs[i].startend_Id[1];

		for (int j = 0; j < HVs[id1].neighbor_vs.size(); j++)
		{
			if (HVs[id1].neighbor_vs[j] == id2)
			{
				havesame = true;
				break;
			}
		}
		if (!havesame)
		{
			Temp_HEs[i].index = E_N++;
			Temp_HEs[i].is_boundary = -1;

			for (int j = 0; j < HVs[Temp_HEs[i].startend_Id[0]].neighbor_Hs.size(); j++)
			{
				int h1 = HVs[Temp_HEs[i].startend_Id[0]].neighbor_Hs[j];
				for (int k = 0; k < HVs[Temp_HEs[i].startend_Id[1]].neighbor_Hs.size(); k++)
				{
					int h2 = HVs[Temp_HEs[i].startend_Id[1]].neighbor_Hs[k];
					if (h1 == h2)
						Temp_HEs[i].neighbor_Hs.push_back(h1);
				}
			}

			HEs.push_back(Temp_HEs[i]);
			HVs[id1].neighbor_Es.push_back(Temp_HEs[i].index);
			HVs[id1].neighbor_vs.push_back(id2);

			HVs[id2].neighbor_Es.push_back(Temp_HEs[i].index);
			HVs[id2].neighbor_vs.push_back(id1);
		}
	}

	for (int i = 0; i < HEs.size(); i++)
	{
		for (int j = 0; j < HVs[HEs[i].startend_Id[0]].neighbor_Hs.size(); j++)
		{
			int nh = HVs[HEs[i].startend_Id[0]].neighbor_Hs[j];
			for (int k = 0; k < HVs[HEs[i].startend_Id[1]].neighbor_Hs.size(); k++)
				if (nh == HVs[HEs[i].startend_Id[1]].neighbor_Hs[k])
					HEs[i].neighbor_Hs.push_back(nh);
		}
	}

	for (int i = 0; i < HEs.size(); i++)
	{
		vector<int> neighborhs;
		for (int j = 0; j < HEs[i].neighbor_Hs.size(); j++)
		{
			bool have = false;
			for (int k = 0; k < neighborhs.size(); k++)
				if (HEs[i].neighbor_Hs[j] == neighborhs[k])
					have = true;
			if (!have)
				neighborhs.push_back(HEs[i].neighbor_Hs[j]);
		}
		HEs[i].neighbor_Hs.clear();
		HEs[i].neighbor_Hs = neighborhs;
	}

	for (int i = 0; i < HVs.size(); i++)
	{
		vector<int> neighborvs;
		for (int j = 0; j < HVs[i].neighbor_vs.size(); j++)
		{
			bool have = false;
			for (int k = 0; k < neighborvs.size(); k++)
				if (HVs[i].neighbor_vs[j] == neighborvs[k])
					have = true;
			if (!have)
				neighborvs.push_back(HVs[i].neighbor_vs[j]);
		}
		HVs[i].neighbor_vs.clear();
		HVs[i].neighbor_vs = neighborvs;
	}
}
void construct_Fs(vector<Hex_V> &HVs, vector<Hex_E> &HEs, vector<Hex_F> &HFs, vector<Hex> &HHs)
{
	HFs.clear();
	vector<Hex_F> Temp_HFs;
	int F_N = 0;
	for (int i = 0; i < HHs.size(); i++)
	{
		Hex_F hf[6];

		hf[0].cv_Ids[0] = HHs[i].V_Ids[0];
		hf[0].cv_Ids[1] = HHs[i].V_Ids[1];
		hf[0].cv_Ids[2] = HHs[i].V_Ids[2];
		hf[0].cv_Ids[3] = HHs[i].V_Ids[3];

		hf[1].cv_Ids[0] = HHs[i].V_Ids[4];
		hf[1].cv_Ids[1] = HHs[i].V_Ids[5];
		hf[1].cv_Ids[2] = HHs[i].V_Ids[6];
		hf[1].cv_Ids[3] = HHs[i].V_Ids[7];

		hf[2].cv_Ids[0] = HHs[i].V_Ids[0];
		hf[2].cv_Ids[1] = HHs[i].V_Ids[1];
		hf[2].cv_Ids[2] = HHs[i].V_Ids[5];
		hf[2].cv_Ids[3] = HHs[i].V_Ids[4];

		hf[3].cv_Ids[0] = HHs[i].V_Ids[0];
		hf[3].cv_Ids[1] = HHs[i].V_Ids[4];
		hf[3].cv_Ids[2] = HHs[i].V_Ids[7];
		hf[3].cv_Ids[3] = HHs[i].V_Ids[3];

		hf[4].cv_Ids[0] = HHs[i].V_Ids[3];
		hf[4].cv_Ids[1] = HHs[i].V_Ids[2];
		hf[4].cv_Ids[2] = HHs[i].V_Ids[6];
		hf[4].cv_Ids[3] = HHs[i].V_Ids[7];

		hf[5].cv_Ids[0] = HHs[i].V_Ids[1];
		hf[5].cv_Ids[1] = HHs[i].V_Ids[2];
		hf[5].cv_Ids[2] = HHs[i].V_Ids[6];
		hf[5].cv_Ids[3] = HHs[i].V_Ids[5];

		for (int j = 0; j < 6; j++)
		{
			bool have = false;
			for (int m = 0; m < 4; m++)
			{
				for (int n = 0; n < HVs[hf[j].cv_Ids[m]].neighbor_Fs.size(); n++)
				{
					int F_id = HVs[hf[j].cv_Ids[m]].neighbor_Fs[n];

					bool all_have = true;
					for (int p = 0; p < 4; p++)
					{
						bool exist_v = false;

						for (int q = 0; q < 4; q++)
							if (hf[j].cv_Ids[p] == HFs[F_id].cv_Ids[q])
								exist_v = true;
						if (!exist_v)
							all_have = false;
					}
					if (all_have)
					{
						have = true;
						HFs[F_id].neighbor_Cs.push_back(i);
					}
				}
			}
			if (!have)
			{
				hf[j].index = F_N++;

				for (int k = 0; k < 4; k++)
				{
					int id1 = hf[j].cv_Ids[k];
					int id2 = hf[j].cv_Ids[(k + 1) % 4];
					bool found = false;
					for (int m = 0; m < HVs[id1].neighbor_Es.size(); m++)
					{
						int edge1 = HVs[id1].neighbor_Es[m];
						for (int n = 0; n < HVs[id2].neighbor_Es.size(); n++)
						{
							int edge2 = HVs[id2].neighbor_Es[n];
							if (edge1 == edge2)
							{
								hf[j].ce_Ids[k] = edge1;
								found = true;
							}
							if (found)
								break;
						}
						if (found)
							break;
					}
				}
				hf[j].frame_boundary = -1;
				HFs.push_back(hf[j]);
				HVs[hf[j].cv_Ids[0]].neighbor_Fs.push_back(hf[j].index);
				HVs[hf[j].cv_Ids[1]].neighbor_Fs.push_back(hf[j].index);
				HVs[hf[j].cv_Ids[2]].neighbor_Fs.push_back(hf[j].index);
				HVs[hf[j].cv_Ids[3]].neighbor_Fs.push_back(hf[j].index);

				HFs[HFs.size() - 1].neighbor_Cs.push_back(i);
			}
		}
	}

	for (int i = 0; i < HFs.size(); i++)
	{
		vector<int> neighborHs = HFs[i].neighbor_Cs;
		HFs[i].neighbor_Cs.clear();

		for (int j = 0; j < neighborHs.size(); j++)
		{
			bool already = false;
			for (int k = 0; k < HFs[i].neighbor_Cs.size(); k++)
			{
				if (neighborHs[j] == HFs[i].neighbor_Cs[k])
				{
					already = true;
					break;
				}
			}
			if (!already)
				HFs[i].neighbor_Cs.push_back(neighborHs[j]);
		}
	}

	for (int i = 0; i < HFs.size(); i++)
	{
		for (int j = 0; j < 4; j++)
		{
			bool found = false;
			for (int k = 0; k < HVs[HFs[i].cv_Ids[j]].neighbor_Es.size(); k++)
			{
				int ne = HVs[HFs[i].cv_Ids[j]].neighbor_Es[k];
				for (int m = 0; m < HVs[HFs[i].cv_Ids[(j + 1) % 4]].neighbor_Es.size(); m++)
				{
					if (ne == HVs[HFs[i].cv_Ids[(j + 1) % 4]].neighbor_Es[m])
					{
						HEs[ne].neighbor_Fs.push_back(i);
						HFs[i].ce_Ids[j] = ne;
						found = true;
						break;
					}
				}
				if (found)
					break;
			}
		}
	}

	for (int i = 0; i < HFs.size(); i++)
	{
		for (int j = 0; j < HFs[i].neighbor_Cs.size(); j++)
			HHs[HFs[i].neighbor_Cs[j]].neighbor_FS.push_back(i);
	}

	for (int i = 0; i < HHs.size(); i++)
	{
		vector<int> f_ids = HHs[i].neighbor_FS;
		HHs[i].neighbor_FS.clear();
		for (int j = 0; j < f_ids.size(); j++)
		{
			bool havesame = false;
			for (int k = j + 1; k < f_ids.size(); k++)
				if (f_ids[j] == f_ids[k])
					havesame = true;
			if (!havesame)
			{
				HHs[i].neighbor_FS.push_back(f_ids[j]);
				HHs[i].F_Ids[HHs[i].neighbor_FS.size() - 1] = HFs[f_ids[j]].index;
			}
		}
	}
}
void determine_boundary_info(vector<Hex_V> &HVs, vector<Hex_E> &HEs, vector<Hex_F> &HFs, vector<Hex> &HHs)
{
	for (int i = 0; i < HFs.size(); i++)
	{
		HFs[i].is_boundary = -1;
		if (HFs[i].neighbor_Cs.size() == 1)
		{
			HFs[i].is_boundary = 1;
			for (int j = 0; j < 4; j++)
			{

				HEs[HFs[i].ce_Ids[j]].is_boundary = 1;
				HVs[HFs[i].cv_Ids[j]].where_location = 1;
			}
		}
	}
	for (int i = 0; i < HHs.size(); i++)
	{
		HHs[i].at_boundary = -1;
		for (int j = 0; j < HHs[i].neighbor_FS.size(); j++)
		{
			int fid = HHs[i].neighbor_FS[j];
			if (HFs[fid].is_boundary == 1)
				HHs[i].at_boundary = 1;
		}
	}
}
float average_len(vector<Hex_V> &HVs, vector<Hex_E> &HEs)
{
	float len_aver = 0;
	for (int i = 0; i < HEs.size(); i++)
	{
		float dis;
		DISTANCE(dis, HVs[HEs[i].startend_Id[0]].v, HVs[HEs[i].startend_Id[1]].v);
		len_aver += dis;
	}
	return len_aver / HEs.size();
}
bool update_one_quadcenter(vector<Hex_V> &HVs, vector<Hex_F> &HFs, int which, int bornot)
{
	if (bornot == 0)
	{
		if (HFs[which].is_boundary != 1)
			return false;
	}
	Hex_V v1, v2, v3, v4;
	v1 = HVs[HFs[which].cv_Ids[0]];
	v2 = HVs[HFs[which].cv_Ids[1]];
	v3 = HVs[HFs[which].cv_Ids[2]];
	v4 = HVs[HFs[which].cv_Ids[3]];

	HFs[which].centroid[0] = (v1.v[0] + v2.v[0] + v3.v[0] + v4.v[0]) / 4;
	HFs[which].centroid[1] = (v1.v[1] + v2.v[1] + v3.v[1] + v4.v[1]) / 4;
	HFs[which].centroid[2] = (v1.v[2] + v2.v[2] + v3.v[2] + v4.v[2]) / 4;

	return true;
}
bool update_one_hex_b_center(vector<Hex_V> &HVs, vector<Hex> &HHs, int which, int bornot)
{
	if (bornot == 0)
	{
		if (HHs[which].at_boundary != 1)
			return false;
	}

	Hex_V v1, v2, v3, v4, v5, v6, v7, v8;
	v1 = HVs[HHs[which].V_Ids[0]];
	v2 = HVs[HHs[which].V_Ids[1]];
	v3 = HVs[HHs[which].V_Ids[2]];
	v4 = HVs[HHs[which].V_Ids[3]];
	v5 = HVs[HHs[which].V_Ids[4]];
	v6 = HVs[HHs[which].V_Ids[5]];
	v7 = HVs[HHs[which].V_Ids[6]];
	v8 = HVs[HHs[which].V_Ids[7]];

	HHs[which].centroid[0] = (v1.v[0] + v2.v[0] + v3.v[0] + v4.v[0] + v5.v[0] + v6.v[0] + v7.v[0] + v8.v[0]) / 8;
	HHs[which].centroid[1] = (v1.v[1] + v2.v[1] + v3.v[1] + v4.v[1] + v5.v[1] + v6.v[1] + v7.v[1] + v8.v[1]) / 8;
	HHs[which].centroid[2] = (v1.v[2] + v2.v[2] + v3.v[2] + v4.v[2] + v5.v[2] + v6.v[2] + v7.v[2] + v8.v[2]) / 8;
	return true;
}
void calculation_surface_centroid(vector<Hex_V> &HVs, vector<Hex_F> &HFs, vector<Hex> &HHs)
{
#pragma omp parallel
	{
#pragma omp for //multi-thread		for (int i = 0; i < HFs.size(); i++)
		{
			update_one_quadcenter(HVs, HFs, i, 0);
		}
#pragma omp for //multi-thread		for (int i = 0; i < HHs.size(); i++)
		{
			update_one_hex_b_center(HVs, HHs, i, 0);
		}
	}
}
void calculation_volume_centroid(vector<Hex_V> &HVs, vector<Hex> &HHs)
{
#pragma omp parallel for //multi-thread	for (int i = 0; i < HHs.size(); i++)
	{
		update_one_hex_b_center(HVs, HHs, i, 1);
	}
}
void calculation_centroid(vector<Hex_V> &HVs, vector<Hex_F> &HFs, vector<Hex> &HHs)
{
#pragma omp parallel
	{
#pragma omp for //multi-thread		for (int i = 0; i < HFs.size(); i++)
		{
			update_one_quadcenter(HVs, HFs, i, 1);
		}
#pragma omp for //multi-thread		for (int i = 0; i < HHs.size(); i++)
		{
			update_one_hex_b_center(HVs, HHs, i, 1);
		}
	}
}

void base_complex_e_onering_fs(int eid, vector<int> &fs)
{
	fs.clear();
	int v1 = hex_mesh.HEs[eid].startend_Id[0], v2 = hex_mesh.HEs[eid].startend_Id[1];
	vector<int> fs_op;
	set_exclusion(hex_mesh.HVs[v1].neighbor_Fs, hex_mesh.HVs[v2].neighbor_Fs, fs_op);
	set_exclusion(hex_mesh.HVs[v1].neighbor_Fs, fs_op, fs);
}
bool insideVectorT(vector<int> Ts, int t)
{
	for (int i = 0; i < Ts.size(); i++)
		if (Ts[i] == t)
			return true;
	return false;
}
;

//basic functions
void set_exclusion(vector<int> large_set, vector<int> small_set, vector<int> &result_set)
{
	for (int i = 0; i < large_set.size(); i++)
	{
		bool inside = false;
		for (int j = 0; j < small_set.size(); j++)
		{
			if (small_set[j] == large_set[i])
			{
				inside = true;
				break;
			}
		}
		if (!inside)
			result_set.push_back(large_set[i]);
	}
}
int set_contain(vector<int> large_set, int elm)
{
	for (int j = 0; j < large_set.size(); j++)
		if (elm == large_set[j])
			return j;

	return -1;
}
;
bool set_contain(vector<int> large_set, vector<int> small_set)
{
	for (int i = 0; i < small_set.size(); i++)
	{
		bool inside = false;
		for (int j = 0; j < large_set.size(); j++)
		{
			if (small_set[i] == large_set[j])
			{
				inside = true;
				break;
			}
		}
		if (!inside)
			return false;
	}
	return true;
}
;
bool set_contain(vector<int> large_set, vector<int> small_set, int num)
{
	for (int i = 0; i < small_set.size(); i++)
	{
		bool inside = false;
		for (int j = 0; j < large_set.size(); j++)
		{
			if (small_set[i] == large_set[j])
			{
				num--;
				break;
			}
		}
	}
	if (num == 0)
		return true;
	return false;
}
void set_cross(vector<int> set1, vector<int> set2, vector<int> &result_set)
{
	result_set.clear();
	for (int i = 0; i < set1.size(); i++)
	{
		bool inside = false;
		for (int j = 0; j < set2.size(); j++)
		{
			if (set2[j] == set1[i])
			{
				inside = true;
				break;
			}
		}
		if (inside)
			result_set.push_back(set1[i]);
	}
}
void set_cross_ids(vector<int> set1, vector<int> set2, vector<int> &result_set)
{
	result_set.clear();
	for (int i = 0; i < set1.size(); i++)
	{
		bool inside = false;
		for (int j = 0; j < set2.size(); j++)
		{
			if (set2[j] == set1[i])
			{
				inside = true;
				break;
			}
		}
		if (inside)
			result_set.push_back(i);
	}
}
void set_redundent_clearn(vector<int> &set)
{
	vector<int> set_copy;
	for (int i = 0; i < set.size(); i++)
	{
		bool have = false;
		for (int j = i + 1; j < set.size(); j++)
			if (set[i] == set[j])
				have = true;
		if (!have)
			set_copy.push_back(set[i]);
	}
	set = set_copy;
}

//math
float matrix_determinant(vector<float> &C1, vector<float> &C2, vector<float> &C3, vector<float> &C4)
{
	Matrix4f M(4, 4);
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			switch (i)
			{
			case 0:
				M(i, j) = C1[j];
				break;
			case 1:
				M(i, j) = C2[j];
				break;
			case 2:
				M(i, j) = C3[j];
				break;
			case 3:
				M(i, j) = C4[j];
				break;
			}
		}
	}
	return M.determinant();
}
float sign_of_value(double x)
{
	if (abs(x) < EPS)
		return 0.0;
	if (x < 0)
		return -1.0;
	return 1.0;
}
//vector
void append_vector(vector<int> &v1, vector<int> v2)
{
	for (int i = 0; i < v2.size(); i++)
		v1.push_back(v2[i]);
}
void interpolation_vector(vector<float> v1, vector<float> v2, vector<float> &v12, float w)
{
	for (int i = 0; i < v1.size(); i++)
	{
		v12.push_back(v1[i] + w * (v2[i] - v1[i]));
	}
}
void reverse_vector(vector<int> &vs)
{
	vector<int> tempvs = vs;
	vs.clear();
	for (int i = tempvs.size() - 1; i >= 0; i--)
		vs.push_back(tempvs[i]);
}
//geometry
void point_line_projection(vector<float> v1, vector<float> v2, vector<float> v, float &t)
{
	vector<float> vv1, v21;
	for (int i = 0; i < v1.size(); i++)
	{
		vv1.push_back(v[i] - v1[i]);
		v21.push_back(v2[i] - v1[i]);
	}
	float nv21_2 = v21[0] * v21[0] + v21[1] * v21[1] + v21[2] * v21[2];
	if (abs(nv21_2) >= EPS)
		t = (vv1[0] * v21[0] + vv1[1] * v21[1] + vv1[2] * v21[2]) / nv21_2;
	else
		t = -1;
}
void point_line_projection(float v1[3], float v2[3], float v[3], float &t)
{
	vector<float> vv1, v21;
	for (int i = 0; i < 3; i++)
	{
		vv1.push_back(v[i] - v1[i]);
		v21.push_back(v2[i] - v1[i]);
	}
	float nv21_2 = v21[0] * v21[0] + v21[1] * v21[1] + v21[2] * v21[2];
	if (abs(nv21_2) >= EPS)
		t = (vv1[0] * v21[0] + vv1[1] * v21[1] + vv1[2] * v21[2]) / nv21_2;
	else
		t = -1;
}
void triangle_coordinates(vector<float> v1, vector<float> v2, vector<float> v3, vector<float> v, vector<float> &ws)
{
	float area, area1, area2, area3;
	triangle_area(v2, v3, v, area1);
	ws.push_back(area1);
	triangle_area(v3, v1, v, area2);
	ws.push_back(area2);
	triangle_area(v1, v2, v, area3);
	ws.push_back(area3);
	area = ws[0] + ws[1] + ws[2];
	if (abs(area) > EPS)
	{
		ws[0] /= area;
		ws[1] /= area;
		ws[2] /= area;
	}
	else
	{
		ws[0] = ws[1] = ws[2] = 1.0 / 3;
	}
}
void triangle_coordinates(float v1[3], float v2[3], float v3[3], float v[3], vector<float> &ws)
{
	for (int i = 0; i < 3; i++)
		printf("v1 %f, v2 %f, v3 %f, v %f\n", v1[i], v2[i], v3[i], v[i]);
	float area, area1, area2, area3;
	triangle_area(v2, v3, v, area1);
	ws.push_back(area1);
	triangle_area(v3, v1, v, area2);
	ws.push_back(area2);
	triangle_area(v1, v2, v, area3);
	ws.push_back(area3);

	area = ws[0] + ws[1] + ws[2];
	printf("a1 %f,a2 %f, a3 %f a %f\n", area1, area2, area3, area);
	ws[0] /= area;
	ws[1] /= area;
	ws[2] /= area;
}
void triangle_area(float v1[3], float v2[3], float v3[3], float &area)
{
	float dis1, dis2, dis3;
	DISTANCE(dis1, v1, v2);
	DISTANCE(dis2, v1, v3);
	DISTANCE(dis3, v3, v2);
	float len = (dis1 + dis2 + dis3) / 2;
	area = sqrt(len * (len - dis1) * (len - dis2) * (len - dis3));
}
void triangle_area(vector<float> v1, vector<float> v2, vector<float> v3, float &area)
{
	float dis1, dis2, dis3;
	DISTANCE(dis1, v1, v2);
	DISTANCE(dis2, v1, v3);
	DISTANCE(dis3, v3, v2);
	float len = (dis1 + dis2 + dis3) / 2;
	area = sqrt(len * (len - dis1) * (len - dis2) * (len - dis3));
}

//topology
int find_opposite_f(int F_f, int F_h)
{
	Hex h = hex_mesh.HHs[F_h];
	Hex_F f_ = hex_mesh.HFs[F_f];
	for (int i = 0; i < h.neighbor_FS.size(); i++)
	{
		int fid_temp = h.neighbor_FS[i];

		Hex_F f = hex_mesh.HFs[fid_temp];

		bool is_this_f = true;
		for (int k = 0; k < 4; k++)
		{
			int vid_temp = hex_mesh.HFs[fid_temp].cv_Ids[k];
			for (int m = 0; m < 4; m++)
			{
				int vid_temp2 = hex_mesh.HFs[F_f].cv_Ids[m];
				if (vid_temp == vid_temp2)
				{
					is_this_f = false;
					break;
				}
			}
			if (!is_this_f)
				break;
		}
		if (is_this_f)
			return fid_temp;
	}
	return -1;
}
int find_opposite_f_frame(int F_f, int F_h)
{
	Frame_H h = FrameI.FHs[F_h];
	Frame_F f_ = FrameI.FFs[F_f];
	for (int i = 0; i < h.neighbor_FS.size(); i++)
	{
		int fid_temp = h.neighbor_FS[i];

		Frame_F f = FrameI.FFs[fid_temp];

		bool is_this_f = true;
		for (int k = 0; k < 4; k++)
		{
			int vid_temp = FrameI.FFs[fid_temp].fv_Ids[k];
			for (int m = 0; m < 4; m++)
			{
				int vid_temp2 = FrameI.FFs[F_f].fv_Ids[m];
				if (vid_temp == vid_temp2)
				{
					is_this_f = false;
					break;
				}
			}
			if (!is_this_f)
				break;
		}
		if (is_this_f)
			return fid_temp;
	}
	return -1;
}
