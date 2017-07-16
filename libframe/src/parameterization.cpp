#include "parameterization.h"

parameterization::parameterization(void)
{
	Escalar = 1.0;
}
void parameterization::parameterization_main(char *path)
{
	if (hex_mesh.average_e_len <= 1.5)
	{
		enlarged = true;
		Escalar = 1000;
		for (int i = 0; i < hex_mesh.HVs.size(); i++)
		{
			hex_mesh.HVs[i].v[0] *= Escalar;
			hex_mesh.HVs[i].v[1] *= Escalar;
			hex_mesh.HVs[i].v[2] *= Escalar;
		}
		hex_mesh.average_e_len *= Escalar;
		printf("Enlarged\n");
	}
	determine_parametric_numbers();

	vector<Parameterization_Cell> PCs;
	for (int i = 0; i < FrameI.FHs.size(); i++)
	{
		Parameterization_Cell pc;
		PCs.push_back(pc);

		Us.push_back(0);
		Vs.push_back(0);
		Ws.push_back(0);
	}

	clock_t start_time = clock();

#pragma omp parallel for //multi-thread	for (int i = 0; i < FrameI.FHs.size(); i++)
	{
		parameterization_cuboid(i, PCs[i]);
		//parameterization_cuboid_tetgen(i,PCs[i]);
	}
	clock_t end_time = clock();
	std::cout << "Running time is: " << static_cast<double>(end_time - start_time) / CLOCKS_PER_SEC * 1000 << "ms" << std::endl;

	printf("produce final hex-mesh\n");
	produce_hex_mesh(PCs, path);
}

void parameterization::determine_parametric_numbers()
{
	//extract edge groups.
	vector<bool> arrayE_test;
	for (int i = 0; i < FrameI.FEs.size(); i++)
	{
		arrayE_test.push_back(false);
		es_Flags.push_back(-1);
		L_para_Ns.push_back(0);
	}
	vector<vector<int>> All_es_areas;
	while (true)
	{
		int start_eid = -1;
		for (int i = 0; i < arrayE_test.size(); i++)
		{
			if (!arrayE_test[i])
			{
				start_eid = i;
				break;
			}
		}
		if (start_eid == -1)
			break;
		All_es_areas.push_back(an_edge_area(start_eid, arrayE_test));
	}
	//solve weights
	vector<float> All_es_Weights;

	for (int i = 0; i < All_es_areas.size(); i++)
	{
		float w = 0;
		for (int j = 0; j < All_es_areas[i].size(); j++)
		{
			//weight calculation 1
			//w+=FrameI.FEs[All_es_areas[i][j]].vs_link.size()-1;
			//weight calculation 2
			for (int k = 0; k < FrameI.FEs[All_es_areas[i][j]].vs_link.size() - 1; k++)
			{
				int v1 = FrameI.FEs[All_es_areas[i][j]].vs_link[k];
				int v2 = FrameI.FEs[All_es_areas[i][j]].vs_link[k + 1];
				float dis;
				DISTANCE(dis, hex_mesh.HVs[v1].v, hex_mesh.HVs[v2].v);
				w += dis;
			}
			es_Flags[All_es_areas[i][j]] = i;
		}
		w /= hex_mesh.average_e_len;

		w /= All_es_areas[i].size();
		All_es_Weights.push_back(w);
	}
	//sum up cuboid hex elements
	float Sum_cuboid_hexes_N = 0;
	for (int i = 0; i < FrameI.FHs.size(); i++)
	{
		vector<int> tags;
		for (int j = 0; j < 12; j++)
		{
			if (!insideVectorT(tags, es_Flags[FrameI.FHs[i].neighbor_ES[j]]))
				tags.push_back(es_Flags[FrameI.FHs[i].neighbor_ES[j]]);
		}
		float cuboid_hexes = 1;
		for (int j = 0; j < tags.size(); j++)
		{
			cuboid_hexes *= All_es_Weights[tags[j]];
		}
		Sum_cuboid_hexes.push_back(cuboid_hexes);
		Sum_cuboid_hexes_N += cuboid_hexes;
	}

	float base_w = pow((double) Para_Total_N / Sum_cuboid_hexes_N, (double) 1.0 / 3.0);

	for (int i = 0; i < L_para_Ns.size(); i++)
	{
		L_para_Ns[i] = int(All_es_Weights[es_Flags[i]] * base_w + 0.5);			//
		if (L_para_Ns[i] <= Para_min_N - 1)
			L_para_Ns[i] = Para_min_N;
// 		else if(L_para_Ns[i]>100)
// 			L_para_Ns[i]=40;
// 		else if(L_para_Ns[i]>40)
// 			L_para_Ns[i]=25;
	}
}
vector<int> parameterization::an_edge_area(int e_id, vector<bool> &arrayE_test)
{
	vector<int> es, candiates, candidates_temp;
	candiates.push_back(e_id);
	bool stillhave = true;
	while (stillhave)
	{
		for (int i = 0; i < candiates.size(); i++)
		{
			e_id = candiates[i];
			if (arrayE_test[e_id])
				continue;
			es.push_back(e_id);

			arrayE_test[e_id] = true;

			for (int j = 0; j < FrameI.FEs[e_id].neighbor_Fs.size(); j++)
			{
				int nf_id = FrameI.FEs[e_id].neighbor_Fs[j];
				int which_e_ind;
				for (int k = 0; k < 4; k++)
				{
					if (FrameI.FFs[nf_id].fe_Ids[k] == e_id)
					{
						which_e_ind = k;
						break;
					}
				}

				int op_eid = FrameI.FFs[nf_id].fe_Ids[(which_e_ind + 2) % 4];
				if (!arrayE_test[op_eid])
					candidates_temp.push_back(op_eid);
			}
		}
		if (candidates_temp.size() == 0)
			stillhave = false;
		else
		{
			candiates.clear();
			candiates = candidates_temp;
			candidates_temp.clear();
		}
	}

	return es;
}
float corners[8][3] =
{
{ 0, 0, 0 },
{ 0, 0, 1 },
{ 0, 1, 1 },
{ 0, 1, 0 },
{ 1, 0, 0 },
{ 1, 0, 1 },
{ 1, 1, 1 },
{ 1, 1, 0 } };

void parameterization::parameterization_cuboid_tetgen(int C_Id, Parameterization_Cell &pc)
{
	vector<int> Vs_arrayIds;
	vector<Vertex> Vs_patch;
	vector<Triangle> Fs_patch;
	cuboid2triangle_mesh(C_Id, Vs_arrayIds, Vs_patch, Fs_patch);

	vector<int> VTs_arrayIds;
	vector<Vertex> VTs_patch;
	vector<Edge> ETs_patch;
	vector<Triangle> FTs_patch;
	vector<Tet> TETs_patch;
	translate_to_tetmesh_tetgen(Vs_arrayIds, Vs_patch, Fs_patch, VTs_arrayIds, VTs_patch, FTs_patch, TETs_patch);

	h_io io;

	for (int i = 0; i < VTs_patch.size(); i++)
	{
		vector<float> uvw;
		pc.UVW_coords.push_back(uvw);
	}

	//point
	//printf("8 Corners UVW\n");
	for (int i = 0; i < 8; i++)
	{
		vector<float> uvw;
		uvw.push_back(corners[i][0]);
		uvw.push_back(corners[i][1]);
		uvw.push_back(corners[i][2]);
		pc.UVW_coords[VTs_arrayIds[FrameI.FVs[FrameI.FHs[C_Id].FV_Ids[i]].index_hex]] = uvw;

		//printf("uvw: %f,%f,%f\n",uvw[0],uvw[1],uvw[2]);
	}
	//curve
	//printf("12 Edges UVW\n");
	for (int i = 0; i < 12; i++)
	{
		int fe = FrameI.FHs[C_Id].neighbor_ES[i];
		float total_len = 0;
		vector<float> ratios;
		for (int j = 0; j < FrameI.FEs[fe].vs_link.size() - 1; j++)
		{
			int v1, v2;
			v1 = FrameI.FEs[fe].vs_link[j];
			v2 = FrameI.FEs[fe].vs_link[j + 1];
			float dis;
			DISTANCE(dis, hex_mesh.HVs[v1].v, hex_mesh.HVs[v2].v);
			total_len += dis;
			ratios.push_back(total_len);
		}
		int sv_id = FrameI.FEs[fe].vs_link[0], ev_id = FrameI.FEs[fe].vs_link[FrameI.FEs[fe].vs_link.size() - 1];
		vector<float> startuvw = pc.UVW_coords[VTs_arrayIds[sv_id]], enduvw = pc.UVW_coords[VTs_arrayIds[ev_id]];

		for (int j = 1; j < FrameI.FEs[fe].vs_link.size() - 1; j++)
		{
			vector<float> uvw;
			uvw.push_back(startuvw[0] + ratios[j - 1] / total_len * (enduvw[0] - startuvw[0]));
			uvw.push_back(startuvw[1] + ratios[j - 1] / total_len * (enduvw[1] - startuvw[1]));
			uvw.push_back(startuvw[2] + ratios[j - 1] / total_len * (enduvw[2] - startuvw[2]));
			pc.UVW_coords[VTs_arrayIds[FrameI.FEs[fe].vs_link[j]]] = uvw;

			//printf("uvw: %f,%f,%f\n",uvw[0],uvw[1],uvw[2]);
		}
	}
	//face
	for (int i = 0; i < 6; i++)
	{
		//printf("Face %d UVW\n",i);
		D2_parameterization_tetgen(FrameI.FHs[C_Id].neighbor_FS[i], pc, Vs_arrayIds, Vs_patch, Fs_patch);
	}
	//volume
	//printf("Component %d UVW\n",C_Id);

	//	clock_t start_time=clock();
	//watch->startTimer();
	D3_parameterization_tetgen(C_Id, pc, VTs_patch, ETs_patch, FTs_patch, TETs_patch);
	//watch->stopTimer();
	// 	clock_t end_time=clock();
	// 	std::cout<< "Running time is: "<<static_cast<double>(end_time-start_time)/CLOCKS_PER_SEC*1000<<"ms"<<std::endl;

	//printf("Parameteric of component %d UVW\n",C_Id);
	parametric_coords(C_Id, pc, VTs_patch, TETs_patch);
}
void parameterization::cuboid2triangle_mesh(int C_Id, vector<int> &Vs_arrayIds, vector<Vertex> &vs, vector<Triangle> &fs)
{
	for (int i = 0; i < hex_mesh.HVs.size(); i++)
		Vs_arrayIds.push_back(-2);

// 	for(int j=0;j<FrameI.FHs[C_Id].hs_net.size();j++)
// 	{
// 		int hid=FrameI.FHs[C_Id].hs_net[j];
// 		for(int k=0;k<8;k++)
// 			Vs_arrayIds[hex_mesh.HHs[hid].V_Ids[k]]=-1;
// 	}
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < FrameI.FFs[FrameI.FHs[C_Id].neighbor_FS[i]].hfs_net_another.size(); j++)
		{
			int fid = FrameI.FFs[FrameI.FHs[C_Id].neighbor_FS[i]].hfs_net_another[j];
			for (int k = 0; k < 4; k++)
				Vs_arrayIds[hex_mesh.HFs[fid].cv_Ids[k]] = -1;
		}
	}
	for (int i = 0; i < Vs_arrayIds.size(); i++)
	{
		if (Vs_arrayIds[i] == -1)
		{
			Vertex v;
			v.v[0] = hex_mesh.HVs[i].v[0];
			v.v[1] = hex_mesh.HVs[i].v[1];
			v.v[2] = hex_mesh.HVs[i].v[2];
			v.index = vs.size();
			v.is_boundary = -1;
			v.on_base_complex = 0;
			v.which_F_face = -1;
			vs.push_back(v);
			Vs_arrayIds[i] = v.index;
		}
	}

	for (int i = 0; i < 12; i++)
	{
		for (int j = 0; j < FrameI.FEs[FrameI.FHs[C_Id].neighbor_ES[i]].vs_link.size() - 1; j++)
		{
			int v1, v2;
			v1 = FrameI.FEs[FrameI.FHs[C_Id].neighbor_ES[i]].vs_link[j];
			v1 = Vs_arrayIds[v1];
			v2 = FrameI.FEs[FrameI.FHs[C_Id].neighbor_ES[i]].vs_link[j + 1];
			v2 = Vs_arrayIds[v2];

			vs[v1].on_base_complex = 1;
			vs[v2].on_base_complex = 1;
		}
	}

	int F_ind = 0;
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < FrameI.FFs[FrameI.FHs[C_Id].neighbor_FS[i]].hfs_net_another.size(); j++)
		{
			int fid = FrameI.FFs[FrameI.FHs[C_Id].neighbor_FS[i]].hfs_net_another[j];
			int v1, v2, v3, v4;
			v1 = Vs_arrayIds[hex_mesh.HFs[fid].cv_Ids[0]];
			v2 = Vs_arrayIds[hex_mesh.HFs[fid].cv_Ids[1]];
			v3 = Vs_arrayIds[hex_mesh.HFs[fid].cv_Ids[2]];
			v4 = Vs_arrayIds[hex_mesh.HFs[fid].cv_Ids[3]];

			Triangle t1, t2;
			t1.index = F_ind++;
			t2.index = F_ind++;
			t1.triangle_v[0] = v1;
			t1.triangle_v[1] = v2;
			t1.triangle_v[2] = v3;
			t2.triangle_v[0] = v1;
			t2.triangle_v[1] = v3;
			t2.triangle_v[2] = v4;
			t1.which_F_face = FrameI.FHs[C_Id].neighbor_FS[i];
			t2.which_F_face = FrameI.FHs[C_Id].neighbor_FS[i];

			vs[v1].is_boundary = 1;
			vs[v2].is_boundary = 1;
			vs[v3].is_boundary = 1;
			vs[v4].is_boundary = 1;
			fs.push_back(t1);
			fs.push_back(t2);
			vs[v1].neighbort.push_back(t1.index);
			vs[v1].neighbort.push_back(t2.index);
			vs[v1].neighborv.push_back(v2);
			vs[v1].neighborv.push_back(v3);
			vs[v1].neighborv.push_back(v4);

			vs[v2].neighbort.push_back(t1.index);
			vs[v2].neighborv.push_back(v1);
			vs[v2].neighborv.push_back(v3);

			vs[v3].neighbort.push_back(t1.index);
			vs[v3].neighbort.push_back(t2.index);
			vs[v3].neighborv.push_back(v1);
			vs[v3].neighborv.push_back(v2);
			vs[v3].neighborv.push_back(v4);

			vs[v4].neighbort.push_back(t2.index);
			vs[v4].neighborv.push_back(v1);
			vs[v4].neighborv.push_back(v3);
		}
	}
	for (int i = 0; i < vs.size(); i++)
		set_redundent_clearn(vs[i].neighborv);
}
void parameterization::translate_to_tetmesh_tetgen(vector<int> &Vs_arrayIds, vector<Vertex> &vs, vector<Triangle> &fs, vector<int> &VTs_arrayIds,
		vector<Vertex> &vts, vector<Triangle> &fts, vector<Tet> &tets)
{
	tetgenio in, out;
	tetgenio::facet *f;
	tetgenio::polygon *p;

	// All indices start from 1.
	in.firstnumber = 0;

	in.numberofpoints = vs.size();
	in.pointlist = new REAL[in.numberofpoints * 3];
	in.pointmarkerlist = new int[in.numberofpoints];
	for (int i = 0; i < in.numberofpoints; i++)
	{
		in.pointlist[3 * i + 0] = vs[i].v[0];
		in.pointlist[3 * i + 1] = vs[i].v[1];
		in.pointlist[3 * i + 2] = vs[i].v[2];
		in.pointmarkerlist[i] = 1;
	}

	in.numberoffacets = fs.size();
	in.facetlist = new tetgenio::facet[in.numberoffacets];
	in.facetmarkerlist = new int[in.numberoffacets];

	for (int i = 0; i < in.numberoffacets; i++)
	{
		// Facet 1. The leftmost facet.
		f = &in.facetlist[i];
		f->numberofpolygons = 1;
		f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
		f->numberofholes = 0;
		f->holelist = NULL;
		p = &f->polygonlist[0];
		p->numberofvertices = 3;
		p->vertexlist = new int[p->numberofvertices];
		p->vertexlist[0] = fs[i].triangle_v[0];
		p->vertexlist[1] = fs[i].triangle_v[1];
		p->vertexlist[2] = fs[i].triangle_v[2];

		// Set 'in.facetmarkerlist'
		in.facetmarkerlist[i] = 1;

		//in.pointmarkerlist[fs[i].triangle_v[0]]=1;
		//in.pointmarkerlist[fs[i].triangle_v[1]]=1;
		//in.pointmarkerlist[fs[i].triangle_v[2]]=1;
	}
	// Tetrahedralize the PLC. Switches are chosen to read a PLC (p),
	//   do quality mesh generation (q) with a specified quality bound
	//   (1.414), and apply a maximum volume constraint (a0.1).

	//tetrahedralize(reinterpret_cast<tetgenbehavior *>("p"), &in, &out);//pq1.414Va0.0001;q1.1
	//tetrahedralize("pYq1.1V", &in, &out);//pq1.414Va0.0001;q1.1
	printf("vs before %d; vs after %d\n", vs.size(), out.numberofpoints);

	for (int i = 0; i < out.numberofpoints; i++)
	{
		Vertex v;
		v.index = i;
		v.v[0] = out.pointlist[3 * i + 0];
		v.v[1] = out.pointlist[3 * i + 1];
		v.v[2] = out.pointlist[3 * i + 2];
		v.is_boundary = vs[i].is_boundary;
		vts.push_back(v);
	}
	for (int i = 0; i < out.numberoftetrahedra; i++)
	{
		Tet tet;
		tet.index = i;
		tet.vs.push_back(out.tetrahedronlist[4 * i + 0]);
		tet.vs.push_back(out.tetrahedronlist[4 * i + 1]);
		tet.vs.push_back(out.tetrahedronlist[4 * i + 2]);
		tet.vs.push_back(out.tetrahedronlist[4 * i + 3]);

		if (vts[tet.vs[0]].is_boundary == 1 && vts[tet.vs[1]].is_boundary == 1 && vts[tet.vs[2]].is_boundary == 1 && vts[tet.vs[3]].is_boundary == 1)
			continue;
		tets.push_back(tet);
	}
	//triangles for a tet
	for (int i = 0; i < tets.size(); i++)
	{
		Triangle t;
		for (int j = 0; j < 4; j++)
		{
			t.triangle_v[0] = tets[i].vs[(j + 1) % 4];
			t.triangle_v[1] = tets[i].vs[(j + 2) % 4];
			t.triangle_v[2] = tets[i].vs[(j + 3) % 4];
			t.index = fts.size();

			vector<int> result1, result2;
			set_cross(vts[t.triangle_v[0]].neighbort, vts[t.triangle_v[1]].neighbort, result1);

			set_cross(result1, vts[t.triangle_v[2]].neighbort, result2);
			if (!result2.size())
			{
				fts.push_back(t);
				vts[t.triangle_v[0]].neighbort.push_back(t.index);
				vts[t.triangle_v[1]].neighbort.push_back(t.index);
				vts[t.triangle_v[2]].neighbort.push_back(t.index);
				tets[i].ts.push_back(t.index);
				vts[tets[i].vs[j]].neighbor_ot.push_back(t.index);
			}
			else
			{
				tets[i].ts.push_back(result2[0]);
				vts[tets[i].vs[j]].neighbor_ot.push_back(result2[0]);
			}
			vts[tets[i].vs[j]].neighborv.push_back(t.triangle_v[0]);
			vts[tets[i].vs[j]].neighborv.push_back(t.triangle_v[1]);
			vts[tets[i].vs[j]].neighborv.push_back(t.triangle_v[2]);

			vts[tets[i].vs[j]].neighbortet.push_back(i);
		}
	}
	for (int i = 0; i < vts.size(); i++)
	{
		set_redundent_clearn(vts[i].neighborv);
	}
	VTs_arrayIds = Vs_arrayIds;
}
struct Lamda
{
	vector<float> lamda_ks;
};
bool parameterization::D2_parameterization_tetgen(int F_Id, Parameterization_Cell &pc, vector<int> &VTs_arrayIds, vector<Vertex> &vts,
		vector<Triangle> &fts)
{
	//re-indexing 2D surface patch
	vector<int> VTs_arrayIds_patch, VTs_arrayIds_patch_reverse;
	vector<Vertex> vts_patch;
	vector<Vertex> vts_patch_Inner;
	vector<Vertex> vts_patch_Outer;
	vector<Triangle> fts_patch;

	for (int i = 0; i < vts.size(); i++)
		VTs_arrayIds_patch.push_back(-2);

	for (int i = 0; i < fts.size(); i++)
	{
		if (fts[i].which_F_face == F_Id)
		{
			Triangle t;
			t.triangle_v[0] = fts[i].triangle_v[0];
			t.triangle_v[1] = fts[i].triangle_v[1];
			t.triangle_v[2] = fts[i].triangle_v[2];
			t.index = fts_patch.size();
			fts_patch.push_back(t);
			for (int j = 0; j < 3; j++)
			{
				if (VTs_arrayIds_patch[fts[i].triangle_v[j]] != -2)
					continue;
				VTs_arrayIds_patch[fts[i].triangle_v[j]] = -1;
				Vertex v;
				v.v[0] = vts[fts[i].triangle_v[j]].v[0];
				v.v[1] = vts[fts[i].triangle_v[j]].v[1];
				v.v[2] = vts[fts[i].triangle_v[j]].v[2];
				v.index = fts[i].triangle_v[j];
				v.on_base_complex = vts[fts[i].triangle_v[j]].on_base_complex;

				if (v.on_base_complex)
					vts_patch_Outer.push_back(v);
				else
					vts_patch_Inner.push_back(v);

				VTs_arrayIds_patch_reverse.push_back(-1);
			}
		}
	}
	if (!vts_patch_Inner.size())
		return false;

	for (int i = 0; i < vts_patch_Inner.size(); i++)
	{
		VTs_arrayIds_patch_reverse[i] = vts_patch_Inner[i].index;
		VTs_arrayIds_patch[vts_patch_Inner[i].index] = i;
		vts_patch_Inner[i].index = i;
		vts_patch.push_back(vts_patch_Inner[i]);
	}
	for (int i = 0; i < vts_patch_Outer.size(); i++)
	{
		VTs_arrayIds_patch_reverse[vts_patch_Inner.size() + i] = vts_patch_Outer[i].index;

		VTs_arrayIds_patch[vts_patch_Outer[i].index] = vts_patch_Inner.size() + i;
		vts_patch_Outer[i].index = vts_patch_Inner.size() + i;
		vts_patch.push_back(vts_patch_Outer[i]);
	}

	for (int i = 0; i < fts_patch.size(); i++)
	{
		int v1, v2, v3;
		v1 = VTs_arrayIds_patch[fts_patch[i].triangle_v[0]];
		v2 = VTs_arrayIds_patch[fts_patch[i].triangle_v[1]];
		v3 = VTs_arrayIds_patch[fts_patch[i].triangle_v[2]];

		fts_patch[i].triangle_v[0] = v1;
		fts_patch[i].triangle_v[1] = v2;
		fts_patch[i].triangle_v[2] = v3;
		vts_patch[v1].neighbort.push_back(i);
		vts_patch[v2].neighbort.push_back(i);
		vts_patch[v3].neighbort.push_back(i);

		vts_patch[v1].neighborv.push_back(v2);
		vts_patch[v1].neighborv.push_back(v3);
		vts_patch[v2].neighborv.push_back(v1);
		vts_patch[v2].neighborv.push_back(v3);
		vts_patch[v3].neighborv.push_back(v1);
		vts_patch[v3].neighborv.push_back(v2);
	}
	for (int i = 0; i < vts_patch.size(); i++)
	{
		set_redundent_clearn(vts_patch[i].neighborv);
		//orient vi's neighborvs
		vector<int> orient_nvs;
		int bef = -1, cur = vts_patch[i].neighborv[0];
		orient_nvs.push_back(cur);
		for (int j = 1; j < vts_patch[i].neighborv.size(); j++)
		{
			for (int k = 0; k < vts_patch[i].neighborv.size(); k++)
			{
				if (set_contain(vts_patch[cur].neighborv, vts_patch[i].neighborv[k]) != -1 && vts_patch[i].neighborv[k] != bef)
				{
					orient_nvs.push_back(vts_patch[i].neighborv[k]);
					bef = cur;
					cur = vts_patch[i].neighborv[k];
					break;
				}
			}
		}
		vts_patch[i].neighborv = orient_nvs;
	}
	//cal_lamdas
	vector<Lamda> LamdasforVS;
	for (int i = 0; i < vts_patch_Inner.size(); i++)
	{
		Lamda ls;
		std::vector<float> angle;
		std::vector<float> weights;
		for (int j = 0; j < vts_patch[i].neighborv.size(); j++)
		{
			float dis;
			DISTANCE(dis, vts_patch[vts_patch[i].neighborv[j]].v, vts_patch[vts_patch[i].neighborv[(j + 1) % vts_patch[i].neighborv.size()]].v);
			float disL;
			DISTANCE(disL, vts_patch[vts_patch[i].neighborv[j]].v, vts_patch[i].v);
			float disR;
			DISTANCE(disR, vts_patch[i].v, vts_patch[vts_patch[i].neighborv[(j + 1) % vts_patch[i].neighborv.size()]].v);

			float dis_minus = (disL * disL + disR * disR - dis * dis);
			float divide = dis_minus / (2 * disL * disR);
			if (divide < -1)
				divide = -1;
			if (divide > 1)
				divide = 1;
			float angle_ = acos(divide);
			angle.push_back(angle_);
			weights.push_back(disL);
		}

		float w_total = 0;
		for (int j = 0; j < angle.size(); j++)
		{
			weights[j] = (tan(angle[(j - 1 + angle.size()) % angle.size()] / 2) + tan(angle[j] / 2)) / weights[j];
			w_total += weights[j];
		}
		for (int j = 0; j < weights.size(); j++)
		{
			ls.lamda_ks.push_back(weights[j] / w_total);
		}
		LamdasforVS.push_back(ls);
	}
	//solve laplacian equation
	int Inner_N = vts_patch_Inner.size();
	MatrixXf A(Inner_N, Inner_N);
	VectorXf X(Inner_N), B(Inner_N), Y(Inner_N), BB(Inner_N), Z(Inner_N), BBB(Inner_N);

	for (int i = 0; i < Inner_N; i++)
	{
		for (int j = 0; j < Inner_N; j++)
		{
			A(i, j) = 0;
		}
		B(i) = 0;
		BB(i) = 0;
		BBB(i) = 0;
		for (int j = 0; j < vts_patch[i].neighborv.size(); j++)
		{
			if (vts_patch[i].neighborv[j] >= Inner_N)
			{
				B(i) += LamdasforVS[i].lamda_ks[j] * pc.UVW_coords[VTs_arrayIds_patch_reverse[vts_patch[i].neighborv[j]]][0];
				BB(i) += LamdasforVS[i].lamda_ks[j] * pc.UVW_coords[VTs_arrayIds_patch_reverse[vts_patch[i].neighborv[j]]][1];
				BBB(i) += LamdasforVS[i].lamda_ks[j] * pc.UVW_coords[VTs_arrayIds_patch_reverse[vts_patch[i].neighborv[j]]][2];
			}
			else
				A(i, vts_patch[i].neighborv[j]) = -LamdasforVS[i].lamda_ks[j];
		}
		A(i, i) = 1;
	}

	//sparse matrix and solver
	std::vector<TT> coefficients;
	for (int i = 0; i < Inner_N; i++)
	{
		for (int j = 0; j < Inner_N; j++)
		{
			if (A(i, j) != 0)
			{
				coefficients.push_back(TT(i, j, A(i, j)));
			}
		}
	}
	SpMat LeftA(Inner_N, Inner_N);
	LeftA.setFromTriplets(coefficients.begin(), coefficients.end());

	Eigen::SparseLU<SpMat> chol(LeftA);
	X = chol.solve(B);
	Y = chol.solve(BB);
	Z = chol.solve(BBB);

	// 	cout<<"A= "<<A<<endl;
	// 	cout<<"B= "<<B<<endl;
	// 	cout<<"BB= "<<BB<<endl;
	// 	cout<<"BBB= "<<BBB<<endl;
	//
	// 	X=A.colPivHouseholderQr().solve(B);
	// 	Y=A.colPivHouseholderQr().solve(BB);
	// 	Z=A.colPivHouseholderQr().solve(BBB);

	for (int i = 0; i < Inner_N; i++)
	{
		vector<float> uvw;
		uvw.push_back(X(i));
		uvw.push_back(Y(i));
		uvw.push_back(Z(i));
		pc.UVW_coords[VTs_arrayIds_patch_reverse[i]] = uvw;
		printf("uvw: %f,%f,%f\n", uvw[0], uvw[1], uvw[2]);
	}

	return true;
}
bool parameterization::D3_parameterization_tetgen(int C_Id, Parameterization_Cell &pc, vector<Vertex> &vts, vector<Edge> &ets, vector<Triangle> &fts,
		vector<Tet> &tets)
{
	//find all 3D 6 surfaces
	vector<int> VTs_arrayIds, VTs_arrayIds_reverse;
	vector<Vertex> vts_Inner;

	for (int i = 0; i < vts.size(); i++)
		VTs_arrayIds.push_back(-1);
	for (int i = 0; i < vts.size(); i++)
	{
		if (vts[i].is_boundary == 1)
			VTs_arrayIds[i] = -2;
	}
	for (int i = 0; i < VTs_arrayIds.size(); i++)
	{
		if (VTs_arrayIds[i] == -1)
		{
			VTs_arrayIds[i] = vts_Inner.size();
			vts_Inner.push_back(vts[i]);
			VTs_arrayIds_reverse.push_back(i);
		}
	}
	if (!vts_Inner.size())
		return false;

	//cal_lamdas
	vector<Lamda> LamdasforVS;
	for (int i = 0; i < vts_Inner.size(); i++)
	{
		Lamda ls;
		for (int j = 0; j < vts_Inner[i].neighborv.size(); j++)
			ls.lamda_ks.push_back(0.0);

		std::vector<float> u_lens;
		std::vector<std::vector<float>> uss;
		bool almost_vj = false;

		for (int j = 0; j < vts_Inner[i].neighborv.size(); j++)
		{
			float dis;
			DISTANCE(dis, vts[vts_Inner[i].neighborv[j]].v, vts_Inner[i].v);
			u_lens.push_back(dis);
			if (dis < EPS)
			{
				ls.lamda_ks[j] = 1.0;
				almost_vj = true;
				break;
			}

			vector<float> minuses;
			minuses.push_back((vts[vts_Inner[i].neighborv[j]].v[0] - vts_Inner[i].v[0]) / dis);
			minuses.push_back((vts[vts_Inner[i].neighborv[j]].v[1] - vts_Inner[i].v[1]) / dis);
			minuses.push_back((vts[vts_Inner[i].neighborv[j]].v[2] - vts_Inner[i].v[2]) / dis);
			uss.push_back(minuses);
		}
		if (almost_vj)
		{
			LamdasforVS.push_back(ls);
			continue;
		}
		bool inside_t = false;
		float Wj_total = 0.0;
		for (int j = 0; j < vts_Inner[i].neighbor_ot.size(); j++)
		{
			int tid = vts_Inner[i].neighbor_ot[j];

			vector<vector<float>> uis;
			vector<float> lis, thetais, cis, sis, diss, wis;
			float h = 0.0;
			for (int k = 0; k < 3; k++)
			{
				int vid = fts[tid].triangle_v[k];
				int whichone = set_contain(vts_Inner[i].neighborv, vid);
				diss.push_back(u_lens[whichone]);
				uis.push_back(uss[whichone]);
			}
			for (int k = 0; k < 3; k++)
			{
				float dis;
				DISTANCE(dis, uis[(k + 1) % 3], uis[(k - 1 + 3) % 3]);
				lis.push_back(dis);
				thetais.push_back(2 * asin(lis[k] / 2));
				h += thetais[k];
			}
			h /= 2;
			if (abs(PI - h) < EPS)
			{
				vector<float> ws;
				for (int k = 0; k < 3; k++)
				{
					ws.push_back(sin(thetais[k]) * diss[(k + 1) % 3] * diss[(k + 2) % 3]);
					Wj_total += ws[k];
				}
				for (int k = 0; k < 3; k++)
				{
					int vid = fts[tid].triangle_v[k];
					int whichone = set_contain(vts_Inner[i].neighborv, vid);
					ls.lamda_ks[whichone] += ws[k];
				}
				inside_t = true;
				break;
			}
			bool outsidet = false;
			for (int k = 0; k < 3; k++)
			{
				cis.push_back(2 * sin(h) * sin(h - thetais[k]) / (sin(thetais[(k + 1) % 3]) * sin(thetais[(k - 1 + 3) % 3])) - 1);
				float deltau = uis[0][0] * (uis[1][1] * uis[2][2] - uis[2][1] * uis[1][2])
						- uis[1][0] * (uis[0][1] * uis[2][2] - uis[2][1] * uis[0][2]) + uis[2][0] * (uis[0][1] * uis[1][2] - uis[1][1] * uis[0][2]);
				//float deltau=uis[0][0]*(uis[1][1]*uis[2][2]-uis[2][1]*uis[1][2])-uis[0][1]*(uis[1][0]*uis[2][2]-uis[1][2]*uis[2][0])+uis[0][2]*(uis[1][0]*uis[2][1]-uis[1][1]*uis[2][0]);
				float minuscisk = 1 - cis[k] * cis[k];
				if (minuscisk < 0)
					minuscisk = 0;
				float si = sqrt(minuscisk);

				if (deltau < 0)
					si = -si;
				sis.push_back(si);
				if (abs(sis[k]) < EPS)
					outsidet = true;
			}

			if (outsidet)
				continue;

			for (int k = 0; k < 3; k++)
			{
				wis.push_back(
						(thetais[k] - cis[(k + 1) % 3] * thetais[(k - 1 + 3) % 3] - cis[(k - 1 + 3) % 3] * thetais[(k + 1) % 3])
								/ (diss[k] * sin(thetais[(k + 1) % 3]) * sis[(k - 1 + 3) % 3]));

				int vid = fts[tid].triangle_v[k];
				int whichone = set_contain(vts_Inner[i].neighborv, vid);

				ls.lamda_ks[whichone] += abs(wis[k]);
				Wj_total += abs(wis[k]);
			}
		}

		for (int j = 0; j < vts_Inner[i].neighborv.size(); j++)
		{
			ls.lamda_ks[j] /= Wj_total;
		}
		LamdasforVS.push_back(ls);
	}
	//solve laplacian equation
	int Inner_N = vts_Inner.size();
	//MatrixXf A(Inner_N,Inner_N);
	VectorXf X(Inner_N), B(Inner_N), Y(Inner_N), BB(Inner_N), Z(Inner_N), BBB(Inner_N);

	std::vector<TT> coefficients;

	for (int i = 0; i < Inner_N; i++)
	{
		// 		for(int j=0;j<Inner_N;j++)
		// 		{
		// 			A(i,j)=0;
		// 		}
		B(i) = 0;
		BB(i) = 0;
		BBB(i) = 0;
		for (int j = 0; j < vts_Inner[i].neighborv.size(); j++)
		{
			if (VTs_arrayIds[vts_Inner[i].neighborv[j]] == -2)
			{
				B(i) += LamdasforVS[i].lamda_ks[j] * pc.UVW_coords[vts_Inner[i].neighborv[j]][0];
				BB(i) += LamdasforVS[i].lamda_ks[j] * pc.UVW_coords[vts_Inner[i].neighborv[j]][1];
				BBB(i) += LamdasforVS[i].lamda_ks[j] * pc.UVW_coords[vts_Inner[i].neighborv[j]][2];
			}
			else
				//A(i,VTs_arrayIds[vts_Inner[i].neighborv[j]])=-LamdasforVS[i].lamda_ks[j];
				coefficients.push_back(TT(i, VTs_arrayIds[vts_Inner[i].neighborv[j]], -LamdasforVS[i].lamda_ks[j]));
		}
		coefficients.push_back(TT(i, i, 1.0));

	}
	//sparse matrix and solver
	SpMat LeftA(Inner_N, Inner_N);
	LeftA.setFromTriplets(coefficients.begin(), coefficients.end());

	Eigen::SparseLU<SpMat> chol(LeftA);
	X = chol.solve(B);
	Y = chol.solve(BB);
	Z = chol.solve(BBB);

	// 	{
	// 		// 	cout<<"A= "<<A<<endl;
	// 		// 	cout<<"B= "<<B<<endl;
	// 		// 	cout<<"BB= "<<BB<<endl;
	// 		// 	cout<<"BBB= "<<BBB<<endl;
	// 		printf("solving Linear system for X coordinates %d\n",Inner_N);
	// 		X=A.colPivHouseholderQr().solve(B);
	// 		printf("solving Linear system for Y coordinates %d\n",Inner_N);
	// 		Y=A.colPivHouseholderQr().solve(BB);
	// 		printf("solving Linear system for Z coordinates %d\n",Inner_N);
	// 		Z=A.colPivHouseholderQr().solve(BBB);
	// 	}

	for (int i = 0; i < VTs_arrayIds_reverse.size(); i++)
	{
		vector<float> uvw;
		uvw.push_back(X(i));
		uvw.push_back(Y(i));
		uvw.push_back(Z(i));
		pc.UVW_coords[VTs_arrayIds_reverse[i]] = uvw;
		//printf("uvw: %f,%f,%f\n",uvw[0],uvw[1],uvw[2]);
	}
	return true;
}

//a hex to 24 tets
void parameterization::parameterization_cuboid(int C_Id, Parameterization_Cell &pc)
{
	h_io io;

	vector<Hex> hs;
	for (int m = 0; m < FrameI.FHs[C_Id].hs_net.size(); m++)
	{
		int hid = FrameI.FHs[C_Id].hs_net[m];
		Hex h = hex_mesh.HHs[hid];
		hs.push_back(h);
	}
	char fname1[300];
	sprintf(fname1, "%s%d%s", "C:/Users/Xifeng_Gao/Desktop/temps/x", C_Id, "_before.mesh");
	//io.write_hex_mesh_mesh(hex_mesh.HVs,hs,fname1);

	vector<int> Vs_arrayIds;
	vector<Hex_V> Vs_patch;
	vector<Hex_E> Es_patch;
	vector<Hex_F> Fs_patch;
	vector<Hex> Hs_patch;
	reindexing_cuboid_mesh(C_Id, Vs_arrayIds, Vs_patch, Es_patch, Fs_patch, Hs_patch);

	vector<int> VTs_arrayIds;
	vector<Vertex> VTs_patch;
	vector<Edge> ETs_patch;
	vector<Triangle> FTs_patch;
	vector<Tet> TETs_patch;
	translate_to_tetmesh(Vs_arrayIds, Vs_patch, Es_patch, Fs_patch, Hs_patch, VTs_arrayIds, VTs_patch, ETs_patch, FTs_patch, TETs_patch);

	//h_io io;
	char fname[300];
	sprintf(fname, "%s%d%s", "C:/Users/Xifeng_Gao/Desktop/temps/cuboid_", C_Id, ".obj");
	//io.write_tet_mesh_obj(VTs_patch,TETs_patch,"C:/Users/Xifeng_Gao/Desktop/temps/cuboid.obj");
	//io.write_triangle_mesh_obj(VTs_patch,FTs_patch,fname);

	for (int i = 0; i < VTs_patch.size(); i++)
	{
		vector<float> uvw;
		pc.UVW_coords.push_back(uvw);
	}

	//point
	//printf("8 Corners UVW\n");
	for (int i = 0; i < 8; i++)
	{
		vector<float> uvw;
		uvw.push_back(corners[i][0]);
		uvw.push_back(corners[i][1]);
		uvw.push_back(corners[i][2]);
		pc.UVW_coords[VTs_arrayIds[FrameI.FVs[FrameI.FHs[C_Id].FV_Ids[i]].index_hex]] = uvw;

		//printf("uvw: %f,%f,%f\n",uvw[0],uvw[1],uvw[2]);
	}
	//curve
	//printf("12 Edges UVW\n");
	for (int i = 0; i < 12; i++)
	{
		int fe = FrameI.FHs[C_Id].neighbor_ES[i];
		float total_len = 0;
		vector<float> ratios;
		for (int j = 0; j < FrameI.FEs[fe].vs_link.size() - 1; j++)
		{
			int v1, v2;
			v1 = FrameI.FEs[fe].vs_link[j];
			v2 = FrameI.FEs[fe].vs_link[j + 1];
			float dis;
			DISTANCE(dis, hex_mesh.HVs[v1].v, hex_mesh.HVs[v2].v);
			total_len += dis;
			ratios.push_back(total_len);
		}
		int sv_id = FrameI.FEs[fe].vs_link[0], ev_id = FrameI.FEs[fe].vs_link[FrameI.FEs[fe].vs_link.size() - 1];
		vector<float> startuvw = pc.UVW_coords[VTs_arrayIds[sv_id]], enduvw = pc.UVW_coords[VTs_arrayIds[ev_id]];

		for (int j = 1; j < FrameI.FEs[fe].vs_link.size() - 1; j++)
		{
			vector<float> uvw;
			uvw.push_back(startuvw[0] + ratios[j - 1] / total_len * (enduvw[0] - startuvw[0]));
			uvw.push_back(startuvw[1] + ratios[j - 1] / total_len * (enduvw[1] - startuvw[1]));
			uvw.push_back(startuvw[2] + ratios[j - 1] / total_len * (enduvw[2] - startuvw[2]));
			pc.UVW_coords[VTs_arrayIds[FrameI.FEs[fe].vs_link[j]]] = uvw;

			//printf("uvw: %f,%f,%f\n",uvw[0],uvw[1],uvw[2]);
		}
	}
	//face
	for (int i = 0; i < 6; i++)
	{
		//printf("Face %d UVW\n",i);
		D2_parameterization(FrameI.FHs[C_Id].neighbor_FS[i], pc, VTs_arrayIds, VTs_patch, ETs_patch, FTs_patch);
	}
	//volume
	//printf("Component %d UVW\n",C_Id);

	clock_t start_time = clock();
	//watch->startTimer();
	D3_parameterization(C_Id, pc, VTs_patch, ETs_patch, FTs_patch, TETs_patch);
	//watch->stopTimer();
	clock_t end_time = clock();
	//	std::cout<< "Running time is: "<<static_cast<double>(end_time-start_time)/CLOCKS_PER_SEC*1000<<"ms"<<std::endl;

	//printf("Parameteric of component %d UVW\n",C_Id);
	parametric_coords(C_Id, pc, VTs_patch, TETs_patch);
}
void parameterization::reindexing_cuboid_mesh(int C_Id, vector<int> &Vs_arrayIds, vector<Hex_V> &vs, vector<Hex_E> &es, vector<Hex_F> &fs,
		vector<Hex> &hs)
{
	for (int i = 0; i < hex_mesh.HVs.size(); i++)
		Vs_arrayIds.push_back(-2);

	for (int j = 0; j < FrameI.FHs[C_Id].hs_net.size(); j++)
	{
		int hid = FrameI.FHs[C_Id].hs_net[j];
		for (int k = 0; k < 8; k++)
			Vs_arrayIds[hex_mesh.HHs[hid].V_Ids[k]] = -1;

		Hex h;
		for (int k = 0; k < 8; k++)
			h.V_Ids[k] = hex_mesh.HHs[hid].V_Ids[k];
		h.index = hs.size();
		hs.push_back(h);
	}

//	printf("hs collection finished!.\n");

	for (int i = 0; i < Vs_arrayIds.size(); i++)
	{
		if (Vs_arrayIds[i] == -1)
		{
			Hex_V v;
			v.v[0] = hex_mesh.HVs[i].v[0];
			v.v[1] = hex_mesh.HVs[i].v[1];
			v.v[2] = hex_mesh.HVs[i].v[2];
			v.index = vs.size();
			v.fixed = false;
			v.Frame_V_id = -1;
			v.where_location = -1;
			v.is_on_base_complex = 0;
			v.which_F_face = -1;
			vs.push_back(v);
			Vs_arrayIds[i] = v.index;
		}
	}
//	printf("Vs_arrayIds collection finished!.\n");
	for (int i = 0; i < hs.size(); i++)
	{
		for (int j = 0; j < 8; j++)
		{
			hs[i].V_Ids[j] = Vs_arrayIds[hs[i].V_Ids[j]];
			vs[hs[i].V_Ids[j]].neighbor_Hs.push_back(i);
		}
	}
//	printf("hs vs replacing finished!.\n");
	construct_Es(vs, es, hs);
//	printf("Es construction finished!.\n");
	construct_Fs(vs, es, fs, hs);
//	printf("Fs construction finished!.\n");
	determine_boundary_info(vs, es, fs, hs);
//	printf("boundary info finished!.\n");
	calculation_centroid(vs, fs, hs);
//	printf("centroid finished!.\n");

	for (int i = 0; i < 8; i++)
		vs[Vs_arrayIds[FrameI.FVs[FrameI.FHs[C_Id].FV_Ids[i]].index_hex]].Frame_V_id = FrameI.FVs[FrameI.FHs[C_Id].FV_Ids[i]].index_own;

	for (int i = 0; i < es.size(); i++)
		es[i].which_F_edge = -1;
	for (int i = 0; i < 12; i++)
	{
		for (int j = 0; j < FrameI.FEs[FrameI.FHs[C_Id].neighbor_ES[i]].vs_link.size() - 1; j++)
		{
			int v1, v2;
			v1 = FrameI.FEs[FrameI.FHs[C_Id].neighbor_ES[i]].vs_link[j];
			v1 = Vs_arrayIds[v1];
			v2 = FrameI.FEs[FrameI.FHs[C_Id].neighbor_ES[i]].vs_link[j + 1];
			v2 = Vs_arrayIds[v2];

			vs[v1].is_on_base_complex = 1;
			vs[v2].is_on_base_complex = 1;
		}
	}

	for (int i = 0; i < fs.size(); i++)
		fs[i].which_F_face = -1;
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < FrameI.FFs[FrameI.FHs[C_Id].neighbor_FS[i]].hfs_net_another.size(); j++)
		{
			int fid = FrameI.FFs[FrameI.FHs[C_Id].neighbor_FS[i]].hfs_net_another[j];
			int v1, v2, v3, v4;
			v1 = Vs_arrayIds[hex_mesh.HFs[fid].cv_Ids[0]];
			v2 = Vs_arrayIds[hex_mesh.HFs[fid].cv_Ids[1]];
			v3 = Vs_arrayIds[hex_mesh.HFs[fid].cv_Ids[2]];
			v4 = Vs_arrayIds[hex_mesh.HFs[fid].cv_Ids[3]];

			vector<int> sharedf12;
			set_cross(vs[v1].neighbor_Fs, vs[v2].neighbor_Fs, sharedf12);
			vector<int> sharedf34;
			set_cross(vs[v3].neighbor_Fs, vs[v4].neighbor_Fs, sharedf34);
			vector<int> sharedf;
			set_cross(sharedf12, sharedf34, sharedf);
			fs[sharedf[0]].which_F_face = FrameI.FHs[C_Id].neighbor_FS[i];

			vs[v1].which_F_face = FrameI.FHs[C_Id].neighbor_FS[i];
			vs[v2].which_F_face = FrameI.FHs[C_Id].neighbor_FS[i];
			vs[v3].which_F_face = FrameI.FHs[C_Id].neighbor_FS[i];
			vs[v4].which_F_face = FrameI.FHs[C_Id].neighbor_FS[i];
		}
	}
}
void parameterization::translate_to_tetmesh(vector<int> &Vs_arrayIds, vector<Hex_V> &vs, vector<Hex_E> &es, vector<Hex_F> &fs, vector<Hex> &hs,
		vector<int> &VTs_arrayIds, vector<Vertex> &vts, vector<Edge> &ets, vector<Triangle> &fts, vector<Tet> &tets)
{
	for (int i = 0; i < vs.size(); i++)
	{
		Vertex v;
		v.index = i;
		v.v[0] = vs[i].v[0];
		v.v[1] = vs[i].v[1];
		v.v[2] = vs[i].v[2];
		v.on_base_complex = vs[i].is_on_base_complex;
		v.which_F_face = vs[i].which_F_face;
		vts.push_back(v);
	}
	vector<Vertex> FVs, HVs;
	//face v
	for (int i = 0; i < fs.size(); i++)
	{
		Vertex v;
		v.index = vts.size();
		v.v[0] = fs[i].centroid[0];
		v.v[1] = fs[i].centroid[1];
		v.v[2] = fs[i].centroid[2];
		v.on_base_complex = false;
		v.which_F_face = fs[i].which_F_face;
		FVs.push_back(v);
		vts.push_back(v);
	}
	//hex v
	for (int i = 0; i < hs.size(); i++)
	{
		Vertex v;
		v.index = vts.size();
		v.v[0] = hs[i].centroid[0];
		v.v[1] = hs[i].centroid[1];
		v.v[2] = hs[i].centroid[2];
		v.on_base_complex = false;
		v.which_F_face = -1;
		HVs.push_back(v);
		vts.push_back(v);
	}
	//tets for a hex
	for (int i = 0; i < hs.size(); i++)
	{
		for (int j = 0; j < 6; j++)
		{
			int fid = hs[i].F_Ids[j];
			for (int k = 0; k < 4; k++)
			{
				int v1 = es[fs[fid].ce_Ids[k]].startend_Id[0], v2 = es[fs[fid].ce_Ids[k]].startend_Id[1];
				Tet t;
				t.vs.push_back(HVs[i].index);
				t.vs.push_back(v1);
				t.vs.push_back(v2);
				t.vs.push_back(FVs[fid].index);
				tets.push_back(t);
			}
		}
	}
	//triangles for a tet
	for (int i = 0; i < tets.size(); i++)
	{
		Triangle t;
		for (int j = 0; j < 4; j++)
		{
			t.triangle_v[0] = tets[i].vs[(j + 1) % 4];
			t.triangle_v[1] = tets[i].vs[(j + 2) % 4];
			t.triangle_v[2] = tets[i].vs[(j + 3) % 4];
			t.index = fts.size();

			vector<int> result1, result2;
			set_cross(vts[t.triangle_v[0]].neighbort, vts[t.triangle_v[1]].neighbort, result1);

			set_cross(result1, vts[t.triangle_v[2]].neighbort, result2);
			if (!result2.size())
			{
				t.which_F_face = -1;
				if (vts[t.triangle_v[0]].which_F_face != -1)
				{
					int which_f = -1;
					vector<int> non_baseedgevs;
					for (int k = 0; k < 3; k++)
					{
						if (!vts[t.triangle_v[k]].on_base_complex)
						{
							which_f = vts[t.triangle_v[k]].which_F_face;
							non_baseedgevs.push_back(t.triangle_v[k]);
						}
					}
					bool allequal = true;
					for (int k = 0; k < non_baseedgevs.size(); k++)
					{
						if (vts[non_baseedgevs[k]].which_F_face != which_f)
							allequal = false;
					}
					if (which_f != -1 && allequal)
						t.which_F_face = which_f;
				}

				fts.push_back(t);
				vts[t.triangle_v[0]].neighbort.push_back(t.index);
				vts[t.triangle_v[1]].neighbort.push_back(t.index);
				vts[t.triangle_v[2]].neighbort.push_back(t.index);
				tets[i].ts.push_back(t.index);
				vts[tets[i].vs[j]].neighbor_ot.push_back(t.index);
			}
			else
			{
				tets[i].ts.push_back(result2[0]);
				vts[tets[i].vs[j]].neighbor_ot.push_back(result2[0]);
			}
			vts[tets[i].vs[j]].neighborv.push_back(t.triangle_v[0]);
			vts[tets[i].vs[j]].neighborv.push_back(t.triangle_v[1]);
			vts[tets[i].vs[j]].neighborv.push_back(t.triangle_v[2]);

			vts[tets[i].vs[j]].neighbortet.push_back(i);
		}
	}
	for (int i = 0; i < vts.size(); i++)
	{
		set_redundent_clearn(vts[i].neighborv);
	}
	//assign which face, which edge.
	VTs_arrayIds = Vs_arrayIds;
}
bool parameterization::D2_parameterization(int F_Id, Parameterization_Cell &pc, vector<int> &VTs_arrayIds, vector<Vertex> &vts, vector<Edge> &ets,
		vector<Triangle> &fts)
{
	//re-indexing 2D surface patch
	vector<int> VTs_arrayIds_patch, VTs_arrayIds_patch_reverse;
	vector<Vertex> vts_patch;
	vector<Vertex> vts_patch_Inner;
	vector<Vertex> vts_patch_Outer;
	vector<Triangle> fts_patch;

	for (int i = 0; i < vts.size(); i++)
		VTs_arrayIds_patch.push_back(-2);

	for (int i = 0; i < fts.size(); i++)
	{
		if (fts[i].which_F_face == F_Id)
		{
			Triangle t;
			t.triangle_v[0] = fts[i].triangle_v[0];
			t.triangle_v[1] = fts[i].triangle_v[1];
			t.triangle_v[2] = fts[i].triangle_v[2];
			t.index = fts_patch.size();
			fts_patch.push_back(t);
			for (int j = 0; j < 3; j++)
			{
				if (VTs_arrayIds_patch[fts[i].triangle_v[j]] != -2)
					continue;
				VTs_arrayIds_patch[fts[i].triangle_v[j]] = -1;
				Vertex v;
				v.v[0] = vts[fts[i].triangle_v[j]].v[0];
				v.v[1] = vts[fts[i].triangle_v[j]].v[1];
				v.v[2] = vts[fts[i].triangle_v[j]].v[2];
				v.index = fts[i].triangle_v[j];
				v.on_base_complex = vts[fts[i].triangle_v[j]].on_base_complex;

				if (v.on_base_complex)
					vts_patch_Outer.push_back(v);
				else
					vts_patch_Inner.push_back(v);

				VTs_arrayIds_patch_reverse.push_back(-1);
			}
		}
	}
	if (!vts_patch_Inner.size())
		return false;

	for (int i = 0; i < vts_patch_Inner.size(); i++)
	{
		VTs_arrayIds_patch_reverse[i] = vts_patch_Inner[i].index;
		VTs_arrayIds_patch[vts_patch_Inner[i].index] = i;
		vts_patch_Inner[i].index = i;
		vts_patch.push_back(vts_patch_Inner[i]);
	}
	for (int i = 0; i < vts_patch_Outer.size(); i++)
	{
		VTs_arrayIds_patch_reverse[vts_patch_Inner.size() + i] = vts_patch_Outer[i].index;

		VTs_arrayIds_patch[vts_patch_Outer[i].index] = vts_patch_Inner.size() + i;
		vts_patch_Outer[i].index = vts_patch_Inner.size() + i;
		vts_patch.push_back(vts_patch_Outer[i]);
	}

	for (int i = 0; i < fts_patch.size(); i++)
	{
		int v1, v2, v3;
		v1 = VTs_arrayIds_patch[fts_patch[i].triangle_v[0]];
		v2 = VTs_arrayIds_patch[fts_patch[i].triangle_v[1]];
		v3 = VTs_arrayIds_patch[fts_patch[i].triangle_v[2]];

		fts_patch[i].triangle_v[0] = v1;
		fts_patch[i].triangle_v[1] = v2;
		fts_patch[i].triangle_v[2] = v3;
		vts_patch[v1].neighbort.push_back(i);
		vts_patch[v2].neighbort.push_back(i);
		vts_patch[v3].neighbort.push_back(i);

		vts_patch[v1].neighborv.push_back(v2);
		vts_patch[v1].neighborv.push_back(v3);
		vts_patch[v2].neighborv.push_back(v1);
		vts_patch[v2].neighborv.push_back(v3);
		vts_patch[v3].neighborv.push_back(v1);
		vts_patch[v3].neighborv.push_back(v2);
	}
	for (int i = 0; i < vts_patch.size(); i++)
	{
		set_redundent_clearn(vts_patch[i].neighborv);
		//orient vi's neighborvs
		vector<int> orient_nvs;
		int bef = -1, cur = vts_patch[i].neighborv[0];
		orient_nvs.push_back(cur);
		for (int j = 1; j < vts_patch[i].neighborv.size(); j++)
		{
			for (int k = 0; k < vts_patch[i].neighborv.size(); k++)
			{
				if (set_contain(vts_patch[cur].neighborv, vts_patch[i].neighborv[k]) != -1 && vts_patch[i].neighborv[k] != bef)
				{
					orient_nvs.push_back(vts_patch[i].neighborv[k]);
					bef = cur;
					cur = vts_patch[i].neighborv[k];
					break;
				}
			}
		}
		vts_patch[i].neighborv = orient_nvs;
	}
	//cal_lamdas
	vector<Lamda> LamdasforVS;
	for (int i = 0; i < vts_patch_Inner.size(); i++)
	{
		Lamda ls;
		std::vector<float> angle;
		std::vector<float> weights;
		for (int j = 0; j < vts_patch[i].neighborv.size(); j++)
		{
			float dis;
			DISTANCE(dis, vts_patch[vts_patch[i].neighborv[j]].v, vts_patch[vts_patch[i].neighborv[(j + 1) % vts_patch[i].neighborv.size()]].v);
			float disL;
			DISTANCE(disL, vts_patch[vts_patch[i].neighborv[j]].v, vts_patch[i].v);
			float disR;
			DISTANCE(disR, vts_patch[i].v, vts_patch[vts_patch[i].neighborv[(j + 1) % vts_patch[i].neighborv.size()]].v);

			float dis_minus = (disL * disL + disR * disR - dis * dis);
			float divide = dis_minus / (2 * disL * disR);
			if (divide < -1)
				divide = -1;
			if (divide > 1)
				divide = 1;
			float angle_ = acos(divide);
			angle.push_back(angle_);
			weights.push_back(disL);
		}

		float w_total = 0;
		for (int j = 0; j < angle.size(); j++)
		{
			weights[j] = (tan(angle[(j - 1 + angle.size()) % angle.size()] / 2) + tan(angle[j] / 2)) / weights[j];
			w_total += weights[j];
		}
		for (int j = 0; j < weights.size(); j++)
		{
			ls.lamda_ks.push_back(weights[j] / w_total);
		}
		LamdasforVS.push_back(ls);
	}
	//solve laplacian equation
	int Inner_N = vts_patch_Inner.size();
	MatrixXf A(Inner_N, Inner_N);
	VectorXf X(Inner_N), B(Inner_N), Y(Inner_N), BB(Inner_N), Z(Inner_N), BBB(Inner_N);

	for (int i = 0; i < Inner_N; i++)
	{
		for (int j = 0; j < Inner_N; j++)
		{
			A(i, j) = 0;
		}
		B(i) = 0;
		BB(i) = 0;
		BBB(i) = 0;
		for (int j = 0; j < vts_patch[i].neighborv.size(); j++)
		{
			if (vts_patch[i].neighborv[j] >= Inner_N)
			{
				B(i) += LamdasforVS[i].lamda_ks[j] * pc.UVW_coords[VTs_arrayIds_patch_reverse[vts_patch[i].neighborv[j]]][0];
				BB(i) += LamdasforVS[i].lamda_ks[j] * pc.UVW_coords[VTs_arrayIds_patch_reverse[vts_patch[i].neighborv[j]]][1];
				BBB(i) += LamdasforVS[i].lamda_ks[j] * pc.UVW_coords[VTs_arrayIds_patch_reverse[vts_patch[i].neighborv[j]]][2];
			}
			else
				A(i, vts_patch[i].neighborv[j]) = -LamdasforVS[i].lamda_ks[j];
		}
		A(i, i) = 1;
	}

	//sparse matrix and solver
	std::vector<TT> coefficients;
	for (int i = 0; i < Inner_N; i++)
	{
		for (int j = 0; j < Inner_N; j++)
		{
			if (A(i, j) != 0)
			{
				coefficients.push_back(TT(i, j, A(i, j)));
			}
		}
	}
	SpMat LeftA(Inner_N, Inner_N);
	LeftA.setFromTriplets(coefficients.begin(), coefficients.end());

	Eigen::SparseLU<SpMat> chol(LeftA);
	X = chol.solve(B);
	Y = chol.solve(BB);
	Z = chol.solve(BBB);

// 	cout<<"A= "<<A<<endl;
// 	cout<<"B= "<<B<<endl;
// 	cout<<"BB= "<<BB<<endl;
// 	cout<<"BBB= "<<BBB<<endl;
//
// 	X=A.colPivHouseholderQr().solve(B);
// 	Y=A.colPivHouseholderQr().solve(BB);
// 	Z=A.colPivHouseholderQr().solve(BBB);

	for (int i = 0; i < Inner_N; i++)
	{
		vector<float> uvw;
		uvw.push_back(X(i));
		uvw.push_back(Y(i));
		uvw.push_back(Z(i));
		pc.UVW_coords[VTs_arrayIds_patch_reverse[i]] = uvw;
		//printf("uvw: %f,%f,%f\n",uvw[0],uvw[1],uvw[2]);
	}

	return true;
}
bool parameterization::D3_parameterization(int C_Id, Parameterization_Cell &pc, vector<Vertex> &vts, vector<Edge> &ets, vector<Triangle> &fts,
		vector<Tet> &tets)
{
	//find all 3D 6 surfaces
	vector<int> VTs_arrayIds, VTs_arrayIds_reverse;
	vector<Vertex> vts_Inner;

	for (int i = 0; i < vts.size(); i++)
		VTs_arrayIds.push_back(-1);
	for (int i = 0; i < fts.size(); i++)
	{
		if (fts[i].which_F_face != -1)
			for (int j = 0; j < 3; j++)
				VTs_arrayIds[fts[i].triangle_v[j]] = -2;
	}
	for (int i = 0; i < VTs_arrayIds.size(); i++)
	{
		if (VTs_arrayIds[i] == -1)
		{
			VTs_arrayIds[i] = vts_Inner.size();
			vts_Inner.push_back(vts[i]);
			VTs_arrayIds_reverse.push_back(i);
		}
	}
	if (!vts_Inner.size())
		return false;

	//cal_lamdas
	vector<Lamda> LamdasforVS;
	for (int i = 0; i < vts_Inner.size(); i++)
	{
		Lamda ls;
		for (int j = 0; j < vts_Inner[i].neighborv.size(); j++)
			ls.lamda_ks.push_back(0.0);

		std::vector<float> u_lens;
		std::vector<std::vector<float>> uss;
		bool almost_vj = false;

		for (int j = 0; j < vts_Inner[i].neighborv.size(); j++)
		{
			float dis;
			DISTANCE(dis, vts[vts_Inner[i].neighborv[j]].v, vts_Inner[i].v);
			u_lens.push_back(dis);
			if (dis < EPS)
			{
				ls.lamda_ks[j] = 1.0;
				almost_vj = true;
				break;
			}

			vector<float> minuses;
			minuses.push_back((vts[vts_Inner[i].neighborv[j]].v[0] - vts_Inner[i].v[0]) / dis);
			minuses.push_back((vts[vts_Inner[i].neighborv[j]].v[1] - vts_Inner[i].v[1]) / dis);
			minuses.push_back((vts[vts_Inner[i].neighborv[j]].v[2] - vts_Inner[i].v[2]) / dis);
			uss.push_back(minuses);
		}
		if (almost_vj)
		{
			LamdasforVS.push_back(ls);
			continue;
		}
		bool inside_t = false;
		float Wj_total = 0.0;
		for (int j = 0; j < vts_Inner[i].neighbor_ot.size(); j++)
		{
			int tid = vts_Inner[i].neighbor_ot[j];

			vector<vector<float>> uis;
			vector<float> lis, thetais, cis, sis, diss, wis;
			float h = 0.0;
			for (int k = 0; k < 3; k++)
			{
				int vid = fts[tid].triangle_v[k];
				int whichone = set_contain(vts_Inner[i].neighborv, vid);
				diss.push_back(u_lens[whichone]);
				uis.push_back(uss[whichone]);
			}
			for (int k = 0; k < 3; k++)
			{
				float dis;
				DISTANCE(dis, uis[(k + 1) % 3], uis[(k - 1 + 3) % 3]);
				lis.push_back(dis);
				thetais.push_back(2 * asin(lis[k] / 2));
				h += thetais[k];
			}
			h /= 2;
			if (abs(PI - h) < EPS)
			{
				vector<float> ws;
				for (int k = 0; k < 3; k++)
				{
					ws.push_back(sin(thetais[k]) * diss[(k + 1) % 3] * diss[(k + 2) % 3]);
					Wj_total += ws[k];
				}
				for (int k = 0; k < 3; k++)
				{
					int vid = fts[tid].triangle_v[k];
					int whichone = set_contain(vts_Inner[i].neighborv, vid);
					ls.lamda_ks[whichone] += ws[k];
				}
				inside_t = true;
				break;
			}
			bool outsidet = false;
			for (int k = 0; k < 3; k++)
			{
				cis.push_back(2 * sin(h) * sin(h - thetais[k]) / (sin(thetais[(k + 1) % 3]) * sin(thetais[(k - 1 + 3) % 3])) - 1);
				float deltau = uis[0][0] * (uis[1][1] * uis[2][2] - uis[2][1] * uis[1][2])
						- uis[1][0] * (uis[0][1] * uis[2][2] - uis[2][1] * uis[0][2]) + uis[2][0] * (uis[0][1] * uis[1][2] - uis[1][1] * uis[0][2]);
				//float deltau=uis[0][0]*(uis[1][1]*uis[2][2]-uis[2][1]*uis[1][2])-uis[0][1]*(uis[1][0]*uis[2][2]-uis[1][2]*uis[2][0])+uis[0][2]*(uis[1][0]*uis[2][1]-uis[1][1]*uis[2][0]);
				float minuscisk = 1 - cis[k] * cis[k];
				if (minuscisk < 0)
					minuscisk = 0;
				float si = sqrt(minuscisk);

				if (deltau < 0)
					si = -si;
				sis.push_back(si);
				if (abs(sis[k]) < EPS)
					outsidet = true;
			}

			if (outsidet)
				continue;

			for (int k = 0; k < 3; k++)
			{
				wis.push_back(
						(thetais[k] - cis[(k + 1) % 3] * thetais[(k - 1 + 3) % 3] - cis[(k - 1 + 3) % 3] * thetais[(k + 1) % 3])
								/ (diss[k] * sin(thetais[(k + 1) % 3]) * sis[(k - 1 + 3) % 3]));

				int vid = fts[tid].triangle_v[k];
				int whichone = set_contain(vts_Inner[i].neighborv, vid);

				ls.lamda_ks[whichone] += abs(wis[k]);
				Wj_total += abs(wis[k]);
			}
		}

		for (int j = 0; j < vts_Inner[i].neighborv.size(); j++)
		{
			ls.lamda_ks[j] /= Wj_total;
		}
		LamdasforVS.push_back(ls);
	}
	//solve laplacian equation
	int Inner_N = vts_Inner.size();
	//MatrixXf A(Inner_N,Inner_N);
	VectorXf X(Inner_N), B(Inner_N), Y(Inner_N), BB(Inner_N), Z(Inner_N), BBB(Inner_N);

	std::vector<TT> coefficients;

	for (int i = 0; i < Inner_N; i++)
	{
// 		for(int j=0;j<Inner_N;j++)
// 		{
// 			A(i,j)=0;
// 		}
		B(i) = 0;
		BB(i) = 0;
		BBB(i) = 0;
		for (int j = 0; j < vts_Inner[i].neighborv.size(); j++)
		{
			if (VTs_arrayIds[vts_Inner[i].neighborv[j]] == -2)
			{
				B(i) += LamdasforVS[i].lamda_ks[j] * pc.UVW_coords[vts_Inner[i].neighborv[j]][0];
				BB(i) += LamdasforVS[i].lamda_ks[j] * pc.UVW_coords[vts_Inner[i].neighborv[j]][1];
				BBB(i) += LamdasforVS[i].lamda_ks[j] * pc.UVW_coords[vts_Inner[i].neighborv[j]][2];
			}
			else
				//A(i,VTs_arrayIds[vts_Inner[i].neighborv[j]])=-LamdasforVS[i].lamda_ks[j];
				coefficients.push_back(TT(i, VTs_arrayIds[vts_Inner[i].neighborv[j]], -LamdasforVS[i].lamda_ks[j]));
		}
		coefficients.push_back(TT(i, i, 1.0));

	}
	//sparse matrix and solver
	SpMat LeftA(Inner_N, Inner_N);
	LeftA.setFromTriplets(coefficients.begin(), coefficients.end());

	Eigen::SparseLU<SpMat> chol(LeftA);
	X = chol.solve(B);
	Y = chol.solve(BB);
	Z = chol.solve(BBB);

// 	{
// 		// 	cout<<"A= "<<A<<endl;
// 		// 	cout<<"B= "<<B<<endl;
// 		// 	cout<<"BB= "<<BB<<endl;
// 		// 	cout<<"BBB= "<<BBB<<endl;
// 		printf("solving Linear system for X coordinates %d\n",Inner_N);
// 		X=A.colPivHouseholderQr().solve(B);
// 		printf("solving Linear system for Y coordinates %d\n",Inner_N);
// 		Y=A.colPivHouseholderQr().solve(BB);
// 		printf("solving Linear system for Z coordinates %d\n",Inner_N);
// 		Z=A.colPivHouseholderQr().solve(BBB);
// 	}

	for (int i = 0; i < VTs_arrayIds_reverse.size(); i++)
	{
		vector<float> uvw;
		uvw.push_back(X(i));
		uvw.push_back(Y(i));
		uvw.push_back(Z(i));
		pc.UVW_coords[VTs_arrayIds_reverse[i]] = uvw;
		//printf("uvw: %f,%f,%f\n",uvw[0],uvw[1],uvw[2]);
	}
	return true;
}
//a hex to 24 tets

void parameterization::parametric_coords(int C_Id, Parameterization_Cell &pc, vector<Vertex> &vts, vector<Tet> &tets)
{
	int U, V, W;
	vector<int> result1;
	set_cross(FrameI.FVs[FrameI.FHs[C_Id].FV_Ids[0]].neighbor_Es, FrameI.FVs[FrameI.FHs[C_Id].FV_Ids[4]].neighbor_Es, result1);
	vector<int> result2;
	set_cross(FrameI.FVs[FrameI.FHs[C_Id].FV_Ids[1]].neighbor_Es, FrameI.FVs[FrameI.FHs[C_Id].FV_Ids[2]].neighbor_Es, result2);
	vector<int> result3;
	set_cross(FrameI.FVs[FrameI.FHs[C_Id].FV_Ids[0]].neighbor_Es, FrameI.FVs[FrameI.FHs[C_Id].FV_Ids[1]].neighbor_Es, result3);

	vector<int> cedges;
	for (int i = 0; i < 12; i++)
		cedges.push_back(FrameI.FHs[C_Id].neighbor_ES[i]);
	if (result1.size() > 1)
	{
		for (int i = 0; i < result1.size(); i++)
		{
			if (set_contain(cedges, result1[i]) != -1)
			{
				int tempid = result1[i];
				result1.clear();
				result1.push_back(tempid);
			}
		}
	}
	if (result2.size() > 1)
	{
		for (int i = 0; i < result2.size(); i++)
		{
			if (set_contain(cedges, result2[i]) != -1)
			{
				int tempid = result2[i];
				result2.clear();
				result2.push_back(tempid);
			}
		}
	}
	if (result3.size() > 1)
	{
		for (int i = 0; i < result3.size(); i++)
		{
			if (set_contain(cedges, result3[i]) != -1)
			{
				int tempid = result3[i];
				result3.clear();
				result3.push_back(tempid);
			}
		}
	}
	U = L_para_Ns[result1[0]];
	V = L_para_Ns[result2[0]];
	W = L_para_Ns[result3[0]];

	Us[C_Id] = U;
	Vs[C_Id] = V;
	Ws[C_Id] = W;

	int ***para_indicator;
	para_indicator = new int **[U + 1];
	for (int i = 0; i < U + 1; i++)
	{
		para_indicator[i] = new int *[V + 1];
		vector<vector<Hex_V>> hhvs;
		for (int j = 0; j < V + 1; j++)
		{
			para_indicator[i][j] = new int[W + 1];
			vector<Hex_V> hvs;
			for (int k = 0; k < W + 1; k++)
			{
				Hex_V hv;
				hv.index = -1;
				hvs.push_back(hv);
				para_indicator[i][j][k] = false;
			}
			hhvs.push_back(hvs);
		}
		pc.Parametric_coords.push_back(hhvs);
	}

	float u_step = 1.0 / U, v_step = 1.0 / V, w_step = 1.0 / W;
	vector<float> us, vs, ws;
	for (int i = 0; i < U + 1; i++)
		us.push_back(i * u_step);
	for (int i = 0; i < V + 1; i++)
		vs.push_back(i * v_step);
	for (int i = 0; i < W + 1; i++)
		ws.push_back(i * w_step);

	for (int i = 0; i < tets.size(); i++)
	{
		vector<vector<float>> vers;
		vers.push_back(pc.UVW_coords[tets[i].vs[0]]);
		vers.push_back(pc.UVW_coords[tets[i].vs[1]]);
		vers.push_back(pc.UVW_coords[tets[i].vs[2]]);
		vers.push_back(pc.UVW_coords[tets[i].vs[3]]);

		//min max
		float min_u = 1, max_u = 0, min_v = 1, max_v = 0, min_w = 1, max_w = 0;
		for (int j = 0; j < 4; j++)
		{
			if (vers[j][0] > max_u)
				max_u = vers[j][0];
			if (vers[j][0] < min_u)
				min_u = vers[j][0];

			if (vers[j][1] > max_v)
				max_v = vers[j][1];
			if (vers[j][1] < min_v)
				min_v = vers[j][1];

			if (vers[j][2] > max_w)
				max_w = vers[j][2];
			if (vers[j][2] < min_w)
				min_w = vers[j][2];
			vers[j].push_back(1.0);
		}
		vector<int> u_ids, v_ids, w_ids;
		for (int j = 0; j < us.size(); j++)
		{
			if (us[j] >= min_u - 100 * EPS && us[j] <= max_u + 100 * EPS)
				u_ids.push_back(j);
		}
		for (int j = 0; j < vs.size(); j++)
		{
			if (vs[j] >= min_v - 100 * EPS && vs[j] <= max_v + 100 * EPS)
				v_ids.push_back(j);
		}
		for (int j = 0; j < ws.size(); j++)
		{
			if (ws[j] >= min_w - 100 * EPS && ws[j] <= max_w + 100 * EPS)
				w_ids.push_back(j);
		}

		if (!u_ids.size() || !v_ids.size() || !w_ids.size())
			continue;

		bool tet_coplanar = false;
		float deter0;
		vector<float> ver1 = vers[0], ver2 = vers[1], ver3 = vers[2], ver4 = vers[3];
		deter0 = matrix_determinant(ver1, ver2, ver3, ver4);
		if (abs(deter0) < EPS)
		{
			printf("Tet %d is a coplanar-four point set\n", i);
			tet_coplanar = true;
			continue;;	//tet is a coplanar-four point set
		}
		for (int j = 0; j < u_ids.size(); j++)
		{
			for (int k = 0; k < v_ids.size(); k++)
			{
				for (int m = 0; m < w_ids.size(); m++)
				{
					if (para_indicator[u_ids[j]][v_ids[k]][w_ids[m]])
						continue;

					vector<float> lamdas;
					if (!tet_coplanar)
					{
						vector<float> deters;
						bool next = false;

						vector<float> cur_uvw;
						cur_uvw.push_back(us[u_ids[j]]);
						cur_uvw.push_back(vs[v_ids[k]]);
						cur_uvw.push_back(ws[w_ids[m]]);
						cur_uvw.push_back(1.0);
						for (int n = 0; n < 4; n++)
						{
							switch (n)
							{
							case 0:
								deters.push_back(matrix_determinant(cur_uvw, vers[1], vers[2], vers[3]));
								break;
							case 1:
								deters.push_back(matrix_determinant(vers[0], cur_uvw, vers[2], vers[3]));
								break;
							case 2:
								deters.push_back(matrix_determinant(vers[0], vers[1], cur_uvw, vers[3]));
								break;
							case 3:
								deters.push_back(matrix_determinant(vers[0], vers[1], vers[2], cur_uvw));
								break;
							}
							if (sign_of_value(deter0) * sign_of_value(deters[n]) < 0)
							{
								next = true;
								break;
							}
						}
						if (next)
							continue;
						for (int n = 0; n < 4; n++)
							lamdas.push_back(deters[n] / deter0);
					}

					Hex_V hv;
					hv.v[0] = hv.v[1] = hv.v[2] = 0;
					for (int n = 0; n < 4; n++)
					{
						hv.v[0] += lamdas[n] * vts[tets[i].vs[n]].v[0];
						hv.v[1] += lamdas[n] * vts[tets[i].vs[n]].v[1];
						hv.v[2] += lamdas[n] * vts[tets[i].vs[n]].v[2];
					}
					vector<float> v_temp;
					v_temp.push_back(-0.06078311056 * Escalar);
					v_temp.push_back(-0.06332313269 * Escalar);
					v_temp.push_back(0.7998717427 * Escalar);
					float dis;
					DISTANCE(dis, hv.v, v_temp);
					if (abs(dis) < 100 * EPS)
					{
						printf("x %f y %f z %f; lamda %f %f %f %f; deter0 %f\n", hv.v[0], hv.v[1], hv.v[2], lamdas[0], lamdas[1], lamdas[2],
								lamdas[3], deter0);

					}
					hv.index = -1;
					pc.Parametric_coords[u_ids[j]][v_ids[k]][w_ids[m]] = hv;
					para_indicator[u_ids[j]][v_ids[k]][w_ids[m]] = true;
				}
			}
		}
	}

	h_io io;
	vector<Hex_V> tempvers;
	for (int i = 0; i < U + 1; i++)
	{
		for (int j = 0; j < V + 1; j++)
		{
			for (int k = 0; k < W + 1; k++)
			{
				//printf("para_indicator: i %d, j %d, k %d: %d; x %f,y %f,z %f\n",i,j,k,para_indicator[i][j][k],pc.Parametric_coords[i][j][k].v[0],pc.Parametric_coords[i][j][k].v[1],pc.Parametric_coords[i][j][k].v[2]);
				if (!para_indicator[i][j][k])
				{
					pc.Parametric_coords[i][j][k].v[0] = 0;
					pc.Parametric_coords[i][j][k].v[1] = 0;
					pc.Parametric_coords[i][j][k].v[2] = 0;
					printf("para_indicator: i %d, j %d, k %d: u %f,v %f,w %f\n", i, j, k, float(i) / U, float(j) / V, float(k) / W);
				}
				Hex_V v = pc.Parametric_coords[i][j][k];
				para_indicator[i][j][k] = i * (V + 1) * (W + 1) + j * (W + 1) + k;
				tempvers.push_back(v);
			}
		}
	}
	vector<Hex> hs;
	for (int j = 0; j < U; j++)
	{
		for (int k = 0; k < V; k++)
		{
			for (int m = 0; m < W; m++)
			{
				Hex h;
				h.V_Ids[0] = para_indicator[j][k][m];
				h.V_Ids[1] = para_indicator[j][k + 1][m];
				h.V_Ids[2] = para_indicator[j][k + 1][m + 1];
				h.V_Ids[3] = para_indicator[j][k][m + 1];

				h.V_Ids[4] = para_indicator[j + 1][k][m];
				h.V_Ids[5] = para_indicator[j + 1][k + 1][m];
				h.V_Ids[6] = para_indicator[j + 1][k + 1][m + 1];
				h.V_Ids[7] = para_indicator[j + 1][k][m + 1];

				hs.push_back(h);
			}
		}
	}
// 	for(int i=0;i<vts.size();i++)
// 	{
// 		Vertex v;
// 		v.v[0]=pc.UVW_coords[i][0];
// 		v.v[1]=pc.UVW_coords[i][1];
// 		v.v[2]=pc.UVW_coords[i][2];
// 		tempvers.push_back(v);
// 	}
	char fname[300];
	sprintf(fname, "%s%d%s", "C:/Users/Xifeng_Gao/Desktop/temps/x", C_Id, ".mesh");
	//io.write_hex_mesh_mesh(tempvers,hs,fname);

}
void parameterization::produce_hex_mesh(vector<Parameterization_Cell> &PCs, char *path)
{
	vector<Hex_V> hvs;
	unsigned int V_Ind = 0, H_Ind = 0;
	vector<Hex> hhs;
	for (int i = 0; i < FrameI.FVs.size(); i++)
	{
		for (int j = 0; j < FrameI.FVs[i].neighbor_Hs.size(); j++)
		{
			int nh = FrameI.FVs[i].neighbor_Hs[j];
			for (int k = 0; k < 8; k++)
			{
				if (FrameI.FHs[nh].FV_Ids[k] == i)
				{
					switch (k)
					{
					case 0:
						PCs[nh].Parametric_coords[0][0][0].index = V_Ind;
						break;
					case 1:
						PCs[nh].Parametric_coords[0][0][Ws[nh]].index = V_Ind;
						break;
					case 2:
						PCs[nh].Parametric_coords[0][Vs[nh]][Ws[nh]].index = V_Ind;
						break;
					case 3:
						PCs[nh].Parametric_coords[0][Vs[nh]][0].index = V_Ind;
						break;
					case 4:
						PCs[nh].Parametric_coords[Us[nh]][0][0].index = V_Ind;
						break;
					case 5:
						PCs[nh].Parametric_coords[Us[nh]][0][Ws[nh]].index = V_Ind;
						break;
					case 6:
						PCs[nh].Parametric_coords[Us[nh]][Vs[nh]][Ws[nh]].index = V_Ind;
						break;
					case 7:
						PCs[nh].Parametric_coords[Us[nh]][Vs[nh]][0].index = V_Ind;
						break;
					}
				}
			}
		}
		V_Ind++;
	}

	for (int i = 0; i < FrameI.FEs.size(); i++)
	{
		int v1 = FrameI.FEs[i].startend_Id[0], v2 = FrameI.FEs[i].startend_Id[1];

		int num = 0;
		for (int j = 0; j < FrameI.FEs[i].neighbor_Hs.size(); j++)
		{
			int nh = FrameI.FEs[i].neighbor_Hs[j];
			int vi, vj;
			for (int k = 0; k < 8; k++)
			{
				if (FrameI.FHs[nh].FV_Ids[k] == v1)
					vi = k;
				if (FrameI.FHs[nh].FV_Ids[k] == v2)
					vj = k;
			}

			int U1, V1, W1, U2, V2, W2;
			if (corners[vi][0] == 0)
				U1 = 0;
			else
				U1 = Us[nh];
			if (corners[vj][0] == 0)
				U2 = 0;
			else
				U2 = Us[nh];
			if (corners[vi][1] == 0)
				V1 = 0;
			else
				V1 = Vs[nh];
			if (corners[vj][1] == 0)
				V2 = 0;
			else
				V2 = Vs[nh];
			if (corners[vi][2] == 0)
				W1 = 0;
			else
				W1 = Ws[nh];
			if (corners[vj][2] == 0)
				W2 = 0;
			else
				W2 = Ws[nh];

			int v_ind = V_Ind;

			if (U1 < U2)
				for (int k = U1 + 1; k < U2; k++)
					PCs[nh].Parametric_coords[k][V1][W1].index = v_ind++;
			else if (U1 > U2)
				for (int k = U1 - 1; k > U2; k--)
					PCs[nh].Parametric_coords[k][V1][W1].index = v_ind++;

			if (V1 < V2)
				for (int k = V1 + 1; k < V2; k++)
					PCs[nh].Parametric_coords[U1][k][W1].index = v_ind++;
			else if (V1 > V2)
				for (int k = V1 - 1; k > V2; k--)
					PCs[nh].Parametric_coords[U1][k][W1].index = v_ind++;

			if (W1 < W2)
				for (int k = W1 + 1; k < W2; k++)
					PCs[nh].Parametric_coords[U1][V1][k].index = v_ind++;
			else if (W1 > W2)
				for (int k = W1 - 1; k > W2; k--)
					PCs[nh].Parametric_coords[U1][V1][k].index = v_ind++;
			num = v_ind;
		}
		V_Ind = num;
	}

	for (int i = 0; i < FrameI.FFs.size(); i++)
	{
		int num = 0;
		for (int j = 0; j < FrameI.FFs[i].neighbor_Cs.size(); j++)
		{
			int nh = FrameI.FFs[i].neighbor_Cs[j];

			vector<int> v_ids;
			for (int m = 0; m < 4; m++)
			{
				for (int k = 0; k < 8; k++)
				{
					if (FrameI.FHs[nh].FV_Ids[k] == FrameI.FFs[i].fv_Ids[m])
						v_ids.push_back(k);
				}
			}
			float min_u = 1, max_u = 0, min_v = 1, max_v = 0, min_w = 1, max_w = 0;
			for (int k = 0; k < 4; k++)
			{
				if (corners[v_ids[k]][0] > max_u)
					max_u = corners[v_ids[k]][0];
				if (corners[v_ids[k]][0] < min_u)
					min_u = corners[v_ids[k]][0];

				if (corners[v_ids[k]][1] > max_v)
					max_v = corners[v_ids[k]][1];
				if (corners[v_ids[k]][1] < min_v)
					min_v = corners[v_ids[k]][1];

				if (corners[v_ids[k]][2] > max_w)
					max_w = corners[v_ids[k]][2];
				if (corners[v_ids[k]][2] < min_w)
					min_w = corners[v_ids[k]][2];
			}
			if (min_u == 1)
				min_u = Us[nh];
			if (max_u == 1)
				max_u = Us[nh];
			if (min_v == 1)
				min_v = Vs[nh];
			if (max_v == 1)
				max_v = Vs[nh];
			if (min_w == 1)
				min_w = Ws[nh];
			if (max_w == 1)
				max_w = Ws[nh];

			int v_ind = V_Ind;

			if (min_u == max_u)
			{
				if (corners[v_ids[0]][1] == 0 && corners[v_ids[0]][2] == 0)
				{
					if (corners[v_ids[0]][1] != corners[v_ids[1]][1])
					{
						for (int k = min_v + 1; k < max_v; k++)
							for (int m = min_w + 1; m < max_w; m++)
								PCs[nh].Parametric_coords[min_u][k][m].index = v_ind++;
					}
					else
					{
						for (int m = min_w + 1; m < max_w; m++)
							for (int k = min_v + 1; k < max_v; k++)
								PCs[nh].Parametric_coords[min_u][k][m].index = v_ind++;
					}
				}
				if (corners[v_ids[0]][1] == 1 && corners[v_ids[0]][2] == 0)
				{
					if (corners[v_ids[0]][1] != corners[v_ids[1]][1])
					{
						for (int k = max_v - 1; k >= min_v + 1; k--)
							for (int m = min_w + 1; m < max_w; m++)
								PCs[nh].Parametric_coords[min_u][k][m].index = v_ind++;
					}
					else
					{
						for (int m = min_w + 1; m < max_w; m++)
							for (int k = max_v - 1; k >= min_v + 1; k--)
								PCs[nh].Parametric_coords[min_u][k][m].index = v_ind++;
					}
				}
				if (corners[v_ids[0]][1] == 1 && corners[v_ids[0]][2] == 1)
				{
					if (corners[v_ids[0]][1] != corners[v_ids[1]][1])
					{
						for (int k = max_v - 1; k >= min_v + 1; k--)
							for (int m = max_w - 1; m >= min_w + 1; m--)
								PCs[nh].Parametric_coords[min_u][k][m].index = v_ind++;
					}
					else
					{
						for (int m = max_w - 1; m >= min_w + 1; m--)
							for (int k = max_v - 1; k >= min_v + 1; k--)
								PCs[nh].Parametric_coords[min_u][k][m].index = v_ind++;
					}
				}
				if (corners[v_ids[0]][1] == 0 && corners[v_ids[0]][2] == 1)
				{
					if (corners[v_ids[0]][1] != corners[v_ids[1]][1])
					{
						for (int k = min_v + 1; k < max_v; k++)
							for (int m = max_w - 1; m >= min_w + 1; m--)
								PCs[nh].Parametric_coords[min_u][k][m].index = v_ind++;
					}
					else
					{
						for (int m = max_w - 1; m >= min_w + 1; m--)
							for (int k = min_v + 1; k < max_v; k++)
								PCs[nh].Parametric_coords[min_u][k][m].index = v_ind++;
					}
				}
			}
			if (min_v == max_v)
			{
				if (corners[v_ids[0]][0] == 0 && corners[v_ids[0]][2] == 0)
				{
					if (corners[v_ids[0]][0] != corners[v_ids[1]][0])
					{
						for (int k = min_u + 1; k < max_u; k++)
							for (int m = min_w + 1; m < max_w; m++)
								PCs[nh].Parametric_coords[k][min_v][m].index = v_ind++;
					}
					else
					{
						for (int m = min_w + 1; m < max_w; m++)
							for (int k = min_u + 1; k < max_u; k++)
								PCs[nh].Parametric_coords[k][min_v][m].index = v_ind++;
					}
				}
				if (corners[v_ids[0]][0] == 1 && corners[v_ids[0]][2] == 0)
				{
					if (corners[v_ids[0]][0] != corners[v_ids[1]][0])
					{
						for (int k = max_u - 1; k >= min_u + 1; k--)
							for (int m = min_w + 1; m < max_w; m++)
								PCs[nh].Parametric_coords[k][min_v][m].index = v_ind++;
					}
					else
					{
						for (int m = min_w + 1; m < max_w; m++)
							for (int k = max_u - 1; k >= min_u + 1; k--)
								PCs[nh].Parametric_coords[k][min_v][m].index = v_ind++;
					}
				}
				if (corners[v_ids[0]][0] == 1 && corners[v_ids[0]][2] == 1)
				{
					if (corners[v_ids[0]][0] != corners[v_ids[1]][0])
					{
						for (int k = max_u - 1; k >= min_u + 1; k--)
							for (int m = max_w - 1; m >= min_w + 1; m--)
								PCs[nh].Parametric_coords[k][min_v][m].index = v_ind++;
					}
					else
					{
						for (int m = max_w - 1; m >= min_w + 1; m--)
							for (int k = max_u - 1; k >= min_u + 1; k--)
								PCs[nh].Parametric_coords[k][min_v][m].index = v_ind++;
					}
				}
				if (corners[v_ids[0]][0] == 0 && corners[v_ids[0]][2] == 1)
				{
					if (corners[v_ids[0]][0] != corners[v_ids[1]][0])
					{
						for (int k = min_u + 1; k < max_u; k++)
							for (int m = max_w - 1; m >= min_w + 1; m--)
								PCs[nh].Parametric_coords[k][min_v][m].index = v_ind++;
					}
					else
					{
						for (int m = max_w - 1; m >= min_w + 1; m--)
							for (int k = min_u + 1; k < max_u; k++)
								PCs[nh].Parametric_coords[k][min_v][m].index = v_ind++;
					}
				}
			}
			if (min_w == max_w)
			{
				if (corners[v_ids[0]][0] == 0 && corners[v_ids[0]][1] == 0)
				{
					if (corners[v_ids[0]][0] != corners[v_ids[1]][0])
					{
						for (int k = min_u + 1; k < max_u; k++)
							for (int m = min_v + 1; m < max_v; m++)
								PCs[nh].Parametric_coords[k][m][min_w].index = v_ind++;
					}
					else
					{
						for (int m = min_v + 1; m < max_v; m++)
							for (int k = min_u + 1; k < max_u; k++)
								PCs[nh].Parametric_coords[k][m][min_w].index = v_ind++;
					}
				}
				if (corners[v_ids[0]][0] == 1 && corners[v_ids[0]][1] == 0)
				{
					if (corners[v_ids[0]][0] != corners[v_ids[1]][0])
					{
						for (int k = max_u - 1; k >= min_u + 1; k--)
							for (int m = min_v + 1; m < max_v; m++)
								PCs[nh].Parametric_coords[k][m][min_w].index = v_ind++;
					}
					else
					{
						for (int m = min_v + 1; m < max_v; m++)
							for (int k = max_u - 1; k >= min_u + 1; k--)
								PCs[nh].Parametric_coords[k][m][min_w].index = v_ind++;
					}
				}
				if (corners[v_ids[0]][0] == 1 && corners[v_ids[0]][1] == 1)
				{
					if (corners[v_ids[0]][0] != corners[v_ids[1]][0])
					{
						for (int k = max_u - 1; k >= min_u + 1; k--)
							for (int m = max_v - 1; m >= min_v + 1; m--)
								PCs[nh].Parametric_coords[k][m][min_w].index = v_ind++;
					}
					else
					{
						for (int m = max_v - 1; m >= min_v + 1; m--)
							for (int k = max_u - 1; k >= min_u + 1; k--)
								PCs[nh].Parametric_coords[k][m][min_w].index = v_ind++;
					}
				}
				if (corners[v_ids[0]][0] == 0 && corners[v_ids[0]][1] == 1)
				{
					if (corners[v_ids[0]][0] != corners[v_ids[1]][0])
					{
						for (int k = min_u + 1; k < max_u; k++)
							for (int m = max_v - 1; m >= min_v + 1; m--)
								PCs[nh].Parametric_coords[k][m][min_w].index = v_ind++;
					}
					else
					{
						for (int m = max_v - 1; m >= min_v + 1; m--)
							for (int k = min_u + 1; k < max_u; k++)
								PCs[nh].Parametric_coords[k][m][min_w].index = v_ind++;
					}
				}
			}
			num = v_ind;
		}
		if (num != 0)
			V_Ind = num;
		else
		{
			printf("i %d: nh size %d\n", i, FrameI.FFs[i].neighbor_Cs.size());
			for (int j = 0; j < FrameI.FFs[i].neighbor_Cs.size(); j++)
				printf("i %d: nh1 %d, nh2 %d\n", i, FrameI.FFs[i].neighbor_Cs[j]);

		}
		if (i == FrameI.FFs.size() - 1)
			continue;
	}

	for (int i = 0; i < FrameI.FHs.size(); i++)
	{
		for (int j = 1; j < Us[i]; j++)
		{
			for (int k = 1; k < Vs[i]; k++)
			{
				for (int m = 1; m < Ws[i]; m++)
					PCs[i].Parametric_coords[j][k][m].index = V_Ind++;

			}
		}
	}
	for (int i = 0; i < V_Ind; i++)
	{
		Hex_V v;
		v.index = i;
		v.where_location = -1;
		hvs.push_back(v);
	}

	for (int i = 0; i < PCs.size(); i++)
	{
		for (int j = 0; j < Us[i] + 1; j++)
		{
			for (int k = 0; k < Vs[i] + 1; k++)
			{
				for (int m = 0; m < Ws[i] + 1; m++)
				{
					hvs[PCs[i].Parametric_coords[j][k][m].index].v[0] = PCs[i].Parametric_coords[j][k][m].v[0];
					hvs[PCs[i].Parametric_coords[j][k][m].index].v[1] = PCs[i].Parametric_coords[j][k][m].v[1];
					hvs[PCs[i].Parametric_coords[j][k][m].index].v[2] = PCs[i].Parametric_coords[j][k][m].v[2];
				}
			}
		}
	}
	for (int i = 0; i < PCs.size(); i++)
	{
		for (int j = 0; j < Us[i]; j++)
		{
			for (int k = 0; k < Vs[i]; k++)
			{
				for (int m = 0; m < Ws[i]; m++)
				{
					Hex h;
					h.index = H_Ind++;
					h.V_Ids[0] = PCs[i].Parametric_coords[j][k][m].index;
					h.V_Ids[1] = PCs[i].Parametric_coords[j][k + 1][m].index;
					h.V_Ids[2] = PCs[i].Parametric_coords[j][k + 1][m + 1].index;
					h.V_Ids[3] = PCs[i].Parametric_coords[j][k][m + 1].index;

					h.V_Ids[4] = PCs[i].Parametric_coords[j + 1][k][m].index;
					h.V_Ids[5] = PCs[i].Parametric_coords[j + 1][k + 1][m].index;
					h.V_Ids[6] = PCs[i].Parametric_coords[j + 1][k + 1][m + 1].index;
					h.V_Ids[7] = PCs[i].Parametric_coords[j + 1][k][m + 1].index;

					hhs.push_back(h);

					hvs[h.V_Ids[0]].neighbor_Hs.push_back(h.index);
					hvs[h.V_Ids[1]].neighbor_Hs.push_back(h.index);
					hvs[h.V_Ids[2]].neighbor_Hs.push_back(h.index);
					hvs[h.V_Ids[3]].neighbor_Hs.push_back(h.index);
					hvs[h.V_Ids[4]].neighbor_Hs.push_back(h.index);
					hvs[h.V_Ids[5]].neighbor_Hs.push_back(h.index);
					hvs[h.V_Ids[6]].neighbor_Hs.push_back(h.index);
					hvs[h.V_Ids[7]].neighbor_Hs.push_back(h.index);
				}
			}
		}
	}

	if (enlarged)
	{
		for (int i = 0; i < hex_mesh.HVs.size(); i++)
		{
			hex_mesh.HVs[i].v[0] /= Escalar;
			hex_mesh.HVs[i].v[1] /= Escalar;
			hex_mesh.HVs[i].v[2] /= Escalar;
		}
		for (int i = 0; i < hvs.size(); i++)
		{
			hvs[i].v[0] /= Escalar;
			hvs[i].v[1] /= Escalar;
			hvs[i].v[2] /= Escalar;
		}
	}

	char fname[300];

	h_io io;
	sprintf(fname, "%s%s", path, ".off");
	io.write_hex_mesh_off(hvs, hhs, fname);
}
