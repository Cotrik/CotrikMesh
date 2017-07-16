#include "parametric_optimization.h"
#include <cmath>

parametric_optimization::parametric_optimization(void)
{
	File_Index = 0;
}
void parametric_optimization::optimization_pipeline(char *fname)
{
	float hex_average_len = hex_mesh.average_e_len;

	bool enlarged = false;
	float Escalar = 1.0;
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
		hex_average_len = hex_mesh.average_e_len;
		printf("Enlarged\n");
	}

	clock_t start_time = clock();

	h_io io;

	int Iteration = 5;
	for (int i = 0; i < Iteration; i++)
	{
		laplacian_smoothing laps;
		laps.Iteration = 1;
		//laps.laplacian_global_pipeline(hex_mesh.HVs,hex_mesh.HFs,hex_mesh.HHs);
		laps.project_surface_global(hex_mesh.HVs);
		//optimization_pipeline_face(i,true);
		//laps.project_surface_global(hex_mesh.HVs);
		optimization_pipeline_edge(i);
		//laps.project_surface_global(hex_mesh.HVs);
		//optimization_pipeline_face(i,false);
		//optimization_pipeline_edgeface(i,fname);
		//laps.project_surface_global(hex_mesh.HVs);

		if (enlarged)
		{
			enlarged = false;
			for (int i = 0; i < hex_mesh.HVs.size(); i++)
			{
				hex_mesh.HVs[i].v[0] /= Escalar;
				hex_mesh.HVs[i].v[1] /= Escalar;
				hex_mesh.HVs[i].v[2] /= Escalar;
			}
			hex_mesh.average_e_len /= Escalar;
			hex_average_len = hex_mesh.average_e_len;
			printf("Shrinked\n");
		}

		if (1)
		{
			char pathH[300];
			sprintf(pathH, "%s%s", fname, "Iter5_temp.off");
			io.write_hex_mesh_off(hex_mesh.HVs, hex_mesh.HHs, pathH);

			initialization_parameters();
			{
				io.read_hex_mesh_off(hex_mesh.HVs, hex_mesh.HHs, pathH);

				construct_Es(hex_mesh.HVs, hex_mesh.HEs, hex_mesh.HHs);
				construct_Fs(hex_mesh.HVs, hex_mesh.HEs, hex_mesh.HFs, hex_mesh.HHs);
				determine_boundary_info(hex_mesh.HVs, hex_mesh.HEs, hex_mesh.HFs, hex_mesh.HHs);

				laps.project_surface_global(hex_mesh.HVs);

				hex_mesh.average_e_len = hex_average_len;
				calculation_centroid(hex_mesh.HVs, hex_mesh.HFs, hex_mesh.HHs);
			}
			frame_of_mesh fom;
			fom.base_complex_extraction();
		}
		if (1)
		{
			printf("start Parmaterization\n");
			parameterization pa;
			Para_min_N = 1;
			sprintf(fname, "%s%d", fname, Iteration);
			pa.parameterization_main(fname);

			initialization_parameters();
			{
				char path[300];
				sprintf(path, "%s%s", fname, ".off");
				io.read_hex_mesh_off(hex_mesh.HVs, hex_mesh.HHs, path);

				construct_Es(hex_mesh.HVs, hex_mesh.HEs, hex_mesh.HHs);
				construct_Fs(hex_mesh.HVs, hex_mesh.HEs, hex_mesh.HFs, hex_mesh.HHs);
				determine_boundary_info(hex_mesh.HVs, hex_mesh.HEs, hex_mesh.HFs, hex_mesh.HHs);

				laps.project_surface_global(hex_mesh.HVs);

				hex_mesh.average_e_len = hex_average_len;
				calculation_centroid(hex_mesh.HVs, hex_mesh.HFs, hex_mesh.HHs);
			}

			frame_of_mesh fom;
			fom.base_complex_extraction();
		}
	}
	clock_t end_time = clock();
	std::cout << "Running time is: " << static_cast<double>(end_time - start_time) / CLOCKS_PER_SEC * 1000 << "ms" << std::endl;

}
void parametric_optimization::optimization_pipelineF(char *fname)
{
	float hex_average_len = hex_mesh.average_e_len;

	clock_t start_time = clock();

	h_io io;

	int Iteration = 2;
	for (int i = 0; i < Iteration; i++)
	{
		laplacian_smoothing laps;
		laps.Iteration = 1;
		//laps.laplacian_global_pipeline(hex_mesh.HVs,hex_mesh.HFs,hex_mesh.HHs);
		//laps.project_surface_global(hex_mesh.HVs);
		optimization_pipeline_face(i, false);
		//optimization_pipeline_edge(i);

		if (1)
		{
			char pathH[300];
			sprintf(pathH, "%s%s", fname, "Iter5_temp.off");
			io.write_hex_mesh_off(hex_mesh.HVs, hex_mesh.HHs, pathH);

			initialization_parameters();
			{
				io.read_hex_mesh_off(hex_mesh.HVs, hex_mesh.HHs, pathH);

				construct_Es(hex_mesh.HVs, hex_mesh.HEs, hex_mesh.HHs);
				construct_Fs(hex_mesh.HVs, hex_mesh.HEs, hex_mesh.HFs, hex_mesh.HHs);
				determine_boundary_info(hex_mesh.HVs, hex_mesh.HEs, hex_mesh.HFs, hex_mesh.HHs);
				hex_mesh.average_e_len = hex_average_len;
				calculation_centroid(hex_mesh.HVs, hex_mesh.HFs, hex_mesh.HHs);
			}
			//laps.project_surface_global(hex_mesh.HVs);

			frame_of_mesh fom;
			fom.base_complex_extraction();
		}
		if (0)
		{
			printf("start Parmaterization\n");
			parameterization pa;
			Para_min_N = 1;
			sprintf(fname, "%s%d", fname, Iteration);
			pa.parameterization_main(fname);

			initialization_parameters();
			{
				char path[300];
				sprintf(path, "%s%s", fname, ".off");
				io.read_hex_mesh_off(hex_mesh.HVs, hex_mesh.HHs, path);

				construct_Es(hex_mesh.HVs, hex_mesh.HEs, hex_mesh.HHs);
				construct_Fs(hex_mesh.HVs, hex_mesh.HEs, hex_mesh.HFs, hex_mesh.HHs);
				determine_boundary_info(hex_mesh.HVs, hex_mesh.HEs, hex_mesh.HFs, hex_mesh.HHs);
				hex_mesh.average_e_len = hex_average_len;
				calculation_centroid(hex_mesh.HVs, hex_mesh.HFs, hex_mesh.HHs);
			}
			//laps.project_surface_global(hex_mesh.HVs);

			frame_of_mesh fom;
			fom.base_complex_extraction();
		}
	}
	clock_t end_time = clock();
	std::cout << "Running time is: " << static_cast<double>(end_time - start_time) / CLOCKS_PER_SEC * 1000 << "ms" << std::endl;

}
void parametric_optimization::optimization_pipeline_edge(int Iter)
{
	//extract edge chains
	ECs.clear();
	edge_chain_extraction();
	ParaCs.clear();
	parameter_index_correspondence_new();
	vector<bool> ecs_Ind;
	initializeVectorT(ecs_Ind, true, ECs.size());

	vector<vector<int>> ECs_idscurS;
	vector<int> ECs_idscur;
	find_ecs_cur(ECs_idscur, ecs_Ind);
	while (ECs_idscur.size())
	{
		ECs_idscurS.push_back(ECs_idscur);
		find_ecs_cur(ECs_idscur, ecs_Ind);
	}

	int which_circle = 1;
	for (int i = 0; i < ECs_idscurS.size(); i++)
	{
		printf("iteration %d circle %d\n", Iter, which_circle++);
#pragma omp parallel for //multi-thread		for (int j = 0; j < ECs_idscurS[i].size(); j++)
		{
			a_chain_parameterization(ECs_idscurS[i][j]);
		}
	}
}
void parametric_optimization::optimization_pipeline_face(int Iter, bool smooth)
{
	//extract face chains
	FCs.clear();
	face_center_extraction();
	ParaFCs.clear();
	face_parameter_index_correspondence();
	vector<bool> fcs_Ind;
	initializeVectorT(fcs_Ind, true, FCs.size());

	vector<vector<int>> FCs_idscurS;
	vector<int> FCs_idscur;
	find_fcs_cur(FCs_idscur, fcs_Ind);
	while (FCs_idscur.size())
	{
		FCs_idscurS.push_back(FCs_idscur);
		find_fcs_cur(FCs_idscur, fcs_Ind);
	}

	int which_circle = 0;
	for (int i = 0; i < FCs_idscurS.size(); i++)
	{
		printf("iteration %d circle %d\n", Iter, which_circle++);
#pragma omp parallel for //multi-thread		for (int j = 0; j < FCs_idscurS[i].size(); j++)
		{
			face_center_parameterization(FCs_idscurS[i][j]);
		}
		//project to surface
		vector<Hex_V> hvs;
		hvs = hex_mesh.HVs;
		for (int j = 0; j < hvs.size(); j++)
			hvs[j].where_location = -1;
		for (int j = 0; j < FCs_idscurS[i].size(); j++)
		{
			for (int k = 0; k < FCs[FCs_idscurS[i][j]].fs_boundary.size(); k++)
			{
				int bfid = FCs[FCs_idscurS[i][j]].fs_boundary[k];
				if (FrameI.FFs[bfid].is_boundary)
				{
					for (int m = 0; m < FrameI.FFs[bfid].hfs_net_another.size(); m++)
					{
						int hffid = FrameI.FFs[bfid].hfs_net_another[m];
						for (int n = 0; n < 4; n++)
						{
							int hffvid = hex_mesh.HFs[hffid].cv_Ids[n];
							hvs[hffvid].where_location = 1;
						}
					}
				}
			}
		}
		laplacian_smoothing laps;
		laps.project_surface_global(hvs);
		for (int j = 0; j < hvs.size(); j++)
		{
			if (hvs[j].where_location == 1)
				hex_mesh.HVs[j].v[0] = hvs[j].v[0];
			hex_mesh.HVs[j].v[1] = hvs[j].v[1];
			hex_mesh.HVs[j].v[2] = hvs[j].v[2];
		}
		//smoothing or not
		if (!smooth)
			continue;

		vector<int> vs_all_cur, vs_all_cur2;
		for (int j = 0; j < FCs_idscurS[i].size(); j++)
		{
			append_vector(vs_all_cur, FCs[FCs_idscurS[i][j]].vs_all);
		}
		//vector<Hex_V> hvs;
		hvs = hex_mesh.HVs;
		for (int j = 0; j < vs_all_cur.size(); j++)
		{
			vs_all_cur[j] = FrameI.FVs[vs_all_cur[j]].index_hex;
			if (hex_mesh.HVs[vs_all_cur[j]].where_location != 1)
				vs_all_cur2.push_back(vs_all_cur[j]);
		}
		for (int j = 0; j < hvs.size(); j++)
			hvs[j].where_location = 1;
		for (int j = 0; j < vs_all_cur2.size(); j++)
			hvs[vs_all_cur2[j]].where_location = -1;
		//laplacian_smoothing laps;
		for (int j = 0; j < 3; j++)
			laps.laplacian_inner_volume(hvs, hex_mesh.HHs, 1);
		for (int j = 0; j < vs_all_cur2.size(); j++)
		{
			hex_mesh.HVs[vs_all_cur2[j]].v[0] = hvs[vs_all_cur2[j]].v[0];
			hex_mesh.HVs[vs_all_cur2[j]].v[1] = hvs[vs_all_cur2[j]].v[1];
			hex_mesh.HVs[vs_all_cur2[j]].v[2] = hvs[vs_all_cur2[j]].v[2];
		}
	}
}
void parametric_optimization::optimization_pipeline_edgeface(int Iter, char *fname)
{
	//extract edge chains
	ECs.clear();
	edge_chain_extraction();
	ParaCs.clear();
	parameter_index_correspondence_new();
	vector<bool> ecs_Ind;
	initializeVectorT(ecs_Ind, true, ECs.size());

	//extract face chains
	FCs.clear();
	face_center_extraction();
	ParaFCs.clear();
	face_parameter_index_correspondence();
	vector<bool> fcs_Ind;
	initializeVectorT(fcs_Ind, true, FCs.size());

	vector<vector<int>> ECs_idscurS, FCs_idscurS;
	vector<int> ECs_idscur, FCs_idscur;
	find_efcs_cur(ECs_idscur, FCs_idscur, ecs_Ind, fcs_Ind);
	while (FCs_idscur.size() || ECs_idscur.size())
	{
		ECs_idscurS.push_back(ECs_idscur);
		FCs_idscurS.push_back(FCs_idscur);
		find_efcs_cur(ECs_idscur, FCs_idscur, ecs_Ind, fcs_Ind);
	}

	h_io io;

	int which_circle = 1;
	for (int i = 0; i < FCs_idscurS.size(); i++)
	{
		////assign edge domain colors
		for (int j = 0; j < hex_mesh.HHs.size(); j++)
			hex_mesh.HHs[j].Color_EF_ID = 0;

		int Ind_C = 0;
		int Base_ = 8;
		for (int m = 0; m < ECs_idscurS[i].size(); m++)
		{
			Ind_C++;
			Ind_C = Ind_C % Base_;

			for (int j = 0; j < ECs[ECs_idscurS[i][m]].bhs.size(); j++)
			{
				for (int k = 0; k < ECs[ECs_idscurS[i][m]].bhs[j].size(); k++)
				{
					int fhid = ECs[ECs_idscurS[i][m]].bhs[j][k];
					for (int n = 0; n < FrameI.FHs[fhid].hs_net.size(); n++)
					{
						int hid = FrameI.FHs[fhid].hs_net[n];
						hex_mesh.HHs[hid].Color_EF_ID = Ind_C;
					}
				}
			}
		}
		for (int m = 0; m < FCs_idscurS[i].size(); m++)
		{
			Ind_C++;
			Ind_C = Ind_C % Base_ + 10;
			vector<int> bhs = FrameI.FFs[FCs[FCs_idscurS[i][m]].f_id].neighbor_Cs;
			for (int j = 0; j < bhs.size(); j++)
			{
				int fhid = bhs[j];
				for (int n = 0; n < FrameI.FHs[fhid].hs_net.size(); n++)
				{
					int hid = FrameI.FHs[fhid].hs_net[n];
					hex_mesh.HHs[hid].Color_EF_ID = Ind_C;
				}
			}
		}

		char path[300];
		sprintf(path, "%s%d%s%d%s", fname, Iter, "_part_", i, "_before.vtk");
		//io.write_VTK(hex_mesh.HVs,hex_mesh.HHs,path);

		printf("iteration %d circle %d\n", Iter, which_circle++);
#pragma omp parallel for //multi-thread		for (int j = 0; j < ECs_idscurS[i].size(); j++)
		{
			a_chain_parameterization(ECs_idscurS[i][j]);
		}
		printf("iteration %d circle %d\n", Iter, which_circle++);
#pragma omp parallel for //multi-thread		for (int j = 0; j < FCs_idscurS[i].size(); j++)
		{
			face_center_parameterization(FCs_idscurS[i][j]);
		}
	}
}
void parametric_optimization::find_efcs_cur(vector<int> &ecs_ids, vector<int> &fcs_ids, vector<bool> &ecs_Ind, vector<bool> &fcs_Ind)
{
	ecs_ids.clear();
	for (int i = 0; i < ECs.size(); i++)
	{
		if (ecs_Ind[i])
		{
			bool no_inter = true;
			for (int j = 0; j < ecs_ids.size(); j++)
			{
				vector<int> vectors1, vectors2, vectors3, vectors4;
				append_vector(vectors1, ECs[i].boundary_fs);
				for (int k = 0; k < ECs[i].bfs.size(); k++)
					append_vector(vectors1, ECs[i].bfs[k]);
				append_vector(vectors2, ECs[ecs_ids[j]].boundary_fs);
				for (int k = 0; k < ECs[ecs_ids[j]].bfs.size(); k++)
					append_vector(vectors2, ECs[ecs_ids[j]].bfs[k]);
				set_exclusion(vectors1, vectors2, vectors3);
				if (vectors1.size() != vectors3.size())
				{
					no_inter = false;
					break;
				}
			}
			if (no_inter)
			{
				ecs_ids.push_back(i);
				ecs_Ind[i] = false;
			}
		}
	}

	fcs_ids.clear();
	for (int i = 0; i < FCs.size(); i++)
	{
		if (fcs_Ind[i])
		{
			bool no_inter = true;
			for (int j = 0; j < fcs_ids.size(); j++)
			{
				vector<int> excfs;
				set_exclusion(FCs[i].fs_boundary, FCs[fcs_ids[j]].fs_boundary, excfs);
				if (FCs[i].fs_boundary.size() != excfs.size())
				{
					no_inter = false;
					break;
				}
			}
			for (int j = 0; j < ecs_ids.size(); j++)
			{
				vector<int> vectors1, vectors2, vectors3, vectors4;
				append_vector(vectors1, FCs[i].fs_boundary);
				append_vector(vectors1, FCs[i].fs_paralell);
				set_redundent_clearn(vectors1);
				append_vector(vectors2, ECs[ecs_ids[j]].boundary_fs);
				for (int k = 0; k < ECs[ecs_ids[j]].bfs.size(); k++)
					append_vector(vectors2, ECs[ecs_ids[j]].bfs[k]);
				set_exclusion(vectors1, vectors2, vectors3);
				if (vectors1.size() != vectors3.size())
				{
					no_inter = false;
					break;
				}
			}

			if (no_inter)
			{
				fcs_ids.push_back(i);
				fcs_Ind[i] = false;
			}
		}
	}
}

//edge
void parametric_optimization::edge_chain_extraction()
{
	vector<bool> bes_Inds;
	initializeVectorT(bes_Inds, false, FrameI.FEs.size());
	bool stillleft = true;
	while (stillleft)
	{
		int starteid = -1;
		for (int i = 0; i < bes_Inds.size(); i++)
		{
			if (!bes_Inds[i])
			{
				starteid = i;
				break;
			}
		}
		if (starteid == -1)
			break;

		bool inside = false;
		for (int i = 0; i < FrameI.FEs[starteid].vs_link.size(); i++)
		{
			int vid = FrameI.FEs[starteid].vs_link[i];
			if (hex_mesh.HVs[vid].where_location == -1)
				inside = true;
		}
		if (!inside)
		{
			bes_Inds[starteid] = true;
			continue;
		}
		// 		if(FrameI.FEs[starteid].neighbor_Hs.size()!=4)
		// 		{
		// 			bes_Inds[starteid]=true;
		// 			continue;
		// 		}

		vector<int> vs_total, es_total, fs_total, hs_total;
		append_vector(hs_total, FrameI.FEs[starteid].neighbor_Hs);
		vector<int> fs_1ring;
		base_complex_e_onering_fs(starteid, fs_1ring);
		append_vector(fs_total, fs_1ring);

		int startv1 = FrameI.FEs[starteid].startend_Id[0];
		int startv2 = FrameI.FEs[starteid].startend_Id[1];
		vs_total.push_back(startv1);
		es_total.push_back(starteid);

		//towards 1
		int startv1_next = -1;
		int next_e1 = -1;
		while (FrameI.FVs[startv1].which_singularity == -1 && hex_mesh.HVs[FrameI.FVs[startv1].index_hex].where_location == -1)
		{
			for (int i = 0; i < FrameI.FVs[startv1].neighbor_Es.size(); i++)
			{
				int neid = FrameI.FVs[startv1].neighbor_Es[i];
				if (bes_Inds[neid])
					continue;
				vector<int> crosshs;
				set_cross(FrameI.FEs[starteid].neighbor_Hs, FrameI.FEs[neid].neighbor_Hs, crosshs);
				if (!crosshs.size())
				{
					next_e1 = neid;
					break;
				}
			}
			if (next_e1 == -1)
				break;

			startv1_next = FrameI.FEs[next_e1].startend_Id[0];
			if (FrameI.FEs[next_e1].startend_Id[0] == startv1)
				startv1_next = FrameI.FEs[next_e1].startend_Id[1];

			if (startv1_next == startv2)
				break;		//going to be form a circle
			vector<int> crosshs;
			set_cross(hs_total, FrameI.FEs[next_e1].neighbor_Hs, crosshs);
			if (crosshs.size())
				break;		//cross hex components

			base_complex_e_onering_fs(next_e1, fs_1ring);
			vector<int> crossfs;
			set_cross(fs_total, fs_1ring, crosshs);
			if (crossfs.size())
				break;		//neighboring hex components

			es_total.push_back(next_e1);
			vs_total.push_back(startv1_next);
			append_vector(fs_total, fs_1ring);
			append_vector(hs_total, FrameI.FEs[next_e1].neighbor_Hs);

			startv1 = startv1_next;
			starteid = next_e1;
			next_e1 = -1;
		}

		//towards 2
		startv1 = vs_total[vs_total.size() - 1];
		starteid = es_total[0];

		vector<int> vs2, es2, fs2, hs2;
		vs2.push_back(startv2);
		int startv2_next = -1;
		int next_e2 = -1;
		while (FrameI.FVs[startv2].which_singularity == -1 && hex_mesh.HVs[FrameI.FVs[startv2].index_hex].where_location == -1)
		{

			for (int i = 0; i < FrameI.FVs[startv2].neighbor_Es.size(); i++)
			{
				int neid = FrameI.FVs[startv2].neighbor_Es[i];
				if (bes_Inds[neid])
					continue;
				vector<int> crosshs;
				set_cross(FrameI.FEs[starteid].neighbor_Hs, FrameI.FEs[neid].neighbor_Hs, crosshs);
				if (!crosshs.size())
				{
					next_e2 = neid;
					break;
				}
			}
			if (next_e2 == -1)
				break;
			startv2_next = FrameI.FEs[next_e2].startend_Id[0];
			if (FrameI.FEs[next_e2].startend_Id[0] == startv2)
				startv2_next = FrameI.FEs[next_e2].startend_Id[1];

			if (startv2_next == startv1)
				break;		//going to be form a circle
			vector<int> crosshs;
			set_cross(hs_total, FrameI.FEs[next_e2].neighbor_Hs, crosshs);
			if (crosshs.size())
				break;		//cross hex components

			base_complex_e_onering_fs(next_e2, fs_1ring);
			vector<int> crossfs;
			set_cross(fs_total, fs_1ring, crosshs);
			if (crossfs.size())
				break;		//neighboring hex components

			es2.push_back(next_e1);
			vs2.push_back(startv2_next);
			append_vector(fs_total, fs_1ring);
			append_vector(hs_total, FrameI.FEs[next_e2].neighbor_Hs);

			startv2 = startv2_next;
			starteid = next_e2;
			next_e2 = -1;
		}
		//merge 1 and 2
		edge_chain ec;
		for (int i = es2.size() - 1; i >= 0; i--)
			ec.bes.push_back(es2[i]);
		for (int i = 0; i < es_total.size(); i++)
			ec.bes.push_back(es_total[i]);
		for (int i = vs2.size() - 1; i >= 0; i--)
			ec.bvs.push_back(vs2[i]);
		for (int i = 0; i < vs_total.size(); i++)
			ec.bvs.push_back(vs_total[i]);
		ECs.push_back(ec);
		for (int i = 0; i < ec.bes.size(); i++)
			bes_Inds[ec.bes[i]] = true;
	}

	for (int i = 0; i < ECs.size(); i++)
	{
		ECs[i].id = i;
		ECs[i].startendIds[0] = ECs[i].bvs[0];
		ECs[i].startendIds[1] = ECs[i].bvs[ECs[i].bvs.size() - 1];

		//valence fs, each layer are not guaranteeed to have the same start f
		for (int j = 0; j < ECs[i].bes.size(); j++)
		{
			vector<int> valence_fs;
			valence_fs.push_back(FrameI.FEs[ECs[i].bes[j]].neighbor_Fs[0]);
			int cur_pos = valence_fs[0];
			bool found = true;
			while (found)
			{
				found = false;
				for (int k = 0; k < FrameI.FEs[ECs[i].bes[j]].neighbor_Fs.size(); k++)
				{
					int nfid = FrameI.FEs[ECs[i].bes[j]].neighbor_Fs[k];
					if (set_contain(valence_fs, nfid) == -1)
					{
						vector<int> cross_hs;
						set_cross(FrameI.FFs[nfid].neighbor_Cs, FrameI.FFs[cur_pos].neighbor_Cs, cross_hs);
						if (!cross_hs.size())
							continue;
						valence_fs.push_back(nfid);
						cur_pos = nfid;
						found = true;
						break;
					}
				}
			}
			ECs[i].bfs.push_back(valence_fs);
		}
		//valence vs, boundary vs
		for (int j = 0; j < ECs[i].bvs.size(); j++)
		{
			int vid = ECs[i].bvs[j];
			int eid = j;
			if (j == ECs[i].bvs.size() - 1)
				eid--;

			vector<int> hnvs;
			for (int k = 0; k < FrameI.FEs[ECs[i].bes[eid]].neighbor_Hs.size(); k++)
			{
				int nhid = FrameI.FEs[ECs[i].bes[eid]].neighbor_Hs[k];
				for (int m = 0; m < 8; m++)
					hnvs.push_back(FrameI.FHs[nhid].FV_Ids[m]);
			}
			set_redundent_clearn(hnvs);
			vector<int> valence_vs;
			//find start vid
			if (j == 0)
			{
				for (int m = 0; m < ECs[i].bfs[eid].size(); m++)
				{
					int val_fid = ECs[i].bfs[eid][m];
					for (int k = 0; k < FrameI.FVs[vid].neighbor_vs.size(); k++)
					{
						int nvid = FrameI.FVs[vid].neighbor_vs[k];
						if (set_contain(hnvs, nvid) != -1 && set_contain(ECs[i].bvs, nvid) == -1)
						{
							vector<int> fvs;
							fvs.push_back(FrameI.FFs[val_fid].fv_Ids[0]);
							fvs.push_back(FrameI.FFs[val_fid].fv_Ids[1]);
							fvs.push_back(FrameI.FFs[val_fid].fv_Ids[2]);
							fvs.push_back(FrameI.FFs[val_fid].fv_Ids[3]);
							if (set_contain(fvs, nvid) != -1)
							{
								valence_vs.push_back(nvid);
								break;
							}
						}
					}
				}
			}
			else
			{
				for (int k = 0; k < ECs[i].bvs_valence_ring[j - 1].size(); k++)
				{
					int bvvrid = ECs[i].bvs_valence_ring[j - 1][k];
					for (int m = 0; m < FrameI.FVs[bvvrid].neighbor_vs.size(); m++)
					{
						int mnvid = FrameI.FVs[bvvrid].neighbor_vs[m];
						if (set_contain(ECs[i].bvs, mnvid) != -1)
							continue;
						if (set_contain(FrameI.FVs[vid].neighbor_vs, mnvid) != -1)
						{
							valence_vs.push_back(mnvid);
							break;
						}
					}
				}
			}
			ECs[i].bvs_valence_ring.push_back(valence_vs);
			//collect 1-ring vs around vi
			vector<int> ring1_vs;
			for (int k = 0; k < valence_vs.size(); k++)
			{
				ring1_vs.push_back(valence_vs[k]);
				vector<int> cross_vs;
				set_cross(FrameI.FVs[valence_vs[k]].neighbor_vs, FrameI.FVs[valence_vs[(k + 1) % valence_vs.size()]].neighbor_vs, cross_vs);
				for (int m = 0; m < cross_vs.size(); m++)
				{
					if (set_contain(hnvs, cross_vs[m]) != -1 && set_contain(ECs[i].bvs, cross_vs[m]) == -1)
					{
						ring1_vs.push_back(cross_vs[m]);
						break;
					}
				}
			}
			ECs[i].bvs_1ring.push_back(ring1_vs);
		}
		for (int j = 0; j < ECs[i].bes.size(); j++)
		{
			int eid = ECs[i].bes[j];
			vector<int> valence_hs;
			for (int m = 0; m < ECs[i].bvs_valence_ring[j].size(); m++)
			{
				int v1 = ECs[i].bvs_valence_ring[j][m], v2 = ECs[i].bvs_valence_ring[j][(m + 1) % ECs[i].bvs_valence_ring[j].size()];
				for (int k = 0; k < FrameI.FEs[eid].neighbor_Hs.size(); k++)
				{
					int nhid = FrameI.FEs[eid].neighbor_Hs[k];
					bool v1_have = false, v2_have = false;
					for (int n = 0; n < 8; n++)
					{
						if (FrameI.FHs[nhid].FV_Ids[n] == v1)
							v1_have = true;
						if (FrameI.FHs[nhid].FV_Ids[n] == v2)
							v2_have = true;
					}
					if (v1_have && v2_have)
					{
						valence_hs.push_back(nhid);
						break;
					}
				}
			}
			ECs[i].bhs.push_back(valence_hs);
		}

		vector<int> fs_all;
		for (int j = 0; j < ECs[i].bhs.size(); j++)
		{
			for (int k = 0; k < ECs[i].bhs[j].size(); k++)
				append_vector(fs_all, FrameI.FHs[ECs[i].bhs[j][k]].neighbor_FS);
		}

		set_redundent_clearn(fs_all);
		ECs[i].all_fs = fs_all;

		if (ECs[i].bvs.size() > 2)
		{
			for (int j = 1; j < ECs[i].bvs.size() - 1; j++)
			{
				vector<int> resultfs;
				set_exclusion(fs_all, FrameI.FVs[ECs[i].bvs[j]].neighbor_Fs, resultfs);
				fs_all = resultfs;
			}
		}
		else
		{
			for (int j = 0; j < ECs[i].bvs.size(); j++)
			{
				vector<int> resultfs;
				set_exclusion(fs_all, FrameI.FVs[ECs[i].bvs[j]].neighbor_Fs, resultfs);
				fs_all = resultfs;
			}

			vector<int> topbottomvalencefs;
			for (int j = 0; j < FrameI.FVs[ECs[i].bvs[0]].neighbor_Fs.size(); j++)
			{
				int nfid = FrameI.FVs[ECs[i].bvs[0]].neighbor_Fs[j];
				vector<int> v4;
				v4.push_back(FrameI.FFs[nfid].fv_Ids[0]);
				v4.push_back(FrameI.FFs[nfid].fv_Ids[1]);
				v4.push_back(FrameI.FFs[nfid].fv_Ids[2]);
				v4.push_back(FrameI.FFs[nfid].fv_Ids[3]);
				vector<int> uplayervs = ECs[i].bvs_1ring[0];
				uplayervs.push_back(ECs[i].bvs[0]);
				if (set_contain(uplayervs, v4))
					topbottomvalencefs.push_back(nfid);
			}
			for (int j = 0; j < FrameI.FVs[ECs[i].bvs[1]].neighbor_Fs.size(); j++)
			{
				int nfid = FrameI.FVs[ECs[i].bvs[1]].neighbor_Fs[j];
				vector<int> v4;
				v4.push_back(FrameI.FFs[nfid].fv_Ids[0]);
				v4.push_back(FrameI.FFs[nfid].fv_Ids[1]);
				v4.push_back(FrameI.FFs[nfid].fv_Ids[2]);
				v4.push_back(FrameI.FFs[nfid].fv_Ids[3]);
				vector<int> bottomlayervs = ECs[i].bvs_1ring[1];
				bottomlayervs.push_back(ECs[i].bvs[1]);
				if (set_contain(bottomlayervs, v4))
					topbottomvalencefs.push_back(nfid);
			}
			append_vector(fs_all, topbottomvalencefs);
		}
		ECs[i].boundary_fs = fs_all;

		for (int j = 0; j < ECs[i].bvs_1ring.size(); j++)
		{
			vector<int> ringes;
			for (int k = 0; k < ECs[i].bvs_1ring[j].size(); k++)
			{
				int v1 = ECs[i].bvs_1ring[j][k];
				int v2 = ECs[i].bvs_1ring[j][(k + 1) % ECs[i].bvs_1ring[j].size()];
				vector<int> crosses;
				set_cross(FrameI.FVs[v1].neighbor_Es, FrameI.FVs[v2].neighbor_Es, crosses);
				ringes.push_back(crosses[0]);
			}
			ECs[i].bes_1ring.push_back(ringes);
		}
		//boundary es except top and bottom valence es
		vector<int> ringes;
		for (int k = 0; k < ECs[i].boundary_fs.size(); k++)
		{
			int bfid = ECs[i].boundary_fs[k];
			for (int j = 0; j < 4; j++)
				ringes.push_back(FrameI.FFs[bfid].fe_Ids[j]);
		}
		set_redundent_clearn(ringes);
		vector<int> valence_edges;
		for (int k = 0; k < ECs[i].bvs_1ring[0].size(); k++)
		{
			int bvid = ECs[i].bvs_1ring[0][k];
			vector<int> cross_e;
			set_cross(FrameI.FVs[ECs[i].bvs[0]].neighbor_Es, FrameI.FVs[bvid].neighbor_Es, cross_e);
			if (cross_e.size())
				valence_edges.push_back(cross_e[0]);
		}
		for (int k = 0; k < ECs[i].bvs_1ring[ECs[i].bvs.size() - 1].size(); k++)
		{
			int bvid = ECs[i].bvs_1ring[ECs[i].bvs.size() - 1][k];
			vector<int> cross_e;
			set_cross(FrameI.FVs[ECs[i].bvs[ECs[i].bvs.size() - 1]].neighbor_Es, FrameI.FVs[bvid].neighbor_Es, cross_e);
			if (cross_e.size())
				valence_edges.push_back(cross_e[0]);
		}
		//ECs[i].boundaryes_excepttb=ringes;
		vector<int> real_ringes;
		set_exclusion(ringes, valence_edges, real_ringes);
		ECs[i].boundaryes_excepttb = real_ringes;
		for (int k = 1; k < ECs[i].bes_1ring.size() - 1; k++)
		{
			set_exclusion(real_ringes, ECs[i].bes_1ring[k], ECs[i].boundaryes_excepttb);
			real_ringes = ECs[i].boundaryes_excepttb;
		}
		ECs[i].boundaryes_excepttb = real_ringes;

		// 		vector<int> alles;
		// 		for(int k=0;k<ECs[i].bhs.size();k++)
		// 		{
		// 			for(int n=0;n<ECs[i].bhs[k].size();n++)
		// 			{
		// 				int bhid=ECs[i].bhs[k][n];
		// 				for(int j=0;j<6;j++)
		// 				{
		// 					int bfid=FrameI.FHs[bhid].neighbor_FS[j];
		// 					for(int m=0;m<4;m++)
		// 					{
		// 						alles.push_back(FrameI.FFs[bfid].fe_Ids[m]);
		// 					}
		// 				}
		// 			}
		// 		}
		// 		set_redundent_clearn(alles);
		// 		ECs[i].boundaryes_excepttb=alles;
	}
}

void parametric_optimization::parameter_index_correspondence()
{
	for (int i = 0; i < ECs.size(); i++)
	{
//assign frame node parameter coords.
		float len = 0.0;
		for (int j = 0; j < ECs[i].bvs_1ring.size(); j++)
		{
			vector<para_coor> onelayer_ring;
			int sizering = ECs[i].bvs_1ring[j].size();
			if (sizering == 8)
			{
				float para[8][2] =
				{
				{ 1, 0 },
				{ 1, 1 },
				{ 0, 1 },
				{ -1, 1 },
				{ -1, 0 },
				{ -1, -1 },
				{ 0, -1 },
				{ 1, -1 } };
				for (int k = 0; k < sizering; k++)
				{
					para_coor pc;
					pc.coords.push_back(para[k][0]);
					pc.coords.push_back(para[k][1]);
					pc.coords.push_back(len);
					onelayer_ring.push_back(pc);
				}
			}
			else
			{
				for (int k = 0; k < sizering; k++)
				{
					para_coor pc;
					pc.coords.push_back(cos(k * 2 * PI / sizering));
					pc.coords.push_back(sin(k * 2 * PI / sizering));
					pc.coords.push_back(len);
					onelayer_ring.push_back(pc);
				}
			}

			ECs[i].ring1_paracoords.push_back(onelayer_ring);
			len += 1;
		}
//assign all parameter coords
		//printf("This is %dth Ec\n",i);

		Parameterization_Chain ParaC;
		int N_fan = ECs[i].bvs_valence_ring[0].size();
		for (int j = 0; j < N_fan; j++)
		{
			//printf("This is %dth Fan\n",j);
			//determine x, y, z dimensions
			int X_dim, Y_dim;
			vector<int> Z_dims;
			vector<int> xeid;
			vector<int> ids_vx, ids_vy, ids_vz;
			vector<int> ids_vzz;

			set_cross(FrameI.FVs[ECs[i].bvs[0]].neighbor_Es, FrameI.FVs[ECs[i].bvs_valence_ring[0][j]].neighbor_Es, xeid);
			X_dim = FrameI.FEs[xeid[0]].vs_link.size() - 1;

			ids_vx = FrameI.FEs[xeid[0]].vs_link;
			if (FrameI.FEs[xeid[0]].vs_link[0] != FrameI.FVs[ECs[i].bvs[0]].index_hex)
				reverse_vector(ids_vx);

			set_cross(FrameI.FVs[ECs[i].bvs[0]].neighbor_Es, FrameI.FVs[ECs[i].bvs_valence_ring[0][(j + 1) % N_fan]].neighbor_Es, xeid);
			Y_dim = FrameI.FEs[xeid[0]].vs_link.size() - 1;

			ids_vy = FrameI.FEs[xeid[0]].vs_link;
			if (FrameI.FEs[xeid[0]].vs_link[0] != FrameI.FVs[ECs[i].bvs[0]].index_hex)
				reverse_vector(ids_vy);

			for (int k = 0; k < ECs[i].bes.size(); k++)
			{
				Z_dims.push_back(FrameI.FEs[ECs[i].bes[k]].vs_link.size() - 1);
				ids_vz = FrameI.FEs[ECs[i].bes[k]].vs_link;
				if (FrameI.FEs[ECs[i].bes[k]].vs_link[0] != FrameI.FVs[ECs[i].bvs[k]].index_hex)
					reverse_vector(ids_vz);
				append_vector(ids_vzz, ids_vz);
			}
			set_redundent_clearn(ids_vzz);

			len = 0.0;
			fan_chain fanc;
			for (int k = 0; k < Z_dims.size(); k++)
			{
				for (int p = 0; p < Z_dims[k] + 1; p++)
				{
					if (k != Z_dims.size() - 1 && p == Z_dims[k])
						continue;

					para_coor vz0, vz1, vz2, vz3;
					vector<float> origin_0;
					origin_0.push_back(0);
					origin_0.push_back(0);
					origin_0.push_back(len);
					vector<float> origin_1;
					origin_1.push_back(0);
					origin_1.push_back(0);
					origin_1.push_back(len + 1);
					interpolation_vector(origin_0, origin_1, vz0.coords, float(p) / Z_dims[k]);
					interpolation_vector(ECs[i].ring1_paracoords[k][2 * j].coords, ECs[i].ring1_paracoords[k + 1][2 * j].coords, vz1.coords,
							float(p) / Z_dims[k]);
					interpolation_vector(ECs[i].ring1_paracoords[k][2 * j + 1].coords, ECs[i].ring1_paracoords[k + 1][2 * j + 1].coords, vz2.coords,
							float(p) / Z_dims[k]);
					interpolation_vector(ECs[i].ring1_paracoords[k][(2 * j + 2) % (N_fan * 2)].coords,
							ECs[i].ring1_paracoords[k + 1][(2 * j + 2) % (N_fan * 2)].coords, vz3.coords, float(p) / Z_dims[k]);

					vector<vector<para_coor>> paracoords_xy;
					for (int n = 0; n < Y_dim + 1; n++)
					{
						para_coor vy0, vy1;
						interpolation_vector(vz0.coords, vz3.coords, vy0.coords, float(n) / Y_dim);
						interpolation_vector(vz1.coords, vz2.coords, vy1.coords, float(n) / Y_dim);

						vector<para_coor> paracoords_x;
						for (int m = 0; m < X_dim + 1; m++)
						{
							para_coor vx0;
							interpolation_vector(vy0.coords, vy1.coords, vx0.coords, float(m) / X_dim);
							paracoords_x.push_back(vx0);
						}
						paracoords_xy.push_back(paracoords_x);
					}
					fanc.UVW_coords.push_back(paracoords_xy);
				}
			}

			//assign all hex-v index
			vector<bool> arrayhv_ids;
			initializeVectorT(arrayhv_ids, true, hex_mesh.HVs.size());
			for (int m = 0; m < ECs[i].bhs.size(); m++)
			{
				for (int n = 0; n < FrameI.FHs[ECs[i].bhs[m][j]].hs_net.size(); n++)
				{
					int hid = FrameI.FHs[ECs[i].bhs[m][j]].hs_net[n];
					for (int p = 0; p < 8; p++)
						arrayhv_ids[hex_mesh.HHs[hid].V_Ids[p]] = false;
				}
			}

			vector<vector<int>> ids_xy;
			ids_xy.push_back(ids_vx);
			for (int n = 0; n < ids_vx.size(); n++)
				arrayhv_ids[ids_vx[n]] = true;
			for (int n = 0; n < ids_vy.size(); n++)
				arrayhv_ids[ids_vy[n]] = true;
			for (int n = 1; n < Y_dim + 1; n++)
			{
				vector<int> ids_vx_temp;
				ids_vx_temp.push_back(ids_vy[n]);
				for (int m = 1; m < X_dim + 1; m++)
				{
					vector<int> cross_vs;
					set_cross(hex_mesh.HVs[ids_vx_temp[m - 1]].neighbor_vs, hex_mesh.HVs[ids_xy[n - 1][m]].neighbor_vs, cross_vs);
					for (int q = 0; q < cross_vs.size(); q++)
					{
						if (!arrayhv_ids[cross_vs[q]])
						{
							ids_vx_temp.push_back(cross_vs[q]);
							arrayhv_ids[cross_vs[q]] = true;
						}
					}
				}
				ids_xy.push_back(ids_vx_temp);
			}
			fanc.hv_ids.push_back(ids_xy);

			for (int p = 1; p < ids_vzz.size(); p++)
			{
				ids_xy.clear();
				for (int n = 0; n < Y_dim + 1; n++)
				{
					vector<int> ids_vx_temp;
					for (int m = 0; m < X_dim + 1; m++)
					{
						if (m == 0 && n == 0)
						{
							ids_vx_temp.push_back(ids_vzz[p]);
							arrayhv_ids[ids_vzz[p]] = true;
						}
						else if (n == 0)
						{
							vector<int> cross_vs;
							set_cross(hex_mesh.HVs[ids_vx_temp[m - 1]].neighbor_vs, hex_mesh.HVs[fanc.hv_ids[p - 1][n][m]].neighbor_vs, cross_vs);
							for (int q = 0; q < cross_vs.size(); q++)
							{
								if (!arrayhv_ids[cross_vs[q]])
								{
									ids_vx_temp.push_back(cross_vs[q]);
									arrayhv_ids[cross_vs[q]] = true;
								}
							}
						}
						else
						{
							vector<int> cross_vs;
							set_cross(hex_mesh.HVs[ids_xy[n - 1][m]].neighbor_vs, hex_mesh.HVs[fanc.hv_ids[p - 1][n][m]].neighbor_vs, cross_vs);
							for (int q = 0; q < cross_vs.size(); q++)
							{
								if (!arrayhv_ids[cross_vs[q]])
								{
									ids_vx_temp.push_back(cross_vs[q]);
									arrayhv_ids[cross_vs[q]] = true;
								}
							}
						}
					}
					ids_xy.push_back(ids_vx_temp);
				}
				fanc.hv_ids.push_back(ids_xy);
			}
			ParaC.fan_cicle.push_back(fanc);
		}
		ParaCs.push_back(ParaC);
	}
}
void parametric_optimization::parameter_index_correspondence_new()
{
	for (int i = 0; i < ECs.size(); i++)
	{
		vector<int> Z_dims, ids_vz, ids_vzz;
		for (int k = 0; k < ECs[i].bes.size(); k++)
		{
			Z_dims.push_back(FrameI.FEs[ECs[i].bes[k]].vs_link.size() - 1);
			ids_vz = FrameI.FEs[ECs[i].bes[k]].vs_link;
			if (FrameI.FEs[ECs[i].bes[k]].vs_link[0] != FrameI.FVs[ECs[i].bvs[k]].index_hex)
				reverse_vector(ids_vz);
			append_vector(ids_vzz, ids_vz);
		}
		set_redundent_clearn(ids_vzz);

		float kai = 4.0 / ECs[i].bvs_valence_ring[0].size();
		//assign all parameter coords
		//printf("This is %dth Ec\n",i);
		Parameterization_Chain ParaC;
		int N_fan = ECs[i].bvs_valence_ring[0].size();
		for (int j = 0; j < N_fan; j++)
		{
			//printf("This is %dth Fan\n",j);
			//determine x, y, z dimensions
			int X_dim, Y_dim;
			vector<int> xeid;
			vector<int> ids_vx, ids_vy;

			set_cross(FrameI.FVs[ECs[i].bvs[0]].neighbor_Es, FrameI.FVs[ECs[i].bvs_valence_ring[0][j]].neighbor_Es, xeid);
			X_dim = FrameI.FEs[xeid[0]].vs_link.size() - 1;

			ids_vx = FrameI.FEs[xeid[0]].vs_link;
			if (FrameI.FEs[xeid[0]].vs_link[0] != FrameI.FVs[ECs[i].bvs[0]].index_hex)
				reverse_vector(ids_vx);

			set_cross(FrameI.FVs[ECs[i].bvs[0]].neighbor_Es, FrameI.FVs[ECs[i].bvs_valence_ring[0][(j + 1) % N_fan]].neighbor_Es, xeid);
			Y_dim = FrameI.FEs[xeid[0]].vs_link.size() - 1;

			ids_vy = FrameI.FEs[xeid[0]].vs_link;
			if (FrameI.FEs[xeid[0]].vs_link[0] != FrameI.FVs[ECs[i].bvs[0]].index_hex)
				reverse_vector(ids_vy);

			float len = 0.0;
			fan_chain fanc;
			for (int k = 0; k < Z_dims.size(); k++)
			{
				for (int p = 0; p < Z_dims[k] + 1; p++)
				{
					if (k != Z_dims.size() - 1 && p == Z_dims[k])
						continue;
					vector<vector<para_coor>> paracoords_xy;
					for (int n = 0; n < Y_dim + 1; n++)
					{
						vector<para_coor> paracoords_x;
						for (int m = 0; m < X_dim + 1; m++)
						{
							para_coor vx0;
							vector<float> polar_coords;
							float rou = sqrt(double(m * m + n * n));
							float rou_powed = pow(rou, kai);
							polar_coords.push_back(rou_powed);

							float alpha = atan2(double(n), double(m));
							float alpha_exponetial = alpha * kai + j * 2 * PI / ECs[i].bvs_valence_ring[0].size();
							polar_coords.push_back(alpha_exponetial);

							vx0.coords.push_back(polar_coords[0] * cos(polar_coords[1]));
							vx0.coords.push_back(polar_coords[0] * sin(polar_coords[1]));
							vx0.coords.push_back(len);
							paracoords_x.push_back(vx0);
						}
						paracoords_xy.push_back(paracoords_x);
					}
					fanc.UVW_coords.push_back(paracoords_xy);
					len++;
				}
			}

			//assign all hex-v index
			vector<bool> arrayhv_ids;
			initializeVectorT(arrayhv_ids, true, hex_mesh.HVs.size());
			for (int m = 0; m < ECs[i].bhs.size(); m++)
			{
				for (int n = 0; n < FrameI.FHs[ECs[i].bhs[m][j]].hs_net.size(); n++)
				{
					int hid = FrameI.FHs[ECs[i].bhs[m][j]].hs_net[n];
					for (int p = 0; p < 8; p++)
						arrayhv_ids[hex_mesh.HHs[hid].V_Ids[p]] = false;
				}
			}

			vector<vector<int>> ids_xy;
			ids_xy.push_back(ids_vx);
			for (int n = 0; n < ids_vx.size(); n++)
				arrayhv_ids[ids_vx[n]] = true;
			for (int n = 0; n < ids_vy.size(); n++)
				arrayhv_ids[ids_vy[n]] = true;
			for (int n = 1; n < Y_dim + 1; n++)
			{
				vector<int> ids_vx_temp;
				ids_vx_temp.push_back(ids_vy[n]);
				for (int m = 1; m < X_dim + 1; m++)
				{
					vector<int> cross_vs;
					set_cross(hex_mesh.HVs[ids_vx_temp[m - 1]].neighbor_vs, hex_mesh.HVs[ids_xy[n - 1][m]].neighbor_vs, cross_vs);
					for (int q = 0; q < cross_vs.size(); q++)
					{
						if (!arrayhv_ids[cross_vs[q]])
						{
							ids_vx_temp.push_back(cross_vs[q]);
							arrayhv_ids[cross_vs[q]] = true;
						}
					}
				}
				ids_xy.push_back(ids_vx_temp);
			}
			fanc.hv_ids.push_back(ids_xy);

			for (int p = 1; p < ids_vzz.size(); p++)
			{
				ids_xy.clear();
				for (int n = 0; n < Y_dim + 1; n++)
				{
					vector<int> ids_vx_temp;
					for (int m = 0; m < X_dim + 1; m++)
					{
						if (m == 0 && n == 0)
						{
							ids_vx_temp.push_back(ids_vzz[p]);
							arrayhv_ids[ids_vzz[p]] = true;
						}
						else if (n == 0)
						{
							vector<int> cross_vs;
							set_cross(hex_mesh.HVs[ids_vx_temp[m - 1]].neighbor_vs, hex_mesh.HVs[fanc.hv_ids[p - 1][n][m]].neighbor_vs, cross_vs);
							for (int q = 0; q < cross_vs.size(); q++)
							{
								if (!arrayhv_ids[cross_vs[q]])
								{
									ids_vx_temp.push_back(cross_vs[q]);
									arrayhv_ids[cross_vs[q]] = true;
								}
							}
						}
						else
						{
							vector<int> cross_vs;
							set_cross(hex_mesh.HVs[ids_xy[n - 1][m]].neighbor_vs, hex_mesh.HVs[fanc.hv_ids[p - 1][n][m]].neighbor_vs, cross_vs);
							for (int q = 0; q < cross_vs.size(); q++)
							{
								if (!arrayhv_ids[cross_vs[q]])
								{
									ids_vx_temp.push_back(cross_vs[q]);
									arrayhv_ids[cross_vs[q]] = true;
								}
							}
						}
					}
					ids_xy.push_back(ids_vx_temp);
				}
				fanc.hv_ids.push_back(ids_xy);
			}
			ParaC.fan_cicle.push_back(fanc);
		}
		ParaCs.push_back(ParaC);
	}
}
void parametric_optimization::find_ecs_cur(vector<int> &ecs_ids, vector<bool> &ecs_Ind)
{
	ecs_ids.clear();
	for (int i = 0; i < ECs.size(); i++)
	{
		if (ecs_Ind[i])
		{
			bool no_inter = true;
			for (int j = 0; j < ecs_ids.size(); j++)
			{
				vector<int> vectors1, vectors2, vectors3, vectors4;
				append_vector(vectors1, ECs[i].boundary_fs);
				for (int k = 0; k < ECs[i].bfs.size(); k++)
					append_vector(vectors1, ECs[i].bfs[k]);
				append_vector(vectors2, ECs[ecs_ids[j]].boundary_fs);
				for (int k = 0; k < ECs[ecs_ids[j]].bfs.size(); k++)
					append_vector(vectors2, ECs[ecs_ids[j]].bfs[k]);
				set_exclusion(vectors1, vectors2, vectors3);
				if (vectors1.size() != vectors3.size())
				{
					no_inter = false;
					break;
				}
			}
			if (no_inter)
			{
				ecs_ids.push_back(i);
				ecs_Ind[i] = false;
			}
		}
	}
}
void parametric_optimization::a_chain_parameterization(int EC_ID)
{
	File_Index = EC_ID;

	h_io io;

	vector<Hex> hs;
	for (int j = 0; j < ECs[EC_ID].bhs.size(); j++)
	{
		for (int k = 0; k < ECs[EC_ID].bhs[j].size(); k++)
		{
			int Fhid = ECs[EC_ID].bhs[j][k];
			for (int m = 0; m < FrameI.FHs[Fhid].hs_net.size(); m++)
			{
				int hid = FrameI.FHs[Fhid].hs_net[m];
				Hex h = hex_mesh.HHs[hid];
				hs.push_back(h);
			}
		}
	}
	char fname1[300];
	sprintf(fname1, "%s%d%s", "C:/Users/Xifeng_Gao/Desktop/temps/x", EC_ID, "_before.mesh");
	//io.write_hex_mesh_mesh(hex_mesh.HVs,hs,fname1);

	vector<int> Vs_arrayIds;
	vector<Hex_V> Vs_patch;
	vector<Hex_E> Es_patch;
	vector<Hex_F> Fs_patch;
	vector<Hex> Hs_patch;
	reindexing_cuboid_mesh(EC_ID, Vs_arrayIds, Vs_patch, Es_patch, Fs_patch, Hs_patch);

	vector<Vertex> VTs_patch;
	vector<Edge> ETs_patch;
	vector<Triangle> FTs_patch;
	vector<Tet> TETs_patch;
	translate_to_tetmesh(Vs_arrayIds, Vs_patch, Es_patch, Fs_patch, Hs_patch, VTs_patch, ETs_patch, FTs_patch, TETs_patch);

	char fname[300];
	sprintf(fname, "%s%d%s", "C:/Users/Xifeng_Gao/Desktop/temps/cuboid_", EC_ID, ".obj");
	//io.write_tet_mesh_obj(VTs_patch,TETs_patch,"C:/Users/Xifeng_Gao/Desktop/temps/cuboid.obj");
	//io.write_triangle_mesh_obj(VTs_patch,FTs_patch,fname);

	Parameterization_Cell pc;
	for (int i = 0; i < VTs_patch.size(); i++)
	{
		vector<float> uvw;
		uvw.push_back(0);
		uvw.push_back(0);
		uvw.push_back(0);
		pc.UVW_coords.push_back(uvw);
	}

	//point
	//printf("POINT:::: \n");
	for (int i = 0; i < ParaCs[EC_ID].fan_cicle.size(); i++)
	{
		int m = ParaCs[EC_ID].fan_cicle[i].UVW_coords[0].size();
		int n = ParaCs[EC_ID].fan_cicle[i].UVW_coords[0][0].size();
		int size_z = ParaCs[EC_ID].fan_cicle[i].UVW_coords.size();
		int size_layer = ECs[EC_ID].bvs_valence_ring.size();
		for (int j = 0; j < 2; j++)
		{
			para_coor ppc = ParaCs[EC_ID].fan_cicle[i].UVW_coords[0][0][n - 1];
			int vid = ECs[EC_ID].bvs_valence_ring[j][i];
			if (j == 1)
			{
				ppc = ParaCs[EC_ID].fan_cicle[i].UVW_coords[size_z - 1][0][n - 1];
				vid = ECs[EC_ID].bvs_valence_ring[size_layer - 1][i];
			}
			vid = FrameI.FVs[vid].index_hex;
			pc.UVW_coords[Vs_arrayIds[vid]] = ppc.coords;

			ppc = ParaCs[EC_ID].fan_cicle[i].UVW_coords[0][m - 1][n - 1];
			vid = ECs[EC_ID].bvs_1ring[j][2 * i + 1];
			if (j == 1)
			{
				ppc = ParaCs[EC_ID].fan_cicle[i].UVW_coords[size_z - 1][m - 1][n - 1];
				vid = ECs[EC_ID].bvs_1ring[size_layer - 1][2 * i + 1];
			}
			vid = FrameI.FVs[vid].index_hex;
			pc.UVW_coords[Vs_arrayIds[vid]] = ppc.coords;
			//printf(" %f %f %f\n",ppc.coords[0],ppc.coords[1],ppc.coords[2]);
		}
	}

// 	int hex_vid=FrameI.FVs[ECs[EC_ID].bvs[0]].index_hex;
// 	//if(hex_mesh.HVs[hex_vid].where_location!=1)
// 	{
// 		//printf("add this.");
// 		vector<float> ppcc;ppcc.push_back(0);
// 		ppcc.push_back(0);ppcc.push_back(0);
// 		pc.UVW_coords[Vs_arrayIds[hex_vid]]=ppcc;
// 	}
// 	hex_vid=FrameI.FVs[ECs[EC_ID].bvs[ECs[EC_ID].bvs.size()-1]].index_hex;
// 	//if(hex_mesh.HVs[hex_vid].where_location!=1)
// 	{
// 		//printf("add this.");
// 		vector<float> ppcc;ppcc.push_back(0);
// 		ppcc.push_back(0);ppcc.push_back(ECs[EC_ID].bvs.size()-1);
// 		pc.UVW_coords[Vs_arrayIds[hex_vid]]=ppcc;
// 	}
	//curve
	//printf("CURVE:::: \n");
	int size_layer = ECs[EC_ID].bvs_1ring.size();
	for (int i = 0; i < 2; i++)
	{
		float kai = 4.0 / ECs[EC_ID].bvs_valence_ring[0].size();

		int size_zlen = 0;
		if (i == 1)
		{
			i = size_layer - 1;
			size_zlen = ParaCs[EC_ID].fan_cicle[0].UVW_coords.size() - 1;
		}
		for (int j = 0; j < ECs[EC_ID].bvs_valence_ring[i].size(); j++)
		{
			int X_dim, Y_dim;
			vector<int> Z_dims;
			vector<int> xeid;
			vector<int> cur_evs[2];
			vector<int> ids_vzz;

			set_cross(FrameI.FVs[ECs[EC_ID].bvs[0]].neighbor_Es, FrameI.FVs[ECs[EC_ID].bvs_valence_ring[0][j]].neighbor_Es, xeid);
			X_dim = FrameI.FEs[xeid[0]].vs_link.size() - 1;

			set_cross(FrameI.FVs[ECs[EC_ID].bvs[0]].neighbor_Es,
					FrameI.FVs[ECs[EC_ID].bvs_valence_ring[0][(j + 1) % ECs[EC_ID].bvs_valence_ring[i].size()]].neighbor_Es, xeid);
			Y_dim = FrameI.FEs[xeid[0]].vs_link.size() - 1;

			set_cross(FrameI.FVs[ECs[EC_ID].bvs_1ring[i][2 * j]].neighbor_Es, FrameI.FVs[ECs[EC_ID].bvs_1ring[i][2 * j + 1]].neighbor_Es, xeid);
			cur_evs[0] = FrameI.FEs[xeid[0]].vs_link;
			if (cur_evs[0][0] != FrameI.FVs[ECs[EC_ID].bvs_1ring[i][2 * j]].index_hex)
				reverse_vector(cur_evs[0]);

			set_cross(FrameI.FVs[ECs[EC_ID].bvs_1ring[i][2 * j + 1]].neighbor_Es,
					FrameI.FVs[ECs[EC_ID].bvs_1ring[i][(2 * j + 2) % ECs[EC_ID].bvs_1ring[0].size()]].neighbor_Es, xeid);
			cur_evs[1] = FrameI.FEs[xeid[0]].vs_link;
			if (cur_evs[1][0] != FrameI.FVs[ECs[EC_ID].bvs_1ring[i][2 * j + 1]].index_hex)
				reverse_vector(cur_evs[1]);

			for (int k = 0; k < 2; k++)
			{
				float total_len = 0;
				vector<float> ratios;
				for (int m = 0; m < cur_evs[k].size() - 1; m++)
				{
					float dis;
					DISTANCE(dis, hex_mesh.HVs[cur_evs[k][m]].v, hex_mesh.HVs[cur_evs[k][m + 1]].v);
					total_len += dis;
					ratios.push_back(total_len);
				}

				vector<float> startuvw, enduvw;
				startuvw.push_back(X_dim);
				startuvw.push_back(0);
				startuvw.push_back(size_zlen);
				enduvw.push_back(X_dim);
				enduvw.push_back(Y_dim);
				enduvw.push_back(size_zlen);
				if (k == 1)
				{
					startuvw = enduvw;
					enduvw.clear();
					enduvw.push_back(0);
					enduvw.push_back(Y_dim);
					enduvw.push_back(size_zlen);
				}
				for (int m = 1; m < cur_evs[k].size() - 1; m++)
				{
					vector<float> uvw;
					uvw.push_back(startuvw[0] + ratios[m - 1] / total_len * (enduvw[0] - startuvw[0]));
					uvw.push_back(startuvw[1] + ratios[m - 1] / total_len * (enduvw[1] - startuvw[1]));
					uvw.push_back(startuvw[2] + ratios[m - 1] / total_len * (enduvw[2] - startuvw[2]));

					para_coor vx0;
					vector<float> polar_coords;
					float rou = sqrt(uvw[0] * uvw[0] + uvw[1] * uvw[1]);
					float rou_powed = pow(rou, kai);
					polar_coords.push_back(rou_powed);

					float alpha = atan2(double(uvw[1]), double(uvw[0])) * kai + j * 2 * PI / ECs[EC_ID].bvs_valence_ring[0].size();
					polar_coords.push_back(alpha);

					vx0.coords.push_back(polar_coords[0] * cos(polar_coords[1]));
					vx0.coords.push_back(polar_coords[0] * sin(polar_coords[1]));
					vx0.coords.push_back(uvw[2]);

					pc.UVW_coords[Vs_arrayIds[cur_evs[k][m]]] = vx0.coords;
					//printf(" %f %f %f\n",uvw[0],uvw[1],uvw[2]);
				}
			}
		}
	}

	for (int j = 0; j < ECs[EC_ID].bvs_1ring[0].size(); j++)
	{
		float kai = 4.0 / ECs[EC_ID].bvs_valence_ring[0].size();

		vector<int> ids_vzz;
		for (int i = 0; i < ECs[EC_ID].bvs_1ring.size() - 1; i++)
		{
			int v1 = ECs[EC_ID].bvs_1ring[i][j], v2 = ECs[EC_ID].bvs_1ring[i + 1][j];

			vector<int> cur_evs, xeid;
			set_cross(FrameI.FVs[v1].neighbor_Es, FrameI.FVs[v2].neighbor_Es, xeid);
			cur_evs = FrameI.FEs[xeid[0]].vs_link;
			if (cur_evs[0] != FrameI.FVs[v1].index_hex)
				reverse_vector(cur_evs);
			append_vector(ids_vzz, cur_evs);
		}

		float total_len = 0;
		vector<float> ratios;
		for (int m = 0; m < ids_vzz.size() - 1; m++)
		{
			float dis;
			DISTANCE(dis, hex_mesh.HVs[ids_vzz[m]].v, hex_mesh.HVs[ids_vzz[m + 1]].v);
			total_len += dis;
			ratios.push_back(total_len);
		}

		int m = ParaCs[EC_ID].fan_cicle[j / 2].UVW_coords[0].size();
		int n = ParaCs[EC_ID].fan_cicle[j / 2].UVW_coords[0][0].size();

		vector<float> startuvw = ParaCs[EC_ID].fan_cicle[j / 2].UVW_coords[0][0][n - 1].coords, enduvw =
				ParaCs[EC_ID].fan_cicle[j / 2].UVW_coords[ids_vzz.size() - 1][0][n - 1].coords;
		if (j % 2 == 1)
		{
			startuvw = ParaCs[EC_ID].fan_cicle[j / 2].UVW_coords[0][m - 1][n - 1].coords, enduvw =
					ParaCs[EC_ID].fan_cicle[j / 2].UVW_coords[ids_vzz.size() - 1][m - 1][n - 1].coords;
		}
		for (int m = 1; m < ids_vzz.size() - 1; m++)
		{
			vector<float> uvw;
			uvw.push_back(startuvw[0] + ratios[m - 1] / total_len * (enduvw[0] - startuvw[0]));
			uvw.push_back(startuvw[1] + ratios[m - 1] / total_len * (enduvw[1] - startuvw[1]));
			uvw.push_back(startuvw[2] + ratios[m - 1] / total_len * (enduvw[2] - startuvw[2]));

			pc.UVW_coords[Vs_arrayIds[ids_vzz[m]]] = uvw;
			//printf(" %f %f %f\n",uvw[0],uvw[1],uvw[2]);
		}
	}

	//face
	//printf("FACE:::: \n");
	D2_parameterization(0, pc, VTs_patch, ETs_patch, FTs_patch);
	//volume
	//	clock_t start_time=clock();
	//watch->startTimer();
	//printf("VOLUME:::: \n");
	D3_parameterization(pc, VTs_patch, ETs_patch, FTs_patch, TETs_patch);
	//watch->stopTimer();
	// 	clock_t end_time=clock();
	// 	std::cout<< "Running time is: "<<static_cast<double>(end_time-start_time)/CLOCKS_PER_SEC*1000<<"ms"<<std::endl;

	//printf("Parameteric of component %d UVW\n",EC_ID);
	parametric_coords(EC_ID, Vs_arrayIds, pc, VTs_patch, FTs_patch, TETs_patch);
}
void parametric_optimization::reindexing_cuboid_mesh(int EC_ID, vector<int> &Vs_arrayIds, vector<Hex_V> &vs, vector<Hex_E> &es, vector<Hex_F> &fs,
		vector<Hex> &hs)
{
	initializeVectorT(Vs_arrayIds, -2, hex_mesh.HVs.size());
	for (int i = 0; i < ECs[EC_ID].bhs.size(); i++)
	{
		int Color_id = 0;
		for (int j = 0; j < ECs[EC_ID].bhs[i].size(); j++)
		{
			Color_id++;
			int Fhid = ECs[EC_ID].bhs[i][j];

			for (int m = 0; m < FrameI.FHs[Fhid].hs_net.size(); m++)
			{
				int hid = FrameI.FHs[Fhid].hs_net[m];
				for (int k = 0; k < 8; k++)
				{
					Vs_arrayIds[hex_mesh.HHs[hid].V_Ids[k]] = -1;
				}
				Hex h;
				for (int k = 0; k < 8; k++)
					h.V_Ids[k] = hex_mesh.HHs[hid].V_Ids[k];
				h.index = hs.size();
				h.Color_ID = Color_id;
				hs.push_back(h);
			}
		}
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

	//on base-complex vs
	int hex_vid = FrameI.FVs[ECs[EC_ID].bvs[0]].index_hex;
	//if(hex_mesh.HVs[hex_vid].where_location!=1)
	//vs[Vs_arrayIds[hex_vid]].is_on_base_complex=1;;//
	hex_vid = FrameI.FVs[ECs[EC_ID].bvs[ECs[EC_ID].bvs.size() - 1]].index_hex;
	//if(hex_mesh.HVs[hex_vid].where_location!=1)
	//vs[Vs_arrayIds[hex_vid]].is_on_base_complex=1;;//

	vector<int> ringes = ECs[EC_ID].boundaryes_excepttb;
	for (int i = 0; i < ringes.size(); i++)
	{
		for (int j = 0; j < FrameI.FEs[ringes[i]].vs_link.size(); j++)
		{
			int vid = FrameI.FEs[ringes[i]].vs_link[j];
			vs[Vs_arrayIds[vid]].is_on_base_complex = 1;
		}
	}
	//on boundary faces vs
	for (int i = 0; i < fs.size(); i++)
		fs[i].which_F_face = -1;
	for (int i = 0; i < ECs[EC_ID].boundary_fs.size(); i++)
	{
		int Ffid = ECs[EC_ID].boundary_fs[i];
		// 	for(int i=0;i<ECs[EC_ID].all_fs.size();i++)
		// 	{
		// 		int Ffid=ECs[EC_ID].all_fs[i];

		for (int j = 0; j < FrameI.FFs[Ffid].hfs_net_another.size(); j++)
		{
			int fid = FrameI.FFs[Ffid].hfs_net_another[j];
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
			fs[sharedf[0]].which_F_face = 0;

			vs[v1].which_F_face = 0;
			vs[v2].which_F_face = 0;
			vs[v3].which_F_face = 0;
			vs[v4].which_F_face = 0;
		}
	}
}
void parametric_optimization::parametric_coords(int EC_ID, vector<int> &Vs_arrayIds, Parameterization_Cell &pc, vector<Vertex> &vts,
		vector<Triangle> &fts, vector<Tet> &tets)
{
	vector<bool> array_hv_ids;
	initializeVectorT(array_hv_ids, false, hex_mesh.HVs.size());
	for (int i = 0; i < ParaCs[EC_ID].fan_cicle.size(); i++)
	{
		for (int j = 0; j < ParaCs[EC_ID].fan_cicle[i].hv_ids.size(); j++)
			for (int k = 0; k < ParaCs[EC_ID].fan_cicle[i].hv_ids[j].size(); k++)
				for (int m = 0; m < ParaCs[EC_ID].fan_cicle[i].hv_ids[j][k].size(); m++)
				{
					int hid = ParaCs[EC_ID].fan_cicle[i].hv_ids[j][k][m];
					array_hv_ids[hid] = true;
				}
	}
// 	if(ECs[EC_ID].bvs_valence_ring[0].size()!=4)
// 	{
// 		int v1=ECs[EC_ID].startendIds[0],v2=ECs[EC_ID].startendIds[1];
// 		int hv1=FrameI.FVs[v1].index_hex,hv2=FrameI.FVs[v2].index_hex;
// 		if(hex_mesh.HVs[hv1].where_location==1&&FrameI.FVs[v1].neighbor_Hs.size()!=4)
// 			array_hv_ids[hv1]=false;
// 		if(hex_mesh.HVs[hv2].where_location==1&&FrameI.FVs[v2].neighbor_Hs.size()!=4)
// 			array_hv_ids[hv2]=false;
// 	}

	int N_fan = ParaCs[EC_ID].fan_cicle.size();
	vector<vector<vector<para_coor>>> Fan_UVW_coords;
	for (int i = 0; i < N_fan; i++)
	{
		Fan_UVW_coords.push_back(ParaCs[EC_ID].fan_cicle[i].UVW_coords[0]);
	}
	int Z_dim = ECs[EC_ID].bes.size();
	vector<float> Z_dims;
	float len = 0;
	for (int i = 0; i < Z_dim; i++)
	{
		int beid = ECs[EC_ID].bes[i];
		int cursize = FrameI.FEs[beid].vs_link.size() - 1;
		for (int j = 0; j < cursize; j++)
		{
			Z_dims.push_back(len++);
		}
	}
	Z_dims.push_back(len);

	float cur_EPS = 100 * EPS;
	for (int i = 0; i < tets.size(); i++)
	{
		vector<vector<float>> vers;
		vers.push_back(pc.UVW_coords[tets[i].vs[0]]);
		vers.push_back(pc.UVW_coords[tets[i].vs[1]]);
		vers.push_back(pc.UVW_coords[tets[i].vs[2]]);
		vers.push_back(pc.UVW_coords[tets[i].vs[3]]);

		//min max
		float min_u = 10000, max_u = -10000, min_v = 10000, max_v = -10000, min_w = Z_dims[Z_dims.size() - 1], max_w = 0;
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
		min_u -= cur_EPS;
		min_v -= cur_EPS;
		min_w -= cur_EPS;
		max_u += cur_EPS;
		max_v += cur_EPS;
		max_w += cur_EPS;

		vector<vector<int>> candidates;
		vector<vector<float>> candidates_coords;

		for (int j = 0; j < N_fan; j++)
		{
			vector<float> z_dims, z_dimsint;
			for (int n = 0; n < Z_dims.size(); n++)
			{
				if (Z_dims[n] >= min_w && Z_dims[n] <= max_w)
				{
					z_dims.push_back(Z_dims[n]);
					z_dimsint.push_back(n);
				}
			}
			if (!z_dims.size())
				break;
			vector<vector<int>> candidates_temp;
			for (int k = 0; k < Fan_UVW_coords[j].size(); k++)
			{
				for (int m = 0; m < Fan_UVW_coords[j][k].size(); m++)
				{
					if (Fan_UVW_coords[j][k][m].coords[0] >= min_u && Fan_UVW_coords[j][k][m].coords[0] <= max_u)
					{
						if (Fan_UVW_coords[j][k][m].coords[1] >= min_v && Fan_UVW_coords[j][k][m].coords[1] <= max_v)
						{
							vector<int> onep;
							onep.push_back(k);
							onep.push_back(m);
							candidates_temp.push_back(onep);
						}
					}
				}
			}
			for (int k = 0; k < candidates_temp.size(); k++)
			{
				for (int m = 0; m < z_dims.size(); m++)
				{
					vector<int> onep = candidates_temp[k];
					onep.push_back(z_dimsint[m]);
					onep.push_back(j);
					int hid = ParaCs[EC_ID].fan_cicle[j].hv_ids[onep[2]][onep[0]][onep[1]];
					if (array_hv_ids[hid])
					{
						candidates.push_back(onep);
						vector<float> onepc;
						onepc = Fan_UVW_coords[j][onep[0]][onep[1]].coords;
						onepc[2] = z_dims[m];
						candidates_coords.push_back(onepc);
					}
				}
			}
		}

		if (!candidates.size())
			continue;

		bool tet_coplanar = false;
		float deter0;
		vector<float> ver1 = vers[0], ver2 = vers[1], ver3 = vers[2], ver4 = vers[3];
		deter0 = matrix_determinant(ver1, ver2, ver3, ver4);
		if (abs(deter0) < EPS)
		{
			printf("Tet %d is a coplanar-four point set\n", i);
			tet_coplanar = true;
			continue;;		//tet is a coplanar-four point set
		}

		for (int j = 0; j < candidates_coords.size(); j++)
		{
			int hvid = ParaCs[EC_ID].fan_cicle[candidates[j][3]].hv_ids[candidates[j][2]][candidates[j][0]][candidates[j][1]];
			if (!array_hv_ids[hvid])
				continue;

			vector<float> lamdas;
			vector<float> deters;
			if (!tet_coplanar)
			{
				bool next = false;

				vector<float> cur_uvw;
				cur_uvw.push_back(candidates_coords[j][0]);
				cur_uvw.push_back(candidates_coords[j][1]);
				cur_uvw.push_back(candidates_coords[j][2]);
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
			else
			{
				for (int n = 0; n < 4; n++)
					lamdas.push_back(0.25);
			}

			if (!isfinite(lamdas[0]) || !isfinite(lamdas[1]) || !isfinite(lamdas[2]) || !isfinite(lamdas[3]))
			{
				printf("volume %d", EC_ID);
				continue;
			}

			hex_mesh.HVs[hvid].v[0] = hex_mesh.HVs[hvid].v[1] = hex_mesh.HVs[hvid].v[2] = 0;
			for (int n = 0; n < 4; n++)
			{
				hex_mesh.HVs[hvid].v[0] += lamdas[n] * vts[tets[i].vs[n]].v[0];
				hex_mesh.HVs[hvid].v[1] += lamdas[n] * vts[tets[i].vs[n]].v[1];
				hex_mesh.HVs[hvid].v[2] += lamdas[n] * vts[tets[i].vs[n]].v[2];
			}

			vector<float> illv;
			illv.push_back(-0.003533);
			illv.push_back(-0.001745);
			illv.push_back(-0.003753);
			float dis = -1;
			DISTANCE(dis, illv, hex_mesh.HVs[hvid].v);
			if (abs(dis) < EPS)
			{
				printf("EC_ID %d, %f %f %f; ", EC_ID, hex_mesh.HVs[hvid].v[0], hex_mesh.HVs[hvid].v[1], hex_mesh.HVs[hvid].v[2]);
				printf("%f %f %f %f\n", lamdas[0], lamdas[1], lamdas[2], lamdas[3]);
			}
			array_hv_ids[hvid] = false;
		}
	}

	laplacian_smoothing laps;
	for (int j = 0; j < N_fan; j++)
	{
		for (int k = 0; k < Fan_UVW_coords[j].size(); k++)
		{
			for (int m = 0; m < Fan_UVW_coords[j][k].size(); m++)
			{
				for (int n = 0; n < Z_dims.size(); n++)
				{
					int hid = ParaCs[EC_ID].fan_cicle[j].hv_ids[n][k][m];
					if (array_hv_ids[hid])
					{
						int thid = Vs_arrayIds[hid];
						vector<float> uvw_cur = ParaCs[EC_ID].fan_cicle[j].UVW_coords[n][k][m].coords;

						if (vts[thid].on_base_complex)
						{
							vector<int> cur_linevs = ParaCs[EC_ID].fan_cicle[j].hv_ids[n][k];
							if ((m == 0 || m == Fan_UVW_coords[j][k].size() - 1) && (k != 0 && k != Fan_UVW_coords[j].size() - 1))
							{
								cur_linevs.clear();
								for (int p = 0; p < Fan_UVW_coords[j].size(); p++)
									cur_linevs.push_back(ParaCs[EC_ID].fan_cicle[j].hv_ids[n][p][m]);
							}
							else if ((m == 0 || m == Fan_UVW_coords[j][k].size() - 1) && (k == 0 || k == Fan_UVW_coords[j].size() - 1))
							{
								cur_linevs.clear();
								for (int p = 0; p < Z_dims.size(); p++)
									cur_linevs.push_back(ParaCs[EC_ID].fan_cicle[j].hv_ids[p][k][m]);
							}
							for (int p = 0; p < cur_linevs.size() - 1; p++)
							{
								int nv1 = Vs_arrayIds[cur_linevs[p]], nv2 = Vs_arrayIds[cur_linevs[p + 1]];
								float kk = -1;
								point_line_projection(pc.UVW_coords[nv1], pc.UVW_coords[nv2], uvw_cur, kk);
								if (kk >= 0 && kk <= 1)
								{
									if (isfinite(kk))
									{
										continue;
										printf("base-complex curve kk %d", EC_ID);
									}
									for (int q = 0; q < 3; q++)
									{
										hex_mesh.HVs[hid].v[q] = vts[nv1].v[q] + kk * (vts[nv2].v[q] - vts[nv1].v[q]);
									}

									vector<float> illv;
									illv.push_back(-0.003533);
									illv.push_back(-0.001745);
									illv.push_back(-0.003753);
									float dis = -1;
									DISTANCE(dis, illv, hex_mesh.HVs[hid].v);
									if (abs(dis) < EPS)
									{
										printf("EC_ID %d, %f %f %f; %f\n", EC_ID, hex_mesh.HVs[hid].v[0], hex_mesh.HVs[hid].v[1],
												hex_mesh.HVs[hid].v[2], kk);
									}

									array_hv_ids[hid] = false;
									break;
								}
							}
						}
						else
						{
							float Min_dis = 10000;
							int which_one = -1;

							for (int p = 0; p < vts.size(); p++)
							{
								if (vts[p].which_F_face != 0)
									continue;
								float dis;
								DISTANCE(dis, pc.UVW_coords[p], uvw_cur);
								if (Min_dis > dis)
								{
									Min_dis = dis;
									which_one = p;
								}
							}
							vector<int> vidnTs = vts[which_one].neighbort;
							for (int p = 0; p < vts[which_one].neighborv.size(); p++)
							{
								int nvid = vts[which_one].neighborv[p];
								append_vector(vidnTs, vts[nvid].neighbort);
							}
							set_redundent_clearn(vidnTs);
							std::vector<Vertex> projected_vs;
							std::vector<float> Diss;
							Min_dis = 100000;
							int which_one_inT = -1, whichtid = -1;
							for (int p = 0; p < vidnTs.size(); p++)
							{
								if (fts[vidnTs[p]].which_F_face != 0)
									continue;

								float x, y, z;
								float dis;
								bool istrue = laps.projected_v(pc.UVW_coords[fts[vidnTs[p]].triangle_v[0]],
										pc.UVW_coords[fts[vidnTs[p]].triangle_v[1]], pc.UVW_coords[fts[vidnTs[p]].triangle_v[2]], uvw_cur, &x, &y, &z,
										&dis);
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
										whichtid = p;
									}
								}
							}
							if (which_one_inT != -1)
							{
								vector<float> weights;
								triangle_coordinates(pc.UVW_coords[fts[vidnTs[whichtid]].triangle_v[0]],
										pc.UVW_coords[fts[vidnTs[whichtid]].triangle_v[1]], pc.UVW_coords[fts[vidnTs[whichtid]].triangle_v[2]],
										uvw_cur, weights);
								if (isfinite(weights[0]) && isfinite(weights[1]) && isfinite(weights[2]))
								{
									hex_mesh.HVs[hid].v[0] = vts[fts[vidnTs[whichtid]].triangle_v[0]].v[0] * weights[0]
											+ vts[fts[vidnTs[whichtid]].triangle_v[1]].v[0] * weights[1]
											+ vts[fts[vidnTs[whichtid]].triangle_v[2]].v[0] * weights[2];
									hex_mesh.HVs[hid].v[1] = vts[fts[vidnTs[whichtid]].triangle_v[0]].v[1] * weights[0]
											+ vts[fts[vidnTs[whichtid]].triangle_v[1]].v[1] * weights[1]
											+ vts[fts[vidnTs[whichtid]].triangle_v[2]].v[1] * weights[2];
									hex_mesh.HVs[hid].v[2] = vts[fts[vidnTs[whichtid]].triangle_v[0]].v[2] * weights[0]
											+ vts[fts[vidnTs[whichtid]].triangle_v[1]].v[2] * weights[1]
											+ vts[fts[vidnTs[whichtid]].triangle_v[2]].v[2] * weights[2];
									array_hv_ids[hid] = false;

									vector<float> illv;
									illv.push_back(-0.003533);
									illv.push_back(-0.001745);
									illv.push_back(-0.003753);
									float dis = -1;
									DISTANCE(dis, illv, hex_mesh.HVs[hid].v);
									if (abs(dis) < EPS)
									{
										printf("EC_ID %d, %f %f %f; %f %f %f\n", EC_ID, hex_mesh.HVs[hid].v[0], hex_mesh.HVs[hid].v[1],
												hex_mesh.HVs[hid].v[2], weights[0], weights[1], weights[2]);
									}

								}
							}
							{
								int which_one_v = -1;
								vector<int> vidnVs = vts[which_one].neighborv;
								for (int p = 0; p < vidnVs.size(); p++)
								{
									if (vts[vidnVs[p]].which_F_face != 0)
										continue;
									float kk = -1;
									point_line_projection(pc.UVW_coords[vidnVs[p]], pc.UVW_coords[which_one], uvw_cur, kk);
									if (kk >= 0 && kk <= 1)
									{
										if (!isfinite(kk))
											continue;
										Vertex v, v_true;
										for (int q = 0; q < 3; q++)
										{
											v.v[q] = pc.UVW_coords[vidnVs[p]][q] + kk * (pc.UVW_coords[which_one][q] - pc.UVW_coords[vidnVs[p]][q]);
											v_true.v[q] = vts[vidnVs[p]].v[q] + kk * (vts[which_one].v[q] - vts[vidnVs[p]].v[q]);
										}
										float dis;
										DISTANCE(dis, uvw_cur, v.v);
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
									hex_mesh.HVs[hid].v[0] = projected_vs[which_one_v].v[0];
									hex_mesh.HVs[hid].v[1] = projected_vs[which_one_v].v[1];
									hex_mesh.HVs[hid].v[2] = projected_vs[which_one_v].v[2];
									array_hv_ids[hid] = false;

									vector<float> illv;
									illv.push_back(-0.003533);
									illv.push_back(-0.001745);
									illv.push_back(-0.003753);
									float dis = -1;
									DISTANCE(dis, illv, hex_mesh.HVs[hid].v);
									if (abs(dis) < EPS)
									{
										printf("EC_ID %d, %f %f %f; %d\n", EC_ID, hex_mesh.HVs[hid].v[0], hex_mesh.HVs[hid].v[1],
												hex_mesh.HVs[hid].v[2], which_one_v);
									}
								}
							}
						}
					}
				}
			}
		}
	}
	for (int i = 0; i < array_hv_ids.size(); i++)
	{
		if (array_hv_ids[i])
			printf("hereerer i %d\n", EC_ID);
	}

	h_io io;
	vector<Hex> hs;
	for (int j = 0; j < ECs[EC_ID].bhs.size(); j++)
	{
		int Color_id = 0;
		for (int k = 0; k < ECs[EC_ID].bhs[j].size(); k++)
		{
			Color_id++;
			int Fhid = ECs[EC_ID].bhs[j][k];
			for (int m = 0; m < FrameI.FHs[Fhid].hs_net.size(); m++)
			{
				int hid = FrameI.FHs[Fhid].hs_net[m];
				Hex h = hex_mesh.HHs[hid];
				h.Color_ID = Color_id;
				hs.push_back(h);
			}
		}
	}
}

//face
void parametric_optimization::face_center_extraction()
{
	for (int i = 0; i < FrameI.FFs.size(); i++)
	{
		face_center fc;
		fc.id = FCs.size();
		fc.f_id = i;
		//vs_all, es_all, fs_boundary
		vector<int> fs_all;
		for (int j = 0; j < FrameI.FFs[i].neighbor_Cs.size(); j++)
		{
			int cid = FrameI.FFs[i].neighbor_Cs[j];
			for (int k = 0; k < 8; k++)
			{
				fc.vs_all.push_back(FrameI.FHs[cid].FV_Ids[k]);
			}
			for (int k = 0; k < 6; k++)
			{
				fs_all.push_back(FrameI.FHs[cid].neighbor_FS[k]);
				for (int m = 0; m < 4; m++)
					fc.es_all.push_back(FrameI.FFs[FrameI.FHs[cid].neighbor_FS[k]].fe_Ids[m]);
			}
		}
		set_redundent_clearn(fc.vs_all);
		set_redundent_clearn(fc.es_all);
		set_redundent_clearn(fc.fs_boundary);
		vector<int> f_temp;
		f_temp.push_back(i);
		set_exclusion(fs_all, f_temp, fc.fs_boundary);
		//conditions for FCi
		if (FrameI.FFs[i].neighbor_Cs.size() == 1)
			;		//continue;
		int boundarynum = 0;
		for (int j = 0; j < fs_all.size(); j++)
		{
			if (FrameI.FFs[fs_all[j]].neighbor_Cs.size() == 1)
				boundarynum++;
		}
		if (boundarynum != 1)
			;		//continue;

		//fs_paralell
		for (int j = 0; j < fc.fs_boundary.size(); j++)
		{
			vector<int> fvs, fcur_vs;
			for (int k = 0; k < 4; k++)
			{
				fvs.push_back(FrameI.FFs[fc.fs_boundary[j]].fv_Ids[k]);
				fcur_vs.push_back(FrameI.FFs[i].fv_Ids[k]);
			}
			if (set_contain(fvs, fcur_vs, 0))		//not touching
			{
				fc.fs_paralell.push_back(fc.fs_boundary[j]);
				if (fc.fs_paralell.size() == 1)
					fc.fs_paralell.push_back(i);
			}
		}
		//vs_paralell es_paralell
		vector<int> vs_layer;
		for (int j = 0; j < fc.fs_paralell.size(); j++)
		{
			if (j == 0)
			{
				for (int k = 0; k < 4; k++)
					vs_layer.push_back(FrameI.FFs[fc.fs_paralell[j]].fv_Ids[k]);
			}
			else
			{
				vector<int> vs_curf;
				for (int k = 0; k < 4; k++)
					vs_curf.push_back(FrameI.FFs[fc.fs_paralell[j]].fv_Ids[k]);
				for (int k = 0; k < 4; k++)
				{
					vector<int> vs_cross;
					set_cross(FrameI.FVs[fc.vs_paralell[j - 1][k]].neighbor_vs, vs_curf, vs_cross);
					vs_layer.push_back(vs_cross[0]);
				}
			}
			fc.vs_paralell.push_back(vs_layer);
			vs_layer.clear();
		}
		for (int j = 0; j < fc.vs_paralell.size(); j++)
		{
			vector<int> es_layer;
			for (int k = 0; k < fc.vs_paralell[j].size(); k++)
			{
				int v1 = fc.vs_paralell[j][k], v2 = fc.vs_paralell[j][(k + 1) % 4];
				vector<int> eid;
				set_cross(FrameI.FVs[v1].neighbor_Es, FrameI.FVs[v2].neighbor_Es, eid);
				es_layer.push_back(eid[0]);
			}
			fc.es_paralell.push_back(es_layer);
		}
		FCs.push_back(fc);
	}
}
void parametric_optimization::face_parameter_index_correspondence()
{
	float para[4][2] =
	{
	{ 0, 0 },
	{ 1, 0 },
	{ 1, 1 },
	{ 0, 1 } };
	for (int i = 0; i < FCs.size(); i++)
	{
		for (int j = 0; j < FCs[i].vs_paralell.size(); j++)
		{
			vector<para_coor> layer_paracoords;
			for (int k = 0; k < 4; k++)
			{
				para_coor pc;
				pc.coords.push_back(para[k][0]);
				pc.coords.push_back(para[k][1]);
				layer_paracoords.push_back(pc);
			}
			FCs[i].layer_paracoords.push_back(layer_paracoords);
		}

		fan_chain fanc;
		//determine x, y, z dimensions
		int X_dim, Y_dim;
		vector<int> Z_dims;
		vector<int> xeid;
		vector<int> ids_vx, ids_vy, ids_vz;
		vector<int> ids_vzz;

		X_dim = FrameI.FEs[FCs[i].es_paralell[0][0]].vs_link.size() - 1;

		ids_vx = FrameI.FEs[FCs[i].es_paralell[0][0]].vs_link;
		if (FrameI.FEs[FCs[i].es_paralell[0][0]].vs_link[0] != FrameI.FVs[FCs[i].vs_paralell[0][0]].index_hex)
			reverse_vector(ids_vx);

		Y_dim = FrameI.FEs[FCs[i].es_paralell[0][3]].vs_link.size() - 1;

		ids_vy = FrameI.FEs[FCs[i].es_paralell[0][3]].vs_link;
		if (FrameI.FEs[FCs[i].es_paralell[0][3]].vs_link[0] != FrameI.FVs[FCs[i].vs_paralell[0][0]].index_hex)
			reverse_vector(ids_vy);

		for (int j = 0; j < FCs[i].vs_paralell.size() - 1; j++)
		{
			int v1 = FCs[i].vs_paralell[j][0], v2 = FCs[i].vs_paralell[j + 1][0];
			vector<int> eid;
			set_cross(FrameI.FVs[v1].neighbor_Es, FrameI.FVs[v2].neighbor_Es, eid);

			Z_dims.push_back(FrameI.FEs[eid[0]].vs_link.size() - 1);
			ids_vz = FrameI.FEs[eid[0]].vs_link;
			if (FrameI.FEs[eid[0]].vs_link[0] != FrameI.FVs[v1].index_hex)
				reverse_vector(ids_vz);
			append_vector(ids_vzz, ids_vz);
		}
		set_redundent_clearn(ids_vzz);

		for (int p = 0; p < 4; p++)
		{
			FCs[i].layer_paracoords[0][p].coords.push_back(0);
		}
		bool HALF_ORNOT = false;
		vector<float> z_dim_para;
		int count_dim = 1;
		for (int k = 0; k < Z_dims.size(); k++)
		{
			count_dim--;
			for (int p = 0; p < Z_dims[k] + 1; p++)
			{
				if (HALF_ORNOT)
					z_dim_para.push_back(k + float(p) / Z_dims[k]);
				else
					z_dim_para.push_back(count_dim++);
			}
			for (int p = 0; p < 4; p++)
			{
				FCs[i].layer_paracoords[k + 1][p].coords.push_back(z_dim_para[z_dim_para.size() - 1]);
			}
		}
		count_dim = 0;
		for (int k = 0; k < Z_dims.size(); k++)
		{
			for (int p = 0; p < Z_dims[k] + 1; p++)
			{
				count_dim++;
				if (k != Z_dims.size() - 1 && p == Z_dims[k])
					continue;

				para_coor vz0, vz1, vz2, vz3;
				vz0 = FCs[i].layer_paracoords[0][0];
				vz0.coords.push_back(z_dim_para[count_dim - 1]);
				vz1 = FCs[i].layer_paracoords[0][1];
				vz1.coords.push_back(z_dim_para[count_dim - 1]);
				vz2 = FCs[i].layer_paracoords[0][2];
				vz2.coords.push_back(z_dim_para[count_dim - 1]);
				vz3 = FCs[i].layer_paracoords[0][3];
				vz3.coords.push_back(z_dim_para[count_dim - 1]);

				vector<vector<para_coor>> paracoords_xy;
				for (int n = 0; n < Y_dim + 1; n++)
				{
					para_coor vy0, vy1;
					interpolation_vector(vz0.coords, vz3.coords, vy0.coords, float(n) / Y_dim);
					interpolation_vector(vz1.coords, vz2.coords, vy1.coords, float(n) / Y_dim);

					vector<para_coor> paracoords_x;
					for (int m = 0; m < X_dim + 1; m++)
					{
						para_coor vx0;
						interpolation_vector(vy0.coords, vy1.coords, vx0.coords, float(m) / X_dim);
						paracoords_x.push_back(vx0);
					}
					paracoords_xy.push_back(paracoords_x);
				}
				fanc.UVW_coords.push_back(paracoords_xy);
			}
		}

		//assign all hex-v index
		vector<bool> arrayhv_ids;
		initializeVectorT(arrayhv_ids, true, hex_mesh.HVs.size());
		for (int m = 0; m < FrameI.FFs[FCs[i].f_id].neighbor_Cs.size(); m++)
		{
			int fhid = FrameI.FFs[FCs[i].f_id].neighbor_Cs[m];
			for (int n = 0; n < FrameI.FHs[fhid].hs_net.size(); n++)
			{
				int hid = FrameI.FHs[fhid].hs_net[n];
				for (int p = 0; p < 8; p++)
					arrayhv_ids[hex_mesh.HHs[hid].V_Ids[p]] = false;
			}
		}

		vector<vector<int>> ids_xy;
		ids_xy.push_back(ids_vx);
		for (int n = 0; n < ids_vx.size(); n++)
			arrayhv_ids[ids_vx[n]] = true;
		for (int n = 0; n < ids_vy.size(); n++)
			arrayhv_ids[ids_vy[n]] = true;
		for (int n = 1; n < Y_dim + 1; n++)
		{
			vector<int> ids_vx_temp;
			ids_vx_temp.push_back(ids_vy[n]);
			for (int m = 1; m < X_dim + 1; m++)
			{
				vector<int> cross_vs;
				set_cross(hex_mesh.HVs[ids_vx_temp[m - 1]].neighbor_vs, hex_mesh.HVs[ids_xy[n - 1][m]].neighbor_vs, cross_vs);
				for (int q = 0; q < cross_vs.size(); q++)
				{
					if (!arrayhv_ids[cross_vs[q]])
					{
						ids_vx_temp.push_back(cross_vs[q]);
						arrayhv_ids[cross_vs[q]] = true;
					}
				}
			}
			ids_xy.push_back(ids_vx_temp);
		}
		fanc.hv_ids.push_back(ids_xy);

		for (int p = 1; p < ids_vzz.size(); p++)
		{
			ids_xy.clear();
			for (int n = 0; n < Y_dim + 1; n++)
			{
				vector<int> ids_vx_temp;
				for (int m = 0; m < X_dim + 1; m++)
				{
					if (m == 0 && n == 0)
					{
						ids_vx_temp.push_back(ids_vzz[p]);
						arrayhv_ids[ids_vzz[p]] = true;
					}
					else if (n == 0)
					{
						vector<int> cross_vs;
						set_cross(hex_mesh.HVs[ids_vx_temp[m - 1]].neighbor_vs, hex_mesh.HVs[fanc.hv_ids[p - 1][n][m]].neighbor_vs, cross_vs);
						for (int q = 0; q < cross_vs.size(); q++)
						{
							if (!arrayhv_ids[cross_vs[q]])
							{
								ids_vx_temp.push_back(cross_vs[q]);
								arrayhv_ids[cross_vs[q]] = true;
							}
						}
					}
					else
					{
						vector<int> cross_vs;
						set_cross(hex_mesh.HVs[ids_xy[n - 1][m]].neighbor_vs, hex_mesh.HVs[fanc.hv_ids[p - 1][n][m]].neighbor_vs, cross_vs);
						for (int q = 0; q < cross_vs.size(); q++)
						{
							if (!arrayhv_ids[cross_vs[q]])
							{
								ids_vx_temp.push_back(cross_vs[q]);
								arrayhv_ids[cross_vs[q]] = true;
							}
						}
					}
				}
				ids_xy.push_back(ids_vx_temp);
			}
			fanc.hv_ids.push_back(ids_xy);
		}
		ParaFCs.push_back(fanc);
	}
}

void parametric_optimization::find_fcs_cur(vector<int> &fcs_ids, vector<bool> &fcs_Ind)
{
	fcs_ids.clear();
	for (int i = 0; i < FCs.size(); i++)
	{
		if (fcs_Ind[i])
		{
			bool no_inter = true;
			for (int j = 0; j < fcs_ids.size(); j++)
			{
				vector<int> excfs;
				set_exclusion(FCs[i].fs_boundary, FCs[fcs_ids[j]].fs_boundary, excfs);
				if (FCs[i].fs_boundary.size() != excfs.size())
				{
					no_inter = false;
					break;
				}
			}
			if (no_inter)
			{
				fcs_ids.push_back(i);
				fcs_Ind[i] = false;
			}
		}
	}
}
void parametric_optimization::face_center_parameterization(int FC_ID)
{
	File_Index = FC_ID;

	h_io io;

	vector<int> Vs_arrayIds;
	vector<Hex_V> Vs_patch;
	vector<Hex_E> Es_patch;
	vector<Hex_F> Fs_patch;
	vector<Hex> Hs_patch;
	face_reindexing_cuboid_mesh(FC_ID, Vs_arrayIds, Vs_patch, Es_patch, Fs_patch, Hs_patch);

	vector<Vertex> VTs_patch;
	vector<Edge> ETs_patch;
	vector<Triangle> FTs_patch;
	vector<Tet> TETs_patch;
	translate_to_tetmesh(Vs_arrayIds, Vs_patch, Es_patch, Fs_patch, Hs_patch, VTs_patch, ETs_patch, FTs_patch, TETs_patch);

	char fname[300];
	sprintf(fname, "%s%d%s", "C:/Users/Xifeng_Gao/Desktop/temps/cuboid_", FC_ID, ".obj");
	//io.write_tet_mesh_obj(VTs_patch,TETs_patch,"C:/Users/Xifeng_Gao/Desktop/temps/cuboid.obj");
	//io.write_triangle_mesh_obj(VTs_patch,FTs_patch,fname);

	Parameterization_Cell pc;
	for (int i = 0; i < VTs_patch.size(); i++)
	{
		vector<float> uvw;
		uvw.push_back(0);
		uvw.push_back(0);
		uvw.push_back(0);
		pc.UVW_coords.push_back(uvw);
	}

	//point
	//printf("POINT:::: \n");
	for (int i = 0; i < FCs[FC_ID].layer_paracoords.size(); i++)
	{
		for (int j = 0; j < FCs[FC_ID].layer_paracoords[i].size(); j++)
		{
			para_coor ppc = FCs[FC_ID].layer_paracoords[i][j];
			int vid = FCs[FC_ID].vs_paralell[i][j];
			vid = FrameI.FVs[vid].index_hex;
			pc.UVW_coords[Vs_arrayIds[vid]] = ppc.coords;
			//printf(" %f %f %f\n",ppc.coords[0],ppc.coords[1],ppc.coords[2]);
		}
	}
	//curve
	//printf("CURVE:::: \n");
	vector<int> ringes = FCs[FC_ID].es_all;
	if (FCs[FC_ID].es_paralell.size() == 3)
	{
		vector<int> bedges_temp;
		set_exclusion(ringes, FCs[FC_ID].es_paralell[1], bedges_temp);
		;	//ringes=bedges_temp;
	}
// 	vector<int> ringes=FCs[FC_ID].es_paralell[0];
// 	append_vector(ringes,FCs[FC_ID].es_paralell[FCs[FC_ID].es_paralell.size()-1]);
	for (int i = 0; i < ringes.size(); i++)
	{
		float total_len = 0;
		vector<float> ratios;
		for (int j = 0; j < FrameI.FEs[ringes[i]].vs_link.size() - 1; j++)
		{
			int v1, v2;
			v1 = FrameI.FEs[ringes[i]].vs_link[j];
			v2 = FrameI.FEs[ringes[i]].vs_link[j + 1];
			float dis;
			DISTANCE(dis, hex_mesh.HVs[v1].v, hex_mesh.HVs[v2].v);
			total_len += dis;
			ratios.push_back(total_len);
		}
		int sv_id = FrameI.FEs[ringes[i]].vs_link[0], ev_id = FrameI.FEs[ringes[i]].vs_link[FrameI.FEs[ringes[i]].vs_link.size() - 1];
		vector<float> startuvw = pc.UVW_coords[Vs_arrayIds[sv_id]], enduvw = pc.UVW_coords[Vs_arrayIds[ev_id]];

		for (int j = 1; j < FrameI.FEs[ringes[i]].vs_link.size() - 1; j++)
		{
			vector<float> uvw;
			uvw.push_back(startuvw[0] + ratios[j - 1] / total_len * (enduvw[0] - startuvw[0]));
			uvw.push_back(startuvw[1] + ratios[j - 1] / total_len * (enduvw[1] - startuvw[1]));
			uvw.push_back(startuvw[2] + ratios[j - 1] / total_len * (enduvw[2] - startuvw[2]));
			pc.UVW_coords[Vs_arrayIds[FrameI.FEs[ringes[i]].vs_link[j]]] = uvw;
			//printf(" %f %f %f\n",uvw[0],uvw[1],uvw[2]);
		}
	}
	if (0)
	//if(FCs[FC_ID].es_paralell.size()==3)
	{
		int z_len = ParaFCs[FC_ID].hv_ids.size();
		vector<int> corners[2];
		corners[0].push_back(0);
		corners[0].push_back(0);
		corners[0].push_back(ParaFCs[FC_ID].hv_ids[0].size() - 1);
		corners[0].push_back(ParaFCs[FC_ID].hv_ids[0].size() - 1);
		corners[1].push_back(0);
		corners[1].push_back(ParaFCs[FC_ID].hv_ids[0][0].size() - 1);
		corners[1].push_back(0);
		corners[1].push_back(ParaFCs[FC_ID].hv_ids[0][0].size() - 1);

		for (int j = 0; j < 4; j++)
		{
			float total_len = 0;
			vector<float> ratios;
			for (int i = 0; i < z_len - 1; i++)
			{
				int v1, v2;
				v1 = ParaFCs[FC_ID].hv_ids[i][corners[0][j]][corners[1][j]];
				v2 = ParaFCs[FC_ID].hv_ids[i + 1][corners[0][j]][corners[1][j]];
				float dis;
				DISTANCE(dis, hex_mesh.HVs[v1].v, hex_mesh.HVs[v2].v);
				total_len += dis;
				ratios.push_back(total_len);
			}
			int sv_id = ParaFCs[FC_ID].hv_ids[0][corners[0][j]][corners[1][j]],
					ev_id = ParaFCs[FC_ID].hv_ids[z_len - 1][corners[0][j]][corners[1][j]];
			;
			vector<float> startuvw = pc.UVW_coords[Vs_arrayIds[sv_id]], enduvw = pc.UVW_coords[Vs_arrayIds[ev_id]];

			for (int i = 1; i < z_len - 1; i++)
			{
				vector<float> uvw;
				uvw.push_back(startuvw[0] + ratios[i - 1] / total_len * (enduvw[0] - startuvw[0]));
				uvw.push_back(startuvw[1] + ratios[i - 1] / total_len * (enduvw[1] - startuvw[1]));
				uvw.push_back(startuvw[2] + ratios[i - 1] / total_len * (enduvw[2] - startuvw[2]));

				int v1 = ParaFCs[FC_ID].hv_ids[i][corners[0][j]][corners[1][j]];
				pc.UVW_coords[Vs_arrayIds[v1]] = uvw;
				//printf(" %f %f %f\n",uvw[0],uvw[1],uvw[2]);
			}
		}
	}

	//face
	//printf("FACE:::: \n");
	D2_parameterization(0, pc, VTs_patch, ETs_patch, FTs_patch);
	//volume
	//	clock_t start_time=clock();
	//watch->startTimer();
	//printf("VOLUME:::: \n");
	D3_parameterization(pc, VTs_patch, ETs_patch, FTs_patch, TETs_patch);
	//watch->stopTimer();
	// 	clock_t end_time=clock();
	// 	std::cout<< "Running time is: "<<static_cast<double>(end_time-start_time)/CLOCKS_PER_SEC*1000<<"ms"<<std::endl;

	//laplacian_smoothing ls;
	//vector<Vertex> VTs_patch_para;
	//for(int i=0;i<VTs_patch.size();i++)
	//{
	//	VTs_patch_para.push_back(VTs_patch[i]);
	//	VTs_patch_para[i].v[0]=pc.UVW_coords[i][0];
	//	VTs_patch_para[i].v[1]=pc.UVW_coords[i][1];
	//	VTs_patch_para[i].v[2]=pc.UVW_coords[i][2];
	//}
	//ls.laplacian_global_pipeline_tet(VTs_patch,FTs_patch,TETs_patch);
	//for(int i=0;i<VTs_patch.size();i++)
	//{
	//	pc.UVW_coords[i][0]=VTs_patch_para[i].v[0];
	//	pc.UVW_coords[i][1]=VTs_patch_para[i].v[1];
	//	pc.UVW_coords[i][2]=VTs_patch_para[i].v[2];
	//}
	//printf("Parameteric of component %d UVW\n",FC_ID);
	face_parametric_coords(FC_ID, pc, VTs_patch, TETs_patch);
}
void parametric_optimization::face_reindexing_cuboid_mesh(int FC_ID, vector<int> &Vs_arrayIds, vector<Hex_V> &vs, vector<Hex_E> &es,
		vector<Hex_F> &fs, vector<Hex> &hs)
{
	int Color_id = 0;
	initializeVectorT(Vs_arrayIds, -2, hex_mesh.HVs.size());
	for (int i = 0; i < FrameI.FFs[FCs[FC_ID].f_id].neighbor_Cs.size(); i++)
	{
		Color_id++;
		int Fhid = FrameI.FFs[FCs[FC_ID].f_id].neighbor_Cs[i];
		for (int m = 0; m < FrameI.FHs[Fhid].hs_net.size(); m++)
		{
			int hid = FrameI.FHs[Fhid].hs_net[m];
			for (int k = 0; k < 8; k++)
			{
				Vs_arrayIds[hex_mesh.HHs[hid].V_Ids[k]] = -1;
			}
			Hex h;
			for (int k = 0; k < 8; k++)
				h.V_Ids[k] = hex_mesh.HHs[hid].V_Ids[k];
			h.index = hs.size();
			h.Color_ID = Color_id;
			hs.push_back(h);
		}
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

	vector<int> bedges = FCs[FC_ID].es_all;
	if (FCs[FC_ID].es_paralell.size() == 3)
	{
		vector<int> bedges_temp;
		set_exclusion(bedges, FCs[FC_ID].es_paralell[1], bedges_temp);
		//bedges=bedges_temp;
		;
	}
	for (int i = 0; i < bedges.size(); i++)
	{
		for (int j = 0; j < FrameI.FEs[bedges[i]].vs_link.size(); j++)
		{
			int vid = FrameI.FEs[bedges[i]].vs_link[j];
			vs[Vs_arrayIds[vid]].is_on_base_complex = 1;
		}
	}
	//on boundary faces vs
	for (int i = 0; i < fs.size(); i++)
		fs[i].which_F_face = -1;
	vector<int> fss = FCs[FC_ID].fs_boundary;
	if (FCs[FC_ID].fs_boundary.size() == 5)
		fss.push_back(FCs[FC_ID].f_id);
	for (int i = 0; i < fss.size(); i++)
	{
		int Ffid = fss[i];
		for (int j = 0; j < FrameI.FFs[Ffid].hfs_net_another.size(); j++)
		{
			int fid = FrameI.FFs[Ffid].hfs_net_another[j];
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
			fs[sharedf[0]].which_F_face = 0;

			vs[v1].which_F_face = 0;
			vs[v2].which_F_face = 0;
			vs[v3].which_F_face = 0;
			vs[v4].which_F_face = 0;
		}
	}
}

void parametric_optimization::translate_to_tetmesh(vector<int> &Vs_arrayIds, vector<Hex_V> &vs, vector<Hex_E> &es, vector<Hex_F> &fs, vector<Hex> &hs,
		vector<Vertex> &vts, vector<Edge> &ets, vector<Triangle> &fts, vector<Tet> &tets)
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
				t.Color_Id = hs[i].Color_ID;
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
				t.which_F_face = 0;
				for (int k = 0; k < 3; k++)
				{
					if (vts[t.triangle_v[k]].which_F_face == -1)
						t.which_F_face = -1;
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
}
struct Lamda
{
	vector<float> lamda_ks;
};
bool parametric_optimization::D2_parameterization(int F_Id, Parameterization_Cell &pc, vector<Vertex> &vts, vector<Edge> &ets, vector<Triangle> &fts)
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
		std::vector<float> weights, weights_temp;
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
			if (divide != divide)
			{
				float v0[3], v1[3];
				v0[0] = vts_patch[i].v[0];
				v0[1] = vts_patch[i].v[1];
				v0[2] = vts_patch[i].v[2];
				v1[0] = vts_patch[vts_patch[i].neighborv[(j + 1) % vts_patch[i].neighborv.size()]].v[0];
				v1[1] = vts_patch[vts_patch[i].neighborv[(j + 1) % vts_patch[i].neighborv.size()]].v[1];
				v1[2] = vts_patch[vts_patch[i].neighborv[(j + 1) % vts_patch[i].neighborv.size()]].v[2];
				if (disL == 0)
				{
					v1[0] = vts_patch[vts_patch[i].neighborv[j]].v[0];
					v1[1] = vts_patch[vts_patch[i].neighborv[j]].v[1];
					v1[2] = vts_patch[vts_patch[i].neighborv[j]].v[2];
					printf("probably QNAN: v0:%d-%f %f %f; v1:%d-%f %f %f\n", i, v0[0], v0[1], v0[2], vts_patch[i].neighborv[j], v1[0], v1[1], v1[2]);
				}
				else
					printf("probably QNAN: v0:%d-%f %f %f; v1:%d-%f %f %f\n", i, v0[0], v0[1], v0[2],
							vts_patch[i].neighborv[(j + 1) % vts_patch[i].neighborv.size()], v1[0], v1[1], v1[2]);
				divide = 1;
			}
			if (divide < -1)
				divide = -1;
			if (divide > 1)
				divide = 1;
			float angle_ = acos(divide);
			angle.push_back(angle_);
			weights_temp.push_back(disL);
		}

		float w_total = 0;
		weights = weights_temp;
		for (int j = 0; j < angle.size(); j++)
		{
			weights[j] = (tan(angle[(j - 1 + angle.size()) % angle.size()] / 2) + tan(angle[j] / 2)) / weights_temp[j];
			if (weights[j] != weights[j])
			{
				printf("probably QNAN\n");
				weights[j] = 1;
			}
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

	h_io io;
	vector<Vertex> tempvers;
	for (int i = 0; i < vts.size(); i++)
	{
		Vertex v;
		v.v[0] = pc.UVW_coords[i][0];
		v.v[1] = pc.UVW_coords[i][1];
		v.v[2] = pc.UVW_coords[i][2];
		tempvers.push_back(v);
	}
	char fname[300];
	sprintf(fname, "%s%d%s", "C:/Users/Xifeng_Gao/Desktop/temps/para_surface", File_Index, ".obj");
	//io.write_triangle_mesh_obj(tempvers,fts,fname);

	return true;
}
bool parametric_optimization::D3_parameterization(Parameterization_Cell &pc, vector<Vertex> &vts, vector<Edge> &ets, vector<Triangle> &fts,
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
					if (ws[ws.size() - 1] != ws[ws.size() - 1])	//self added, judge if NAN value
						ws[ws.size() - 1] = 1;

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
				if (wis[wis.size() - 1] != wis[wis.size() - 1])	//self added, judge if NAN value
					wis[wis.size() - 1] = 1;

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

	h_io io;
	vector<Vertex> tempvers;
	for (int i = 0; i < vts.size(); i++)
	{
		Vertex v;
		v.v[0] = pc.UVW_coords[i][0];
		v.v[1] = pc.UVW_coords[i][1];
		v.v[2] = pc.UVW_coords[i][2];
		tempvers.push_back(v);
	}
	return true;
}
void parametric_optimization::face_parametric_coords(int FC_ID, Parameterization_Cell &pc, vector<Vertex> &vts, vector<Tet> &tets)
{
	vector<bool> array_hv_ids;
	initializeVectorT(array_hv_ids, false, hex_mesh.HVs.size());
	for (int i = 0; i < ParaFCs[FC_ID].hv_ids.size(); i++)
	{
		for (int k = 0; k < ParaFCs[FC_ID].hv_ids[i].size(); k++)
			for (int m = 0; m < ParaFCs[FC_ID].hv_ids[i][k].size(); m++)
			{
				int hid = ParaFCs[FC_ID].hv_ids[i][k][m];
				array_hv_ids[hid] = true;
			}
	}

	vector<int> Z_dims_int;
	for (int j = 0; j < FCs[FC_ID].vs_paralell.size() - 1; j++)
	{
		int v1 = FCs[FC_ID].vs_paralell[j][0], v2 = FCs[FC_ID].vs_paralell[j + 1][0];
		vector<int> eid;
		set_cross(FrameI.FVs[v1].neighbor_Es, FrameI.FVs[v2].neighbor_Es, eid);
		Z_dims_int.push_back(FrameI.FEs[eid[0]].vs_link.size() - 1);
	}

	bool HALF_ORNOT = false;
	int count_dim = 0;
	vector<float> Z_dims;
	for (int i = 0; i < Z_dims_int.size(); i++)
	{
		int cursize = Z_dims_int[i];
		for (int j = 0; j < cursize; j++)
		{
			if (HALF_ORNOT)
				Z_dims.push_back(i + (float(j) / cursize));
			else
				Z_dims.push_back(count_dim++);
		}
	}
	if (HALF_ORNOT)
		Z_dims.push_back(Z_dims_int.size());
	else
		Z_dims.push_back(count_dim);

	float cur_EPS = 100 * EPS;
	for (int i = 0; i < tets.size(); i++)
	{
		vector<vector<float>> vers;
		vers.push_back(pc.UVW_coords[tets[i].vs[0]]);
		vers.push_back(pc.UVW_coords[tets[i].vs[1]]);
		vers.push_back(pc.UVW_coords[tets[i].vs[2]]);
		vers.push_back(pc.UVW_coords[tets[i].vs[3]]);

		//min max
		float min_u = 1, max_u = 0, min_v = 1, max_v = 0, min_w = Z_dims[Z_dims.size() - 1], max_w = 0;
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
		min_u -= cur_EPS;
		min_v -= cur_EPS;
		min_w -= cur_EPS;
		max_u += cur_EPS;
		max_v += cur_EPS;
		max_w += cur_EPS;

		vector<vector<int>> candidates;
		vector<vector<float>> candidates_coords;

		vector<float> z_dims, z_dimsint;
		for (int n = 0; n < Z_dims.size(); n++)
		{
			if (Z_dims[n] >= min_w && Z_dims[n] <= max_w)
			{
				z_dims.push_back(Z_dims[n]);
				z_dimsint.push_back(n);
			}
		}
		if (!z_dims.size())
			continue;

		vector<vector<int>> candidates_temp;
		for (int k = 0; k < ParaFCs[FC_ID].UVW_coords[0].size(); k++)
		{
			for (int m = 0; m < ParaFCs[FC_ID].UVW_coords[0][k].size(); m++)
			{
				if (ParaFCs[FC_ID].UVW_coords[0][k][m].coords[0] >= min_u && ParaFCs[FC_ID].UVW_coords[0][k][m].coords[0] <= max_u)
				{
					if (ParaFCs[FC_ID].UVW_coords[0][k][m].coords[1] >= min_v && ParaFCs[FC_ID].UVW_coords[0][k][m].coords[1] <= max_v)
					{
						vector<int> onep;
						onep.push_back(k);
						onep.push_back(m);
						candidates_temp.push_back(onep);
					}
				}
			}
		}

		for (int k = 0; k < candidates_temp.size(); k++)
		{
			for (int m = 0; m < z_dims.size(); m++)
			{
				vector<int> onep = candidates_temp[k];
				onep.push_back(z_dimsint[m]);
				int hid = ParaFCs[FC_ID].hv_ids[onep[2]][onep[0]][onep[1]];
				if (array_hv_ids[hid])
				{
					candidates.push_back(onep);
					vector<float> onepc;
					onepc = ParaFCs[FC_ID].UVW_coords[0][onep[0]][onep[1]].coords;
					onepc[2] = z_dims[m];
					candidates_coords.push_back(onepc);
				}
			}
		}

		if (!candidates.size())
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

		for (int j = 0; j < candidates_coords.size(); j++)
		{
			int hvid = ParaFCs[FC_ID].hv_ids[candidates[j][2]][candidates[j][0]][candidates[j][1]];
			if (!array_hv_ids[hvid])
				continue;

			vector<float> lamdas;
			vector<float> deters;
			if (!tet_coplanar)
			{
				bool next = false;

				vector<float> cur_uvw;
				cur_uvw.push_back(candidates_coords[j][0]);
				cur_uvw.push_back(candidates_coords[j][1]);
				cur_uvw.push_back(candidates_coords[j][2]);
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
			else
			{
				for (int n = 0; n < 4; n++)
					lamdas.push_back(0.25);
			}

			hex_mesh.HVs[hvid].v[0] = hex_mesh.HVs[hvid].v[1] = hex_mesh.HVs[hvid].v[2] = 0;
			for (int n = 0; n < 4; n++)
			{
				hex_mesh.HVs[hvid].v[0] += lamdas[n] * vts[tets[i].vs[n]].v[0];
				hex_mesh.HVs[hvid].v[1] += lamdas[n] * vts[tets[i].vs[n]].v[1];
				hex_mesh.HVs[hvid].v[2] += lamdas[n] * vts[tets[i].vs[n]].v[2];
			}
			array_hv_ids[hvid] = false;
		}
	}

	for (int i = 0; i < array_hv_ids.size(); i++)
		if (array_hv_ids[i])
			printf("i %d\n", i);
//
	h_io io;
	vector<Hex> hs;
	int Color_id = 0;

	for (int i = 0; i < FrameI.FFs[FCs[FC_ID].f_id].neighbor_Cs.size(); i++)
	{
		Color_id++;
		int Fhid = FrameI.FFs[FCs[FC_ID].f_id].neighbor_Cs[i];
		for (int m = 0; m < FrameI.FHs[Fhid].hs_net.size(); m++)
		{
			int hid = FrameI.FHs[Fhid].hs_net[m];
			Hex h = hex_mesh.HHs[hid];
			h.Color_ID = Color_id;
			hs.push_back(h);
		}
	}
}
