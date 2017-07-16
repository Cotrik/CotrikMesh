#include "frame_of_mesh.h"

frame_of_mesh::frame_of_mesh(void)
{
}

//extract base-complex
void frame_of_mesh::base_complex_extraction()
{
	extract_singular_node_edge();
	singularity_structure_connectivity();

	extract_base_complex_node_edge();
	extract_base_complex_face();
	assign_hex_mesh_component();

	assign_singular_edge_composedBSEs();
}
void frame_of_mesh::extract_singular_node_edge()
{
	vector<int> arraye_test;
	for (int i = 0; i < hex_mesh.HEs.size(); i++)
		arraye_test.push_back(-1);

	int SV_Ind = 0, SE_Ind = 0;
	vector<Singular_E> circle_ses;
	for (int i = 0; i < hex_mesh.HEs.size(); i++)
	{
		if (arraye_test[i] != -1)
			continue;
		if ((hex_mesh.HEs[i].is_boundary == -1 && hex_mesh.HEs[i].neighbor_Hs.size() != 4)
				|| (hex_mesh.HEs[i].is_boundary == 1 && hex_mesh.HEs[i].neighbor_Hs.size() != 2))
		{
			bool is_circle = false;
			//marching left of the edge
			vector<int> es_left, vs_left;
			int current_e = i;
			es_left.push_back(current_e);
			arraye_test[current_e] = 1;

			int current_v = hex_mesh.HEs[current_e].startend_Id[0];
			vs_left.push_back(current_v);
			int next_e = -1;
			while (!is_v_singular(current_v, current_e, next_e))
			{
				es_left.push_back(next_e);
				current_e = next_e;
				arraye_test[current_e] = 1;
				if (hex_mesh.HEs[current_e].startend_Id[0] == current_v)
					current_v = hex_mesh.HEs[current_e].startend_Id[1];
				else
					current_v = hex_mesh.HEs[current_e].startend_Id[0];
				vs_left.push_back(current_v);
				if (current_e == i)
				{
					is_circle = true;
					break;
				}
			}

			if (is_circle)
			{
				Singular_E se;

				for (int j = 0; j < es_left.size() - 1; j++)
					se.es_link.push_back(es_left[j]);
				for (int j = 0; j < vs_left.size(); j++)
					se.vs_link.push_back(vs_left[j]);
				se.is_boundary = hex_mesh.HEs[es_left[0]].is_boundary;
				circle_ses.push_back(se);
				continue;
			}
			//test if the singular v added or not
			int which_left_sv = -1;
			bool addedornot = false;
			for (int j = 0; j < SingularityI.SVs.size(); j++)
			{
				if (SingularityI.SVs[j].index_hex == current_v)
				{
					addedornot = true;
					which_left_sv = j;
				}
			}
			if (!addedornot)
			{
				Singular_V sv;
				sv.index_own = SV_Ind++;
				sv.index_hex = current_v;
				sv.is_boundary = hex_mesh.HVs[current_v].where_location;
				SingularityI.SVs.push_back(sv);
				which_left_sv = sv.index_own;
			}

			//marching right of the edge
			vector<int> es_right, vs_right;
			current_e = i;
			current_v = hex_mesh.HEs[current_e].startend_Id[1];
			vs_right.push_back(current_v);
			next_e = -1;
			while (!is_v_singular(current_v, current_e, next_e))
			{
				es_right.push_back(next_e);
				current_e = next_e;
				arraye_test[current_e] = 1;
				if (hex_mesh.HEs[current_e].startend_Id[0] == current_v)
					current_v = hex_mesh.HEs[current_e].startend_Id[1];
				else
					current_v = hex_mesh.HEs[current_e].startend_Id[0];
				vs_right.push_back(current_v);
			}
			//test if the singular v added or not
			int which_right_sv = -1;
			addedornot = false;
			for (int j = 0; j < SingularityI.SVs.size(); j++)
			{
				if (SingularityI.SVs[j].index_hex == current_v)
				{
					addedornot = true;
					which_right_sv = j;
				}
			}
			if (!addedornot)
			{
				Singular_V sv;
				sv.index_own = SV_Ind++;
				sv.index_hex = current_v;
				sv.is_boundary = hex_mesh.HVs[current_v].where_location;
				SingularityI.SVs.push_back(sv);
				which_right_sv = sv.index_own;
			}

			Singular_E se;
			se.startend_Id[0] = which_left_sv;
			se.startend_Id[1] = which_right_sv;

			for (int j = es_left.size() - 1; j >= 0; j--)
				se.es_link.push_back(es_left[j]);
			for (int j = vs_left.size() - 1; j >= 0; j--)
				se.vs_link.push_back(vs_left[j]);
			for (int j = 0; j < es_right.size(); j++)
				se.es_link.push_back(es_right[j]);
			for (int j = 0; j < vs_right.size(); j++)
				se.vs_link.push_back(vs_right[j]);
			se.is_boundary = hex_mesh.HEs[es_left[0]].is_boundary;

			se.index_own = SE_Ind++;
			se.edge_type = 0;
			SingularityI.SEs.push_back(se);
		}
	}
	for (int i = 0; i < SingularityI.SVs.size(); i++)
		SingularityI.SVs[i].fake = false;

	vector<int> interior_Ids, interior_ses, boundary_Ids, boundary_ses;
	for (int i = 0; i < SingularityI.SEs.size(); i++)
	{
		vector<int> ses;
		if (SingularityI.SEs[i].is_boundary == -1)
		{
			int eid = SingularityI.SEs[i].es_link[0];
			int whichone = set_contain(interior_Ids, hex_mesh.HEs[eid].neighbor_Hs.size());
			if (whichone != -1)
			{
				interior_ses[whichone]++;
			}
			else
			{
				interior_ses.push_back(1);
				interior_Ids.push_back(hex_mesh.HEs[eid].neighbor_Hs.size());
			}
		}
		else
		{
			int eid = SingularityI.SEs[i].es_link[0];
			int whichone = set_contain(boundary_Ids, hex_mesh.HEs[eid].neighbor_Hs.size());
			if (whichone != -1)
			{
				boundary_ses[whichone]++;
			}
			else
			{
				boundary_ses.push_back(1);
				boundary_Ids.push_back(hex_mesh.HEs[eid].neighbor_Hs.size());
			}
		}
	}
	dealwith_circlesingularedge(circle_ses);
}
void frame_of_mesh::singularity_structure_connectivity()
{
	for (int i = 0; i < SingularityI.SEs.size(); i++)
	{
		int v1 = SingularityI.SEs[i].startend_Id[0], v2 = SingularityI.SEs[i].startend_Id[1];
		SingularityI.SVs[v1].neighbor_Vs.push_back(v2);
		SingularityI.SVs[v2].neighbor_Vs.push_back(v1);

		SingularityI.SVs[v1].neighbor_Es.push_back(i);
		SingularityI.SVs[v2].neighbor_Es.push_back(i);
	}
}
void frame_of_mesh::dealwith_circlesingularedge(vector<Singular_E> &circle_ses)
{
	if (circle_ses.size() > 0)
	{
		vector<int> arrayv_test;
		for (int i = 0; i < hex_mesh.HVs.size(); i++)
			arrayv_test.push_back(-2);
		for (int i = 0; i < SingularityI.SEs.size(); i++)
		{
			for (int j = 0; j < SingularityI.SEs[i].vs_link.size(); j++)
			{
				int vid = SingularityI.SEs[i].vs_link[j];
				arrayv_test[vid] = -1;
			}
		}
		for (int i = 0; i < SingularityI.SVs.size(); i++)
		{
			int vid = SingularityI.SVs[i].index_hex;
			arrayv_test[vid] = -3;
		}
		for (int i = 0; i < circle_ses.size(); i++)
		{
			for (int j = 0; j < circle_ses[i].vs_link.size(); j++)
			{
				int vid = circle_ses[i].vs_link[j];
				arrayv_test[vid] = i;
			}
		}
		//test intersection with current singular nodes
		vector<Singular_E> parallel_circle_ses;
		vector<vector<int> > Vfinds, Whichones;
		for (int i = 0; i < circle_ses.size(); i++)
		{
			vector<int> vfind;
			vector<int> which_ones;
			for (int j = 0; j < circle_ses[i].vs_link.size() - 1; j++)
			{
				int thisv = circle_ses[i].vs_link[j];

				Hex_V hv = hex_mesh.HVs[thisv];
				for (int k = 0; k < hex_mesh.HVs[thisv].neighbor_Es.size(); k++)
				{
					bool found = false;
					int currente = hex_mesh.HVs[thisv].neighbor_Es[k];
					if (currente == circle_ses[i].es_link[(1 + j) % circle_ses[i].es_link.size()] || currente == circle_ses[i].es_link[j])
						continue;

					int currentv = thisv;
					int nextv = hex_mesh.HEs[currente].startend_Id[0];
					if (nextv == thisv)
						nextv = hex_mesh.HEs[currente].startend_Id[1];
					int cout_step = 0;
					vector<int> asery;
					asery.push_back(currentv);

					while (true)
					{
						cout_step++;
						asery.push_back(nextv);
						if (arrayv_test[nextv] == -3)
						{
							vfind.push_back(thisv);
							which_ones.push_back(j);
							found = true;
							break;
						}
						bool foundnextv = false;
						for (int m = 0; m < hex_mesh.HVs[nextv].neighbor_vs.size(); m++)
						{
							int nv = hex_mesh.HVs[nextv].neighbor_vs[m];

							if (nv != currentv)
							{
								if (straight_line_test(nv, nextv, currentv) && arrayv_test[nv] != -1)
								{
									currentv = nextv;
									nextv = nv;
									foundnextv = true;
									break;
								}
							}
						}
						if (nextv == thisv || !foundnextv)
							break;
					}
					if (found)
						break;
				}
			}

			Vfinds.push_back(vfind);
			Whichones.push_back(which_ones);
		}

		//test intersection with no circle singular edge
		for (int i = 0; i < circle_ses.size(); i++)
		{
			if (Vfinds[i].size() != 0)
				continue;

			vector<int> vfind;
			vector<int> which_ones;
			for (int j = 0; j < circle_ses[i].vs_link.size() - 1; j++)
			{
				int thisv = circle_ses[i].vs_link[j];

				Hex_V hv = hex_mesh.HVs[thisv];
				for (int k = 0; k < hex_mesh.HVs[thisv].neighbor_Es.size(); k++)
				{
					bool found = false;
					int currente = hex_mesh.HVs[thisv].neighbor_Es[k];
					if (currente == circle_ses[i].es_link[(1 + j) % circle_ses[i].es_link.size()] || currente == circle_ses[i].es_link[j])
						continue;

					int currentv = thisv;
					int nextv = hex_mesh.HEs[currente].startend_Id[0];
					if (nextv == thisv)
						nextv = hex_mesh.HEs[currente].startend_Id[1];
					int cout_step = 0;
					vector<int> asery;
					asery.push_back(currentv);

					while (true)
					{
						cout_step++;
						asery.push_back(nextv);
						if (arrayv_test[nextv] == -1)
						{
							vfind.push_back(thisv);
							which_ones.push_back(j);
							found = true;
							break;
						}
						bool foundnextv = false;
						for (int m = 0; m < hex_mesh.HVs[nextv].neighbor_vs.size(); m++)
						{
							int nv = hex_mesh.HVs[nextv].neighbor_vs[m];
							if (nv != currentv)
							{
								if (straight_line_test(nv, nextv, currentv))
								{
									currentv = nextv;
									nextv = nv;
									foundnextv = true;
									break;
								}
							}
						}
						if (nextv == thisv || !foundnextv)
							break;
					}
					if (found)
						break;
				}
			}

			Vfinds[i] = vfind;
			Whichones[i] = which_ones;
		}
		//if all the rest circle singular edge are the same length
		int length_circle = -1;
		int first_circle = -1;
		for (int i = 0; i < circle_ses.size(); i++)
		{
			if (Vfinds[i].size() == 0)
			{
				length_circle = circle_ses[i].es_link.size();
				first_circle = i;
				break;
			}
		}
		if (length_circle != -1)
		{
			for (int i = 0; i < circle_ses.size(); i++)
			{
				if (Vfinds[i].size() == 0)
				{
					if (length_circle != circle_ses[i].es_link.size())
						; //exit(0);
				}
			}
		}
		//test intersection with circle singular edge
		if (first_circle != -1)
		{
			vector<int> circlevfind, which_ones;
			int onesection = circle_ses[first_circle].es_link.size() / 3; //3 is the minimum length of circle could be!
			for (int j = 0; j < 3; j++)
			{
				int thisv = circle_ses[first_circle].vs_link[j * onesection];
				circlevfind.push_back(thisv);
				which_ones.push_back(first_circle);

				for (int k = 0; k < hex_mesh.HVs[thisv].neighbor_Es.size(); k++)
				{
					int currente = hex_mesh.HVs[thisv].neighbor_Es[k];
					if (currente == circle_ses[first_circle].es_link[(1 + j * onesection) % circle_ses[first_circle].es_link.size()]
							|| currente == circle_ses[first_circle].es_link[j * onesection])
						continue;

					int currentv = thisv;
					int nextv = hex_mesh.HEs[currente].startend_Id[0];
					if (nextv == thisv)
						nextv = hex_mesh.HEs[currente].startend_Id[1];
					while (true)
					{
						if (arrayv_test[nextv] >= 0 && arrayv_test[nextv] != first_circle)
						{
							bool added = false;
							for (int m = 0; m < which_ones.size(); m++)
								if (which_ones[m] == arrayv_test[nextv])
									added = true;
							if (!added)
							{
								circlevfind.push_back(nextv);
								which_ones.push_back(arrayv_test[nextv]);
							}
						}
						bool foundnextv = false;
						for (int m = 0; m < hex_mesh.HVs[nextv].neighbor_vs.size(); m++)
						{
							int nv = hex_mesh.HVs[nextv].neighbor_vs[m];
							if (nv != currentv)
							{
								if (straight_line_test(nv, nextv, currentv))
								{
									currentv = nextv;
									nextv = nv;
									foundnextv = true;
									break;
								}
							}
						}
						if (nextv == thisv || !foundnextv)
							break;
					}
				}

				for (int k = 0; k < which_ones.size(); k++)
				{
					Vfinds[which_ones[k]].push_back(circlevfind[k]);
				}
			}
		}

// 		Singular_V sv;
// 		sv.index_hex=136;
// 		sv.is_boundary=-1;
// 		sv.index_own=SingularityI.SVs.size();
// 		SingularityI.SVs.push_back(sv);
// 		sv.index_hex=139;
// 		sv.is_boundary=-1;
// 		sv.index_own=SingularityI.SVs.size();
// 		SingularityI.SVs.push_back(sv);

// 		for(int i=0;i<circle_ses.size();i++)
// 			SingularityI.SEs.push_back(circle_ses[i]);

		for (int i = 0; i < Vfinds.size(); i++)
		{
			for (int j = 0; j < Vfinds[i].size(); j++)
			{
				//continue;;

				Singular_V sv;
				sv.index_hex = Vfinds[i][j];
				sv.is_boundary = hex_mesh.HVs[Vfinds[i][j]].where_location;
				sv.fake = true;
				sv.index_own = SingularityI.SVs.size();
				SingularityI.SVs.push_back(sv);

				//* for later use
				Singular_E se;
				se.edge_type = i + 1;
				se.startend_Id[0] = sv.index_own;
				se.startend_Id[1] = sv.index_own + 1;
				if (j == Vfinds[i].size() - 1)
					se.startend_Id[1] = sv.index_own - Vfinds[i].size() + 1;

				int start_p = -1, end_p = -1;
				for (int k = 0; k < circle_ses[i].vs_link.size() - 1; k++)
				{
					if (circle_ses[i].vs_link[k] == Vfinds[i][j])
						start_p = k;
					if (circle_ses[i].vs_link[k] == Vfinds[i][(j + 1) % Vfinds[i].size()])
						end_p = k;
				}
				int length_singulare = (end_p + circle_ses[i].vs_link.size() - 1 - start_p + 1) % (circle_ses[i].vs_link.size() - 1);

				for (int k = 0; k < length_singulare; k++)
					se.vs_link.push_back(circle_ses[i].vs_link[(k + start_p) % (circle_ses[i].vs_link.size() - 1)]);

				int start_e = -1, end_e = -1;
				for (int k = 0; k < circle_ses[i].es_link.size(); k++)
				{
					int e = circle_ses[i].es_link[k];
					if ((circle_ses[i].vs_link[start_p] == hex_mesh.HEs[e].startend_Id[0]
							&& circle_ses[i].vs_link[(start_p + 1) % (circle_ses[i].vs_link.size() - 1)] == hex_mesh.HEs[e].startend_Id[1])
							|| (circle_ses[i].vs_link[start_p] == hex_mesh.HEs[e].startend_Id[1]
									&& circle_ses[i].vs_link[(start_p + 1) % (circle_ses[i].vs_link.size() - 1)] == hex_mesh.HEs[e].startend_Id[0]))
						start_e = k;

					if ((circle_ses[i].vs_link[end_p] == hex_mesh.HEs[e].startend_Id[0]
							&& circle_ses[i].vs_link[(end_p - 1 + circle_ses[i].vs_link.size() - 1) % (circle_ses[i].vs_link.size() - 1)]
									== hex_mesh.HEs[e].startend_Id[1])
							|| (circle_ses[i].vs_link[end_p] == hex_mesh.HEs[e].startend_Id[1]
									&& circle_ses[i].vs_link[(end_p - 1 + circle_ses[i].vs_link.size() - 1) % (circle_ses[i].vs_link.size() - 1)]
											== hex_mesh.HEs[e].startend_Id[0]))
						end_e = k;
				}

				length_singulare = (end_e - start_e + circle_ses[i].es_link.size() + 1) % circle_ses[i].es_link.size();
				for (int k = 0; k < length_singulare; k++)
					se.es_link.push_back(circle_ses[i].es_link[(k + start_e) % (circle_ses[i].es_link.size())]);

				SingularityI.SEs.push_back(se);
				//*/
			}
		}
	}
}
void frame_of_mesh::extract_base_complex_node_edge()
{
	for (int i = 0; i < hex_mesh.HVs.size(); i++)
	{
		hex_mesh.HVs[i].Frame_V_id = -1;
		hex_mesh.HVs[i].is_Extra_ordinary = false;
		hex_mesh.HVs[i].on_whichFrame_edge = -1;
	}
	for (int i = 0; i < SingularityI.SEs.size(); i++)
	{
		for (int j = 0; j < SingularityI.SEs[i].es_link.size(); j++)
		{
			hex_mesh.HVs[hex_mesh.HEs[SingularityI.SEs[i].es_link[j]].startend_Id[0]].is_Extra_ordinary = true;
			hex_mesh.HVs[hex_mesh.HEs[SingularityI.SEs[i].es_link[j]].startend_Id[1]].is_Extra_ordinary = true;
		}
		for (int j = 0; j < SingularityI.SEs[i].vs_link.size(); j++)
			hex_mesh.HVs[SingularityI.SEs[i].vs_link[j]].neighbor_SEs.push_back(i);
	}

	vector<int> frame_nodes; //initialize frame nodes
	int FV_Ind = 0;
	for (int i = 0; i < SingularityI.SVs.size(); i++)
	{
		frame_nodes.push_back(SingularityI.SVs[i].index_hex);
		Frame_V fv;
		fv.index_hex = SingularityI.SVs[i].index_hex;
		fv.what_type = 1;
		fv.index_own = FV_Ind++;
		fv.which_singularity = i;
		hex_mesh.HVs[SingularityI.SVs[i].index_hex].Frame_V_id = fv.index_own;

		FrameI.FVs.push_back(fv);
	}
	// processing
	int FE_Ind = 0;
	while (frame_nodes.size() > 0)
	{
		int id = frame_nodes[frame_nodes.size() - 1];
		frame_nodes.pop_back();

		for (int i = 0; i < hex_mesh.HVs[id].neighbor_vs.size(); i++)
		{
			int next_id = hex_mesh.HVs[id].neighbor_vs[i];

			bool should_continue = true;
			for (int j = 0; j < FrameI.FVs[hex_mesh.HVs[id].Frame_V_id].neighbor_Es.size(); j++)
			{
				int eid = FrameI.FVs[hex_mesh.HVs[id].Frame_V_id].neighbor_Es[j];
				if (FrameI.FEs[eid].startend_Id[0] == hex_mesh.HVs[id].Frame_V_id)
					if (next_id == FrameI.FEs[eid].vs_link[1])
						should_continue = false;
				if (FrameI.FEs[eid].startend_Id[1] == hex_mesh.HVs[id].Frame_V_id)
					if (next_id == FrameI.FEs[eid].vs_link[FrameI.FEs[eid].vs_link.size() - 2])
						should_continue = false;

				if (!should_continue)
					break;
			}
			if (!should_continue)
				continue;

			Frame_E fe;
			fe.index_own = FE_Ind++;
			fe.startend_Id[0] = hex_mesh.HVs[id].Frame_V_id;
			fe.vs_link.push_back(id);
			fe.vs_link.push_back(next_id);

			FrameI.FVs[hex_mesh.HVs[id].Frame_V_id].neighbor_Es.push_back(fe.index_own);

			int pre_id = id;
			bool find_next = true;
			while (true)
			{
				vector<int> sharedse;
				set_cross(hex_mesh.HVs[pre_id].neighbor_SEs, hex_mesh.HVs[next_id].neighbor_SEs, sharedse);
				set_redundent_clearn(sharedse);

				//judge if new frame_node
				if (hex_mesh.HVs[next_id].Frame_V_id >= 0 && find_next)
				{ //frame node
					fe.startend_Id[1] = hex_mesh.HVs[next_id].Frame_V_id;
					FrameI.FEs.push_back(fe);

					FrameI.FVs[hex_mesh.HVs[next_id].Frame_V_id].neighbor_Es.push_back(fe.index_own);
					break;
				}
				else if (hex_mesh.HVs[next_id].on_whichFrame_edge >= 0 && find_next)
				{ //not frame node or extra-ordinary node, but on frame edge
					int on_edge = hex_mesh.HVs[next_id].on_whichFrame_edge;

					frame_nodes.push_back(next_id);
					Frame_V fv;
					fv.index_hex = next_id;

					fv.what_type = 1;

					fv.index_own = FV_Ind++;
					hex_mesh.HVs[next_id].Frame_V_id = fv.index_own;
					FrameI.FVs.push_back(fv);

					fe.startend_Id[1] = hex_mesh.HVs[next_id].Frame_V_id;

					FrameI.FVs[hex_mesh.HVs[next_id].Frame_V_id].neighbor_Es.push_back(fe.index_own);
					FrameI.FEs.push_back(fe);
					//splitting the meeting edge into two edges
					Frame_E fe_2;
					int edge1_start = FrameI.FEs[on_edge].startend_Id[0];
					int edge2_end = FrameI.FEs[on_edge].startend_Id[1];
					FrameI.FEs[on_edge].startend_Id[1] = fv.index_own;

					FrameI.FVs[fv.index_own].neighbor_Es.push_back(on_edge);

					fe_2.startend_Id[0] = fv.index_own;
					fe_2.startend_Id[1] = edge2_end;
					fe_2.index_own = FE_Ind++;

					vector<int> vlink = FrameI.FEs[on_edge].vs_link;
					FrameI.FEs[on_edge].vs_link.clear();

					bool which_edge = true;
					for (int j = 0; j < vlink.size(); j++)
					{
						if (vlink[j] == next_id)
							which_edge = false;
						if (which_edge)
							FrameI.FEs[on_edge].vs_link.push_back(vlink[j]);
						else
						{
							if (hex_mesh.HVs[vlink[j]].Frame_V_id == -1)
								hex_mesh.HVs[vlink[j]].on_whichFrame_edge = fe_2.index_own;
							fe_2.vs_link.push_back(vlink[j]);
						}
					}
					FrameI.FEs[on_edge].vs_link.push_back(next_id);
					FrameI.FEs.push_back(fe_2);

					vector<int> edge2_end_nedges = FrameI.FVs[edge2_end].neighbor_Es;
					FrameI.FVs[edge2_end].neighbor_Es.clear();

					for (int j = 0; j < edge2_end_nedges.size(); j++)
					{
						if (edge2_end != edge1_start)
						{
							if (edge2_end_nedges[j] == on_edge)
								continue;
						}
						FrameI.FVs[edge2_end].neighbor_Es.push_back(edge2_end_nedges[j]);
					}
					FrameI.FVs[edge2_end].neighbor_Es.push_back(fe_2.index_own);
					FrameI.FVs[fv.index_own].neighbor_Es.push_back(fe_2.index_own);
					break;
				}
				else if (hex_mesh.HVs[next_id].is_Extra_ordinary && sharedse.size() == 0 && find_next)
				{ //not frame node,but extra-ordinary node
					frame_nodes.push_back(next_id);
					Frame_V fv;
					fv.index_hex = next_id;
					fv.what_type = 1;
					fv.index_own = FV_Ind++;
					hex_mesh.HVs[next_id].Frame_V_id = fv.index_own;
					FrameI.FVs.push_back(fv);

					fe.startend_Id[1] = hex_mesh.HVs[next_id].Frame_V_id;
					FrameI.FEs.push_back(fe);
					FrameI.FVs[fv.index_own].neighbor_Es.push_back(fe.index_own);
					break;

				}
				else if (!find_next)
				{ //reached the boundary
					frame_nodes.push_back(next_id);
					Frame_V fv;
					fv.index_hex = next_id;
					fv.what_type = 1;
					fv.index_own = FV_Ind++;
					hex_mesh.HVs[next_id].Frame_V_id = fv.index_own;
					FrameI.FVs.push_back(fv);

					fe.startend_Id[1] = fv.index_own;
					FrameI.FEs.push_back(fe);

					FrameI.FVs[fv.index_own].neighbor_Es.push_back(fe.index_own);
					break;
				}

				find_next = false;
				//if not break above, continue increasing edge
				hex_mesh.HVs[next_id].on_whichFrame_edge = fe.index_own;

				for (int j = 0; j < hex_mesh.HVs[next_id].neighbor_vs.size(); j++)
				{
					int nid = hex_mesh.HVs[next_id].neighbor_vs[j];
					if (nid != pre_id)
					{
						int shared_edge = -1;
						bool found = false;
						for (int m = 0; m < hex_mesh.HVs[pre_id].neighbor_Es.size(); m++)
						{
							int pre_ne = hex_mesh.HVs[pre_id].neighbor_Es[m];
							for (int n = 0; n < hex_mesh.HVs[next_id].neighbor_Es.size(); n++)
							{
								int next_ne = hex_mesh.HVs[next_id].neighbor_Es[n];
								if (pre_ne == next_ne)
								{
									shared_edge = pre_ne;
									found = true;
								}
								if (found)
									break;
							}
							if (found)
								break;
						}

						bool is_invalid_nid = false;
						for (int m = 0; m < hex_mesh.HEs[shared_edge].neighbor_Hs.size(); m++)
						{
							int shared_nh = hex_mesh.HEs[shared_edge].neighbor_Hs[m];
							for (int n = 0; n < 8; n++)
							{
								if (hex_mesh.HHs[shared_nh].V_Ids[n] == nid)
									is_invalid_nid = true;
								if (is_invalid_nid)
									break;
							}
							if (is_invalid_nid)
								break;
						}

						if (!is_invalid_nid)
						{
							pre_id = next_id;
							next_id = nid;
							fe.vs_link.push_back(next_id);
							find_next = true;
							break;
						}
					}
				}
			}
		}
	}
	for (int i = 0; i < FrameI.FVs.size(); i++)
	{
		FrameI.FVs[i].neighbor_vs.clear();
		FrameI.FVs[i].neighbor_Es.clear();
	}

	{
		int E_N = 0;
		vector<Frame_E> tempFes = FrameI.FEs;
		FrameI.FEs.clear();
		for (int i = 0; i < tempFes.size(); i++)
		{
			bool havesame = false;

			int id1 = tempFes[i].startend_Id[0];
			int id2 = tempFes[i].startend_Id[1];

			for (int j = 0; j < FrameI.FVs[id1].neighbor_vs.size(); j++)
			{
				if (FrameI.FVs[id1].neighbor_vs[j] == id2)
				{
					//test if they are the same edge or not
					int nfe = FrameI.FVs[id1].neighbor_Es[j];
					if (tempFes[i].vs_link.size() == tempFes[nfe].vs_link.size() && tempFes[nfe].vs_link.size() != 2)
					{
						int vstart_2 = tempFes[i].vs_link[1];
						int whichone = set_contain(tempFes[nfe].vs_link, vstart_2);
						if (whichone != -1)
						{
							havesame = true;
							break;
						}
					}
					else if (tempFes[i].vs_link.size() == tempFes[nfe].vs_link.size() && tempFes[nfe].vs_link.size() == 2)
					{
						havesame = true;
						break;
					}
				}
			}
			if (!havesame)
			{
				tempFes[i].index_own = E_N++;

				FrameI.FEs.push_back(tempFes[i]);
				FrameI.FVs[id1].neighbor_Es.push_back(tempFes[i].index_own);
				FrameI.FVs[id1].neighbor_vs.push_back(id2);

				FrameI.FVs[id2].neighbor_Es.push_back(tempFes[i].index_own);
				FrameI.FVs[id2].neighbor_vs.push_back(id1);
			}
		}
	}

	for (int i = 0; i < FrameI.FVs.size(); i++)
	{
		for (int j = 0; j < FrameI.FVs[i].neighbor_Es.size(); j++)
			if (FrameI.FEs[FrameI.FVs[i].neighbor_Es[j]].startend_Id[0] == FrameI.FVs[i].index_own)
				FrameI.FVs[i].neighbor_vs.push_back(FrameI.FEs[FrameI.FVs[i].neighbor_Es[j]].startend_Id[1]);
			else if (FrameI.FEs[FrameI.FVs[i].neighbor_Es[j]].startend_Id[1] == FrameI.FVs[i].index_own)
				FrameI.FVs[i].neighbor_vs.push_back(FrameI.FEs[FrameI.FVs[i].neighbor_Es[j]].startend_Id[0]);
	}

	for (int i = 0; i < FrameI.FVs.size(); i++)
	{
		vector<int> temp_neighbor;
		for (int j = 0; j < FrameI.FVs[i].neighbor_Es.size(); j++)
		{
			bool already = false;
			for (int k = 0; k < temp_neighbor.size(); k++)
				if (FrameI.FVs[i].neighbor_Es[j] == temp_neighbor[k])
					already = true;
			if (!already)
				temp_neighbor.push_back(FrameI.FVs[i].neighbor_Es[j]);
		}
		FrameI.FVs[i].neighbor_Es = temp_neighbor;

		temp_neighbor.clear();

		for (int j = 0; j < FrameI.FVs[i].neighbor_vs.size(); j++)
		{
			bool already = false;
			for (int k = 0; k < temp_neighbor.size(); k++)
				if (FrameI.FVs[i].neighbor_vs[j] == temp_neighbor[k])
					already = true;
			if (!already)
				temp_neighbor.push_back(FrameI.FVs[i].neighbor_vs[j]);
		}
		FrameI.FVs[i].neighbor_vs = temp_neighbor;
	}

}
void frame_of_mesh::extract_base_complex_face()
{
	for (int i = 0; i < FrameI.FEs.size(); i++)
	{ //find opposite edges for edge e
		vector<int> V1_neighbors, V2_neighbors;

		for (int j = 0; j < FrameI.FVs[FrameI.FEs[i].startend_Id[0]].neighbor_vs.size(); j++)
			if (FrameI.FVs[FrameI.FEs[i].startend_Id[0]].neighbor_vs[j] != FrameI.FEs[i].startend_Id[1])
				V1_neighbors.push_back(FrameI.FVs[FrameI.FEs[i].startend_Id[0]].neighbor_vs[j]);

		for (int j = 0; j < FrameI.FVs[FrameI.FEs[i].startend_Id[1]].neighbor_vs.size(); j++)
			if (FrameI.FVs[FrameI.FEs[i].startend_Id[1]].neighbor_vs[j] != FrameI.FEs[i].startend_Id[0])
				V2_neighbors.push_back(FrameI.FVs[FrameI.FEs[i].startend_Id[1]].neighbor_vs[j]);

		for (int j = 0; j < V1_neighbors.size(); j++)
		{
			for (int k = 0; k < FrameI.FVs[V1_neighbors[j]].neighbor_Es.size(); k++)
			{
				int which_e = FrameI.FVs[V1_neighbors[j]].neighbor_Es[k];
				for (int m = 0; m < V2_neighbors.size(); m++)
				{
					for (int n = 0; n < FrameI.FVs[V2_neighbors[m]].neighbor_Es.size(); n++)
					{
						if (which_e == FrameI.FVs[V2_neighbors[m]].neighbor_Es[n])
							FrameI.FEs[i].opposite_Es.push_back(which_e);
					}
				}
			}
		}
	}
//test!!!
// 	Frame_V fv1,fv2;
// 	for(int i=0;i<FrameI.FVs.size();i++)
// 	{
// 		if(FrameI.FVs[i].index_hex==3963)
// 			fv1=FrameI.FVs[i];
// 		if(FrameI.FVs[i].index_hex==3956)
// 			fv2=FrameI.FVs[i];
// 	}
// 	vector<int> sharedfe;
// 	vector<Hex_E> bes;
// 	set_cross(fv1.neighbor_Es,fv2.neighbor_Es,sharedfe);
// 	Frame_E fe=FrameI.FEs[sharedfe[0]];
// 	for(int i=0;i<fe.opposite_Es.size();i++)
// 	{
// 		int ofeid=fe.opposite_Es[i];
// 		for(int m=1;m<FrameI.FEs[ofeid].vs_link.size();m++)
// 		{
// 			int v1=FrameI.FEs[ofeid].vs_link[m-1];
// 			int v2=FrameI.FEs[ofeid].vs_link[m];
// 			vector<int> sharede;
// 			set_cross(hex_mesh.HVs[v1].neighbor_Es,hex_mesh.HVs[v2].neighbor_Es,sharede);
// 			bes.push_back(hex_mesh.HEs[sharede[0]]);
// 		}
// 	}
// 	char fname[300]="D:/xgao/optimization/simplification/data/t-h/cube_low/cube_low_S162_fe_opposite_es.vtk";
// 	vector<float> Vs_properties;
// 	h_io io;
// 	io.write_VTK_V(hex_mesh.HVs,bes,Vs_properties,fname);

	int Face_Count = 0;
	for (int i = 0; i < FrameI.FEs.size(); i++)
	{ //find all faces for edge e
		for (int j = 0; j < FrameI.FEs[i].opposite_Es.size(); j++)
		{ //for each opposit edge
			int corner1 = FrameI.FEs[i].startend_Id[0];
			int corner2 = FrameI.FEs[i].startend_Id[1];
			int corner3, corner4;
			for (int k = 0; k < FrameI.FVs[corner2].neighbor_vs.size(); k++)
			{
				if (FrameI.FEs[i].index_own != FrameI.FEs[i].opposite_Es[j])
				{
					if (FrameI.FVs[corner2].neighbor_vs[k] == FrameI.FEs[FrameI.FEs[i].opposite_Es[j]].startend_Id[0])
					{
						corner3 = FrameI.FEs[FrameI.FEs[i].opposite_Es[j]].startend_Id[0];
						corner4 = FrameI.FEs[FrameI.FEs[i].opposite_Es[j]].startend_Id[1];
						break;
					}
					else if (FrameI.FVs[corner2].neighbor_vs[k] == FrameI.FEs[FrameI.FEs[i].opposite_Es[j]].startend_Id[1])
					{
						corner3 = FrameI.FEs[FrameI.FEs[i].opposite_Es[j]].startend_Id[1];
						corner4 = FrameI.FEs[FrameI.FEs[i].opposite_Es[j]].startend_Id[0];
						break;
					}
				}
				else
				{ //added for dealing with torus
					if (corner2 == FrameI.FEs[FrameI.FEs[i].opposite_Es[j]].startend_Id[0])
					{
						corner3 = FrameI.FEs[FrameI.FEs[i].opposite_Es[j]].startend_Id[0];
						corner4 = FrameI.FEs[FrameI.FEs[i].opposite_Es[j]].startend_Id[1];
						break;
					}
					else if (corner2 == FrameI.FEs[FrameI.FEs[i].opposite_Es[j]].startend_Id[1])
					{
						corner3 = FrameI.FEs[FrameI.FEs[i].opposite_Es[j]].startend_Id[1];
						corner4 = FrameI.FEs[FrameI.FEs[i].opposite_Es[j]].startend_Id[0];
						break;
					}
				}
			}

			int edge1 = FrameI.FEs[i].index_own;
			int edge2, edge3 = FrameI.FEs[i].opposite_Es[j], edge4;
			for (int k = 0; k < FrameI.FVs[corner2].neighbor_Es.size(); k++)
			{ //edge2
				int temp_e = FrameI.FVs[corner2].neighbor_Es[k];
				if (FrameI.FEs[temp_e].startend_Id[0] == corner2 && FrameI.FEs[temp_e].startend_Id[1] == corner3)
				{
					edge2 = temp_e;
					break;
				}
				else if (FrameI.FEs[temp_e].startend_Id[1] == corner2 && FrameI.FEs[temp_e].startend_Id[0] == corner3)
				{
					edge2 = temp_e;
					break;
				}
			}
			for (int k = 0; k < FrameI.FVs[corner4].neighbor_Es.size(); k++)
			{ //edge4
				int temp_e = FrameI.FVs[corner4].neighbor_Es[k];
				if (FrameI.FEs[temp_e].startend_Id[0] == corner4 && FrameI.FEs[temp_e].startend_Id[1] == corner1)
				{
					edge4 = temp_e;
					break;
				}
				else if (FrameI.FEs[temp_e].startend_Id[1] == corner4 && FrameI.FEs[temp_e].startend_Id[0] == corner1)
				{
					edge4 = temp_e;
					break;
				}
			}

			//judge formed face before or not
			bool formedbefore = false;
			for (int k = 0; k < FrameI.FEs[edge1].neighbor_Fs.size(); k++)
			{
				int fid = FrameI.FEs[edge1].neighbor_Fs[k];
				bool thisone = true;
				for (int m = 0; m < 4; m++)
					if (FrameI.FFs[fid].fe_Ids[m] != edge1 && FrameI.FFs[fid].fe_Ids[m] != edge2 && FrameI.FFs[fid].fe_Ids[m] != edge3
							&& FrameI.FFs[fid].fe_Ids[m] != edge4)
						thisone = false;
				if (thisone)
					formedbefore = true;
			}
			if (formedbefore)
				continue;

			if (FrameI.FEs[edge1].vs_link.size() != FrameI.FEs[edge3].vs_link.size()
					|| FrameI.FEs[edge2].vs_link.size() != FrameI.FEs[edge4].vs_link.size())
				continue;

			//start form
			Frame_F ff;
			ff.fv_Ids[0] = corner1;
			ff.fv_Ids[1] = corner2;
			ff.fv_Ids[2] = corner3;
			ff.fv_Ids[3] = corner4;
			ff.fe_Ids[0] = edge1;
			ff.fe_Ids[1] = edge2;
			ff.fe_Ids[2] = edge3;
			ff.fe_Ids[3] = edge4;
			if (FrameI.FEs[edge1].vs_link[0] != FrameI.FVs[corner1].index_hex)
			{
				vector<int> reverse_vslist = FrameI.FEs[edge1].vs_link;
				FrameI.FEs[edge1].vs_link.clear();
				for (int k = reverse_vslist.size() - 1; k >= 0; k--)
					FrameI.FEs[edge1].vs_link.push_back(reverse_vslist[k]);
			}
			if (FrameI.FEs[edge4].vs_link[0] != FrameI.FVs[corner1].index_hex)
			{
				vector<int> reverse_vslist = FrameI.FEs[edge4].vs_link;
				FrameI.FEs[edge4].vs_link.clear();
				for (int k = reverse_vslist.size() - 1; k >= 0; k--)
					FrameI.FEs[edge4].vs_link.push_back(reverse_vslist[k]);
			}
			if (FrameI.FEs[edge3].vs_link[0] != FrameI.FVs[corner4].index_hex)
			{
				vector<int> reverse_vslist = FrameI.FEs[edge3].vs_link;
				FrameI.FEs[edge3].vs_link.clear();
				for (int k = reverse_vslist.size() - 1; k >= 0; k--)
					FrameI.FEs[edge3].vs_link.push_back(reverse_vslist[k]);
			}

			bool valid_face = true;
			//scan
			ff.U = FrameI.FEs[edge4].vs_link.size() - 1;
			ff.V = FrameI.FEs[edge1].vs_link.size() - 1;
			vector<int> edge1_vslink = FrameI.FEs[edge1].vs_link;

			for (int m = 0; m < FrameI.FEs[edge4].vs_link.size() - 1; m++)
			{
				int v_m1 = FrameI.FEs[edge4].vs_link[m];
				int v_m2 = FrameI.FEs[edge4].vs_link[m + 1];
				vector<int> next_edge1_vslink;
				next_edge1_vslink.push_back(v_m2);
				vector<int> oneline_face;
				for (int n = 1; n < edge1_vslink.size(); n++)
				{
					int v_n1 = edge1_vslink[n];

					int shared_quad = -1;
					bool found = false;
					for (int o = 0; o < hex_mesh.HVs[v_m1].neighbor_Fs.size(); o++)
					{
						for (int p = 0; p < hex_mesh.HVs[v_m2].neighbor_Fs.size(); p++)
						{
							for (int q = 0; q < hex_mesh.HVs[v_n1].neighbor_Fs.size(); q++)
							{
								if (hex_mesh.HVs[v_m1].neighbor_Fs[o] == hex_mesh.HVs[v_m2].neighbor_Fs[p]
										&& hex_mesh.HVs[v_m2].neighbor_Fs[p] == hex_mesh.HVs[v_n1].neighbor_Fs[q])
								{
									shared_quad = hex_mesh.HVs[v_m1].neighbor_Fs[o];
									found = true;
									if (m == FrameI.FEs[edge4].vs_link.size() - 2) //temporarily solve two quads share three points.
									{
										bool contain_fourthv = false;
										for (int r = 0; r < 4; r++)
										{
											if (hex_mesh.HFs[shared_quad].cv_Ids[r] == FrameI.FEs[edge3].vs_link[n])
											{
												contain_fourthv = true;
											}
										}
										if (!contain_fourthv)
											found = false;
									}
								}
								if (found)
									break;
							}
							if (found)
								break;
						}
						if (found)
							break;
					}

					if (!found || shared_quad == -1)
					{
						valid_face = false;
						break;
					}

					oneline_face.push_back(shared_quad);
					int v_n2;
					for (int o = 0; o < 4; o++)
						if (hex_mesh.HFs[shared_quad].cv_Ids[o] == v_m1 || hex_mesh.HFs[shared_quad].cv_Ids[o] == v_m2
								|| hex_mesh.HFs[shared_quad].cv_Ids[o] == v_n1)
							continue;
						else
							v_n2 = hex_mesh.HFs[shared_quad].cv_Ids[o];

// 					if(hex_mesh.HVs[v_n2].Frame_V_id>=0&&!(m==FrameI.FEs[edge4].vs_link.size()-1&&n==edge1_vslink.size()-1))
// 					{
// 						valid_face=false;
// 						break;
// 					}

					next_edge1_vslink.push_back(v_n2);

					v_m1 = v_n1;
					v_m2 = v_n2;
				}

				if (!valid_face)
					break;

				ff.hfs_net.push_back(oneline_face);
				edge1_vslink.clear();
				edge1_vslink = next_edge1_vslink;
			}

			if (!valid_face)
				continue;

			ff.index_own = Face_Count++;
			ff.Color_ID = 0;

			for (int m = 0; m < ff.hfs_net.size(); m++)
			{
				for (int n = 0; n < ff.hfs_net[m].size(); n++)
					ff.hfs_net_another.push_back(ff.hfs_net[m][n]);
			}

			FrameI.FFs.push_back(ff);

			FrameI.FVs[corner1].neighbor_Fs.push_back(ff.index_own);
			FrameI.FVs[corner2].neighbor_Fs.push_back(ff.index_own);
			FrameI.FVs[corner3].neighbor_Fs.push_back(ff.index_own);
			FrameI.FVs[corner4].neighbor_Fs.push_back(ff.index_own);

			FrameI.FEs[edge1].neighbor_Fs.push_back(ff.index_own);
			FrameI.FEs[edge2].neighbor_Fs.push_back(ff.index_own);
			FrameI.FEs[edge3].neighbor_Fs.push_back(ff.index_own);
			FrameI.FEs[edge4].neighbor_Fs.push_back(ff.index_own);
		}
	}
	for (int i = 0; i < FrameI.FEs.size(); i++)
	{
		int sv = FrameI.FEs[i].startend_Id[0];
		if (FrameI.FVs[sv].index_hex != FrameI.FEs[i].vs_link[0])
		{
			vector<int> vslink = FrameI.FEs[i].vs_link;
			FrameI.FEs[i].vs_link.clear();
			for (int j = vslink.size() - 1; j >= 0; j--)
			{
				FrameI.FEs[i].vs_link.push_back(vslink[j]);
			}
		}
	}
}

bool frame_of_mesh::is_v_singular(int hvid, int cur_e, int &next_e)
{
	int num1 = 0, num2 = 0;
	for (int i = 0; i < hex_mesh.HVs[hvid].neighbor_Es.size(); i++)
	{
		int ne = hex_mesh.HVs[hvid].neighbor_Es[i];
		if (ne == cur_e)
			continue;
		if (hex_mesh.HEs[ne].neighbor_Hs.size() == hex_mesh.HEs[cur_e].neighbor_Hs.size())
		{
			num1++;
			next_e = ne;
		}
		else if ((hex_mesh.HEs[ne].is_boundary == -1 && hex_mesh.HEs[ne].neighbor_Hs.size() != 4)
				|| (hex_mesh.HEs[ne].is_boundary == 1 && hex_mesh.HEs[ne].neighbor_Hs.size() != 2))
			num2++;

	}
	if (num1 == 1 && num2 == 0)
		return false;

	return true;
}
bool frame_of_mesh::straight_line_test(int h_v1, int h_v2, int h_v3)
{
	vector<int> nei = hex_mesh.HVs[h_v1].neighbor_Fs;
	vector<int> nei2 = hex_mesh.HVs[h_v2].neighbor_Fs;
	vector<int> nei3 = hex_mesh.HVs[h_v3].neighbor_Fs;

	for (int i = 0; i < hex_mesh.HVs[h_v1].neighbor_Fs.size(); i++)
	{
		for (int j = 0; j < hex_mesh.HVs[h_v3].neighbor_Fs.size(); j++)
		{
			if (hex_mesh.HVs[h_v1].neighbor_Fs[i] == hex_mesh.HVs[h_v3].neighbor_Fs[j])
				return false;
		}
	}
	return true;
}

void frame_of_mesh::assign_singular_edge_composedBSEs()
{
	for (int i = 0; i < hex_mesh.HVs.size(); i++)
	{
		hex_mesh.HVs[i].on_whichFrame_edge = -1;
	}
	for (int i = 0; i < FrameI.FEs.size(); i++)
	{
		FrameI.FEs[i].singular_Id = -1;
	}
	for (int i = 0; i < SingularityI.SEs.size(); i++)
	{
		SingularityI.SEs[i].composed_BS_Es.clear();

		for (int j = 0; j < SingularityI.SEs[i].vs_link.size(); j++)
			hex_mesh.HVs[SingularityI.SEs[i].vs_link[j]].on_whichFrame_edge = 1;

		for (int j = 0; j < FrameI.FEs.size(); j++)
		{
			bool all = true;
			for (int k = 0; k < FrameI.FEs[j].vs_link.size(); k++)
			{
				int id = FrameI.FEs[j].vs_link[k];
				if (hex_mesh.HVs[id].on_whichFrame_edge != 1)
				{
					all = false;
					break;
				}
			}
			if (all)
			{
				SingularityI.SEs[i].composed_BS_Es.push_back(j);
				FrameI.FEs[j].singular_Id = i;
			}
		}
		for (int j = 0; j < SingularityI.SEs[i].vs_link.size(); j++)
			hex_mesh.HVs[SingularityI.SEs[i].vs_link[j]].on_whichFrame_edge = -1;

		if (SingularityI.SEs[i].edge_type == 0) //non circle
		{
			int startv = SingularityI.SVs[SingularityI.SEs[i].startend_Id[0]].index_hex;
			vector<int> orderedBSEs, recorder;
			for (int j = 0; j < SingularityI.SEs[i].composed_BS_Es.size(); j++)
				recorder.push_back(0);
			for (int j = 0; j < SingularityI.SEs[i].composed_BS_Es.size(); j++)
			{
				for (int k = 0; k < SingularityI.SEs[i].composed_BS_Es.size(); k++)
				{
					if (recorder[k] == 0)
					{
						int eid = SingularityI.SEs[i].composed_BS_Es[k];
						int v1, v2;
						v1 = FrameI.FVs[FrameI.FEs[eid].startend_Id[0]].index_hex;
						v2 = FrameI.FVs[FrameI.FEs[eid].startend_Id[1]].index_hex;
						if (v1 == startv)
						{
							recorder[k] = 1;
							orderedBSEs.push_back(eid);
							startv = v2;
							break;
						}
						else if (v2 == startv)
						{
							recorder[k] = 1;
							orderedBSEs.push_back(eid);
							startv = v1;
							break;
						}

					}
				}
			}
			SingularityI.SEs[i].composed_BS_Es = orderedBSEs;
		}
		else
		{

		}
	}
	for (int i = 0; i < hex_mesh.HEs.size(); i++)
	{
		hex_mesh.HEs[i].which_F_edge = -1;
	}
	for (int i = 0; i < FrameI.FEs.size(); i++)
	{
		for (int j = 0; j < FrameI.FEs[i].vs_link.size(); j++)
		{
			vector<int> neighbore;
			int v1 = FrameI.FEs[i].vs_link[j], v2 = FrameI.FEs[i].vs_link[(j + 1) % FrameI.FEs[i].vs_link.size()];
			set_cross(hex_mesh.HVs[v1].neighbor_Es, hex_mesh.HVs[v2].neighbor_Es, neighbore);
			if (neighbore.size())
				hex_mesh.HEs[neighbore[0]].which_F_edge = i;
		}
	}
}
//assign hex mesh component
void frame_of_mesh::assign_hex_mesh_component()
{
	for (int i = 0; i < hex_mesh.HFs.size(); i++)
	{
		hex_mesh.HFs[i].frame_boundary = -1;
	}
	vector<vector<int> > contain_samefs;
	for (int i = 0; i < FrameI.FFs.size(); i++)
	{
		for (int j = 0; j < FrameI.FFs[i].hfs_net_another.size(); j++)
		{
			if (hex_mesh.HFs[FrameI.FFs[i].hfs_net_another[j]].frame_boundary != -1)
			{
				int cur_fb = hex_mesh.HFs[FrameI.FFs[i].hfs_net_another[j]].frame_boundary;
				bool found = false;
				for (int m = 0; m < contain_samefs.size(); m++)
				{
					if (contain_samefs[m][0] == cur_fb)
					{
						contain_samefs[m].push_back(i);
						found = true;
					}
				}
				if (!found)
				{
					vector<int> newsery;
					newsery.push_back(cur_fb);
					newsery.push_back(i);
					contain_samefs.push_back(newsery);
				}
			}
			else
				hex_mesh.HFs[FrameI.FFs[i].hfs_net_another[j]].frame_boundary = FrameI.FFs[i].index_own;
		}
	}

	for (int i = 0; i < hex_mesh.HHs.size(); i++)
		hex_mesh.HHs[i].frame_component = -1;

	assign_component();
	assign_color_hexfaces_component();
}
void frame_of_mesh::assign_component()
{
	vector<bool> hexes_Inds;
	initializeVectorT(hexes_Inds, true, hex_mesh.HHs.size());
	extract_useless_hexes(hexes_Inds);

	int component_Id = -1;
	bool loop = true;

	vector<int> seeds;
	int start_hex = -1;
	for (int i = 0; i < hex_mesh.HHs.size(); i++)
	{
		if (hexes_Inds[i])
		{
			start_hex = i;
			seeds.push_back(start_hex);
			break;
		}
	}
	while (loop)
	{
		component_Id++;

		while (seeds.size() > 0)
		{
			int h_id = seeds[seeds.size() - 1];
			seeds.pop_back();
			hex_mesh.HHs[h_id].frame_component = component_Id;

			for (int i = 0; i < 6; i++)
			{
				int fid = hex_mesh.HHs[h_id].F_Ids[i];

				if (hex_mesh.HFs[fid].frame_boundary > -1)
					continue;

				for (int j = 0; j < hex_mesh.HFs[fid].neighbor_Cs.size(); j++)
				{
					int neighbor_hid = hex_mesh.HFs[fid].neighbor_Cs[j];
					if (neighbor_hid != h_id && hex_mesh.HHs[neighbor_hid].frame_component == -1 && hexes_Inds[neighbor_hid])
						seeds.push_back(neighbor_hid);
				}
			}
		}

		start_hex = -1;
		for (int i = 0; i < hex_mesh.HHs.size(); i++)
		{
			if (hex_mesh.HHs[i].frame_component == -1 && hexes_Inds[i])
			{
				start_hex = i;
				seeds.push_back(start_hex);
				break;
			}
		}
		if (start_hex == -1)
			loop = false;
	}

	for (int i = 0; i < component_Id + 1; i++)
	{
		Frame_H fh;
		fh.Color_ID = -1;
		fh.index_own = -1;
		FrameI.FHs.push_back(fh);
	}

	//collect frameH's hexes
	for (int i = 0; i < hex_mesh.HHs.size(); i++)
	{
		int comp_id = hex_mesh.HHs[i].frame_component;
		if (comp_id == -1)
			continue;
		FrameI.FHs[comp_id].hs_net.push_back(hex_mesh.HHs[i].index);

	}
	//collect frameH's FrameFs/6
	for (int i = 0; i < FrameI.FHs.size(); i++)
	{
		for (int j = 0; j < FrameI.FHs[i].hs_net.size(); j++)
		{
			int hid = FrameI.FHs[i].hs_net[j];

			for (int k = 0; k < 6; k++)
			{
				int fid = hex_mesh.HHs[hid].F_Ids[k];

				if (hex_mesh.HFs[fid].frame_boundary > -1)
				{
					bool have_already = false;
					for (int m = 0; m < FrameI.FHs[i].neighbor_FS.size(); m++)
					{
						if (FrameI.FHs[i].neighbor_FS[m] == hex_mesh.HFs[fid].frame_boundary)
							have_already = true;
					}
					if (!have_already)
					{
						FrameI.FHs[i].neighbor_FS.push_back(hex_mesh.HFs[fid].frame_boundary);
						FrameI.FFs[hex_mesh.HFs[fid].frame_boundary].neighbor_Cs.push_back(i);
					}
				}
			}
		}

		vector<int> neighbor_es;
		for (int j = 0; j < FrameI.FHs[i].neighbor_FS.size(); j++)
		{
			for (int k = 0; k < 4; k++)
			{
				int eid = FrameI.FFs[FrameI.FHs[i].neighbor_FS[j]].fe_Ids[k];
				bool added = false;
				for (int m = 0; m < neighbor_es.size(); m++)
					if (neighbor_es[m] == eid)
						added = true;
				if (!added)
					neighbor_es.push_back(eid);
			}
		}
		FrameI.FHs[i].neighbor_ES = neighbor_es;
		for (int j = 0; j < neighbor_es.size(); j++)
			FrameI.FEs[neighbor_es[j]].neighbor_Hs.push_back(i);

		vector<int> fh_8vs;

		int f1 = FrameI.FHs[i].neighbor_FS[0], f2 = -1;
		vector<int> vidsf1;
		vidsf1.push_back(FrameI.FFs[f1].fv_Ids[0]);
		vidsf1.push_back(FrameI.FFs[f1].fv_Ids[1]);
		vidsf1.push_back(FrameI.FFs[f1].fv_Ids[2]);
		vidsf1.push_back(FrameI.FFs[f1].fv_Ids[3]);
		for (int j = 1; j < FrameI.FHs[i].neighbor_FS.size(); j++)
		{
			f2 = FrameI.FHs[i].neighbor_FS[j];
			vector<int> results1;
			vector<int> vidsf2;
			vidsf2.push_back(FrameI.FFs[f2].fv_Ids[0]);
			vidsf2.push_back(FrameI.FFs[f2].fv_Ids[1]);
			vidsf2.push_back(FrameI.FFs[f2].fv_Ids[2]);
			vidsf2.push_back(FrameI.FFs[f2].fv_Ids[3]);
			set_cross(vidsf1, vidsf2, results1);
			if (!results1.size())
				break;
		}
		for (int j = 0; j < 8; j++)
			fh_8vs.push_back(-1);
		vector<int> vids;
		vids.push_back(FrameI.FFs[f2].fv_Ids[0]);
		vids.push_back(FrameI.FFs[f2].fv_Ids[1]);
		vids.push_back(FrameI.FFs[f2].fv_Ids[2]);
		vids.push_back(FrameI.FFs[f2].fv_Ids[3]);

		for (int j = 0; j < 4; j++)
		{
			fh_8vs[j] = FrameI.FFs[f1].fv_Ids[j];
			for (int k = 0; k < FrameI.FVs[FrameI.FFs[f1].fv_Ids[j]].neighbor_vs.size(); k++)
			{
				int which_vid = set_contain(vids, FrameI.FVs[FrameI.FFs[f1].fv_Ids[j]].neighbor_vs[k]);
				if (which_vid != -1)
				{
					fh_8vs[j + 4] = vids[which_vid];
					break;
				}
			}
		}

		for (int j = 0; j < fh_8vs.size(); j++)
		{
			FrameI.FHs[i].FV_Ids[j] = fh_8vs[j];
			FrameI.FVs[fh_8vs[j]].neighbor_Hs.push_back(i);
		}
	}

	for (int i = 0; i < FrameI.FVs.size(); i++)
	{
		set_redundent_clearn(FrameI.FVs[i].neighbor_vs);
		set_redundent_clearn(FrameI.FVs[i].neighbor_Es);
		set_redundent_clearn(FrameI.FVs[i].neighbor_Fs);
		set_redundent_clearn(FrameI.FVs[i].neighbor_Hs);
	}
	for (int i = 0; i < FrameI.FEs.size(); i++)
	{
		FrameI.FEs[i].is_boundary = -1;
		set_redundent_clearn(FrameI.FEs[i].neighbor_Fs);
		set_redundent_clearn(FrameI.FEs[i].neighbor_Hs);
	}
	for (int i = 0; i < FrameI.FFs.size(); i++)
	{
		set_redundent_clearn(FrameI.FFs[i].neighbor_Cs);
		FrameI.FFs[i].is_boundary = 0;
		if (FrameI.FFs[i].neighbor_Cs.size() == 1)
		{
			FrameI.FFs[i].is_boundary = 1;
			for (int j = 0; j < 4; j++)
			{
				int feid = FrameI.FFs[i].fe_Ids[j];
				FrameI.FEs[feid].is_boundary = 1;
			}
		}
	}
	for (int i = 0; i < FrameI.FHs.size(); i++)
	{
		set_redundent_clearn(FrameI.FHs[i].neighbor_ES);
		set_redundent_clearn(FrameI.FHs[i].neighbor_FS);
	}
}
void frame_of_mesh::extract_useless_hexes(vector<bool> &hexes_Inds)
{
	bool loop = true;

	vector<int> seeds;
	int start_hex = -1;
	for (int i = 0; i < hex_mesh.HFs.size(); i++)
	{
		if (hex_mesh.HFs[i].is_boundary == 1 && hex_mesh.HFs[i].frame_boundary == -1)
		{
			start_hex = hex_mesh.HFs[i].neighbor_Cs[0];
			break;
		}
	}
	if (start_hex == -1)
		loop = false;
	else
		seeds.push_back(start_hex);
	while (loop)
	{
		while (seeds.size() > 0)
		{
			int h_id = seeds[seeds.size() - 1];
			seeds.pop_back();
			hexes_Inds[h_id] = false;

			for (int i = 0; i < 6; i++)
			{
				int fid = hex_mesh.HHs[h_id].F_Ids[i];

				if (hex_mesh.HFs[fid].frame_boundary > -1)
					continue;

				for (int j = 0; j < hex_mesh.HFs[fid].neighbor_Cs.size(); j++)
				{
					int neighbor_hid = hex_mesh.HFs[fid].neighbor_Cs[j];
					if (neighbor_hid != h_id && hexes_Inds[neighbor_hid])
						seeds.push_back(neighbor_hid);
				}
			}
		}

		start_hex = -1;
		for (int i = 0; i < hex_mesh.HFs.size(); i++)
		{
			if (hex_mesh.HFs[i].is_boundary == 1 && hex_mesh.HFs[i].frame_boundary == -1 && hexes_Inds[hex_mesh.HFs[i].neighbor_Cs[0]])
			{
				start_hex = hex_mesh.HFs[i].neighbor_Cs[0];
				break;
			}
		}
		if (start_hex == -1)
			loop = false;
		else
			seeds.push_back(start_hex);
	}
}

void frame_of_mesh::assign_color_hexfaces_component()
{
	int Ind_C = 0;
	int Base_ = 8;
	for (int i = 0; i < FrameI.FHs.size(); i++)
	{
		if (FrameI.FHs[i].Color_ID == -1)
		{
			vector<int> ids;
			for (int j = 0; j < FrameI.FHs[i].neighbor_FS.size(); j++)
			{
				int fid = FrameI.FHs[i].neighbor_FS[j];
				for (int k = 0; k < FrameI.FFs[fid].neighbor_Cs.size(); k++)
				{
					int hid = FrameI.FFs[fid].neighbor_Cs[k];
					if (hid != i && FrameI.FHs[hid].Color_ID != -1)
					{
						ids.push_back(FrameI.FHs[hid].Color_ID);
					}
				}
			}
			while (true)
			{
				if (set_contain(ids, Ind_C) != -1)
				{
					Ind_C++;
					Ind_C = Ind_C % Base_;
				}
				else
					break;
			}
			FrameI.FHs[i].Color_ID = Ind_C;
		}
	}
	for (int i = 0; i < FrameI.FHs.size(); i++)
	{
		for (int j = 0; j < FrameI.FHs[i].hs_net.size(); j++)
		{
			int hid = FrameI.FHs[i].hs_net[j];
			hex_mesh.HHs[hid].Color_ID = FrameI.FHs[i].Color_ID;
		}
	}
}
frame_of_mesh::~frame_of_mesh(void)
{
}
