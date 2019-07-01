/*
 * ExtractHomotopyLoops.cpp
 *
 *  Created on: Dec 8, 2017
 *      Author: cotrik
 */

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "DualMesh.h"
#include "MST.h"
#include "BaseComplexQuad.h"
#include "BaseComplexSheetQuad.h"
#include "BaseComplexChord.h"
#include "ArgumentManager.h"
#include <iostream>
#include <fstream>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <stack>
#include <algorithm>
#include <iterator>

std::unordered_map<size_t, size_t> get_key_edgeid(const Mesh& mesh) {
	std::unordered_map<size_t, size_t> key_edgeid;
	for (auto& e : mesh.E) {
		key_edgeid[(e.Vids[0] << 32) | e.Vids[1]] = e.id;
		key_edgeid[(e.Vids[1] << 32) | e.Vids[0]] = e.id;
	}
	return key_edgeid;
}

////////////////////////////////////////////////////////////////
// C++ program to print all the cycles in an undirected graph 

std::vector<std::vector<size_t>> extractCycleVids(const Mesh& mesh, const std::vector<size_t>& eids) {
	CycleExtractor cycleExtractor(mesh, eids);
	return cycleExtractor.cycleVids;
}

std::vector<size_t> extractLoopEids(const Mesh& mesh, const std::vector<size_t>& eids) {
    std::vector<size_t> res = eids;
	Graph g;
	g.prune(mesh, res);
	return res;
}

std::vector<std::vector<size_t>> extractLoopEids(const Mesh& mesh, const std::vector<size_t>& mst_eids,
	const std::vector<size_t>& reversed_mst_eids) {
	std::vector<std::vector<size_t>> res;
	auto eids = mst_eids;
	for (auto eid : reversed_mst_eids) {
		eids.push_back(eid);
		res.push_back(extractLoopEids(mesh, eids));
		eids.pop_back();
	}
	std::sort(res.begin(), res.end(), [&](const std::vector<size_t>& a, const std::vector<size_t>& b) {
		return a.size() > b.size();
	});
	return res;
}

void GetMST(const Mesh& mesh, std::unordered_map<size_t, size_t>& key_edgeid,
	std::vector<size_t>& mst_eids, std::vector<size_t>& reversed_mst_eids) {
	Graph g(mesh.V.size(), mesh.E.size());
	for (auto& e : mesh.E)
		//if (!e.isBoundary)
		g.addEdge(e.Vids[0], e.Vids[1], 1);
	auto mst_wt = g.kruskalMST();
	// std::vector<size_t> mst_eids;
	for (auto& p : g.mst_edges) mst_eids.push_back(key_edgeid[(p.first << 32) + p.second]);
	std::set<size_t> mst_eids_set(mst_eids.begin(), mst_eids.end());
	//std::vector<size_t> reversed_mst_eids;
	for (auto& e : mesh.E)
		if (mst_eids_set.find(e.id) == mst_eids_set.end()) reversed_mst_eids.push_back(e.id);
}

std::vector<size_t> get_cut_graph_eids(const Mesh& mesh, const Mesh& dualMesh, const std::vector<size_t>& reversed_mst_eids) {
	std::vector<size_t> cut_graph_eids;
	for (auto eid : reversed_mst_eids) {
		auto& e = dualMesh.E.at(eid);
		size_t share_eid = MAXID;
		std::map<size_t, size_t> eid_ocunt;
		for (auto fid : e.Vids) {
			const auto& f = mesh.F.at(fid);
			for (auto eid : f.Eids)
				eid_ocunt[eid]++;
		}
		for (auto& item : eid_ocunt)
			if (item.second > 1) share_eid = item.first;
		if (share_eid == MAXID)	std::cerr << "share_eid Error\n";
		cut_graph_eids.push_back(share_eid);
	}
	return cut_graph_eids;
}

void get_cut_graph_mst_eids(const Mesh& mesh, std::unordered_map<size_t, size_t>& key_edgeIds, 
	const std::vector<size_t>& cut_graph_eids, std::vector<size_t>& cut_mst_eids, std::vector<size_t>& cut_reversed_mst_eids) {
	Graph cut_g(mesh.V.size(), cut_graph_eids.size());
	for (auto& eid : cut_graph_eids) {
		auto& e = mesh.E.at(eid);
		cut_g.addEdge(e.Vids[0], e.Vids[1], 1);
	}
	auto cut_mst_wt = cut_g.kruskalMST();
	// std::vector<size_t> cut_mst_eids;
	for (auto& p : cut_g.mst_edges) cut_mst_eids.push_back(key_edgeIds[(p.first << 32) + p.second]);

	std::set<size_t> cut_mst_eids_set(cut_mst_eids.begin(), cut_mst_eids.end());
	// std::vector<size_t> cut_reversed_mst_eids;
	for (auto& eid : cut_graph_eids)
		if (cut_mst_eids_set.find(eid) == cut_mst_eids_set.end()) cut_reversed_mst_eids.push_back(eid);
}

std::vector<size_t> GetLinkVids(const Mesh& mesh, const std::vector<size_t>& vids) {
	auto Vids = vids;
	std::vector<size_t> ringVids;
	ringVids.push_back(Vids.back());
	Vids.pop_back();
	while (ringVids.size() < vids.size()) {
		auto last_vid = ringVids.back();
		int i = 0;
		for (auto next_vid : Vids) {
			auto& last_v = mesh.V.at(last_vid);
			for (auto nvid : last_v.N_Vids)
				if (nvid == next_vid) {
					ringVids.push_back(next_vid);
					Vids.erase(Vids.begin() + i);
					break;
				}
			++i;
		}
	}
	return ringVids;
}

std::vector<size_t> GetLinkEVids(const Mesh& mesh, const std::vector<size_t>& eids) {
	std::set<size_t> vids;
	for (auto eid : eids) {
		const auto& e = mesh.E.at(eid);
		vids.insert(e.Vids.begin(), e.Vids.end());
	}
	std::vector<size_t> res;
	for (auto vid : vids) res.push_back(vid);
	return res;
}

std::vector<size_t> GetLinkVids_(const Mesh& mesh, const std::vector<size_t>& eids) {
	std::set<size_t> eids_set(eids.begin(), eids.end());
	const auto vids = GetLinkEVids(mesh, eids);
	auto Vids = vids;
	std::vector<size_t> ringVids;
	ringVids.push_back(Vids.back());
	Vids.pop_back();
	while (ringVids.size() <= vids.size()) {
		auto last_vid = ringVids.back();
		int i = 0;
		auto& last_v = mesh.V.at(last_vid);
		auto next_eid = MAXID;
		for (auto neid : last_v.N_Eids) {
			if (eids_set.find(neid) == eids_set.end()) continue;
			next_eid = neid;
			break;
		}
		if (next_eid == MAXID) {
			std::cerr << "Err in GetLinkVids_ ";
			break;
		}
		auto& ne = mesh.E.at(next_eid);
		auto next_vid = ne.Vids[0] == last_vid ? ne.Vids[1] : ne.Vids[0];
		ringVids.push_back(next_vid);
		eids_set.erase(next_eid);
		++i;
	}
	return ringVids;
}

std::vector<std::vector<size_t>> GetLinkVids(const Mesh& mesh, const std::vector<size_t>& cut_mst_eids, 
	const std::vector<size_t>& cut_reversed_mst_eids) {
	std::vector<std::vector<size_t>> res;
	auto eids = cut_mst_eids;
	for (auto eid : cut_reversed_mst_eids) {
		eids.push_back(eid);
		auto vids = GetLinkEVids(mesh, eids);
		res.push_back(GetLinkVids(mesh, vids));
		eids.pop_back();
	}
	return res;
}

std::vector<std::vector<size_t>> GetLinkEids(const Mesh& mesh, const std::vector<size_t>& cut_mst_eids,
	const std::vector<size_t>& cut_reversed_mst_eids) {
	std::vector<std::vector<size_t>> res;
	auto eids = cut_mst_eids;
	for (auto eid : cut_reversed_mst_eids) {
		eids.push_back(eid);
		res.push_back(eids);
		eids.pop_back();
	}
	return res;
}


void AddVertices(const Mesh& mesh, std::vector<std::vector<size_t>>& baseLinkVids) {
	adjacency_list_t adjacency_list(mesh.V.size());
	for (auto& v : mesh.V) {
		if (v.isSingularity) continue;
		for (auto nvid : v.N_Vids) {
			if (mesh.V.at(nvid).isSingularity) continue;
			adjacency_list[v.id].push_back(neighbor(nvid, 1));
		}
	}
	auto linkid = 0;
	for (auto& link : baseLinkVids) {
		std::vector<size_t> new_link;
		while (link.size() > 2 && link.front() != link.back()) link.erase(link.begin());
		if (link.front() != link.back()) {
			std::cerr << "Err in AddVertices, linkid = " << linkid++ << "\n";
			continue;
			//break;
		}
		new_link.push_back(link.front());
		for (int i = 1; i < link.size(); ++i) {
			auto& v = mesh.V.at(link[i]);
			std::set<size_t> nvids(v.N_Vids.begin(), v.N_Vids.end());
			if (nvids.find(link[i - 1]) == nvids.end()) {
				auto src = link[i - 1];
				auto dest = link[i];
				std::vector<weight_t> min_distance;
				std::vector<vertex_t> previous;
				Graph g;
				g.DijkstraComputePaths(src, adjacency_list, min_distance, previous);
				std::list<vertex_t> path = g.DijkstraGetShortestPathTo(dest, previous);
				std::vector<size_t> path_vids(path.begin(), path.end());
				for (int j = 1; j < path_vids.size(); ++j) {
					new_link.push_back(path_vids[j]);
				}
			} else {
				new_link.push_back(v.id);
			}
		}
		// new_link.push_back(link.back());
		link = new_link;
	}
}

std::vector<size_t> GetShortestPath(const adjacency_list_t& adjacency_list, size_t src, size_t dest) {
	std::vector<weight_t> min_distance;
	std::vector<vertex_t> previous;
	Graph g;
	g.DijkstraComputePaths(src, adjacency_list, min_distance, previous);
	std::list<vertex_t> path = g.DijkstraGetShortestPathTo(dest, previous);
	std::vector<size_t> path_vids(path.begin(), path.end());
	return path_vids;
}

void RemoveSingularities(const Mesh& mesh, std::vector<std::vector<size_t>>& baseLinkVids) {
	adjacency_list_t adjacency_list(mesh.V.size());
	for (auto& v : mesh.V) {
		if (v.isSingularity) continue;
		for (auto nvid : v.N_Vids) {
			if (mesh.V.at(nvid).isSingularity) continue;
			adjacency_list[v.id].push_back(neighbor(nvid, 1));
		}
	}
	auto linkid = 0;
	for (auto& link : baseLinkVids) {
		std::vector<size_t> new_link;
		while (link.size() > 2 && link.front() != link.back()) link.erase(link.begin());
		if (link.front() != link.back()) {
			std::cerr << "Err in RemoveSingularities, linkid = " << linkid++ << "\n";
			continue;
			//break;
		}
		new_link.push_back(link.front());
		for (int i = 1; i < link.size() - 1; ++i) {
			auto& v = mesh.V.at(link[i]);
			if (v.isSingularity) {
				auto src = link[i - 1];
				while (i + 1 < link.size() && mesh.V.at(link[i + 1]).isSingularity) {
					++i;
				}
				if (i + 1 < link.size()) {
					auto dest = link[i + 1];
					std::vector<size_t> path_vids = GetShortestPath(adjacency_list, src, dest);
					for (int j = 1; j < path_vids.size() - 1; ++j) {
						new_link.push_back(path_vids[j]);
					}
				}
			} else {
				new_link.push_back(v.id);
			}
		}
		new_link.push_back(link.back());

		auto& v = mesh.V.at(new_link.front());
		if (v.isSingularity) {
			new_link.pop_back();
			new_link.erase(new_link.begin());

			auto src = new_link.back();
			auto dest = new_link.front();
			std::vector<size_t> path_vids = GetShortestPath(adjacency_list, src, dest);
			for (int j = 1; j < path_vids.size() - 1; ++j) {
				new_link.push_back(path_vids[j]);
			}

			new_link.push_back(dest);
		}

		link = new_link;
	}
}

int GetHolonomyGroup(const Mesh& dualMesh, const std::vector<size_t>& baseLinkVids) {
	int res = 0;
	auto& fids = baseLinkVids;
	auto eid0 = dualMesh.F.at(fids.front()).Eids[0];
	{
		auto& prev_f = dualMesh.F.at(fids.at(0));
		auto& curr_f = dualMesh.F.at(fids.at(1));
		std::set<size_t> eids1(prev_f.Eids.begin(), prev_f.Eids.end());
		std::set<size_t> eids2(curr_f.Eids.begin(), curr_f.Eids.end());
		std::set<size_t> eids = Util::get_intersect(eids1, eids2);
		if (eids.size() != 1) {
			std::cerr << "Err in GetHolonomyGroup i = " << 0 << std::endl;
		}
		eid0 = *eids.begin();
	}
	auto parallel_eid = eid0;
	for (auto i = 1; i < fids.size(); ++i) {
		auto& prev_f = dualMesh.F.at(fids.at(i - 1));
		auto& curr_f = dualMesh.F.at(fids.at(i));
		std::set<size_t> eids1(prev_f.Eids.begin(), prev_f.Eids.end());
		std::set<size_t> eids2(curr_f.Eids.begin(), curr_f.Eids.end());
		std::set<size_t> eids = Util::get_intersect(eids1, eids2);
		if (eids.size() != 1) {
			std::cerr << "Err in GetHolonomyGroup i = " << i << std::endl;
			continue;
		}
		if (curr_f.Vids.size() != 4) {
			std::cerr << "Err in GetHolonomyGroup curr_f.Vids.size() = " << curr_f.Vids.size() << " fid = " << curr_f.id << std::endl;
		}
		auto share_eid = *eids.begin();
		auto& share_e = dualMesh.E.at(share_eid);
		if (share_eid == parallel_eid) {
			std::set<size_t> eids1(share_e.parallelEids.begin(), share_e.parallelEids.end());
			std::set<size_t> eids2(curr_f.Eids.begin(), curr_f.Eids.end());
			std::set<size_t> eids = Util::get_intersect(eids1, eids2);
			if (eids.size() != 1) {
				std::cerr << "Err in GetHolonomyGroup parallel_eids, share_eid = " << share_eid << std::endl;
				//continue;
			}
			parallel_eid = *eids.begin();
		} else if (share_e.parallelEids[0] == parallel_eid || share_e.parallelEids[1] == parallel_eid) {
			parallel_eid = share_eid;
		} else {
			auto& parallel_e = dualMesh.E.at(parallel_eid);
			//std::set<size_t> eids1(parallel_e.consecutiveEids.begin(), parallel_e.consecutiveEids.end());
			//std::set<size_t> eids2(curr_f.Eids.begin(), curr_f.Eids.end());
			//std::set<size_t> eids = Util::get_intersect(eids1, eids2);
			//if (eids.size() != 1) {
			//	std::cerr << "Err in GetHolonomyGroup consecutive_eids, share_eid = " << share_eid << std::endl;
			//	continue;
			//}
			//parallel_eid = *eids.begin();

			const auto parallel_eid_old = parallel_eid;
			for (auto eid : curr_f.Eids) {
				if (eid == share_eid) continue;
				auto& e = dualMesh.E.at(eid);
				if (e.Vids[0] == parallel_e.Vids[0] || e.Vids[1] == parallel_e.Vids[0] ||
					e.Vids[0] == parallel_e.Vids[1] || e.Vids[1] == parallel_e.Vids[1]) {
					parallel_eid = eid;
					break;
				}
			}
			if (parallel_eid_old == parallel_eid) {
				std::cerr << "Err in GetHolonomyGroup consecutive_eids, share_eid = " << share_eid
					<< " (" << share_e.Vids[0] << ", " << share_e.Vids[1]
					<< ") share_e.parallel_eids = (" << share_e.parallelEids[0] << ", " << share_e.parallelEids[1]
					<< ") parallel_eid = " << parallel_eid << " (" << parallel_e.Vids[0] << ", " << parallel_e.Vids[1] << ") i = " << i
					<< " prev_fid = " << prev_f.id << " curr_fid = " << curr_f.id << std::endl;
				//continue;
				return -1;
			}
		}
		// std::cout << "i = " << i << " parallel_eid = " << parallel_eid << std::endl;
	}
	auto& e0 = dualMesh.E.at(eid0);
	if (parallel_eid == eid0) res = 0;
	else if (e0.parallelEids[0] == parallel_eid || e0.parallelEids[1] == parallel_eid) res = 2;
	else {
		auto& parallel_e = dualMesh.E.at(parallel_eid);
		std::set<size_t> vids1(parallel_e.Vids.begin(), parallel_e.Vids.end());
		std::set<size_t> vids2(e0.Vids.begin(), e0.Vids.end());
		std::set<size_t> vids = Util::get_intersect(vids1, vids2);
		if (vids.size() != 1) {
			std::cerr << "Err in GetHolonomyGroup vids" << std::endl;
		}
		auto vid = *vids.begin();
		auto& v = dualMesh.V.at(vid);
		auto vid0 = e0.Vids[0] == vid ? e0.Vids[1] : e0.Vids[0];
		auto vid1 = parallel_e.Vids[0] == vid ? parallel_e.Vids[1] : parallel_e.Vids[0];
		auto& v0 = dualMesh.V.at(vid0);
		auto& v1 = dualMesh.V.at(vid1);
		auto dir = glm::cross(v1 - v, v0 - v);

		auto& f0 = dualMesh.F.at(fids.front());
		auto& v_ = dualMesh.V.at(f0.Vids[1]);
		auto& v_0 = dualMesh.V.at(f0.Vids[0]);
		auto& v_1 = dualMesh.V.at(f0.Vids[2]);
		auto dir0 = glm::cross(v_1 - v_, v_0 - v_);

		auto sign = glm::dot(dir, dir0);
		if (sign > 0) res = 1;
		else res = 3;
	}
	return res;
}

std::set<int> GetHolonomyGroup(const Mesh& dualMesh, const std::vector<std::vector<size_t>>& baseLinkVids) {
	// get_key_dualEdgeId(dualMesh);
	std::set<int> res;
	int linkid = 0;
	for (auto& linkVids : baseLinkVids) {
		std::cout << "Processing link " << linkid++ << std::endl;
		auto holonomy = GetHolonomyGroup(dualMesh, linkVids);
		res.insert(holonomy);
	}
	return res;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: ExtractHomotopyLoops <file> <out.vtk>";
        return -1;
    }
    ArgumentManager argumentManager(argc, argv);
	MeshFileReader reader(argv[1]);
	auto mesh = reader.GetMesh();
	mesh.BuildAllConnectivities();
	mesh.ExtractBoundary();
	mesh.ExtractSingularities();
	mesh.BuildParallelE();
	mesh.BuildConsecutiveE();
	mesh.BuildOrthogonalE();

	DualMesh dualHoleMesh;
	dualHoleMesh.Build(mesh);
	std::vector<size_t> quadCids;
	for (auto& c : dualHoleMesh.C)
		if (c.Vids.size() == 4) quadCids.push_back(c.id);
	std::vector<Cell> newC(quadCids.size());
	for (auto i = 0; i < newC.size(); ++i) {
		newC[i].id = i;
		newC[i].Vids = dualHoleMesh.C.at(quadCids[i]).Vids;
		newC[i].cellType = VTK_QUAD;
	}
	dualHoleMesh.C = newC;
	{
		MeshFileWriter writer(dualHoleMesh, "holes.vtk");
		writer.WriteFile();
	}
	dualHoleMesh.BuildAllConnectivities();
	dualHoleMesh.ExtractBoundary();
	dualHoleMesh.ExtractSingularities();


	DualMesh dualMesh;
	dualMesh.Build(mesh);
	dualMesh.BuildAllConnectivities();
	dualMesh.ExtractBoundary();
	dualMesh.ExtractSingularities();
	dualMesh.BuildParallelE();

	{
		MeshFileWriter writer(dualMesh, "dual_edges.vtk");
		writer.WriteEdgesVtk();
	}

	//auto& dualMesh = dualHoleMesh;

	auto key_edgeIds = get_key_edgeid(mesh);
	auto key_dualEdgeIds = get_key_edgeid(dualMesh);
	std::vector<size_t> mst_eids;
	std::vector<size_t> reversed_mst_eids;
	GetMST(dualMesh, key_dualEdgeIds, mst_eids, reversed_mst_eids);
	{
		MeshFileWriter writer(dualMesh, "dual_mst.vtk");
		writer.WriteEdgesVtk(mst_eids);
	}
	{
		MeshFileWriter writer(dualMesh, "dual_rmst.vtk");
		writer.WriteEdgesVtk(reversed_mst_eids);
	}
	auto cut_graph_eids = get_cut_graph_eids(mesh, dualMesh, reversed_mst_eids);
	auto cut_graph_eids_without_prune = cut_graph_eids;
	{
		MeshFileWriter writer(mesh, "CutGraph0.vtk");
		writer.WriteEdgesVtk(cut_graph_eids_without_prune);
	}
	Graph g;
	g.prune(mesh, cut_graph_eids);
	{
		std::cout << "Writing CutGraph.vtk\n";
		MeshFileWriter writer(mesh, "CutGraph.vtk");
		writer.WriteEdgesVtk(cut_graph_eids);
	}

	std::vector<size_t> cut_mst_eids;
	std::vector<size_t> cut_reversed_mst_eids;
	get_cut_graph_mst_eids(mesh, key_edgeIds, cut_graph_eids, cut_mst_eids, cut_reversed_mst_eids);
	{
		std::cout << "Writing CutGraphMST.vtk\n";
		MeshFileWriter writer(mesh, "CutGraphMST.vtk");
		writer.WriteEdgesVtk(cut_mst_eids);
	}
	{
		std::cout << "Writing CutGraphRMST.vtk\n";
		MeshFileWriter writer(mesh, "CutGraphRMST.vtk");
		writer.WriteEdgesVtk(cut_reversed_mst_eids);
	}

	{
		auto linkEids = extractLoopEids(mesh, cut_mst_eids, cut_reversed_mst_eids);
		MeshFileWriter writer(mesh, argv[2]);
		writer.WriteEdgesVtk(linkEids);
		std::vector<std::vector<size_t>> linkVids;
		for (auto& eids : linkEids)
			linkVids.push_back(GetLinkVids_(mesh, eids));
		{
			MeshFileWriter writer(mesh, "baseLinkVids.vtk");
			writer.WriteLinksVtk(linkVids);
		}

		//------------------------------------------
		BaseComplexQuad baseComplex(mesh);
		baseComplex.Build();
		BaseComplexSheetQuad baseComplexSheets(baseComplex);
		baseComplexSheets.Extract();

		RefinedDualQuad dual(mesh);
		dual.Build();

		//std::vector<std::unordered_set<size_t>> all_dualEdgeIds;
		//size_t numOfLines = 0;
		//for (size_t i = 0; i < baseComplexSheets.sheets_componentFaceIds.size(); ++i) {
		//	std::unordered_set<size_t> dualEdgeIds = baseComplexSheets.GetDualEdgeIds(i);
		//	all_dualEdgeIds.push_back(dualEdgeIds);
		//	numOfLines += dualEdgeIds.size();
		//}
		//for (auto x = 0; x < baseComplexSheets.sheets_componentFaceIds.size(); ++x) 

		std::vector<std::vector<size_t>> intersects;
		auto x = 0;
		{
			std::cout << "\n************************************\n";
			std::cout << "*************** x = " << x << std::endl;
			std::unordered_set<size_t> dualEdgeIds = baseComplexSheets.GetDualEdgeIds(x);
			std::set<size_t> dualVids;
			for (auto eid : dualEdgeIds) {
				auto& e = dual.E.at(eid);
				dualVids.insert(e.Vids.begin(), e.Vids.end());
			}

			std::set<size_t> dualFids;
			for (auto vid : dualVids) {
				if (vid >= mesh.V.size() + mesh.E.size()) dualFids.insert(vid - mesh.V.size() - mesh.E.size());
			}
			//std::vector<size_t> dualVids_;
			//std::set<size_t> dualEids;
			//for (auto vid : dualVids) {
			//	if (vid < mesh.V.size() + mesh.E.size() && vid >= mesh.V.size()) {
			//		auto eid = vid - mesh.V.size();
			//		auto& e = mesh.E.at(eid);
			//		dualEids.insert(eid);
			//		if (dualVids_.empty()) {
			//			dualVids_.push_back(e.Vids[0]);
			//		} else {
			//			auto last_vid = dualVids_.back();
			//			auto& v0 = mesh.V.at(e.Vids[0]);
			//			auto& v1 = mesh.V.at(e.Vids[1]);
			//			for (auto nvid : v0.N_Vids) {
			//				if (nvid == last_vid) {
			//					dualVids_.push_back(e.Vids[0]);
			//					break;
			//				}
			//			}
			//			if (last_vid == dualVids_.back())
			//			for (auto nvid : v1.N_Vids) {
			//				if (nvid == last_vid) {
			//					dualVids_.push_back(e.Vids[1]);
			//					break;
			//				}
			//			}
			//		}
			//	}
			//	dualVids.clear();
			//	dualVids.insert(dualVids_.begin(), dualVids_.end());
			//}
			std::cout << "\n-----------------------------\n";
			std::cout << "dualVids size = " << dualVids.size() << "\n";
			auto loopid = 0;
			for (auto& vids : linkVids) {
				std::set<size_t> intersect;
				std::set<size_t> vids_set(vids.begin(), vids.end());
				for (auto fid : dualFids) {
					for (auto vid : mesh.F.at(fid).Vids) {
						if (vids_set.find(vid) != vids_set.end()) {
							intersect.insert(vid);
							break;
						}
					}
				}
				if (!intersect.empty()) {
					std::vector<size_t> y;
					for (auto z : intersect)
						y.push_back(z);
					intersects.push_back(y);
				}
				std::cout << " base loop " << loopid++ << " : " << intersect.size() << std::endl;
				//auto intersect = Util::get_intersect(dualVids, vids_set);
				//std::cout << "loop size = " << vids_set.size() << " base loop " << loopid++ << " : " << 
				//	intersect.size() << std::endl;
			}
			std::cout << "\n-----------------------------\n";
		}
		{
			MeshFileWriter writer(mesh, "intersects.vtk");
			writer.WriteVerticesVtk(intersects);
		}
		//------------------------------------------

		RemoveSingularities(mesh, linkVids);
		// AddVertices(mesh, linkVids);
		g.PruneLinkVids(mesh, linkVids);
		{
			MeshFileWriter writer(mesh, "cleanBaseLinkVids.vtk");
			writer.WriteLinksVtk(linkVids);
		}
		auto holonomyGroup = GetHolonomyGroup(dualMesh, linkVids);
		std::cout << "HolonomyGroup:";
		for (auto holonomy : holonomyGroup)
			std::cout << " " << holonomy;
		std::cout << std::endl;
		{
			MeshFileWriter writer(dualMesh, "face_path.vtk");
			writer.WriteFacesVtk(linkVids);
		}
	}


    return 0;
}
