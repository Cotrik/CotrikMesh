/*
* ExtractCoverSpace2.cpp
*
*  Created on: Dec 18, 2017
*      Author: cotrik
*/

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "DualMesh.h"
#include "MST.h"
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
	std::set<int> res;
	int linkid = 0;
	for (auto& linkVids : baseLinkVids) {
		std::cout << "Processing link " << linkid++ << std::endl;
		auto holonomy = GetHolonomyGroup(dualMesh, linkVids);
		res.insert(holonomy);
	}
	return res;
}

const int connections[4][4][2] = {
{ { 0, 0 },{ 1, 1 },{ 2, 2 },{ 3, 3 } },
{ { 0, 3 },{ 1, 0 },{ 2, 1 },{ 3, 2 } },
{ { 0, 2 },{ 1, 3 },{ 2, 0 },{ 3, 1 } },
{ { 0, 1 },{ 1, 2 },{ 2, 3 },{ 3, 0 } }
};

void initial(const Mesh& mesh, std::vector<std::vector<Vertex>>& V, std::vector<std::vector<Cell>>& F,
	std::unordered_map<size_t, size_t>& g_coverVid_origVid, std::unordered_map<size_t, std::set<size_t>>& g_coverVid_coverFids) {
	for (auto i = 0; i < 4; ++i) {
		for (auto& f : mesh.F) {
			auto fid = f.id + i * mesh.F.size();
			auto& fi = F[i][f.id];
			fi.id = fid;
			fi.Vids.resize(4);
			glm::dvec3 center;
			for (auto vid : f.Vids)
				center += mesh.V.at(vid).xyz();
			center *= 0.25;
			for (auto j = 0; j < 4; ++j) {
				fi.Vids[j] = 4 * fid + j;
				g_coverVid_origVid[fi.Vids[j]] = f.Vids[j];
				g_coverVid_coverFids[fi.Vids[j]].insert(fi.id);

				auto& vi = V[i][f.id * 4 + j];
				auto& v = mesh.V.at(f.Vids[j]);
				vi = (v.xyz() - center) * (0.5 + 0.1 * i) + center;
				vi.id = fi.Vids[j];
			}
		}
	}
}

void initial_1(const Mesh& mesh, std::vector<std::vector<Vertex>>& V, std::vector<std::vector<Cell>>& F,
	std::unordered_map<size_t, size_t>& g_coverVid_origVid, std::unordered_map<size_t, std::set<size_t>>& g_coverVid_coverFids) {

	//glm::dvec3 center;
	//for (auto& v : mesh.V)
	//	center += v.xyz();
	//center /= mesh.V.size();
	auto& e = mesh.E.at(0);
	auto& v0 = mesh.V.at(e.Vids[0]);
	auto& v1 = mesh.V.at(e.Vids[1]);
	auto length = glm::length(v0 - v1);
	for (auto i = 0; i < 4; ++i) {
		for (auto& f : mesh.F) {
			auto fid = f.id + i * mesh.F.size();
			auto& fi = F[i][f.id];
			fi.id = fid;
			fi.Vids.resize(4);
			for (auto j = 0; j < 4; ++j) {
				fi.Vids[j] = 4 * fid + j;
				g_coverVid_origVid[fi.Vids[j]] = f.Vids[j];
				g_coverVid_coverFids[fi.Vids[j]].insert(fi.id);
				glm::dvec3 center;
				for (auto vid : f.Vids)
					center += mesh.V.at(vid).xyz();
				center *= 0.25;
				auto& vi = V[i][f.id * 4 + j];
				auto& v = mesh.V.at(f.Vids[j]);
				vi = (v.xyz() - center) * (0.5 + 0.1 * i) + center + f.normal * length * (1.0 - 0.05 * i);
				vi.id = fi.Vids[j];
			}
		}
	}
}

size_t get_cover_vid(const std::unordered_map<size_t, size_t>& g_coverVid_origVid, const Cell& cover_f, const size_t orig_vid) {
	auto cover_vid = MAXID;
	for (auto vid : cover_f.Vids) {
		auto& m = g_coverVid_origVid;
		auto iter = m.find(vid);
		if (iter != m.end() && iter->second == orig_vid) { 
			cover_vid = vid; 
			break; 
		}
	}
	return cover_vid;
}

void update(const Mesh& mesh, std::vector<std::vector<Cell>>& F, size_t cover_f0_share_vid0, size_t cover_f1_share_vid0,
	std::unordered_map<size_t, std::set<size_t>>& g_coverVid_coverFids) {
	auto& nfids0 = g_coverVid_coverFids[cover_f0_share_vid0];
	for (auto fid : nfids0) {
		auto copy_ = fid / mesh.F.size();
		auto id = fid % mesh.F.size();
		auto& f0 = F[copy_][id];
		for (auto& vid : f0.Vids)
			if (vid == cover_f0_share_vid0 && vid != cover_f1_share_vid0) {
				vid = cover_f1_share_vid0;
				g_coverVid_coverFids[vid].insert(nfids0.begin(), nfids0.end());
				break;
			}
	}
}

Mesh GetBranchCover0(const Mesh& mesh, const DualMesh& dualMesh) {
	std::vector<std::vector<Vertex>> V(4, std::vector<Vertex>(4 * mesh.F.size())); // 4 covers, each cover has C faces, each cell has 4 vertices;
	std::vector<std::vector<Cell>>   F(4, std::vector<Cell>(mesh.F.size()));

	std::unordered_map<size_t, size_t> g_coverVid_origVid;
	std::unordered_map<size_t, std::set<size_t>> g_coverVid_coverFids;

	initial(mesh, V, F, g_coverVid_origVid, g_coverVid_coverFids);

	std::vector<std::vector<size_t>> loops;
	std::vector<int> holonomys;

	size_t orig = 0;
	adjacency_list_t adjacency_list(dualMesh.V.size());
	for (auto& v : dualMesh.V)
		for (auto nvid : v.N_Vids)
			adjacency_list[v.id].push_back(neighbor(nvid, 1));
	for (auto& e : mesh.E) {
		auto vid0 = e.N_Fids[0];
		auto vid1 = e.N_Fids[1];
		std::vector<size_t> linkVids = GetShortestPath(adjacency_list, orig, vid0);
		std::vector<size_t> linkVids1 = GetShortestPath(adjacency_list, vid1, orig);
		std::copy(linkVids1.begin(), linkVids1.end(), std::back_inserter(linkVids));
		loops.push_back(linkVids);
		auto holonomy = GetHolonomyGroup(mesh, linkVids);
		holonomys.push_back(holonomy);
		for (auto i = 0; i < 4; ++i) {
			auto copy_0 = connections[holonomy][i][0];
			auto copy_1 = connections[holonomy][i][1];
			auto& cover_f0 = F[copy_0][e.N_Fids[0]];
			auto& cover_f1 = F[copy_1][e.N_Fids[1]];

			auto& orig_v0 = mesh.E.at(e.Vids[0]);
			auto& orig_v1 = mesh.E.at(e.Vids[1]);

			auto cover_f0_share_vid0 = get_cover_vid(g_coverVid_origVid, cover_f0, orig_v0.id);
			auto cover_f0_share_vid1 = get_cover_vid(g_coverVid_origVid, cover_f0, orig_v1.id);
			auto cover_f1_share_vid0 = get_cover_vid(g_coverVid_origVid, cover_f1, orig_v0.id);
			auto cover_f1_share_vid1 = get_cover_vid(g_coverVid_origVid, cover_f1, orig_v1.id);

			if (cover_f0_share_vid0 == MAXID || cover_f0_share_vid1 == MAXID ||
				cover_f1_share_vid0 == MAXID || cover_f1_share_vid1 == MAXID) {
				std::cerr << "Err in GetBranchCover eid = " << e.id << " holonomy = " << holonomy
					<< " ~(" << copy_0 << ", " << copy_1 << ")\n";
			}

			update(mesh, F, cover_f0_share_vid0, cover_f1_share_vid0, g_coverVid_coverFids);
			update(mesh, F, cover_f0_share_vid1, cover_f1_share_vid1, g_coverVid_coverFids);
		}
		auto base = mesh.E.size() / 10;
		if (e.id % base == 0) std::cout << "Done " << float(e.id) * 100 / mesh.E.size() << "%\n";
	}
	{
		MeshFileWriter writer(mesh, "FacePaths.vtk");
		writer.WriteFacesVtk(loops);
		std::vector<int> newholonomys;
		for (auto i = 0; i < loops.size(); ++i) {
			auto& loop = loops[i];
			for (auto fid : loop) newholonomys.push_back(holonomys[i]);
		}
		writer.WriteCellData(newholonomys, "holonomy", false);
	}
	{
		MeshFileWriter writer(dualMesh, "DualLoops.vtk");
		writer.WriteLinksVtk(loops);
		writer.WriteCellData(holonomys, "holonomy", false);
	}
	{
		MeshFileWriter writer(mesh, "EdgeHolonomy.vtk");
		writer.WriteEdgesVtk();
		writer.WriteCellData(holonomys, "holonomy", true);
	}
	std::cout << "Building CoverMesh V\n";
	std::vector<Vertex> coverV(4 * mesh.F.size() * 4);
	std::vector<Cell>   coverC(4 * mesh.F.size());
	for (auto i = 0; i < 4; ++i) {
		for (size_t j = 0; j < V[i].size(); ++j) {
			auto& v = V[i][j];
			auto& coverv = coverV[v.id];
			coverv.id = v.id;
			coverv = v;
		}
	}
	std::cout << "Building CoverMesh C\n";
	for (auto i = 0; i < 4; ++i) {
		for (size_t j = 0; j < F[i].size(); ++j) {
			auto& c = F[i][j];
			auto& coverc = coverC[c.id];
			coverc.id = c.id;
			coverc.Vids = c.Vids;
		}
	}
	std::cout << "Finish Building CoverMesh\n";
	Mesh coverMesh(coverV, coverC, QUAD);
	return coverMesh;
}

Mesh GetBranchCover1(const Mesh& mesh, const DualMesh& dualMesh, const std::set<size_t>& scut) {
	std::vector<std::vector<Vertex>> V(4, std::vector<Vertex>(4 * mesh.F.size())); // 4 covers, each cover has C faces, each cell has 4 vertices;
	std::vector<std::vector<Cell>>   F(4, std::vector<Cell>(mesh.F.size()));

	std::unordered_map<size_t, size_t> g_coverVid_origVid;
	std::unordered_map<size_t, std::set<size_t>> g_coverVid_coverFids;

	initial(mesh, V, F, g_coverVid_origVid, g_coverVid_coverFids);

	std::vector<std::vector<size_t>> loops;
	std::vector<int> holonomys;

	std::unordered_map<size_t, size_t> scut_key_eids;
	for (auto eid : scut) {
		auto& e = mesh.E.at(eid);
		scut_key_eids[(e.N_Fids[0] << 32) | e.N_Fids[1]] = eid;
		scut_key_eids[(e.N_Fids[1] << 32) | e.N_Fids[0]] = eid;
	}
	std::vector<size_t> dualEids;
	for (auto& v : dualMesh.V)
		for (auto neid : v.N_Eids) {
			auto& ne = dualMesh.E.at(neid);
			auto nvid = ne.Vids[0] == v.id ? ne.Vids[1] : ne.Vids[0];
			auto key = (nvid << 32) | v.id;
			auto iter = scut_key_eids.find(key);
			if (iter == scut_key_eids.end()) continue;
			dualEids.push_back(neid);
		}
	{
		std::cout << "Writing SCutGraph.vtk\n";
		MeshFileWriter writer(dualMesh, "SCutGraph.vtk");
		writer.WriteEdgesVtk(dualEids);
	}
	size_t orig = 0;
	adjacency_list_t adjacency_list(dualMesh.V.size());
	for (auto& v : dualMesh.V)
		for (auto nvid : v.N_Vids) {
			auto key = (nvid << 32) | v.id;
			auto iter = scut_key_eids.find(key);
			if (iter != scut_key_eids.end()) continue;
			adjacency_list[v.id].push_back(neighbor(nvid, 1));
		}
	for (auto& e : mesh.E) {
		auto vid0 = e.N_Fids[0];
		auto vid1 = e.N_Fids[1];
		std::vector<size_t> linkVids = GetShortestPath(adjacency_list, orig, vid0);
		std::vector<size_t> linkVids1 = GetShortestPath(adjacency_list, vid1, orig);
		std::copy(linkVids1.begin(), linkVids1.end(), std::back_inserter(linkVids));
		loops.push_back(linkVids);
		auto holonomy = GetHolonomyGroup(mesh, linkVids);
		holonomys.push_back(holonomy);
		for (auto i = 0; i < 4; ++i) {
			auto copy_0 = connections[holonomy][i][0];
			auto copy_1 = connections[holonomy][i][1];
			auto& cover_f0 = F[copy_0][e.N_Fids[0]];
			auto& cover_f1 = F[copy_1][e.N_Fids[1]];

			auto& orig_v0 = mesh.E.at(e.Vids[0]);
			auto& orig_v1 = mesh.E.at(e.Vids[1]);

			auto cover_f0_share_vid0 = get_cover_vid(g_coverVid_origVid, cover_f0, orig_v0.id);
			auto cover_f0_share_vid1 = get_cover_vid(g_coverVid_origVid, cover_f0, orig_v1.id);
			auto cover_f1_share_vid0 = get_cover_vid(g_coverVid_origVid, cover_f1, orig_v0.id);
			auto cover_f1_share_vid1 = get_cover_vid(g_coverVid_origVid, cover_f1, orig_v1.id);

			if (cover_f0_share_vid0 == MAXID || cover_f0_share_vid1 == MAXID ||
				cover_f1_share_vid0 == MAXID || cover_f1_share_vid1 == MAXID) {
				std::cerr << "Err in GetBranchCover eid = " << e.id << " holonomy = " << holonomy
					<< " ~(" << copy_0 << ", " << copy_1 << ")\n";
			}

			update(mesh, F, cover_f0_share_vid0, cover_f1_share_vid0, g_coverVid_coverFids);
			update(mesh, F, cover_f0_share_vid1, cover_f1_share_vid1, g_coverVid_coverFids);
		}
		auto base = mesh.E.size() / 10;
		if (e.id % base == 0) std::cout << "Done " << float(e.id) * 100 / mesh.E.size() << "%\n";
	}
	{
		MeshFileWriter writer(mesh, "FacePaths.vtk");
		writer.WriteFacesVtk(loops);
		std::vector<int> newholonomys;
		for (auto i = 0; i < loops.size(); ++i) {
			auto& loop = loops[i];
			for (auto fid : loop) newholonomys.push_back(holonomys[i]);
		}
		writer.WriteCellData(newholonomys, "holonomy", false);
	}
	{
		MeshFileWriter writer(dualMesh, "DualLoops.vtk");
		writer.WriteLinksVtk(loops);
		writer.WriteCellData(holonomys, "holonomy", false);
	}
	{
		MeshFileWriter writer(mesh, "EdgeHolonomy.vtk");
		writer.WriteEdgesVtk();
		writer.WriteCellData(holonomys, "holonomy", true);
	}
	std::cout << "Building CoverMesh V\n";
	std::vector<Vertex> coverV(4 * mesh.F.size() * 4);
	std::vector<Cell>   coverC(4 * mesh.F.size());
	for (auto i = 0; i < 4; ++i) {
		for (size_t j = 0; j < V[i].size(); ++j) {
			auto& v = V[i][j];
			auto& coverv = coverV[v.id];
			coverv.id = v.id;
			coverv = v;
		}
	}
	std::cout << "Building CoverMesh C\n";
	for (auto i = 0; i < 4; ++i) {
		for (size_t j = 0; j < F[i].size(); ++j) {
			auto& c = F[i][j];
			auto& coverc = coverC[c.id];
			coverc.id = c.id;
			coverc.Vids = c.Vids;
		}
	}
	std::cout << "Finish Building CoverMesh\n";
	Mesh coverMesh(coverV, coverC, QUAD);
	return coverMesh;
}

std::vector<std::vector<size_t>> get_shortest_paths(const Mesh& mesh, const DualMesh& dualMesh) {
	size_t orig = 0;
	adjacency_list_t adjacency_list(dualMesh.V.size());
	for (auto& v : dualMesh.V)
		for (auto nvid : v.N_Vids) {
			adjacency_list[v.id].push_back(neighbor(nvid, 1));
		}
	std::vector<std::vector<size_t>> shortest_paths;
	for (auto& c : mesh.C)
		shortest_paths.push_back(GetShortestPath(adjacency_list, orig, c.id));

	return shortest_paths;
}

std::vector<std::vector<size_t>> get_loops(const Mesh& mesh, const DualMesh& dualMesh) {
	std::vector<std::vector<size_t>> loops;
	size_t orig = 0;
	auto shortest_paths = get_shortest_paths(mesh, dualMesh);
	for (auto& f : mesh.F) {
		if (f.isBoundary) continue;
		auto vid0 = f.N_Cids[0];
		auto vid1 = f.N_Cids[1];
		std::vector<size_t> linkVids = shortest_paths[vid0];
		std::vector<size_t> linkVids1 = shortest_paths[vid1];
		std::reverse(linkVids1.begin(), linkVids1.end());
		std::copy(linkVids1.begin(), linkVids1.end(), std::back_inserter(linkVids));
		loops.push_back(linkVids);
	}
	return loops;
}

Mesh GetBranchCover(const Mesh& mesh, const DualMesh& dualMesh) {
	std::vector<std::vector<Vertex>> V(4, std::vector<Vertex>(4 * mesh.F.size())); // 4 covers, each cover has C faces, each cell has 4 vertices;
	std::vector<std::vector<Cell>>   F(4, std::vector<Cell>(mesh.F.size()));

	std::unordered_map<size_t, size_t> g_coverVid_origVid;
	std::unordered_map<size_t, std::set<size_t>> g_coverVid_coverFids;

	initial(mesh, V, F, g_coverVid_origVid, g_coverVid_coverFids);

	std::vector<std::vector<size_t>> loops;
	std::vector<int> holonomys;

	size_t orig = 0;
	adjacency_list_t adjacency_list(dualMesh.V.size());
	for (auto& v : dualMesh.V)
		for (auto nvid : v.N_Vids) {
			adjacency_list[v.id].push_back(neighbor(nvid, 1));
		}
	std::vector<std::vector<size_t>> shortest_paths;
	for (auto& f : mesh.F)
		shortest_paths.push_back(GetShortestPath(adjacency_list, orig, f.id));
	for (auto& e : mesh.E) {
		auto vid0 = e.N_Fids[0];
		auto vid1 = e.N_Fids[1];
		std::vector<size_t> linkVids = shortest_paths[vid0];
		std::vector<size_t> linkVids1 = shortest_paths[vid1];
		std::reverse(linkVids1.begin(), linkVids1.end());
		std::copy(linkVids1.begin(), linkVids1.end(), std::back_inserter(linkVids));
		loops.push_back(linkVids);
		auto holonomy = GetHolonomyGroup(mesh, linkVids);
		holonomys.push_back(holonomy);
		for (auto i = 0; i < 4; ++i) {
			auto copy_0 = connections[holonomy][i][0];
			auto copy_1 = connections[holonomy][i][1];
			auto& cover_f0 = F[copy_0][e.N_Fids[0]];
			auto& cover_f1 = F[copy_1][e.N_Fids[1]];

			auto& orig_v0 = mesh.E.at(e.Vids[0]);
			auto& orig_v1 = mesh.E.at(e.Vids[1]);

			auto cover_f0_share_vid0 = get_cover_vid(g_coverVid_origVid, cover_f0, orig_v0.id);
			auto cover_f0_share_vid1 = get_cover_vid(g_coverVid_origVid, cover_f0, orig_v1.id);
			auto cover_f1_share_vid0 = get_cover_vid(g_coverVid_origVid, cover_f1, orig_v0.id);
			auto cover_f1_share_vid1 = get_cover_vid(g_coverVid_origVid, cover_f1, orig_v1.id);

			if (cover_f0_share_vid0 == MAXID || cover_f0_share_vid1 == MAXID ||
				cover_f1_share_vid0 == MAXID || cover_f1_share_vid1 == MAXID) {
				std::cerr << "Err in GetBranchCover eid = " << e.id << " holonomy = " << holonomy
					<< " ~(" << copy_0 << ", " << copy_1 << ")\n";
			}

			update(mesh, F, cover_f0_share_vid0, cover_f1_share_vid0, g_coverVid_coverFids);
			update(mesh, F, cover_f0_share_vid1, cover_f1_share_vid1, g_coverVid_coverFids);
		}
		auto base = mesh.E.size() / 10;
		if (e.id % base == 0) std::cout << "Done " << float(e.id) * 100 / mesh.E.size() << "%\n";
	}
	{
		MeshFileWriter writer(mesh, "FacePaths.vtk");
		writer.WriteFacesVtk(loops);
		std::vector<int> newholonomys;
		for (auto i = 0; i < loops.size(); ++i) {
			auto& loop = loops[i];
			for (auto fid : loop) newholonomys.push_back(holonomys[i]);
		}
		writer.WriteCellData(newholonomys, "holonomy", false);
	}
	{
		MeshFileWriter writer(dualMesh, "DualLoops.vtk");
		writer.WriteLinksVtk(loops);
		writer.WriteCellData(holonomys, "holonomy", false);
	}
	{
		MeshFileWriter writer(mesh, "EdgeHolonomy.vtk");
		writer.WriteEdgesVtk();
		writer.WriteCellData(holonomys, "holonomy", true);
	}
	std::cout << "Building CoverMesh V\n";
	std::vector<Vertex> coverV(4 * mesh.F.size() * 4);
	std::vector<Cell>   coverC(4 * mesh.F.size());
	for (auto i = 0; i < 4; ++i) {
		for (size_t j = 0; j < V[i].size(); ++j) {
			auto& v = V[i][j];
			auto& coverv = coverV[v.id];
			coverv.id = v.id;
			coverv = v;
		}
	}
	std::cout << "Building CoverMesh C\n";
	for (auto i = 0; i < 4; ++i) {
		for (size_t j = 0; j < F[i].size(); ++j) {
			auto& c = F[i][j];
			auto& coverc = coverC[c.id];
			coverc.id = c.id;
			coverc.Vids = c.Vids;
			coverc.cellType = VTK_QUAD;
		}
	}
	std::cout << "Finish Building CoverMesh\n";
	Mesh coverMesh(coverV, coverC, QUAD);
	return coverMesh;
}

void TestBranchCover(const Mesh& mesh, const DualMesh& dualMesh) {
	std::cout << "TestBranchCover\n";
	std::vector<std::vector<Vertex>> V(4, std::vector<Vertex>(4 * mesh.F.size())); // 4 covers, each cover has C faces, each cell has 4 vertices;
	std::vector<std::vector<Cell>>   F(4, std::vector<Cell>(mesh.F.size()));

	std::unordered_map<size_t, size_t> g_coverVid_origVid;
	std::unordered_map<size_t, std::set<size_t>> g_coverVid_coverFids;

	initial(mesh, V, F, g_coverVid_origVid, g_coverVid_coverFids);

	std::vector<std::vector<size_t>> loops;
	std::vector<int> holonomys;

	size_t orig = 0;
	adjacency_list_t adjacency_list(dualMesh.V.size());
	for (auto& v : dualMesh.V)
		for (auto nvid : v.N_Vids)
			adjacency_list[v.id].push_back(neighbor(nvid, 1));
	for (auto& e : mesh.E) {
		auto vid0 = e.N_Fids[0];
		auto vid1 = e.N_Fids[1];
		std::vector<size_t> linkVids = GetShortestPath(adjacency_list, orig, vid0);
		std::vector<size_t> linkVids1 = GetShortestPath(adjacency_list, vid1, orig);
		std::copy(linkVids1.begin(), linkVids1.end(), std::back_inserter(linkVids));
		auto holonomy = GetHolonomyGroup(mesh, linkVids);
		if (holonomy != 0) 
			std::cout << "Err holomomy =  " << holonomy << "\n";
		auto base = mesh.E.size() / 10;
		if (e.id % base == 0) std::cout << "Done " << float(e.id) * 100 / mesh.E.size() << "%\n";
	}
}

bool IsManifold(const Mesh& mesh) {
	for (auto& e : mesh.E)
		if (e.N_Fids.size() > 2) return false;
	return true;
}
void test(Mesh& mesh) {
	if (!IsManifold(mesh)) std::cerr << "Not Manifold!\n";
	DualMesh dualMesh;
	std::cout << "dualMesh.Build(mesh);\n";
	//dualMesh.Build_skip_singularities(mesh);
	dualMesh.Build(mesh);
	std::cout << "dualMesh.BuildAllConnectivities();\n";
	dualMesh.BuildAllConnectivities();
	std::cout << "dualMesh.ExtractBoundary();\n";
	dualMesh.ExtractBoundary();
	std::cout << "dualMesh.ExtractSingularities();\n";
	dualMesh.ExtractSingularities();
	std::cout << "dualMesh.BuildParallelE();\n";
	dualMesh.BuildParallelE();
	//std::cout << "dualMesh.FixOrientation(mesh);\n";
	//dualMesh.FixOrientation(mesh);
	//std::cout << "dualMesh.FixOrientation();\n";
	//dualMesh.FixOrientation();
	//std::cout << "TestBranchCover\n";
	TestBranchCover(mesh, dualMesh);
}

std::vector<std::vector<size_t>> GetSingularityShortestPathVidsToCut(const Mesh& mesh, const std::vector<size_t>& cut_graph_eids) {
	auto cut_graph_vids = GetLinkEVids(mesh, cut_graph_eids);
	std::vector<std::vector<size_t>> res;
	adjacency_list_t adjacency_list(mesh.V.size());
	for (auto& v : mesh.V)
		for (auto nvid : v.N_Vids)
			adjacency_list[v.id].push_back(neighbor(nvid, 1));
	for (auto& v : mesh.V) {
		if (!v.isSingularity) continue;
		std::vector<size_t> shortestPath(mesh.V.size());
		for (auto vid : cut_graph_vids) {
			auto linkVids = GetShortestPath(adjacency_list, v.id, vid);
			if (linkVids.size() < shortestPath.size()) shortestPath = linkVids;
		}
		res.push_back(shortestPath);
	}
	return res;
}

std::set<size_t> GetSingularityShortestPathEidsToCut(const std::vector<std::vector<size_t>>& linkVids, 
	const std::unordered_map<size_t, size_t>& key_edgeid) {
	std::set<size_t> res;
	for (auto& vids : linkVids) {
		for (auto i = 1; i < vids.size(); ++i) {
			auto iter = key_edgeid.find((vids[i - 1] << 32) | vids[i]);
			if (iter != key_edgeid.end()) res.insert(iter->second);
		}
	}
	return res;
}

int main(int argc, char* argv[]) {
	if (argc < 3) {
		std::cerr << "Usage: ExtractCoverSpace2 <file> <out.vtk>";
		return -1;
	}
	ArgumentManager argumentManager(argc, argv);
	MeshFileReader reader(argv[1]);
	auto mesh = reader.GetMesh();
	mesh.BuildAllConnectivities();
	mesh.ExtractBoundary();
	mesh.ExtractSingularities();
	mesh.BuildParallelE();
	mesh.FixOrientation();
	mesh.GetNormalOfSurfaceFaces();

	DualMesh dualMesh;
	dualMesh.Build(mesh);
	dualMesh.BuildConnection(mesh);
	dualMesh.ExtractBoundary();
	dualMesh.ExtractSingularities();
	dualMesh.BuildParallelE();
	dualMesh.FixOrientation(mesh);
	//dualMesh.FixOrientation();

	auto loops = get_loops(mesh, dualMesh);
	{
		std::vector<int> newholonomys;
		std::vector<size_t> cellIds;
		//for (auto i = 0; i < loops.size(); ++i) {
		//	auto& loop = loops[i];
		//	for (auto fid : loop) newholonomys.push_back(holonomys[i]);
		//}
		for (auto i = 0; i < loops.size(); ++i) {
			auto& loop = loops[i];
			for (auto cid : loop) cellIds.push_back(cid);
		}
		for (auto i = 0; i < loops.size(); ++i) {
			auto& loop = loops[i];
			for (auto cid : loop) newholonomys.push_back(i);
		}
		MeshFileWriter writer(mesh, "CellPaths.vtu");
		writer.WriteCellsVtu(cellIds, newholonomys, "id");
		//writer.WriteCellData(newholonomys, "id", false);
	}
	//auto coverMesh = GetBranchCover(mesh, dualMesh);
	//coverMesh.RemoveUselessVertices();
	//{
	//	std::cout << "Writing CoverMesh!\n";
	//	MeshFileWriter writer(coverMesh, argv[2]);
	//	writer.WriteFile();
	//	std::vector<int> layers(coverMesh.C.size());
	//	for (auto i = 0; i < layers.size(); ++i)
	//		layers[i] = i / mesh.F.size();
	//	writer.WriteCellData(layers, "layer", true);
	//}
	//coverMesh.BuildAllConnectivities();	
	//coverMesh.ExtractBoundary();
	//coverMesh.ExtractSingularities();
	//coverMesh.BuildParallelE();
	//{
	//	std::vector<size_t> singularVids;
	//	for (auto& v : coverMesh.V)
	//		if (v.isSingularity) singularVids.push_back(v.id);
	//	MeshFileWriter writer(coverMesh, "CoverMeshSingularities.vtk");
	//	writer.WriteVerticesVtk(singularVids);
	//}
	//test(coverMesh);
	return 0;
}
