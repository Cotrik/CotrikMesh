/*
 * ExtractHomotopyLoops.cpp
 *
 *  Created on: Jan 7, 2019
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

std::string get_facekey(const std::vector<size_t>& Vids) {
	std::set<size_t> vids(Vids.begin(), Vids.end());
	std::string s;
	for (auto vid : vids)
		s += std::to_string(vid) + "@";
	return s;
}

std::string get_facekey(const Face& f) {
	return get_facekey(f.Vids);
}

std::unordered_map<std::string, size_t> get_key_faceId(const Mesh& mesh) {
	std::unordered_map<std::string, size_t> key_faceId;
	for (size_t i = 0; i < mesh.F.size(); ++i) {
		const auto& f = mesh.F.at(i);
		std::string s = get_facekey(f);
		key_faceId[s] = i;
	}
	return key_faceId;
}

std::map<size_t, size_t> get_dualEid_meshFid(const Mesh& mesh, const Mesh& dualMesh) {
	std::map<size_t, size_t> dualEid_meshFid;
	auto key_faceId = get_key_faceId(mesh);
	for (auto& e : dualMesh.E) {
		if (e.isBoundary) continue;
		auto key = get_facekey(e.N_Cids);
		auto iter = key_faceId.find(key);
		if (iter == key_faceId.end()) {
			std::cerr << "Err in get_dualEid_meshFid, eid = " << e.id << std::endl;
			continue;
		}
		dualEid_meshFid[e.id] = iter->second;
	}
	return dualEid_meshFid;
}

std::vector<size_t> get_cut_graph_fids(const Mesh& mesh, const Mesh& dualMesh, const std::vector<size_t>& reversed_mst_eids) {
	auto dualEid_meshFid = get_dualEid_meshFid(mesh, dualMesh);
	std::vector<size_t> cut_graph_fids;
	for (auto eid : reversed_mst_eids) {
		auto& e = dualMesh.E.at(eid);
		if (e.isBoundary) continue;
		auto iter = dualEid_meshFid.find(e.id);
		if (iter == dualEid_meshFid.end())	std::cerr << "get_cut_graph_fids Error\n";
		cut_graph_fids.push_back(iter->second);
	}
	return cut_graph_fids;
}

bool check(const Mesh& mesh, const std::unordered_set<size_t>& vids, const std::unordered_set<size_t>& eids,
	const std::unordered_set<size_t>& fids, const std::unordered_set<size_t>& cids) {
	static std::vector<size_t> v_visited(mesh.V.size(), false);
	static std::vector<size_t> e_visited(mesh.E.size(), false);
	for (auto vid : vids) {
		if (v_visited[vid]) continue;
		auto& v = mesh.V.at(vid);
		std::vector<size_t> ncids;
		for (auto ncid : v.N_Cids)
			if (cids.find(ncid) != cids.end()) {
				ncids.push_back(ncid);
				if (ncids.size() > 2) {
					v_visited[vid] = true;
					break;
				}
			}
		if (ncids.size() == 2) {
			auto& c0 = mesh.C.at(ncids[0]);
			auto& c1 = mesh.C.at(ncids[1]);
			std::set<size_t> s(c0.Vids.begin(), c0.Vids.end());
			s.insert(c1.Vids.begin(), c1.Vids.end());
			if (s.size() > 12) return false;
		}
	}
	for (auto eid : eids) {
		if (e_visited[eid]) continue;
		auto& e = mesh.E.at(eid);
		std::vector<size_t> ncids;
		for (auto ncid : e.N_Cids)
			if (cids.find(ncid) != cids.end()) {
				ncids.push_back(ncid);
				if (ncids.size() > 2) {
					e_visited[eid] = true;
					break;
				}
			}
		if (ncids.size() == 2) {
			auto& c0 = mesh.C.at(ncids[0]);
			auto& c1 = mesh.C.at(ncids[1]);
			std::set<size_t> s(c0.Vids.begin(), c0.Vids.end());
			s.insert(c1.Vids.begin(), c1.Vids.end());
			if (s.size() > 12) return false;
		}
	}
	return true;
}

int get_genus(const std::unordered_set<size_t>& vids, const std::unordered_set<size_t>& eids,
	const std::unordered_set<size_t>& fids, const std::unordered_set<size_t>& cids) {
	auto g = 1 - ((long long)vids.size() - eids.size() + fids.size() - cids.size());
	return g;
}

bool check(const Mesh& mesh, const std::set<size_t>& cids, const std::vector<bool> vids, const std::vector<bool> eids) {
	for (auto& v : mesh.V) {
		if (!vids[v.id]) continue;
		std::vector<size_t> ncids;
		for (auto ncid : v.N_Cids)
			if (cids.find(ncid) != cids.end()) ncids.push_back(ncid);
		if (ncids.size() == 2) {
			auto& c0 = mesh.C.at(ncids[0]);
			auto& c1 = mesh.C.at(ncids[1]);
			std::set<size_t> s(c0.Vids.begin(), c0.Vids.end());
			s.insert(c1.Vids.begin(), c1.Vids.end());
			if (s.size() > 12) return false;
		}
	}

	for (auto& e : mesh.E) {
		if (!eids[e.id]) continue;
		std::vector<size_t> ncids;
		for (auto ncid : e.N_Cids)
			if (cids.find(ncid) != cids.end()) ncids.push_back(ncid);
		if (ncids.size() == 2) {
			auto& c0 = mesh.C.at(ncids[0]);
			auto& c1 = mesh.C.at(ncids[1]);
			std::set<size_t> s(c0.Vids.begin(), c0.Vids.end());
			s.insert(c1.Vids.begin(), c1.Vids.end());
			if (s.size() > 12) return false;
		}
	}
	return true;
}

int get_genus(const Mesh& mesh, const std::vector<size_t>& Cids) {
	std::set<size_t> cids(Cids.begin(), Cids.end());
	std::vector<bool> vids(mesh.V.size(), false);
	std::vector<bool> eids(mesh.E.size(), false);
	std::vector<bool> fids(mesh.F.size(), false);
	for (auto cid : cids) {
		auto& c = mesh.C.at(cid);
		for (auto id : c.Vids) vids[id] = true;
		for (auto id : c.Eids) eids[id] = true;
		for (auto id : c.Fids) fids[id] = true;
		//for (auto id : Cids) cids[id] = true;
	}



	long long V = 0, E = 0, F = 0, C = cids.size();
	for (auto id : vids) if (id) ++V;
	for (auto id : eids) if (id) ++E;
	for (auto id : fids) 
		if (id) 
			++F;

	auto g = 1 - (V - E + F - C);
	if (g == 0 && !check(mesh, cids, vids, eids)) return -1;
	return g;
}

void RegionGrowCell(const Mesh& mesh, const Cell& start_cell, std::vector<bool>& is_cell_visited, std::vector<size_t>& cids) {
	is_cell_visited[start_cell.id] = true;
	cids.push_back(start_cell.id);
	if (cids.size() % 100 == 0) std::cout << "RegionGrowCell cids.size() = " << cids.size() << std::endl;
	for (auto fid : start_cell.Fids) {
		auto& f = mesh.F.at(fid);
		for (auto ncid : f.N_Cids) {
			if (is_cell_visited[ncid]) continue;
			// cids.push_back(ncid);
			if (get_genus(mesh, cids) != 0) cids.pop_back();
			else RegionGrowCell(mesh, mesh.C.at(ncid), is_cell_visited, cids);
		}
	}
}

void RegionGrowCell(const Mesh& mesh) {
	std::vector<bool> is_cell_visited(mesh.C.size(), false);
	std::vector<size_t> cids;
	RegionGrowCell(mesh, mesh.C.at(0), is_cell_visited, cids);
	{
		std::set<size_t> Cids(cids.begin(), cids.end());
		cids.clear();
		std::copy(Cids.begin(), Cids.end(), std::back_inserter(cids));
		{
			std::cout << "Writing cut_sphere.vtk\n";
			MeshFileWriter writer(mesh, "cut_sphere.vtk");
			writer.WriteCellsVtk(cids);
		}
		std::vector<size_t> cut_cids;
		for (auto & c : mesh.C)
			if (Cids.find(c.id) == Cids.end()) cut_cids.push_back(c.id);
		{
			std::cout << "Writing cut_cells.vtk\n";
			MeshFileWriter writer(mesh, "cut_cells.vtk");
			writer.WriteCellsVtk(cut_cids);
		}
	}
}
std::set<size_t> bfs(const Mesh& mesh, std::unordered_set<size_t>& vids, std::unordered_set<size_t>& eids,
	std::unordered_set<size_t>& fids, std::unordered_set<size_t>& cids) {
	{
		auto& c0 = mesh.C.at(0);
		vids.insert(c0.Vids.begin(), c0.Vids.end());
		eids.insert(c0.Eids.begin(), c0.Eids.end());
		fids.insert(c0.Fids.begin(), c0.Fids.end());
		cids.insert(0);
	}
	std::set<size_t> res;
	//std::set<size_t> cut_cids;
	while (true) {
		auto n = cids.size();
		for (auto vid : vids) {
			auto& v = mesh.V.at(vid);
			for (auto ncid : v.N_Cids) {
				if (cids.find(ncid) != cids.end()) continue;
				auto& nc = mesh.C.at(ncid);
				std::vector<size_t> insert_vids, insert_eids, insert_fids;
				for (auto id : nc.Vids)
					if (vids.find(id) == vids.end()) {
						insert_vids.push_back(id);
						vids.insert(id);
					}

				for (auto id : nc.Eids)
					if (eids.find(id) == eids.end()) {
						insert_eids.push_back(id);
						eids.insert(id);
					}

				for (auto id : nc.Fids)
					if (fids.find(id) == fids.end()) {
						insert_fids.push_back(id);
						fids.insert(id);
					}

				cids.insert(ncid);
				auto& V = vids;
				auto& E = eids;
				auto& F = fids;
				auto& C = cids;
				//auto V = vids;
				//auto E = eids;
				//auto F = fids;
				//auto C = cids;
				//V.insert(nc.Vids.begin(), nc.Vids.end());
				//E.insert(nc.Eids.begin(), nc.Eids.end());
				//F.insert(nc.Fids.begin(), nc.Fids.end());
				//C.insert(nc.id);
				auto g = get_genus(V, E, F, C);
				if (g == 0 && check(mesh, V, E, F, C)) {
					//auto& c0 = nc;
					//vids.insert(c0.Vids.begin(), c0.Vids.end());
					//eids.insert(c0.Eids.begin(), c0.Eids.end());
					//fids.insert(c0.Fids.begin(), c0.Fids.end());
					//cids.insert(c0.id);
				} else {
					//cut_cids.insert(ncid);
					for (auto id : insert_vids) vids.erase(id);
					for (auto id : insert_eids) eids.erase(id);
					for (auto id : insert_fids) fids.erase(id);
					cids.erase(ncid);
				}
			}
		}
		if (cids.size() == n) break;
		std::cout << "Done " << double(cids.size()) * 100 / mesh.C.size() << "%" << std::endl;
	}
	//std::vector<size_t> C;
	//for (auto cid : cids)
	//	C.push_back(cid);
	//{
	//	std::cout << "Writing cut_cells.vtu\n";
	//	MeshFileWriter writer(mesh, "cut_cells.vtu");
	//	writer.WriteCellsVtu(C);
	//}
	std::set<size_t> Cids(cids.begin(), cids.end());
	{
		std::cout << "Writing cut_sphere.vtk\n";
		MeshFileWriter writer(mesh, "cut_sphere.vtk");
		writer.WriteCellsVtk(Cids);
	}
	std::vector<size_t> cut_cids;
	for (auto & c : mesh.C)
		if (Cids.find(c.id) == Cids.end()) cut_cids.push_back(c.id);
	{
		std::cout << "Writing cut_cells.vtk\n";
		MeshFileWriter writer(mesh, "cut_cells.vtk");
		writer.WriteCellsVtk(cut_cids);
	}
	return res;
}
std::set<size_t> dfs_without_hole(const Mesh& mesh, std::unordered_set<size_t>& vids, std::unordered_set<size_t>& eids,
	std::unordered_set<size_t>& fids, std::unordered_set<size_t>& cids) {
	std::set<size_t> res;
	std::set<size_t> cut_cids;
	std::stack<size_t> st; // stack of cells;
	st.push(0);
	while (!st.empty()) {
		auto t = st.top();
		st.pop();
		auto& c = mesh.C.at(t);
		{
			vids.insert(c.Vids.begin(), c.Vids.end());
			eids.insert(c.Eids.begin(), c.Eids.end());
			fids.insert(c.Fids.begin(), c.Fids.end());
			cids.insert(c.id);
		}
		for (auto fid : c.Fids) {
			// if (fids.find(fid) != fids.end()) continue;
			auto& f = mesh.F.at(fid);
			for (auto ncid : f.N_Cids) {
				if (cids.find(ncid) != cids.end()) continue;
				auto& nc = mesh.C.at(ncid);
				std::vector<size_t> insert_vids, insert_eids, insert_fids;
				for (auto id : nc.Vids) 
					if (vids.find(id) == vids.end()) {
						insert_vids.push_back(id);
						vids.insert(id);
					}

				for (auto id : nc.Eids)
					if (eids.find(id) == eids.end()) {
						insert_eids.push_back(id);
						eids.insert(id);
					}

				for (auto id : nc.Fids)
					if (fids.find(id) == fids.end()) {
						insert_fids.push_back(id);
						fids.insert(id);
					}

				cids.insert(ncid);
				auto& V = vids;
				auto& E = eids;
				auto& F = fids;
				auto& C = cids;
				//V.insert(nc.Vids.begin(), nc.Vids.end());
				//E.insert(nc.Eids.begin(), nc.Eids.end());
				//F.insert(nc.Fids.begin(), nc.Fids.end());
				//C.insert(nc.id);
				auto g = get_genus(V, E, F, C);
				if (g == 0 && check(mesh, V, E, F, C)) {
					st.push(ncid);
				} else {
					cut_cids.insert(ncid);
					res.insert(fid);
				}
			}
		}
	}
	std::vector<size_t> C;
	for (auto cid : cids)
		C.push_back(cid);
	{
		std::cout << "Writing cut_cells.vtu\n";
		MeshFileWriter writer(mesh, "cut_cells.vtu");
		writer.WriteCellsVtu(C);
	}
	return res;
}

std::set<size_t> bfs_without_hole(const Mesh& mesh, std::unordered_set<size_t>& vids, std::unordered_set<size_t>& eids,
	std::unordered_set<size_t>& fids, std::unordered_set<size_t>& cids) {
	std::set<size_t> res;
	std::set<size_t> cut_cids;
	std::queue<size_t> q; // stack of cells;
	q.push(0);
	while (!q.empty()) {
		auto n = q.size();
		//for (auto i = 0; i < n; ++i) 
		{
			auto t = q.front();
			q.pop();
			auto& c = mesh.C.at(t);
			{
				vids.insert(c.Vids.begin(), c.Vids.end());
				eids.insert(c.Eids.begin(), c.Eids.end());
				fids.insert(c.Fids.begin(), c.Fids.end());
				cids.insert(c.id);
			}
			for (auto fid : c.Fids) {
				// if (fids.find(fid) != fids.end()) continue;
				auto& f = mesh.F.at(fid);
				for (auto ncid : f.N_Cids) {
					if (cids.find(ncid) != cids.end()) continue;
					auto& nc = mesh.C.at(ncid);
					auto V = vids;
					auto E = eids;
					auto F = fids;
					auto C = cids;
					V.insert(nc.Vids.begin(), nc.Vids.end());
					E.insert(nc.Eids.begin(), nc.Eids.end());
					F.insert(nc.Fids.begin(), nc.Fids.end());
					C.insert(nc.id);
					auto g = get_genus(V, E, F, C);
					if (g == 0) q.push(ncid);
					else {
						cut_cids.insert(ncid);
						res.insert(fid);
					}
				}
			}
		}
	}
	std::vector<size_t> C;
	for (auto cid : cids)
		C.push_back(cid);
	{
		std::cout << "Writing cut_cells.vtu\n";
		MeshFileWriter writer(mesh, "cut_cells.vtu");
		writer.WriteCellsVtu(C);
	}
	return res;
}

std::set<size_t> get_cut_faces(const Mesh& mesh) {
	std::set<size_t> res;
	std::unordered_set<size_t> vids, eids, fids, cids;
	return bfs(mesh, vids, eids, fids, cids);
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: ExtractHomotopyBases <file> <out.vtk>";
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

	if (mesh.m_cellType == POLYGON || mesh.m_cellType == TRIANGLE || mesh.m_cellType == QUAD) {
		DualMesh dualMesh;
		dualMesh.Build(mesh);
		if (dualMesh.m_cellType == POLYGON) dualMesh.BuildAllConnectivities();
		else if (dualMesh.m_cellType == POLYHEDRA) dualMesh.BuildConnection(mesh);
		dualMesh.ExtractBoundary();
		dualMesh.ExtractSingularities();
		dualMesh.BuildParallelE();
		{
			std::cout << "Writing dual_edges.vtk\n";
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
			std::cout << "Writing dual_mst.vtk\n";
			MeshFileWriter writer(dualMesh, "dual_mst.vtk");
			writer.WriteEdgesVtk(mst_eids);
		}
		{
			std::cout << "Writing dual_rmst.vtk\n";
			MeshFileWriter writer(dualMesh, "dual_rmst.vtk");
			writer.WriteEdgesVtk(reversed_mst_eids);
		}
		auto cut_graph_eids = get_cut_graph_eids(mesh, dualMesh, reversed_mst_eids);
		auto cut_graph_eids_without_prune = cut_graph_eids;
		{
			std::cout << "Writing CutGraph0.vtk\n";
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
				std::cout << "Writing baseLinkVids.vtk\n";
				MeshFileWriter writer(mesh, "baseLinkVids.vtk");
				writer.WriteLinksVtk(linkVids);
			}

			RemoveSingularities(mesh, linkVids);
			// AddVertices(mesh, linkVids);
			g.PruneLinkVids(mesh, linkVids);
			{
				std::cout << "Writing cleanBaseLinkVids.vtk\n";
				MeshFileWriter writer(mesh, "cleanBaseLinkVids.vtk");
				writer.WriteLinksVtk(linkVids);
			}
			auto holonomyGroup = GetHolonomyGroup(dualMesh, linkVids);
			std::cout << "HolonomyGroup:";
			for (auto holonomy : holonomyGroup)
				std::cout << " " << holonomy;
			std::cout << std::endl;
			{
				std::cout << "Writing face_path.vtk\n";
				MeshFileWriter writer(dualMesh, "face_path.vtk");
				writer.WriteFacesVtk(linkVids);
			}
		}
	} else if (mesh.m_cellType == HEXAHEDRA || mesh.m_cellType == TETRAHEDRA || mesh.m_cellType == POLYHEDRA) {
		/*std::map<size_t, size_t> dualEid_meshFid;
		auto cut_graph_fids = get_cut_graph_fids(mesh, dualMesh, reversed_mst_eids);
		{
			std::cout << "Writing CutGraphF0.vtk\n";
			MeshFileWriter writer(mesh, "CutGraphF0.vtk");
			writer.WriteFacesVtk(cut_graph_fids);
		}*/
		auto cut_graph_fids = get_cut_faces(mesh);
		//Graph g;
		//g.prune(mesh, cut_graph_fids);
		{
			std::cout << "Writing CutGraphF0.vtk\n";
			MeshFileWriter writer(mesh, "CutGraphF0.vtk");
			writer.WriteFacesVtk(cut_graph_fids);
		}
		// RegionGrowCell(mesh);
	}

    return 0;
}
