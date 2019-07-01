/*
* ExtractBaseHomology.cpp
*
*  Created on: Oct 20, 2018
*      Author: cotrik
*/

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "DualMesh.h"
#include "MST.h"
#include "ArgumentManager.h"
#include <iostream>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <stack>
#include <algorithm>
#include <list>
#include <limits> // for numeric_limits
#include <utility> // for pair
#include <iterator>

struct VertexUV : public Vertex {
	size_t father;
	glm::dvec2 uv;
};

void read(const char* filename, std::vector<VertexUV>& V, std::vector<Face>& F, std::vector<Edge>& E) {
	std::ifstream ifs(filename);
	std::string line;
	while (getline(ifs, line)) {
		std::stringstream ss(line);
		std::string keyword;
		VertexUV v;
		ss >> keyword;
		if (keyword == "Vertex") {
			ss >> v.id >> v.x >> v.y >> v.z;
			--v.id;
			V.push_back(v);
		} else if (keyword == "Face") {
			Face f;
			ss >> f.id;
			--f.id;
			size_t vid;
			while (ss >> vid)
				f.Vids.push_back(vid - 1);

			F.push_back(f);
		} else if (keyword == "Edge") {
			Edge e(2);
			ss >> e.Vids[0] >> e.Vids[1];
			--e.Vids[0];
			--e.Vids[1];
			E.push_back(e);
		}
	}
}

void write(const char* filename, const std::vector<VertexUV>& V, const std::vector<Face>& F) {
	std::vector<Vertex> triV;
	Vertex x;
	for (auto& v : V) {
		x.id = v.id;
		x.x = v.x;
		x.y = v.y;
		x.z = v.z;
		triV.push_back(x);
	}
	std::vector<Cell> C;
	Cell c;
	for (auto& f : F) {
		c.id = f.id;
		c.Vids = f.Vids;
		C.push_back(c);
	}
	MeshFileWriter writer(triV, C, filename, F.front().Vids.size() == 3 ? TRIANGLE : QUAD);
	writer.WriteFile();
}

void generate_tri_mesh(const std::vector<VertexUV>& V, const std::vector<Face>& F, Mesh& triMesh) {
	std::vector<Vertex> triV;
	Vertex x;
	for (auto& v : V) {
		x.id = v.id;
		x.x = v.x;
		x.y = v.y;
		x.z = v.z;
		triV.push_back(x);
	}
	std::vector<Cell> C;
	Cell c;
	for (auto& f : F) {
		c.id = f.id;
		c.Vids = f.Vids;
		C.push_back(c);
	}

	triMesh.V = triV;
	triMesh.C = C;
	triMesh.m_cellType = TRIANGLE;
	triMesh.BuildAllConnectivities();
	triMesh.ExtractBoundary();
}

void generate_tri_mesh(const std::vector<Vertex>& V, const std::vector<Face>& F, Mesh& triMesh) {
	std::vector<Vertex> triV;
	Vertex x;
	for (auto& v : V) {
		x.id = v.id;
		x.x = v.x;
		x.y = v.y;
		x.z = v.z;
		triV.push_back(x);
	}
	std::vector<Cell> C;
	Cell c;
	for (auto& f : F) {
		c.id = f.id;
		c.Vids = f.Vids;
		C.push_back(c);
	}

	triMesh.V = triV;
	triMesh.C = C;
	triMesh.m_cellType = TRIANGLE;
	triMesh.BuildAllConnectivities();
	triMesh.ExtractBoundary();
}

void generate_quad_mesh(const std::vector<Vertex>& V, const std::vector<Face>& F, Mesh& quadMesh) {
	std::vector<Cell> C;
	Cell c;
	for (auto& f : F) {
		c.id = f.id;
		c.Vids = f.Vids;
		C.push_back(c);
	}

	quadMesh.V = V;
	quadMesh.C = C;
	quadMesh.m_cellType = QUAD;
	quadMesh.BuildAllConnectivities();
	quadMesh.ExtractBoundary();
}

void write(const char* filename, const std::vector<Vertex>& V, const std::vector<Face>& F) {
	std::vector<Vertex> triV;
	Vertex x;
	for (auto& v : V) {
		x.id = v.id;
		x.x = v.x;
		x.y = v.y;
		x.z = v.z;
		triV.push_back(x);
	}
	std::vector<Cell> C;
	Cell c;
	for (auto& f : F) {
		c.id = f.id;
		c.Vids = f.Vids;
		C.push_back(c);
	}
	MeshFileWriter writer(triV, C, filename, F.front().Vids.size() == 3 ? TRIANGLE : QUAD);
	writer.WriteFile();
}

std::vector<std::vector<size_t>> GetBaseEdgeIds(const Mesh& mesh, const size_t numOfBase) {
	std::unordered_map<size_t, size_t> key_edgeId;
	for (auto& e : mesh.E) {
		key_edgeId[(e.Vids[0] << 32) | e.Vids[1]] = e.id;
		key_edgeId[(e.Vids[1] << 32) | e.Vids[0]] = e.id;
	}

	std::vector<std::vector<size_t>> baseEdgeIds;
	for (auto i = 0; i < numOfBase; ++i) {
		std::vector<VertexUV> V;
		std::vector<Face> F;
		std::vector<Edge> E;
		std::string filename = std::string("base_") + std::to_string(i) + ".m";
		read(filename.c_str(), V, F, E);
		std::vector<size_t> eids;
		for (auto& e : E)
			eids.push_back(key_edgeId[(e.Vids[0] << 32) | e.Vids[1]]);
		baseEdgeIds.push_back(eids);
	}
	return baseEdgeIds;
}

std::vector<std::vector<size_t>> GetBaseLinkVids(const Mesh& mesh, const std::vector<std::vector<size_t>>& baseEdgeIds) {
	std::vector<std::vector<size_t>> baseLinkVids;
	for (auto& eids_ : baseEdgeIds) {
		auto eids = eids_;
		std::vector<size_t> ringEids;
		ringEids.push_back(eids.back());
		eids.pop_back();
		while (ringEids.size() < eids_.size()) {
			auto last_eid = ringEids.back();
			int i = 0;
			for (auto next_eid : eids) {
				auto& last_e = mesh.E.at(last_eid);
				std::set<size_t> vids(last_e.Vids.begin(), last_e.Vids.end());
				auto& next_e = mesh.E.at(next_eid);
				vids.insert(next_e.Vids.begin(), next_e.Vids.end());
				if (vids.size() == 3) {
					ringEids.push_back(next_eid);
					eids.erase(eids.begin() + i);
					break;
				}
				++i;
			}
		}
		std::vector<size_t> ringVids;
		auto& front_e = mesh.E.at(ringEids.front());
		ringVids.push_back(front_e.Vids[0]);
		for (auto eid : ringEids) {
			auto last_vid = ringVids.back();
			auto& e = mesh.E.at(eid);
			auto next_vid = e.Vids[0] == last_vid ? e.Vids[1] : e.Vids[0];
			ringVids.push_back(next_vid);
		}
		baseLinkVids.push_back(ringVids);
	}
	return baseLinkVids;
}

void RemoveVertices(const Mesh& mesh, std::vector<std::vector<size_t>>& baseLinkVids) {
	for (auto& link : baseLinkVids) {
		while (true) {
			bool remove = false;
			for (int i = 1; i < link.size() - 1; ++i) {
				auto& vid_curr = link[i];
				auto& vid_prev = link[i - 1];
				auto& vid_next = link[i + 1];
				bool isAllInOneQuad = false;
				for (auto nfid : mesh.V.at(vid_curr).N_Fids) {
					int count = 0;
					for (auto nfvid : mesh.F.at(nfid).Vids) {
						if (nfvid == vid_prev) ++count;
						else if (nfvid == vid_next) ++count;
					}
					if (count == 2) {
						isAllInOneQuad = true;
						break;
					}
				}
				if (isAllInOneQuad) {
					link.erase(link.begin() + i);
					remove = true;
					break;
				}
			}
			if (!remove) break;
		}
		if (link.front() != link.back()) link.push_back(link.front());
	}
}

void ChangeVertices(const Mesh& mesh, std::vector<std::vector<size_t>>& baseLinkVids) {
	for (auto& link : baseLinkVids) {
		int iters = link.size();
		while (true) {
			bool change = false;
			for (int i = 1; i < link.size() - 1; ++i) {
				auto& vid_curr = link[i];
				auto& vid_prev = link[i - 1];
				auto& vid_next = link[i + 1];
				bool found_prev = false;
				bool found_next = false;
				for (auto nvid : mesh.V.at(vid_curr).N_Vids) {
					if (nvid == vid_prev) found_prev = true;
					else if (nvid == vid_next) found_next = true;
				}
				if (!found_prev && !found_next) {
					for (auto neid : mesh.V.at(vid_curr).N_Eids) {
						auto& e = mesh.E.at(neid);
						std::set<size_t> vids;
						for (auto nfid : e.N_Fids) {
							auto& nf = mesh.F.at(nfid);
							vids.insert(nf.Vids.begin(), nf.Vids.end());
						}
						if (vids.find(vid_prev) != vids.end() && vids.find(vid_next) != vids.end()) {
							auto vid_change = e.Vids[0] == vid_curr ? e.Vids[1] : e.Vids[0];
							if (link.back() == vid_curr) link.back() = vid_change;
							if (link.front() == vid_curr) link.front() = vid_change;
							vid_curr = vid_change;
							change = true;
							break;
						}
					}
					if (change) break;
				}
			}
			if (!change) break;
		}
		if (link.front() != link.back()) link.push_back(link.front());
	}
}

void InsertVertices(const Mesh& mesh, std::vector<std::vector<size_t>>& baseLinkVids) {
	std::unordered_map<size_t, size_t> key_edgeId;
	for (auto& e : mesh.E) {
		key_edgeId[(e.Vids[0] << 32) | e.Vids[1]] = e.id;
		key_edgeId[(e.Vids[1] << 32) | e.Vids[0]] = e.id;
	}
	for (auto& link : baseLinkVids) {
		std::vector<size_t> new_link;
		new_link.push_back(link.front());
		for (int i = 1; i < link.size(); ++i) {
			auto& vid_curr = link[i];
			auto& vid_prev = link[i - 1];
			size_t key = (vid_curr << 32) | vid_prev;
			if (key_edgeId.find(key) == key_edgeId.end()) {
				for (auto nvid : mesh.V.at(vid_curr).N_Vids) {
					bool found = false;
					for (auto nnvid : mesh.V.at(nvid).N_Vids) {
						if (nnvid == vid_prev) {
							found == true;
							new_link.push_back(nvid);
							vid_prev = nvid;
							break;
						}
					}
					if (found) break;
				}
			}
			new_link.push_back(vid_curr);
		}
		link = new_link;
		if (link.front() != link.back()) link.push_back(link.front());
	}
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
					std::vector<weight_t> min_distance;
					std::vector<vertex_t> previous;
					Graph g;
					g.DijkstraComputePaths(src, adjacency_list, min_distance, previous);
					std::list<vertex_t> path = g.DijkstraGetShortestPathTo(dest, previous);
					std::vector<size_t> path_vids(path.begin(), path.end());
					for (int j = 1; j < path_vids.size() - 1; ++j) {
						new_link.push_back(path_vids[j]);
					}
				}
			} else {
				new_link.push_back(v.id);
			}
		}
		new_link.push_back(link.back());
		link = new_link;
	}
}

std::unordered_map<size_t, size_t> key_dualEdgeId;
void get_key_dualEdgeId(const DualMesh& dualMesh) {
	for (auto& e : dualMesh.E) {
		key_dualEdgeId[(e.Vids[0] << 32) | e.Vids[1]] = e.id;
		key_dualEdgeId[(e.Vids[1] << 32) | e.Vids[0]] = e.id;
	}
}

std::set<size_t> get_intersect(const std::set<size_t>& s1, const std::set<size_t>& s2) {
    std::set<size_t> intersect;
	set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), std::inserter(intersect, intersect.begin()));
	return intersect;
}

int GetHolonomyGroup(const DualMesh& dualMesh, const std::vector<size_t>& baseLinkVids) {
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
		auto share_eid = *eids.begin();
		auto& share_e = dualMesh.E.at(share_eid);
		if (share_eid == parallel_eid) {
			std::set<size_t> eids1(share_e.parallelEids.begin(), share_e.parallelEids.end());
			std::set<size_t> eids2(curr_f.Eids.begin(), curr_f.Eids.end());
			std::set<size_t> eids = Util::get_intersect(eids1, eids2);
			if (eids.size() != 1) {
				std::cerr << "Err in GetHolonomyGroup parallel_eids, share_eid = " << share_eid << std::endl;
				continue;
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
				std::cerr << "Err in GetHolonomyGroup consecutive_eids, share_eid = " << share_eid << std::endl;
				continue;
			}
		}
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

std::set<int> GetHolonomyGroup(const DualMesh& dualMesh, const std::vector<std::vector<size_t>>& baseLinkVids) {
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


void Clean(std::vector<std::vector<size_t>>& baseLinkVids, const char* filename) {
	MeshFileReader reader(filename);
	auto mesh = reader.GetMesh();
	mesh.BuildAllConnectivities();
	mesh.ExtractBoundary();
	mesh.ExtractSingularities();

	DualMesh dualMesh;
	dualMesh.Build(mesh);
	dualMesh.BuildAllConnectivities();	
	dualMesh.BuildParallelE();
	//dualMesh.BuildConsecutiveE();
	//dualMesh.BuildOrthogonalE();
	//{
	//	MeshFileWriter writer(dualMesh, "FacePath0.vtk");
	//	writer.WriteFacesVtk(baseLinkVids);
	//}
	AddVertices(mesh, baseLinkVids);
	//{
	//	MeshFileWriter writer(dualMesh, "AddVertices.vtk");
	//	writer.WriteFacesVtk(baseLinkVids);
	//}
	RemoveVertices(mesh, baseLinkVids); 
	//{
	//	MeshFileWriter writer(dualMesh, "RemoveVerticesFacePath.vtk");
	//	writer.WriteFacesVtk(baseLinkVids);
	//}
	ChangeVertices(mesh, baseLinkVids); 
	//{
	//	MeshFileWriter writer(dualMesh, "ChangeVerticesFacePath.vtk");
	//	writer.WriteFacesVtk(baseLinkVids);
	//}
	InsertVertices(mesh, baseLinkVids); 
	//{
	//	MeshFileWriter writer(dualMesh, "InsertVerticesFacePath.vtk");
	//	writer.WriteFacesVtk(baseLinkVids);
	//}
	RemoveSingularities(mesh, baseLinkVids);
	AddVertices(mesh, baseLinkVids);
	MeshFileWriter writer(dualMesh, "FacePath.vtk");
	writer.WriteFacesVtk(baseLinkVids);

	auto holonomyGroup = GetHolonomyGroup(dualMesh, baseLinkVids);
	std::cout << "HolonomyGroup:";
	for (auto holonomy : holonomyGroup)
		std::cout << " " << holonomy;
	std::cout << std::endl;
}

int main(int argc, char* argv[]) {
	if (argc < 4) {
		std::cout << "Usage: ExtractBaseHomologyQuad input.vtk output.vtk num_of_base quad=<quad.vtk>\n";
		return -1;
	}
	ArgumentManager am(argc, argv);
	std::string strQuadMeshFileName = am.get("quad");
	MeshFileReader reader(argv[1]);
	auto mesh = reader.GetMesh();
	mesh.BuildAllConnectivities();

	size_t numOfBase = std::stoi(argv[3]);
	
	std::vector<std::vector<size_t>> baseEdgeIds = GetBaseEdgeIds(mesh, numOfBase);
	std::vector<std::vector<size_t>> baseLinkVids = GetBaseLinkVids(mesh, baseEdgeIds);
	if (!strQuadMeshFileName.empty())
		Clean(baseLinkVids, strQuadMeshFileName.c_str());
	MeshFileWriter writer(mesh, argv[2]);
	writer.WriteLinksVtk(baseLinkVids);
	return 0;
}
