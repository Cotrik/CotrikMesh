/*
 * ExtractCutGraph.cpp
 *
 *  Created on: Dec 8, 2017
 *      Author: cotrik
 */

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "DualMesh.h"
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
#include "MST.h"

std::unordered_map<size_t, size_t> get_key_edgeid(const Mesh& mesh) {
	std::unordered_map<size_t, size_t> key_edgeid;
	for (auto& e : mesh.E) {
		key_edgeid[(e.Vids[0] << 32) | e.Vids[1]] = e.id;
		key_edgeid[(e.Vids[1] << 32) | e.Vids[0]] = e.id;
	}
	return key_edgeid;
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

void prune(const Mesh& mesh, std::vector<size_t>& cut_graph_eids) {
	std::set<size_t> eids(cut_graph_eids.begin(), cut_graph_eids.end());
	while (true) {
		std::map<size_t, std::set<size_t>> vid_eids;
		for (auto eid : eids) {
			auto& e = mesh.E.at(eid);
			vid_eids[e.Vids[0]].insert(eid);
			vid_eids[e.Vids[1]].insert(eid);
		}
		bool needPrune = false;
		for (auto& item : vid_eids)
			if (item.second.size() == 1) {
				needPrune = true;
				auto eid = *item.second.begin();
				auto iter = eids.find(eid);
				if (iter != eids.end()) eids.erase(eid);
			}
		if (!needPrune) break;
	}

	cut_graph_eids.clear();
	for (auto eid : eids)
		cut_graph_eids.push_back(eid);
}

std::set<size_t> get_cut_graph_intersect_vids(const Mesh& mesh, const std::vector<size_t>& cut_graph_eids) {
	std::set<size_t> cut_graph_eids_set(cut_graph_eids.begin(), cut_graph_eids.end());
	std::set<size_t> cut_graph_intersect_vids;
	for (auto& v : mesh.V) {
		int count = 0;
		for (auto neid : v.N_Eids)
			if (cut_graph_eids_set.find(neid) != cut_graph_eids_set.end()) ++count;
		if (count > 2) cut_graph_intersect_vids.insert(v.id);
	}
	return cut_graph_intersect_vids;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: ExtractCutGraph <file> <out.vtk>";
        return -1;
    }
    ArgumentManager argumentManager(argc, argv);
	MeshFileReader reader(argv[1]);
	auto mesh = reader.GetMesh();
	mesh.BuildAllConnectivities();
	mesh.ExtractBoundary();
	mesh.ExtractSingularities();

	DualMesh dualMesh;
	dualMesh.Build(mesh);
	//std::vector<size_t> quadCids;
	//for (auto& c : dualMesh.C)
	//	if (c.Vids.size() == 4) quadCids.push_back(c.id);
	//std::vector<Cell> newC(quadCids.size());
	//for (auto i = 0; i < newC.size(); ++i) {
	//	newC[i].id = i;
	//	newC[i].Vids = dualMesh.C.at(quadCids[i]).Vids;
	//	newC[i].cellType = VTK_QUAD;
	//}
	//dualMesh.C = newC;
	//{
	//	MeshFileWriter writer(dualMesh, "holes.vtk");
	//	writer.WriteFile();
	//}
	dualMesh.BuildAllConnectivities();
	dualMesh.ExtractBoundary();
	dualMesh.ExtractSingularities();
	//dualMesh.BuildParallelE();
	//dualMesh.FixOrientation(mesh);


	auto key_edgeIds = get_key_edgeid(mesh);
	auto key_dualEdgeIds = get_key_edgeid(dualMesh);
	std::vector<size_t> mst_eids;
	std::vector<size_t> reversed_mst_eids;
	GetMST(dualMesh, key_dualEdgeIds, mst_eids, reversed_mst_eids);

	auto cut_graph_eids = get_cut_graph_eids(mesh, dualMesh, reversed_mst_eids);
	auto cut_graph_eids_without_prune = cut_graph_eids;

	Graph g;
	g.prune(mesh, cut_graph_eids);
	{
		std::cout << "Writing " << argv[2] << "\n";
		MeshFileWriter writer(mesh, argv[2]);
		writer.WriteEdgesVtk(cut_graph_eids);
	}
	auto cut_graph_intersect_vids = get_cut_graph_intersect_vids(mesh, cut_graph_eids);
	{
		std::cout << "Writing cut_graph_intersect_vids\n";
		MeshFileWriter writer(mesh, "cut_graph_intersect_vids.vtk");
		writer.WriteVerticesVtk(cut_graph_intersect_vids);
	}
    return 0;
}
