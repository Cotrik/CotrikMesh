/*
* ExtractLoops.cpp
*
*  Created on: Nove 27, 2018
*      Author: cotrik
*/

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "ArgumentManager.h"
#include "DualMesh.h"
#include "MST.h"
#include <iostream>
#include <fstream>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <stack>
#include <algorithm>
#include <list> 
#include <limits.h> 


std::vector<std::vector<size_t>> extractCycleVids(const Mesh& mesh, const std::vector<size_t>& eids) {
	CycleExtractor cycleExtractor(mesh, eids);
	return cycleExtractor.cycleVids;
}

////////////////////////////////////////////////////////////////
int main(int argc, char* argv[]) {
	if (argc < 3) {
		std::cerr << "Usage: ExtractDualMesh <file> <loops.vtk>";
		return -1;
	}

	ArgumentManager argumentManager(argc, argv);
	MeshFileReader reader(argv[1]);
	auto mesh = reader.GetMesh();
	mesh.BuildAllConnectivities();
	mesh.ExtractBoundary();
	mesh.ExtractSingularities();
	//mesh.BuildParallelE();
	//mesh.BuildConsecutiveE();
	//mesh.BuildOrthogonalE();

	Graph g(mesh.V.size(), mesh.E.size());
	for (auto& e : mesh.E) g.addEdge(e.Vids[0], e.Vids[1], 1);
	auto mst_wt = g.kruskalMST();

	std::unordered_map<size_t, size_t> key_edgeid;
	for (auto& e : mesh.E) {
		key_edgeid[(e.Vids[0] << 32) | e.Vids[1]] = e.id;
		key_edgeid[(e.Vids[1] << 32) | e.Vids[0]] = e.id;
	}

	std::vector<size_t> mst_eids;
	for (auto& p : g.mst_edges) mst_eids.push_back(key_edgeid[(p.first << 32) + p.second]);
	std::set<size_t> mst_eids_set(mst_eids.begin(), mst_eids.end());
	std::vector<size_t> reversed_mst_eids;
	for (auto& e : mesh.E)
		if (mst_eids_set.find(e.id) == mst_eids_set.end()) reversed_mst_eids.push_back(e.id);

	std::vector<std::vector<size_t>> total_cycleVids;
	std::vector<std::set<size_t>> total_cycleEids;
	std::unordered_map<size_t, size_t> key_edgeIds;
	for (auto edgeid : mst_eids) {
		auto& e = mesh.E.at(edgeid);
		key_edgeIds[(e.Vids[0] << 32) | e.Vids[1]] = e.id;
		key_edgeIds[(e.Vids[1] << 32) | e.Vids[0]] = e.id;
	}
	for (auto eid : reversed_mst_eids) {
		auto eids = mst_eids;
		eids.push_back(eid);
		std::cout << "process eid = " << eid << " n = " << mesh.E.size() << std::endl;
		auto cycleVids = extractCycleVids(mesh, eids);
		if (cycleVids.size() != 1) {
			std::cerr << "Err!\n";
			continue;
		}
		if (cycleVids.size() == 1) {
			total_cycleVids.push_back(cycleVids.front());
			auto& e = mesh.E.at(eid);
			key_edgeIds[(e.Vids[0] << 32) | e.Vids[1]] = e.id;
			key_edgeIds[(e.Vids[1] << 32) | e.Vids[0]] = e.id;
			std::set<size_t> vids(cycleVids.front().begin(), cycleVids.front().end());
			std::set<size_t> eids;
			for (auto vid : cycleVids.front()) {
				for (auto nvid : mesh.V.at(vid).N_Vids) {
					if (vids.find(nvid) == vids.end()) continue;
					auto key = (vid << 32) | nvid;
					auto iter = key_edgeIds.find(key);
					if (iter != key_edgeIds.end()) eids.insert(iter->second);
				}
			}
			total_cycleEids.push_back(eids);
			key_edgeIds.erase((e.Vids[0] << 32) | e.Vids[1]);
			key_edgeIds.erase((e.Vids[1] << 32) | e.Vids[0]);
		}
	}


	//for (auto& cycleVids : total_cycleVids) {
	//	std::set<size_t> vids(cycleVids.begin(), cycleVids.end());
	//	std::set<size_t> eids;
	//	for (auto vid : cycleVids) {
	//		for (auto nvid : mesh.V.at(vid).N_Vids) {
	//			if (vids.find(nvid) == vids.end()) continue;
	//			auto key = (vid << 32) | nvid;
	//			auto iter = key_edgeIds.find(key);
	//			if (iter != key_edgeIds.end()) eids.insert(iter->second);
	//		}
	//	}
	//	total_cycleEids.push_back(eids);
	//}

	std::vector<std::set<size_t>> long_cycleEids;
	for (auto& eids : total_cycleEids)
		if (eids.size() >= 10) long_cycleEids.push_back(eids);
	total_cycleEids = long_cycleEids;

	std::sort(total_cycleEids.begin(), total_cycleEids.end(), [&](const std::set<size_t>& a, const std::set<size_t>& b) {
		return a.size() > b.size();
	});

	MeshFileWriter writer(mesh, argv[2]);
	writer.WriteEdgesVtk(total_cycleEids);
	return 0;
}
