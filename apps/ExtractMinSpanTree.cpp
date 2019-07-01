/*
 * ExtractMinSpanTree.cpp
 *
 *  Created on: Nove 27, 2017
 *      Author: cotrik
 */

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "ArgumentManager.h"
#include "BaseComplex.h"
#include "BaseComplexQuad.h"
#include "BaseComplexSheetQuad.h"
#include "RefinedDual.h"
#include <iostream>
#include <fstream>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <stack>
#include <algorithm>
using namespace std;

// C++ program for Kruskal's algorithm to find Minimum Spanning Tree of a given connected, undirected and weighted graph 

// Creating shortcut for an integer pair 
typedef pair<size_t, size_t> iPair;
// To represent Disjoint Sets 
struct DisjointSets {
	size_t *parent, *rnk;
	size_t n;

	// Constructor. 
	DisjointSets(int n) {
		// Allocate memory 
		this->n = n;
		parent = new size_t[n + 1];
		rnk = new size_t[n + 1];

		// Initially, all vertices are in different sets and have rank 0. 
		for (size_t i = 0; i <= n; i++) {
			rnk[i] = 0;
			parent[i] = i; // every element is parent of itself 
		}
	}

	// Find the parent of a node 'u' 
	// Path Compression 
	size_t find(size_t u) {
		/* Make the parent of the nodes in the path
		from u--> parent[u] point to parent[u] */
		if (u != parent[u])	parent[u] = find(parent[u]);
		return parent[u];
	}

	// Union by rank 
	void merge(size_t x, size_t y) {
		x = find(x), y = find(y);
		/* Make tree with smaller height a subtree of the other tree */
		if (rnk[x] > rnk[y]) parent[y] = x;
		else parent[x] = y;
		if (rnk[x] == rnk[y]) rnk[y]++;
	}
};
// Structure to represent a graph 
struct Graph {
	size_t V, E;
	vector< pair<size_t, iPair> > edges;
	std::vector<iPair> mst_edges;

	// Constructor 
	Graph(size_t V, size_t E) {
		this->V = V;
		this->E = E;
	}

	// Utility function to add an edge 
	void addEdge(size_t u, size_t v, size_t w) {
		edges.push_back({ w, { u, v } });
	}

	size_t kruskalMST() // Function to find MST using Kruskal's MST algorithm 
	{
		size_t mst_wt = 0;
		sort(edges.begin(), edges.end());

		DisjointSets ds(V);
		for (auto it = edges.begin(); it != edges.end(); ++it) {
			size_t u = it->second.first;
			size_t v = it->second.second;

			size_t set_u = ds.find(u);
			size_t set_v = ds.find(v);

			// Check if the selected edge is creating a cycle or not (Cycle is created if u and v belong to same set) 
			if (set_u != set_v) {
				// Current edge will be in the MST so print it 
				//cout << u << " - " << v << endl;
				mst_edges.push_back(it->second);
				mst_wt += it->first;
				ds.merge(set_u, set_v);
			}
		}

		return mst_wt;
	}
};


int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: ExtractMinSpanTree <file> <out.vtk>";
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

	MeshFileWriter writer(mesh, argv[2]);
	writer.WriteEdgesVtk(mst_eids);

	std::set<size_t> mst_eids_set(mst_eids.begin(), mst_eids.end());
	std::vector<size_t> reversed_mst_eids;
	for (auto& e : mesh.E)
		if (mst_eids_set.find(e.id) == mst_eids_set.end()) reversed_mst_eids.push_back(e.id);
	MeshFileWriter writer_(mesh, "reversed_mst_eids.vtk");
	writer_.WriteEdgesVtk(reversed_mst_eids);

    return 0;
}
