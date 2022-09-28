#include <iostream>
#include <fstream>
#include <list>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <stack>
#include <algorithm>
#include "Mesh.h"
// C++ program for Kruskal's algorithm to find Minimum Spanning Tree of a given connected, undirected and weighted graph 
typedef std::pair<size_t, size_t> iPair;

typedef size_t vertex_t;
typedef double weight_t;

const weight_t max_weight = std::numeric_limits<double>::infinity();

struct DisjointSets {
	size_t *parent, *rnk;
	size_t n;

	DisjointSets(int n);
	size_t find(size_t u);
	void merge(size_t x, size_t y);
};

struct neighbor {
	vertex_t target;
	weight_t weight;
	neighbor(vertex_t arg_target, weight_t arg_weight)
		: target(arg_target), weight(arg_weight) {}
};

typedef std::vector<std::vector<neighbor> > adjacency_list_t;

struct Graph {
	size_t V, E;
	std::vector<std::pair<size_t, iPair> > edges;
	std::vector<iPair> mst_edges;
	Graph();
	Graph(size_t V, size_t E);
	void addEdge(size_t u, size_t v, size_t w);
	size_t kruskalMST();
	void prune(const Mesh& mesh, std::vector<size_t>& cut_graph_eids);
	void prune(const Mesh& mesh, std::set<size_t>& cut_graph_fids);
	void PruneLinkVids(const Mesh& mesh, std::vector<size_t>& linkVids);
	void PruneLinkVids(const Mesh& mesh, std::vector<std::vector<size_t>>& linkVids);
	void DijkstraComputePaths(vertex_t source, const adjacency_list_t &adjacency_list,
		std::vector<weight_t> &min_distance, std::vector<vertex_t> &previous);
	std::list<vertex_t> DijkstraGetShortestPathTo(vertex_t vertex, const std::vector<vertex_t> &previous);
	std::list<vertex_t> DijkstraGetShortestPath(const Mesh& mesh, const size_t src, const size_t dest);
	std::vector<size_t> GetShortestPath(const adjacency_list_t& adjacency_list, size_t src, size_t dest);
};

struct CycleExtractor {
	std::vector<std::vector<size_t>> graph;
	std::vector<std::vector<size_t>> cycles;
	std::vector<std::vector<size_t>> cycleVids;

	CycleExtractor(const Mesh& mesh, const std::vector<size_t>& eids);
	~CycleExtractor();
	void addEdge(size_t u, size_t v);
	// Function to mark the vertex with different colors for different cycles 
	void dfs_cycle(size_t u, size_t p, std::vector<size_t>& color, std::vector<size_t>& mark, std::vector<size_t>& par, 
		size_t& cyclenumber);
	std::vector<std::vector<size_t>> getCycleVids(size_t edges, std::vector<size_t>& mark, size_t& cyclenumber);
};
