#include "MST.h"

DisjointSets::DisjointSets(int n) {
	this->n = n;
	parent = new size_t[n + 1];
	rnk = new size_t[n + 1];
	for (size_t i = 0; i <= n; i++) {
		rnk[i] = 0;
		parent[i] = i; // every element is parent of itself 
	}
}

// Find the parent of a node 'u' Path Compression 
size_t DisjointSets::find(size_t u) {
	/* Make the parent of the nodes in the path	from u--> parent[u] point to parent[u] */
	if (u != parent[u])	parent[u] = find(parent[u]);
	return parent[u];
}

// Union by rank 
void DisjointSets::merge(size_t x, size_t y) {
	x = find(x), y = find(y);
	/* Make tree with smaller height a subtree of the other tree */
	if (rnk[x] > rnk[y]) parent[y] = x;
	else parent[x] = y;
	if (rnk[x] == rnk[y]) rnk[y]++;
}


Graph::Graph() {

}

Graph::Graph(size_t V, size_t E) {
	this->V = V;
	this->E = E;
}

// Utility function to add an edge 
void Graph::addEdge(size_t u, size_t v, size_t w) {
	edges.push_back({ w,{ u, v } });
}

size_t Graph::kruskalMST() {
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

void Graph::prune(const Mesh& mesh, std::vector<size_t>& cut_graph_eids) {
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

void Graph::prune(const Mesh& mesh, std::set<size_t>& cut_graph_fids) {
	auto fids = cut_graph_fids;
	while (true) {
		std::map<size_t, std::set<size_t>> eid_fids;
		for (auto fid : fids) {
			auto& f = mesh.F.at(fid);
			for (auto eid : f.Eids)
				eid_fids[eid].insert(fid);
		}
		bool needPrune = false;
		//for (auto& item : eid_fids) {
		//	auto& e = mesh.E.at(item.first);
		//	if (!e.isBoundary && item.second.size() == 1) {
		//		needPrune = true;
		//		auto fid = *item.second.begin();
		//		auto iter = fids.find(fid);
		//		if (iter != fids.end()) fids.erase(fid);
		//	}
		//}
		std::vector<size_t> erase_fids;
		for (auto fid : fids) {
			auto count = 0;
			auto& f = mesh.F.at(fid);
			for (auto eid : f.Eids) {
				auto& e = mesh.E.at(eid);
				if (!e.isBoundary && eid_fids[eid].size() == 1) ++count;
			}
			if (count >= 1) {
				erase_fids.push_back(fid);
				needPrune = true;
			}
		}
		if (!needPrune) break;
		for (auto fid : erase_fids)
			fids.erase(fid);
	}

	cut_graph_fids = fids;
}

static std::vector<size_t> GetEids(const Mesh& mesh, const std::vector<size_t>& linkVids) {
	std::vector<size_t> eids;
	for (auto i = 1; i < linkVids.size(); ++i) {
		auto curr_vid = linkVids[i];
		auto prev_vid = linkVids[i - 1];
		auto& curr_v = mesh.V.at(curr_vid);
		auto& prev_v = mesh.V.at(prev_vid);
		for (auto eid : curr_v.N_Eids) {
			auto& e = mesh.E.at(eid);
			auto nvid = e.Vids[0] == curr_vid ? e.Vids[1] : e.Vids[0];
			if (nvid == prev_vid) {
				eids.push_back(eid);
				break;
			}
		}
	}
	return eids;
}

static std::vector<size_t> GetLinkEVids(const Mesh& mesh, const std::vector<size_t>& eids) {
	std::set<size_t> vids;
	for (auto eid : eids) {
		const auto& e = mesh.E.at(eid);
		vids.insert(e.Vids.begin(), e.Vids.end());
	}
	std::vector<size_t> res;
	for (auto vid : vids) res.push_back(vid);
	return res;
}

static std::vector<size_t> GetRingVids(const Mesh& mesh, const std::vector<size_t>& eids) {
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
			std::cerr << "Err in GetRingVids\n";
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

void Graph::PruneLinkVids(const Mesh& mesh, std::vector<size_t>& linkVids) {
	auto eids = GetEids(mesh, linkVids);
	prune(mesh, eids);
	linkVids = GetRingVids(mesh, eids);
}

void Graph::PruneLinkVids(const Mesh& mesh, std::vector<std::vector<size_t>>& linkVids) {
	for (auto& vids : linkVids)
		PruneLinkVids(mesh, vids);
}

void Graph::DijkstraComputePaths(vertex_t source,
	const adjacency_list_t &adjacency_list,
	std::vector<weight_t> &min_distance,
	std::vector<vertex_t> &previous) {
	int n = adjacency_list.size();
	min_distance.clear();
	min_distance.resize(n, max_weight);
	min_distance[source] = 0;
	previous.clear();
	previous.resize(n, -1);
	std::set<std::pair<weight_t, vertex_t> > vertex_queue;
	vertex_queue.insert(std::make_pair(min_distance[source], source));

	while (!vertex_queue.empty()) {
		weight_t dist = vertex_queue.begin()->first;
		vertex_t u = vertex_queue.begin()->second;
		vertex_queue.erase(vertex_queue.begin());

		// Visit each edge exiting u
		const std::vector<neighbor> &neighbors = adjacency_list[u];
		for (std::vector<neighbor>::const_iterator neighbor_iter = neighbors.begin();
			neighbor_iter != neighbors.end();
			neighbor_iter++) {
			vertex_t v = neighbor_iter->target;
			weight_t weight = neighbor_iter->weight;
			weight_t distance_through_u = dist + weight;
			if (distance_through_u < min_distance[v]) {
				vertex_queue.erase(std::make_pair(min_distance[v], v));

				min_distance[v] = distance_through_u;
				previous[v] = u;
				vertex_queue.insert(std::make_pair(min_distance[v], v));

			}

		}
	}
}

std::list<vertex_t> Graph::DijkstraGetShortestPathTo(
	vertex_t vertex, const std::vector<vertex_t> &previous) {
	std::list<vertex_t> path;
	for (; vertex != -1; vertex = previous[vertex])
		path.push_front(vertex);
	return path;
}

std::list<vertex_t> Graph::DijkstraGetShortestPath(const Mesh& mesh, const size_t src, const size_t dest) {
	adjacency_list_t adjacency_list(mesh.V.size());
	for (auto& v : mesh.V)
		for (auto nvid : v.N_Vids)
			adjacency_list[v.id].push_back(neighbor(nvid, /*1*/glm::length(v - mesh.V.at(nvid))));

	std::vector<weight_t> min_distance;
	std::vector<vertex_t> previous;
	DijkstraComputePaths(src, adjacency_list, min_distance, previous);
	std::list<vertex_t> path = DijkstraGetShortestPathTo(dest, previous);
	return path;
}

std::vector<size_t> Graph::GetShortestPath(const adjacency_list_t& adjacency_list, size_t src, size_t dest) {
    std::vector<weight_t> min_distance;
    std::vector<vertex_t> previous;
    DijkstraComputePaths(src, adjacency_list, min_distance, previous);
    std::list<vertex_t> path = DijkstraGetShortestPathTo(dest, previous);
    std::vector<size_t> path_vids(path.begin(), path.end());
    return path_vids;
}

CycleExtractor::CycleExtractor(const Mesh& mesh, const std::vector<size_t>& eids) {
	const auto N = mesh.V.size() > eids.size() ? mesh.V.size() + 1 : eids.size() + 1;
	graph.resize(N);
	cycles.resize(N);
	for (auto eid : eids)
		addEdge(mesh.E[eid].Vids[0] + 1, mesh.E[eid].Vids[1] + 1);

	std::vector<size_t> color(N, 0); // arrays required to color the graph, 
	std::vector<size_t> par(N, 0);   // store the parent of node 
	std::vector<size_t> mark(N, 0);  // mark with unique numbers 

	size_t cyclenumber = 0;                   // store the numbers of cycle 
	auto start_vid = mesh.E.at(eids.front()).Vids.front();
	dfs_cycle(start_vid, 0, color, mark, par, cyclenumber); // call DFS to mark the cycles 
	std::cout << "cyclenumber = " << cyclenumber << std::endl;
	auto res = getCycleVids(eids.size() - 1, mark, cyclenumber);
	graph.clear();
	cycles.clear();
}

CycleExtractor::~CycleExtractor() {
	graph.clear();
	cycles.clear();
	cycleVids.clear();
}

void CycleExtractor::addEdge(size_t u, size_t v) {
	graph[u].push_back(v);
	graph[v].push_back(u);
}

// Function to mark the vertex with different colors for different cycles 
void CycleExtractor::dfs_cycle(size_t u, size_t p, std::vector<size_t>& color, std::vector<size_t>& mark, std::vector<size_t>& par, 
	size_t& cyclenumber) {
	// already (completely) visited vertex. 
	if (color[u] == 2) return;

	// seen vertex, but was not completely visited -> cycle detected. backtrack based on parents to find the complete cycle. 
	if (color[u] == 1) {
		cyclenumber++;
		int cur = p;
		mark[cur] = cyclenumber;
		// backtrack the vertex which are in the current cycle thats found 
		while (cur != u) {
			cur = par[cur];
			mark[cur] = cyclenumber;
		}
		return;
	}
	par[u] = p;
	color[u] = 1; // partially visited. 

					// simple dfs on graph 
	for (auto v : graph[u]) {
		// if it has not been visited previously 
		if (v == par[u]) continue;
		dfs_cycle(v, u, color, mark, par, cyclenumber);
	}

	color[u] = 2; // completely visited. 
}

std::vector<std::vector<size_t>> CycleExtractor::getCycleVids(size_t edges, std::vector<size_t>& mark, size_t& cyclenumber) {
	for (size_t i = 1; i <= edges; ++i)
		if (mark[i] != 0) cycles[mark[i]].push_back(i);
	std::vector<std::vector<size_t>> res(cyclenumber);
	for (int i = 1; i <= cyclenumber; ++i) {
		res[i - 1] = cycles[i];
		for (auto& x : res[i - 1]) --x;
	}
	return res;
}
