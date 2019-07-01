/*
 * ConvexCut.cpp
 *
 *  Created on: Dec 8, 2017
 *      Author: cotrik
 */

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "DualMesh.h"
#include "ArgumentManager.h"
#include "Patches.h"
#include "Simplifier.h"
#include "MST.h"
#include <math.h>

const double angle = 155;

std::vector<std::vector<size_t>> get_boundary_link_vids(const Mesh& mesh);
std::vector<size_t> get_boundary_link_vids(const Mesh& mesh, const Vertex& start_v, std::vector<bool>& visited);
std::vector<std::vector<size_t>> get_corner_link_vids(const Mesh& mesh, const std::vector<std::vector<size_t>>& boundary_links);
std::vector<std::pair<size_t, size_t>> get_concave_corner_pairs(const Mesh& mesh, std::vector<std::vector<size_t>>& corner_link_vids);
int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: ConvexCut <file> <out.vtk>";
        return -1;
    }
    ArgumentManager argumentManager(argc, argv);
	MeshFileReader reader(argv[1]);
	auto mesh = reader.GetMesh();
	mesh.BuildAllConnectivities();
	mesh.ExtractBoundary();
	mesh.ExtractSingularities();
	Simplifier simplifier(mesh);
	simplifier.get_feature();

	auto boundary_links = get_boundary_link_vids(mesh);
	auto corner_links = get_corner_link_vids(mesh, boundary_links);
	auto pairs = get_concave_corner_pairs(mesh, corner_links);

    adjacency_list_t adjacency_list(mesh.V.size());
    for (auto& v : mesh.V)
        for (auto nvid : v.N_Vids) {
            adjacency_list[v.id].push_back(neighbor(nvid, /*1*/glm::length(v - mesh.V.at(nvid))));
        }
    Graph g;
	std::vector<std::vector<size_t>> cut_links;
	for (auto& p : pairs)
	    cut_links.push_back(g.GetShortestPath(adjacency_list, p.first, p.second));

	MeshFileWriter writer(mesh, argv[2]);
	writer.WriteLinksVtk(cut_links);
    return 0;
}

std::vector<std::vector<size_t>> get_boundary_link_vids(const Mesh& mesh) {
    std::vector<std::vector<size_t>> res;
    std::vector<bool> visited(mesh.V.size(), false);
    for (auto& v : mesh.V) {
        if (!v.isBoundary || visited[v.id]) continue;
        auto link = get_boundary_link_vids(mesh, v, visited);
        res.push_back(link);
    }
    return res;
}

std::vector<size_t> get_boundary_link_vids(const Mesh& mesh, const Vertex& start_v, std::vector<bool>& visited) {
    std::vector<size_t> res;
    visited[start_v.id] = true;
    res.push_back(start_v.id);
    auto next_vid = start_v.id;
    while (true) {
        auto old_size = res.size();
        auto& curr_v = mesh.V.at(next_vid);
        for (auto nvid : curr_v.N_Vids) {
            auto& v = mesh.V.at(nvid);
            if (!v.isBoundary || visited[v.id]) continue;
            visited[v.id] = true;
            res.push_back(v.id);
            next_vid = v.id;
            break;
        }
        if (old_size == res.size()) break;
    }
    return res;
}

std::vector<size_t> get_corner_link_vids(const Mesh& mesh, const std::vector<size_t>& boundary_link) {
    std::vector<size_t> corner_link;
    for (auto vid : boundary_link)
        if (mesh.V.at(vid).isCorner) corner_link.push_back(vid);
    return corner_link;
}

std::vector<std::vector<size_t>> get_corner_link_vids(const Mesh& mesh, const std::vector<std::vector<size_t>>& boundary_links) {
    std::vector<std::vector<size_t>> res;
    for (auto& boundary_link : boundary_links) {
        auto corner_link = get_corner_link_vids(mesh, boundary_link);
        if (!corner_link.empty()) res.push_back(corner_link);
    }
    return res;
}

std::vector<std::pair<size_t, size_t>> get_concave_corner_pairs(const Mesh& mesh, const std::vector<size_t>& corner_link) {
    std::vector<std::pair<size_t, size_t>> res;
    auto link = corner_link;
    if (link.size() > 4) {
        link.push_back(link.at(0));
        link.push_back(link.at(1));
        link.push_back(link.at(2));
    }
    for (auto i = 3; i < link.size(); ++i) {
        auto vid = link.at(i);
        auto& v = mesh.V.at(vid);
        auto prev_vid = link.at(i - 3);
        auto& prev_v = mesh.V.at(prev_vid);
        if (v.idealValence > 2 && prev_v.idealValence > 2) res.push_back(std::make_pair(vid, prev_vid));
    }
    return res;
}

std::vector<std::pair<size_t, size_t>> get_concave_corner_pairs(const Mesh& mesh, std::vector<std::vector<size_t>>& corner_links) {
    std::vector<std::pair<size_t, size_t>> res;
    for (auto& link : corner_links) {
        auto pairs = get_concave_corner_pairs(mesh, link);
        if (!pairs.empty()) std::copy(pairs.begin(), pairs.end(), std::back_inserter(res));
    }
    return res;
}
