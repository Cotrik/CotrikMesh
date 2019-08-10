/*
 * MeshMap.cpp
 *
 *  Created on: Aug 05, 2019
 *      Author: cotrik
 */

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include <iostream>
#include <vector>
#include "ArgumentManager.h"
#include <math.h>

std::vector<glm::dvec2> GetPolarCoordinates(const Vertex& v0, const Vertex& v1, const std::vector<Vertex>& V) {
	std::vector<glm::dvec2> res(V.size());

	glm::dvec3 d0n = glm::normalize(v1 - v0);
	res[0] = glm::dvec2(0, 0);
	size_t i = 0;
	for (auto& v : V) {
		if (i != 0) {
			glm::dvec3 d1 = v - v0;
			auto d1n = glm::normalize(d1);
			res[i].x = glm::length(d1);                                              // radius
			res[i].y = atan2(glm::length(glm::cross(d0n, d1n)), glm::dot(d0n, d1n)); // theta
		}
		++i;
	}
	
	return res;
}

std::vector<glm::dvec2> GetXYCoordinates(const std::vector<glm::dvec2>& polarCoordinates, const std::vector<Vertex>& V) {
	std::vector<glm::dvec2> res(V.size());
	res[0] = glm::dvec2(0, 0);
	size_t i = 0;
	for (auto& v : V) {
		if (i != 0) {
			auto r = polarCoordinates[i].x;
			auto theta = polarCoordinates[i].y;
			res[i].x = r * cos(theta);
			res[i].y = r * sin(theta);
		}
		++i;
	}

	return res;
}

//std::vector<std::pair<std::pair<size_t, size_t>, glm::dvec2>> GetLocalPolarCoordinates(const Vertex& v0, const Vertex& v1, const std::vector<Vertex>& V) {
//	std::vector<std::pair<std::pair<size_t, size_t>, glm::dvec2>> res;//(V.size());
//	res.reserve(V.size());
//	std::vector<bool> visitedV(V.size());	           
//	std::queue<std::vector<size_t>> q({ {v0.id, v0.id, v1.id} }); // {vid, refVid, orientVid}
//	while (!q.empty()) {
//		auto n = q.size();
//		while (n--) {
//			auto vid_refVid_orientVid = q.front();
//			q.pop();
//			auto vid = vid_refVid_orientVid[0];
//			auto refVid = vid_refVid_orientVid[1];
//			auto orientVid = vid_refVid_orientVid[2];
//			auto& v = V.at(vid);
//			auto& v0 = V.at(refVid);
//			auto& v1 = V.at(orientVid);
//			if (!visitedV[vid]) {
//				visitedV[vid] = true;
//				for (auto nvid : v.N_Vids) {
//					if (visitedV[nvid]) continue;
//					q.push({ nvid, vid, v.N_Vids[0]});
//				}
//				glm::dvec3 d0n = glm::normalize(v1 - v0);
//				glm::dvec3 d1 = v - v0;
//				auto d1n = glm::normalize(d1);
//				std::pair<std::pair<size_t, size_t>, glm::dvec2> p;
//				p.second.x = glm::length(d1);                                              // radius
//				p.second.y = atan2(glm::length(glm::cross(d0n, d1n)), glm::dot(d0n, d1n)); // theta
//				p.first.first = vid;
//				p.first.second = refVid;
//				res.push_back(p);
//			}	
//		}
//	}
//
//	res[0] = { {0, 0}, glm::dvec2(0, 0) };
//
//	return res;
//}

std::vector<std::pair<std::vector<size_t>, glm::dvec2>> GetLocalPolarCoordinates(const Vertex& v0, const Vertex& v1, const std::vector<Vertex>& V) {
	std::cout << "v0.id = " << v0.id << " v1.id = " << v1.id << std::endl;
	std::vector<std::pair<std::vector<size_t>, glm::dvec2>> res;//(V.size());
	res.reserve(V.size());
	std::vector<bool> visitedV(V.size());
	auto d = glm::length(v1 - v0);
	res.push_back({ {0, 0}, glm::dvec2(0, 0) });
	res.push_back({ {v1.id, 0}, glm::dvec2(d, 0) });
	visitedV[v0.id] = true;
	visitedV[v1.id] = true;
	std::queue<std::vector<size_t>> q({ {v0.N_Vids[1], v0.id, v1.id} }); // {vid, refVid, orientVid}
	std::ofstream ofs("VisitedVertices.txt");
	while (!q.empty()) {
		auto n = q.size();
		while (n--) {
			auto vid_refVid_orientVid = q.front();
			q.pop();
			auto vid = vid_refVid_orientVid[0];
			auto refVid = vid_refVid_orientVid[1];
			auto orientVid = vid_refVid_orientVid[2];
			auto& v = V.at(vid);
			auto& v0 = V.at(refVid);
			auto& v1 = V.at(orientVid);
			if (!visitedV[vid] && visitedV[refVid] && visitedV[orientVid]) {
				std::cout << vid << "," << refVid << std::endl;
				ofs << 1 << " " << vid << std::endl;
				visitedV[vid] = true;
				for (auto nvid : v.N_Vids) {
					if (visitedV[nvid]/* || !visitedV[v.N_Vids[0]]*/) continue;
					for (auto nnvid : v.N_Vids)
						if (visitedV[nnvid] && nnvid != vid) {
							q.push({ nvid, vid, nnvid });
							break;
						}
					//q.push({ nvid, vid, v.N_Vids[0] });
				}
				glm::dvec3 d0n = glm::normalize(v1 - v0);
				glm::dvec3 d1 = v - v0;
				auto d1n = glm::normalize(d1);
				std::pair<std::vector<size_t>, glm::dvec2> p;
				p.second.x = glm::length(d1);                                              // radius
				p.second.y = atan2(glm::length(glm::cross(d0n, d1n)), glm::dot(d0n, d1n)); // theta
				p.first.resize(3);
				p.first[0] = vid;
				p.first[1] = refVid;
				p.first[2] = orientVid;
				res.push_back(p);
			}
			else if (!visitedV[vid] && (!visitedV[refVid] || !visitedV[orientVid])) {
				q.push(vid_refVid_orientVid);
			}
		}
	}
	if (res.size() != V.size()) {
		std::cerr << "Error in GetLocalPolarCoordinates\n";
		std::cerr << "res.size() = " << res.size() << " V.size() = " << V.size() << std::endl;
	}

	return res;
}

std::vector<glm::dvec2> GetLocalXYCoordinates(const std::vector<std::pair<std::vector<size_t>, glm::dvec2>>& polarCoordinates, const std::vector<Vertex>& V) {
	std::vector<glm::dvec2> res(V.size());
	res[0] = glm::dvec2(0, 0);
	size_t i = 0;
	const double PI = 3.1415926535;
	for (auto& polarCoordinate : polarCoordinates) {
		//if (i != 0) 
		{
			auto vid = polarCoordinate.first[0];
			auto refVid = polarCoordinate.first[1];
			auto orientVid = polarCoordinate.first[2];
			if (res[orientVid] == glm::dvec2(0, 0) || res[refVid] == glm::dvec2(0, 0)) {
				std::cerr << "refVid = " << refVid << " orientVid = " << orientVid << std::endl;
			}
			auto d0 = res[orientVid] - res[refVid];
			auto theta0 = atan2(d0.y, d0.x) + 2 * PI;
			auto r = polarCoordinate.second.x;
			auto theta = polarCoordinate.second.y;
			res[vid].x = r * cos(theta + theta0) + res[refVid].x;
			res[vid].y = r * sin(theta + theta0) + res[refVid].y;
		}
		++i;
	}

	return res;
}

std::vector<Vertex> Get2DVertices(const Mesh& mesh) {
	std::vector<Vertex> res(mesh.V.size());

	auto& v0 = mesh.V[0];
	auto& v1 = mesh.V[v0.N_Vids[0]];

	auto polarCoordinates = GetLocalPolarCoordinates(v0, v1, mesh.V);
	auto xyCoordinates = GetLocalXYCoordinates(polarCoordinates, mesh.V);
	size_t i = 0;
	for (auto& v : mesh.V) {
		res[i] = glm::dvec3(xyCoordinates[i].x, xyCoordinates[i].y, 0);
		++i;
	}
	return res;
}
int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: MeshMap open_surface_patch_file 2d_patch_file\n";
        return -1;
    }

    MeshFileReader reader(argv[1]);
    //reader.GetPointsScalarFields();
    Mesh& mesh = (Mesh&)reader.GetMesh();
	mesh.RemoveUselessVertices();
    mesh.BuildAllConnectivities();
    mesh.ExtractBoundary();
	mesh.GetAvgEdgeLength();

	auto newV = Get2DVertices(mesh);

	MeshFileWriter writer(newV, mesh.C, argv[2], mesh.m_cellType);
	writer.WriteFile();
    return 0;
}


