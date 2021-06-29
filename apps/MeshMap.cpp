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

const double PI = 3.1415926535;
size_t _v0id = 0;

glm::dvec3 GetNormal(const Vertex& v, const Vertex& v0, const Vertex& v1) {
	const glm::dvec3 d0 = v1.xyz() - v0.xyz();
	const glm::dvec3 d1 = v.xyz() - v0.xyz();
	return glm::normalize(glm::cross(d0, d1));
}

std::vector<glm::dvec2> GetPolarCoordinates(const Vertex& v0, const Vertex& v1, const std::vector<Vertex>& V) {
	std::vector<glm::dvec2> res(V.size());

	glm::dvec3 d0n = glm::normalize(v1 - v0);
	res[v0.id] = glm::dvec2(0, 0);
	for (auto& v : V) {
		glm::dvec3 d1 = v - v0;
		if (glm::length(d1) < 1e-8 || std::isnan(glm::length(d1))) continue;
		auto d1n = glm::normalize(d1);
		auto det = glm::dot(v.normal, glm::cross(d0n, d1n));
		auto dot = glm::dot(d0n, d1n);
		res[v.id].x = glm::length(d1); // radius
		res[v.id].y = atan2(det, dot); // theta
	}
	
	return res;
}

std::vector<glm::dvec2> GetXYCoordinates(const std::vector<glm::dvec2>& polarCoordinates, const std::vector<Vertex>& V) {
	std::vector<glm::dvec2> res(V.size());

	for (auto& v : V) {
		auto i = v.id;
		auto r = polarCoordinates[i].x;
		auto theta = polarCoordinates[i].y;
		res[i].x = r * cos(theta);
		res[i].y = r * sin(theta);
	}

	return res;
}

std::vector<std::pair<std::vector<size_t>, glm::dvec2>> GetLocalPolarCoordinates(const Vertex& v0, const Vertex& v1, const std::vector<Vertex>& V) {
	std::cout << "v0.id = " << v0.id << " v1.id = " << v1.id << std::endl;
	std::vector<std::pair<std::vector<size_t>, glm::dvec2>> res;//(V.size());
	res.reserve(V.size());
	std::vector<bool> visitedV(V.size());
	auto d = glm::length(v1 - v0);
	res.push_back({ {v0.id, v0.id, v1.id}, glm::dvec2(0, 0) });
	res.push_back({ {v1.id, v0.id, v1.id}, glm::dvec2(d, 0) });
	visitedV[v0.id] = true;
	visitedV[v1.id] = true;
	std::queue<std::vector<size_t>> q({ {v0.N_Vids[1], v0.id, v1.id} }); // {vid, refVid, orientVid}
	// std::ofstream ofs("VisitedVertices.txt");
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
				size_t nOrientvid = MAXID;
				for (auto nvid : v.N_Vids)
					if (visitedV[nvid]) {
						nOrientvid = nvid;
						break;
					}
				if (nOrientvid != MAXID) {
					for (auto nvid : v.N_Vids) {
						if (visitedV[nvid]) continue;
						q.push({ nvid, vid, nOrientvid });
					}
				}

				// https://math.stackexchange.com/questions/878785/how-to-find-an-angle-in-range0-360-between-2-vectors
				visitedV[vid] = true;
				glm::dvec3 d0n = glm::normalize(v1 - v0);
				glm::dvec3 d1 = v - v0;
				if (glm::length(d1) < 1e-8 || std::isnan(glm::length(d1))) continue;
				auto d1n = glm::normalize(d1);
				auto n = v.normal;
				// auto n = GetNormal(v, v0, v1);
				auto det = glm::dot(n, glm::cross(d0n, d1n));
				auto dot = glm::dot(d0n, d1n);
				std::pair<std::vector<size_t>, glm::dvec2> p;
				p.first = vid_refVid_orientVid;
				p.second.x = glm::length(d1);                // radius
				p.second.y = atan2(det, dot);                // theta
				res.push_back(p);

				//std::cout << "("<< vid << "," << refVid << "," << orientVid << "), (" 
				//	<< p.second.x << ", " << p.second.y * 180 / PI << "), ("
				//	<< v.normal.x << ", " << v.normal.y << ", " << v.normal.z << ")" << std::endl;
				//ofs << 1 << " " << vid << std::endl;
			}
		}
	}
	if (res.size() != V.size()) {
		std::cerr << "Error in GetLocalPolarCoordinates\n";
		std::cerr << "res.size() = " << res.size() << " V.size() = " << V.size() << std::endl;
	}

	return res;
}

std::vector<int> colors;

std::vector<std::pair<std::vector<size_t>, glm::dvec2>> GetLocalPolarCoordinates(const Vertex& v0, const Vertex& v1, const std::vector<Vertex>& V, const std::vector<Face>& F) {
	std::cout << "v0.id = " << v0.id << " v1.id = " << v1.id << std::endl;
	std::vector<std::pair<std::vector<size_t>, glm::dvec2>> res;//(V.size());
	res.reserve(V.size());
	std::vector<bool> visitedV(V.size());
	auto d = glm::length(v1 - v0);
	res.push_back({ { v0.id, v0.id, v1.id }, glm::dvec2(0, 0) });
	res.push_back({ { v1.id, v0.id, v1.id }, glm::dvec2(d, 0) });
	visitedV[v0.id] = true;
	visitedV[v1.id] = true;
	std::queue<std::vector<size_t>> q({ { v0.N_Vids[1], v0.id, v1.id } }); // {vid, refVid, orientVid}
																		   // std::ofstream ofs("VisitedVertices.txt");
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
				visitedV[vid] = true;
				for (auto nfid : v.N_Fids) {
					auto& f = F.at(nfid);
					size_t nTargetVid = MAXID;
					size_t nOrientVid = MAXID;
					int count = 0;
					for (auto nvid : f.Vids) {
						if (visitedV[nvid]) {
							++count;
							if (nvid != vid) nOrientVid = nvid;
						} else nTargetVid = nvid;
					}
					if (count == 2 && nTargetVid != MAXID && nOrientVid != MAXID) {
						q.push({ nTargetVid, vid, nOrientVid });
					}
				}

				// https://math.stackexchange.com/questions/878785/how-to-find-an-angle-in-range0-360-between-2-vectors
				glm::dvec3 d0n = glm::normalize(v1 - v0);
				glm::dvec3 d1 = v - v0;
				if (glm::length(d1) < 1e-8 || std::isnan(glm::length(d1))) continue;
				auto d1n = glm::normalize(d1);
				auto n = v.normal;
				// auto n = GetNormal(v, v0, v1);
				auto det = glm::dot(n, glm::cross(d0n, d1n));
				auto dot = glm::dot(d0n, d1n);
				std::pair<std::vector<size_t>, glm::dvec2> p;
				p.first = vid_refVid_orientVid;
				p.second.x = glm::length(d1);                // radius
				p.second.y = atan2(det, dot);                // theta
				res.push_back(p);

				//std::cout << "("<< vid << "," << refVid << "," << orientVid << "), (" 
				//	<< p.second.x << ", " << p.second.y * 180 / PI << "), ("
				//	<< v.normal.x << ", " << v.normal.y << ", " << v.normal.z << ")" << std::endl;
				//ofs << 1 << " " << vid << std::endl;
			}
		}
	}
	if (res.size() != V.size()) {
		std::cerr << "Error in GetLocalPolarCoordinates\n";
		std::cerr << "res.size() = " << res.size() << " V.size() = " << V.size() << std::endl;
	}

	return res;
}

void Get2DPlaneCoordinate(const Mesh& mesh, const std::pair<std::vector<size_t>, glm::dvec2>& polarCoordinate, std::vector<Vertex>& res) {
	auto vid = polarCoordinate.first[0];
	auto refVid = polarCoordinate.first[1];
	auto orientVid = polarCoordinate.first[2];
	auto d0 = res[orientVid] - res[refVid];
	auto theta0 = atan2(d0.y, d0.x);
	auto r = polarCoordinate.second.x;
	auto theta = polarCoordinate.second.y;
	res[vid].x = r * cos(theta + theta0) + res[refVid].x;
	res[vid].y = r * sin(theta + theta0) + res[refVid].y;
}

std::vector<Vertex> Map3DTo2DPlane_bak(const Mesh& mesh, const Vertex& _v0, const Vertex& _v1) {
	int colorid = 0;
	colors.resize(mesh.V.size());
	auto& V = mesh.V;
	auto& F = mesh.F;
	std::vector<Vertex> newV(V.size());
	std::cout << "v0.id = " << _v0.id << " v1.id = " << _v1.id << std::endl;
	std::vector<std::pair<std::vector<size_t>, glm::dvec2>> res;//(V.size());
	res.reserve(V.size());
	std::vector<bool> visitedV(V.size(), false);
	auto d = glm::length(_v1 - _v0);
	res.push_back({ { _v0.id, _v0.id, _v1.id }, glm::dvec2(0, 0) });
	res.push_back({ { _v1.id, _v0.id, _v1.id }, glm::dvec2(d, 0) });
	visitedV[_v0.id] = true;
	visitedV[_v1.id] = true;
	newV[_v0.id] = glm::dvec3(0, 0, 0);
	newV[_v1.id] = glm::dvec3(d, 0, 0);
	colors[_v0.id] = colorid++;
	colors[_v1.id] = colorid++;

	std::queue<std::vector<size_t>> q({ { _v0.N_Vids[1], _v0.id, _v1.id } }); // {vid, refVid, orientVid}
																			  // std::ofstream ofs("VisitedVertices.txt");
																			  //for (auto nfid : _v0.N_Fids) {
																			  //	auto& f = F.at(nfid);
																			  //	size_t nTargetVid = MAXID;
																			  //	size_t nOrientVid = MAXID;
																			  //	int count = 0;
																			  //	for (auto nvid : f.Vids) {
																			  //		if (visitedV[nvid]) {
																			  //			++count;
																			  //			if (nvid != _v0.id) nOrientVid = nvid;
																			  //		} else nTargetVid = nvid;
																			  //	}
																			  //	if (count == 2 && nTargetVid != MAXID && nOrientVid != MAXID) {
																			  //		q.push({ nTargetVid, _v0.id, nOrientVid });
																			  //		visitedV[nTargetVid] = true;
																			  //	}
																			  //}

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
				std::cout << "vid = " << vid << std::endl;
				visitedV[vid] = true;
				for (auto nfid : v.N_Fids) {
					auto& f = F.at(nfid);
					size_t nTargetVid = MAXID;
					size_t nOrientVid = MAXID;
					int count = 0;
					for (auto nvid : f.Vids) {
						if (visitedV[nvid]) {
							++count;
							if (nvid != vid) nOrientVid = nvid;
						} else nTargetVid = nvid;
					}
					if (count == 2 && nTargetVid != MAXID && nOrientVid != MAXID) {
						q.push({ nTargetVid, vid, nOrientVid });
					}
				}

				// https://math.stackexchange.com/questions/878785/how-to-find-an-angle-in-range0-360-between-2-vectors
				glm::dvec3 d0n = glm::normalize(v1 - v0);
				glm::dvec3 d1 = v - v0;
				if (glm::length(d1) < 1e-8 || std::isnan(glm::length(d1))) continue;
				auto d1n = glm::normalize(d1);
				auto n = v.normal;
				// auto n = GetNormal(v, v0, v1);
				auto det = glm::dot(n, glm::cross(d0n, d1n));
				auto dot = glm::dot(d0n, d1n);
				std::pair<std::vector<size_t>, glm::dvec2> p;
				p.first = vid_refVid_orientVid;
				p.second.x = glm::length(d1);                // radius
				p.second.y = atan2(det, dot);                // theta
				Get2DPlaneCoordinate(mesh, p, newV);
				auto distance_v_v0 = glm::length(v - v0);
				auto distance_newv_newv0 = glm::length(newV[vid] - newV[refVid]);
				if (fabs(distance_v_v0 - distance_newv_newv0) < 1e-12) {
					res.push_back(p);
					colors[v1.id] = colorid++;
				} else {
					std::cerr << "Parametrization is not correct!!" << "(" << vid << "," << refVid << "," << orientVid << "), ("
						<< p.second.x << ", " << p.second.y * 180 / PI << "), ("
						<< v.normal.x << ", " << v.normal.y << ", " << v.normal.z << ")" << std::endl;
					visitedV[vid] = false;
					q.push(vid_refVid_orientVid);
				}
			}
		}
	}
	if (res.size() != V.size()) {
		std::cerr << "Error in Map3DTo2DPlane\n";
		std::cerr << "res.size() = " << res.size() << " V.size() = " << V.size() << std::endl;
	}

	return newV;
}

std::vector<Vertex> Map3DTo2DPlane(Mesh& mesh, const Vertex& _v0, const Vertex& _v1) {
	int colorid = 0;
	colors.resize(mesh.V.size());
	auto& V = mesh.V;
	auto& F = mesh.F;
	std::vector<Vertex> newV(V.size());
	std::cout << "v0.id = " << _v0.id << " v1.id = " << _v1.id << std::endl;
	std::vector<std::pair<std::vector<size_t>, glm::dvec2>> res;//(V.size());
	res.reserve(V.size());
	std::vector<bool> visitedV(V.size(), false);
	auto d = glm::length(_v1 - _v0);
	res.push_back({ { _v0.id, _v0.id, _v1.id }, glm::dvec2(0, 0) });
	res.push_back({ { _v1.id, _v0.id, _v1.id }, glm::dvec2(d, 0) });
	visitedV[_v0.id] = true;
	visitedV[_v1.id] = true;
	newV[_v0.id] = glm::dvec3(0, 0, 0);
	newV[_v1.id] = glm::dvec3(d, 0, 0);
	colors[_v0.id] = colorid++;
	colors[_v1.id] = colorid++;

	auto avg_edge_length = mesh.GetAvgEdgeLength();

	std::queue<std::vector<size_t>> q;// ({ { _v0.N_Vids[1], _v0.id, _v1.id } }); // {vid, refVid, orientVid}
																		   // std::ofstream ofs("VisitedVertices.txt");
	for (auto nfid : _v0.N_Fids) {
		auto& f = F.at(nfid);
		size_t nTargetVid = MAXID;
		size_t nOrientVid = MAXID;
		int count = 0;
		for (auto nvid : f.Vids) {
			if (visitedV[nvid]) {
				++count;
				if (nvid != _v0.id) nOrientVid = nvid;
			} else nTargetVid = nvid;
		}
		if (count == 2 && nTargetVid != MAXID && nOrientVid != MAXID) {
			q.push({ nTargetVid, _v0.id, nOrientVid });
			//visitedV[nTargetVid] = true;
		}
	}

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
				// std::cout << "vid = " << vid << std::endl;
				visitedV[vid] = true;
				for (auto nfid : v.N_Fids) {
					auto& f = F.at(nfid);
					size_t nTargetVid = MAXID;
					size_t nOrientVid = MAXID;
					int count = 0;
					for (auto nvid : f.Vids) {
						if (visitedV[nvid]) {
							++count;
							if (nvid != vid) nOrientVid = nvid;
						} else nTargetVid = nvid;
					}
					if (count == 2 && nTargetVid != MAXID && nOrientVid != MAXID) {
						q.push({ nTargetVid, vid, nOrientVid });
					}
				}

				// https://math.stackexchange.com/questions/878785/how-to-find-an-angle-in-range0-360-between-2-vectors
				glm::dvec3 d0n = glm::normalize(v1 - v0);
				glm::dvec3 d1 = v - v0;
				if (glm::length(d1) < 1e-8 || std::isnan(glm::length(d1))) {
					std::cerr << "glm::length(d1) < 1e-8 || std::isnan(glm::length(d1))\n";
					continue;
				}
				auto d1n = glm::normalize(d1);
				auto n = v.normal;
				// auto n = GetNormal(v, v0, v1);
				auto det = glm::dot(n, glm::cross(d0n, d1n));
				auto dot = glm::dot(d0n, d1n);
				std::pair<std::vector<size_t>, glm::dvec2> p;
				p.first = vid_refVid_orientVid;
				p.second.x = glm::length(d1);                // radius
				p.second.y = atan2(det, dot);                // theta
				Get2DPlaneCoordinate(mesh, p, newV);
				auto distance_v_v0 = glm::length(v - v0);
				auto distance_newv_newv0 = glm::length(newV[vid] - newV[refVid]);
				auto distance_v_v1 = glm::length(v - v1);
				auto distance_newv_newv1 = glm::length(newV[vid] - newV[orientVid]);
				if (fabs(distance_v_v0 - distance_newv_newv0) < 5e-1 * avg_edge_length && fabs(distance_v_v1 - distance_newv_newv1) < 5e-1 * avg_edge_length) {
					res.push_back(p);
					colors[v.id] = colorid++;
				} else {
					std::cerr << "Parametrization is not correct!!" << "(" << vid << "," << refVid << "," << orientVid << "), ("
						<< p.second.x << ", " << p.second.y * 180 / PI << "), ("
						<< v.normal.x << ", " << v.normal.y << ", " << v.normal.z << ")" << std::endl;
					visitedV[vid] = false;
					// q.push(vid_refVid_orientVid);
				}
			}
		}
	}
	if (res.size() != V.size()) {
		std::cerr << "Error in Map3DTo2DPlane\n";
		std::cerr << "res.size() = " << res.size() << " V.size() = " << V.size() << std::endl;
	}

	return newV;
}

std::vector<glm::dvec2> GetLocalXYCoordinates(const std::vector<std::pair<std::vector<size_t>, glm::dvec2>>& polarCoordinates, const std::vector<Vertex>& V) {
	std::vector<glm::dvec2> res(V.size());
	res[0] = glm::dvec2(0, 0);
	size_t i = 0;
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
			auto theta0 = atan2(d0.y, d0.x);
			auto r = polarCoordinate.second.x;
			auto theta = polarCoordinate.second.y;
			res[vid].x = r * cos(theta + theta0) + res[refVid].x;
			res[vid].y = r * sin(theta + theta0) + res[refVid].y;
		}
		++i;
	}

	return res;
}

std::vector<Vertex> Get2DVertices(Mesh& mesh, const bool useLocalPolarCoordinate = true) {
	std::vector<Vertex> res(mesh.V.size());

	auto& v0 = mesh.V[_v0id];
	auto& v1 = mesh.V[v0.N_Vids[0]];

	if (useLocalPolarCoordinate) {
		/*auto polarCoordinates = GetLocalPolarCoordinates(v0, v1, mesh.V, mesh.F);
		auto xyCoordinates = GetLocalXYCoordinates(polarCoordinates, mesh.V);
		for (auto& v : mesh.V)
			res[v.id] = glm::dvec3(xyCoordinates[v.id].x, xyCoordinates[v.id].y, 0);*/
		res = Map3DTo2DPlane(mesh, v0, v1);
	} else {
		auto polarCoordinates = GetPolarCoordinates(v0, v1, mesh.V);
		auto xyCoordinates = GetXYCoordinates(polarCoordinates, mesh.V);
		for (auto& v : mesh.V)
			res[v.id] = glm::dvec3(xyCoordinates[v.id].x, xyCoordinates[v.id].y, 0);
	}

	return res;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: MeshMap open_surface_patch_file 2d_patch_file local=false v0id=0\n";
        return -1;
    }
	ArgumentManager am(argc, argv);
    MeshFileReader reader(argv[1]);
    //reader.GetPointsScalarFields();
    Mesh& mesh = (Mesh&)reader.GetMesh();
	mesh.RemoveUselessVertices();
    mesh.BuildAllConnectivities();
    mesh.ExtractBoundary();
	mesh.GetAvgEdgeLength();
	mesh.GetNormalOfSurfaceFaces();
	mesh.unifyOrientation();
	mesh.GetNormalOfSurfaceVertices();

	bool useLocalPolarCoordinate = false;
	am.get("local", useLocalPolarCoordinate);
	am.get("v0id", _v0id);
	auto newV = Get2DVertices(mesh, useLocalPolarCoordinate);
	{
		MeshFileWriter writer(newV, mesh.C, argv[2], mesh.m_cellType);
		writer.WriteFile();
		writer.WritePointData(colors, argv[2]);
	}
	{
		MeshFileWriter writer(mesh, "colors.vtk");
		writer.WriteFile();
		writer.WritePointData(colors, "colors.vtk");
	}
	{
		MeshFileWriter writer(mesh, "vertex_colors.vtk");
		writer.WriteVerticesVtk();
		writer.WritePointData(colors, "vertex_colors.vtk");
	}
    return 0;
}


