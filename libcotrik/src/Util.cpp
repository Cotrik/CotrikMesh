/*
 * Util.cpp
 *
 *  Created on: Feb 20, 2019
 *      Author: cotrik
 */

#include "Util.h"

const double PI = 3.1415926535898;

Util::Util() {
    // TODO Auto-generated constructor stub

}

Util::~Util() {
    // TODO Auto-generated destructor stub
}

std::set<size_t> Util::get_intersect(const std::set<size_t>& s1, const std::set<size_t>& s2) {
    std::set<size_t> intersect;
    std::set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), std::inserter(intersect, intersect.begin()));
    return intersect;
}

std::vector<std::vector<size_t>> Util::combine(int n, int k) {
    std::vector<std::vector<size_t>> res;
    std::vector<size_t> com;
    combine(com, res, n, k, 0);
    return res;
}

void Util::combine(std::vector<size_t>& com, std::vector<std::vector<size_t>> &res, int n, int k, int start) {
    if (k == com.size()) {
        res.push_back(com);
        return;
    }
    for (int i = start; i < n; ++i) {
        com.push_back(i);
        combine(com, res, n, k, i + 1);
        com.pop_back();
    }
}

FACE_TYPE Util::GetFaceType(const glm::dvec3& normal) {
    const float l = glm::length(normal);
    if (fabs(normal.x / l) > 0.8) return FACE_X;
    else if (fabs(normal.y / l) > 0.8) return FACE_Y;
    else if (fabs(normal.z / l) > 0.8) return FACE_Z;

    FACE_TYPE faceType = FACE_X;
    double max = fabs(normal.x / l);
    if (fabs(normal.y / l) > max) {
        max = fabs(normal.y / l);
        faceType = FACE_Y;
    }
    if (fabs(normal.z / l) > max) {
        max = fabs(normal.z / l);
        faceType = FACE_Z;
    }
    return faceType;
}

double Util::GetAngle(const Vertex& v, const Vertex& v0, const Vertex& v1) {
	auto d0 = v0.xyz() - v.xyz();
	auto d1 = v1.xyz() - v.xyz();
	auto cosangle = glm::dot(glm::normalize(d0), glm::normalize(d1));
	return acos(cosangle) * 180.0 / PI;
}

std::vector<size_t> Util::GetNeighborVids(const Mesh& mesh, const Vertex& v, size_t fid) {
	std::vector<size_t> res;
	auto& f = mesh.F.at(fid);
	for (auto eid : f.Eids) {
		auto& e = mesh.E.at(eid);
		if (e.Vids[0] == v.id || e.Vids[1] == v.id) {
			auto vid = e.Vids[0] == v.id ? e.Vids[1] : e.Vids[0];
			res.push_back(vid);
		}
	}
	return res;
}

double Util::GetAngle(const Mesh& mesh, const Vertex& v, const std::vector<size_t>& fids) {
	double total_angle = 0;
	for (auto fid : fids) {
		auto vids = GetNeighborVids(mesh, v, fid);
		total_angle += GetAngle(v, mesh.V.at(vids[0]), mesh.V.at(vids[1]));
	}
	return total_angle;
}


void Util::set_redundent_clearn(std::vector<size_t>& set) {
	std::vector<size_t> set_copy;
	for (int i = 0; i < set.size(); i++) {
		bool have = false;
		for (size_t j = i + 1; j < set.size(); j++)
			if (set[i] == set[j]) have = true;
		if (!have) set_copy.push_back(set[i]);
	}
	set = set_copy;
}

bool Util::set_contain(std::vector<size_t>& large_set, size_t element) {
	for (size_t j = 0; j < large_set.size(); j++)
		if (element == large_set[j]) return true;

	return false;
}

void Util::set_exclusion(std::vector<size_t>& large_set, std::vector<size_t>& small_set, std::vector<size_t> &result_set) {
	result_set.clear();
	for (size_t i = 0; i < large_set.size(); i++) {
		bool inside = false;
		for (size_t j = 0; j < small_set.size(); j++) {
			if (small_set[j] == large_set[i]) {
				inside = true;
				break;
			}
		}
		if (!inside) result_set.push_back(large_set[i]);
	}
}

bool Util::IsOverlap(const Face& f1, const Face& f2) {
	bool bRet = false;
	for (size_t i = 0; i < f1.Vids.size(); i++) {
		for (size_t j = 0; j < f2.Vids.size(); j++) {
			if (f1.Vids.at(i) == f2.Vids.at(j)) {
				bRet = true;
				break;
			}
		}
	}
	return bRet;
}

bool Util::Find(const std::vector<size_t>& Ids, const size_t targetId) {
	bool isFound = false;
	for (size_t i = 0; i < Ids.size(); i++) {
		if (targetId == Ids.at(i)) {
			isFound = true;
			break;
		}
	}
	return isFound;
}

size_t Util::GetOppositeFaceId(const Mesh& mesh, const size_t cellId, const size_t faceId) {
	size_t oppositeFaceId = MAXID;
	const Face& face = mesh.F.at(faceId);
	const Cell& cell = mesh.C.at(cellId);
	for (size_t i = 0; i < cell.Fids.size(); i++) {
		const Face& cellFace = mesh.F.at(cell.Fids.at(i));
		if (!IsOverlap(face, cellFace)) {
			oppositeFaceId = cell.Fids.at(i);
			break;
		}
	}
	return oppositeFaceId;
}

bool Util::IsEdgeInCell(const Mesh& mesh, const size_t cellId, const size_t edgeId) {
	const Cell& cell = mesh.C.at(cellId);
	for (size_t i = 0; i < cell.Eids.size(); i++)
		if (edgeId == cell.Eids.at(i)) return true;
	return false;
}

void Util::GetBoundingBox(const std::vector<Vertex>& V, glm::dvec3& Max, glm::dvec3& Min) {
	Min = glm::dvec3(-1e+10);
	Max = glm::dvec3(1e+10);
	for (int i = 0; i < V.size(); i++) {
		for (int j = 0; j < 3; j++) {
			if (V[i][j] > Max[j]) Max[j] = V[i][j];
			if (V[i][j] < Min[j]) Min[j] = V[i][j];
		}
	}
}

glm::dvec3 Util::GetCenter(const std::vector<Vertex>& V) {
	glm::dvec3 sum(0.0, 0.0, 0.0);
	for (size_t i = 0; i < V.size(); i++)
		sum += V.at(i).xyz();
	return glm::dvec3(sum.x / V.size(), sum.y / V.size(), sum.z / V.size());
}
