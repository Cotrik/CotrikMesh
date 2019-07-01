/*
 * Util.cpp
 *
 *  Created on: Feb 20, 2019
 *      Author: cotrik
 */

#include "Util.h"

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

