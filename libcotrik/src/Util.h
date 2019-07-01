/*
 * Util.h
 *
 *  Created on: Feb 20, 2019
 *      Author: cotrik
 */

#ifndef LIBCOTRIK_SRC_UTIL_H_
#define LIBCOTRIK_SRC_UTIL_H_

#include <unordered_map>
#include <unordered_set>
#include <set>
#include <map>
#include <vector>
#include <stack>
#include <queue>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <math.h>
#include <glm/glm.hpp>

enum FACE_TYPE {
    FACE_UNKNOWN = 0,
    FACE_X = 1,
    FACE_Y = 2,
    FACE_Z = 3,
};

class Util {
public:
    Util();
    virtual ~Util();
public:
    static std::set<size_t> get_intersect(const std::set<size_t>& s1, const std::set<size_t>& s2);
    static std::vector<std::vector<size_t>> combine(int n, int k);

    static FACE_TYPE GetFaceType(const glm::dvec3& normal);
private:
    static void combine(std::vector<size_t>& com, std::vector<std::vector<size_t>> &res, int n, int k, int start);
};

#endif /* LIBCOTRIK_SRC_UTIL_H_ */
