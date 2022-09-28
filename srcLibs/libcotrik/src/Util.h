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
#include "Mesh.h"

class Vertex;
class Edge;
class Face;
class Cell;
class Mesh;

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
	static double GetAngle(const Vertex& v, const Vertex& v0, const Vertex& v1);
	static std::vector<size_t> GetNeighborVids(const Mesh& mesh, const Vertex& v, size_t fid);
	static double GetAngle(const Mesh& mesh, const Vertex& v, const std::vector<size_t>& fids);

	static void set_redundent_clearn(std::vector<size_t>& set);
	static bool set_contain(std::vector<size_t>& large_set, size_t element);
	static void set_exclusion(std::vector<size_t>& large_set, std::vector<size_t>& small_set, std::vector<size_t> &result_set);
	static bool IsOverlap(const Face& f1, const Face& f2);
	static bool Find(const std::vector<size_t>& Ids, const size_t targetId);
	static size_t GetOppositeFaceId(const Mesh& mesh, const size_t cellId, const size_t faceId);
	static bool IsEdgeInCell(const Mesh& mesh, const size_t cellId, const size_t edgeId);
	static glm::dvec3 GetCenter(const std::vector<Vertex>& V);
	static void GetBoundingBox(const std::vector<Vertex>& V, glm::dvec3& Max, glm::dvec3& Min);

private:
    static void combine(std::vector<size_t>& com, std::vector<std::vector<size_t>> &res, int n, int k, int start);
};

#endif /* LIBCOTRIK_SRC_UTIL_H_ */
