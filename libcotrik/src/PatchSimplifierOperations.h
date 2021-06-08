/*
 * PatchSimplifierOperations.h
 *
 *  Created on: June 5, 2021
 *      Author: Naeem
 */

#ifndef LIBCOTRIK_SRC_PATCH_SIMPLIFIER_OPERATIONS_H_
#define LIBCOTRIK_SRC_PATCH_SIMPLIFIER_OPERATIONS_H_

#include "Mesh.h"

#include <unordered_set>
#include <set>
#include <iostream>

struct SimplificationOperation {
    std::string type;
    double profitability = 0;
    std::set<size_t> canceledFids;
    std::vector<Face> newFaces;
    std::vector<size_t> updateVertexIds;
    std::vector<glm::dvec3> updatedVertexPos;
};

#endif /* LIBCOTRIK_SRC_PATCH_SIMPLIFIER_OPERATIONS_H_ */
