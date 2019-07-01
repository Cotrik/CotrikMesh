/*
 * TriangleSimplifier.cpp
 *
 *  Created on: Feb 05, 2019
 *      Author: cotrik
 */

#include "TriangleSimplifier.h"

TriangleSimplifier::TriangleSimplifier(Mesh& mesh) : Simplifier(mesh) {
    // TODO Auto-generated constructor stub

}

TriangleSimplifier::~TriangleSimplifier() {
    // TODO Auto-generated destructor stub
}

bool TriangleSimplifier::AllNeighboringVerticesRegular(const Vertex& v) {
    bool res = true;
    for (auto nvid : v.N_Vids) {
        auto& nv = mesh.V.at(nvid);
        if ((nv.isBoundary && nv.N_Fids.size() != 2) || (!nv.isBoundary && nv.N_Fids.size() != 4)) {
            res = false;
            break;
        }
    }
    return res;
}

bool TriangleSimplifier::AllNeighboringVerticesInterior(const Vertex& v) {
    bool res = true;
    for (auto nvid : v.N_Vids) {
        auto& nv = mesh.V.at(nvid);
        if (nv.isBoundary) {
            res = false;
            break;
        }
    }
    return res;
}

bool TriangleSimplifier::AllDiagonalVerticesSingular(const Vertex& v) {
    bool res = true;
    for (auto nfid : v.N_Fids) {
        auto nvid = get_diagnal_vid(v.id, nfid);
        auto& nv = mesh.V.at(nvid);
        if ((nv.isBoundary && nv.N_Fids.size() == 2) || (!nv.isBoundary && nv.N_Fids.size() == 4)) {
            res = false;
            break;
        }
    }
    return res;
}

void TriangleSimplifier::Run(std::set<size_t>& canceledFids) {
    for (auto& v : mesh.V) {
        if (v.isBoundary || v.N_Fids.size() != 3) continue;
        if (!AllNeighboringVerticesRegular(v) || !AllNeighboringVerticesInterior(v) || !AllDiagonalVerticesSingular(v)) continue;
        auto valence = GetCollapseValence(v);
        if (valence < minValence || valence > maxValence) continue;
        canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
        Collapse(v);
        break;
    }
}

void TriangleSimplifier::Collapse(const Vertex& v) {
    collapse_vids(v.N_Vids, v.id);
}

int TriangleSimplifier::GetCollapseValence(const Vertex& v) {
    int sum = 0;
    for (auto nfid : v.N_Fids) {
        auto nvid = get_diagnal_vid(v.id, nfid);
        auto& nv = mesh.V.at(nvid);
        sum += nv.N_Fids.size() - 4;
    }
    return sum + 3;
}
