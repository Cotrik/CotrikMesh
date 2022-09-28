/*
 * DoubletSimplifier.cpp
 *
 *  Created on: Jan 20, 2019
 *      Author: cotrik
 */

#include "DoubletSimplifier.h"

DoubletSimplifier::DoubletSimplifier(Mesh& mesh) : Simplifier(mesh) {
    // TODO Auto-generated constructor stub
}

DoubletSimplifier::~DoubletSimplifier() {
    // TODO Auto-generated destructor stub
}

void DoubletSimplifier::Run(std::set<size_t>& canceledFids) {
    for (auto& v : mesh.V) {
        if (v.N_Fids.size() != 2 || v.isBoundary) continue;
        auto& v0 = mesh.V.at(v.N_Vids[0]);
        auto& v1 = mesh.V.at(v.N_Vids[1]);

        if (checkCorner && v0.isCorner && v0.idealValence >= 3 && v0.N_Fids.size() <= 3) continue;
        if (checkCorner && v1.isCorner && v1.idealValence >= 3 && v1.N_Fids.size() <= 3) continue;

        auto& f0 = mesh.F.at(v.N_Fids.front());
        auto& f1 = mesh.F.at(v.N_Fids.back());
        std::set<size_t> eids_set(f0.Eids.begin(), f0.Eids.end());
        eids_set.insert(f1.Eids.begin(), f1.Eids.end());
        for (auto eid : v.N_Eids)
            eids_set.erase(eid);
        std::vector<size_t> eids;
        std::copy(eids_set.begin(), eids_set.end(), std::back_inserter(eids));
        auto linkVids = GetLinkVidsFromEids(eids);
        if (linkVids.empty()) continue;
        linkVids.pop_back();
        Face newf;
        newf.id = mesh.F.size();
        newf.Vids = linkVids;
        mesh.F.push_back(newf);
        canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
        return;
    }

    for (auto& v : mesh.V) {
        if (v.N_Fids.size() != 1 || !v.isBoundary) continue;
        auto& v0 = mesh.V.at(v.N_Vids[0]);
        auto& v1 = mesh.V.at(v.N_Vids[1]);
        //if (v0.type != FEATURE || v1.type != FEATURE || v0.label != v.label || v1.label != v.label) continue;
        auto target_vid = MAXID;
        if (v0.type == FEATURE && v1.type == FEATURE && v.type == FEATURE && (v0.label && v.label && v1.label == v.label)) {
            target_vid = v.id;
        } else if (v0.type == FEATURE && v1.type == FEATURE && v.type == CORNER && v.N_Fids.size() != 1 &&
            (v.labels.find(v0.label) != v.labels.end() && v.labels.find(v1.label) != v.labels.end())) {
            target_vid = v.id;
        } else if (v0.type == CORNER && v1.type == FEATURE && v.type == FEATURE &&
            (v0.labels.find(v.label) != v0.labels.end() && v1.label == v.label)) {
            target_vid = v0.id;
        } else if (v0.type == FEATURE && v1.type == CORNER && v.type == FEATURE &&
            (v1.labels.find(v.label) != v1.labels.end() && v0.label == v.label)) {
            target_vid = v1.id;
        }
        if (target_vid == MAXID) continue;
        auto vids = { v.id, v0.id, v1.id };
        for (auto vid : vids)
            Collapse(vid, target_vid);
        canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
        return;
    }
}


void DoubletSimplifier::RunCollective(std::set<size_t>& canceledFids) {
    for (auto& v : mesh.V) {
        if (v.N_Fids.size() != 2 || v.isBoundary) continue;
        auto& v0 = mesh.V.at(v.N_Vids[0]);
        auto& v1 = mesh.V.at(v.N_Vids[1]);

        if (checkCorner && v0.isCorner && v0.idealValence >= 3 && v0.N_Fids.size() <= 3) continue;
        if (checkCorner && v1.isCorner && v1.idealValence >= 3 && v1.N_Fids.size() <= 3) continue;

        auto& f0 = mesh.F.at(v.N_Fids.front());
        auto& f1 = mesh.F.at(v.N_Fids.back());
        std::set<size_t> eids_set(f0.Eids.begin(), f0.Eids.end());
        eids_set.insert(f1.Eids.begin(), f1.Eids.end());
        for (auto eid : v.N_Eids)
            eids_set.erase(eid);
        std::vector<size_t> eids;
        std::copy(eids_set.begin(), eids_set.end(), std::back_inserter(eids));
        auto linkVids = GetLinkVidsFromEids(eids);
        if (linkVids.empty()) continue;
        linkVids.pop_back();
        Face newf;
        newf.id = mesh.F.size();
        newf.Vids = linkVids;
        mesh.F.push_back(newf);
        canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
        // return;
    }
    return;
    for (auto& v : mesh.V) {
        if (v.N_Fids.size() != 1 || !v.isBoundary) continue;
        auto& v0 = mesh.V.at(v.N_Vids[0]);
        auto& v1 = mesh.V.at(v.N_Vids[1]);
        //if (v0.type != FEATURE || v1.type != FEATURE || v0.label != v.label || v1.label != v.label) continue;
        auto target_vid = MAXID;
        if (v0.type == FEATURE && v1.type == FEATURE && v.type == FEATURE && (v0.label && v.label && v1.label == v.label)) {
            target_vid = v.id;
        } else if (v0.type == FEATURE && v1.type == FEATURE && v.type == CORNER && v.N_Fids.size() != 1 &&
            (v.labels.find(v0.label) != v.labels.end() && v.labels.find(v1.label) != v.labels.end())) {
            target_vid = v.id;
        } else if (v0.type == CORNER && v1.type == FEATURE && v.type == FEATURE &&
            (v0.labels.find(v.label) != v0.labels.end() && v1.label == v.label)) {
            target_vid = v0.id;
        } else if (v0.type == FEATURE && v1.type == CORNER && v.type == FEATURE &&
            (v1.labels.find(v.label) != v1.labels.end() && v0.label == v.label)) {
            target_vid = v1.id;
        }
        if (target_vid == MAXID) continue;
        auto vids = { v.id, v0.id, v1.id };
        for (auto vid : vids)
            Collapse(vid, target_vid);
        canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
        return;
    }
}
