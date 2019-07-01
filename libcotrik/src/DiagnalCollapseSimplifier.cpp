/*
 * DiagnalCollapseSimplifier.cpp
 *
 *  Created on: Jan 20, 2019
 *      Author: cotrik
 */

#include "DiagnalCollapseSimplifier.h"

DiagnalCollapseSimplifier::DiagnalCollapseSimplifier(Mesh& mesh) : Simplifier(mesh) {
    // TODO Auto-generated constructor stub

}

DiagnalCollapseSimplifier::~DiagnalCollapseSimplifier() {
    // TODO Auto-generated destructor stub
}

void DiagnalCollapseSimplifier::Run3(std::set<size_t>& canceledFids) {
    for (auto& f : mesh.F) {
        auto& v0 = mesh.V.at(f.Vids[0]);
        auto& v1 = mesh.V.at(f.Vids[1]);
        auto& v2 = mesh.V.at(f.Vids[2]);
        auto& v3 = mesh.V.at(f.Vids[3]);
        if (v0.type == CORNER || v1.type == CORNER || v2.type == CORNER || v3.type == CORNER) continue;

        if (v0.type == FEATURE && v1.type == FEATURE && v2.type == FEATURE && v0.label == v1.label && v1.label == v2.label) {
            canceledFids.insert(f.id);
            Collapse(v0.id, v2.id);
            return;
        } else if (v1.type == FEATURE && v2.type == FEATURE && v3.type == FEATURE && v1.label == v2.label && v2.label == v3.label) {
            canceledFids.insert(f.id);
            Collapse(v1.id, v3.id);
            return;
        } else if (v2.type == FEATURE && v3.type == FEATURE && v0.type == FEATURE && v2.label == v3.label && v3.label == v0.label) {
            canceledFids.insert(f.id);
            Collapse(v2.id, v0.id);
            return;
        } else if (v3.type == FEATURE && v0.type == FEATURE && v1.type == FEATURE && v3.label == v0.label && v0.label == v1.label) {
            canceledFids.insert(f.id);
            Collapse(v3.id, v1.id);
            return;
        }
    }
}

void DiagnalCollapseSimplifier::Run(std::set<size_t>& canceledFids) {
    for (auto valence = Simplifier::maxValence - 1; valence > 4; --valence) {
        for (auto& f : mesh.F) {
            auto& v0 = mesh.V.at(f.Vids[0]);
            auto& v1 = mesh.V.at(f.Vids[1]);
            auto& v2 = mesh.V.at(f.Vids[2]);
            auto& v3 = mesh.V.at(f.Vids[3]);
            //if (v0.type == CORNER || v1.type == CORNER || v2.type == CORNER || v3.type == CORNER) continue;
            if (v1.N_Fids.size() >= valence && v3.N_Fids.size() >= valence && !v1.isBoundary && !v3.isBoundary) {
                if (v0.type == CORNER || v2.type == CORNER) continue;
                if (CanCollapseDiagnal(v0, v2)) {
                    Collapse(v0.id, v2.id);
                    canceledFids.insert(f.id);
                    return;
                } else if (CanCollapseDiagnal(v2, v0)) {
                    Collapse(v2.id, v0.id);
                    canceledFids.insert(f.id);
                    return;
                }
            }
            //else
                if (v0.N_Fids.size() >= valence && v2.N_Fids.size() >= valence && !v0.isBoundary && !v2.isBoundary) {
                if (v1.type == CORNER || v3.type == CORNER) continue;
                if (CanCollapseDiagnal(v1, v3)) {
                    Collapse(v1.id, v3.id);
                    canceledFids.insert(f.id);
                    return;
                } else if (CanCollapseDiagnal(v3, v1)) {
                    Collapse(v3.id, v1.id);
                    canceledFids.insert(f.id);
                    return;
                }
            }
        }
    }
}

void DiagnalCollapseSimplifier::CollapseDiagnal(const Vertex& v0, const Vertex& v2) {
    auto target_vid = v0.type == CORNER ? v0.id : v2.id;
    auto vid = v0.type == CORNER ? v2.id : v0.id;
    Collapse(vid, target_vid);
}

//bool DiagnalCollapseSimplifier::CanCollapseDiagnal(const Vertex& v, const Vertex& target_v) {
//    return !v.isBoundary && v.N_Fids.size() == 3 && target_v.N_Fids.size() == 3;
//}
bool DiagnalCollapseSimplifier::CanCollapseDiagnal(const Vertex& v, const Vertex& target_v) {
    return !v.isBoundary && v.N_Fids.size() == 3 && (target_v.N_Fids.size() >= Simplifier::minValence && target_v.N_Fids.size() <= Simplifier::maxValence);
}
