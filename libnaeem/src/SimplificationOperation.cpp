#include "SimplificationOperation.h"

SimplificationOperation::SimplificationOperation() {}

SimplificationOperation::SimplificationOperation(Mesh& mesh_, MeshUtil& mu_) : mesh(mesh_), mu(mu_) {}

SimplificationOperation::~SimplificationOperation() {}

void SimplificationOperation::CheckValidity() {
    if (mesh.V.size() == 0 || mesh.F.size() == 0 || mesh.C.size() == 0) {
        std::cout << "No mesh to use for Simplification Operation." << std::endl;
        exit(0);
    }
}

void SimplificationOperation::SetMesh(Mesh& mesh_) {
    mesh = mesh_;
}

std::vector<size_t> SimplificationOperation::GetDifference(std::vector<size_t>& a, std::vector<size_t>& b) {
    return mu.GetDifference(a, b);
}

std::vector<size_t> SimplificationOperation::GetUnion(std::vector<size_t>& a, std::vector<size_t>& b) {
    return mu.GetUnion(a, b);
}

std::vector<size_t> SimplificationOperation::GetIntersection(std::vector<size_t>& a, std::vector<size_t>& b) {
    return mu.GetIntersection(a, b);
}

void SimplificationOperation::AddContents(std::vector<size_t>& a, std::vector<size_t>& b) {
    mu.AddContents(a, b);
}

void SimplificationOperation::UpdateContents(std::vector<size_t>& a, std::vector<size_t>& b) {
    mu.UpdateContents(a, b);
}

bool SimplificationOperation::IsCollapsable(size_t vid1, size_t vid2) {
    CheckValidity();

    bool res = true;

    auto& v0 = mesh.V.at(vid1);
    auto& v1 = mesh.V.at(vid2);
    int count = 0;
    if (v0.isCorner) ++count;
    if (v1.isCorner) ++count;
    if (count > 1) res = false;
    std::set<size_t> labels;
    if (!v0.isCorner && v0.label != MAXID) labels.insert(v0.label);
    if (!v1.isCorner && v1.label != MAXID) labels.insert(v1.label);
    if (!v0.isBoundary && !v1.isBoundary && v0.N_Fids.size() + v1.N_Fids.size() - 4 > 6)
        res = false;
    if (!v0.isBoundary && !v1.isBoundary && v0.N_Fids.size() + v1.N_Fids.size() - 4 < 2/*Simplifier::minValence*/)
        res = false;
    if (v0.type == FEATURE && v1.type == FEATURE && v0.N_Fids.size() + v1.N_Fids.size() != 8/*Simplifier::minValence*/)
        res = false;
    if (v0.type == CORNER && v1.type == FEATURE && v0.N_Fids.size() + v1.N_Fids.size() != 2 + v0.idealValence/*Simplifier::minValence*/)
        res = false;
    if (v1.type == CORNER && v0.type == FEATURE && v0.N_Fids.size() + v1.N_Fids.size() != 2 + v1.idealValence/*Simplifier::minValence*/)
        res = false;

    // if (v0.isBoundary && !v1.isBoundary && v0.N_Fids.size() + v1.N_Fids.size() - 4 > 5)
    //     res = false;
    // if (!v0.isBoundary && v1.isBoundary && v0.N_Fids.size() + v1.N_Fids.size() - 4 > 5)
    //     res = false;
    // if (v0.isBoundary && !v1.isBoundary && v0.N_Fids.size() + v1.N_Fids.size() - 4 < 3)
    //     res = false;
    // if (!v0.isBoundary && v1.isBoundary && v0.N_Fids.size() + v1.N_Fids.size() - 4 < 3)
    //     res = false;

    // if ((v0.idealValence >= 3 || v1.idealValence >= 3) && v0.isBoundary && v1.isBoundary &&
    //         v0.N_Fids.size() + v1.N_Fids.size() - 2 < 3)
    //     res = false;

    if (labels.size() == 1 && v0.type == CORNER && v1.type == FEATURE && v0.labels.find(v1.label) == v0.labels.end())
        res = false;
    if (labels.size() == 1 && v0.type == FEATURE && v1.type == CORNER && v1.labels.find(v0.label) == v1.labels.end())
        res = false;
    if (labels.size() > 2)
        res = false;
    if (labels.size() == 2 && (v0.isCorner || v1.isCorner))
        res = false;

    return res;
}


void SimplificationOperation::FixDoublet(size_t vid) {
    CheckValidity();

    Vertex& v = mesh.V.at(vid);
    if (v.N_Vids.size() != 2) return;
    for (auto fid: v.N_Fids) {
        if (mesh.F.at(fid).Vids.empty()) return;
    }

    std::cout << "Removing doublet at: " << v.id << std::endl;

    Face& faceToKeep = mesh.F.at(v.N_Fids.at(0));
    Face& faceToRemove = mesh.F.at(v.N_Fids.at(1));

    size_t vertexToAdd = faceToRemove.Vids.at((std::distance(faceToRemove.Vids.begin(), std::find(faceToRemove.Vids.begin(), faceToRemove.Vids.end(), vid)) + 2) % faceToRemove.Vids.size());
    faceToKeep.Vids.at(std::distance(faceToKeep.Vids.begin(), std::find(faceToKeep.Vids.begin(), faceToKeep.Vids.end(), vid))) = vertexToAdd;
    for (auto fvid: faceToRemove.Vids) {
        if (fvid == vid) continue;
        Vertex& fv = mesh.V.at(fvid);
        UpdateContents(fv.N_Vids, std::vector<size_t>{vid});
        UpdateContents(fv.N_Fids, std::vector<size_t>{faceToRemove.id});
        AddContents(fv.N_Fids, std::vector<size_t>{faceToKeep.id});
    }
    for (auto feid: faceToRemove.Eids) {
        Edge& fe = mesh.E.at(feid);
        if (std::find(fe.Vids.begin(), fe.Vids.end(), vid) != fe.Vids.end()) {
            size_t evid = fe.Vids.at(0) == vid ? fe.Vids.at(1) : fe.Vids.at(0);
            UpdateContents(mesh.V.at(evid).N_Eids, std::vector<size_t>{feid});
            UpdateContents(faceToKeep.Eids, std::vector<size_t>{feid});
            fe.Vids.clear();
            fe.N_Fids.clear();
        } else {
            AddContents(faceToKeep.Eids, std::vector<size_t>{feid});
            UpdateContents(fe.N_Fids, std::vector<size_t>{faceToRemove.id});
            AddContents(fe.N_Fids, std::vector<size_t>{faceToKeep.id});
        }
    }
    v.N_Vids.clear();
    v.N_Eids.clear();
    v.N_Fids.clear();
    faceToRemove.Vids.clear();
    faceToRemove.N_Fids.clear();
    for (auto fvid: faceToKeep.Vids) {
        SetSingularity(fvid);
    }
    for (auto fvid: faceToKeep.Vids) {
        FixDoublet(fvid);
    }
}

void SimplificationOperation::SetSingularity(size_t vid) {
    CheckValidity();

    auto& v = mesh.V.at(vid);
    v.isSingularity = v.N_Vids.size() == 4 ? false : true;
}

