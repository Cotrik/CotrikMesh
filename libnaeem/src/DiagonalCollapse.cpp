#include "DiagonalCollapse.h"

DiagonalCollapse::DiagonalCollapse() : SimplificationOperation() {}
DiagonalCollapse::DiagonalCollapse(Mesh& mesh_, size_t f, size_t d_idx1_, size_t d_idx2_) : SimplificationOperation(mesh_) {
    fId = f;
    d_idx1 = d_idx1_;
    d_idx2 = d_idx2_;
}

DiagonalCollapse::~DiagonalCollapse() {}

void DiagonalCollapse::SetRanking(MeshUtil& mu) {
    CheckValidity();

    Face& f = mesh.F.at(fId);

    Vertex& diag_v1 = mesh.V.at(f.Vids.at(d_idx1));
    Vertex& diag_v2 = mesh.V.at(f.Vids.at(d_idx2));
    Vertex& v3 = mesh.V.at(f.Vids.at((d_idx1 + 1) % f.Vids.size()));
    Vertex& v4 = mesh.V.at(f.Vids.at((d_idx2 + 1) % f.Vids.size()));

    double min = mu.GetVertexEnergy(diag_v1.id) + mu.GetVertexEnergy(diag_v2.id);
    double max = mu.GetVertexEnergy(v3.id) + mu.GetVertexEnergy(v4.id);
    double normalized_area = mu.GetMeshArea() / mu.GetFaceArea(fId);

    ranking = min / (max * normalized_area);
}

bool DiagonalCollapse::IsOperationValid() {
    if (mesh.F.at(fId).N_Fids.size() == 0) return false;
    return true;
}

void DiagonalCollapse::PerformOperation() {
    CheckValidity();

    if (!IsOperationValid()) return;

    Face& f = mesh.F.at(fId);

    Vertex& target = mesh.V.at(f.Vids.at(d_idx1));
    Vertex& source = mesh.V.at(f.Vids.at(d_idx2));
    Vertex& v3 = mesh.V.at(f.Vids.at((d_idx1 + 1) % f.Vids.size()));
    Vertex& v4 = mesh.V.at(f.Vids.at((d_idx2 + 1) % f.Vids.size()));

    // step 1: Set target location halfway between source and target
    target = 0.5 * (target.xyz() + source.xyz());

    // step 2: add source's neighboring vertices, edges and faces to target's neighbors
    UpdateNeighborInfo(target, source);

}

void DiagonalCollapse::UpdateNeighborInfo(Vertex& target, Vertex& source) {
    Face& face = mesh.F.at(fId);
    std::vector<size_t> verticesToRemove{source.id};
    std::vector<size_t> edgesToRemove;
    std::vector<size_t> edgesToKeep;
    std::vector<size_t> facesToRemove{fId};
    for (auto eid: face.Eids) {
        Edge& e = mesh.E.at(eid);
        if (std::find(e.Vids.begin(), e.Vids.end(), source.id) != e.Vids.end()) {
            edgesToRemove.push_back(eid);
        } else {
            edgesToKeep.push_back(eid);
        }
    }
    
    std::vector<size_t> diffSourceVertices = GetDifference(source.N_Vids, target.N_Vids);
    std::vector<size_t> diffSourceEdges = GetDifference(source.N_Eids, edgesToRemove);
    std::vector<size_t> diffSourceFaces = GetDifference(source.N_Fids, facesToRemove);

    std::vector<size_t> diffTargetVertices = GetDifference(target.N_Vids, source.N_Vids);
    std::vector<size_t> diffTargetEdges = GetDifference(target.N_Eids, edgesToKeep);
    std::vector<size_t> diffTargetFaces = GetDifference(source.N_Fids, facesToRemove);

    std::cout << "target: ";
    for (auto id: target.N_Vids) {
        std::cout << id << " ";
    }
    std::cout << std::endl;

    std::cout << "source: ";
    for (auto id: source.N_Vids) {
        std::cout << id << " ";
    }
    std::cout << std::endl;

    std::cout << "diffTargetVertices: ";
    for (auto id: diffTargetVertices) {
        std::cout << id << " ";
    }
    std::cout << std::endl;

    std::cout << "diffSourceVertices: ";
    for (auto id: diffSourceVertices) {
        std::cout << id << " ";
    }
    std::cout << std::endl;

    for (auto vid: source.N_Vids) {
        Vertex& v = mesh.V.at(vid);
        for (size_t i = 0; i < v.N_Vids.size(); i++) {
            if (v.N_Vids.at(i) == source.id) v.N_Vids.at(i) = target.id;
        }
    }

    for (auto eid: diffSourceEdges) {
        Edge& e = mesh.E.at(eid);
        if (e.Vids.at(0) == source.id) e.Vids.at(0) = target.id;
        if (e.Vids.at(1) == source.id) e.Vids.at(1) = target.id;
        e.N_Vids.insert(e.N_Vids.end(), diffTargetVertices.begin(), diffTargetVertices.end());
        e.N_Eids.insert(e.N_Eids.end(), target.N_Eids.begin(), target.N_Eids.end());
        UpdateContents(e.N_Eids, edgesToRemove);
    }

    for (auto fid: diffSourceFaces) {
        Face& f = mesh.F.at(fid);
        for (size_t i = 0; i < f.Vids.size(); i++) {
            if (f.Vids.at(i) == source.id) {
                f.Vids.at(i) = target.id;
                break;
            }
        }
        f.N_Vids.insert(f.N_Vids.end(), diffTargetVertices.begin(), diffTargetVertices.end());
        f.N_Eids.insert(f.N_Eids.end(), diffTargetEdges.begin(), diffTargetEdges.end());
        UpdateContents(f.N_Eids, edgesToRemove);
        for (auto nfid: diffTargetFaces) {
            if (std::find(f.N_Fids.begin(), f.N_Fids.end(), nfid) != f.N_Fids.end()) continue;
            f.N_Fids.push_back(nfid);
        }
        UpdateContents(f.N_Fids, facesToRemove);
    }

    target.N_Vids.insert(target.N_Vids.end(), diffSourceVertices.begin(), diffSourceVertices.end());
    target.N_Eids.insert(target.N_Eids.end(), diffSourceEdges.begin(), diffSourceEdges.end());
    target.N_Fids.insert(target.N_Fids.end(), diffSourceFaces.begin(), diffSourceFaces.end());
    UpdateContents(target.N_Fids, facesToRemove);
    
    for (auto eid: target.N_Eids) {
        Edge& e = mesh.E.at(eid);
        e.N_Vids.insert(e.N_Vids.end(), diffSourceVertices.begin(), diffSourceVertices.end());
        e.N_Eids.insert(e.N_Eids.end(), diffSourceEdges.begin(), diffSourceEdges.end());
        UpdateContents(e.N_Eids, edgesToRemove);
    }

    for (auto fid: target.N_Fids) {
        if (fid == face.id) continue; 
        Face& f = mesh.F.at(fid);
        f.N_Vids.insert(f.N_Vids.end(), diffSourceVertices.begin(), diffSourceVertices.end());
        f.N_Eids.insert(f.N_Eids.end(), diffSourceEdges.begin(), diffSourceEdges.end());
        UpdateContents(f.N_Eids, edgesToRemove);
        f.N_Fids.insert(f.N_Fids.end(), diffSourceFaces.begin(), diffSourceFaces.end());
        UpdateContents(f.N_Fids, facesToRemove);
    }

    for (auto vid: face.Vids) {
        if (vid == target.id || vid == source.id) continue;
        Vertex& v = mesh.V.at(vid);
        size_t e1, e2;
        for (auto eid: face.Eids) {
            Edge& e = mesh.E.at(eid);
            if (std::find(e.Vids.begin(), e.Vids.end(), target.id) != e.Vids.end() &&
                std::find(e.Vids.begin(), e.Vids.end(), v.id) != e.Vids.end()) {
                    e1 = e.id;
            } else if (std::find(e.Vids.begin(), e.Vids.end(), source.id) != e.Vids.end() && 
                std::find(e.Vids.begin(), e.Vids.end(), v.id) != e.Vids.end()) {
                    e2 = e.id;
            }
        }
        Edge& edgeToKeep = mesh.E.at(e1);
        Edge& edgeToRemove = mesh.E.at(e2);
        Face& faceTokeep = edgeToKeep.N_Fids.at(0) == face.id ? mesh.F.at(edgeToKeep.N_Fids.at(1)) : mesh.F.at(edgeToKeep.N_Fids.at(0));
        Face& faceToChange = edgeToRemove.N_Fids.at(0) == face.id ? mesh.F.at(edgeToRemove.N_Fids.at(1)) : mesh.F.at(edgeToRemove.N_Fids.at(0));
        for (int i = 0; i < faceToChange.Eids.size(); i++) {
            if (faceToChange.Eids.at(i) == edgeToRemove.id) {
                faceToChange.Eids.at(i) = edgeToKeep.id;
                break;
            }
        }
        for (int i = 0; i < edgeToKeep.N_Fids.size(); i++) {
            if (edgeToKeep.N_Fids.at(i) == face.id) {
                edgeToKeep.N_Fids.at(i) = faceToChange.id;
            }
        }
        UpdateContents(v.N_Vids, verticesToRemove);
        UpdateContents(v.N_Eids, edgesToRemove);
        UpdateContents(v.N_Fids, facesToRemove);
        for (auto neid: v.N_Eids) {
            Edge& e = mesh.E.at(neid);
            UpdateContents(e.N_Vids, verticesToRemove);
            UpdateContents(e.N_Eids, edgesToRemove);
        }
        for (auto nfid: v.N_Fids) {
            Face& f = mesh.F.at(nfid);
            UpdateContents(f.N_Vids, verticesToRemove);
            UpdateContents(f.N_Eids, edgesToRemove);
            UpdateContents(f.N_Fids, facesToRemove);
        }
        edgeToRemove.N_Vids.clear();
        edgeToRemove.N_Eids.clear();
        edgeToRemove.N_Fids.clear();
    }
    source.N_Vids.clear();
    source.N_Eids.clear();
    source.N_Fids.clear();
    face.N_Vids.clear();
    face.N_Eids.clear();
    face.N_Fids.clear();
}

std::vector<size_t> DiagonalCollapse::GetDifference(std::vector<size_t>& a, std::vector<size_t>& b) {
    std::vector<size_t> diff;
    std::set_difference(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(diff));
    return diff;
}

void DiagonalCollapse::UpdateContents(std::vector<size_t>& a, std::vector<size_t>& b) {
    std::vector<size_t> temp = GetDifference(a, b);
    a.clear();
    a.insert(a.begin(), temp.begin(), temp.end());
}

// std::cout << "";
// for (auto id: ) {
//     std::cout << id << " ";
// }
// std::cout << std::endl;