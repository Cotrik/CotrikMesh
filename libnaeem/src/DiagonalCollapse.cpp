#include "DiagonalCollapse.h"

DiagonalCollapse::DiagonalCollapse() : SimplificationOperation() {}
DiagonalCollapse::DiagonalCollapse(Mesh& mesh_, MeshUtil& mu_, Smoother& smoother_, int f, size_t d_idx1_, size_t d_idx2_) : SimplificationOperation(mesh_, mu_, smoother_) {
    fId = f;
    d_idx1 = d_idx1_;
    d_idx2 = d_idx2_;
}

DiagonalCollapse::~DiagonalCollapse() {}

void DiagonalCollapse::SetFaceInfo() {
    fId = GetIntersection(mesh.V.at(d_idx1).N_Fids, mesh.V.at(d_idx2).N_Fids)[0];
    auto& f = mesh.F.at(fId);
    d_idx1 = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), d_idx1));
    d_idx2 = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), d_idx2));
}

void DiagonalCollapse::SetRanking(glm::dvec3 d) {
    CheckValidity();

    Face& f = mesh.F.at(fId);

    Vertex& diag_v1 = mesh.V.at(f.Vids.at(d_idx1));
    Vertex& diag_v2 = mesh.V.at(f.Vids.at(d_idx2));
    Vertex& v3 = mesh.V.at(f.Vids.at((d_idx1 + 1) % f.Vids.size()));
    Vertex& v4 = mesh.V.at(f.Vids.at((d_idx2 + 1) % f.Vids.size()));

    double min = mu.GetVertexEnergy(diag_v1.id) + mu.GetVertexEnergy(diag_v2.id);
    double max = mu.GetVertexEnergy(v3.id) + mu.GetVertexEnergy(v4.id);
    double normalized_area = mu.GetFaceArea(fId) / mu.GetMeshArea();

    ranking = min / (max * normalized_area);
}

bool DiagonalCollapse::IsOperationValid() {
    CheckValidity();
    auto& f = mesh.F.at(fId);
    if (f.N_Fids.size() == 0 || f.Vids.empty()) return false;
    // if (mesh.V.at(f.Vids.at(d_idx1)).N_Fids.size() != 3 || mesh.V.at(f.Vids.at(d_idx2)).N_Fids.size() != 3) return false;
    return true;
}

void DiagonalCollapse::PerformOperation() {
    CheckValidity();
    if (fId == -1) {
        SetFaceInfo();
    } else {
        if (!IsOperationValid()) return;
    }
    Face& f = mesh.F.at(fId);

    Vertex& target = mesh.V.at(f.Vids.at(d_idx1));
    Vertex& source = mesh.V.at(f.Vids.at(d_idx2));
    Vertex& v3 = mesh.V.at(f.Vids.at((d_idx1 + 1) % f.Vids.size()));
    Vertex& v4 = mesh.V.at(f.Vids.at((d_idx2 + 1) % f.Vids.size()));

    // step 1: Set target location halfway between source and target
    target = 0.5 * (target.xyz() + source.xyz());

    // step 2: add source's neighboring vertices, edges and faces to target's neighbors
    // UpdateNeighborInfo(target, source);
    UpdateNeighborInfo(target, source, fId);
    ToSmooth.push_back(target.id);
    AddContents(ToSmooth, target.N_Vids);
    Smooth();
}

/*void DiagonalCollapse::UpdateNeighborInfo(Vertex& target, Vertex& source) {
    Face& face = mesh.F.at(fId);
    std::vector<size_t> verticesToRemove{source.id};
    std::vector<size_t> edgesToRemove;
    std::vector<size_t> edgesToKeep;
    std::vector<size_t> facesToRemove{(size_t) fId};
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
    std::vector<size_t> diffTargetFaces = GetDifference(target.N_Fids, facesToRemove);

    // std::cout << "target: " << target.id << std::endl;
    // std::cout << "source: " << source.id << std::endl;

    // for (auto fvid: face.Vids) {
    //     std::cout << "face v: " << fvid << std::endl;
    // }

    // for (auto eid: face.Eids) {
    //     Edge& e = mesh.E.at(eid);
    //     std::cout << "face e: " << e.Vids.at(0) << " " << e.Vids.at(1) << std::endl;
    // }

    // std::cout << "verticesToRemove: " << verticesToRemove.size() << std::endl;
    // std::cout << "edgesToRemove: " << edgesToRemove.size() << std::endl;
    // std::cout << "edgesToKeep: " << edgesToKeep.size() << std::endl;
    // std::cout << "diffSourceVertices: " << diffSourceVertices.size() << std::endl;
    // std::cout << "diffSourceEdges: " << diffSourceEdges.size() << std::endl;
    // std::cout << "diffSourceFaces: " << diffSourceFaces.size() << std::endl;
    // std::cout << "diffTargetVertices: " << diffTargetVertices.size() << std::endl;
    // std::cout << "diffTargetEdges: " << diffTargetEdges.size() << std::endl;
    // std::cout << "diffTargetFaces: " << diffTargetFaces.size() << std::endl;
    
    for (auto vid: source.N_Vids) {
        Vertex& v = mesh.V.at(vid);
        AddContents(v.N_Vids, std::vector<size_t>{target.id});
        UpdateContents(v.N_Vids, verticesToRemove);
        UpdateContents(v.N_Eids, edgesToRemove);
        UpdateContents(v.N_Fids, facesToRemove);
    }

    for (auto eid: diffSourceEdges) {
        Edge& e = mesh.E.at(eid);
        if (e.Vids.at(0) == source.id) e.Vids.at(0) = target.id;
        if (e.Vids.at(1) == source.id) e.Vids.at(1) = target.id;
    }

    for (auto fid: diffSourceFaces) {
        Face& f = mesh.F.at(fid);
        for (size_t i = 0; i < f.Vids.size(); i++) {
            if (f.Vids.at(i) == source.id) {
                f.Vids.at(i) = target.id;
                break;
            }
        }
        AddContents(f.N_Fids, diffTargetFaces);
        UpdateContents(f.N_Fids, facesToRemove);
    }

    AddContents(target.N_Vids, diffSourceVertices);
    AddContents(target.N_Eids, diffSourceEdges);
    AddContents(target.N_Fids, diffSourceFaces);
    UpdateContents(target.N_Fids, facesToRemove);
    
    for (auto fid: diffTargetFaces) { 
        Face& f = mesh.F.at(fid);
        AddContents(f.N_Fids, diffSourceFaces);
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
        for (auto nfid: v.N_Fids) {
            Face& f = mesh.F.at(nfid);
            UpdateContents(f.N_Fids, facesToRemove);
        }
    }
    source.N_Vids.clear();
    source.N_Eids.clear();
    source.N_Fids.clear();
    for (auto id: face.Vids) {
        FixDoublet(id);
    }
    for (auto id: face.Vids) {
        SetSingularity(id);
    }
    face.N_Fids.clear();
    // face.Vids.clear();
    // face.Eids.clear();
}*/
// std::cout << "";
// for (auto id: ) {
//     std::cout << id << " ";
// }
// std::cout << std::endl;