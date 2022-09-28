#include "SimplificationOperation.h"

SimplificationOperation::SimplificationOperation() {}

SimplificationOperation::SimplificationOperation(Mesh& mesh_, MeshUtil& mu_, Smoother& smoother_) : mesh(&mesh_), mu(&mu_), smoother(&smoother_) {}

SimplificationOperation::~SimplificationOperation() {}

void SimplificationOperation::CheckValidity() {
    if (mesh == NULL) {
        std::cout << "No mesh to use for Simplification Operation." << std::endl;
        exit(0);
    }
    if (mu == NULL) {
        std::cout << "MeshUtil is not initialized for Simplification Operation." << std::endl;
        exit(0);
    }
    if (smoother == NULL) {
        std::cout << "Smoother is not intitalized for Simplification Operation." << std::endl;
        exit(0);
    }
    if (mesh->V.size() == 0 || mesh->F.size() == 0 || mesh->C.size() == 0) {
        std::cout << "No mesh to use for Simplification Operation." << std::endl;
        exit(0);
    }
}

void SimplificationOperation::SetMembers(Mesh& mesh_, MeshUtil& mu_, Smoother& smoother_) {
    mesh = &mesh_;
    mu = &mu_;
    smoother = &smoother_;
}

std::vector<size_t> SimplificationOperation::GetDifference(std::vector<size_t>& a, std::vector<size_t>& b) {
    return mu->GetDifference(a, b);
}

std::vector<size_t> SimplificationOperation::GetUnion(std::vector<size_t>& a, std::vector<size_t>& b) {
    return mu->GetUnion(a, b);
}

std::vector<size_t> SimplificationOperation::GetIntersection(std::vector<size_t>& a, std::vector<size_t>& b) {
    return mu->GetIntersection(a, b);
}

void SimplificationOperation::AddContents(std::vector<size_t>& a, std::vector<size_t>& b) {
    mu->AddContents(a, b);
}

void SimplificationOperation::UpdateContents(std::vector<size_t>& a, std::vector<size_t>& b) {
    mu->UpdateContents(a, b);
}

bool SimplificationOperation::IsCollapsable(size_t vid1, size_t vid2) {
    CheckValidity();

    bool res = true;

    auto& v0 = mesh->V.at(vid1);
    auto& v1 = mesh->V.at(vid2);
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

void SimplificationOperation::UpdateNeighborInfo(Vertex& target, Vertex& source, int fId) {
    CheckValidity();

    Face& face = mesh->F.at(fId);
    std::vector<size_t> verticesToRemove{source.id};
    std::vector<size_t> edgesToRemove;
    std::vector<size_t> edgesToKeep;
    std::vector<size_t> facesToRemove{(size_t) fId};
    for (auto eid: face.Eids) {
        Edge& e = mesh->E.at(eid);
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
    //     Edge& e = mesh->E.at(eid);
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
        Vertex& v = mesh->V.at(vid);
        AddContents(v.N_Vids, std::vector<size_t>{target.id});
        UpdateContents(v.N_Vids, verticesToRemove);
        UpdateContents(v.N_Eids, edgesToRemove);
        UpdateContents(v.N_Fids, facesToRemove);
    }
    // std::cout << "After updating source nvids" << std::endl;

    for (auto eid: diffSourceEdges) {
        Edge& e = mesh->E.at(eid);
        if (e.Vids.at(0) == source.id) e.Vids.at(0) = target.id;
        if (e.Vids.at(1) == source.id) e.Vids.at(1) = target.id;
    }
    // std::cout << "After updating source edges" << std::endl;

    for (auto fid: diffSourceFaces) {
        Face& f = mesh->F.at(fid);
        for (size_t i = 0; i < f.Vids.size(); i++) {
            if (f.Vids.at(i) == source.id) {
                f.Vids.at(i) = target.id;
                break;
            }
        }
        AddContents(f.N_Fids, diffTargetFaces);
        UpdateContents(f.N_Fids, facesToRemove);
    }
    // std::cout << "After updating source faces" << std::endl;
    
    AddContents(target.N_Vids, diffSourceVertices);
    AddContents(target.N_Eids, diffSourceEdges);
    AddContents(target.N_Fids, diffSourceFaces);
    UpdateContents(target.N_Fids, facesToRemove);
    // std::cout << "After updating target info" << std::endl;

    for (auto fid: diffTargetFaces) { 
        Face& f = mesh->F.at(fid);
        AddContents(f.N_Fids, diffSourceFaces);
        UpdateContents(f.N_Fids, facesToRemove);
    }
    // std::cout << "After updating target faces" << std::endl;


    for (auto vid: face.Vids) {
        if (vid == target.id || vid == source.id) continue;
        Vertex& v = mesh->V.at(vid);
        size_t e1, e2;
        for (auto eid: face.Eids) {
            Edge& e = mesh->E.at(eid);

            if (std::find(e.Vids.begin(), e.Vids.end(), target.id) != e.Vids.end() &&
                std::find(e.Vids.begin(), e.Vids.end(), v.id) != e.Vids.end()) {
                    e1 = e.id;
            } else if (std::find(e.Vids.begin(), e.Vids.end(), source.id) != e.Vids.end() && 
                std::find(e.Vids.begin(), e.Vids.end(), v.id) != e.Vids.end()) {
                    e2 = e.id;
            }
        }
        Edge& edgeToKeep = mesh->E.at(e1);
        Edge& edgeToRemove = mesh->E.at(e2);
        // std::cout << edgeToKeep.N_Fids.size() << " " << edgeToRemove.N_Fids.size() << std::endl;
        if (edgeToRemove.N_Fids.size() < 2) continue;
        // Face& faceTokeep = edgeToKeep.N_Fids.at(0) == face.id ? mesh->F.at(edgeToKeep.N_Fids.at(1)) : mesh->F.at(edgeToKeep.N_Fids.at(0));
        Face& faceToChange = edgeToRemove.N_Fids.at(0) == face.id ? mesh->F.at(edgeToRemove.N_Fids.at(1)) : mesh->F.at(edgeToRemove.N_Fids.at(0));
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
            Face& f = mesh->F.at(nfid);
            UpdateContents(f.N_Fids, facesToRemove);
        }
    }
    source.N_Vids.clear();
    source.N_Eids.clear();
    source.N_Fids.clear();
    // std::cout << "After updating face info" << std::endl;
    for (auto id: face.Vids) {
        FixDoublet(id);
    }
    for (auto id: face.Vids) {
        SetSingularity(id);
        auto& fv = mesh->V.at(id);
        AddContents(ToSmooth, std::vector<size_t>{id});
        AddContents(ToSmooth, fv.N_Vids);
        // ToSmooth.push_back(id);
        // ToSmooth.insert(ToSmooth.end(), fv.N_Vids.begin(), fv.N_Vids.end());
        // ToSmooth.insert(ToSmooth.end(), fv.N_Vids.begin(), fv.N_Vids.end());
    }
    
    face.N_Fids.clear();
    // face.Vids.clear();
    // face.Eids.clear();
    // Smooth();
}


void SimplificationOperation::FixDoublet(size_t vid) {
    CheckValidity();

    Vertex& v = mesh->V.at(vid);
    if (v.N_Vids.size() == 1) {
        v.N_Fids.clear();
        v.N_Vids.clear();
        v.N_Eids.clear();
        return;
    }
    if (v.N_Vids.size() != 2) return;
    if (v.type == FEATURE || v.isBoundary) return;
    for (auto fid: v.N_Fids) {
        if (mesh->F.at(fid).Vids.empty()) return;
    }

    // std::cout << "Removing doublet at: " << v.id << std::endl;

    Face& faceToRemove = mesh->F.at(v.N_Fids.at(0));
    size_t targetId = faceToRemove.Vids.at((std::distance(faceToRemove.Vids.begin(), std::find(faceToRemove.Vids.begin(), faceToRemove.Vids.end(), vid)) + 2) % faceToRemove.Vids.size());
    Vertex& target = mesh->V.at(targetId);
    Vertex& source = mesh->V.at(vid);
    UpdateNeighborInfo(target, source, faceToRemove.id);
    
    // for (auto fvid: faceToRemove.Vids) {
    //     FixDoublet(fvid);
    // }

    /*Face& faceToKeep = mesh->F.at(v.N_Fids.at(0));
    Face& faceToRemove = mesh->F.at(v.N_Fids.at(1));

    size_t vertexToAdd = faceToRemove.Vids.at((std::distance(faceToRemove.Vids.begin(), std::find(faceToRemove.Vids.begin(), faceToRemove.Vids.end(), vid)) + 2) % faceToRemove.Vids.size());
    faceToKeep.Vids.at(std::distance(faceToKeep.Vids.begin(), std::find(faceToKeep.Vids.begin(), faceToKeep.Vids.end(), vid))) = vertexToAdd;
    for (auto fvid: faceToRemove.Vids) {
        if (fvid == vid) continue;
        Vertex& fv = mesh->V.at(fvid);
        UpdateContents(fv.N_Vids, std::vector<size_t>{vid});
        UpdateContents(fv.N_Fids, std::vector<size_t>{faceToRemove.id});
        AddContents(fv.N_Fids, std::vector<size_t>{faceToKeep.id});
    }
    for (auto feid: faceToRemove.Eids) {
        Edge& fe = mesh->E.at(feid);
        if (std::find(fe.Vids.begin(), fe.Vids.end(), vid) != fe.Vids.end()) {
            size_t evid = fe.Vids.at(0) == vid ? fe.Vids.at(1) : fe.Vids.at(0);
            UpdateContents(mesh->V.at(evid).N_Eids, std::vector<size_t>{feid});
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
    faceToRemove.Eids.clear();
    faceToRemove.N_Fids.clear();
    for (auto fvid: faceToKeep.Vids) {
        SetSingularity(fvid);
    }
    for (auto fvid: faceToKeep.Vids) {
        FixDoublet(fvid);
    }*/
}

void SimplificationOperation::SetSingularity(size_t vid) {
    CheckValidity();

    auto& v = mesh->V.at(vid);
    v.isSingularity = v.N_Vids.size() == 4 ? false : true;
}

void SimplificationOperation::Smooth() {
    // return;
    int n = ToSmooth.size();
    for (int i = 0; i < n; i++) {
        auto& v = mesh->V.at(ToSmooth.at(i));
        AddContents(ToSmooth, v.N_Vids);
        // ToSmooth.insert(ToSmooth.end(), v.N_Vids.begin(), v.N_Vids.end());
    }
    smoother->Smooth(ToSmooth);
    return;
    int it = 0;
    while (it < 10) {
        std::vector<glm::dvec3> delta(ToSmooth.size());        
        for (int i = 0; i < ToSmooth.size(); i++) {
            auto& v = mesh->V.at(ToSmooth.at(i));
            if (v.N_Vids.empty() || v.N_Fids.empty() || v.type == FEATURE || v.isBoundary) continue;
            double n = 0.0;
            glm::dvec3 center(0.0, 0.0, 0.0);
            double w_agg = 0.0;
            for (auto nvid: v.N_Vids) {
                auto& nv = mesh->V.at(nvid);
                w_agg += glm::length(nv.xyz() - v.xyz());
            }
            for (auto nvid: v.N_Vids) {
                auto& nv = mesh->V.at(nvid);
                double w = 1.0 - (glm::length(nv.xyz() - v.xyz()) / w_agg);
                center += (w * nv.xyz());
                n += w;
            }
            delta.at(i) = center / n;
        }
        for (int i = 0; i < ToSmooth.size(); i++) {
            auto& v = mesh->V.at(ToSmooth.at(i));
            if (v.N_Vids.empty() || v.N_Fids.empty() || v.type == FEATURE || v.isBoundary) continue;
            v.xyz(delta.at(i));
        }
        it++;
    }
}

void SimplificationOperation::SetCoords(size_t vid, glm::dvec3& c) {
    CheckValidity();

    mesh->V.at(vid).x = c.x;
    mesh->V.at(vid).y = c.y;
    mesh->V.at(vid).z = c.z;
}

