#include "VertexSplit.h"

VertexSplit::VertexSplit() {}
VertexSplit::VertexSplit(Mesh& mesh_, MeshUtil& mu_, Smoother& smoother_, size_t vid_, std::vector<size_t> edgesToSplit_) : SimplificationOperation(mesh_, mu_, smoother_) {
    vid = vid_;
    splitEdges = edgesToSplit_;
}
VertexSplit::~VertexSplit() {}

bool VertexSplit::IsOperationValid() {
    CheckValidity();

    auto& v = mesh->V.at(vid);
    int count = 0;
    if (v.type == FEATURE) {
        for (auto nvid: v.N_Vids) {
            if (mesh->V.at(nvid).type == FEATURE) count += 1;
        }
        if (count < 2 || count > 2) return false;
    }
    if (v.isBoundary) {
        for (auto nvid: v.N_Vids) {
            if (mesh->V.at(nvid).isBoundary) count += 1;
        }
        if (count < 2 || count > 2) return false;
    }
    return true;
}

void VertexSplit::SetRanking(glm::dvec3 d) {
    CheckValidity();

}

void VertexSplit::PerformOperation() {
    CheckValidity();
    if (!IsOperationValid()) return;

    if (mesh->V.at(vid).type != FEATURE && !mesh->V.at(vid).isBoundary && mesh->V.at(vid).N_Vids.size() == 5) {
        FiveVertexSplit();
        return;
    }

    if (splitEdges.empty()) {
       splitEdges = SelectEdgesToSplit();
    }

    Edge& e1 = mesh->E.at(splitEdges.at(0));
    Edge& e2 = mesh->E.at(splitEdges.at(1));
    
    Vertex new_V1(0.5 * (mesh->V.at(e1.Vids.at(0)).xyz() + mesh->V.at(e1.Vids.at(1)).xyz()));
    Vertex new_V2(0.5 * (mesh->V.at(e2.Vids.at(0)).xyz() + mesh->V.at(e2.Vids.at(1)).xyz()));
    new_V1.id = mesh->V.size();
    mesh->V.push_back(new_V1);
    new_V2.id = mesh->V.size();
    mesh->V.push_back(new_V2);

    auto& v = mesh->V.at(vid);
    auto& newV1 = mesh->V.at(new_V1.id);
    auto& newV2 = mesh->V.at(new_V2.id);
    newV1.type = v.type;
    newV2.type = v.type;
    newV1.isBoundary = v.isBoundary;
    newV2.isBoundary = v.isBoundary;
    newV1.idealValence = v.idealValence;
    newV2.idealValence = v.idealValence;

    std::vector<size_t> facesToInspect;
    AddContents(facesToInspect, v.N_Fids);

    std::vector<size_t> edgesToChange = {};
    std::vector<size_t> facesToRemove = GetDifference(v.N_Fids, e1.N_Fids);
    for (auto fid: e1.N_Fids) {
        auto& f = mesh->F.at(fid);
        for (int i = 0; i < f.Vids.size(); i++) {
            if (f.Vids.at(i) == vid) {
                f.Vids.at(i) = newV1.id;
                break;    
            }
        }
        AddContents(edgesToChange, GetIntersection(v.N_Eids, f.Eids));
        UpdateContents(f.N_Fids, facesToRemove);
        AddContents(newV1.N_Fids, std::vector<size_t>{f.id});
    }

    std::vector<size_t> edgesToAvoid = GetDifference(v.N_Eids, edgesToChange);
    std::vector<size_t> facesToAvoid = {};
    for (auto eid: edgesToChange) {
        auto& e = mesh->E.at(eid);
        std::vector<size_t> edgeV;
        e.Vids.at(0) == vid ? edgeV.push_back(e.Vids.at(1)) : edgeV.push_back(e.Vids.at(0));
        e.Vids.at(0) == vid ? e.Vids.at(0) = newV1.id : e.Vids.at(1) = newV1.id;
        
        AddContents(facesToAvoid, GetDifference(e.N_Fids, e1.N_Fids));
        AddContents(newV1.N_Eids, std::vector<size_t>{eid});
        AddContents(newV1.N_Vids, edgeV);
        UpdateContents(mesh->V.at(edgeV.at(0)).N_Vids, std::vector<size_t>{vid});
        AddContents(mesh->V.at(edgeV.at(0)).N_Vids, std::vector<size_t>{newV1.id});
    }
    // std::cout << "faces to avoid: " << facesToAvoid.size() << std::endl;
    // return;

    std::vector<size_t> facesToChange = {};
    for (auto eid: edgesToAvoid) {
        auto& e = mesh->E.at(eid);
        std::vector<size_t> edgeV;
        e.Vids.at(0) == vid ? edgeV.push_back(e.Vids.at(1)) : edgeV.push_back(e.Vids.at(0));
        e.Vids.at(0) == vid ? e.Vids.at(0) = newV2.id : e.Vids.at(1) = newV2.id;
        AddContents(facesToChange, GetDifference(e.N_Fids, facesToAvoid));
        AddContents(newV2.N_Eids, std::vector<size_t>{eid});
        AddContents(newV2.N_Vids, edgeV);
        UpdateContents(mesh->V.at(edgeV.at(0)).N_Vids, std::vector<size_t>{vid});
        AddContents(mesh->V.at(edgeV.at(0)).N_Vids, std::vector<size_t>{newV2.id});
    }

    facesToRemove = GetDifference(v.N_Fids, facesToChange);
    for (auto fid: facesToChange) {
        auto& f = mesh->F.at(fid);
        for (int i = 0; i < f.Vids.size(); i++) {
            if (f.Vids.at(i) == vid) {
                f.Vids.at(i) = newV2.id;
                break;    
            }
        }
        UpdateContents(f.N_Fids, facesToRemove);
        AddContents(newV2.N_Fids, std::vector<size_t>{f.id});
    }

    int offset = 0;
    std::vector<Face> newFaces;
    // std::cout << facesToAvoid.size() << std::endl;
    for (auto fid: facesToAvoid) {
        auto& f = mesh->F.at(fid);
        for (auto eid: f.Eids) {
            auto& e = mesh->E.at(eid);
            if (std::find(v.N_Eids.begin(), v.N_Eids.end(), eid) != v.N_Eids.end()) {
                if (e.Vids.at(0) == newV1.id || e.Vids.at(1) == newV1.id) {
                    size_t evid = e.Vids.at(0) == newV1.id ? e.Vids.at(1) : e.Vids.at(0);
                    int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), evid));
                    if (f.Vids.at((idx+1)%f.Vids.size()) == vid) {
                        std::vector<size_t> newVids = {evid, newV1.id, vid, f.Vids.at((idx+3)%f.Vids.size())};
                        Face newF(newVids);
                        newFaces.push_back(newF);
                    } else if (f.Vids.at((idx+3)%f.Vids.size()) == vid) {
                        std::vector<size_t> newVids = {evid, f.Vids.at((idx+1)%f.Vids.size()), vid, newV1.id};
                        Face newF(newVids);
                        newFaces.push_back(newF);
                    }
                } else if (e.Vids.at(0) == newV2.id || e.Vids.at(1) == newV2.id) {
                    size_t evid = e.Vids.at(0) == newV2.id ? e.Vids.at(1) : e.Vids.at(0);
                    int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), evid));
                    if (f.Vids.at((idx+1)%f.Vids.size()) == vid) {
                        std::vector<size_t> newVids = {evid, newV2.id, vid, f.Vids.at((idx+3)%f.Vids.size())};
                        Face newF(newVids);
                        newFaces.push_back(newF);
                    } else if (f.Vids.at((idx+3)%f.Vids.size()) == vid) {
                        std::vector<size_t> newVids = {evid, f.Vids.at((idx+1)%f.Vids.size()), vid, newV2.id};
                        Face newF(newVids);
                        newFaces.push_back(newF);
                    }
                }
            }
            UpdateContents(e.N_Fids, std::vector<size_t>{fid});
        }
    }
    std::vector<Edge> newEdges;
    newEdges.emplace_back(std::vector<size_t>{vid, newV1.id});
    newEdges.emplace_back(std::vector<size_t>{vid, newV2.id});
    for (auto fid: facesToAvoid) {
        auto& f = mesh->F.at(fid);
        int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), vid));
        newEdges.emplace_back(std::vector<size_t>{vid, f.Vids.at((idx+2)%f.Vids.size())});
        for (auto fvid: f.Vids) {
            UpdateContents(mesh->V.at(fvid).N_Fids, std::vector<size_t>{f.id});
        }
    }

    v.N_Vids.clear();
    v.N_Eids.clear();
    v.N_Fids.clear();
    std::vector<size_t> edgesToUpdate;
    for (auto& e: newEdges) {
        e.id = mesh->E.size();
        mesh->E.push_back(e);
        edgesToUpdate.push_back(e.id);
    }
    std::vector<size_t> facesToUpdate;
    for (auto& f: newFaces) {
        f.id = mesh->F.size();
        mesh->F.push_back(f);
        facesToUpdate.push_back(f.id);
    }

    for (auto eid: edgesToUpdate) {
        auto& e = mesh->E.at(eid);
        AddContents(mesh->V.at(e.Vids.at(0)).N_Vids, std::vector<size_t>{e.Vids.at(1)});
        AddContents(mesh->V.at(e.Vids.at(1)).N_Vids, std::vector<size_t>{e.Vids.at(0)});
        AddContents(mesh->V.at(e.Vids.at(0)).N_Eids, std::vector<size_t>{e.id});
        AddContents(mesh->V.at(e.Vids.at(1)).N_Eids, std::vector<size_t>{e.id});
        AddContents(v.N_Eids, std::vector<size_t>{e.id});
    }
    for (auto fid: facesToUpdate) {
        auto& f = mesh->F.at(fid);
        for (auto fid: facesToAvoid) {
            auto& nf = mesh->F.at(fid);
            for (auto eid: nf.Eids) {
                auto& e = mesh->E.at(eid);
                if (std::find(f.Vids.begin(), f.Vids.end(), e.Vids.at(0)) != f.Vids.end() && 
                    std::find(f.Vids.begin(), f.Vids.end(), e.Vids.at(1)) != f.Vids.end()) {
                        AddContents(f.Eids, std::vector<size_t>{e.id});
                        AddContents(e.N_Fids, std::vector<size_t>{f.id});
                }
            }
        }
        for (auto eid: edgesToUpdate) {
        // for (auto& e: newEdges) {
            auto& e = mesh->E.at(eid);
            if (std::find(f.Vids.begin(), f.Vids.end(), e.Vids.at(0)) != f.Vids.end() && 
                std::find(f.Vids.begin(), f.Vids.end(), e.Vids.at(1)) != f.Vids.end()) {
                    AddContents(f.Eids, std::vector<size_t>{e.id});
                    AddContents(e.N_Fids, std::vector<size_t>{f.id});
            }
        }
        for (auto fvid: f.Vids) {
            auto& fv = mesh->V.at(fvid);
            UpdateContents(fv.N_Fids, facesToAvoid);
            AddContents(fv.N_Fids, std::vector<size_t>{f.id});
        }
        AddContents(v.N_Fids, std::vector<size_t>{f.id});
        for (auto fvid: f.Vids) {
            AddContents(f.N_Fids, GetDifference(mesh->V.at(fvid).N_Fids, std::vector<size_t>{f.id}));
        }
    }
    AddContents(facesToInspect, facesToUpdate);
    UpdateContents(facesToInspect, facesToRemove);
    std::vector<size_t> verticesToCheck;
    for (auto fid: facesToInspect) {
        auto& f = mesh->F.at(fid);
        AddContents(verticesToCheck, std::vector<size_t>(f.Vids.begin(), f.Vids.end()));
    }
    for (auto id: verticesToCheck) {
        SetSingularity(id);
    }
    for (auto fid: facesToAvoid) {
        auto& f = mesh->F.at(fid);
        f.Vids.clear();
        f.Eids.clear();
        f.N_Fids.clear();
    }
    AddContents(ToSmooth, verticesToCheck);
    /*std::cout << "checking validity of nvids and faces" << std::endl;
    std::cout << "v nvids #: " << v.N_Vids.size() << std::endl;
    std::cout << "v neids #: " << v.N_Eids.size() << std::endl;
    std::cout << "v nfids #: " << v.N_Fids.size() << std::endl;
    for (auto fid: v.N_Fids) {
        auto& f = mesh->F.at(fid);
        std::vector<size_t> a(f.Vids.begin(), f.Vids.end());
        std::vector<size_t> b(v.N_Vids.begin(), v.N_Vids.end());
        std::vector<size_t> comv = GetIntersection(a, b);
        if (comv.size() != 2) {
            std::cout << "nvids and face vids are not valid" << std::endl;
        }
    }

    std::cout << "checking validity of new faces " << std::endl;
    for (auto eid: v.N_Eids) {
        auto& e = mesh->E.at(eid);
        size_t nvid = e.Vids.at(0) != v.id ? e.Vids.at(0) : e.Vids.at(1);
        if (std::find(v.N_Vids.begin(), v.N_Vids.end(), nvid) == v.N_Vids.end()) {
            std::cout << "edge of center vertex is wrong" << std::endl;
        }
        for (auto fid: e.N_Fids) {
            auto& f = mesh->F.at(fid);
            std::vector<size_t> a(f.Vids.begin(), f.Vids.end());
            std::vector<size_t> b(e.Vids.begin(), e.Vids.end());
            std::vector<size_t> fvs = GetDifference(a, b);
            if (fvs.size() != 2) {
                std::cout << "difference between face vids and edge vids is not 2" << std::endl;
            }
            std::vector<size_t> c(v.N_Vids.begin(), v.N_Vids.end());
            std::vector<size_t> comv = GetIntersection(fvs, c);
            if (comv.size() != 1) {
                std::cout << "cant find common v with face v" << std::endl;
            }
        }
    }

    std::cout << "checking validity of newv1 nvids and faces" << std::endl;
    std::cout << "newv1 nvids #: " << newV1.N_Vids.size() << std::endl;
    std::cout << "newv1 neids #: " << newV1.N_Eids.size() << std::endl;
    std::cout << "newv1 nfids #: " << newV1.N_Fids.size() << std::endl;
    for (auto fid: newV1.N_Fids) {
        auto& f = mesh->F.at(fid);
        std::vector<size_t> a(f.Vids.begin(), f.Vids.end());
        std::vector<size_t> b(newV1.N_Vids.begin(), newV1.N_Vids.end());
        std::vector<size_t> comv = GetIntersection(a, b);
        if (comv.size() != 2) {
            std::cout << "nvids and face vids are not valid" << std::endl;
        }
    }

    std::cout << "checking validity of new faces " << std::endl;
    for (auto eid: newV1.N_Eids) {
        auto& e = mesh->E.at(eid);
        size_t nvid = e.Vids.at(0) != newV1.id ? e.Vids.at(0) : e.Vids.at(1);
        if (std::find(newV1.N_Vids.begin(), newV1.N_Vids.end(), nvid) == newV1.N_Vids.end()) {
            std::cout << "edge of center vertex is wrong" << std::endl;
        }
        for (auto fid: e.N_Fids) {
            auto& f = mesh->F.at(fid);
            std::vector<size_t> a(f.Vids.begin(), f.Vids.end());
            std::vector<size_t> b(e.Vids.begin(), e.Vids.end());
            std::vector<size_t> fvs = GetDifference(a, b);
            if (fvs.size() != 2) {
                std::cout << "difference between face vids and edge vids is not 2" << std::endl;
            }
            std::vector<size_t> c(newV1.N_Vids.begin(), newV1.N_Vids.end());
            std::vector<size_t> comv = GetIntersection(fvs, c);
            if (comv.size() != 1) {
                std::cout << "cant find common v with face v" << std::endl;
            }
        }
    }

    std::cout << "checking validity of newv2 nvids and faces" << std::endl;
    std::cout << "newv2 nvids #: " << newV2.N_Vids.size() << std::endl;
    std::cout << "newv2 neids #: " << newV2.N_Eids.size() << std::endl;
    std::cout << "newv2 nfids #: " << newV2.N_Fids.size() << std::endl;
    for (auto fid: newV2.N_Fids) {
        auto& f = mesh->F.at(fid);
        std::vector<size_t> a(f.Vids.begin(), f.Vids.end());
        std::vector<size_t> b(newV2.N_Vids.begin(), newV2.N_Vids.end());
        std::vector<size_t> comv = GetIntersection(a, b);
        if (comv.size() != 2) {
            std::cout << "nvids and face vids are not valid" << std::endl;
        }
    }

    std::cout << "checking validity of new faces " << std::endl;
    for (auto eid: newV2.N_Eids) {
        auto& e = mesh->E.at(eid);
        size_t nvid = e.Vids.at(0) != newV2.id ? e.Vids.at(0) : e.Vids.at(1);
        if (std::find(newV2.N_Vids.begin(), newV2.N_Vids.end(), nvid) == newV2.N_Vids.end()) {
            std::cout << "edge of center vertex is wrong" << std::endl;
        }
        for (auto fid: e.N_Fids) {
            auto& f = mesh->F.at(fid);
            std::vector<size_t> a(f.Vids.begin(), f.Vids.end());
            std::vector<size_t> b(e.Vids.begin(), e.Vids.end());
            std::vector<size_t> fvs = GetDifference(a, b);
            if (fvs.size() != 2) {
                std::cout << "difference between face vids and edge vids is not 2" << std::endl;
            }
            std::vector<size_t> c(newV2.N_Vids.begin(), newV2.N_Vids.end());
            std::vector<size_t> comv = GetIntersection(fvs, c);
            if (comv.size() != 1) {
                std::cout << "cant find common v with face v" << std::endl;
            }
        }
    }*/

    /*std::cout << "checking edges of center vertex: " << v.id << std::endl;
    for (auto eid: v.N_Eids) {
        auto& e = mesh->E.at(eid);
        if (e.N_Fids.size() != 2) {
            std::cout << "center vertex has wrong edge " << e.N_Fids.size() << std::endl;
        }
    }

    std::cout << "checking edges of new vertex 1: " << newV1.id << std::endl;
    for (auto eid: newV1.N_Eids) {
        auto& e = mesh->E.at(eid);
        if (e.N_Fids.size() != 2) {
            std::cout << "center vertex has wrong edge "<< e.N_Fids.size()  << std::endl;
        }
    }

    std::cout << "checking edges of new vertex 2: " << newV2.id << std::endl;
    for (auto eid: newV2.N_Eids) {
        auto& e = mesh->E.at(eid);
        if (e.N_Fids.size() != 2) {
            std::cout << "center vertex has wrong edge "<< e.N_Fids.size()  << std::endl;
        }
    }

    std::cout << "checking faces of center vertex: " << v.id << std::endl;
    for (auto fid: v.N_Fids) {
        auto& f = mesh->F.at(fid);
        std::cout << f.Vids.size() << " " << f.Eids.size() << std::endl;
        if (f.Vids.size() != 4) {
            std::cout << "center vertex has wrong face with vertices " << f.Vids.size() << std::endl;
        }
        if (f.Eids.size() != 4) {
            std::cout << "center vertex has wrong face with edges " << f.Eids.size() << std::endl;
        }
    }

    std::cout << "checking faces of new vertex 1: " << newV1.id << std::endl;
    for (auto fid: newV1.N_Fids) {
        auto& f = mesh->F.at(fid);
        std::cout << f.Vids.size() << " " << f.Eids.size() << std::endl;
        if (f.Vids.size() != 4) {
            std::cout << "center vertex has wrong face with vertices " << f.Vids.size() << std::endl;
        }
        if (f.Eids.size() != 4) {
            std::cout << "center vertex has wrong face with edges " << f.Eids.size() << std::endl;
        }
    }

    std::cout << "checking faces of new vertex 2: " << newV2.id << std::endl;
    for (auto fid: newV2.N_Fids) {
        auto& f = mesh->F.at(fid);
        std::cout << f.Vids.size() << " " << f.Eids.size() << std::endl;
        if (f.Vids.size() != 4) {
            std::cout << "center vertex has wrong face with vertices " << f.Vids.size() << std::endl;
        }
        if (f.Eids.size() != 4) {
            std::cout << "center vertex has wrong face with edges " << f.Eids.size() << std::endl;
        }
    }*/

    for (auto fid: facesToAvoid) {
        mesh->F.at(fid).N_Fids.clear();
    }
    Smooth();   
}

void VertexSplit::FiveVertexSplit() {
    if (splitEdges.empty()) {
       splitEdges = SelectEdgesToSplit();
    }
    Edge* e1 = &mesh->E.at(splitEdges.at(0));
    Edge* e2 = &mesh->E.at(splitEdges.at(1));
    
    Vertex new_V1(0.5 * (mesh->V.at(e1->Vids.at(0)).xyz() + mesh->V.at(e1->Vids.at(1)).xyz()));
    Vertex new_V2(0.5 * (mesh->V.at(e2->Vids.at(0)).xyz() + mesh->V.at(e2->Vids.at(1)).xyz()));
    new_V1.id = mesh->V.size();
    mesh->V.push_back(new_V1);
    new_V2.id = mesh->V.size();
    mesh->V.push_back(new_V2);

    Vertex* v = &mesh->V.at(vid);
    Vertex* newV1 = &mesh->V.at(new_V1.id);
    Vertex* newV2 = &mesh->V.at(new_V2.id);
    newV1->type = v->type;
    newV2->type = v->type;
    newV1->isBoundary = v->isBoundary;
    newV2->isBoundary = v->isBoundary;
    newV1->idealValence = v->idealValence;
    newV2->idealValence = v->idealValence;

    std::vector<size_t> temp_a = GetIntersection(v->N_Eids, mesh->F.at(e1->N_Fids.at(0)).Eids);
    AddContents(temp_a, GetIntersection(v->N_Eids, mesh->F.at(e1->N_Fids.at(1)).Eids));
    
    std::vector<size_t> temp_b = GetIntersection(v->N_Eids, mesh->F.at(e2->N_Fids.at(0)).Eids);
    AddContents(temp_b, GetIntersection(v->N_Eids, mesh->F.at(e2->N_Fids.at(1)).Eids));

    std::vector<size_t> diffFaces = e1->N_Fids;
    AddContents(diffFaces, e2->N_Fids);

    size_t edgeToRemoveId = GetIntersection(temp_a, temp_b).at(0);
    size_t faceToRemoveId = GetDifference(v->N_Fids, diffFaces).at(0);
    Edge* edgeToRemove = &mesh->E.at(edgeToRemoveId);
    Face* faceToRemove = &mesh->F.at(faceToRemoveId);
    
    // std::cout << "faceToRemove: " << faceToRemove->id << std::endl;
    
    // std::cout << "faceToRemove Vids: ";
    // for (auto id: faceToRemove->Vids) std::cout << id << " ";
    // std::cout << std::endl;

    // std::cout << "faceToRemove Eids: ";
    // for (auto id: faceToRemove->Eids) std::cout << id << " ";
    // std::cout << std::endl;

    Edge newEdge_1;
    newEdge_1.Vids = std::vector<size_t>{newV1->id, GetDifference(edgeToRemove->Vids, std::vector<size_t>{v->id}).at(0)};
    Edge newEdge_2;
    newEdge_2.Vids = std::vector<size_t>{newV2->id, GetDifference(edgeToRemove->Vids, std::vector<size_t>{v->id}).at(0)};
    Edge newEdge_3;
    newEdge_3.Vids = std::vector<size_t>{faceToRemove->Vids.at((std::distance(faceToRemove->Vids.begin(), std::find(faceToRemove->Vids.begin(), faceToRemove->Vids.end(), v->id))+2)%faceToRemove->Vids.size()), v->id};
    Edge newEdge_4;
    newEdge_4.Vids = std::vector<size_t>{v->id, newV1->id};
    Edge newEdge_5;
    newEdge_5.Vids = std::vector<size_t>{v->id, newV2->id};

    newEdge_1.id = mesh->E.size();
    mesh->E.push_back(newEdge_1);
    newEdge_2.id = mesh->E.size();
    mesh->E.push_back(newEdge_2);
    newEdge_3.id = mesh->E.size();
    mesh->E.push_back(newEdge_3);
    newEdge_4.id = mesh->E.size();
    mesh->E.push_back(newEdge_4);
    newEdge_5.id = mesh->E.size();
    mesh->E.push_back(newEdge_5);

    edgeToRemove = &mesh->E.at(edgeToRemoveId);
    e1 = &mesh->E.at(splitEdges.at(0));
    e2 = &mesh->E.at(splitEdges.at(1));

    Edge* newEdge1 = &mesh->E.at(newEdge_1.id);
    Edge* newEdge2 = &mesh->E.at(newEdge_2.id);
    Edge* newEdge3 = &mesh->E.at(newEdge_3.id);
    Edge* newEdge4 = &mesh->E.at(newEdge_4.id);
    Edge* newEdge5 = &mesh->E.at(newEdge_5.id);

    Face newF_1;
    size_t evid = edgeToRemove->Vids.at(0) == v->id ? edgeToRemove->Vids.at(1) : edgeToRemove->Vids.at(0);
    newF_1.Vids.push_back(evid);
    for (auto fid: edgeToRemove->N_Fids) {
        Face* f = &mesh->F.at(fid);
        int idx = std::distance(f->Vids.begin(), std::find(f->Vids.begin(), f->Vids.end(), v->id));
        if (f->Vids.at((idx+1)%f->Vids.size()) == evid && std::find(f->Eids.begin(), f->Eids.end(), e1->id) != f->Eids.end()) {
            newF_1.Vids.push_back(newV1->id);
            newF_1.Vids.push_back(v->id);
            newF_1.Vids.push_back(newV2->id);
            break;
        } else if (f->Vids.at((idx+1)%f->Vids.size()) == evid && std::find(f->Eids.begin(), f->Eids.end(), e2->id) != f->Eids.end()) {
            newF_1.Vids.push_back(newV2->id);
            newF_1.Vids.push_back(v->id);
            newF_1.Vids.push_back(newV1->id);
            break;
        }    
    }

    Face newF_2;
    Face newF_3;
    int idx = std::distance(newF_1.Vids.begin(), std::find(newF_1.Vids.begin(), newF_1.Vids.end(), v->id));
    int idx2 = std::distance(faceToRemove->Vids.begin(), std::find(faceToRemove->Vids.begin(), faceToRemove->Vids.end(), v->id));
    if (newF_1.Vids.at((idx+1)%newF_1.Vids.size()) == newV1->id) {
        newF_2.Vids = {newV1->id, v->id, faceToRemove->Vids.at((idx2+2)%faceToRemove->Vids.size()), faceToRemove->Vids.at((idx2+3)%faceToRemove->Vids.size())};
        newF_3.Vids = {v->id, newV2->id, faceToRemove->Vids.at((idx2+1)%faceToRemove->Vids.size()), faceToRemove->Vids.at((idx2+2)%faceToRemove->Vids.size())};
    } else if (newF_1.Vids.at((idx+1)%newF_1.Vids.size()) == newV2->id) {
        newF_2.Vids = {newV2->id, v->id, faceToRemove->Vids.at((idx2+2)%faceToRemove->Vids.size()), faceToRemove->Vids.at((idx2+3)%faceToRemove->Vids.size())};
        newF_3.Vids = {v->id, newV1->id, faceToRemove->Vids.at((idx2+1)%faceToRemove->Vids.size()), faceToRemove->Vids.at((idx2+2)%faceToRemove->Vids.size())};
    }


    newF_1.id = mesh->F.size();
    mesh->F.push_back(newF_1);
    newF_2.id = mesh->F.size();
    mesh->F.push_back(newF_2);
    newF_3.id = mesh->F.size();
    mesh->F.push_back(newF_3);
    mu->delta += 3;

    faceToRemove = &mesh->F.at(faceToRemoveId);
    // std::cout << "faceToRemove: " << faceToRemove->id << std::endl;
    
    // std::cout << "faceToRemove Vids: ";
    // for (auto id: faceToRemove->Vids) std::cout << id << " ";
    // std::cout << std::endl;

    // std::cout << "faceToRemove Eids: ";
    // for (auto id: faceToRemove->Eids) std::cout << id << " ";
    // std::cout << std::endl;
    

    Face* newF1 = &mesh->F.at(newF_1.id);
    Face* newF2 = &mesh->F.at(newF_2.id);
    Face* newF3 = &mesh->F.at(newF_3.id);
    // std::cout << "newV1: " << newV1->id << " newV2: " << newV2->id << std::endl; 
    // std::cout << "edgeToRemove: " << edgeToRemove->id << " " << edgeToRemove->Vids.at(0) << " " << edgeToRemove->Vids.at(1) << std::endl;
    // std::cout << "e1: " << e1->id << " " << e1->Vids.at(0) << " " << e1->Vids.at(1) << std::endl;
    // std::cout << "e2: " << e2->id << " " << e2->Vids.at(0) << " " << e2->Vids.at(1) << std::endl;
    // std::cout << "newEdge1: " << newEdge1->id << " " << newEdge1->Vids.at(0) << " " << newEdge1->Vids.at(1) << std::endl;
    // std::cout << "newEdge2: " << newEdge2->id << " " << newEdge2->Vids.at(0) << " " << newEdge2->Vids.at(1) << std::endl;
    // std::cout << "newEdge3: " << newEdge3->id << " " << newEdge3->Vids.at(0) << " " << newEdge3->Vids.at(1) << std::endl;
    // std::cout << "newEdge4: " << newEdge4->id << " " << newEdge4->Vids.at(0) << " " << newEdge4->Vids.at(1) << std::endl;
    // std::cout << "newEdge5: " << newEdge5->id << " " << newEdge5->Vids.at(0) << " " << newEdge5->Vids.at(1) << std::endl;

    // std::cout << "faceToRemove Vids: ";
    // for (auto id: faceToRemove->Vids) std::cout << id << " ";
    // std::cout << std::endl;

    // std::cout << "faceToRemove Eids: ";
    // for (auto id: faceToRemove->Eids) std::cout << id << " ";
    // std::cout << std::endl;

    // std::cout << "edgeToRemove NFids: " << edgeToRemove->N_Fids.at(0) << " " << edgeToRemove->N_Fids.at(1) << std::endl;
    // std::cout << "e1 NFids: " << e1->N_Fids.at(0) << " " << e1->N_Fids.at(1) << std::endl;
    // std::cout << "e2 NFids: " << e2->N_Fids.at(0) << " " << e2->N_Fids.at(1) << std::endl;
    // std::cout << "faceToRemove: " << faceToRemove->id << std::endl;
    // std::cout << "newF1 Vids: " << newF1->Vids.size() << " " << newF1->Vids.at(0) << " " << newF1->Vids.at(1) << " " << newF1->Vids.at(2) << " " << newF1->Vids.at(3) << " " << std::endl;
    // std::cout << "newF2: Vids: " << newF2->Vids.size() << " " << newF2->Vids.at(0) << " " << newF2->Vids.at(1) << " " << newF2->Vids.at(2) << " " << newF2->Vids.at(3) << " " <<  std::endl;
    // std::cout << "newF3: Vids: " << newF3->Vids.size() << " " << newF3->Vids.at(0) << " " << newF3->Vids.at(1) << " " << newF3->Vids.at(2) << " " << newF3->Vids.at(3) << " " <<  std::endl;

    // e1->Vids.at(0) == v->id ? newV1->N_Vids.push_back(e1->Vids.at(1)) : newV1->N_Vids.push_back(e1->Vids.at(0));
    // newV1->N_Eids.push_back(e1->id);
    e1->Vids.at(0) == v->id ? mesh->V.at(e1->Vids.at(1)).N_Vids = GetDifference(mesh->V.at(e1->Vids.at(1)).N_Vids, std::vector<size_t>{v->id}) : mesh->V.at(e1->Vids.at(0)).N_Vids = GetDifference(mesh->V.at(e1->Vids.at(0)).N_Vids, std::vector<size_t>{v->id});
    e1->Vids.at(std::distance(e1->Vids.begin(), std::find(e1->Vids.begin(), e1->Vids.end(), v->id))) = newV1->id;
    
    // e2->Vids.at(0) == v->id ? newV2->N_Vids.push_back(e2->Vids.at(1)) : newV2->N_Vids.push_back(e2->Vids.at(0));
    // newV2->N_Eids.push_back(e2->id);
    e2->Vids.at(0) == v->id ? mesh->V.at(e2->Vids.at(1)).N_Vids = GetDifference(mesh->V.at(e2->Vids.at(1)).N_Vids, std::vector<size_t>{v->id}) : mesh->V.at(e2->Vids.at(0)).N_Vids = GetDifference(mesh->V.at(e2->Vids.at(0)).N_Vids, std::vector<size_t>{v->id});
    e2->Vids.at(std::distance(e2->Vids.begin(), std::find(e2->Vids.begin(), e2->Vids.end(), v->id))) = newV2->id;
    for (auto fid: e1->N_Fids) {
        auto& f = mesh->F.at(fid);
        for (auto eid: f.Eids) {
            auto& e = mesh->E.at(eid);
            if (eid == edgeToRemove->id || eid == e1->id || (e.Vids.at(0) != v->id && e.Vids.at(1) != v->id)) continue;
            e.Vids.at(0) == v->id ? newV1->N_Vids.push_back(e.Vids.at(1)) : newV1->N_Vids.push_back(e.Vids.at(0));
            e.Vids.at(0) == v->id ? mesh->V.at(e.Vids.at(1)).N_Vids.at(std::distance(mesh->V.at(e.Vids.at(1)).N_Vids.begin(), std::find(mesh->V.at(e.Vids.at(1)).N_Vids.begin(), mesh->V.at(e.Vids.at(1)).N_Vids.end(), v->id))) = newV1->id : mesh->V.at(e.Vids.at(0)).N_Vids.at(std::distance(mesh->V.at(e.Vids.at(0)).N_Vids.begin(), std::find(mesh->V.at(e.Vids.at(0)).N_Vids.begin(), mesh->V.at(e.Vids.at(0)).N_Vids.end(), v->id))) = newV1->id;
            newV1->N_Eids.push_back(e.id);
            e.Vids.at(std::distance(e.Vids.begin(), std::find(e.Vids.begin(), e.Vids.end(), v->id))) = newV1->id;
        }
    }
    for (auto fid: e2->N_Fids) {
        auto& f = mesh->F.at(fid);
        for (auto eid: f.Eids) {
            auto& e = mesh->E.at(eid);
            if (eid == edgeToRemove->id || eid == e2->id || (e.Vids.at(0) != v->id && e.Vids.at(1) != v->id)) continue;
            e.Vids.at(0) == v->id ? newV2->N_Vids.push_back(e.Vids.at(1)) : newV2->N_Vids.push_back(e.Vids.at(0));
            e.Vids.at(0) == v->id ? mesh->V.at(e.Vids.at(1)).N_Vids.at(std::distance(mesh->V.at(e.Vids.at(1)).N_Vids.begin(), std::find(mesh->V.at(e.Vids.at(1)).N_Vids.begin(), mesh->V.at(e.Vids.at(1)).N_Vids.end(), v->id))) = newV2->id : mesh->V.at(e.Vids.at(0)).N_Vids.at(std::distance(mesh->V.at(e.Vids.at(0)).N_Vids.begin(), std::find(mesh->V.at(e.Vids.at(0)).N_Vids.begin(), mesh->V.at(e.Vids.at(0)).N_Vids.end(), v->id))) = newV2->id;
            newV2->N_Eids.push_back(e.id);
            e.Vids.at(std::distance(e.Vids.begin(), std::find(e.Vids.begin(), e.Vids.end(), v->id))) = newV2->id;
        }
    }
    mesh->F.at(e1->N_Fids.at(0)).Vids.at(std::distance(mesh->F.at(e1->N_Fids.at(0)).Vids.begin(), std::find(mesh->F.at(e1->N_Fids.at(0)).Vids.begin(), mesh->F.at(e1->N_Fids.at(0)).Vids.end(), v->id))) = newV1->id;
    mesh->F.at(e1->N_Fids.at(1)).Vids.at(std::distance(mesh->F.at(e1->N_Fids.at(1)).Vids.begin(), std::find(mesh->F.at(e1->N_Fids.at(1)).Vids.begin(), mesh->F.at(e1->N_Fids.at(1)).Vids.end(), v->id))) = newV1->id;
    mesh->F.at(e2->N_Fids.at(0)).Vids.at(std::distance(mesh->F.at(e2->N_Fids.at(0)).Vids.begin(), std::find(mesh->F.at(e2->N_Fids.at(0)).Vids.begin(), mesh->F.at(e2->N_Fids.at(0)).Vids.end(), v->id))) = newV2->id;
    mesh->F.at(e2->N_Fids.at(1)).Vids.at(std::distance(mesh->F.at(e2->N_Fids.at(1)).Vids.begin(), std::find(mesh->F.at(e2->N_Fids.at(1)).Vids.begin(), mesh->F.at(e2->N_Fids.at(1)).Vids.end(), v->id))) = newV2->id;

    newEdge1->N_Fids.push_back(newF1->id);
    newEdge2->N_Fids.push_back(newF1->id);
    newEdge3->N_Fids = std::vector<size_t>{newF2->id, newF3->id};
    newEdge4->N_Fids.push_back(newF1->id);
    newEdge5->N_Fids.push_back(newF1->id);
    for (auto fid: edgeToRemove->N_Fids) {
        auto& f = mesh->F.at(fid);
        if (std::find(f.Vids.begin(), f.Vids.end(), newV1->id) != f.Vids.end()) {
            f.Eids.at(std::distance(f.Eids.begin(), std::find(f.Eids.begin(), f.Eids.end(), edgeToRemove->id))) = newEdge1->id;
            newEdge1->N_Fids.push_back(f.id); 
        }
        if (std::find(f.Vids.begin(), f.Vids.end(), newV2->id) != f.Vids.end()) {
            f.Eids.at(std::distance(f.Eids.begin(), std::find(f.Eids.begin(), f.Eids.end(), edgeToRemove->id))) = newEdge2->id;
            newEdge2->N_Fids.push_back(f.id);
        }
    }
    AddContents(newF1->Eids, std::vector<size_t>{newEdge1->id, newEdge2->id, newEdge4->id, newEdge5->id});
    newF2->Eids.push_back(newEdge3->id);
    newF3->Eids.push_back(newEdge3->id);
    if (std::find(newF2->Vids.begin(), newF2->Vids.end(), new_V1.id) != newF2->Vids.end()) {
        newF2->Eids.push_back(newEdge4->id);
        newEdge4->N_Fids.push_back(newF2->id);
        newF3->Eids.push_back(newEdge5->id);
        newEdge5->N_Fids.push_back(newF3->id);
    } else {
        newF2->Eids.push_back(newEdge5->id);
        newEdge5->N_Fids.push_back(newF2->id);
        newF3->Eids.push_back(newEdge4->id);
        newEdge4->N_Fids.push_back(newF3->id);
    }
    for (auto eid: faceToRemove->Eids) {
        auto& e = mesh->E.at(eid);
        if (std::find(newF2->Vids.begin(), newF2->Vids.end(), e.Vids.at(0)) != newF2->Vids.end() && std::find(newF2->Vids.begin(), newF2->Vids.end(), e.Vids.at(1)) != newF2->Vids.end()) {
            newF2->Eids.push_back(e.id);
            e.N_Fids.at(0) == faceToRemove->id ? e.N_Fids.at(0) = newF2->id : e.N_Fids.at(1) = newF2->id;
        } else if (std::find(newF3->Vids.begin(), newF3->Vids.end(), e.Vids.at(0)) != newF3->Vids.end() && std::find(newF3->Vids.begin(), newF3->Vids.end(), e.Vids.at(1)) != newF3->Vids.end()) {
            newF3->Eids.push_back(e.id);
            e.N_Fids.at(0) == faceToRemove->id ? e.N_Fids.at(0) = newF3->id : e.N_Fids.at(1) = newF3->id;
        }
    }
    std::vector<size_t> v_nvids = v->N_Vids;
    std::vector<size_t> v_neids = v->N_Eids;
    std::vector<size_t> v_nfids = v->N_Fids;
    for (auto eid: std::vector<size_t>{e1->id, e2->id, newEdge1->id, newEdge2->id, newEdge3->id, newEdge4->id, newEdge5->id}) {
        auto& e = mesh->E.at(eid);
        AddContents(mesh->V.at(e.Vids.at(0)).N_Vids, std::vector<size_t>{e.Vids.at(1)});
        AddContents(mesh->V.at(e.Vids.at(1)).N_Vids, std::vector<size_t>{e.Vids.at(0)});
        AddContents(mesh->V.at(e.Vids.at(0)).N_Eids, std::vector<size_t>{e.id});
        AddContents(mesh->V.at(e.Vids.at(1)).N_Eids, std::vector<size_t>{e.id});
        AddContents(mesh->V.at(e.Vids.at(0)).N_Fids, e.N_Fids);
        AddContents(mesh->V.at(e.Vids.at(1)).N_Fids, e.N_Fids);
    }
    if (edgeToRemove->Vids.at(0) == v->id) {
        mesh->V.at(edgeToRemove->Vids.at(1)).N_Eids = GetDifference(mesh->V.at(edgeToRemove->Vids.at(1)).N_Eids, std::vector<size_t>{edgeToRemove->id});
        mesh->V.at(edgeToRemove->Vids.at(1)).N_Vids = GetDifference(mesh->V.at(edgeToRemove->Vids.at(1)).N_Vids, std::vector<size_t>{v->id});
    } else {
        mesh->V.at(edgeToRemove->Vids.at(0)).N_Eids = GetDifference(mesh->V.at(edgeToRemove->Vids.at(0)).N_Eids, std::vector<size_t>{edgeToRemove->id});
        mesh->V.at(edgeToRemove->Vids.at(0)).N_Vids = GetDifference(mesh->V.at(edgeToRemove->Vids.at(0)).N_Vids, std::vector<size_t>{v->id});
    }
    v->N_Vids = GetDifference(v->N_Vids, v_nvids);
    v->N_Eids = GetDifference(v->N_Eids, v_neids);
    v->N_Fids = GetDifference(v->N_Fids, v_nfids);
    for (auto fvid: faceToRemove->Vids) {
        auto& fv = mesh->V.at(fvid);
        fv.N_Fids = GetDifference(fv.N_Fids, std::vector<size_t>{faceToRemove->id});
        
    }
    for (auto nfid: faceToRemove->N_Fids) {
        mesh->F.at(nfid).N_Fids = GetDifference(mesh->F.at(nfid).N_Fids, std::vector<size_t>{faceToRemove->id});
    }
    for (auto fid: std::vector<size_t>{newF1->id, newF2->id, newF3->id}) {
        auto& newF = mesh->F.at(fid);
        for (auto fvid: newF.Vids) {
            AddContents(ToSmooth, std::vector<size_t>{fvid});
            AddContents(ToSmooth, mesh->V.at(fvid).N_Vids);
            AddContents(newF.N_Fids, mesh->V.at(fvid).N_Fids);
            AddContents(mesh->V.at(fvid).N_Fids, std::vector<size_t>{newF.id});
        }
    }

    // std::cout << "newF1 Vids: ";
    // for (auto id: newF1->Vids) std::cout << id << " ";
    // std::cout << std::endl;

    // std::cout << "newF1 Eids: ";
    // for (auto id: newF1->Eids) std::cout << id << " ";
    // std::cout << std::endl;

    // std::cout << "newF2 Vids: ";
    // for (auto id: newF2->Vids) std::cout << id << " ";
    // std::cout << std::endl;

    // std::cout << "newF2 Eids: ";
    // for (auto id: newF2->Eids) std::cout << id << " ";
    // std::cout << std::endl;

    // std::cout << "newF3 Vids: ";
    // for (auto id: newF3->Vids) std::cout << id << " ";
    // std::cout << std::endl;

    // std::cout << "newF3 Eids: ";
    // for (auto id: newF3->Eids) std::cout << id << " ";
    // std::cout << std::endl;

    // std::cout << "faceToRemove: " << faceToRemove->id << std::endl;
    
    // std::cout << "faceToRemove Vids: ";
    // for (auto id: faceToRemove->Vids) std::cout << id << " ";
    // std::cout << std::endl;

    // std::cout << "faceToRemove Eids: ";
    // for (auto id: faceToRemove->Eids) std::cout << id << " ";
    // std::cout << std::endl;

    edgeToRemove->Vids.clear();
    edgeToRemove->N_Fids.clear();
    faceToRemove->Vids.clear();
    faceToRemove->Eids.clear();
    faceToRemove->N_Fids.clear();

    // std::cout << "faceToRemove: " << faceToRemove->id << std::endl;
    
    // std::cout << "faceToRemove Vids: ";
    // for (auto id: faceToRemove->Vids) std::cout << id << " ";
    // std::cout << std::endl;

    // std::cout << "faceToRemove Eids: ";
    // for (auto id: faceToRemove->Eids) std::cout << id << " ";
    // std::cout << std::endl;
    // std::cout << "Before smoothing" << std::endl;
    Smooth();
    // std::cout << "After smoothing" << std::endl;
}

std::vector<size_t> VertexSplit::SelectEdgesToSplit() {
    std::vector<size_t> splitEdges;
    std::vector<size_t> EdgesToAvoid;

    auto& v = mesh->V.at(vid);
    if (v.type == FEATURE) {
        for (auto eid: v.N_Eids) {
            auto& e = mesh->E.at(eid);
            if (mesh->V.at(e.Vids.at(0)).type == FEATURE && mesh->V.at(e.Vids.at(1)).type == FEATURE) AddContents(splitEdges, std::vector<size_t>{e.id});
        }
        return splitEdges;
    }
    if (v.isBoundary) {
        for (auto eid: v.N_Eids) {
            auto& e = mesh->E.at(eid);
            if (mesh->V.at(e.Vids.at(0)).isBoundary && mesh->V.at(e.Vids.at(1)).isBoundary) AddContents(splitEdges, std::vector<size_t>{e.id});
        }
        return splitEdges;
    }
    auto& e = mesh->E.at(v.N_Eids.at(0));
    AddContents(splitEdges, std::vector<size_t>{e.id});
    for (auto fid: e.N_Fids) {
        AddContents(EdgesToAvoid, GetIntersection(v.N_Eids, mesh->F.at(fid).Eids));
    }
    
    std::vector<size_t> RemainingEdges = GetDifference(v.N_Eids, EdgesToAvoid);
    AddContents(splitEdges, std::vector<size_t>{RemainingEdges.at(0)});
    
    return splitEdges;
}