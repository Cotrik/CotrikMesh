#include "VertexSplit.h"

VertexSplit::VertexSplit() {}
VertexSplit::VertexSplit(Mesh& mesh_, MeshUtil& mu_, Smoother& smoother_, size_t vid_, std::vector<size_t> edgesToSplit_) : SimplificationOperation(mesh_, mu_, smoother_) {
    vid = vid_;
    splitEdges = edgesToSplit_;
}
VertexSplit::~VertexSplit() {}

bool VertexSplit::IsOperationValid() {
    CheckValidity();

    auto& v = mesh.V.at(vid);
    int count = 0;
    if (v.type == FEATURE) {
        for (auto nvid: v.N_Vids) {
            if (mesh.V.at(nvid).type == FEATURE) count += 1;
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


    if (splitEdges.empty()) {
       splitEdges = SelectEdgesToSplit();
    }
    // std::cout << "Chose edges: " << splitEdges.size() << std::endl;

    Edge& e1 = mesh.E.at(splitEdges.at(0));
    Edge& e2 = mesh.E.at(splitEdges.at(1));
    
    Vertex new_V1(0.5 * (mesh.V.at(e1.Vids.at(0)).xyz() + mesh.V.at(e1.Vids.at(1)).xyz()));
    Vertex new_V2(0.5 * (mesh.V.at(e2.Vids.at(0)).xyz() + mesh.V.at(e2.Vids.at(1)).xyz()));
    new_V1.id = mesh.V.size();
    mesh.V.push_back(new_V1);
    new_V2.id = mesh.V.size();
    mesh.V.push_back(new_V2);

    auto& v = mesh.V.at(vid);
    auto& newV1 = mesh.V.at(new_V1.id);
    auto& newV2 = mesh.V.at(new_V2.id);
    newV1.type = v.type;
    newV2.type = v.type;

    std::vector<size_t> facesToInspect;
    AddContents(facesToInspect, v.N_Fids);
    
    // std::cout << v.id << std::endl;

    // std::cout << "Added 2 new vertices" << std::endl;

    std::vector<size_t> edgesToChange = {};
    std::vector<size_t> facesToRemove = GetDifference(v.N_Fids, e1.N_Fids);
    for (auto fid: e1.N_Fids) {
        auto& f = mesh.F.at(fid);
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
    // std::cout << "edges to avoid: " << edgesToChange.size() << std::endl;
    // for (auto eid: edgesToChange) {
    //     std::cout << eid << std::endl;
    // }
    // std::cout << "edges to avoid: " << edgesToAvoid.size() << std::endl;
    // for (auto eid: edgesToAvoid) {
    //     std::cout << eid << std::endl;
    // }
    std::vector<size_t> facesToAvoid = {};
    for (auto eid: edgesToChange) {
        auto& e = mesh.E.at(eid);
        std::vector<size_t> edgeV;
        e.Vids.at(0) == vid ? edgeV.push_back(e.Vids.at(1)) : edgeV.push_back(e.Vids.at(0));
        e.Vids.at(0) == vid ? e.Vids.at(0) = newV1.id : e.Vids.at(1) = newV1.id;
        
        AddContents(facesToAvoid, GetDifference(e.N_Fids, e1.N_Fids));
        AddContents(newV1.N_Eids, std::vector<size_t>{eid});
        AddContents(newV1.N_Vids, edgeV);
        UpdateContents(mesh.V.at(edgeV.at(0)).N_Vids, std::vector<size_t>{vid});
        AddContents(mesh.V.at(edgeV.at(0)).N_Vids, std::vector<size_t>{newV1.id});
    }
    // std::cout << "faces to avoid: " << facesToAvoid.size() << std::endl;
    // return;

    std::vector<size_t> facesToChange = {};
    for (auto eid: edgesToAvoid) {
        auto& e = mesh.E.at(eid);
        std::vector<size_t> edgeV;
        e.Vids.at(0) == vid ? edgeV.push_back(e.Vids.at(1)) : edgeV.push_back(e.Vids.at(0));
        e.Vids.at(0) == vid ? e.Vids.at(0) = newV2.id : e.Vids.at(1) = newV2.id;
        AddContents(facesToChange, GetDifference(e.N_Fids, facesToAvoid));
        AddContents(newV2.N_Eids, std::vector<size_t>{eid});
        AddContents(newV2.N_Vids, edgeV);
        UpdateContents(mesh.V.at(edgeV.at(0)).N_Vids, std::vector<size_t>{vid});
        AddContents(mesh.V.at(edgeV.at(0)).N_Vids, std::vector<size_t>{newV2.id});
    }

    facesToRemove = GetDifference(v.N_Fids, facesToChange);
    for (auto fid: facesToChange) {
        auto& f = mesh.F.at(fid);
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
    for (auto fid: facesToAvoid) {
        auto& f = mesh.F.at(fid);
        for (auto eid: f.Eids) {
            auto& e = mesh.E.at(eid);
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
        auto& f = mesh.F.at(fid);
        int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), vid));
        newEdges.emplace_back(std::vector<size_t>{vid, f.Vids.at((idx+2)%f.Vids.size())});
        for (auto fvid: f.Vids) {
            UpdateContents(mesh.V.at(fvid).N_Fids, std::vector<size_t>{f.id});
        }
    }

    v.N_Vids.clear();
    v.N_Eids.clear();
    v.N_Fids.clear();

    std::vector<size_t> edgesToUpdate;
    for (auto& e: newEdges) {
        e.id = mesh.E.size();
        mesh.E.push_back(e);
        edgesToUpdate.push_back(e.id);
    }

    std::vector<size_t> facesToUpdate;
    for (auto& f: newFaces) {
        f.id = mesh.F.size();
        mesh.F.push_back(f);
        facesToUpdate.push_back(f.id);
    }

    for (auto eid: edgesToUpdate) {
        auto& e = mesh.E.at(eid);
        AddContents(mesh.V.at(e.Vids.at(0)).N_Vids, std::vector<size_t>{e.Vids.at(1)});
        AddContents(mesh.V.at(e.Vids.at(1)).N_Vids, std::vector<size_t>{e.Vids.at(0)});
        AddContents(mesh.V.at(e.Vids.at(0)).N_Eids, std::vector<size_t>{e.id});
        AddContents(mesh.V.at(e.Vids.at(1)).N_Eids, std::vector<size_t>{e.id});
        AddContents(v.N_Eids, std::vector<size_t>{e.id});
    }

    for (auto fid: facesToUpdate) {
        auto& f = mesh.F.at(fid);
        for (auto fid: facesToAvoid) {
            auto& nf = mesh.F.at(fid);
            for (auto eid: nf.Eids) {
                auto& e = mesh.E.at(eid);
                if (std::find(f.Vids.begin(), f.Vids.end(), e.Vids.at(0)) != f.Vids.end() && 
                    std::find(f.Vids.begin(), f.Vids.end(), e.Vids.at(1)) != f.Vids.end()) {
                        AddContents(f.Eids, std::vector<size_t>{e.id});
                        AddContents(e.N_Fids, std::vector<size_t>{f.id});
                }
            }
        }
        for (auto eid: edgesToUpdate) {
        // for (auto& e: newEdges) {
            auto& e = mesh.E.at(eid);
            if (std::find(f.Vids.begin(), f.Vids.end(), e.Vids.at(0)) != f.Vids.end() && 
                std::find(f.Vids.begin(), f.Vids.end(), e.Vids.at(1)) != f.Vids.end()) {
                    AddContents(f.Eids, std::vector<size_t>{e.id});
                    AddContents(e.N_Fids, std::vector<size_t>{f.id});
            }
        }
        for (auto fvid: f.Vids) {
            auto& fv = mesh.V.at(fvid);
            UpdateContents(fv.N_Fids, facesToAvoid);
            AddContents(fv.N_Fids, std::vector<size_t>{f.id});
        }
        AddContents(v.N_Fids, std::vector<size_t>{f.id});
        for (auto fvid: f.Vids) {
            AddContents(f.N_Fids, GetDifference(mesh.F.at(fvid).N_Fids, std::vector<size_t>{f.id}));
        }
    }

    AddContents(facesToInspect, facesToUpdate);
    UpdateContents(facesToInspect, facesToRemove);
    std::vector<size_t> verticesToCheck;
    for (auto fid: facesToInspect) {
        auto& f = mesh.F.at(fid);
        AddContents(verticesToCheck, std::vector<size_t>(f.Vids.begin(), f.Vids.end()));
    }
    for (auto id: verticesToCheck) {
        SetSingularity(id);
    }
    AddContents(ToSmooth, verticesToCheck);
    /*std::cout << "checking validity of nvids and faces" << std::endl;
    std::cout << "v nvids #: " << v.N_Vids.size() << std::endl;
    std::cout << "v neids #: " << v.N_Eids.size() << std::endl;
    std::cout << "v nfids #: " << v.N_Fids.size() << std::endl;
    for (auto fid: v.N_Fids) {
        auto& f = mesh.F.at(fid);
        std::vector<size_t> a(f.Vids.begin(), f.Vids.end());
        std::vector<size_t> b(v.N_Vids.begin(), v.N_Vids.end());
        std::vector<size_t> comv = GetIntersection(a, b);
        if (comv.size() != 2) {
            std::cout << "nvids and face vids are not valid" << std::endl;
        }
    }

    std::cout << "checking validity of new faces " << std::endl;
    for (auto eid: v.N_Eids) {
        auto& e = mesh.E.at(eid);
        size_t nvid = e.Vids.at(0) != v.id ? e.Vids.at(0) : e.Vids.at(1);
        if (std::find(v.N_Vids.begin(), v.N_Vids.end(), nvid) == v.N_Vids.end()) {
            std::cout << "edge of center vertex is wrong" << std::endl;
        }
        for (auto fid: e.N_Fids) {
            auto& f = mesh.F.at(fid);
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
        auto& f = mesh.F.at(fid);
        std::vector<size_t> a(f.Vids.begin(), f.Vids.end());
        std::vector<size_t> b(newV1.N_Vids.begin(), newV1.N_Vids.end());
        std::vector<size_t> comv = GetIntersection(a, b);
        if (comv.size() != 2) {
            std::cout << "nvids and face vids are not valid" << std::endl;
        }
    }

    std::cout << "checking validity of new faces " << std::endl;
    for (auto eid: newV1.N_Eids) {
        auto& e = mesh.E.at(eid);
        size_t nvid = e.Vids.at(0) != newV1.id ? e.Vids.at(0) : e.Vids.at(1);
        if (std::find(newV1.N_Vids.begin(), newV1.N_Vids.end(), nvid) == newV1.N_Vids.end()) {
            std::cout << "edge of center vertex is wrong" << std::endl;
        }
        for (auto fid: e.N_Fids) {
            auto& f = mesh.F.at(fid);
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
        auto& f = mesh.F.at(fid);
        std::vector<size_t> a(f.Vids.begin(), f.Vids.end());
        std::vector<size_t> b(newV2.N_Vids.begin(), newV2.N_Vids.end());
        std::vector<size_t> comv = GetIntersection(a, b);
        if (comv.size() != 2) {
            std::cout << "nvids and face vids are not valid" << std::endl;
        }
    }

    std::cout << "checking validity of new faces " << std::endl;
    for (auto eid: newV2.N_Eids) {
        auto& e = mesh.E.at(eid);
        size_t nvid = e.Vids.at(0) != newV2.id ? e.Vids.at(0) : e.Vids.at(1);
        if (std::find(newV2.N_Vids.begin(), newV2.N_Vids.end(), nvid) == newV2.N_Vids.end()) {
            std::cout << "edge of center vertex is wrong" << std::endl;
        }
        for (auto fid: e.N_Fids) {
            auto& f = mesh.F.at(fid);
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
        auto& e = mesh.E.at(eid);
        if (e.N_Fids.size() != 2) {
            std::cout << "center vertex has wrong edge " << e.N_Fids.size() << std::endl;
        }
    }

    std::cout << "checking edges of new vertex 1: " << newV1.id << std::endl;
    for (auto eid: newV1.N_Eids) {
        auto& e = mesh.E.at(eid);
        if (e.N_Fids.size() != 2) {
            std::cout << "center vertex has wrong edge "<< e.N_Fids.size()  << std::endl;
        }
    }

    std::cout << "checking edges of new vertex 2: " << newV2.id << std::endl;
    for (auto eid: newV2.N_Eids) {
        auto& e = mesh.E.at(eid);
        if (e.N_Fids.size() != 2) {
            std::cout << "center vertex has wrong edge "<< e.N_Fids.size()  << std::endl;
        }
    }

    std::cout << "checking faces of center vertex: " << v.id << std::endl;
    for (auto fid: v.N_Fids) {
        auto& f = mesh.F.at(fid);
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
        auto& f = mesh.F.at(fid);
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
        auto& f = mesh.F.at(fid);
        std::cout << f.Vids.size() << " " << f.Eids.size() << std::endl;
        if (f.Vids.size() != 4) {
            std::cout << "center vertex has wrong face with vertices " << f.Vids.size() << std::endl;
        }
        if (f.Eids.size() != 4) {
            std::cout << "center vertex has wrong face with edges " << f.Eids.size() << std::endl;
        }
    }*/

    for (auto fid: facesToAvoid) {
        mesh.F.at(fid).N_Fids.clear();
    }
    Smooth();   
}

std::vector<size_t> VertexSplit::SelectEdgesToSplit() {
    std::vector<size_t> splitEdges;
    std::vector<size_t> EdgesToAvoid;

    auto& v = mesh.V.at(vid);
    if (v.type == FEATURE) {
        for (auto eid: v.N_Eids) {
            auto& e = mesh.E.at(eid);
            if (mesh.V.at(e.Vids.at(0)).type == FEATURE && mesh.V.at(e.Vids.at(1)).type == FEATURE) AddContents(splitEdges, std::vector<size_t>{e.id});
        }
        return splitEdges;
    }
    auto& e = mesh.E.at(v.N_Eids.at(0));
    AddContents(splitEdges, std::vector<size_t>{e.id});
    for (auto fid: e.N_Fids) {
        AddContents(EdgesToAvoid, GetIntersection(v.N_Eids, mesh.F.at(fid).Eids));
    }
    
    std::vector<size_t> RemainingEdges = GetDifference(v.N_Eids, EdgesToAvoid);
    AddContents(splitEdges, std::vector<size_t>{RemainingEdges.at(0)});
    
    return splitEdges;
}