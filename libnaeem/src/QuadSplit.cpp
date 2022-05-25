#include "QuadSplit.h"

QuadSplit::QuadSplit() {}
QuadSplit::QuadSplit(Mesh& mesh_, MeshUtil& mu_, Smoother& smoother_, size_t vid_, std::vector<size_t> verticesToSplit_, std::vector<size_t> verticesToChange_) : SimplificationOperation(mesh_, mu_, smoother_) {
    vid = vid_;
    verticesToSplit = verticesToSplit_;
    verticesToChange = verticesToChange_;
}
QuadSplit::~QuadSplit() {}

void QuadSplit::SetRanking(glm::dvec3 d) {
    CheckValidity();
}

bool QuadSplit::IsOperationValid() {
    CheckValidity();

    return true;
}

void QuadSplit::PerformOperation() {
    CheckValidity();


    // glm::dvec3 newCoords = mesh.V.at(vid).xyz();
    glm::dvec3 newCoords(0.0, 0.0, 0.0);
    // for (auto id: verticesToSplit) {
    //     newCoords += mesh.V.at(id).xyz();
    // }
    for (auto id: verticesToChange) {
        newCoords += mesh.V.at(id).xyz();
    }
    // newCoords /= (verticesToSplit.size()+verticesToChange.size()+1);
    newCoords /= (verticesToChange.size());

    Vertex newV_(newCoords);
    newV_.id = mesh.V.size();
    mesh.V.push_back(newV_);

    auto& newV = mesh.V.at(mesh.V.size()-1);
    auto& v = mesh.V.at(vid);

    // std::cout << "splitting at vertex: " << v.id << std::endl;
    // std::cout << "vertex neighbors: " << v.N_Vids.size() << " ";
    // for (auto id: v.N_Vids) std::cout << id << " ";
    // std::cout << std::endl;
    // std::cout << "vertices to split: " << verticesToSplit.size() << " ";
    // for (auto id: verticesToSplit) std::cout << id << " ";
    // std::cout << std::endl;
    // std::cout << "vertices to change: " << verticesToChange.size() << " ";
    // for (auto id: verticesToChange) std::cout << id << " ";
    // std::cout << std::endl;

    std::vector<size_t> newEdges;
    for (auto id: verticesToSplit) {
        Edge edge(std::vector<size_t>{id, newV.id});
        edge.id = mesh.E.size();
        mesh.E.push_back(edge);
        newEdges.push_back(edge.id);
        AddContents(mesh.V.at(id).N_Vids, std::vector<size_t>{newV.id});
        AddContents(mesh.V.at(id).N_Eids, std::vector<size_t>{edge.id});
    }

    std::vector<size_t> edgesToChange;
    for (auto id: verticesToChange) {
        for (auto eid: v.N_Eids) {
            auto& e = mesh.E.at(eid);
            if (e.Vids.at(0) == id || e.Vids.at(1) == id) {
                edgesToChange.push_back(e.id);
            }
        }
        UpdateContents(mesh.V.at(id).N_Vids, std::vector<size_t>{v.id});
        AddContents(mesh.V.at(id).N_Vids, std::vector<size_t>{newV.id});
    }

    std::vector<size_t> edgesToSplit;
    for (auto id: verticesToSplit) {
        for (auto eid: v.N_Eids) {
            auto& e = mesh.E.at(eid);
            if (e.Vids.at(0) == id || e.Vids.at(1) == id) {
                edgesToSplit.push_back(e.id);
            }
        }
    }

    // std::cout << "edge neighbors: " << v.N_Eids.size() << " ";
    // for (auto id: v.N_Eids) std::cout << id << " ";
    // std::cout << std::endl;
    // std::cout << "edges to split: " << edgesToSplit.size() << std::endl;
    // for (auto id: edgesToSplit) {
    //     auto& e = mesh.E.at(id);
    //     std::cout << id << " " << e.Vids.at(0) << " " << e.Vids.at(1) << std::endl;
    // };
    // std::cout << "edges to change: " << edgesToChange.size() << std::endl;
    // for (auto id: edgesToChange) {
    //     auto& e = mesh.E.at(id);
    //     std::cout << id << " " << e.Vids.at(0) << " " << e.Vids.at(1) << std::endl;
    // }


    std::vector<size_t> facesToChange;
    for (auto id: edgesToChange) {
        auto& e = mesh.E.at(id);
        AddContents(facesToChange, e.N_Fids);
    }

    // std::cout << "face neighbors: " << v.N_Fids.size() << " ";
    // for (auto id: v.N_Fids) std::cout << id << " ";
    // std::cout << std::endl;
    // std::cout << "faces to change: " << facesToChange.size() << std::endl;
    // for (auto id: facesToChange) {
    //     auto& f = mesh.F.at(id);
    //     std::cout << id << " " << f.Vids.at(0) << " " << f.Vids.at(1) << " " << f.Vids.at(2) << " " << f.Vids.at(3) << std::endl;
    // }

    std::vector<size_t> newFaceVids = {v.id};
    // std::cout << "newFaceVids: " << newFaceVids.at(0) << std::endl;
    for (auto id: facesToChange) {
        auto& f = mesh.F.at(id);
        int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
        size_t svid = f.Vids.at((idx+1)%f.Vids.size());
        if (std::find(verticesToSplit.begin(), verticesToSplit.end(), svid) != verticesToSplit.end()) {
            newFaceVids.push_back(svid);
            break;
        }
        std::vector<size_t> neighborFaces;
        for (auto fvid: f.Vids) {
            if (fvid == v.id) continue;
            AddContents(neighborFaces, mesh.V.at(fvid).N_Fids);
        }
        UpdateContents(f.N_Fids, GetDifference(v.N_Fids, neighborFaces));
    }
    // std::cout << "newFaceVids: " << newFaceVids.at(0) << " " << newFaceVids.at(1) << std::endl;

    std::vector<size_t> facesToExclude = GetDifference(v.N_Fids, facesToChange);
    for (auto id: facesToExclude) {
        auto& f = mesh.F.at(id);
        std::vector<size_t> neighborFaces;
        for (auto fvid: f.Vids) {
            if (fvid == v.id) continue;
            AddContents(neighborFaces, mesh.V.at(fvid).N_Fids);
        }
        UpdateContents(f.N_Fids, GetDifference(v.N_Fids, neighborFaces));
    }

    newFaceVids.push_back(newV.id);
    // std::cout << "newFaceVids: " << newFaceVids.at(0) << " " << newFaceVids.at(1) << " " << newFaceVids.at(2) << std::endl;
    size_t lastVid = GetDifference(verticesToSplit, std::vector<size_t>(newFaceVids.begin(), newFaceVids.end())).at(0);
    newFaceVids.push_back(lastVid);
    // std::cout << "newFaceVids: " << newFaceVids.at(0) << " " << newFaceVids.at(1) << " " << newFaceVids.at(2) << " " << newFaceVids.at(3) << std::endl;
    Face newF_(newFaceVids);
    newF_.id = mesh.F.size();
    mesh.F.push_back(newF_);
    // std::cout << "newV: " << newV.id << std::endl;
    // std::cout << "newF: " << newF_.id << " ";
    // for (auto fvid: newF_.Vids) {
    //     std::cout << fvid << " ";
    // }
    // std::cout  << std::endl;
    auto& newF = mesh.F.at(mesh.F.size()-1);
    AddContents(mesh.F.at(newF.id).Eids, edgesToSplit);
    AddContents(mesh.F.at(newF.id).Eids, newEdges);
    for (auto id: newF.Vids) {
        AddContents(mesh.F.at(newF.id).N_Fids, mesh.V.at(id).N_Fids);
    }
    AddContents(newV.N_Vids, verticesToSplit);
    AddContents(newV.N_Vids, verticesToChange);
    AddContents(newV.N_Eids, newEdges);
    AddContents(newV.N_Eids, edgesToChange);
    AddContents(newV.N_Fids, facesToChange);
    AddContents(newV.N_Fids, std::vector<size_t>{newF.id});

    for (auto id: facesToChange) {
        auto& f = mesh.F.at(id);
        int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
        f.Vids.at(idx) = newV.id;
    }

    for (auto id: edgesToChange) {
        auto& e = mesh.E.at(id);
        int idx = std::distance(e.Vids.begin(), std::find(e.Vids.begin(), e.Vids.end(), v.id));
        e.Vids.at(idx) = newV.id;
    }
    
    for (auto eid: edgesToSplit) {
        auto& e = mesh.E.at(eid);
        AddContents(e.N_Fids, std::vector<size_t>{newF.id});
        for (auto fid: facesToChange) {
            auto& f = mesh.F.at(fid);
            if (std::find(f.Eids.begin(), f.Eids.end(), e.id) != f.Eids.end()) {
                UpdateContents(f.Eids, std::vector<size_t>{e.id});
                UpdateContents(e.N_Fids, std::vector<size_t>{f.id});
                break;
            }
        }
    }
    for (auto eid: newEdges) {
        auto& e = mesh.E.at(eid);
        AddContents(e.N_Fids, std::vector<size_t>{newF.id});
        size_t nvid = e.Vids.at(0) == v.id ? e.Vids.at(1) : e.Vids.at(0);
        for (auto fid: facesToChange) {
            auto& f = mesh.F.at(fid);
            if (std::find(f.Vids.begin(), f.Vids.end(), nvid) != f.Vids.end()) {
                AddContents(f.Eids, std::vector<size_t>{e.id});
                AddContents(e.N_Fids, std::vector<size_t>{f.id});
                break;
            }
        }
    }
    UpdateContents(v.N_Vids, verticesToChange);
    UpdateContents(v.N_Eids, edgesToChange);
    UpdateContents(v.N_Fids, facesToChange);
    for (auto id: newF.Vids) {
        AddContents(mesh.V.at(id).N_Fids, std::vector<size_t>{newF.id});
    }
    AddContents(v.N_Fids, std::vector<size_t>{newF.id});
    for (auto id: newF.Vids) {
        SetSingularity(id);
        auto& fv = mesh.V.at(id);
        // ToSmooth.push_back(id);
        AddContents(ToSmooth, std::vector<size_t>{id});
        AddContents(ToSmooth, fv.N_Vids);
        // ToSmooth.insert(ToSmooth.end(), fv.N_Vids.begin(), fv.N_Vids.end());
    }
    // std::cout << "new face vertices neighbors" << std::endl;
    // for (auto id: mesh.F.at(newF_.id).Vids) {
    //     std::cout << "vertex: " << id << std::endl;
    //     std::cout << mesh.V.at(id).N_Vids.size() << std::endl;
    //     std::cout << mesh.V.at(id).N_Eids.size() << std::endl;
    //     std::cout << mesh.V.at(id).N_Fids.size() << std::endl;
    // }
    // std::cout << "new face edges " << mesh.F.at(newF_.id).Eids.size() << std::endl;
    // std::cout << "new face edges neighbor faces" << std::endl;
    // for (auto id: mesh.F.at(newF_.id).Eids) {
    //     std::cout << mesh.E.at(id).N_Fids.size() << " " << mesh.E.at(id).N_Fids.at(0) << " " << mesh.E.at(id).N_Fids.at(1) << std::endl;
    // } 
    // std::cout << "new face neighbor faces " << mesh.F.at(newF_.id).N_Fids.size() << std::endl;
    Smooth();
}