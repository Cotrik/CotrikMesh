/*
 * QuadMeshLocalSimplify.cpp
 *
 *  Created on: Nov 6, 2018
 *      Author: cotrik
 */
#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "BaseComplexQuad.h"
#include "BaseComplexSheetQuad.h"
#include "Patches.h"
#include "ArgumentManager.h"

#include <unordered_map>
#include <unordered_set>
#include <set>
#include <queue>
#include <iostream>
#include <algorithm>

int maxValence = 6;
int minValence = 3;
void get_parallel_edgeids(const Mesh& mesh, size_t start_edge_id, size_t start_face_id,
        std::set<size_t>& parallel_edgeids, std::set<size_t>& parallel_faceids) {
    if (parallel_edgeids.find(start_edge_id) != parallel_edgeids.end()) return;
    parallel_edgeids.insert(start_edge_id);
    parallel_faceids.insert(start_face_id);
    size_t next_edge_id;
    auto& start_edge = mesh.E.at(start_edge_id);
    auto& v0 = mesh.V.at(start_edge.Vids[0]);
    auto& v1 = mesh.V.at(start_edge.Vids[1]);
    for (auto edgeid : mesh.F.at(start_face_id).Eids) {
        if (edgeid == start_edge_id) continue;
        auto& e = mesh.E.at(edgeid);
        auto& v_0 = mesh.V.at(e.Vids[0]);
        auto& v_1 = mesh.V.at(e.Vids[1]);
        if (v_0.id != v0.id && v_0.id != v1.id && v_1.id != v0.id && v_1.id != v1.id) {
            next_edge_id = edgeid;
            break;
        }
    }
    auto& next_edge = mesh.E.at(next_edge_id);
    size_t next_face_id = next_edge.N_Fids[0] == start_face_id ? next_edge.N_Fids[1] : next_edge.N_Fids[0];
    get_parallel_edgeids(mesh, next_edge_id, next_face_id, parallel_edgeids, parallel_faceids);
}

size_t get_faceid(const Mesh& mesh, size_t vid, size_t exclude_vid) {
    auto& v = mesh.V.at(vid);
    size_t res;
    for (auto fid : v.N_Fids) {
        auto& f = mesh.F.at(fid);
        bool found_exclude_vid = false;
        for (auto fvid : f.Vids)
            if (fvid == exclude_vid) {
                found_exclude_vid = true;
                break;
            }
        if (!found_exclude_vid) {
            res = fid;
            break;
        }
    }
    return res;
}

size_t get_diagnal_vid(const Mesh& mesh, size_t vid, size_t fid) {
    auto& f = mesh.F.at(fid);
    for (int i = 0; i < 4; ++i) {
        if (f.Vids[i] == vid) return f.Vids.at((i + 2) % 4);
    }
    std::cerr << "ERROR get_diagnal_vid\n";
    return MAXID;
}

std::vector<size_t> get_collapse_vids(Mesh& mesh, size_t vid, size_t eid) {
    std::vector<size_t> p;
    auto& v = mesh.V.at(vid);
    auto& e = mesh.E.at(eid);
    for (auto fid : e.N_Fids) {
        auto& f = mesh.F.at(fid);
        for (auto feid : f.Eids) {
            auto& fe = mesh.E.at(feid);
            std::set<int> vids(fe.Vids.begin(), fe.Vids.end());
            vids.insert(e.Vids.begin(), e.Vids.end());
            if (vids.size() == 4) {
                for (auto nvid : v.N_Vids) {
                    if (nvid == fe.Vids[0] || nvid == fe.Vids[1]) {
                        p.push_back(nvid);
                        break;
                    }
                }
                break;
            }
        }
    }
    return p;
}

bool can_collapse_vids(Mesh& mesh, const std::vector<size_t>& vids, size_t target_vid) {
    return ((mesh.V.at(vids[0]).N_Fids.size() + mesh.V.at(vids[1]).N_Fids.size() < 10) &&
            (mesh.V.at(vids[0]).N_Fids.size() + mesh.V.at(vids[1]).N_Fids.size() > 6));
}

void collapse_vids(Mesh& mesh, const std::vector<size_t>& vids, size_t target_vid) {
    for (auto vid : vids) {
        auto& v = mesh.V.at(vid);
        for (auto fid : v.N_Fids) {
            auto& f = mesh.F.at(fid);
            for (auto& fvid : f.Vids) {
                if (fvid == vid) {
                    fvid = target_vid;
                    break;
                }
            }
        }
    }
}

bool can_collapse(Mesh& mesh, const std::vector<size_t>& linkVids, const std::vector<size_t>& linkEids) {
    std::set<size_t> vids;
    for (size_t i = 0; i < linkEids.size(); ++i) {
        auto vid = linkVids[i];
        auto eid = linkEids[i];
        auto p = get_collapse_vids(mesh, vid, eid);
        if (!can_collapse_vids(mesh, p, vid)) return false;
        vids.insert(p.begin(), p.end());
    }
    auto vid = linkVids.back();
    auto eid = linkEids.back();
    auto p = get_collapse_vids(mesh, vid, eid);
    if (!can_collapse_vids(mesh, p, vid)) return false;
    vids.insert(p.begin(), p.end());

    if (vids.size() < 2 * linkVids.size()) return false; // tangent;
    return true;
}

void collapse(Mesh& mesh, const std::vector<size_t>& linkVids, const std::vector<size_t>& linkEids) {
    for (size_t i = 0; i < linkEids.size(); ++i) {
        auto vid = linkVids[i];
        auto eid = linkEids[i];
        auto p = get_collapse_vids(mesh, vid, eid);
        collapse_vids(mesh, p, vid);
    }
    auto vid = linkVids.back();
    auto eid = linkEids.back();
    auto p = get_collapse_vids(mesh, vid, eid);
    collapse_vids(mesh, p, vid);
}

bool can_collapse_vids_with_feature_preserved(Mesh& mesh, const std::vector<size_t>& vids, size_t target_vid) {
    int count  = 0;
    if (mesh.V.at(vids[0]).isCorner) ++count;
    if (mesh.V.at(vids[1]).isCorner) ++count;
    if (mesh.V.at(target_vid).isCorner) ++count;
	const int max_mergeValence = maxValence * 2 + 1;
    return ((mesh.V.at(vids[0]).N_Fids.size() + mesh.V.at(vids[1]).N_Fids.size() < max_mergeValence) &&
            (mesh.V.at(vids[0]).N_Fids.size() + mesh.V.at(vids[1]).N_Fids.size() > 6) &&
            count <= 1);
}

void collapse_vids_with_feature_preserved(Mesh& mesh, const std::vector<size_t>& vids, size_t target_vid) {
    auto target_vid_copy = target_vid;
    if (mesh.V.at(vids[0]).type == FEATURE) target_vid = vids[0];
    if (mesh.V.at(vids[1]).type == FEATURE) target_vid = vids[1];
    if (mesh.V.at(target_vid_copy).type == FEATURE) target_vid = target_vid_copy;
    if (mesh.V.at(vids[0]).type == CORNER) target_vid = vids[0];
    if (mesh.V.at(vids[1]).type == CORNER) target_vid = vids[1];
    if (mesh.V.at(target_vid_copy).type == CORNER) target_vid = target_vid_copy;
    for (auto vid : vids) {
        auto& v = mesh.V.at(vid);
        for (auto fid : v.N_Fids) {
            auto& f = mesh.F.at(fid);
            for (auto& fvid : f.Vids) {
                if (fvid == vid) {
                    fvid = target_vid;
                    break;
                }
            }
        }
    }
}


bool can_collapse_with_feature_preserved(Mesh& mesh, const std::vector<size_t>& linkVids, const std::vector<size_t>& linkEids) {
    std::set<size_t> vids;
    for (size_t i = 0; i < linkEids.size(); ++i) {
        auto vid = linkVids[i];
        auto eid = linkEids[i];
        auto p = get_collapse_vids(mesh, vid, eid);
        if (!can_collapse_vids_with_feature_preserved(mesh, p, vid)) return false;
        vids.insert(p.begin(), p.end());
    }
    auto vid = linkVids.back();
    auto eid = linkEids.back();
    auto p = get_collapse_vids(mesh, vid, eid);
    if (!can_collapse_vids_with_feature_preserved(mesh, p, vid)) return false;
    vids.insert(p.begin(), p.end());

    if (vids.size() < 2 * linkVids.size()) return false; // tangent;
    return true;
}

void collapse_with_feature_preserved(Mesh& mesh, const std::vector<size_t>& linkVids, const std::vector<size_t>& linkEids) {
    for (size_t i = 0; i < linkEids.size(); ++i) {
        auto vid = linkVids[i];
        auto eid = linkEids[i];
        auto p = get_collapse_vids(mesh, vid, eid);
        collapse_vids_with_feature_preserved(mesh, p, vid);
    }
    auto vid = linkVids.back();
    auto eid = linkEids.back();
    auto p = get_collapse_vids(mesh, vid, eid);
    collapse_vids_with_feature_preserved(mesh, p, vid);
}


//void get_collapse_vids2(Mesh& mesh, size_t vid, size_t eid, std::vector<std::vector<size_t>>& twolinkVids) {
//    std::vector<size_t> p;
//    auto& v = mesh.V.at(vid);
//    auto& e = mesh.E.at(eid);
//    for (auto fid : e.N_Fids) {
//        auto& f = mesh.F.at(fid);
//        for (auto feid : f.Eids) {
//            auto& fe = mesh.E.at(feid);
//            std::set<int> vids(fe.Vids.begin(), fe.Vids.end());
//            vids.insert(e.Vids.begin(), e.Vids.end());
//            if (vids.size() == 4) {
//                for (auto nvid : v.N_Vids) {
//                    if (nvid == fe.Vids[0] || nvid == fe.Vids[1]) {
//                        p.push_back(nvid);
//                        break;
//                    }
//                }
//                break;
//            }
//        }
//    }
//    auto& v0 = mesh.V.at(twolinkVids[0].back());
//    auto& v1 = mesh.V.at(twolinkVids[1].back());
//    for (auto nvid : v0.N_Vids)
//        for (auto pvid : p) {
//            if (nvid == pvid) {
//                twolinkVids[0].push_back(pvid);
//                break;
//            }
//        }
//
//    for (auto nvid : v1.N_Vids)
//        for (auto pvid : p) {
//            if (nvid == pvid) {
//                twolinkVids[1].push_back(pvid);
//                break;
//            }
//        }
//    //return p;
//}
//
//void collapse_with_feature_preserved(Mesh& mesh, std::vector<std::vector<size_t>>& threelinkVids) {
//    bool haveFeature[3] = {false, false, false};
//    bool haveCorners[3] = {false, false, false};
//
//    for (auto i = 0; i < 3; ++i) {
//        for (auto vid : threelinkVids[i]) {
//            if (mesh.V.at(vid).type == FEATURE) {
//                haveFeature[i] = true;
//                break;
//            }
//        }
//    }
//
//    for (auto i = 0; i < 3; ++i) {
//        for (auto vid : threelinkVids[i]) {
//            if (mesh.V.at(vid).type == CORNER) {
//                haveCorners[i] = true;
//                break;
//            }
//        }
//    }
//
//}
//void collapse_with_feature_preserved2(Mesh& mesh, const vector<size_t>& linkVids, const vector<size_t>& linkEids) {
//    std::vector<std::vector<size_t>> twolinkVids(2);
//    {
//        auto vid = linkVids.front();
//        auto eid = linkEids.front();
//        auto p = get_collapse_vids(mesh, vid, eid);
//        twolinkVids[0].push_back(p[0]);
//        twolinkVids[1].push_back(p[1]);
//        //collapse_vids_with_feature_preserved(mesh, p, vid);
//    }
//    for (size_t i = 1; i < linkEids.size(); ++i) {
//        auto vid = linkVids[i];
//        auto eid = linkEids[i];
//        get_collapse_vids2(mesh, vid, eid, twolinkVids);
//        //collapse_vids_with_feature_preserved(mesh, p, vid);
//    }
//    {
//        auto vid = linkVids.back();
//        auto eid = linkEids.back();
//        get_collapse_vids2(mesh, vid, eid, twolinkVids);
//        // collapse_vids_with_feature_preserved(mesh, p, vid);
//    }
//    twolinkVids.push_back(linkVids);
//
//    collapse_with_feature_preserved(mesh, twolinkVids);
//
//
//}

void get_feature(Mesh& mesh) {
    const double cosangle = 0.866025404;
    mesh.SetCosAngleThreshold(cosangle);
//    mesh.LabelSurface();
//    mesh.LabelSharpEdges(true);
    mesh.GetNormalOfSurfaceFaces();
    mesh.GetNormalOfSurfaceVertices();

    Patches patches(mesh);
    patches.Extract();
    std::vector<size_t> edgeIds;
    for (size_t i = 0; i < patches.patches.size(); i++) {
        const Patch& patch = patches.patches.at(i);
        std::copy(patch.edgeIds.begin(), patch.edgeIds.end(), back_inserter(edgeIds));
    }
    std::set<size_t> sharpEdgeVertexIds;
    for (auto& eid : edgeIds) {
        auto& e = mesh.E.at(eid);
        sharpEdgeVertexIds.insert(e.Vids.begin(), e.Vids.end());
    }
//    {
//        MeshFileWriter writer(mesh, "sharpEdgeVertex.vtk");
//        writer.WriteVerticesVtk(sharpEdgeVertexIds);
//    }

    std::vector<size_t> roundVertexIds;
    std::map<size_t, size_t> roundVertexIds_count;
    for (const Patch& patch : patches.patches) {
        for (auto eid : patch.edgeIds) {
            auto& e = mesh.E.at(eid);
            ++roundVertexIds_count[e.Vids[0]];
            ++roundVertexIds_count[e.Vids[1]];
//            mesh.V.at(e.Vids[0]).type = FEATURE;
//            mesh.V.at(e.Vids[1]).type = FEATURE;
        }
    }
    for (auto& v : mesh.V)
        if (v.type == CORNER || v.type == FEATURE) roundVertexIds.push_back(v.id);
    for (auto& p : roundVertexIds_count)
        if (p.second >= 6) {
            roundVertexIds.push_back(p.first);
            mesh.V.at(p.first).isCorner = true;
            mesh.V.at(p.first).type = CORNER;
        }
    for (const Patch& patch : patches.patches) {
        for (auto eid : patch.edgeIds) {
            auto& e = mesh.E.at(eid);
            if (mesh.V.at(e.Vids[0]).type != CORNER) mesh.V.at(e.Vids[0]).type = FEATURE;
            if (mesh.V.at(e.Vids[1]).type != CORNER) mesh.V.at(e.Vids[1]).type = FEATURE;
        }
    }
    {
        MeshFileWriter writer(mesh, "cornerVertex.vtk");
        writer.WriteVerticesVtk(roundVertexIds);
    }
}

void collapse_one_face_with_sharp_feature_preserved(Mesh& mesh, std::set<size_t>& canceledFids) {
    bool all_singular = true;
    for (auto& f : mesh.F) {
        for (auto vid : f.Vids)
            if (!mesh.V.at(vid).isSingularity) {
                all_singular = false;
                break;
            }
        if (all_singular) {
            auto& v0 = mesh.V.at(f.Vids[0]);
            auto& v1 = mesh.V.at(f.Vids[1]);
            auto& v2 = mesh.V.at(f.Vids[2]);
            auto& v3 = mesh.V.at(f.Vids[3]);
            size_t numOfCorners = 0;
            for (auto vid : f.Vids)
                if (mesh.V.at(vid).isCorner) ++numOfCorners;
            if (numOfCorners <= 1)
                if ((v0.N_Fids.size() == 3 && v1.N_Fids.size() == 5 && v2.N_Fids.size() == 3 && v3.N_Fids.size() == 5)
                        || (v0.N_Fids.size() == 5 && v1.N_Fids.size() == 3 && v2.N_Fids.size() == 5 && v3.N_Fids.size() == 3)) {
                    auto center_vid = v0.N_Fids.size() == 3 ? v0.id : v1.id;
                    if (numOfCorners == 1) {
                        for (auto vid : f.Vids)
                            if (mesh.V.at(vid).isCorner) {
                                center_vid = vid;
                                break;
                            }
                    }
                    for (auto vid : f.Vids) {
                        auto& v = mesh.V[vid];
                        for (auto n_fid : v.N_Fids) {
                            if (n_fid == f.id) continue;
                            auto& n_f = mesh.F.at(n_fid);
                            for (auto& n_vid : n_f.Vids) {
                                if (n_vid == vid) {
                                    n_vid = center_vid;
                                    break;
                                }
                            }
                        }
                    }
                }
            canceledFids.insert(f.id);
            break;
        }
    }
}

void update(Mesh& mesh, std::set<size_t>& canceledFids) {
    std::vector<size_t> FaceIds;
    FaceIds.reserve(mesh.F.size());
    for (auto& f : mesh.F)
        if (canceledFids.find(f.id) == canceledFids.end())
            FaceIds.push_back(f.id);
    std::vector<Vertex> newV(mesh.V.size());
    std::vector<Face> newF(FaceIds.size());
    std::vector<Cell> newC(FaceIds.size());
    for (size_t i = 0; i < mesh.V.size(); ++i) {
        newV.at(i).id = i;
        newV.at(i) = mesh.V.at(i).xyz();
        newV.at(i).type = mesh.V.at(i).type;
        newV.at(i).isCorner = mesh.V.at(i).isCorner;
    }
    for (size_t i = 0; i < FaceIds.size(); ++i) {
        newF.at(i).id = i;
        newF.at(i).Vids = mesh.F.at(FaceIds[i]).Vids;
    }
    for (size_t i = 0; i < FaceIds.size(); ++i) {
        newC.at(i).id = i;
        newC.at(i).Vids = mesh.F.at(FaceIds[i]).Vids;
        newC.at(i).cellType = VTK_QUAD;
    }
    mesh.V.clear();
    mesh.E.clear();
    mesh.F.clear();
    mesh.C.clear();

    mesh.V = newV;
    mesh.F = newF;
    mesh.C = newC;
    canceledFids.clear();
}
int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: QuadMeshLocalSimplify quad.vtk simplified.vtk iters=<1> maxValence=<6> featurePreserved=true" << std::endl;
        return -1;
    }
    ArgumentManager argumentManager(argc, argv);
    std::string input = argv[1];
    std::string output = argv[2];
    int iters = 1;
	bool featurePreserved = true;
    auto strIters = argumentManager.get("iters");
    if (!strIters.empty()) iters = std::stoi(strIters);
	auto strmaxValence = argumentManager.get("maxValence");
	if (!strmaxValence.empty()) maxValence = std::stoi(strmaxValence);
	auto strfeaturePreserved = argumentManager.get("featurePreserved");
	if (!strfeaturePreserved.empty()) featurePreserved = strfeaturePreserved == "false" ? false : true;
    std::cout << "---------------------------------------" << std::endl;
    std::cout << "input  = " << input << std::endl;
    std::cout << "output = " << output << std::endl;
    std::cout << "---------------------------------------" << std::endl;
    std::set<size_t> canceledFids;
    int iter = 0;
    MeshFileReader reader(input.c_str());
    Mesh& mesh = (Mesh&) reader.GetMesh();
    while (iters--) {
        mesh.CompressWithFeaturePreserved();
        mesh.BuildAllConnectivities();
        mesh.ExtractBoundary();
        mesh.ExtractSingularities();
        if (iter == 0 && featurePreserved)
            get_feature(mesh);
        BaseComplexQuad baseComplex(mesh);
        baseComplex.Build();
        for (auto t = 0; t < 2; ++t) {
            size_t id = 0;
            for (const auto& link : baseComplex.separatedVertexIdsLink) {
                auto& v_front = mesh.V.at(link.front());
                auto& v_back = mesh.V.at(link.back());
                if ((v_front.N_Fids.size() <= 5 && v_back.N_Fids.size() <= 5) && v_front.N_Fids.size() != v_back.N_Fids.size()) {
                    ; //ofs << 0 << std::endl;
                }
                else if (v_front.N_Fids.size() == 3 && v_back.N_Fids.size() == 3) {
                    // ofs << 1 << std::endl;
                    auto v_front_fid = get_faceid(mesh, v_front.id, link[1]);
                    auto v_front_fvid = get_diagnal_vid(mesh, v_front.id, v_front_fid);

                    auto v_back_fid = get_faceid(mesh, v_back.id, link[link.size() - 2]);
                    auto v_back_fvid = get_diagnal_vid(mesh, v_back.id, v_back_fid);

                    bool condition = false;
                    if ((t == 0 && mesh.V.at(v_front_fvid).N_Fids.size() == 5 && mesh.V.at(v_back_fvid).N_Fids.size() == 5)
                            || (t == 1 && (mesh.V.at(v_front_fvid).N_Fids.size() > 3 && mesh.V.at(v_back_fvid).N_Fids.size() > 3)))
                        condition = true;
                    if (condition) {
                        if (can_collapse_with_feature_preserved(mesh, baseComplex.separatedVertexIdsLink.at(id),
                                baseComplex.separatedEdgeIdsLink.at(id))) {
                            for (auto vid : link) {
                                auto& v = mesh.V.at(vid);
                                canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
                            }
                            collapse_with_feature_preserved(mesh, baseComplex.separatedVertexIdsLink.at(id),
                                    baseComplex.separatedEdgeIdsLink.at(id));

                            break;
                        }
                    }
                } else if (v_front.N_Fids.size() == 5 && v_back.N_Fids.size() == 5) {
                    ; //ofs << 2 << std::endl;
                } else {
                    ; //ofs << 3 << std::endl;
                }
                ++id;
            }
            if (!canceledFids.empty()) break;
        }
        std::cout << "iter = " << iter++ << std::endl;
        if (canceledFids.empty()) collapse_one_face_with_sharp_feature_preserved(mesh, canceledFids);
        if (canceledFids.empty()) break;
        update(mesh, canceledFids);
    }

    MeshFileWriter writer(mesh, output.c_str());
    writer.WriteFile();
    return 0;
}

int main_(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: QuadMeshLocalSimplify quad.vtk simplified.vtk iters=100" << std::endl;
        return -1;
    }
    ArgumentManager argumentManager(argc, argv);
    std::string input = argv[1];
    std::string output = argv[2];
    int iters = 1;
    auto strIters = argumentManager.get("iters");
    if (!strIters.empty()) iters = std::stoi(strIters);
    std::cout << "---------------------------------------" << std::endl;
    std::cout << "input  = " << input << std::endl;
    std::cout << "output = " << output << std::endl;
    std::cout << "---------------------------------------" << std::endl;
    std::set<size_t> canceledFids;
    int iter = 0;
    while (iters--) {
        MeshFileReader reader(input.c_str());
        Mesh& mesh = (Mesh&) reader.GetMesh();
        mesh.CompressWithFeaturePreserved();
        mesh.BuildAllConnectivities();
        mesh.ExtractBoundary();
        mesh.ExtractSingularities();
        //if (iter == 0)
            get_feature(mesh);
        BaseComplexQuad baseComplex(mesh);
        baseComplex.Build();
        for (auto t = 0; t < 2; ++t) {
            size_t id = 0;
            for (const auto& link : baseComplex.separatedVertexIdsLink) {
                auto& v_front = mesh.V.at(link.front());
                auto& v_back = mesh.V.at(link.back());
                if ((v_front.N_Fids.size() <= 5 && v_back.N_Fids.size() <= 5) &&
                        v_front.N_Fids.size() != v_back.N_Fids.size()) {
                    ; //ofs << 0 << std::endl;
                }
                else if (v_front.N_Fids.size() == 3 && v_back.N_Fids.size() == 3) {
                    // ofs << 1 << std::endl;
                    auto v_front_fid = get_faceid(mesh, v_front.id, link[1]);
                    auto v_front_fvid = get_diagnal_vid(mesh, v_front.id, v_front_fid);

                    auto v_back_fid = get_faceid(mesh, v_back.id, link[link.size() - 2]);
                    auto v_back_fvid = get_diagnal_vid(mesh, v_back.id, v_back_fid);

                    bool condition = false;
                    if ((t == 0 && mesh.V.at(v_front_fvid).N_Fids.size() == 5 && mesh.V.at(v_back_fvid).N_Fids.size() == 5) ||
                            (t == 1 && (mesh.V.at(v_front_fvid).N_Fids.size() > 3 && mesh.V.at(v_back_fvid).N_Fids.size() > 3)) )
                        condition = true;
                    if (condition) {
                        if (can_collapse_with_feature_preserved(mesh, baseComplex.separatedVertexIdsLink.at(id), baseComplex.separatedEdgeIdsLink.at(id))) {
                            for (auto vid : link) {
                                auto& v = mesh.V.at(vid);
                                canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
                            }
                            collapse_with_feature_preserved(mesh, baseComplex.separatedVertexIdsLink.at(id), baseComplex.separatedEdgeIdsLink.at(id));

                            break;
                        }
                    }
                }
                else if (v_front.N_Fids.size() == 5 && v_back.N_Fids.size() == 5) {
                    ; //ofs << 2 << std::endl;
                } else {
                    ; //ofs << 3 << std::endl;
                }
                ++id;
            }
            if (!canceledFids.empty()) break;
        }
        std::cout << "iter = " << iter++ << std::endl;
        if (canceledFids.empty()) {
            bool all_singular = true;
            for (auto& f : mesh.F) {
                for (auto vid : f.Vids)
                    if (!mesh.V.at(vid).isSingularity) {
                        all_singular = false;
                        break;
                    }
                if (all_singular) {
                    auto& v0 = mesh.V.at(f.Vids[0]);
                    auto& v1 = mesh.V.at(f.Vids[1]);
                    auto& v2 = mesh.V.at(f.Vids[2]);
                    auto& v3 = mesh.V.at(f.Vids[3]);
                    if ((v0.N_Fids.size() == 3 && v1.N_Fids.size() == 5 && v2.N_Fids.size() == 3 && v3.N_Fids.size() == 5) ||
                        (v0.N_Fids.size() == 5 && v1.N_Fids.size() == 3 && v2.N_Fids.size() == 5 && v3.N_Fids.size() == 3)) {
                        auto center_vid = v0.N_Fids.size() == 3 ? v0.id : v1.id;
                        for (auto vid : f.Vids) {
                            auto& v = mesh.V[vid];
                            for (auto n_fid : v.N_Fids) {
                                if (n_fid == f.id) continue;
                                auto& n_f = mesh.F.at(n_fid);
                                for (auto& n_vid : n_f.Vids) {
                                    if (n_vid == vid) {
                                        n_vid = center_vid;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    canceledFids.insert(f.id);
                    break;
                }
            }
        }
        if (canceledFids.empty()) break;
        std::vector<size_t> FaceIds;
        FaceIds.reserve(mesh.F.size());
        for (auto& f : mesh.F)
            if (canceledFids.find(f.id) == canceledFids.end())
                FaceIds.push_back(f.id);
//        std::vector<Vertex> newV(mesh.V.size());
//        std::vector<Face> newF(FaceIds.size());
//        std::vector<Cell> newC(FaceIds.size());
//        for (size_t i = 0; i < mesh.V.size(); ++i) {
//            newV.at(i).id = i;
//            newV.at(i) = mesh.V.at(i).xyz();
//        }
//        for (size_t i = 0; i < FaceIds.size(); ++i) {
//            newF.at(i).id = i;
//            newF.at(i).Vids = mesh.F.at(FaceIds[i]).Vids;
//        }
//        for (size_t i = 0; i < FaceIds.size(); ++i) {
//            newC.at(i).id = i;
//            newC.at(i).Vids = mesh.F.at(FaceIds[i]).Vids;
//            newC.at(i).cellType = VTK_QUAD;
//        }
//        mesh.V.clear();
//        mesh.E.clear();
//        mesh.F.clear();
//        mesh.C.clear();
//
//        mesh.V = newV;
//        mesh.F = newF;
//        mesh.C = newC;

        canceledFids.clear();
        MeshFileWriter writer(mesh, output.c_str());
        writer.WriteFacesVtk(FaceIds);
        input = output;
    }

//    MeshFileWriter writer(mesh, output.c_str());
//    writer.WriteFile();
    return 0;
}


int main1(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: QuadMeshLocalSimplify quad.vtk simplified.vtk iters=100" << std::endl;
        return -1;
    }
    ArgumentManager argumentManager(argc, argv);
    std::string input = argv[1];
    std::string output = argv[2];
    int iters = 1;
    auto strIters = argumentManager.get("iters");
    if (!strIters.empty()) iters = std::stoi(strIters);
    std::cout << "---------------------------------------" << std::endl;
    std::cout << "input  = " << input << std::endl;
    std::cout << "output = " << output << std::endl;
    std::cout << "---------------------------------------" << std::endl;
    int iter = 0;
    std::set<size_t> canceledFids;
    while (iters--) {
        MeshFileReader reader(input.c_str());
        Mesh& mesh = (Mesh&) reader.GetMesh();
        mesh.RemoveUselessVertices();
        mesh.BuildAllConnectivities();
        mesh.ExtractBoundary();
        mesh.ExtractSingularities();

//    size_t start_edge_id = 62418;
//    size_t start_face_id = 29752;
//    std::set<size_t> parallel_edgeids;
//    std::set<size_t> parallel_faceids;
//    get_parallel_edgeids(mesh, start_edge_id, start_face_id, parallel_edgeids, parallel_faceids);
//        {
//            MeshFileWriter writer(mesh, output.c_str());
//            writer.WriteFacesVtk(parallel_faceids);
//        }
//    return 0;

        for (auto& f : mesh.F) {
            auto& v0 = mesh.V.at(f.Vids[0]);
            auto& v1 = mesh.V.at(f.Vids[1]);
            auto& v2 = mesh.V.at(f.Vids[2]);
            auto& v3 = mesh.V.at(f.Vids[3]);
            if ((v0.N_Fids.size() == 3 && v1.N_Fids.size() == 4 && v2.N_Fids.size() == 5 && v3.N_Fids.size() == 4) ||
                    (v0.N_Fids.size() == 5 && v1.N_Fids.size() == 4 && v2.N_Fids.size() == 3 && v3.N_Fids.size() == 4) ||
                    (v0.N_Fids.size() == 4 && v1.N_Fids.size() == 3 && v2.N_Fids.size() == 4 && v3.N_Fids.size() == 5) ||
                    (v0.N_Fids.size() == 4 && v1.N_Fids.size() == 5 && v2.N_Fids.size() == 4 && v3.N_Fids.size() == 3)) {
                bool found = false;
                for (auto vid : f.Vids) {
                    for (auto fid : mesh.V.at(vid).N_Fids) {
                        for (auto vvid : mesh.F.at(fid).Vids) {
//                            if (mesh.V.at(vvid).isSingularity) {
//                                found = true;
//                                break;
//                            }
                            for (auto n_fid : mesh.V.at(vvid).N_Fids)
                                if (canceledFids.find(n_fid) != canceledFids.end()) {
                                    found = true;
                                    break;
                                }
                            if (found) break;
                        }
                        if (found) break;
                    }
                    if (found) break;
                }
                if (!found) {
                    std::set<size_t> fids;
                    for (auto vid : f.Vids)
                        fids.insert(mesh.V.at(vid).N_Fids.begin(), mesh.V.at(vid).N_Fids.end());
                    //fids.erase(f.id);
                    std::set<size_t> vids;
                    for (auto fid : fids)
                        vids.insert(mesh.F.at(fid).Vids.begin(), mesh.F.at(fid).Vids.end());
                    for (auto vid : f.Vids)
                        vids.erase(vid);
                    for (auto vid : vids)
                        if (mesh.V.at(vid).isSingularity) {
                            found = true;
                            break;
                        }
                }
                if (!found)
                for (auto vid : f.Vids)
                    canceledFids.insert(mesh.V.at(vid).N_Fids.begin(), mesh.V.at(vid).N_Fids.end());
            }
        }
        if (canceledFids.empty()) continue;
        std::vector<size_t> FaceIds;
        FaceIds.reserve(mesh.F.size());
        for (auto& f : mesh.F)
            if (canceledFids.find(f.id) == canceledFids.end())
                FaceIds.push_back(f.id);
        MeshFileWriter writer(mesh, output.c_str());
        writer.WriteFacesVtk(FaceIds);
    }
    return 0;
}
