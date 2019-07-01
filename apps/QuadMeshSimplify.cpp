/*
 * QuadMeshSimplify.cpp
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
#include <iostream>
#include <algorithm>

Mesh Refine(const Mesh& quad_mesh, int clockwise);
static std::unordered_map<size_t, size_t> get_key_edgeId(const Mesh& mesh);
static std::unordered_map<std::string, size_t> get_key_faceId(const Mesh& mesh);
static std::string get_facekey(const Face& f) {
    std::set<size_t> vids(f.Vids.begin(), f.Vids.end());
    std::string s;
    for (auto vid : vids)
        s += std::to_string(vid) + "@";
    return s;
}

std::set<size_t> get_canceledEdgeIds(const BaseComplexSheetQuad& baseComplexSheets, std::map<size_t, size_t>& canceledFaceIds, size_t sheetId) {
    std::set<size_t> canceledEdgeIds;
    for (auto baseComplexEdgeId : baseComplexSheets.sheets_componentEdgeIds[sheetId]) {
        for (auto edgeId : baseComplexSheets.baseComplex.componentE[baseComplexEdgeId].eids_link) {
            auto& e = baseComplexSheets.baseComplex.mesh.E[edgeId];
            for (auto n_fid : e.N_Fids)
                ++canceledFaceIds[n_fid];
            canceledEdgeIds.insert(edgeId);
        }
    }
    return canceledEdgeIds;
}

bool can_collapse(const BaseComplexSheetQuad& baseComplexSheets, std::map<size_t, size_t>& canceledFaceIds, size_t sheetId) {
    const auto& baseComplexFIds = baseComplexSheets.sheets_componentFaceIds[sheetId];
    std::set<size_t> baseComplexFaceIds_set(baseComplexFIds.begin(), baseComplexFIds.end());
    std::set<size_t> canceledEdgeIds;
    std::set<size_t> baseComplexVIds;
    for (auto baseComplexFId : baseComplexFIds)
        for (auto baseComplexVId : baseComplexSheets.baseComplex.componentF[baseComplexFId].Vids)
            baseComplexVIds.insert(baseComplexVId);

    for (auto baseComplexVId : baseComplexVIds) {
        bool all_in = true;
        const auto& bCFids = baseComplexSheets.baseComplex.componentV[baseComplexVId].N_Fids;
        for (auto bCFid : bCFids)
            if (baseComplexFaceIds_set.find(bCFid) == baseComplexFaceIds_set.end()) {
                all_in = false;
                break;
            }
        if (all_in)
            return false;
    }

    const auto& mesh = baseComplexSheets.baseComplex.mesh;
    for (auto& item : canceledFaceIds) {
        if (item.second >= 4) {
            auto& f = mesh.F.at(item.first);
            for (auto vid : f.Vids)
                if (vid >= mesh.V.size()) return false;
        }
    }
    for (auto edgeId : canceledEdgeIds) {
        auto& e = mesh.E[edgeId];
        for (auto vid : e.Vids)
            if (vid >= mesh.V.size()) return false;
    }

    return true;
}

void collapse(Mesh& mesh, std::unordered_map<size_t, size_t>& key_edgeId, std::unordered_map<std::string, size_t>& key_faceId,
        std::map<size_t, size_t>& canceledFaceIds, std::set<size_t>& canceledEdgeIds) {
    for (auto& item : canceledFaceIds) {
        if (item.second >= 4) {
            auto& f  = mesh.F.at(item.first);
            auto key = get_facekey(f);
            auto centerVid = mesh.V.size() + mesh.E.size() + key_faceId[key];
            for (auto vid : f.Vids) {
                // if (vid >= mesh.V.size()) continue;
                auto& v = mesh.V.at(vid);
                for (auto n_fid : v.N_Fids) {
                    auto& n_f = mesh.F.at(n_fid);
                    for (auto& n_vid : n_f.Vids)
                        if (n_vid == vid) n_vid = centerVid;
                }
            }
            for (auto eid : f.Eids)
                if (canceledEdgeIds.find(eid) != canceledEdgeIds.end()) canceledEdgeIds.erase(eid);
        }
    }
    for (auto edgeId : canceledEdgeIds) {
        auto& e = mesh.E[edgeId];
        auto& v0 = mesh.V[e.Vids[0]];
        auto& v1 = mesh.V[e.Vids[1]];
        auto key = (e.Vids[0] << 32) | e.Vids[1];
        if (key_edgeId.find(key) == key_edgeId.end()) std::cout << "Edge search Error !" << std::endl;
        auto centerVid = mesh.V.size() + key_edgeId[key];
        for (auto vid : e.Vids) {
            auto& v = mesh.V.at(vid);
            for (auto n_fid : v.N_Fids) {
                auto& n_f = mesh.F.at(n_fid);
                for (auto& n_vid : n_f.Vids)
                    if (n_vid == vid) n_vid = centerVid;
            }
        }
    }
}

bool can_collapse_with_feature_preserved(const BaseComplexSheetQuad& baseComplexSheets, std::map<size_t, size_t>& canceledFaceIds, size_t sheetId) {
    const auto& baseComplexFIds = baseComplexSheets.sheets_componentFaceIds[sheetId];
    for (auto baseComplexFId : baseComplexFIds) {
        int count = 0;
        for (auto baseComplexVId : baseComplexSheets.baseComplex.componentF[baseComplexFId].Vids) {
            auto vid = baseComplexSheets.baseComplex.Vids.at(baseComplexVId);
            auto& v = baseComplexSheets.baseComplex.mesh.V.at(vid);
            if (v.type == FEATURE || v.type == CORNER) ++count;
        }
        if (count > 3) return false;
    }

    std::set<size_t> baseComplexFaceIds_set(baseComplexFIds.begin(), baseComplexFIds.end());
    std::set<size_t> canceledEdgeIds;
    std::set<size_t> baseComplexVIds;
    for (auto baseComplexFId : baseComplexFIds)
        for (auto baseComplexVId : baseComplexSheets.baseComplex.componentF[baseComplexFId].Vids)
            baseComplexVIds.insert(baseComplexVId);

    for (auto baseComplexVId : baseComplexVIds) {
//        bool all_in = true;
//        const auto& bCFids = baseComplexSheets.baseComplex.componentV[baseComplexVId].N_Fids;
//        for (auto bCFid : bCFids)
//            if (baseComplexFaceIds_set.find(bCFid) == baseComplexFaceIds_set.end()) {
//                all_in = false;
//                break;
//            }
//        if (all_in)
//            return false;
        int count = 0;
        const auto& bCFids = baseComplexSheets.baseComplex.componentV[baseComplexVId].N_Fids;
        for (auto bCFid : bCFids)
            if (baseComplexFaceIds_set.find(bCFid) != baseComplexFaceIds_set.end()) {
                ++count;
            }
        if (count >= 4)
            return false;
    }

    const auto& mesh = baseComplexSheets.baseComplex.mesh;
    for (auto& item : canceledFaceIds) {
        if (item.second >= 3) {
            auto& f = mesh.F.at(item.first);
            for (auto vid : f.Vids)
                if (vid >= mesh.V.size()) return false;
            size_t count = 0;
            for (auto vid : f.Vids) {
                auto& v = mesh.V.at(vid);
                if (v.type == CORNER || v.type == FEATURE) ++count;
            }
            if (count > 1) return false;
        }
    }
    for (auto edgeId : canceledEdgeIds) {
        auto& e = mesh.E[edgeId];
        for (auto vid : e.Vids)
            if (vid >= mesh.V.size()) return false;
        size_t count = 0;
        for (auto vid : e.Vids) {
            auto& v = mesh.V.at(vid);
            if (v.type == CORNER || v.type == FEATURE) ++count;
        }
        if (count > 1) return false;
    }

    return true;
}

void collapse_with_feature_preserved(Mesh& mesh, std::unordered_map<size_t, size_t>& key_edgeId, std::unordered_map<std::string, size_t>& key_faceId,
        std::map<size_t, size_t>& canceledFaceIds, std::set<size_t>& canceledEdgeIds) {
    for (auto& item : canceledFaceIds) {
        if (item.second >= 4) {
            auto& f  = mesh.F.at(item.first);
            auto key = get_facekey(f);
            auto centerVid = mesh.V.size() + mesh.E.size() + key_faceId[key];
            size_t featureType = 0;
            for (auto vid : f.Vids) {
                auto& v = mesh.V.at(vid);
                if (v.type == FEATURE && featureType == 0) centerVid = vid;
                if (v.type == CORNER) centerVid = vid;
            }
            for (auto vid : f.Vids) {
                auto& v = mesh.V.at(vid);
                for (auto n_fid : v.N_Fids) {
                    auto& n_f = mesh.F.at(n_fid);
                    for (auto& n_vid : n_f.Vids)
                        if (n_vid == vid) n_vid = centerVid;
                }
            }
            for (auto eid : f.Eids)
                if (canceledEdgeIds.find(eid) != canceledEdgeIds.end()) canceledEdgeIds.erase(eid);
        }
    }
    for (auto edgeId : canceledEdgeIds) {
        auto& e = mesh.E[edgeId];
        auto& v0 = mesh.V[e.Vids[0]];
        auto& v1 = mesh.V[e.Vids[1]];
        auto key = (e.Vids[0] << 32) | e.Vids[1];
        if (key_edgeId.find(key) == key_edgeId.end()) std::cout << "Edge search Error !" << std::endl;
        auto centerVid = mesh.V.size() + key_edgeId[key];
        size_t featureType = 0;
        for (auto vid : e.Vids) {
            auto& v = mesh.V.at(vid);
            //if (v.type == FEATURE && featureType == 0) centerVid = vid;
            //if (v.type == CORNER) centerVid = vid;
            if (v.type != MAXID && v.type > featureType) {
                featureType = v.type;
                centerVid = vid;
            }
        }
        for (auto vid : e.Vids) {
            auto& v = mesh.V.at(vid);
            for (auto n_fid : v.N_Fids) {
                auto& n_f = mesh.F.at(n_fid);
                for (auto& n_vid : n_f.Vids)
                    if (n_vid == vid) n_vid = centerVid;
            }
        }
    }
}

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
            mesh.V.at(e.Vids[0]).type = FEATURE;
            mesh.V.at(e.Vids[1]).type = FEATURE;
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
    {
        MeshFileWriter writer(mesh, "cornerVertex.vtk");
        writer.WriteVerticesVtk(roundVertexIds);
    }
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: QuadMeshSimplify quad.vtk simplified.vtk iters=100" << std::endl;
        return -1;
    }
    ArgumentManager argumentManager(argc, argv);
    std::string input = argv[1];
    std::string output = argv[2];
    int iters = 100;
    auto strIters = argumentManager.get("iters");
    if (!strIters.empty()) iters = std::stoi(strIters);
    std::cout << "---------------------------------------" << std::endl;
    std::cout << "input  = " << input << std::endl;
    std::cout << "output = " << output << std::endl;
    std::cout << "---------------------------------------" << std::endl;

    std::vector<glm::dvec3> corners;
    int iter = 0;
    while (iters--) {
        MeshFileReader reader(input.c_str());
        Mesh& mesh = (Mesh&) reader.GetMesh();
        mesh.RemoveUselessVertices();
        mesh.BuildAllConnectivities();
        mesh.ExtractBoundary();
        mesh.ExtractSingularities();
        if (!mesh.hasBoundary)
            get_feature(mesh);
        else
        {
            if (iter == 0) {
                for (auto& v : mesh.V)
                    if (v.isBoundary && v.isSingularity) {
                        v.type = CORNER;
                        corners.push_back(v.xyz());
                    }
            } else {
                for (auto& v : mesh.V) {
                    if (!v.isBoundary) continue;
                    for (auto& cv : corners) {
                        if (glm::length(v - cv) < 1e-5) {
                            v.type = CORNER;
                            break;
                        }
                    }
                }
            }

        }
        BaseComplexQuad baseComplex(mesh);
        baseComplex.Build();

        BaseComplexSheetQuad baseComplexSheets(baseComplex);
        baseComplexSheets.Extract();

        auto dualMesh = Refine(mesh, 0);
        auto key_edgeId = get_key_edgeId(mesh);
        auto key_faceId = get_key_faceId(mesh);


        std::map<size_t, size_t> total_canceledFaceIds;
        for (int sheetId = 0; sheetId < baseComplexSheets.sheets_componentEdgeIds.size(); ++sheetId) {
            std::map<size_t, size_t> canceledFaceIds;
            std::set<size_t> canceledEdgeIds = get_canceledEdgeIds(baseComplexSheets, canceledFaceIds, sheetId);
            if (can_collapse_with_feature_preserved(baseComplexSheets, canceledFaceIds, sheetId)) {
                auto savedMesh = mesh;
                std::cout << "collapse sheet " << sheetId << "\n";
                collapse_with_feature_preserved(mesh, key_edgeId, key_faceId, canceledFaceIds, canceledEdgeIds);
                total_canceledFaceIds.insert(canceledFaceIds.begin(), canceledFaceIds.end());

                if (total_canceledFaceIds.empty()) continue;
                std::vector<size_t> FaceIds;
                FaceIds.reserve(mesh.F.size());
                for (auto& f : mesh.F)
                    if (total_canceledFaceIds.find(f.id) == total_canceledFaceIds.end())
                        FaceIds.push_back(f.id);

                std::vector<Vertex> quadV;
                std::vector<Cell> quadC;
                Cell face;
                size_t index = 0;
                for (auto fid : FaceIds) {
                    face.id = index++;
                    face.Vids = mesh.F.at(fid).Vids;
                    quadC.push_back(face);
                }
                Vertex vertex;
                index = 0;
                for (auto& v : dualMesh.V) {
                    vertex.id = index++;
                    vertex = v;
                    quadV.push_back(vertex);
                }
                Mesh newMesh(quadV, quadC, QUAD);
//                {
//                    MeshFileWriter writer(dualMesh.V, mesh.F, "temp.vtk");
//                    writer.WriteFacesVtk(FaceIds);
//                }
//                MeshFileReader newreader("temp.vtk");
//                Mesh& newMesh = (Mesh&) newreader.GetMesh();
                newMesh.RemoveUselessVertices();
                newMesh.BuildAllConnectivities();
                newMesh.ExtractBoundary();
                newMesh.ExtractSingularities();
                bool isValid = true;
                for (auto& v: newMesh.V) {
                    if (!v.isBoundary && v.N_Fids.size() < 3) {
                        isValid = false;
                        break;
                    }
//                    if (v.isBoundary && v.isCorner && v.N_Fids.size() < 2) {
//                        isValid = false;
//                        break;
//                    }
                }
                if (!isValid) {
                    std::cout << "recover previous mesh\n";
                    mesh = savedMesh;
                    total_canceledFaceIds.clear();
                } else {
                    MeshFileWriter writer(dualMesh.V, mesh.F, output.c_str());
                    writer.WriteFacesVtk(FaceIds);
                    break;
                }
            }
        }
        std::cout << "iter = " << ++iter << "\n";
        if (total_canceledFaceIds.empty()) break;
        input = output;
    }
}

//int main(int argc, char* argv[]) {
//    if (argc < 2) {
//        std::cout << "Usage: QuadMeshSimplify quad.vtk simplified.vtk" << std::endl;
//        return -1;
//    }
//    ArgumentManager argumentManager(argc, argv);
//    std::string input = argv[1];
//    std::string output = argv[2];
//    std::cout << "---------------------------------------" << std::endl;
//    std::cout << "input  = " << input << std::endl;
//    std::cout << "output = " << output << std::endl;
//    std::cout << "---------------------------------------" << std::endl;
//
//    MeshFileReader reader(input.c_str());
//    Mesh& mesh = (Mesh&) reader.GetMesh();
//    mesh.RemoveUselessVertices();
//    mesh.BuildAllConnectivities();
//    mesh.ExtractBoundary();
//    mesh.ExtractSingularities();
//
//    BaseComplexQuad baseComplex(mesh);
//    baseComplex.Build();
//
//    BaseComplexSheetQuad baseComplexSheets(baseComplex);
//    baseComplexSheets.Extract();
//
//    auto dualMesh = Refine(mesh, 0);
//    auto key_edgeId = get_key_edgeId(mesh);
//    auto key_faceId = get_key_faceId(mesh);
//
//    std::map<size_t, size_t> canceledFaceIds;
//    std::set<size_t> canceledEdgeIds;
//    //for (auto& baseComplexSheetEdgeIds : baseComplexSheets.sheets_componentEdgeIds)
//    for (auto baseComplexEdgeId : baseComplexSheets.sheets_componentEdgeIds[17]) {
//        for (auto edgeId : baseComplex.componentE[baseComplexEdgeId].eids_link) {
//            auto& e = mesh.E[edgeId];
//            for (auto n_fid : e.N_Fids)
//                ++canceledFaceIds[n_fid];
//            canceledEdgeIds.insert(edgeId);
//        }
//    }
//    for (auto& item : canceledFaceIds) {
//        if (item.second >= 4) {
//            auto& f  = mesh.F.at(item.first);
//            auto key = get_facekey(f);
//            auto centerVid = key_faceId[key];
//            for (auto vid : f.Vids) {
//                if (vid >= mesh.V.size()) continue;
//                auto& v = mesh.V.at(vid);
//                for (auto n_fid : v.N_Fids) {
//                    auto& n_f = mesh.F.at(n_fid);
//                    for (auto& n_vid : n_f.Vids)
//                        if (n_vid == vid) n_vid = mesh.V.size() + mesh.E.size() + centerVid;
//                }
//            }
//            for (auto eid : f.Eids)
//                if (canceledEdgeIds.find(eid) != canceledEdgeIds.end()) canceledEdgeIds.erase(eid);
//        }
//    }
//    for (auto edgeId : canceledEdgeIds) {
//        auto& e = mesh.E[edgeId];
//        auto& v0 = mesh.V[e.Vids[0]];
//        auto& v1 = mesh.V[e.Vids[1]];
//        auto key = (e.Vids[0] << 32) | e.Vids[1];
//        if (key_edgeId.find(key) == key_edgeId.end()) std::cout << "Edge search Error !" << std::endl;
//        auto centerVid = key_edgeId[key];
//        for (auto vid : e.Vids) {
//            auto& v = mesh.V.at(vid);
//            for (auto n_fid : v.N_Fids) {
//                auto& n_f = mesh.F.at(n_fid);
//                for (auto& n_vid : n_f.Vids)
//                    if (n_vid == vid) n_vid = mesh.V.size() + centerVid;
//            }
//        }
//    }
//
//    std::vector<size_t> FaceIds;
//    FaceIds.reserve(mesh.F.size());
//    for (auto& f : mesh.F)
//        if (canceledFaceIds.find(f.id) == canceledFaceIds.end())
//            FaceIds.push_back(f.id);
//    MeshFileWriter writer(dualMesh.V, mesh.F, output.c_str());
//    writer.WriteFacesVtk(FaceIds);
//}

std::unordered_map<size_t, size_t> get_key_edgeId(const Mesh& mesh) {
    std::unordered_map<size_t, size_t> key_edgeId;
    for (size_t i = 0; i < mesh.E.size(); ++i) {
        const auto& e = mesh.E.at(i);
        key_edgeId[(e.Vids[0] << 32) | e.Vids[1]] = i;
        key_edgeId[(e.Vids[1] << 32) | e.Vids[0]] = i;
    }
    return key_edgeId;
}

std::unordered_map<std::string, size_t> get_key_faceId(const Mesh& mesh) {
    std::unordered_map<std::string, size_t> key_faceId;
    for (size_t i = 0; i < mesh.F.size(); ++i) {
        const auto& f = mesh.F.at(i);
        std::set<size_t> vids(f.Vids.begin(), f.Vids.end());
        std::string s;
        for (auto vid : vids)
            s += std::to_string(vid) + "@";
        key_faceId[s] = i;
    }
    return key_faceId;
}

const int QuadRefine[4][4] = {
    0, 4, 8, 7,
    1, 5, 8, 4,
    2, 6, 8, 5,
    3, 7, 8, 6
};

const int TriRefine[3][4] = {
    0, 3, 6, 5,
    1, 4, 6, 3,
    2, 5, 6, 4
};

Mesh Refine(const Mesh& hex_mesh, int clockwise) {
    const Mesh& new_mesh = hex_mesh;
    ////////////////////////////////////////////////////////////////////////////
    // add vertices
    std::vector<Vertex> new_vertex(new_mesh.V.size() + new_mesh.E.size() + new_mesh.C.size());
    for (size_t i = 0; i < new_mesh.V.size(); i++)
        new_vertex.at(i) = new_mesh.V.at(i);
    size_t offset = new_mesh.V.size();
    for (size_t i = 0; i < new_mesh.E.size(); i++) {
        const Edge& e = new_mesh.E.at(i);
        const Vertex& v0 = new_mesh.V.at(e.Vids[0]);
        const Vertex& v1 = new_mesh.V.at(e.Vids[1]);
        new_vertex.at(offset + i) = 0.5 * (v0.xyz() + v1.xyz());
    }
    offset = new_mesh.V.size() + new_mesh.E.size();
    size_t numOfTri = 0, numOfQuad = 0;
    for (size_t i = 0; i < new_mesh.C.size(); i++) {
        const auto& f = new_mesh.C.at(i);
        if (f.Vids.size() == 4) {
            const Vertex& v0 = new_mesh.V.at(f.Vids[0]);
            const Vertex& v1 = new_mesh.V.at(f.Vids[1]);
            const Vertex& v2 = new_mesh.V.at(f.Vids[2]);
            const Vertex& v3 = new_mesh.V.at(f.Vids[3]);
            new_vertex.at(offset + i) = 0.25 * (v0.xyz() + v1.xyz() + v2.xyz() + v3.xyz());
            ++numOfQuad;
        } else  if (f.Vids.size() == 3) {
            const Vertex& v0 = new_mesh.V.at(f.Vids[0]);
            const Vertex& v1 = new_mesh.V.at(f.Vids[1]);
            const Vertex& v2 = new_mesh.V.at(f.Vids[2]);
            new_vertex.at(offset + i) = 0.3333333 * (v0.xyz() + v1.xyz() + v2.xyz());
            ++numOfTri;
        }
    }
    auto key_edgeId = get_key_edgeId(new_mesh);
    //auto key_faceId = get_key_faceId(new_mesh);
    Cell cell(4);
    std::vector<Cell> new_cells(numOfTri * 3 + numOfQuad * 4, cell);
    int count = 0;
    for (size_t i = 0; i < new_mesh.C.size(); i++) {
        unsigned long v_index[9];
        const auto & f = new_mesh.C.at(i);
        for (auto j = 0; j < f.Vids.size(); j++)
            v_index[j] = f.Vids.at(j);
//        if (clockwise != 0) {
//            std::swap(v_index[1], v_index[3]);
//            std::swap(v_index[5], v_index[7]);
//        }
        if (f.Vids.size() == 4) {
            for (unsigned long j = 0; j < 4; j++) {
                const Edge e( { f.Vids.at(QuadEdge[j][0]), f.Vids.at(QuadEdge[j][1]) });
                unsigned long long key = (e.Vids[0] << 32) | e.Vids[1];
                if (key_edgeId.find(key) == key_edgeId.end()) std::cout << "Edge search Error !" << std::endl;
                auto e_index = key_edgeId[key];
                v_index[4 + j] = new_mesh.V.size() + e_index;
            }
            v_index[8] = new_mesh.V.size() + new_mesh.E.size() + i;
            for (int k = 0; k < 4; k++, count++)
                for (int j = 0; j < 4; j++)
                    new_cells[count].Vids[j] = v_index[QuadRefine[k][j]];
        } else if (f.Vids.size() == 3) {
            for (unsigned long j = 0; j < 3; j++) {
                const Edge e( { f.Vids.at(TriEdge[j][0]), f.Vids.at(TriEdge[j][1]) });
                unsigned long long key = (e.Vids[0] << 32) | e.Vids[1];
                if (key_edgeId.find(key) == key_edgeId.end()) std::cout << "Edge search Error !" << std::endl;
                auto e_index = key_edgeId[key];
                v_index[3 + j] = new_mesh.V.size() + e_index;
            }
            v_index[6] = new_mesh.V.size() + new_mesh.E.size() + i;
            for (int k = 0; k < 3; k++, count++)
                for (int j = 0; j < 4; j++)
                    new_cells[count].Vids[j] = v_index[TriRefine[k][j]];
        }
    }
    Mesh mesh(new_vertex, new_cells, QUAD);
    return mesh;
}
