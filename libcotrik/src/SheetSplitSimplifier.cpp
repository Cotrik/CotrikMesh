/*
 * SheetSimplifier.cpp
 *
 *  Created on: Jan 20, 2019
 *      Author: cotrik
 */

#include "SheetSplitSimplifier.h"

SheetSplitSimplifier::SheetSplitSimplifier(Mesh& mesh) : Simplifier(mesh) {
    // TODO Auto-generated constructor stub

}

SheetSplitSimplifier::~SheetSplitSimplifier() {
    // TODO Auto-generated destructor stub
}

void SheetSplitSimplifier::Run(std::set<size_t>& canceledFids) {
    for (auto& v : mesh.V) {
        // if (v.N_Fids.size() != 2 || v.isBoundary) continue;
        if ((v.N_Fids.size() == 2 && !v.isBoundary) || (v.N_Fids.size() == 1 && v.type == FEATURE))
        {
            auto paralell_eids = GetAllParallelEdgeIdsFromDoubletVertex(v);
            canceledFids = GetCanceledFids(paralell_eids);
            Split(v, paralell_eids, canceledFids);
            return;
        }
    }
}

void SheetSplitSimplifier::Split(const Vertex& v, const std::set<size_t>& paralell_eids, std::set<size_t>& canceledFids) {
    auto overlapFids = GetOverlapFids(paralell_eids, canceledFids);
    auto overlapFid_newVid = GetInsertFaceVertices(overlapFids);

    auto paralellEid_newVid = GetInsertEdgeVertices(paralell_eids);
    auto doubletFid_newVid = GetInsertFaceVertices(v);

    auto nonoverlapFids = canceledFids;
    for (auto fid : v.N_Fids)
        nonoverlapFids.erase(fid);
    for (auto fid : overlapFids)
        nonoverlapFids.erase(fid);

    Insert2Faces(paralellEid_newVid, nonoverlapFids);                  // insert 2 faces for nonoverlap faces
    Insert3Faces(paralellEid_newVid, doubletFid_newVid, v);            // insert 3 faces for doublet faces
    Insert4Faces(paralellEid_newVid, overlapFid_newVid, overlapFids);  // insert 4 faces for overlap faces
}

std::map<size_t, size_t> SheetSplitSimplifier::GetInsertEdgeVertices(const std::set<size_t>& paralell_eids) {
    auto orig_patch_id = MAXID;
    for (auto& v : origMesh.V) {
        if (!v.isBoundary && v.patch_id != MAXID) {
            orig_patch_id = v.patch_id;
            break;
        }
    }
    std::map<size_t, size_t> res;
    for (auto eid : paralell_eids) {
        auto& e = mesh.E.at(eid);
        auto& v0 = mesh.V.at(e.Vids[0]);
        auto& v1 = mesh.V.at(e.Vids[1]);
        Vertex v = 0.5 * (v0 + v1);
        v.id = mesh.V.size();
        v.patch_id = GetPatchid(v0, v1);
        if (v.patch_id == MAXID || v.patch_ids.empty()) v.patch_ids.insert(orig_patch_id);
		for (auto fvid : e.Vids) {
			auto& fv = mesh.V.at(fvid);
			v.patch_ids.insert(fv.patch_ids.begin(), fv.patch_ids.end());
		}
        v.label = GetLabel(v0, v1, v);
        mesh.V.push_back(v);
        res[eid] = v.id;
    }
    return res;
}

std::map<size_t, size_t> SheetSplitSimplifier::GetInsertFaceVertices(const std::set<size_t>& overlapFids) {
        auto orig_patch_id = MAXID;
        for (auto& v : origMesh.V) {
            if (v.patch_id != MAXID) {
                orig_patch_id = v.patch_id;
                break;
            }
        }
    std::map<size_t, size_t> res;
    for (auto fid : overlapFids) {
        auto& f = mesh.F.at(fid);
        auto& v0 = mesh.V.at(f.Vids[0]);
        auto& v1 = mesh.V.at(f.Vids[1]);
        auto& v2 = mesh.V.at(f.Vids[2]);
        auto& v3 = mesh.V.at(f.Vids[3]);
        Vertex v = 0.25 * (v0 + v1 + v2 + v3);
        v.id = mesh.V.size();
        v.patch_id = GetPatchid(v0, v1, v2, v3);
        if (v.patch_id == MAXID) v.patch_id = orig_patch_id;
		for (auto fvid : f.Vids) {
			auto& fv = mesh.V.at(fvid);
			v.patch_ids.insert(fv.patch_ids.begin(), fv.patch_ids.end());
		}
        mesh.V.push_back(v);
        res[fid] = v.id;
    }
    return res;
}

std::map<size_t, size_t> SheetSplitSimplifier::GetInsertFaceVertices(const Vertex& v) {
    std::set<size_t> fids(v.N_Fids.begin(), v.N_Fids.end());
    return GetInsertFaceVertices(fids);
}

size_t SheetSplitSimplifier::find(const std::map<size_t, size_t>& m, const size_t key, const char* mapName/* = NULL*/) {
    auto iter = m.find(key);
    if (iter == m.end()) {
        std::cerr << " Err in " << mapName << ".find(" << key << ")\n";
        MeshFileWriter writer(mesh, "E.vtk");
        writer.WriteEdgesVtk();
        return MAXID;
    }
    return iter->second;
}

void SheetSplitSimplifier::Insert2Faces(const std::map<size_t, size_t>& paralellEid_newVid, const std::set<size_t>& fids) {
    for (auto fid : fids) {
        auto& f = mesh.F.at(fid);
        auto paralellEids = GetParallelEids(f, paralellEid_newVid);
        auto vid0 = find(paralellEid_newVid, paralellEids.front(), "paralellEid_newVid");
        auto vid3 = find(paralellEid_newVid, paralellEids.back(), "paralellEid_newVid");
        auto& e0 = mesh.E.at(paralellEids.front());
        auto& e3 = mesh.E.at(paralellEids.back());
        auto perpendicularEids = GetPerpendicularEids(f, paralellEid_newVid);
        for (auto eid : perpendicularEids) {
            auto& e = mesh.E.at(eid);
            auto vid1 = GetSharedVid(e, e0);
            auto vid2 = GetSharedVid(e, e3);
            Face newf;
            newf.id = mesh.F.size();
            newf.Vids = {vid0, vid1, vid2, vid3};
            mesh.F.push_back(newf);
        }
    }
}

void SheetSplitSimplifier::Insert3Faces(const std::map<size_t, size_t>& paralellEid_newVid,
        const std::map<size_t, size_t>& doubletFid_newVid, const Vertex& v) {
    for (auto fid : v.N_Fids) {
        auto diagnalVid = get_diagnal_vid(v.id, fid);
        auto vid0 = find(doubletFid_newVid, fid, "doubletFid_newVid");
        auto& v0 = mesh.V.at(vid0);
        auto& f = mesh.F.at(fid);
        for (auto vid : f.Vids) {
            if (vid == v.id) continue;
            if (vid == diagnalVid) {
                auto& v2 = mesh.V.at(vid);
                auto eids = GetNeignborEids(v2, f);
                if (eids.size() != 2) std::cerr << " Err in eids.size() v = " << vid << "\n";
                auto vid1 = find(paralellEid_newVid, eids.front(), "paralellEid_newVid");
                auto vid3 = find(paralellEid_newVid, eids.back(), "paralellEid_newVid");
                auto& v1 = mesh.V.at(vid1);
                auto& v3 = mesh.V.at(vid3);
                Face newf;
                newf.id = mesh.F.size();
                newf.Vids = {v0.id, v1.id, v2.id, v3.id};
                mesh.F.push_back(newf);
            } else {
                auto& v2 = mesh.V.at(vid);
                auto eids = GetNeignborEids(v2, f);
                if (eids.size() != 2) std::cerr << " Err in eids.size() v = " << vid << "\n";
                std::set<size_t> eids_set(eids.begin(), eids.end());
                std::set<size_t> veids_set(v.N_Eids.begin(), v.N_Eids.end());
                auto intersections = Util::get_intersect(eids_set, veids_set);
                if (intersections.size() != 1) std::cerr << " Err in intersections.size() v = " << vid << "\n";
                auto eid = *intersections.begin() == eids.front() ? eids.back() : eids.front();

                auto vid1 = find(paralellEid_newVid, eid, "paralellEid_newVid");
                auto vid3 = v.id;
                auto& v1 = mesh.V.at(vid1);
                auto& v3 = mesh.V.at(vid3);
                Face newf;
                newf.id = mesh.F.size();
                newf.Vids = {v0.id, v1.id, v2.id, v3.id};
                mesh.F.push_back(newf);
            }
        }
    }
}

void SheetSplitSimplifier::Insert4Faces(const std::map<size_t, size_t>& paralellEid_newVid,
        const std::map<size_t, size_t>& overlapFid_newVid, const std::set<size_t>& fids) {
    for (auto fid : fids) {
        auto vid0 = find(overlapFid_newVid, fid, "overlapFid_newVid");
        auto& v0 = mesh.V.at(vid0);
        auto& f = mesh.F.at(fid);
        for (auto vid : f.Vids) {
            auto& v2 = mesh.V.at(vid);
            auto eids = GetNeignborEids(v2, f);
            if (eids.size() != 2) std::cerr << " Err in eids.size() v = " << vid << "\n";
            auto vid1 = find(paralellEid_newVid, eids.front(), "paralellEid_newVid");
            auto vid3 = find(paralellEid_newVid, eids.back(), "paralellEid_newVid");
            auto& v1 = mesh.V.at(vid1);
            auto& v3 = mesh.V.at(vid3);
            Face newf;
            newf.id = mesh.F.size();
            newf.Vids = {v0.id, v1.id, v2.id, v3.id};
            mesh.F.push_back(newf);
        }
    }
}

std::vector<size_t> SheetSplitSimplifier::GetNeignborEids(const Vertex& v, const Face& f) {
    std::vector<size_t> res;
    std::set<size_t> eids(f.Eids.begin(), f.Eids.end());
    for (auto eid : v.N_Eids)
        if (eids.find(eid) != eids.end()) res.push_back(eid);
    return res;
}

std::vector<size_t> SheetSplitSimplifier::GetParallelEids(const Face& f, const std::map<size_t, size_t>& paralellEid_newVid) {
    std::vector<size_t> res;
    for (auto eid : f.Eids)
        if (paralellEid_newVid.find(eid) != paralellEid_newVid.end()) res.push_back(eid);
    return res;
}

std::vector<size_t> SheetSplitSimplifier::GetPerpendicularEids(const Face& f, const std::map<size_t, size_t>& paralellEid_newVid) {
    std::vector<size_t> res;
    for (auto eid : f.Eids)
        if (paralellEid_newVid.find(eid) == paralellEid_newVid.end()) res.push_back(eid);
    return res;
}

size_t SheetSplitSimplifier::GetPatchid(const Vertex& v0, const Vertex& v1, const Vertex& v2, const Vertex& v3) {
    auto patch_id = MAXID;
    if (v0.isBoundary && v1.isBoundary && v2.isBoundary && v3.isBoundary) {
        auto intersections1 = Util::get_intersect(v0.patch_ids, v1.patch_ids);
        auto intersections2 = Util::get_intersect(v2.patch_ids, v3.patch_ids);
        auto intersections = Util::get_intersect(intersections1, intersections2);
		if (intersections.size() != 1) {
			std::cerr << "Err in SheetSplitSimplifier::GetPatchid!\n";
			{
				static int num = 0;
				std::string fname = std::string("ErrGetPatchid") + std::to_string(num++) + ".vtk";
				std::cout << "writing " << fname << std::endl;
				MeshFileWriter writer(mesh, fname.c_str());
				writer.WriteFile();
			}
		}
        else patch_id = *intersections.begin();
    } else if (!v0.isBoundary) patch_id = v0.patch_id;
    else if (!v1.isBoundary) patch_id = v1.patch_id;
    else if (!v2.isBoundary) patch_id = v2.patch_id;
    else if (!v3.isBoundary) patch_id = v3.patch_id;
	else {
		std::cerr << "Err in SheetSplitSimplifier::GetPatchid!\n";
		{
			static int num = 0;
			std::string fname = std::string("ErrGetPatchid_") + std::to_string(num++) + ".vtk";
			std::cout << "writing " << fname << std::endl;
			MeshFileWriter writer(mesh, fname.c_str());
			writer.WriteFile();
		}
	}
    return patch_id;
}

size_t SheetSplitSimplifier::GetPatchid(const Vertex& v0, const Vertex& v1) {
    auto patch_id = MAXID;
    if (v0.isBoundary && v1.isBoundary) {
        auto intersections = Util::get_intersect(v0.patch_ids, v1.patch_ids);
		if (intersections.size() != 1) {
			std::cerr << "Err in SheetSplitSimplifier::GetPatchid!\n";
			{
				static int num = 0;
				std::string fname = std::string("ErrGetPatchid__") + std::to_string(num++) + ".vtk";
				std::cout << "writing " << fname << std::endl;
				MeshFileWriter writer(mesh, fname.c_str());
				writer.WriteFile();
			}
		}
        else patch_id = *intersections.begin();
    } else patch_id = v0.isBoundary ? v1.patch_id : v0.patch_id;
    return patch_id;
}

//size_t SheetSplitSimplifier::GetPatchid(const Vertex& v0, const Vertex& v1, Vertex& v) {
//    auto orig_patch_id = MAXID;
//    for (auto& v : origMesh.V) {
//        if (v.patch_id != MAXID) {
//            orig_patch_id = v.patch_id;
//            break;
//        }
//    }
//    auto patch_id = MAXID;
//    if (v0.isBoundary && v1.isBoundary) {
//        if (v0.type == FEATURE && v1.type == FEATURE) {
//            v.patch_ids.insert(orig_patch_id);
//        } else if (v0.type == CORNER && v1.type == FEATURE) {
//            ;
//        } else if (v1.type == CORNER && v0.type == FEATURE && v1.labels.find(v0.label) != v1.labels.end()) {
//            ;
//        } else if (v1.type == CORNER && v0.type == CORNER) {
//            auto intersections = Util::get_intersect(v0.labels, v1.labels);
//            if (intersections.size() == 1) {
//                label = *intersections.begin();
//                v.label = label;
//                v.type = FEATURE;
//            }
//        } else if (v0.type == FEATURE && v1.type == FEATURE && v0.label != v1.label) {
//            ;
//        }
//    }
//    return label;
//}

size_t SheetSplitSimplifier::GetLabel(const Vertex& v0, const Vertex& v1, Vertex& v) {
    auto label = MAXID;
    if (v0.isBoundary && v1.isBoundary) {
        if (v0.type == FEATURE && v1.type == FEATURE && v0.label == v1.label) {
            label = v0.label;
            v.label = label;
            v.type = FEATURE;
        } else if (v0.type == CORNER && v1.type == FEATURE && v0.labels.find(v1.label) != v0.labels.end()) {
            label = v1.label;
            v.label = label;
            v.type = FEATURE;
        } else if (v1.type == CORNER && v0.type == FEATURE && v1.labels.find(v0.label) != v1.labels.end()) {
            label = v0.label;
            v.label = label;
            v.type = FEATURE;
        } else if (v1.type == CORNER && v0.type == CORNER) {
            auto intersections = Util::get_intersect(v0.labels, v1.labels);
            if (intersections.size() == 1) {
                label = *intersections.begin();
                v.label = label;
                v.type = FEATURE;
            }
        } else if (v0.type == FEATURE && v1.type == FEATURE && v0.label != v1.label) {
            ;
        }
    }
    return label;
}

size_t SheetSplitSimplifier::GetSharedVid(const Edge& e0, const Edge& e1) {
    size_t res = MAXID;
    if (e0.Vids[0] == e1.Vids[0] || e0.Vids[0] == e1.Vids[1]) res = e0.Vids[0];
    else if (e0.Vids[1] == e1.Vids[0] || e0.Vids[1] == e1.Vids[1]) res = e0.Vids[1];
    if (res == MAXID)
        std::cerr << "Err in SheetSplitSimplifier::GetSharedVid()\n";
    return res;
}

std::set<size_t> SheetSplitSimplifier::GetAllParallelEdgeIdsFromDoubletVertex(const Vertex& v) {
    auto& v0 = mesh.V.at(v.N_Vids[0]);
    auto& v1 = mesh.V.at(v.N_Vids[1]);
    auto& f0 = mesh.F.at(v.N_Fids.front());
    auto& f1 = mesh.F.at(v.N_Fids.back());
    std::set<size_t> eids_set(f0.Eids.begin(), f0.Eids.end());
    eids_set.insert(f1.Eids.begin(), f1.Eids.end());
    for (auto eid : v.N_Eids)
        eids_set.erase(eid);
    std::set<size_t> exclude_eids(v.N_Eids.begin(), v.N_Eids.end());
    std::set<size_t> totalParalellEids;
    for (auto eid : eids_set) {
        auto paralellEids = GetAllParallelEdgeIds(eid, exclude_eids);
        totalParalellEids.insert(paralellEids.begin(), paralellEids.end());
    }
    return totalParalellEids;
}

std::set<size_t> SheetSplitSimplifier::GetAllParallelEdgeIds(const size_t eid, const std::set<size_t>& exclude_eids) {
    std::set<size_t> res;
    res.insert(eid);
    std::queue<size_t> q;
    q.push(eid);
    while (!q.empty()) {
        auto n = q.size();
        for (auto i = 0; i < n; ++i) {
            auto t = q.front();
            q.pop();
            auto& e = mesh.E.at(t);
            for (auto peid : e.parallelEids)
                if (res.find(peid) == res.end() && exclude_eids.find(peid) == exclude_eids.end()) {
                    res.insert(peid);
                    q.push(peid);
                }
        }
    }
    return res;
}

std::set<size_t> SheetSplitSimplifier::GetCanceledFids(const std::set<size_t>& paralell_eids) {
    std::set<size_t> canceledFids;
    for (auto eid : paralell_eids) {
        auto& e = mesh.E.at(eid);
        canceledFids.insert(e.N_Fids.begin(), e.N_Fids.end());
    }
    return canceledFids;
}

std::set<size_t> SheetSplitSimplifier::GetOverlapFids(const std::set<size_t>& paralell_eids, const std::set<size_t>& canceledFids) {
    std::set<size_t> overlapFids;
    for (auto fid : canceledFids) {
        auto& f = mesh.F.at(fid);
        bool foundAll = true;
        for (auto eid : f.Eids)
            if (paralell_eids.find(eid) == paralell_eids.end()) {
                foundAll = false;
                break;
            }
        if (foundAll) overlapFids.insert(fid);
    }
    return overlapFids;
}
