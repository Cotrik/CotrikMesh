/*
 * ExtractLocalRegion.cpp
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

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: ExtractLocalRegion quad.vtk local.quad.vtk restricted=<false>" << std::endl;
        return -1;
    }
    ArgumentManager argumentManager(argc, argv);
    std::string input = argv[1];
    std::string output = argv[2];
	bool restricted = false;
    int iters = 1;
    auto strIters = argumentManager.get("iters");
    if (!strIters.empty()) iters = std::stoi(strIters);
	auto strrestricted = argumentManager.get("restricted");
	if (!strrestricted.empty()) restricted = strrestricted != "false" ? true : false;
    std::cout << "---------------------------------------" << std::endl;
    std::cout << "input  = " << input << std::endl;
    std::cout << "output = " << output << std::endl;
    std::cout << "---------------------------------------" << std::endl;
    std::vector<std::set<size_t>> localFids;
    MeshFileReader reader(input.c_str());
    Mesh& mesh = (Mesh&) reader.GetMesh();
    mesh.CompressWithFeaturePreserved();
    mesh.BuildAllConnectivities();
    mesh.ExtractBoundary();
    mesh.ExtractSingularities();
    get_feature(mesh);
    BaseComplexQuad baseComplex(mesh);
    baseComplex.Build();
    for (const auto& link : baseComplex.separatedVertexIdsLink) {
        auto& v_front = mesh.V.at(link.front());
        auto& v_back = mesh.V.at(link.back());
        if ((v_front.N_Fids.size() <= 5 && v_back.N_Fids.size() <= 5) && v_front.N_Fids.size() != v_back.N_Fids.size()) {
            ; //ofs << 0 << std::endl;
        } else if (v_front.N_Fids.size() == 3 && v_back.N_Fids.size() == 3) {
            std::set<size_t> canceledFids;
			if (!restricted) {
				for (auto vid : link) {
					auto& v = mesh.V.at(vid);
					canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
				}
			} else {
				auto v_front_fid = get_faceid(mesh, v_front.id, link[1]);
				auto v_front_fvid = get_diagnal_vid(mesh, v_front.id, v_front_fid);

				auto v_back_fid = get_faceid(mesh, v_back.id, link[link.size() - 2]);
				auto v_back_fvid = get_diagnal_vid(mesh, v_back.id, v_back_fid);

				if (mesh.V.at(v_front_fvid).N_Fids.size() == 5 && mesh.V.at(v_back_fvid).N_Fids.size() == 5) {
					for (auto vid : link) {
						auto& v = mesh.V.at(vid);
						canceledFids.insert(v.N_Fids.begin(), v.N_Fids.end());
					}
				}
			}
            localFids.push_back(canceledFids);
        } else if (v_front.N_Fids.size() == 5 && v_back.N_Fids.size() == 5) {
            ; //ofs << 2 << std::endl;
        } else {
            ; //ofs << 3 << std::endl;
        }
    }

    MeshFileWriter writer(mesh, output.c_str());
    writer.WriteFacesVtk(localFids);
    return 0;
}

