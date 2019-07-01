/*
 * KillHighValenceSingularities.cpp
 *
 *  Created on: Nov 6, 2018
 *      Author: cotrik
 */
#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "ArgumentManager.h"

#include <unordered_map>
#include <unordered_set>
#include <set>
#include <iostream>
#include <set>
#include <algorithm>

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

void KillHighValenceSingularities(Mesh& mesh) {
	size_t numOfSingularV = 0;
	for (const auto& v : mesh.V)
		if (v.N_Fids.size() == 6) ++numOfSingularV;

	while (numOfSingularV--) {
		bool found_valence6 = false;
		for (const auto& v : mesh.V) {
			if (v.N_Fids.size() != 6) continue;
			auto centerV = v;
			found_valence6 = true;
			std::set<size_t> neighborFVids;
			for (auto nfid : v.N_Fids) {
				auto& nf = mesh.F.at(nfid);
				neighborFVids.insert(nf.Vids.begin(), nf.Vids.end());
			}
			neighborFVids.erase(v.id);
			std::vector<size_t> linkVids;
			for (auto nfvid : neighborFVids)
				if (mesh.V.at(nfvid).N_Fids.size() == 4) {
					linkVids.push_back(nfvid);
					neighborFVids.erase(nfvid);
					break;
				}

			auto start_vid = linkVids.back();
			while (!neighborFVids.empty()) {
				auto& start_v = mesh.V.at(start_vid);
				for (auto nvid : start_v.N_Vids) {
					if (neighborFVids.find(nvid) != neighborFVids.end()) {
						linkVids.push_back(nvid);
						neighborFVids.erase(nvid);
						start_vid = nvid;
						break;
					}
				}
			}

			Vertex new_v;
			std::vector<size_t> ids = { 2, 3, 4, 8, 9, 10 };
			auto xyz = centerV.xyz();
			for (auto id : ids) {
				new_v = 0.5 * (xyz + mesh.V.at(linkVids[id]).xyz());
				new_v.id = mesh.V.size();
				mesh.V.push_back(new_v);
				linkVids.push_back(new_v.id);
			}
			//// 12
			//new_v.xyz = 0.5 * (v + mesh.V.at(linkVids[2]));
			//new_v.id = mesh.V.size();
			//mesh.V.push_back(new_v);
			//linkVids.push_back(new_v.id);
			//// 13
			//new_v.xyz = 0.5 * (v + mesh.V.at(linkVids[3]));
			//new_v.id = mesh.V.size();
			//mesh.V.push_back(new_v);
			//linkVids.push_back(new_v.id);
			//// 14
			//new_v.xyz = 0.5 * (v + mesh.V.at(linkVids[4]));
			//new_v.id = mesh.V.size();
			//mesh.V.push_back(new_v);

			//// 15
			//new_v.xyz = 0.5 * (v + mesh.V.at(linkVids[8]));
			//new_v.id = mesh.V.size();
			//mesh.V.push_back(new_v);
			//linkVids.push_back(new_v.id);
			//// 16
			//new_v.xyz = 0.5 * (v + mesh.V.at(linkVids[9]));
			//new_v.id = mesh.V.size();
			//mesh.V.push_back(new_v);
			//linkVids.push_back(new_v.id);
			//// 17
			//new_v.xyz = 0.5 * (v + mesh.V.at(linkVids[10]));
			//new_v.id = mesh.V.size();
			//mesh.V.push_back(new_v);
			//linkVids.push_back(new_v.id);

			const std::vector<std::vector<size_t>> new_Fvids = {
				{ linkVids[0], linkVids[1], linkVids[2], linkVids[12] },
				{ linkVids[2], linkVids[3], linkVids[13], linkVids[12] },
				{ linkVids[3], linkVids[4], linkVids[14], linkVids[13] },
			    { linkVids[4], linkVids[5], linkVids[6], linkVids[14] },
				{ linkVids[6], linkVids[7], linkVids[8], linkVids[15] },
				{ linkVids[8], linkVids[9], linkVids[16], linkVids[15] },
				{ linkVids[9], linkVids[10], linkVids[17], linkVids[16] },
				{ linkVids[10], linkVids[11], linkVids[0], linkVids[17] },
				{ centerV.id, linkVids[17], linkVids[0], linkVids[12] },
				{ centerV.id, linkVids[12], linkVids[13], linkVids[14] },
				{ centerV.id, linkVids[14], linkVids[6], linkVids[15] },
				{ centerV.id, linkVids[15], linkVids[16], linkVids[17] }
			};
			Face new_f;
			for (auto& vids : new_Fvids) {
				new_f.Vids = vids;
				new_f.id = mesh.F.size();
				mesh.F.push_back(new_f);
			}

			std::set<size_t> canceledFids(centerV.N_Fids.begin(), centerV.N_Fids.end());
			update(mesh, canceledFids);
			mesh.BuildAllConnectivities();
			mesh.ExtractBoundary();
			mesh.ExtractSingularities();
			break;
		}
		if (!found_valence6) break;
	}
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: KillHighValenceSingularities quad.vtk out.quad.vtk" << std::endl;
        return -1;
    }
    ArgumentManager argumentManager(argc, argv);
    std::string input = argv[1];
    std::string output = argv[2];
    std::cout << "---------------------------------------" << std::endl;
    std::cout << "input  = " << input << std::endl;
    std::cout << "output = " << output << std::endl;
    std::cout << "---------------------------------------" << std::endl;

    MeshFileReader reader(input.c_str());
    Mesh& mesh = (Mesh&) reader.GetMesh();
    mesh.BuildAllConnectivities();
    mesh.ExtractBoundary();
    mesh.ExtractSingularities();
//    mesh.BuildParallelE();
//    mesh.BuildConsecutiveE();
//    mesh.BuildOrthogonalE();

    KillHighValenceSingularities(mesh);
    MeshFileWriter writer(mesh, output.c_str());
    writer.WriteFile();
}
