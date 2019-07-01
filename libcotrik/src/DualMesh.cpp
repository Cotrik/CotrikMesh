#include "DualMesh.h"
#include "MeshFileWriter.h"
#include <iostream>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>

DualMesh::DualMesh() {}


DualMesh::~DualMesh() {}

std::vector<size_t> get_link_fids(const Mesh& mesh, const Vertex& v, std::vector<size_t>& ringFids) {
	auto fids = v.N_Fids;
	fids.pop_back();
	while (ringFids.size() < v.N_Fids.size()) {
		auto last_fid = ringFids.back();
		int j = 0;
		for (auto next_fid : fids) {
			auto& last_f = mesh.F.at(last_fid);
			std::set<size_t> vids(last_f.Vids.begin(), last_f.Vids.end());
			auto& next_f = mesh.F.at(next_fid);
			vids.insert(next_f.Vids.begin(), next_f.Vids.end());
			if (vids.size() == last_f.Vids.size() + next_f.Vids.size() - 2) {
				ringFids.push_back(next_fid);
				fids.erase(fids.begin() + j);
				break;
			}
			++j;
		}
	}
	return ringFids;
}

std::vector<size_t> get_link_fids1(const Mesh& mesh, const Vertex& v, std::vector<size_t>& ringFids) {
	auto fids = v.N_Fids;
	for (auto iter = fids.begin(); iter != fids.end(); ++iter)
		if (*iter == ringFids.back()) {
			fids.erase(iter);
			break;
		}

	while (ringFids.size() < v.N_Fids.size()) {
		auto last_fid = ringFids.back();
		int j = 0;
		for (auto next_fid : fids) {
			auto& last_f = mesh.F.at(last_fid);
			std::set<size_t> vids(last_f.Vids.begin(), last_f.Vids.end());
			auto& next_f = mesh.F.at(next_fid);
			vids.insert(next_f.Vids.begin(), next_f.Vids.end());
			if (vids.size() == last_f.Vids.size() + next_f.Vids.size() - 2) {
				ringFids.push_back(next_fid);
				fids.erase(fids.begin() + j);
				break;
			}
			++j;
		}
	}
	return ringFids;
}


std::vector<size_t> get_link_fids_boundary(const Mesh& mesh, const Vertex& v) {
	std::vector<size_t> ringFids;
	std::vector<size_t> N_Fids_boundary;
	for (auto nfid : v.N_Fids)
		if (mesh.F.at(nfid).isBoundary) N_Fids_boundary.push_back(nfid);

	std::vector<size_t> fids = N_Fids_boundary;
	ringFids.push_back(N_Fids_boundary.back());
	fids.pop_back();
	while (ringFids.size() < N_Fids_boundary.size()) {
		auto last_fid = ringFids.back();
		int j = 0;
		for (auto next_fid : fids) {
			auto& last_f = mesh.F.at(last_fid);
			std::set<size_t> vids(last_f.Vids.begin(), last_f.Vids.end());
			auto& next_f = mesh.F.at(next_fid);
			vids.insert(next_f.Vids.begin(), next_f.Vids.end());
			if (vids.size() == last_f.Vids.size() + next_f.Vids.size() - 2) {
				ringFids.push_back(next_fid);
				fids.erase(fids.begin() + j);
				break;
			}
			++j;
		}
	}
	return ringFids;
}

std::vector<size_t> get_link_fids(const Mesh& mesh, const Vertex& v) {
	if (v.isBoundary) return {};
	std::vector<size_t> ringFids;
	ringFids.push_back(v.N_Fids.back());
	return get_link_fids(mesh, v, ringFids);
}

std::vector<size_t> get_link_fids(const Mesh& mesh, const Vertex& v, std::map<size_t, size_t>& boundaryEid_dualVid) {
	std::vector<size_t> ringFids;
	std::vector<size_t> boundaryEids;
	for (auto eid : v.N_Eids) {
		auto& e = mesh.E.at(eid);
		if (e.isBoundary) boundaryEids.push_back(e.id);
	}
	ringFids.push_back(mesh.E.at(boundaryEids.front()).N_Fids.front());
	get_link_fids1(mesh, v, ringFids);
	ringFids.push_back(boundaryEid_dualVid[boundaryEids.back()]);
	ringFids.push_back(boundaryEid_dualVid[boundaryEids.front()]);
	return ringFids;
}

std::vector<size_t> get_link_cids(const Mesh& mesh, const Edge& e) {
	if (e.isBoundary) return {};
	auto cids = e.N_Cids;
	std::vector<size_t> ringCids;
	ringCids.push_back(cids.back());
	cids.pop_back();
	while (ringCids.size() < e.N_Cids.size()) {
		auto last_cid = ringCids.back();
		int j = 0;
		for (auto next_cid : cids) {
			auto& last_c = mesh.C.at(last_cid);
			std::set<size_t> vids(last_c.Vids.begin(), last_c.Vids.end());
			auto& next_c = mesh.C.at(next_cid);
			vids.insert(next_c.Vids.begin(), next_c.Vids.end());
			if (vids.size() == last_c.Vids.size() + next_c.Vids.size() - 4) {
				ringCids.push_back(next_cid);
				cids.erase(cids.begin() + j);
				break;
			}
			++j;
		}
	}
	return ringCids;
}

std::vector<size_t> get_link_cids(const Mesh& mesh, const Edge& e, std::vector<size_t>& ringCids) {
	if (!e.isBoundary) return {};
	auto cids = e.N_Cids;
	for (auto iter = cids.begin(); iter != cids.end(); ++iter)
		if (*iter == ringCids.back()) {
			cids.erase(iter);
			break;
		}
	while (ringCids.size() < e.N_Cids.size()) {
		auto last_cid = ringCids.back();
		int j = 0;
		for (auto next_cid : cids) {
			auto& last_c = mesh.C.at(last_cid);
			std::set<size_t> vids(last_c.Vids.begin(), last_c.Vids.end());
			auto& next_c = mesh.C.at(next_cid);
			vids.insert(next_c.Vids.begin(), next_c.Vids.end());
			if (vids.size() == last_c.Vids.size() + next_c.Vids.size() - 4) {
				ringCids.push_back(next_cid);
				cids.erase(cids.begin() + j);
				break;
			}
			++j;
		}
	}
	return ringCids;
}

std::vector<size_t> get_link_cids(const Mesh& mesh, const Edge& e, std::map<size_t, size_t>& boundaryFid_dualVid) {
	std::vector<size_t> ringCids;
	std::vector<size_t> boundaryFids;
	for (auto fid : e.N_Fids) {
		auto& f = mesh.F.at(fid);
		if (f.isBoundary) boundaryFids.push_back(f.id);
	}
	ringCids.push_back(mesh.F.at(boundaryFids.front()).N_Cids.front());
	get_link_cids(mesh, e, ringCids);
	ringCids.push_back(boundaryFid_dualVid[boundaryFids.back()]);
	ringCids.push_back(boundaryFid_dualVid[boundaryFids.front()]);
	return ringCids;
}

void DualMesh::Build(const Mesh& mesh) {
	if (mesh.m_cellType == TRIANGLE || mesh.m_cellType == QUAD || mesh.m_cellType == POLYGON) {
		m_cellType = POLYGON;
		if (mesh.hasBoundary) {
			size_t numOfBoundaryE = 0;
			for (auto& e : mesh.E)
				if (e.isBoundary) ++numOfBoundaryE;
			V.resize(mesh.F.size() + numOfBoundaryE);
			BuildBoundaryV(mesh);
		} else V.resize(mesh.F.size());
		BuildV(mesh);

		C.resize(mesh.V.size());
		BuildC(mesh);
	} else if (mesh.m_cellType == TETRAHEDRA || mesh.m_cellType == HEXAHEDRA) {
		m_cellType = POLYHEDRA;
		size_t numOfBoundaryF = 0;	
		for (auto& f : mesh.F) if (f.isBoundary) ++numOfBoundaryF;
		std::map<size_t, size_t> boundaryFid_dualVid;
		std::map<size_t, size_t> vid_dualFid;
		V.resize(mesh.C.size() + numOfBoundaryF);
		BuildV(mesh);
		BuildF(mesh, boundaryFid_dualVid, vid_dualFid);
		C.reserve(mesh.V.size());
		BuildC(mesh, boundaryFid_dualVid, vid_dualFid);
		MeshFileWriter writer(V, F, "F.vtk", POLYGON);
		writer.WriteFacesVtk();
	}
	//} else if (mesh.m_cellType == TETRAHEDRA || mesh.m_cellType == HEXAHEDRA) {
	//	m_cellType = POLYHEDRA;
	//	V.resize(mesh.C.size());
	//	for (auto i = 0; i < mesh.C.size(); ++i) {
	//		auto& v = V.at(i);
	//		auto& c = mesh.C.at(i);
	//		v.id = i;
	//		for (auto vid : c.Vids) v += mesh.V.at(vid).xyz();
	//		v /= c.Vids.size();
	//	}

	//	C.reserve(mesh.V.size());
	//	Cell c;
	//	auto id = 0;
	//	for (auto i = 0; i < mesh.V.size(); ++i) {
	//		auto& v = mesh.V.at(i);
	//		if (v.isBoundary) continue;
	//		//if (v.isSingularity) continue;
	//		auto eid0 = v.N_Eids.front();
	//		bool found_singularity = false;
	//		for (auto eid : v.N_Eids) {
	//			auto& e = mesh.E.at(eid);
	//			if (e.isSingularity) {
	//				found_singularity = true;
	//				eid0 = eid;
	//				break;
	//			}
	//		}
	//		//if (found_singularity) continue;
	//		auto& e0 = mesh.E.at(eid0);
	//		c.Vids = e0.N_Cids;
	//		auto cids = e0.N_Cids;
	//		std::vector<size_t> ringCids;
	//		ringCids.push_back(cids.back());
	//		cids.pop_back();
	//		while (ringCids.size() < c.Vids.size()) {
	//			auto last_cid = ringCids.back();
	//			int j = 0;
	//			for (auto next_cid : cids) {
	//				auto& last_c = mesh.C.at(last_cid);
	//				std::set<size_t> vids(last_c.Vids.begin(), last_c.Vids.end());
	//				auto& next_c = mesh.C.at(next_cid);
	//				vids.insert(next_c.Vids.begin(), next_c.Vids.end());
	//				if (vids.size() == last_c.Vids.size() + next_c.Vids.size() - 4) {
	//					ringCids.push_back(next_cid);
	//					cids.erase(cids.begin() + j);
	//					break;
	//				}
	//				++j;
	//			}
	//		}
	//		std::set<size_t> ringCids_set(ringCids.begin(), ringCids.end());
	//		//for (auto cid : v.N_Cids) {
	//		//	if (ringCids_set.find(cid) != ringCids_set.end()) continue;
	//		//	auto& cc = mesh.C.at(cid);
	//		//	for (auto ncid : ringCids_set) {
	//		//		std::set<size_t> vids(cc.Vids.begin(), cc.Vids.end());
	//		//		auto& nc = mesh.C.at(ncid);
	//		//		vids.insert(nc.Vids.begin(), nc.Vids.end());
	//		//		if (vids.size() == cc.Vids.size() + nc.Vids.size() - 4) {
	//		//			ringCids.push_back(cid);
	//		//			break;
	//		//		}
	//		//	}
	//		//}

	//		auto x = ringCids;
	//		for (auto ncid : x) {
	//			auto& nc = mesh.C.at(ncid);
	//			for (auto cid : v.N_Cids) {
	//				if (ringCids_set.find(cid) != ringCids_set.end()) continue;
	//				auto& cc = mesh.C.at(cid);
	//				std::set<size_t> vids(cc.Vids.begin(), cc.Vids.end());
	//				vids.insert(nc.Vids.begin(), nc.Vids.end());
	//				if (vids.size() == cc.Vids.size() + nc.Vids.size() - 4) {
	//					ringCids.push_back(cid);
	//					break;
	//				}
	//			}
	//		}

	//		c.Vids = ringCids; 
	//		c.id = id++;
	//		if (c.Vids.size() == 6) c.cellType = VTK_WEDGE;
	//		else if (c.Vids.size() == 8) c.cellType = VTK_HEXAHEDRON;
	//		else if (c.Vids.size() == 10) c.cellType = VTK_PENTAGONAL_PRISM;
	//		C.push_back(c);
	//	}
	//	m_cellTypes.resize(C.size());
	//	for (auto& c : C) m_cellTypes[c.id] = c.cellType;
	//}
}

void DualMesh::BuildV(const Mesh& mesh) {
	for (auto i = 0; i < mesh.C.size(); ++i) {
		auto& v = V.at(i);
		auto& c = mesh.C.at(i);
		v.id = i;
		for (auto vid : c.Vids) v += mesh.V.at(vid).xyz();
		v /= c.Vids.size();
	}
}

void DualMesh::BuildBoundaryV(const Mesh& mesh) {
	size_t id = 0;
	for (auto& e : mesh.E) {
		if (!e.isBoundary) continue;
		auto& v = V.at(mesh.F.size() + id);
		v.id = mesh.F.size() + id++;
		for (auto vid : e.Vids) v += mesh.V.at(vid).xyz();
		v /= e.Vids.size();
	}
}

void DualMesh::BuildF(const Mesh& mesh, std::map<size_t, size_t>& boundaryFid_dualVid, std::map<size_t, size_t>& vid_dualFid) {
	size_t numOfBoundaryF = 0;
	//std::map<size_t, size_t> boundaryFid_dualVid;
	for (auto& f : mesh.F)
		if (f.isBoundary) boundaryFid_dualVid[f.id] = mesh.C.size() + numOfBoundaryF++;
	std::map<size_t, std::vector<size_t>> eid_dualCids;
	for (auto& e : mesh.E) {
		if (!e.isBoundary) eid_dualCids[e.id] = get_link_cids(mesh, e);
		else eid_dualCids[e.id] = get_link_cids(mesh, e, boundaryFid_dualVid);
	}
	size_t numOfBoundaryV = 0;
	std::map<size_t, std::vector<size_t>> vid_dualFids;
	//std::map<size_t, size_t> vid_dualFid;
	for (auto& v : mesh.V)
		if (v.isBoundary) {
			vid_dualFids[v.id] = get_link_fids_boundary(mesh, v);
			vid_dualFid[v.id] = mesh.E.size() + numOfBoundaryV++;
		}
	size_t id = 0;
	for (auto& f : mesh.F) {
		if (!f.isBoundary) continue;
		auto& v = V.at(mesh.C.size() + id);
		v.id = mesh.C.size() + id++;
		for (auto vid : f.Vids) v += mesh.V.at(vid).xyz();
		v /= f.Vids.size();
	}
	F.reserve(mesh.E.size() + numOfBoundaryV);
	Face f;
	for (auto& item : eid_dualCids) {
		f.id = item.first;
		f.Vids = item.second;
		F.push_back(f);
	}
	id = mesh.E.size();
	for (auto& item : vid_dualFids) {
		f.id = id++;
		f.Vids = item.second;
		for (auto& vid : f.Vids)
			vid = boundaryFid_dualVid[vid];
		F.push_back(f);
	}
}

void DualMesh::BuildC(const Mesh& mesh) {
	if (mesh.m_cellType == TRIANGLE || mesh.m_cellType == QUAD || mesh.m_cellType == POLYGON) {
/*		for (auto i = 0; i < mesh.V.size(); ++i) {
			auto& v = mesh.V.at(i);
			auto& c = C.at(i);
			c.id = i;
			c.Vids = get_link_fids(mesh, v);
		}*/			
		size_t numOfBoundaryE = 0;
		std::map<size_t, size_t> boundaryEid_dualVid;
		for (auto& e : mesh.E)
			if (e.isBoundary) boundaryEid_dualVid[e.id] = mesh.F.size() + numOfBoundaryE++;
		for (auto& v : mesh.V) {
			auto& c = C.at(v.id);
			c.id = v.id;
			c.cellType = VTK_POLYGON;
			if (!v.isBoundary) c.Vids = get_link_fids(mesh, v);
			else c.Vids = get_link_fids(mesh, v, boundaryEid_dualVid);
		}
	}
}

void DualMesh::BuildC(const Mesh& mesh, std::map<size_t, size_t>& boundaryFid_dualVid, std::map<size_t, size_t>& vid_dualFid) {
	Cell c;
	for (auto& v : mesh.V) {
		c.id = v.id;
		c.cellType = VTK_POLYHEDRON;
		std::set<size_t> vids;
		c.Fids.clear();
		for (auto eid : v.N_Eids) {
			auto& e = mesh.E.at(eid);
			std::vector<size_t> linkVids;
			if (!e.isBoundary) linkVids = get_link_cids(mesh, e);
			else linkVids = get_link_cids(mesh, e, boundaryFid_dualVid);
			vids.insert(linkVids.begin(), linkVids.end());
			c.Fids.push_back(e.id);
		}
		c.Vids.clear();
		std::copy(vids.begin(), vids.end(), std::back_inserter(c.Vids));
		if (v.isBoundary) c.Fids.push_back(vid_dualFid[v.id]);
		C.push_back(c);
	}
}

void DualMesh::Build_skip_singularities(const Mesh& mesh) {
	m_cellType = POLYGON;
	V.resize(mesh.F.size());
	for (auto i = 0; i < mesh.F.size(); ++i) {
		auto& v = V.at(i);
		auto& c = mesh.F.at(i);
		v.id = i;
		for (auto vid : c.Vids) v += mesh.V.at(vid).xyz();
		v *= 0.25;
	}

	C.resize(mesh.V.size());
	for (auto i = 0; i < mesh.V.size(); ++i) {
		auto& v = mesh.V.at(i);
		auto& c = C.at(i);
		if (v.isSingularity) continue;
		c.id = i;
		c.Vids = v.N_Fids;
		auto fids = v.N_Fids;
		std::vector<size_t> ringFids;
		ringFids.push_back(fids.back());
		fids.pop_back();
		while (ringFids.size() < c.Vids.size()) {
			auto last_fid = ringFids.back();
			int j = 0;
			for (auto next_fid : fids) {
				auto& last_f = mesh.F.at(last_fid);
				std::set<size_t> vids(last_f.Vids.begin(), last_f.Vids.end());
				auto& next_f = mesh.F.at(next_fid);
				vids.insert(next_f.Vids.begin(), next_f.Vids.end());
				if (vids.size() == last_f.Vids.size() + next_f.Vids.size() - 2) {
					ringFids.push_back(next_fid);
					fids.erase(fids.begin() + j);
					break;
				}
				++j;
			}
		}
		c.Vids = ringFids;
	}
}

//void DualMesh::FixOrientation(const Mesh& mesh) {
//	if (m_cellType == POLYGON || m_cellType == TRIANGLE || m_cellType == QUAD) {
//		for (auto& f : F) {
//			auto f0 = mesh.F.at(f.Vids[0]);
//			auto f0v0 = mesh.V[f0.Vids[0]].xyz();
//			auto f0v1 = mesh.V[f0.Vids[1]].xyz();
//			auto f0v2 = mesh.V[f0.Vids[2]].xyz();
//			auto vec_o = glm::cross(f0v2 - f0v1, f0v0 - f0v1);
//
//			auto fv0 = this->V[f.Vids[0]].xyz();
//			auto fv1 = this->V[f.Vids[1]].xyz();
//			auto fv2 = this->V[f.Vids[2]].xyz();
//			auto vec_f = glm::cross(fv2 - fv1, fv0 - fv1);
//			auto sign = glm::dot(vec_f, vec_o);
//			if (sign < 0) {
//				std::cout << " std::reverse(f.Vids.begin(), f.Vids.end()); fid = " << f.id << std::endl;
//				std::reverse(f.Vids.begin(), f.Vids.end());
//			}
//		}
//	}
//}

static void Print(Mesh& mesh, const Face& f) {
	std::cout << "vids = (";
	for (auto vid : f.Vids)
		std::cout << " " << vid;
	std::cout << ")\n";
}

void DualMesh::FixOrientation(Mesh& mesh) {
	mesh.GetNormalOfSurfaceFaces();
	mesh.GetNormalOfSurfaceVertices();

	this->GetNormalOfSurfaceFaces();
	this->GetNormalOfSurfaceVertices();

	for (auto& f : F) {
		auto& f0 = mesh.F.at(f.Vids[0]);
		auto sign = glm::dot(f.normal, f0.normal);
		if (sign < 0) {
			// Print(*this, f);
			std::reverse(f.Vids.begin(), f.Vids.end());
			// Print(*this, f);
			if (m_cellType == POLYGON || m_cellType == TRIANGLE || m_cellType == QUAD) {
				std::reverse(C[f.id].Vids.begin(), C[f.id].Vids.end());
			}
		}
	}
}

void DualMesh::FixOrientation() {
	if (m_cellType == POLYGON || m_cellType == TRIANGLE || m_cellType == QUAD)
		for (auto& e : E) {
			auto& f0 = F[e.N_Fids[0]];
			auto& f1 = F[e.N_Fids[1]];
			auto feid = MAXID;
			for (auto eid : f1.Eids)
				if (eid == e.id) {
					feid = eid;
					break;
				}
			auto fvid0 = MAXID;
			auto fvid1 = MAXID;
			for (auto i = 0; i < f1.Vids.size(); ++i) {
				auto fvid0_ = f1.Vids.at(i);
				auto fvid1_ = f1.Vids.at((i + 1) % f1.Vids.size());
				if ((fvid0_ == e.Vids[0] && fvid1_ == e.Vids[1]) || (fvid0_ == e.Vids[1] && fvid1_ == e.Vids[0])) {
					fvid0 = fvid0_;
					fvid1 = fvid1_;
					break;
				}
			}
			auto& v0 = V[fvid0];
			auto& v1 = V[fvid1];
			auto f0center = 0.25 * (V[f0.Vids[0]].xyz() + V[f0.Vids[1]].xyz() + V[f0.Vids[2]].xyz() + V[f0.Vids[3]].xyz());
			auto f1center = 0.25 * (V[f1.Vids[0]].xyz() + V[f1.Vids[1]].xyz() + V[f1.Vids[2]].xyz() + V[f1.Vids[3]].xyz());
			auto vec_f = glm::cross(f1center - f0center, v1.xyz() - v0.xyz());
			auto vec_o = glm::cross(V[f0.Vids[2]].xyz() - V[f0.Vids[1]].xyz(), V[f0.Vids[0]].xyz() - V[f0.Vids[1]].xyz());
			auto sign = glm::dot(vec_f, vec_o);
			if (sign < 0) {
				std::swap(e.N_Fids[0], e.N_Fids[1]);
				//std::cout << "std::swap(e.N_Fids[0], e.N_Fids[1]);\n";
			}
		}
}

void DualMesh::BuildConnection() {
	if (m_cellType == TRIANGLE || m_cellType == QUAD || m_cellType == POLYGON)
		BuildAllConnectivities();
}

void DualMesh::BuildConnection(const Mesh& mesh) {
	if (m_cellType == HEXAHEDRA || m_cellType == TETRAHEDRA || m_cellType == POLYHEDRA) {
		BuildF_C(mesh);

		BuildE(mesh);
		BuildF_E(mesh);
		BuildC_E(mesh);

		BuildE_F(mesh);
		BuildE_C(mesh);

		BuildV_V(mesh);
		BuildV_E(mesh);
		BuildV_F(mesh);
		BuildV_C(mesh);
	}
}

void DualMesh::BuildE(const Mesh& mesh) {
	std::unordered_map<size_t, Edge> key_edges;
	for (auto& f : F) {
		Edge e(2);
		for (auto i = 0; i < f.Vids.size(); ++i) {
			e.Vids.at(0) = f.Vids.at(i);
			e.Vids.at(1) = f.Vids.at((i + 1) % f.Vids.size());
			//edges.insert(e);
			auto key0 = (e.Vids.at(0) << 32) | e.Vids.at(1);
			auto key1 = (e.Vids.at(1) << 32) | e.Vids.at(0);
			auto iter0 = key_edges.find(key0);
			auto iter1 = key_edges.find(key1);
			if (iter0 == key_edges.end() && iter1 == key_edges.end()) key_edges[key0] = e;
		}
	}
	E.resize(key_edges.size());
	size_t id = 0;
	for (auto& item : key_edges) {
		E[id].Vids = item.second.Vids;
		E[id].id = id++;
	}
}
void DualMesh::BuildV_V(const Mesh& mesh) {
	for (auto& e : E) {
		V.at(e.Vids[0]).N_Vids.push_back(e.Vids[1]);
		V.at(e.Vids[1]).N_Vids.push_back(e.Vids[0]);
	}
}
void DualMesh::BuildV_E(const Mesh& mesh) {
	for (auto& e : E) {
		V.at(e.Vids[0]).N_Eids.push_back(e.id);
		V.at(e.Vids[1]).N_Eids.push_back(e.id);
	}
}
void DualMesh::BuildV_F(const Mesh& mesh) {
	for (auto& f : F)
		for (auto vid : f.Vids)
			V.at(vid).N_Fids.push_back(f.id);
}
void DualMesh::BuildV_C(const Mesh& mesh) {
	for (auto& c : C)
		for (auto vid : c.Vids)
			V.at(vid).N_Cids.push_back(c.id);
}
void DualMesh::BuildE_V(const Mesh& mesh) {}
void DualMesh::BuildE_E(const Mesh& mesh) {}
void DualMesh::BuildE_F(const Mesh& mesh) {
	for (auto& f : F)
		for (auto eid : f.Eids)
			E.at(eid).N_Fids.push_back(f.id);
}
void DualMesh::BuildE_C(const Mesh& mesh) {
	for (auto& e : E) {
		std::set<size_t> cids;
		for (auto fid : e.N_Fids)
			for (auto cid : F.at(fid).N_Cids)
				cids.insert(cid);
		e.N_Cids.resize(cids.size());
		size_t id = 0;
		for (auto cid : cids)
			e.N_Cids[id++] = cid;
	}
}
void DualMesh::BuildF_V(const Mesh& mesh) {

}

size_t get_eid(const std::vector<Edge>& E, const Edge& edge) {
	for (auto& e : E)
		if (edge == e) return e.id;
	std::cerr << "Cannot get_eid\n";
	return MAXID;
}

void DualMesh::BuildF_E(const Mesh& mesh) {
	for (auto& f : F) {
		f.Eids.resize(f.Vids.size());
		size_t id = 0;
		Edge e(2);
		for (auto i = 0; i < f.Vids.size(); ++i) {
			e.Vids.at(0) = f.Vids.at(i);
			e.Vids.at(1) = f.Vids.at((i + 1) % f.Vids.size());
			f.Eids[id++] = get_eid(E, e);
		}
	}
}
void DualMesh::BuildF_F(const Mesh& mesh) {}
void DualMesh::BuildF_C(const Mesh& mesh) {
	for (auto& c : C) {
		for (auto fid : c.Fids) {
			auto& f = F.at(fid);
			f.N_Cids.push_back(c.id);
		}
	}
}
void DualMesh::BuildC_V(const Mesh& mesh) {}
void DualMesh::BuildC_E(const Mesh& mesh) {
	for (auto& c : C) {
		std::set<size_t> eids;
		for (auto fid : c.Fids) {
			auto& f = F.at(fid);
			eids.insert(f.Eids.begin(), f.Eids.end());
		}
		std::copy(eids.begin(), eids.end(), std::back_inserter(c.Eids));
	}
}
void DualMesh::BuildC_F(const Mesh& mesh) {}
void DualMesh::BuildC_C(const Mesh& mesh) {}