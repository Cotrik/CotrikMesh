#include "SeparatrixCollapse.h"

SeparatrixCollapse::SeparatrixCollapse() : SimplificationOperation() {}
SeparatrixCollapse::SeparatrixCollapse(Mesh& mesh_, MeshUtil& mu_,  Smoother& smoother_, std::vector<size_t> linkV, std::vector<size_t> linkE, bool half_) : SimplificationOperation(mesh_, mu_, smoother_) {
    half = half_;
    BuildSeparatrix(linkV, linkE);
}

SeparatrixCollapse::~SeparatrixCollapse() {}


void SeparatrixCollapse::BuildSeparatrix(std::vector<size_t> linkV, std::vector<size_t> linkE) {
    for (size_t i = 0; i < linkV.size(); i++) {
        auto vid = linkV.at(i);
        auto eid = i == linkV.size()-1 ? linkE.back() : linkE.at(i);
        // std::cout << eid << ": " << mesh.E.at(eid).N_Fids.size() << std::endl;
        auto collapseVids = GetCollapseVids(vid, eid);
        // std::cout << "collapse vids: " << collapseVids.size() << std::endl;
        auto& v0 = mesh.V.at(collapseVids[0]);
        auto& v1 = mesh.V.at(collapseVids[1]);
        auto& v = mesh.V.at(vid);
        if (v0.type == FEATURE && v.type != FEATURE) {
            std::swap(collapseVids[0], vid);
        } else if (v1.type == FEATURE && v.type != FEATURE) {
            std::swap(collapseVids[1], vid);
            std::swap(collapseVids[0], collapseVids[1]);
        }
        // if (v0.type == CORNER && v0.isBoundary) std::swap(collapseVids[0], vid);
        // else if (v1.type == CORNER && v1.isBoundary) std::swap(collapseVids[1], vid);
        target.push_back(vid);
        collapse.push_back(collapseVids);
    }
    SetUpdateElements();
}

bool SeparatrixCollapse::IsOperationValid() {
    bool isValid = true;
    int i = 0;
    // int count = 0;
    for (auto c: collapse) {
        if (target.at(i) == c.at(0) || target.at(i) == c.at(1) || c.at(0) == c.at(1)) isValid = false;
        if (mesh.V.at(target.at(i)).N_Fids.empty()) isValid = false;
        if (mesh.V.at(c.at(0)).N_Fids.empty()) isValid = false;
        if (mesh.V.at(c.at(1)).N_Fids.empty()) isValid = false;
        if (!IsCollapsable(target.at(i), c.at(0)) || !IsCollapsable(target.at(i), c.at(1))) isValid = false;
        i += 1;
        if (!isValid) break;
    }
    // if (count > 2) isValid = false;
    if (ranking != GetRanking()) isValid = false;
    if (!isValid) toUpdate.clear();
    return isValid;
}

void SeparatrixCollapse::SetRanking(glm::dvec3 d) {
    CheckValidity();

    ranking = GetRanking();
}

double SeparatrixCollapse::GetRanking() {
    double min = 0.0;
    double max = 0.0;
    int i = 0;
    for (auto c: collapse) {
        // min += mesh.V.at(c.at(0)).N_Vids.size() + mesh.V.at(c.at(1)).N_Vids.size();
        size_t tId = target.at(i);
        i += 1;

        std::vector<size_t> commonFs = GetIntersection(mesh.V.at(c.at(0)).N_Fids, mesh.V.at(c.at(1)).N_Fids);
        if (commonFs.empty()) continue;
        Face& f = mesh.F.at(commonFs.at(0));
        int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), tId));
        auto& v = mesh.V.at(f.Vids.at((idx+2)%f.Vids.size()));
        min += v.N_Vids.size() > 4 ? v.N_Vids.size() : v.N_Vids.size() * 4;

    }
    return min;
}

void SeparatrixCollapse::PerformOperation() {
    CheckValidity();
    if (!IsOperationValid()) return;
    std::cout << "Performing Separatrix Collapse: " << ranking << std::endl;
    std::vector<size_t> linkFaces;
    for (int i = 0; i < target.size(); i++) {
        auto& v = mesh.V.at(target.at(i));
        auto& v1 = mesh.V.at(collapse.at(i)[0]);
        auto& v2 = mesh.V.at(collapse.at(i)[1]);
        // AddContents(linkFaces, v.N_Fids);
        AddContents(linkFaces, GetIntersection(v.N_Fids, v1.N_Fids));
        AddContents(linkFaces, GetIntersection(v.N_Fids, v2.N_Fids));
        AddContents(linkFaces, GetIntersection(v1.N_Fids, v2.N_Fids));
    }
    for (int i = 0; i < target.size(); i++) {
        auto& v = mesh.V.at(target.at(i));
        auto& v1 = mesh.V.at(collapse.at(i)[0]);
        auto& v2 = mesh.V.at(collapse.at(i)[1]);
        if (v.type != FEATURE) {
            v = 0.5 * (v1.xyz() + v2.xyz());
        }
        Collapse(target.at(i), collapse.at(i)[0]);
        Collapse(target.at(i), collapse.at(i)[1]);
        for (int j = i + 1; j < target.size(); j++) {
            if (collapse.at(j)[0] == collapse.at(i)[0] || collapse.at(j)[0] == collapse.at(i)[1]) {
                collapse.at(j)[0] = target.at(i);
                break;
            }
            if (collapse.at(j)[1] == collapse.at(i)[0] || collapse.at(j)[1] == collapse.at(i)[1]) {
                collapse.at(j)[1] = target.at(i);
                break;
            }
        }
    }
    for (auto fid: linkFaces) {
        auto& f = mesh.F.at(fid);
        f.N_Fids.clear();
        for (auto vid: f.Vids) {
            UpdateContents(mesh.V.at(vid).N_Fids, linkFaces);
            SetSingularity(vid);
        }
        f.Vids.clear();
        f.Eids.clear();
    }
    for (auto tid: target) {
        auto& v = mesh.V.at(tid);
        SetSingularity(v.id);
    }
    // std::cout << "Finished Separatrix Collapse Operation" << std::endl;
}

void SeparatrixCollapse::Collapse(size_t targetId, size_t collapseId) {

    auto& targetV = mesh.V.at(targetId);
    auto& collapseV = mesh.V.at(collapseId);
    std::vector<size_t> vDiff = GetDifference(collapseV.N_Vids, std::vector<size_t>{targetV.id});
    for (auto vid: vDiff) {
        auto& v = mesh.V.at(vid);
        UpdateContents(v.N_Vids, std::vector<size_t>{collapseV.id});
        AddContents(v.N_Vids, std::vector<size_t>{targetV.id});
    }
    UpdateContents(targetV.N_Vids, std::vector<size_t>{collapseV.id});
    AddContents(targetV.N_Vids, vDiff);

    std::vector<size_t> eDiff = GetDifference(collapseV.N_Eids, targetV.N_Eids);
    std::vector<size_t> eComm = GetIntersection(collapseV.N_Eids, targetV.N_Eids);
    std::vector<size_t> edgesToRemove;
    for (auto eid: eDiff) {
        auto& e = mesh.E.at(eid);
        int idx = std::distance(e.Vids.begin(), std::find(e.Vids.begin(), e.Vids.end(), collapseId));
        e.Vids.at(idx) = targetV.id;
        for (auto t_eid: targetV.N_Eids) {
            auto& t_e = mesh.E.at(t_eid);
            if (std::find(e.Vids.begin(), e.Vids.end(), t_e.Vids.at(0)) != e.Vids.end() &&
                std::find(e.Vids.begin(), e.Vids.end(), t_e.Vids.at(1)) != e.Vids.end()) {
                    auto& edgeV = mesh.V.at(e.Vids.at((idx+1)%e.Vids.size()));
                    if (GetIntersection(e.N_Fids, t_e.N_Fids).empty()) continue;
                    UpdateContents(edgeV.N_Eids, std::vector<size_t>{e.id});
                    auto& faceToRemove = mesh.F.at(GetIntersection(e.N_Fids, t_e.N_Fids).at(0));
                    auto& faceToKeep = mesh.F.at(GetDifference(e.N_Fids, t_e.N_Fids).at(0));
                    UpdateContents(t_e.N_Fids, std::vector<size_t>{faceToRemove.id});
                    AddContents(t_e.N_Fids, std::vector<size_t>{faceToKeep.id});
                    UpdateContents(faceToKeep.Eids, std::vector<size_t>{e.id});
                    AddContents(faceToKeep.Eids, std::vector<size_t>{t_e.id});
                    edgesToRemove.push_back(e.id);
                    break;
            }
        }
    }
    UpdateContents(eDiff, edgesToRemove);
    AddContents(eComm, edgesToRemove);
    // UpdateContents(targetV.N_Eids, edgesToRemove);
    UpdateContents(targetV.N_Eids, GetIntersection(targetV.N_Eids, collapseV.N_Eids));
    AddContents(targetV.N_Eids, eDiff);

    std::vector<size_t> fDiff = GetDifference(collapseV.N_Fids, targetV.N_Fids);
    std::vector<size_t> fComm = GetIntersection(targetV.N_Fids, collapseV.N_Fids);
    for (auto fid: fDiff) {
        auto& f = mesh.F.at(fid);
        int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), collapseId));
        f.Vids.at(idx) = targetV.id;
        UpdateContents(f.N_Fids, targetV.N_Fids);
    }
    UpdateContents(targetV.N_Fids, fComm);
    AddContents(targetV.N_Fids, fDiff);

    for (auto eid: eComm) {
        mesh.E.at(eid).Vids.clear();
        mesh.E.at(eid).N_Fids.clear();
    }
    collapseV.N_Vids.clear();
    collapseV.N_Eids.clear();
    collapseV.N_Fids.clear();
}

void SeparatrixCollapse::SetUpdateElements() {
    std::vector<size_t> outerFaces;
    for (int i = 0; i < collapse.size(); i++) {
        AddContents(outerFaces, GetDifference(mesh.V.at(collapse.at(i).at(0)).N_Fids, mesh.V.at(target.at(i)).N_Fids));
        AddContents(outerFaces, GetDifference(mesh.V.at(collapse.at(i).at(1)).N_Fids, mesh.V.at(target.at(i)).N_Fids));
        size_t tId = target.at(i);
        for (auto fid: GetIntersection(mesh.V.at(collapse.at(i).at(0)).N_Fids, mesh.V.at(collapse.at(i).at(1)).N_Fids)) {
            auto& f = mesh.F.at(fid);
            int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), tId));
            auto& v = mesh.V.at(f.Vids.at((idx+2)%f.Vids.size()));
            AddContents(outerFaces, GetDifference(v.N_Fids, mesh.V.at(target.at(i)).N_Fids));
        }
    }
    for (auto fid: outerFaces) {
        auto& f = mesh.F.at(fid);
        for (auto vid: f.Vids) {
            if (mesh.V.at(vid).N_Vids.size() == 3) AddContents(toUpdate, std::vector<size_t>{vid});
        }
    }
}

std::vector<size_t> SeparatrixCollapse::GetCollapseVids(size_t vid, size_t eid) {
    std::vector<size_t> res;
	auto& v = mesh.V.at(vid);
	auto& e = mesh.E.at(eid);
	for (auto fid : e.N_Fids) {
		auto& f = mesh.F.at(fid);
        // std::cout << "f eids: " << f.Eids.size() << std::endl;
		for (auto feid : f.Eids) {
			auto& fe = mesh.E.at(feid);
			std::set<int> vids(fe.Vids.begin(), fe.Vids.end());
			vids.insert(e.Vids.begin(), e.Vids.end());
            // std::cout << "e vids: " << e.id << " " << e.Vids.at(0) << " " << e.Vids.at(1) << std::endl;
            // std::cout << "num fe vids: " << fe.Vids.size() << std::endl;
            // std::cout << "fe vids: " << feid << " " << fe.Vids.at(0) << " " << fe.Vids.at(1) << std::endl;
            // std::cout << "vids: " << vids.size() << std::endl; 
			if (vids.size() == 4) {
				for (auto nvid : v.N_Vids) {
					if (nvid == fe.Vids[0] || nvid == fe.Vids[1]) {
						res.push_back(nvid);
						break;
					}
				}
				break;
			}
		}
	}
	return res;
}