#include "SeparatrixCollapse.h"

SeparatrixCollapse::SeparatrixCollapse() : SimplificationOperation() {}
SeparatrixCollapse::SeparatrixCollapse(Mesh& mesh_, MeshUtil& mu_, std::vector<size_t> linkV, std::vector<size_t> linkE) : SimplificationOperation(mesh_, mu_) {

}

SeparatrixCollapse::~SeparatrixCollapse() {}

bool SeparatrixCollapse::IsOperationValid() {
    return false;
}

void SeparatrixCollapse::BuildSeparatrix(std::vector<size_t> linkV, std::vector<size_t> linkE) {
    for (size_t i = 0; i < linkV.size(); i++) {
        auto vid = linkV.at(i);
        auto eid = i == linkV.size()-1 ? linkE.back() : linkE.at(i);
        auto collapseVids = GetCollapseVids(vid, eid);
        auto& v0 = mesh.V.at(collapseVids[0]);
        auto& v1 = mesh.V.at(collapseVids[1]);
        auto& v = mesh.V.at(vid);
        if (v0.type == CORNER && v0.isBoundary) std::swap(collapseVids[0], vid);
        else if (v1.type == CORNER && v1.isBoundary) std::swap(collapseVids[1], vid);
        target.push_back(vid);
        collapse.push_back(collapseVids);
    }
}

void SeparatrixCollapse::SetRanking(glm::dvec3 d) {
    CheckValidity();
}

void SeparatrixCollapse::PerformOperation() {
    CheckValidity();
    if (!IsOperationValid()) return;
}

void SeparatrixCollapse::Collapse(size_t targetId, std::vector<size_t> collapseIds) {
    auto& target = mesh.V.at(targetId);
    auto& c1 = mesh.V.at(collapseIds.at(0));
    auto& c2 = mesh.V.at(collapseIds.at(1));

}

std::vector<size_t> SeparatrixCollapse::GetCollapseVids(size_t vid, size_t eid) {
    std::vector<size_t> res;
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