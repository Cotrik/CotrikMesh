#include <algorithm>
#include "SemiGlobalSimplifier.h"

SemiGlobalSimplifier::SemiGlobalSimplifier() {}

SemiGlobalSimplifier::SemiGlobalSimplifier(Mesh& mesh_, MeshUtil& mu_, Smoother& smoother_) : mesh(mesh_), mu(mu_), smoother(smoother_) {
    // mu.SetMesh(mesh);
    // smoother.SetMesh(mesh);
}

SemiGlobalSimplifier::~SemiGlobalSimplifier() {}

void SemiGlobalSimplifier::CheckValidity() {
    if (mesh.V.size() == 0 || mesh.F.size() == 0 || mesh.C.size() == 0) {
        std::cout << "No mesh to use for Semi Global Simplifier." << std::endl;
        exit(0);
    }
}

void SemiGlobalSimplifier::SetMesh(Mesh& mesh_) {
    mesh = mesh_;
    mu.SetMesh(mesh);
    smoother.SetMesh(mesh);
}

void SemiGlobalSimplifier::SetSimplificationOperations() {
    CheckValidity();
    
    // SetDiagonalCollapseOperations();
    SetDirectSeparatrixOperations();
}

void SemiGlobalSimplifier::SetDiagonalCollapseOperations() {
    CheckValidity();

    for (auto& f: mesh.F) {
        std::unique_ptr<SimplificationOperation> dc1 = std::make_unique<DiagonalCollapse>(mesh, mu, f.id, 0, 2);
        dc1->SetRanking();
        Ops.push_back(std::move(dc1));

        std::unique_ptr<SimplificationOperation> dc2 = std::make_unique<DiagonalCollapse>(mesh, mu, f.id, 1, 3);
        dc2->SetRanking();
        Ops.push_back(std::move(dc2));
    }

    int i = 0;
    for (auto& op: Ops) {
        // op->PerformOperation();
        std::cout << op->ranking << std::endl;
        // if (i == 10) {
        //     break;
        // }
        // i += 1;
    }
    std::cout << Ops.size() << std::endl;

}

void SemiGlobalSimplifier::SetDirectSeparatrixOperations() {
    CheckValidity();

    Op_Q.setMaxQueueOn();
    for (auto& v: mesh.V) {
        if (v.N_Vids.size() != 4 || v.type == FEATURE) continue;
        std::vector<size_t> c1, c2;
        for (auto vid: v.N_Vids) mesh.V.at(vid).N_Vids.size() == 3 ? c1.push_back(vid) : c2.push_back(vid);
        if (c1.size() < 2) continue;
        auto& s3_v1 = mesh.V.at(c1.at(0));
        auto& s3_v2 = mesh.V.at(c1.at(1));
        if (mu.GetDifference(s3_v1.N_Fids, s3_v2.N_Fids).size() != s3_v1.N_Fids.size()) continue;

        std::shared_ptr<SimplificationOperation> ds = std::make_shared<DirectSeparatrixCollapse>(mesh, mu, v.id, c1, c2);
        ds->SetRanking();
        if (ds->ranking < 0) continue;
        Op_Q.insert(ds->ranking, v.id, ds);
    }
    int i = 0;
    while (!Op_Q.empty()) {
        auto& op = Op_Q.pop();
        op->PerformOperation();
        // smoother.Smooth(op->smoothV);
        for (auto key: op->toUpdate) {
            auto nop = Op_Q.getByKey(key);
            if (!nop) continue;
            nop->SetRanking(op->GetLocation());
            Op_Q.update(nop->ranking, key);
        }
    }

    // std::cout << "Smoothing mesh" << std::endl;
    // std::vector<size_t> smoothv;
    // for (auto& v: mesh.V) smoothv.push_back(v.id);
    // smoother.Smooth(smoothv);

}

void SemiGlobalSimplifier::SetSeparatrixOperations() {
    CheckValidity();

    BaseComplexQuad bc(mesh);
    for (auto& v: mesh.V) {
        if (!v.isSingularity) continue;
        TraceSingularityLinks(v, bc);
    }
    for (int i = 0; i < bc.separatedVertexIdsLink.size(); i++) {
        auto& linkV = bc.separatedVertexIdsLink.at(i);
        auto& linkE = bc.separatedEdgeIdsLink.at(i);

        std::cout << "s1: " << linkV.front() << " valence: " << mesh.V.at(linkV.front()).N_Vids.size() << " s2: " << linkV.back() << " valence: " << mesh.V.at(linkV.back()).N_Vids.size() << std::endl; 
        std::cout << "link V: " << linkV.size() << " link e: " << linkE.size() << std::endl;
        std::cout << "*******************************" << std::endl;
    }
    std::cout << "printed all singularities links " << bc.separatedVertexIdsLink.size() << std::endl;
}

void SemiGlobalSimplifier::TraceSingularityLinks(Vertex& v, BaseComplexQuad& bc) {
    std::vector<bool> is_mesh_edge_visited(mesh.E.size(), false);
    for (auto edgeid : v.N_Eids) {
        const Edge& edge = mesh.E.at(edgeid);
        if (!is_mesh_edge_visited[edgeid]) {
            std::vector<size_t> link_vids;
            std::vector<size_t> link_eids;
            bc.TraceAlongEdge(v, edge, is_mesh_edge_visited, link_vids, link_eids);

            auto& v_front = mesh.V.at(link_vids.front());
            auto& v_back = mesh.V.at(link_vids.back());
            if (v_front.isBoundary || v_back.isBoundary) continue;
            if (v_front.N_Fids.size() != 3 || v_back.N_Fids.size() != 3) continue;
            if (mesh.V.at(GetDiagonalV(v_front.id, GetFaceId(v_front.id, link_vids.at(1)))).N_Fids.size() <= 4 || mesh.V.at(GetDiagonalV(v_back.id, GetFaceId(v_back.id, link_vids.at(link_vids.size()-2)))).N_Fids.size() <= 4) continue;

            bc.separatedVertexIdsLink.push_back(link_vids);
            bc.separatedEdgeIdsLink.push_back(link_eids);
        }
    }
}

void SemiGlobalSimplifier::PerformGlobalOperations() {
    CheckValidity();
}

size_t SemiGlobalSimplifier::GetFaceId(size_t vid, size_t exclude_vid) {
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

size_t SemiGlobalSimplifier::GetDiagonalV(size_t vid, size_t fid) {
	auto& f = mesh.F.at(fid);
    int index = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), vid));
	return f.Vids.at((index+2)%f.Vids.size());
}