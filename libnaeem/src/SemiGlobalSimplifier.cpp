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

void SemiGlobalSimplifier::SetIters(int iters_) {
    iters = iters_;
}

void SemiGlobalSimplifier::SetSimplificationOperations() {
    CheckValidity();
    
    // SetDiagonalCollapseOperations();
    SetDirectSeparatrixOperations();
}

void SemiGlobalSimplifier::FixBoundary() {
    CheckValidity();

    // int i = 0;
    for (int i = 0; i < mesh.V.size(); i++) {
        // if (i >= iters) break;
        auto& v = mesh.V.at(i);
        int valence = v.N_Fids.size();
        bool failed = false;
        if (v.type == FEATURE && valence > 4) {
            int featureCount = 0;
            for (auto vid: v.N_Vids) {
                if (mesh.V.at(vid).type == FEATURE) featureCount += 1;
            }
            if (featureCount > 2) continue;
            bool performVertexSplit = true;
            for (int j = 0; j < v.N_Eids.size(); j++) {
                auto& e = mesh.E.at(v.N_Eids.at(j));
                if (mesh.V.at(e.Vids.at(0)).type == FEATURE && mesh.V.at(e.Vids.at(1)).type == FEATURE) continue;
                int count = 0;
                for (auto fid: e.N_Fids) {
                    auto& f = mesh.F.at(fid);
                    for (auto vid: f.Vids) {
                        if (mesh.V.at(vid).type == FEATURE) count += 1;
                    }
                }
                if (count == 4) performVertexSplit = false;
            }
            for (int j = 0; j < v.N_Fids.size(); j++) {
                int count = 0;
                auto& f = mesh.F.at(v.N_Fids.at(j));
                for (auto vid: f.Vids) {
                    if (mesh.V.at(vid).type == FEATURE) count += 1;
                }
                if (count > 2) performVertexSplit = false;
            }
            if (performVertexSplit && v.N_Vids.size() > 5) {
                std::shared_ptr<SimplificationOperation> s = std::make_shared<VertexSplit>(mesh, mu, v.id);
                s->PerformOperation();
                continue;
            }
            while (valence > 4) {
                bool breakLoop = false;
                for (int j = 0; j < v.N_Eids.size(); j++) {
                    auto& e = mesh.E.at(v.N_Eids.at(j));
                    if (mesh.V.at(e.Vids.at(0)).type == FEATURE && mesh.V.at(e.Vids.at(1)).type == FEATURE) continue;
                    int count = 0;
                    for (auto fid: e.N_Fids) {
                        auto& f = mesh.F.at(fid);
                        for (auto vid: f.Vids) {
                            if (mesh.V.at(vid).type == FEATURE) count += 1;
                        }
                    }
                    if (count == 4) continue;
                    bool clockwise = false;
                    for (auto fid: e.N_Fids) {
                        auto& f = mesh.F.at(fid);
                        int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
                        if (mesh.V.at(f.Vids.at((idx+1)%f.Vids.size())).type == FEATURE) {
                            clockwise = true;
                        }
                    }
                    std::shared_ptr<SimplificationOperation> s = std::make_shared<EdgeRotation>(mesh, mu, e.id, clockwise);
                    s->PerformOperation();
                }
                if (v.N_Fids.size() == valence) {
                    failed = true;
                    break;
                }
                valence = v.N_Fids.size();
            }
            // i += 1;
            if (failed) break;
        }
        if (v.isBoundary && v.N_Fids.size() != 2) {

        }
    }
}

void SemiGlobalSimplifier::SetDiagonalCollapseOperations() {
    CheckValidity();

    for (auto& f: mesh.F) {
        if (f.N_Fids.size() == 0 || f.Vids.empty()) continue;
        for (int i = 0; i < f.Vids.size(); i++) {
            if (mesh.V.at(f.Vids.at(i)).N_Fids.size() == 3 && mesh.V.at(f.Vids.at((i+2)%f.Vids.size())).N_Fids.size() == 3) {
                    std::shared_ptr<SimplificationOperation> dc1 = std::make_shared<DiagonalCollapse>(mesh, mu, f.id, i, (i+2)%f.Vids.size());
                    dc1->SetRanking();
                    Ops.push_back(dc1);
                break;
            }
        }
        
        // std::unique_ptr<SimplificationOperation> dc2 = std::make_unique<DiagonalCollapse>(mesh, mu, f.id, 1, 3);
        // dc2->SetRanking();
        // Ops.push_back(std::move(dc2));
    }
    std::cout << Ops.size() << std::endl;
    int i = 0;
    for (auto& op: Ops) {
        op->PerformOperation();
        i += 1;
        if (i >= iters) break;        
    }
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
    std::cout << Op_Q.size() << " direct separatrix operations" << std::endl;
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
    Op_Q.clear();
    Op_Q.setMinQueueOn();
    // Op_Q.setMaxQueueOn();
    for (auto& v: mesh.V) {
        if (!v.isSingularity || v.N_Fids.empty()) continue;
        std::vector<SingularityLink> links = TraceSingularityLinks(v, bc);
        for (int i = 0; i < links.size(); i++) {
            auto& linkV = links.at(i).linkVids;
            auto& linkE = links.at(i).linkEids;
            auto& v_front = mesh.V.at(linkV.front());
            auto& v_back = mesh.V.at(linkV.back());
            if (v_front.isBoundary || v_back.isBoundary) continue;
            if (v_front.N_Fids.size() != 3 || v_back.N_Fids.size() != 3) continue;
            if (mesh.V.at(GetDiagonalV(v_front.id, GetFaceId(v_front.id, linkV.at(1)))).N_Fids.size() <= 4 || mesh.V.at(GetDiagonalV(v_back.id, GetFaceId(v_back.id, linkV.at(linkV.size()-2)))).N_Fids.size() <= 4) continue;

            std::shared_ptr<SimplificationOperation> s = std::make_shared<SeparatrixCollapse>(mesh, mu, linkV, linkE);
            s->SetRanking();
            Op_Q.insert(s->ranking, s->GetCenterId(), s);
        }
    }
    std::cout << Op_Q.size() << " separatrix operations" << std::endl;
    int i = 0;
    while (!Op_Q.empty()) {
        auto& op = Op_Q.pop();
        op->PerformOperation();
        // std::cout << "Getting to update operations" << std::endl;
        for (auto vid: op->toUpdate) {
            if (!mesh.V.at(vid).isSingularity || mesh.V.at(vid).N_Fids.empty()) continue;
            // std::cout << mesh.V.at(vid).N_Vids.size() << " " << mesh.V.at(vid).N_Eids.size() << " " << mesh.V.at(vid).N_Fids.size() << std::endl;
            std::vector<SingularityLink> links = TraceSingularityLinks(mesh.V.at(vid), bc);
            // std::cout << "links " << links.size() << std::endl;
            for (int i = 0; i < links.size(); i++) {
                auto& linkV = links.at(i).linkVids;
                auto& linkE = links.at(i).linkEids;
                auto& v_front = mesh.V.at(linkV.front());
                auto& v_back = mesh.V.at(linkV.back());
                if (v_front.isBoundary || v_back.isBoundary) continue;
                if (v_front.N_Fids.size() != 3 || v_back.N_Fids.size() != 3) continue;
                if (mesh.V.at(GetDiagonalV(v_front.id, GetFaceId(v_front.id, linkV.at(1)))).N_Fids.size() <= 4 || mesh.V.at(GetDiagonalV(v_back.id, GetFaceId(v_back.id, linkV.at(linkV.size()-2)))).N_Fids.size() <= 4) continue;

                // std::cout << "making a separatrix operation" << std::endl;
                std::shared_ptr<SimplificationOperation> s = std::make_shared<SeparatrixCollapse>(mesh, mu, linkV, linkE);
                // std::cout << "made a separatrix operation" << std::endl;
                s->SetRanking();
                Op_Q.insert(s->ranking, s->GetCenterId(), s);
            }
            // std::cout << "got updated operation" << std::endl;
        }
        i += 1;
        // if (i > iters) break;
    }
    // std::cout << "Finished all separatrix collapse operations" << std::endl;
    // std::cout << "Smoothing mesh" << std::endl;
    // std::vector<size_t> smoothv;
    // for (auto& v: mesh.V) {
    //     if (v.N_Vids.empty() || v.N_Eids.empty() || v.N_Fids.empty()) continue;
    //     smoothv.push_back(v.id);
    // }
    // smoother.Smooth(smoothv);
}

void SemiGlobalSimplifier::SetHalfSeparatrixOperations() {
    CheckValidity();

    BaseComplexQuad bc(mesh);
    Op_Q.clear();
    Op_Q.setMaxQueueOn();
    for (auto& v: mesh.V) {
        if (!v.isSingularity) continue;
        std::vector<SingularityLink> links = TraceSingularityLinks(v, bc);
        for (int i = 0; i < links.size(); i++) {
            auto& linkV = links.at(i).linkVids;
            auto& linkE = links.at(i).linkEids;
            auto& v_front = mesh.V.at(linkV.front());
            auto& v_back = mesh.V.at(linkV.back());
                
            if (!(v_front.isBoundary ^ v_back.isBoundary)) continue;
            if (v_front.isBoundary && !v_back.isBoundary) {
                std::reverse(linkV.begin(), linkV.end());
                std::reverse(linkE.begin(), linkE.end());
                auto& v_front = mesh.V.at(linkV.front());
                auto& v_back = mesh.V.at(linkV.back());
            }
            if (!(v_front.N_Fids.size() == 3 && v_back.N_Fids.size() == 2)) continue;
            std::shared_ptr<SimplificationOperation> s = std::make_shared<SeparatrixCollapse>(mesh, mu, linkV, linkE, true);
            s->SetRanking();
            Op_Q.insert(s->ranking, s->GetCenterId(), s);
        }
    }
    std::cout << Op_Q.size() << " half separatrix operations" << std::endl;
    int i = 0;
    while (!Op_Q.empty()) {
        auto& op = Op_Q.pop();
        op->PerformOperation();
        for (auto vid: op->toUpdate) {
            std::vector<SingularityLink> links = TraceSingularityLinks(mesh.V.at(vid), bc);
            for (int i = 0; i < links.size(); i++) {
                auto& linkV = links.at(i).linkVids;
                auto& linkE = links.at(i).linkEids;
                auto& v_front = mesh.V.at(linkV.front());
                auto& v_back = mesh.V.at(linkV.back());
                
                if (!(v_front.isBoundary ^ v_back.isBoundary)) continue;
                if (v_front.isBoundary && !v_back.isBoundary) {
                    std::reverse(linkV.begin(), linkV.end());
                    std::reverse(linkE.begin(), linkE.end());
                    auto& v_front = mesh.V.at(linkV.front());
                    auto& v_back = mesh.V.at(linkV.back());
                }
                if (!(v_front.N_Fids.size() == 3 && v_back.N_Fids.size() == 2)) continue;
                std::shared_ptr<SimplificationOperation> s = std::make_shared<SeparatrixCollapse>(mesh, mu, linkV, linkE, true);
                s->SetRanking();
                Op_Q.insert(s->ranking, s->GetCenterId(), s);
            }
        }
        i += 1;
        // if (i >= iters) break;
    }
    // std::cout << "Smoothing mesh" << std::endl;
    // std::vector<size_t> smoothv;
    // for (auto& v: mesh.V) {
    //     if (v.N_Vids.empty() || v.N_Eids.empty() || v.N_Fids.empty()) continue;
    //     smoothv.push_back(v.id);
    // }
    // smoother.Smooth(smoothv);
}

void SemiGlobalSimplifier::SetChordCollapseOperations() {
    CheckValidity();

    ChordExtractor ce(mesh);
    ce.Extract();
    // std::vector<size_t> chords = ce.SelectChords();
    // std::cout << chords.size() << std::endl;
    // ce.Write(std::vector<size_t>{chords.at(0)});
    // ce.Write(chords);
    // return;

    for (auto& chord: ce.Chords) {
        std::shared_ptr<SimplificationOperation> s = std::make_shared<ChordCollapse>(mesh, mu, ce, chord.id);
        s->PerformOperation();
    }
    
    // std::cout << "Smoothing mesh" << std::endl;
    // std::vector<size_t> smoothv;
    // for (auto& v: mesh.V) {
    //     if (v.N_Vids.empty() || v.N_Eids.empty() || v.N_Fids.empty()) continue;
    //     smoothv.push_back(v.id);
    // }
    // smoother.Smooth(smoothv);
}

void SemiGlobalSimplifier::SetEdgeRotationOperations() {
    CheckValidity();

    int i = 0;
    for (auto& e: mesh.E) {
        std::shared_ptr<SimplificationOperation> s = std::make_shared<EdgeRotation>(mesh, mu, e.id, true);
        s->PerformOperation();
        i += 1;
        // if (i > iters) break;
    }
    // std::cout << "Smoothing mesh" << std::endl;
    // std::vector<size_t> smoothv;
    // for (auto& v: mesh.V) {
    //     if (v.N_Vids.empty() || v.N_Eids.empty() || v.N_Fids.empty()) continue;
    //     smoothv.push_back(v.id);
    // }
    // smoother.Smooth(smoothv);
}

void SemiGlobalSimplifier::SetVertexRotationOperations() {
    CheckValidity();

    int i = 0;
    for (auto& v: mesh.V) {
        std::shared_ptr<SimplificationOperation> s = std::make_shared<VertexRotation>(mesh, mu, v.id);
        s->PerformOperation();
        // i += 1;
        // if (i > iters) break;
    }
    // std::cout << "Smoothing mesh" << std::endl;
    // std::vector<size_t> smoothv;
    // for (auto& v: mesh.V) {
    //     if (v.N_Vids.empty() || v.N_Eids.empty() || v.N_Fids.empty()) continue;
    //     smoothv.push_back(v.id);
    // }
    // smoother.Smooth(smoothv);
}

void SemiGlobalSimplifier::SetVertexSplitOperations() {
    CheckValidity();

    int i = 0;
    for (int j = 0; j < mesh.V.size(); j++) {
        auto& v = mesh.V.at(j);
        if (v.type != FEATURE) continue;
        if (v.N_Vids.size() < 6) continue;
        std::shared_ptr<SimplificationOperation> s = std::make_shared<VertexSplit>(mesh, mu, v.id);
        s->PerformOperation();
        // i += 1;
        // if (i > iters) break;
        // break;
    }
}

void SemiGlobalSimplifier::SetEdgeCollapseOperations() {
    CheckValidity();

    Edge& e = mesh.E.at(0);
    size_t source_id = e.Vids.at(0);
    size_t targte_id = e.Vids.at(1);
    std::shared_ptr<SimplificationOperation> vr = std::make_shared<VertexRotation>(mesh, mu, targte_id);
    std::shared_ptr<SimplificationOperation> dc = std::make_shared<DiagonalCollapse>(mesh, mu, -1, targte_id, source_id);
    std::shared_ptr<SimplificationOperation> s = std::make_shared<EdgeCollapse>(mesh, mu, vr, dc);
    s->PerformOperation();
}


std::vector<SingularityLink> SemiGlobalSimplifier::TraceSingularityLinks(Vertex& v, BaseComplexQuad& bc) {
    std::vector<bool> is_mesh_edge_visited(mesh.E.size(), false);
    std::vector<SingularityLink> links;
    for (auto edgeid : v.N_Eids) {
        const Edge& edge = mesh.E.at(edgeid);
        if (!is_mesh_edge_visited[edgeid]) {
            SingularityLink l;
            std::vector<size_t> link_vids;
            std::vector<size_t> link_eids;
            bc.TraceAlongEdge(v, edge, is_mesh_edge_visited, l.linkVids, l.linkEids);
            links.push_back(l);
        }
    }
    return links;
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

void SemiGlobalSimplifier::Smooth() {
    std::cout << "Smoothing mesh" << std::endl;
    // std::vector<size_t> smoothv;
    // for (auto& v: mesh.V) {
    //     if (v.N_Vids.empty() || v.N_Eids.empty() || v.N_Fids.empty()) continue;
    //     smoothv.push_back(v.id);
    // }
    // smoother.Smooth(smoothv);
    smoother.Smooth(std::vector<size_t>{});
}

// std::cout << "Writing output file" << std::endl;
// std::string outputf = "links.vtk";
// std::ofstream ofs(outputf.c_str());
// ofs << "# vtk DataFile Version 3.0\n"
//     << outputf.c_str() << "\n"
//     << "ASCII\n\n"
//     << "DATASET UNSTRUCTURED_GRID\n";
// ofs << "POINTS " << mesh.V.size() << " double\n";
// std::vector<size_t> c_indices;
// for (auto vid: target) {
//     c_indices.push_back(vid);
// }
// for (auto c: collapse) {
//     c_indices.push_back(c.at(0));
//     c_indices.push_back(c.at(1));
// }
// // std::vector<size_t> c_indices = {12, 296};
// // std::cout << c_indices.size() << std::endl;
// for (size_t i = 0; i < mesh.V.size(); i++) {
//     ofs << std::fixed << std::setprecision(7) <<  mesh.V.at(i).x << " " <<  mesh.V.at(i).y << " " <<  mesh.V.at(i).z << "\n";
// }
// ofs << "CELLS " << c_indices.size() << " " << 2 * c_indices.size() << std::endl;
// for (size_t i = 0; i < c_indices.size(); i++) {
//     ofs << "1 " << c_indices.at(i) << std::endl;
// }
// ofs << "CELL_TYPES " << c_indices.size() << "\n";
// for (size_t i = 0; i < c_indices.size(); i++) {
//     ofs << "1" << std::endl;
// }