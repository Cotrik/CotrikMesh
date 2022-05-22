#include "ChordCollapse.h"

ChordCollapse::ChordCollapse() {}
ChordCollapse::ChordCollapse(Mesh& mesh_, MeshUtil& mu_, Smoother& smoother_, ChordExtractor& ce_, size_t chordId_) : SimplificationOperation(mesh_, mu_, smoother_), ce(ce_) {
    chord = ce.Chords.at(chordId_);
}
ChordCollapse::~ChordCollapse() {}

void ChordCollapse::SetRanking(glm::dvec3 d) {
    CheckValidity();
}

bool ChordCollapse::IsOperationValid() {
    CheckValidity();

    // return ce.isChordValid(chord.id);
    bool isValid = false;
    for (auto eId: chord.Edges) {
        Edge& e = ce.Edges.at(eId);
        for (auto eid: e.Vids) {
            Edge& meshEdge = mesh.E.at(eid);
            if (!IsCollapsable(meshEdge.Vids[0], meshEdge.Vids[1])) {
                isValid = false;
                break;
            }
        }
        if (!isValid) break;
    }
    return isValid;
}

void ChordCollapse::PerformOperation() {
    CheckValidity();

    if (!IsOperationValid()) return;
    std::cout << "Performing Chord Collapse" << std::endl;
    int i = 0;
    for (auto& eid: chord.Edges) {
        Edge& e = ce.Edges.at(eid);
        if (e.Vids.empty()) continue;
        Edge& e1 = mesh.E.at(e.Vids.at(0));
        Edge& e2 = mesh.E.at(e.Vids.at(1));
        size_t fid = GetIntersection(e1.N_Fids, e2.N_Fids).at(0);
        std::vector<size_t> diffEdges = GetDifference(mesh.F.at(fid).Eids, std::vector<size_t>{e1.id, e2.id});
        size_t edgeToKeepId = diffEdges.at(0);
        size_t edgeToRemoveId = diffEdges.at(1);
        if (mesh.E.at(edgeToRemoveId).isBoundary) {
            edgeToKeepId = diffEdges.at(1);
            edgeToRemoveId = diffEdges.at(0);
        }
        
        CollapseEdge(e1, mesh.E.at(edgeToRemoveId));
        CollapseEdge(e2, mesh.E.at(edgeToRemoveId));
        UpdateEdges(fid, mesh.E.at(edgeToKeepId), mesh.E.at(edgeToRemoveId));
    }
    for (auto& eid: chord.Edges) {
        Edge& e = ce.Edges.at(eid);
        if (e.Vids.empty()) continue;
        Edge& e1 = mesh.E.at(e.Vids.at(0));
        Edge& e2 = mesh.E.at(e.Vids.at(1));
        e1.N_Fids.clear();
        e2.N_Fids.clear();
    }
}

void ChordCollapse::CollapseEdge(Edge& e, Edge& edgeToRemove) {
    if (e.Vids.empty()) return;
    int idx = std::distance(edgeToRemove.Vids.begin(), std::find(edgeToRemove.Vids.begin(), edgeToRemove.Vids.end(), e.Vids.at(0)));
    size_t source_id = e.Vids.at(0);
    size_t target_id = e.Vids.at(1);
    if (idx > 1) {
        source_id = e.Vids.at(1);
        target_id = e.Vids.at(0);
    }
    Vertex& source = mesh.V.at(source_id);
    Vertex& target = mesh.V.at(target_id);
    if (!(target.isBoundary ^ source.isBoundary)) target = 0.5 * (source.xyz() + target.xyz());
    
    UpdateContents(target.N_Vids, std::vector<size_t>{source.id});
    AddContents(target.N_Vids, GetDifference(source.N_Vids, std::vector<size_t>{target.id}));

    UpdateContents(target.N_Eids, std::vector<size_t>{e.id});
    AddContents(target.N_Eids, GetDifference(source.N_Eids, std::vector<size_t>{e.id}));

    UpdateContents(target.N_Fids, e.N_Fids);
    AddContents(target.N_Fids, GetDifference(source.N_Fids, e.N_Fids));

    for (auto nvid: source.N_Vids) {
        if (nvid == target.id) continue;
        auto& v = mesh.V.at(nvid);
        UpdateContents(v.N_Vids, std::vector<size_t>{source.id});
        AddContents(v.N_Vids, std::vector<size_t>{target.id});
    }

    std::vector<size_t> eids = GetDifference(source.N_Eids, std::vector<size_t>{e.id});
    for (auto neid: eids) {
        Edge& ne = mesh.E.at(neid);
        ne.Vids.at(0) == source.id ? ne.Vids.at(0) = target.id : ne.Vids.at(1) = target.id;
    }

    std::vector<size_t> fids = GetDifference(source.N_Fids, e.N_Fids);
    std::vector<size_t> tfids = GetDifference(target.N_Fids, e.N_Fids);
    for (auto nfid: fids) {
        Face& nf = mesh.F.at(nfid);
        int idx = std::distance(nf.Vids.begin(), std::find(nf.Vids.begin(), nf.Vids.end(), source.id));
        nf.Vids.at(idx) = target.id;
        UpdateContents(nf.N_Fids, e.N_Fids);
        AddContents(nf.N_Fids, tfids);
    }
    for (auto nfid: tfids) {
        Face& nf = mesh.F.at(nfid);
        UpdateContents(nf.N_Fids, e.N_Fids);
        AddContents(nf.N_Fids, fids);
    }

    source.N_Vids.clear();
    source.N_Eids.clear();
    source.N_Fids.clear();

    e.Vids.clear();
}

void ChordCollapse::UpdateEdges(size_t fid, Edge& edgeToKeep, Edge& edgeToRemove) {

    Face& f = mesh.F.at(fid);
    Vertex& vertexToKeep = ce.Vertices.at(edgeToKeep.id);
    Vertex& vertexToRemove = ce.Vertices.at(edgeToRemove.id);

    std::vector<size_t> commonEdges = GetIntersection(vertexToKeep.N_Eids, vertexToRemove.N_Eids);
    std::vector<size_t> diffEdges = GetDifference(vertexToRemove.N_Eids, vertexToKeep.N_Eids);
    ce.Edges.at(commonEdges.at(0)).Vids.clear();
    UpdateContents(vertexToKeep.N_Eids, commonEdges);
    AddContents(vertexToKeep.N_Eids, diffEdges);
    UpdateContents(edgeToKeep.N_Fids, std::vector<size_t>{fid});
    std::vector<size_t> facesToUpdate = GetDifference(edgeToRemove.N_Fids, std::vector<size_t>{fid});
    if (!facesToUpdate.empty()) {
        AddContents(edgeToKeep.N_Fids, facesToUpdate);
        Face& faceToUpdate = mesh.F.at(facesToUpdate.at(0));
        UpdateContents(faceToUpdate.Eids, std::vector<size_t>{edgeToRemove.id});
        AddContents(faceToUpdate.Eids, std::vector<size_t>{edgeToKeep.id});

        Edge& edgeToupdate = ce.Edges.at(diffEdges.at(0));
        if (!edgeToupdate.Vids.empty()) {
            edgeToupdate.Vids.at(0) == vertexToRemove.id ? edgeToupdate.Vids.at(0) = vertexToKeep.id : edgeToupdate.Vids.at(1) = vertexToKeep.id;
        }
    }

    f.N_Fids.clear();
}