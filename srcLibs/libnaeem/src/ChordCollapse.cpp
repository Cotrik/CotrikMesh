#include "ChordCollapse.h"
#include <cstdlib>

ChordCollapse::ChordCollapse() {}
ChordCollapse::ChordCollapse(Mesh& mesh_, MeshUtil& mu_, Smoother& smoother_, ChordExtractor& ce_, size_t chordId_) : SimplificationOperation(mesh_, mu_, smoother_), ce(&ce_) {
    chord = &ce->Chords.at(chordId_);
}
ChordCollapse::~ChordCollapse() {}

void ChordCollapse::SetRanking(glm::dvec3 d) {
    CheckValidity();
    if (ce == NULL) {
        std::cout << "No Chord Extractor initialized for Chord Collapse." << std::endl;
        exit(0);
    }
}

bool ChordCollapse::IsOperationValid() {
    CheckValidity();

    // return ce->isChordValid(chord->id);
    bool isValid = true;
    for (auto eId: chord->Edges) {
        Edge& e = ce->Edges.at(eId);
        for (auto eid: e.Vids) {
            Edge& meshEdge = mesh->E.at(eid);
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
    for (auto& eid: chord->Edges) {
        Edge& e = ce->Edges.at(eid);
        if (e.Vids.empty()) continue;
        Edge& e1 = mesh->E.at(e.Vids.at(0));
        Edge& e2 = mesh->E.at(e.Vids.at(1));
        size_t fid = GetIntersection(e1.N_Fids, e2.N_Fids).at(0);
        std::vector<size_t> diffEdges = GetDifference(mesh->F.at(fid).Eids, std::vector<size_t>{e1.id, e2.id});
        size_t edgeToKeepId = diffEdges.at(0);
        size_t edgeToRemoveId = diffEdges.at(1);
        if (mesh->E.at(edgeToRemoveId).isBoundary) {
            edgeToKeepId = diffEdges.at(1);
            edgeToRemoveId = diffEdges.at(0);
        }
        
        CollapseEdge(e1, mesh->E.at(edgeToRemoveId));
        CollapseEdge(e2, mesh->E.at(edgeToRemoveId));
        UpdateEdges(fid, mesh->E.at(edgeToKeepId), mesh->E.at(edgeToRemoveId));
    }
    for (auto& eid: chord->Edges) {
        Edge& e = ce->Edges.at(eid);
        if (e.Vids.empty()) continue;
        Edge& e1 = mesh->E.at(e.Vids.at(0));
        Edge& e2 = mesh->E.at(e.Vids.at(1));
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
    Vertex& source = mesh->V.at(source_id);
    Vertex& target = mesh->V.at(target_id);
    // if (!(target.isBoundary ^ source.isBoundary)) {
    //     target = 0.5 * (source.xyz() + target.xyz());
    // } else if (!(target.type == FEATURE ^ source.type == FEATURE)) {
    //     target = 0.5 * (source.xyz() + target.xyz());
    // }
    
    UpdateContents(target.N_Vids, std::vector<size_t>{source.id});
    AddContents(target.N_Vids, GetDifference(source.N_Vids, std::vector<size_t>{target.id}));

    UpdateContents(target.N_Eids, std::vector<size_t>{e.id});
    AddContents(target.N_Eids, GetDifference(source.N_Eids, std::vector<size_t>{e.id}));

    UpdateContents(target.N_Fids, e.N_Fids);
    AddContents(target.N_Fids, GetDifference(source.N_Fids, e.N_Fids));

    for (auto nvid: source.N_Vids) {
        if (nvid == target.id) continue;
        auto& v = mesh->V.at(nvid);
        UpdateContents(v.N_Vids, std::vector<size_t>{source.id});
        AddContents(v.N_Vids, std::vector<size_t>{target.id});
    }

    std::vector<size_t> eids = GetDifference(source.N_Eids, std::vector<size_t>{e.id});
    for (auto neid: eids) {
        Edge& ne = mesh->E.at(neid);
        ne.Vids.at(0) == source.id ? ne.Vids.at(0) = target.id : ne.Vids.at(1) = target.id;
    }

    std::vector<size_t> fids = GetDifference(source.N_Fids, e.N_Fids);
    std::vector<size_t> tfids = GetDifference(target.N_Fids, e.N_Fids);
    for (auto nfid: fids) {
        Face& nf = mesh->F.at(nfid);
        int idx = std::distance(nf.Vids.begin(), std::find(nf.Vids.begin(), nf.Vids.end(), source.id));
        nf.Vids.at(idx) = target.id;
        UpdateContents(nf.N_Fids, e.N_Fids);
        AddContents(nf.N_Fids, tfids);
    }
    for (auto nfid: tfids) {
        Face& nf = mesh->F.at(nfid);
        UpdateContents(nf.N_Fids, e.N_Fids);
        AddContents(nf.N_Fids, fids);
    }

    source.N_Vids.clear();
    source.N_Eids.clear();
    source.N_Fids.clear();

    e.Vids.clear();
}

void ChordCollapse::UpdateEdges(size_t fid, Edge& edgeToKeep, Edge& edgeToRemove) {

    Face& f = mesh->F.at(fid);
    Vertex& vertexToKeep = ce->Vertices.at(edgeToKeep.id);
    Vertex& vertexToRemove = ce->Vertices.at(edgeToRemove.id);

    std::vector<size_t> commonEdges = GetIntersection(vertexToKeep.N_Eids, vertexToRemove.N_Eids);
    std::vector<size_t> diffEdges = GetDifference(vertexToRemove.N_Eids, vertexToKeep.N_Eids);
    ce->Edges.at(commonEdges.at(0)).Vids.clear();
    UpdateContents(vertexToKeep.N_Eids, commonEdges);
    AddContents(vertexToKeep.N_Eids, diffEdges);
    UpdateContents(edgeToKeep.N_Fids, std::vector<size_t>{fid});
    std::vector<size_t> facesToUpdate = GetDifference(edgeToRemove.N_Fids, std::vector<size_t>{fid});
    if (!facesToUpdate.empty()) {
        AddContents(edgeToKeep.N_Fids, facesToUpdate);
        Face& faceToUpdate = mesh->F.at(facesToUpdate.at(0));
        UpdateContents(faceToUpdate.Eids, std::vector<size_t>{edgeToRemove.id});
        AddContents(faceToUpdate.Eids, std::vector<size_t>{edgeToKeep.id});

        Edge& edgeToupdate = ce->Edges.at(diffEdges.at(0));
        if (!edgeToupdate.Vids.empty()) {
            edgeToupdate.Vids.at(0) == vertexToRemove.id ? edgeToupdate.Vids.at(0) = vertexToKeep.id : edgeToupdate.Vids.at(1) = vertexToKeep.id;
        }
    }

    f.N_Fids.clear();
}

double ChordCollapse::CalculateRanking() {
    double alpha_q = 0.05;
    double alpha_d = 0.05;
    double alpha_v = 0.9;
    double E_q = 0.0;
    double E_d = 0.0;
    double E_v = 0.0;

    std::set<size_t> edges;
    for (auto& eid: chord->Edges) {
        Edge& e = ce->Edges.at(eid);
        edges.insert(e.Vids.at(0));
        edges.insert(e.Vids.at(1));
    }
    double it = 0.0;
    double max_qem = -1.0;
    double max_distance = -1.0;
    double max_val = -1.0;
    double avg_val = 0.0;
    for (auto eid: edges) {
        Edge& e = mesh->E.at(eid);
        glm::dvec4 newV = CalculateQEM(e.Vids.at(0), e.Vids.at(1));
        if (newV[3] > max_qem) max_qem = newV[3];
        double d = CalculateAreaDistance(e.Vids.at(0), e.Vids.at(1));
        if (d > max_distance) max_distance = d;
        double val = CalculateValenceTerm(e.Vids.at(0), e.Vids.at(1));
        if (val > max_val) max_val = val;
        avg_val += val;
        it += 1.0;
    }
    E_q = max_qem;
    E_d = max_distance;
    E_v = max_val + (avg_val / it);
    double r = (alpha_q * (1 - exp(-E_q))) + (alpha_d * (1 - exp(-E_d))) + (alpha_v * (1 - exp(-E_v)));
    std::cout << "ranking: " << r << std::endl;
    return r;
}

glm::dvec4 ChordCollapse::CalculateQEM(size_t v1_id, size_t v2_id) {
    auto& v1 = mesh->V.at(v1_id);
    auto& v2 = mesh->V.at(v2_id);

    // std::cout << "Vertex v1: " << v1.id << " " << v1.x << " " << v1.y << " " << v1.z << std::endl;
    // std::cout << "Vertex v2: " << v2.id << " " << v2.x << " " << v2.y << " " << v2.z << std::endl;

    glm::dmat4 Q1(0.0);
    // std::cout << "Calculating Q1" << std::endl;
    for (auto fid: v1.N_Fids) {
        auto& f = mesh->F.at(fid);
        // std::cout << "neighbor face: " << f.id << " " << f.Vids.at(0) << " " << f.Vids.at(1) << " " << f.Vids.at(2) << " " << f.Vids.at(3) << std::endl;
        int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v1.id));
        auto& v_a = mesh->V.at(f.Vids.at((idx+1)%f.Vids.size()));
        auto& v_b = mesh->V.at(f.Vids.at((idx+3)%f.Vids.size()));

        // std::cout << "Vertex v_a: " << v_a.id << " " << v_a.x << " " << v_a.y << " " << v_a.z << std::endl;
        // std::cout << "Vertex v_b: " << v_b.id << " " << v_b.x << " " << v_b.y << " " << v_b.z << std::endl;
        
        glm::dvec3 A = v_a.xyz() - v1.xyz();
        glm::dvec3 B = v_b.xyz() - v1.xyz();
        // std::cout << "A: " << A.x << " " << A.y << " " << A.z << std::endl;
        // std::cout << "B: " << B.x << " " << B.y << " " << B.z << std::endl;
        
        glm::dvec3 n = glm::normalize(glm::cross(A, B));
        // std::cout << "normal: " << n.x << " " << n.y << " " << n.z << std::endl;

        glm::dvec4 p(n, glm::dot(n,v1.xyz()));
        // std::cout << "plane: " << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << std::endl;
        
        glm::dmat4 Q;
        Q[0] = glm::dvec4(p[0]*p[0], p[0]*p[1], p[0]*p[2], p[0]*p[3]);
        Q[1] = glm::dvec4(p[1]*p[0], p[1]*p[1], p[1]*p[2], p[1]*p[3]);
        Q[2] = glm::dvec4(p[2]*p[0], p[2]*p[1], p[2]*p[2], p[2]*p[3]);
        Q[3] = glm::dvec4(p[3]*p[0], p[3]*p[1], p[3]*p[2], p[3]*p[3]);

        // std::cout << "Q: " << std::endl;
        // std::cout << Q[0][0] << " " << Q[1][0] << " " << Q[2][0] << " " << Q[3][0] << std::endl;
        // std::cout << Q[0][1] << " " << Q[1][1] << " " << Q[2][1] << " " << Q[3][1] << std::endl;
        // std::cout << Q[0][2] << " " << Q[1][2] << " " << Q[2][2] << " " << Q[3][2] << std::endl;
        // std::cout << Q[0][3] << " " << Q[1][3] << " " << Q[2][3] << " " << Q[3][3] << std::endl;

        Q1 += Q;

        // std::cout << "Q1: " << std::endl;
        // std::cout << Q1[0][0] << " " << Q1[1][0] << " " << Q1[2][0] << " " << Q1[3][0] << std::endl;
        // std::cout << Q1[0][1] << " " << Q1[1][1] << " " << Q1[2][1] << " " << Q1[3][1] << std::endl;
        // std::cout << Q1[0][2] << " " << Q1[1][2] << " " << Q1[2][2] << " " << Q1[3][2] << std::endl;
        // std::cout << Q1[0][3] << " " << Q1[1][3] << " " << Q1[2][3] << " " << Q1[3][3] << std::endl;
        // std::cout << "***********************************" << std::endl;
    }
    // std::cout << "Q1: " << std::endl;
    // std::cout << Q1[0][0] << " " << Q1[1][0] << " " << Q1[2][0] << " " << Q1[3][0] << std::endl;
    // std::cout << Q1[0][1] << " " << Q1[1][1] << " " << Q1[2][1] << " " << Q1[3][1] << std::endl;
    // std::cout << Q1[0][2] << " " << Q1[1][2] << " " << Q1[2][2] << " " << Q1[3][2] << std::endl;
    // std::cout << Q1[0][3] << " " << Q1[1][3] << " " << Q1[2][3] << " " << Q1[3][3] << std::endl;


    glm::dmat4 Q2(0.0);
    // std::cout << "Calculating Q2" << std::endl;
    for (auto fid: v2.N_Fids) {
        auto& f = mesh->F.at(fid);
        // std::cout << "neighbor face: " << f.id << " " << f.Vids.at(0) << " " << f.Vids.at(1) << " " << f.Vids.at(2) << " " << f.Vids.at(3) << std::endl;
        int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v2.id));
        auto& v_a = mesh->V.at(f.Vids.at((idx+1)%f.Vids.size()));
        auto& v_b = mesh->V.at(f.Vids.at((idx+3)%f.Vids.size()));

        // std::cout << "Vertex v_a: " << v_a.id << " " << v_a.x << " " << v_a.y << " " << v_a.z << std::endl;
        // std::cout << "Vertex v_b: " << v_b.id << " " << v_b.x << " " << v_b.y << " " << v_b.z << std::endl;
        
        glm::dvec3 A = v_a.xyz() - v2.xyz();
        glm::dvec3 B = v_b.xyz() - v2.xyz();
        // std::cout << "A: " << A.x << " " << A.y << " " << A.z << std::endl;
        // std::cout << "B: " << B.x << " " << B.y << " " << B.z << std::endl;
        
        glm::dvec3 n = glm::normalize(glm::cross(A, B));
        // std::cout << "normal: " << n.x << " " << n.y << " " << n.z << std::endl;
        glm::dvec4 p(n, glm::dot(n,v2.xyz()));
        // std::cout << "plane: " << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << std::endl;
        
        glm::dmat4 Q;
        Q[0] = glm::dvec4(p[0]*p[0], p[0]*p[1], p[0]*p[2], p[0]*p[3]);
        Q[1] = glm::dvec4(p[1]*p[0], p[1]*p[1], p[1]*p[2], p[1]*p[3]);
        Q[2] = glm::dvec4(p[2]*p[0], p[2]*p[1], p[2]*p[2], p[2]*p[3]);
        Q[3] = glm::dvec4(p[3]*p[0], p[3]*p[1], p[3]*p[2], p[3]*p[3]);

        // std::cout << "Q: " << std::endl;
        // std::cout << Q[0][0] << " " << Q[1][0] << " " << Q[2][0] << " " << Q[3][0] << std::endl;
        // std::cout << Q[0][1] << " " << Q[1][1] << " " << Q[2][1] << " " << Q[3][1] << std::endl;
        // std::cout << Q[0][2] << " " << Q[1][2] << " " << Q[2][2] << " " << Q[3][2] << std::endl;
        // std::cout << Q[0][3] << " " << Q[1][3] << " " << Q[2][3] << " " << Q[3][3] << std::endl;

        Q2 += Q;
        // std::cout << "Q2: " << std::endl;
        // std::cout << Q2[0][0] << " " << Q2[1][0] << " " << Q2[2][0] << " " << Q2[3][0] << std::endl;
        // std::cout << Q2[0][1] << " " << Q2[1][1] << " " << Q2[2][1] << " " << Q2[3][1] << std::endl;
        // std::cout << Q2[0][2] << " " << Q2[1][2] << " " << Q2[2][2] << " " << Q2[3][2] << std::endl;
        // std::cout << Q2[0][3] << " " << Q2[1][3] << " " << Q2[2][3] << " " << Q2[3][3] << std::endl;
        // std::cout << "***********************************" << std::endl;
    }

    // std::cout << "Q2: " << std::endl;
    // std::cout << Q2[0][0] << " " << Q2[1][0] << " " << Q2[2][0] << " " << Q2[3][0] << std::endl;
    // std::cout << Q2[0][1] << " " << Q2[1][1] << " " << Q2[2][1] << " " << Q2[3][1] << std::endl;
    // std::cout << Q2[0][2] << " " << Q2[1][2] << " " << Q2[2][2] << " " << Q2[3][2] << std::endl;
    // std::cout << Q2[0][3] << " " << Q2[1][3] << " " << Q2[2][3] << " " << Q2[3][3] << std::endl;

    glm::dmat4 Q_ = Q1+Q2;
    // std::cout << "Q_: " << std::endl;
    // std::cout << Q_[0][0] << " " << Q_[1][0] << " " << Q_[2][0] << " " << Q_[3][0] << std::endl;
    // std::cout << Q_[0][1] << " " << Q_[1][1] << " " << Q_[2][1] << " " << Q_[3][1] << std::endl;
    // std::cout << Q_[0][2] << " " << Q_[1][2] << " " << Q_[2][2] << " " << Q_[3][2] << std::endl;
    // std::cout << Q_[0][3] << " " << Q_[1][3] << " " << Q_[2][3] << " " << Q_[3][3] << std::endl;

    glm::dmat4 Q_v = Q1+Q2;
    Q_v[0][3] = 0.0;
    Q_v[1][3] = 0.0;
    Q_v[2][3] = 0.0;
    Q_v[3][3] = 1.0;
    // std::cout << "Q_v: " << std::endl;
    // std::cout << Q_v[0][0] << " " << Q_v[1][0] << " " << Q_v[2][0] << " " << Q_v[3][0] << std::endl;
    // std::cout << Q_v[0][1] << " " << Q_v[1][1] << " " << Q_v[2][1] << " " << Q_v[3][1] << std::endl;
    // std::cout << Q_v[0][2] << " " << Q_v[1][2] << " " << Q_v[2][2] << " " << Q_v[3][2] << std::endl;
    // std::cout << Q_v[0][3] << " " << Q_v[1][3] << " " << Q_v[2][3] << " " << Q_v[3][3] << std::endl;

    glm::dvec4 v_in(0.0, 0.0, 0.0, 1.0);

    // std::cout << "determinant: " << glm::determinant(Q_v) << std::endl; 
    double determinant = glm::determinant(Q_v);
    glm::dvec4 v_(0.5 * (v1.xyz() + v2.xyz()), 1.0);
    if (fabs(determinant) > 0.0) {
        glm::dmat4 inv_Q_v = glm::inverse(Q_v);
        // std::cout << "inv_Q_v: " << std::endl;
        // std::cout << inv_Q_v[0][0] << " " << inv_Q_v[1][0] << " " << inv_Q_v[2][0] << " " << inv_Q_v[3][0] << std::endl;
        // std::cout << inv_Q_v[0][1] << " " << inv_Q_v[1][1] << " " << inv_Q_v[2][1] << " " << inv_Q_v[3][1] << std::endl;
        // std::cout << inv_Q_v[0][2] << " " << inv_Q_v[1][2] << " " << inv_Q_v[2][2] << " " << inv_Q_v[3][2] << std::endl;
        // std::cout << inv_Q_v[0][3] << " " << inv_Q_v[1][3] << " " << inv_Q_v[2][3] << " " << inv_Q_v[3][3] << std::endl;

        v_ = inv_Q_v*v_in;
    }

    glm::dvec4 v_int = Q_*v_;
    // std::cout << "v_: " << v_[0] << " " << v_[1] << " " << v_[2] << " " << v_[3] << std::endl;
    // std::cout << "v_int: " << v_int[0] << " " << v_int[1] << " " << v_int[2] << " " << v_int[3] << std::endl;

    v_[3] = glm::dot(v_, v_int);
    return v_;
}

double ChordCollapse::CalculateAreaDistance(size_t v1_id, size_t v2_id) {
    return glm::distance(mesh->V.at(v1_id).xyz(), mesh->V.at(v2_id).xyz());
}

double ChordCollapse::CalculateValenceTerm(size_t v1_id, size_t v2_id) {
    auto& v1 = mesh->V.at(v1_id);
    auto& v2 = mesh->V.at(v2_id);
    int target_v = v1.N_Vids.size() + v2.N_Vids.size() - 4;

    int val_v1 = v1.N_Vids.size() >= 4 ? v1.N_Vids.size() - 4 : 4 - v1.N_Vids.size();
    int val_v2 = v2.N_Vids.size() >= 4 ? v2.N_Vids.size() - 4 : 4 - v2.N_Vids.size();
    int target_valence = target_v >= 4 ? target_v - 4 : 4 - target_v;
    double w_val = 0;
    w_val += target_valence > val_v1 ? target_valence - val_v1 : 0.0;
    w_val += target_valence > val_v2 ? target_valence - val_v2 : 0.0;

    return w_val;
}
