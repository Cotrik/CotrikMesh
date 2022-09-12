#include "ChordExtractor.h"

ChordExtractor::ChordExtractor() {}
ChordExtractor::ChordExtractor(Mesh& mesh_) : mesh(mesh_) {}
ChordExtractor::ChordExtractor(const ChordExtractor& r) {
    std::cout << "Setting extractor" << std::endl;
    Vertices = r.Vertices;
    Edges = r.Edges;
    Chords = r.Chords;
    mesh = r.mesh;
}
ChordExtractor::~ChordExtractor() {}

void ChordExtractor::Extract() {
    SetVertices();
    InititateVisitedEdges();
    for (auto& e: mesh.E) {
        if (isEdgeVisited[e.id] || e.Vids.empty()) continue;
        BuildChord(e);
    }
    std::cout << "Extracted " << Chords.size() << " chords" << std::endl;
}

void ChordExtractor::SetVertices() {
    for (auto& e: mesh.E) {
        if (e.Vids.empty()) {
            Vertex v = glm::dvec3(0.0, 0.0, 0.0);
            Vertices.push_back(v);
            continue;
        }
        glm::dvec3 coords = 0.5 * (mesh.V.at(e.Vids.at(0)).xyz() + mesh.V.at(e.Vids.at(1)).xyz());
        Vertex v = coords;
        v.id = e.id;
        Vertices.push_back(v);
    }
}

void ChordExtractor::BuildChord(Edge& start) {
    Chord chord;
    std::queue<size_t> q;
    q.push(start.id);
    while (!q.empty()) {
        size_t id = q.front();
        q.pop();
        if (isEdgeVisited[id]) continue;
        Edge& currentEdge = mesh.E.at(id);
        isEdgeVisited[id] = true;
        std::vector<size_t> parallelEs = GetParallelEdges(currentEdge);
        if (parallelEs.empty()) continue;
        for (auto eid: parallelEs) {
            q.push(eid);
            Edge& parallelEdge = mesh.E.at(eid);
            Edge newEdge(std::vector<size_t>{currentEdge.id, parallelEdge.id});
            newEdge.id = Edges.size();
            Edges.push_back(newEdge);
            Vertices.at(currentEdge.id).N_Eids.push_back(newEdge.id);
            Vertices.at(parallelEdge.id).N_Eids.push_back(newEdge.id);
            chord.Edges.push_back(newEdge.id);
        }
    }
    if (!chord.Edges.empty()) {
        chord.id = Chords.size();
        Chords.push_back(chord);
    }
}

std::vector<size_t> ChordExtractor::SelectChords() {
    std::vector<size_t> res;
    for (int i = 0; i < Chords.size(); i++) {
        // if (isChordValid(i)) res.push_back(i);
        res.push_back(i);
    }
    return res;
}

bool ChordExtractor::isChordValid(size_t id) {
    bool res = true;
    auto& chord = Chords.at(id);
    for (auto eId: chord.Edges) {
        Edge& e = Edges.at(eId);
        for (auto eid: e.Vids) {
            Edge& meshEdge = mesh.E.at(eid);
            auto& v0 = mesh.V.at(meshEdge.Vids[0]);
            auto& v1 = mesh.V.at(meshEdge.Vids[1]);
            int count = 0;
            if (v0.isCorner) ++count;
            if (v1.isCorner) ++count;
            if (count > 1) res = false;
            std::set<size_t> labels;
            if (!v0.isCorner && v0.label != MAXID) labels.insert(v0.label);
            if (!v1.isCorner && v1.label != MAXID) labels.insert(v1.label);
            if (!v0.isBoundary && !v1.isBoundary && v0.N_Fids.size() + v1.N_Fids.size() - 4 > 6)
                res = false;
            if (!v0.isBoundary && !v1.isBoundary && v0.N_Fids.size() + v1.N_Fids.size() - 4 < 2/*Simplifier::minValence*/)
                res = false;
            if (v0.type == FEATURE && v1.type == FEATURE && v0.N_Fids.size() + v1.N_Fids.size() != 4/*Simplifier::minValence*/)
                res = false;
            if (v0.type == CORNER && v1.type == FEATURE && v0.N_Fids.size() + v1.N_Fids.size() != 2 + v0.idealValence/*Simplifier::minValence*/)
                res = false;
            if (v1.type == CORNER && v0.type == FEATURE && v0.N_Fids.size() + v1.N_Fids.size() != 2 + v1.idealValence/*Simplifier::minValence*/)
                res = false;

            // if (v0.isBoundary && !v1.isBoundary && v0.N_Fids.size() + v1.N_Fids.size() - 4 > 5)
            //     res = false;
            // if (!v0.isBoundary && v1.isBoundary && v0.N_Fids.size() + v1.N_Fids.size() - 4 > 5)
            //     res = false;
            // if (v0.isBoundary && !v1.isBoundary && v0.N_Fids.size() + v1.N_Fids.size() - 4 < 3)
            //     res = false;
            // if (!v0.isBoundary && v1.isBoundary && v0.N_Fids.size() + v1.N_Fids.size() - 4 < 3)
            //     res = false;

            // if ((v0.idealValence >= 3 || v1.idealValence >= 3) && v0.isBoundary && v1.isBoundary &&
            //         v0.N_Fids.size() + v1.N_Fids.size() - 2 < 3)
            //     res = false;

            if (labels.size() == 1 && v0.type == CORNER && v1.type == FEATURE && v0.labels.find(v1.label) == v0.labels.end())
                res = false;
            if (labels.size() == 1 && v0.type == FEATURE && v1.type == CORNER && v1.labels.find(v0.label) == v1.labels.end())
                res = false;
            if (labels.size() > 2)
                res = false;
            if (labels.size() == 2 && (v0.isCorner || v1.isCorner))
                res = false;

            if (!res) break;
        }
        if (!res) break;
    }
    return res;
}

std::vector<size_t> ChordExtractor::GetParallelEdges(Edge& e) {
    std::vector<size_t> res;
    for (auto nfid: e.N_Fids) {
        Face& f = mesh.F.at(nfid);
        for (auto faceEid: f.Eids) {
            if (faceEid == e.id) continue;
            if (isEdgeVisited[faceEid]) continue;
            Edge& faceEdge = mesh.E.at(faceEid);
            if (std::find(faceEdge.Vids.begin(), faceEdge.Vids.end(), e.Vids.at(0)) == faceEdge.Vids.end() &&
                std::find(faceEdge.Vids.begin(), faceEdge.Vids.end(), e.Vids.at(1)) == faceEdge.Vids.end()) {
                    res.push_back(faceEid);
            }
        }
    }
    return res;
}

void ChordExtractor::InititateVisitedEdges() {
    isEdgeVisited.clear();
    isEdgeVisited.resize(mesh.E.size(), false);
}

void ChordExtractor::Write(std::vector<size_t> chords) {
    std::cout << "Writing output file" << std::endl;
    std::string outputf = "chords.vtk";
    std::ofstream ofs(outputf.c_str());
    ofs << "# vtk DataFile Version 3.0\n"
        << outputf.c_str() << "\n"
        << "ASCII\n\n"
        << "DATASET UNSTRUCTURED_GRID\n";
    ofs << "POINTS " << Vertices.size() << " double\n";
    std::vector<Edge> c_indices;
    if (chords.empty()) {
        std::cout << "writing all extracted chords" << std::endl;
        for (auto& chord: Chords) {
            for (auto& eid: chord.Edges) {
                c_indices.push_back(Edges.at(eid));
            }
        }
    } else {
        std::cout << "writing " << chords.size() << " chords" << std::endl;
        for (auto id: chords) {
            auto& chord = Chords.at(id);
            for (auto& eid: chord.Edges) {
                c_indices.push_back(Edges.at(eid));
            }
        }
    }
    // std::vector<size_t> c_indices = {12, 296};
    // std::cout << c_indices.size() << std::endl;
    for (size_t i = 0; i < Vertices.size(); i++) {
        ofs << std::fixed << std::setprecision(7) <<  Vertices.at(i).x << " " <<  Vertices.at(i).y << " " <<  Vertices.at(i).z << "\n";
    }
    ofs << "CELLS " << c_indices.size() << " " << 3 * c_indices.size() << std::endl;
    for (size_t i = 0; i < c_indices.size(); i++) {
        ofs << "2 " << c_indices.at(i).Vids.at(0) << " " << c_indices.at(i).Vids.at(1) << std::endl;
    }
    ofs << "CELL_TYPES " << c_indices.size() << "\n";
    for (size_t i = 0; i < c_indices.size(); i++) {
        ofs << "3" << std::endl;
    }
}