/*
* ChordExtractor.h
*
*  Created on: January 18, 2022
*      Author: https://github.com/naeem014
*/

#ifndef _CHORD_EXTRACTOR_H_
#define _CHORD_EXTRACTOR_H_

#include <glm/glm.hpp>

#include "Mesh.h"
#include "MeshUtil.h"

struct Chord {
    int id;
    std::vector<size_t> Edges;
    std::unordered_map<size_t, std::vector<size_t>> eIds; 
};

class ChordExtractor {
    public:
        ChordExtractor();
        ChordExtractor(Mesh& mesh_);
        ChordExtractor(const ChordExtractor& r);
        ~ChordExtractor();

        void Extract();
        void SetVertices();
        void BuildChord(Edge& start);
        std::vector<size_t> GetParallelEdges(Edge& e);
        std::vector<size_t> SelectChords();
        void InititateVisitedEdges();
        bool isChordValid(size_t id);
        void Write(std::vector<size_t> chords = {});

        std::vector<Chord> Chords;
        std::vector<Vertex> Vertices;
        std::vector<Edge> Edges;
    private:
        Mesh& mesh = Mesh();
        std::unordered_map<size_t, size_t> edgeMap;
        std::vector<bool> isEdgeVisited;
};

#endif