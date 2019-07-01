/*
 * BaseComplexQuad.h
 *
 *  Created on: Jan 11, 2018
 *      Author: cotrik
 */

#ifndef LIBCOTRIK_SRC_BASECOMPLEXQUAD_H_
#define LIBCOTRIK_SRC_BASECOMPLEXQUAD_H_

#include "BaseComplex.h"

class BaseComplexQuad : public BaseComplex
{
public:
    BaseComplexQuad(Mesh& mesh);
    BaseComplexQuad(const BaseComplexQuad& component);
    virtual ~BaseComplexQuad();
private:
    BaseComplexQuad();

public:
    const Mesh& GetMesh() const;
    void WriteSingularities_VTK(const char *filename) const;
    void WriteSingularV_VTK(const char* filename) const;
    void WriteSingularE_VTK(const char *filename) const;
    void WriteBaseComplex_VTK(const char *filename) const;
    void WriteBaseComplexQuadVTK(const char *filename) const;
    void WriteBaseComplex_ColorVerticesVTK(const char *filename) const;
    void WriteBaseComplex_ColorEdgesVTK(const char *filename) const;
    void WriteBaseComplex_ColorFacesVTK(const char *filename) const;
    void WriteBaseComplexComponentsVTK(const char *filename) const;
    void WriteBaseComplexAllComponentsVTK(const char *filename_prefix) const;
    void WriteBaseComplexComponentsVTK(const char *filename_prefix, const size_t id) const;
    void WriteBaseComplexComponentsWithoutSingularitiesVTK(const char *filename) const;
    void WriteBaseComplexComponentsWithSingularitiesVTK(const char *filename) const;
    void WriteBaseComplexSeparatedEdgeLinksVTK(const char *filename) const;
    // Write Neighboring information
    void WriteComponentEdge_NeighborComponentFaces_VTK(const char *filename) const;
    void WriteComponentEdge_NeighborComponentCells_VTK(const char *filename) const;
    void WriteComponentFace_NeighborComponentCells_VTK(const char *filename) const;
    void WriteSingularEdge_NeighborSeparatedFacePatches_VTK(const char *filename) const;
    void WriteSingularEdge_NeighborSeparatedComponentFacePatches_VTK(const char *filename) const;
    void WriteSingularEdge_NeighborSeparatedComponentFacePatches_VTK(const char *filename, const int singularEdgeId) const;
    void WriteAllSingularEdge_NeighborSeparatedComponentFacePatches_VTK(const char *filename_prefix) const;

    void Build();

public:
    virtual void BuildV();
    virtual void BuildE();
    virtual void BuildF();

    virtual void BuildComponentV();
    virtual void BuildComponentE();
    virtual void BuildComponentF();
    virtual void BuildComponentC();

    virtual void BuildComponentColor();
    virtual void BuildComponentConnectivities();
    virtual void BuildSigularEdge_separatedPatches();
    virtual void BuildSigularEdge_separatedFacePatches();
//    virtual void BuildSigularEdge_separatedEdgePatches();
    virtual void BuildSigularEdge_separatedComponentFacePatches();
    virtual void BuildSigularEdge_separatedComponentEdgePatches();
    std::vector<char> GetColorsOfNeighborComponents(const ComponentCell& component);


public:
    size_t GetNextSingularEdgeId(const size_t vid, const size_t cur_edge_id);
    // vid is the starting point, eid tells the direction; circular edge link returns true;
    bool TraceEdge(const size_t vid, const size_t eid, std::vector<size_t>& Vids, std::vector<size_t>& Eids, std::vector<bool>& is_edge_visited);
    void AddSingularV(const size_t current_vid, const size_t singularVid, size_t& ending_singular_vid);
    void AddSingularE(const std::vector<size_t>& leftVids, const std::vector<size_t>& leftEids,
            const std::vector<size_t>& rightVids, const std::vector<size_t>& rightEids,
            const size_t left_singular_vid, const size_t right_singular_vid, const size_t singularEid);
    void AddcircularSingularE(const std::vector<size_t>& leftVids, const std::vector<size_t>& leftEids, const size_t singularEid);
    bool straight_line_test(const size_t h_v1, const size_t h_v2, const size_t h_v3);
    void AddcircularSingularE(const std::vector<SingularE> &circle_ses);
    void ExtractSingularVandE();
    void BuildSingularityConnectivity();
    void TraceEdge(const std::vector<size_t>& mesh_edge_ids_on_singular_edge, const Face& start_face, std::vector<bool>& is_mesh_edge_visited);
    void TraceFace(const std::vector<size_t>& mesh_edge_ids_on_singular_edge, const Face& face, std::vector<bool>& is_mesh_face_visited);
    void TraceFace(const Face& start_face, std::vector<bool>& is_mesh_face_visited);
    void TraceFace(const Face& start_face, std::vector<bool>& is_mesh_face_visited, std::vector<size_t>& fids);
    bool IsFaceMetSingularNode(const Face& face);
    bool IsEdgeMetSingularNode(const Edge& edge);
    bool IsEdgeMetSingularEdge(const Edge& edge);
    bool IsEdgeMetBaseComplexVertex(const Edge& edge);
    bool IsEdgeOnBaseComplexEdge(const Edge& edge);
    bool IsFaceOnBaseComplexFace(const Face& face);
    inline bool IsOnBaseComplex(const Vertex& vertex) const;
    inline bool IsOnBaseComplex(const Edge& edge) const;
    inline bool IsOnBaseComplex(const Face& face) const;
    size_t GetOppositEdgeId(const Face& face, const size_t edgeid);
    size_t GetOppositFaceId(const Cell& cell, const size_t faceid);

    size_t TraceVertex(const Vertex& start_vertex, const Edge& start_edge,
            std::vector<bool>& is_edge_visited, std::vector<size_t>& vids_link, std::vector<size_t> &eids_link); // return end_mesh_vertex_id
    size_t TraceEdge(const std::vector<size_t>& start_eids_link, const Face& start_face,
            std::vector<bool>& is_face_visited, std::vector<size_t>& fids_link); // return end_component_edge_id
    size_t TraceEdge(const Edge& start_edge, const Face& start_face, std::vector<bool>& is_face_visited, std::vector<size_t>& fids_patch);
    size_t TraceFace(const std::vector<size_t>& start_fids_patch, const Cell& start_cell,
            std::vector<bool>& is_cell_visited, std::vector<size_t>& cids_patch); // return end_component_face_id
    size_t TraceFace(const Face& start_face, const Cell& start_cell, std::vector<bool>& is_cell_visited, std::vector<size_t>& cids_patch); // return end_face_id
    Face* GetNextFace(const Edge* currentEdge, const Face* currentFace, const std::vector<bool>& is_mesh_edge_visited);
    Face* GetNextFace(const Edge& currentEdge, const Face& currentFace);
    const Face& GetNextFace(const Edge& currentEdge, const Face& currentFace, const Edge& nextEdge);
    const Cell& GetNextCell(const Face& currentFace, const Cell& currentCell, const Face& nextFace);

    Edge* GetNextEdge(const Vertex& currentVertex, const Edge& currentEdge);

    size_t RegionGrowCell(const Cell& start_cell, std::vector<bool>& is_cell_visited, std::vector<size_t>& cids);
    size_t RegionGrowFace(const Face& start_face, std::vector<bool>& is_face_visited, std::vector<size_t>& fids);
    size_t RegionGrowEdge(const Edge& start_edge, std::vector<bool>& is_edge_visited, std::vector<size_t>& eids);

    std::vector<size_t> GetComponentFids(const ComponentCell& component) const;
    std::vector<size_t> GetComponentEids(const ComponentCell& component) const;
    std::vector<size_t> GetComponentVids(const ComponentCell& component) const;
    std::vector<size_t> GetComponentEids(const ComponentFace& componentFace) const;
    std::vector<size_t> GetComponentVids(const ComponentFace& componentFace) const;

    std::vector<size_t> GetNeighborComponentCids(const ComponentFace& componentFace) const;
    std::vector<size_t> GetNeighborComponentCids(const ComponentEdge& componentEdge) const;
    std::vector<size_t> GetNeighborComponentFids(const ComponentEdge& componentEdge) const;

    const std::vector<size_t>& GetComponentEidsLink(SingularE& singularEdge) const;
    const std::vector<std::vector<size_t> >& GetNeighborComponentCidsGroups(SingularE& singularEdge) const;
    const std::vector<std::vector<size_t> >& GetNeighborComponentFidsGroups(SingularE& singularEdge) const;

    std::vector<size_t> GetNeighborComponentCids(const SingularE& singularEdge) const;
    std::vector<size_t> GetNeighborComponentFids(const SingularE& singularEdge) const;

    void BuildComponentEdgeLink(ComponentEdge& componentEdge);
    size_t TraceVertex(ComponentEdge& componentEdge, const Vertex& start_vertex, const Edge& start_edge,
            std::vector<bool>& is_edge_visited, std::vector<size_t>& vids_link, std::vector<size_t> &eids_link); // return end_mesh_vertex_id

    void AlignComponentE();
public:
    Edge* GetNextContinuousEdge(const Vertex& vertex, const Edge& edge) const;
    size_t TraceAlongEdge(const Vertex& start_vertex, const Edge& start_edge,
            std::vector<bool>& is_edge_visited, std::vector<size_t>& vids_link, std::vector<size_t> &eids_link); // return end_mesh_vertex_id
    std::vector<char> GetColorsOfNeighborComponents(const ComponentFace& component);
public:
//    Mesh& mesh;
//    std::vector<size_t> Vids;
//    std::vector<size_t> Eids;
//    std::vector<size_t> Fids;
//    std::vector<size_t> Cids;
//
//    std::unordered_set<size_t> hashVids;
//    std::unordered_set<size_t> hashEids;
//    std::unordered_set<size_t> hashFids;
//    std::unordered_set<size_t> hashCids;
//
//    std::vector<std::vector<size_t> > separatedFacePatches;
//    std::vector<std::vector<size_t> > separatedEdgePatches;
    std::vector<std::vector<size_t> > separatedEdgeIdsLink;
    std::vector<std::vector<size_t> > separatedVertexIdsLink;
//    std::vector<std::vector<size_t> > separatedComponentFacePatches;
//    std::vector<std::vector<size_t> > separatedComponentEdgePatches;
//
//    std::vector<ComponentVertex> componentV;
//    std::vector<ComponentEdge> componentE;
//    std::vector<ComponentFace> componentF;
//    std::vector<ComponentCell> componentC;
//    std::vector<Component> components;
//    SingularityInfo SingularityI;
};

#endif /* LIBCOTRIK_SRC_BASECOMPLEX_H_ */
