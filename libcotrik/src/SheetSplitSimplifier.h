/*
 * SheetSplitSimplifier.h
 *
 *  Created on: Jan 20, 2019
 *      Author: cotrik
 */

#ifndef LIBCOTRIK_SRC_SHEET_SPLIT_SIMPLIFIER_H_
#define LIBCOTRIK_SRC_SHEET_SPLIT_SIMPLIFIER_H_

#include "Simplifier.h"
class SheetSplitSimplifier : public Simplifier {
public:
    SheetSplitSimplifier(Mesh& mesh);
    virtual ~SheetSplitSimplifier();
private:
    SheetSplitSimplifier();
    SheetSplitSimplifier(const SheetSplitSimplifier&);
    SheetSplitSimplifier& operator = (const SheetSplitSimplifier&);
public:
    void Run(std::set<size_t>& canceledFids);
    std::set<size_t> GetAllParallelEdgeIds(const size_t eid, const std::set<size_t>& exclude_eids);
    std::set<size_t> GetAllParallelEdgeIdsFromDoubletVertex(const Vertex& v);
    std::set<size_t> GetCanceledFids(const std::set<size_t>& paralell_eids);
    std::set<size_t> GetOverlapFids(const std::set<size_t>& paralell_eids, const std::set<size_t>& canceledFids);
    void Split(const Vertex& v, const std::set<size_t>& paralell_eids, std::set<size_t>& canceledFids);
    std::map<size_t, size_t> GetInsertEdgeVertices(const std::set<size_t>& paralell_eids);
    std::map<size_t, size_t> GetInsertFaceVertices(const std::set<size_t>& overlapFids);
    std::map<size_t, size_t> GetInsertFaceVertices(const Vertex& v);
    void Insert2Faces(const std::map<size_t, size_t>& paralellEid_newVid, const std::set<size_t>& fids);
    void Insert3Faces(const std::map<size_t, size_t>& paralellEid_newVid, const std::map<size_t, size_t>& doubletFid_newVid,
            const Vertex& v);
    void Insert4Faces(const std::map<size_t, size_t>& paralellEid_newVid, const std::map<size_t, size_t>& overlapFid_newVid,
            const std::set<size_t>& fids);
    void GetOverlapFaces(const std::set<size_t>& paralell_eids, const std::set<size_t>& canceledFids);
private:
    std::vector<size_t> GetNeignborEids(const Vertex& v, const Face& f);
    size_t GetPatchid(const Vertex& v0, const Vertex& v1, const Vertex& v2, const Vertex& v3);
    size_t GetPatchid(const Vertex& v0, const Vertex& v1);
    size_t GetPatchid(const Vertex& v0, const Vertex& v1, Vertex& v);
    size_t GetLabel(const Vertex& v0, const Vertex& v1, Vertex& v);
    std::vector<size_t> GetParallelEids(const Face& f, const std::map<size_t, size_t>& paralellEid_newVid);
    std::vector<size_t> GetPerpendicularEids(const Face& f, const std::map<size_t, size_t>& paralellEid_newVid);
    size_t GetSharedVid(const Edge& e0, const Edge& e1);
    size_t find(const std::map<size_t, size_t>& m, const size_t key, const char* mapName = NULL);
};

#endif /* LIBCOTRIK_SRC_SHEET_SPLIT_SIMPLIFIER_H_ */
