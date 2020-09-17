/*
 * EdgeRotateSimplifier.h
 *
 *  Created on: Jan 20, 2019
 *      Author: cotrik
 */

#ifndef LIBCOTRIK_SRC_EDGE_ROTATE_SIMPLIFIER_H_
#define LIBCOTRIK_SRC_EDGE_ROTATE_SIMPLIFIER_H_

#include "Simplifier.h"
class EdgeRotateSimplifier : public Simplifier {
public:
    EdgeRotateSimplifier(Mesh& mesh);
    virtual ~EdgeRotateSimplifier();
private:
    EdgeRotateSimplifier();
    EdgeRotateSimplifier(const EdgeRotateSimplifier&);
    EdgeRotateSimplifier& operator = (const EdgeRotateSimplifier&);
public:
    void Run();
    void Run(std::set<size_t>& canceledFids);
    void RunCollective(std::set<size_t>& canceledFids);

    void Combine(std::vector<size_t>& com, std::vector<std::vector<size_t>> &res, int n, int k, int start);
    std::vector<std::vector<size_t>> Combine(int n, int k);
    size_t GetDiagnalVid(size_t vid, size_t fid);
    double GetAngle(const Vertex& v, const Vertex& v0, const Vertex& v1);
    bool IsConvex(const Vertex& v, const std::vector<size_t>& fids);
    std::map<size_t, std::set<size_t>> GetPatchid_Fids();
    std::map<size_t, std::set<size_t>> GetPatchid_Vids(const std::map<size_t, std::set<size_t>>& patchid_fids);
    std::vector<size_t> GetNeighborFids(const Vertex& v, const std::set<size_t>& patch_fids);
    std::vector<size_t> GetNeighborVids(const Vertex& v, size_t fid);
    std::vector<size_t> GetBoundaryEids(const Face& f0, const Face& f1, const Edge& exclude_e);
    std::vector<size_t> GetRotateVids(const std::vector<size_t>& boundary_eids, size_t diag_vid0, size_t diag_vid1);
    std::set<size_t> GetRotateEids(const Vertex& v, const std::vector<size_t>& fids);
    std::set<size_t> GetRotateEids();
    std::set<size_t> GetRotateFids();
    void Rotate(std::set<size_t>& canceledFids);
    void Rotate(const Edge& e, const Vertex& v, std::set<size_t>& canceledFids);
};

#endif /* LIBCOTRIK_SRC_EDGE_ROTATE_SIMPLIFIER_H_ */
