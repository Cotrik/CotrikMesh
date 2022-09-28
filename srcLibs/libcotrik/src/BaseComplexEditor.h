/*
 * BaseComplexEditor.h
 *
 *  Created on: Jul 10, 2017
 *      Author: cotrik
 */

#ifndef LIBCOTRIK_SRC_BASECOMPLEXEDITOR_H_
#define LIBCOTRIK_SRC_BASECOMPLEXEDITOR_H_

#include "BaseComplex.h"
#include <unordered_map>
#include <unordered_set>

class BaseComplexEditor
{
public:
    BaseComplexEditor(Mesh& mesh, BaseComplex& baseComplex);
    BaseComplexEditor(const BaseComplexEditor& rhs);
    virtual ~BaseComplexEditor();
private:
    BaseComplexEditor();
public:
    void Run();
    /*
     * MoveSingularity() will move a singularity to one of its neighboring parallel base Complex Edges.
     * At the mean time the separated face patches which are connected to other neighboring parallel base Complex Edges
     * will be removed. This operation will used to align singularities
     * */
    void MoveSingularity(SingularE& singularEdge, const size_t neigborId);
    void EraseSeparatedFaces(SingularE& singularEdge);
    void Remove(SingularE& singularEdge);
    void RemoveOnePadding();
    void CollapseComponent();
    void CollapseComponent(SingularE& singularEdge);
    std::vector<size_t> GetTargetNeighborComponentCidesGroup(const SingularE& singularEdge, const std::vector<size_t>& separatedFacePatchIds) const;
    std::vector<std::vector<std::vector<size_t>>> GetStructureVids(const SingularE& singularEdge,
            const std::vector<size_t>& separatedFacePatchIds, const std::vector<size_t>& targetNeighborComponentCidsGroup) const;
    void GetStructureIdSets(const std::vector<size_t>& targetNeighborComponentCidsGroup,
            std::unordered_set<size_t>& vids, std::unordered_set<size_t>& eids, std::unordered_set<size_t>& fids, std::unordered_set<size_t>& cids) const;
    std::vector<std::vector<std::vector<size_t>>> GetStructureVids(const SingularE& singularEdge,
            std::unordered_set<size_t>& vids, std::unordered_set<size_t>& eids, std::unordered_set<size_t>& fids, std::unordered_set<size_t>& cids) const;
    bool CanNeighborCompoenntsCollapsed(SingularE& singularEdge, std::vector<size_t>& separatedFacePatchIds);
    bool GetNext(size_t& vid, int keys[], int& val, std::unordered_set<size_t>& searchVids,
            std::unordered_map<size_t, size_t>& key_val, std::vector<bool>& visited);
    void ModifyIdOfSeparatedFacePatch(SingularE& singularEdge, std::vector<size_t>& separatedFacePatchIds, std::unordered_map<size_t, size_t>& key_val);
    void WriteBaseComplex_ColorFacesVTK(const char *filename) const;
    void WriteBaseComplex(const char *filename) const;
    void WriteMeshVTK(const char *filename) const;
    Mesh& GetMesh() const;
    BaseComplex& GetBaseComplex() const;
private:
    bool IsOnBoundary(const SingularE& singularEdge) const;
    bool IsConnectedToBoundary(const SingularE& singularEdge) const;
    size_t GetOppositeComponentEdgeId(const ComponentCell& component, size_t componentFaceId1, size_t componentFaceId2);
    bool isAligned(SingularE& singularEdge, const std::vector<size_t>& targetNeighborComponentCidsGroup, size_t oppositeEdgeFrontVid, size_t oppositeEdgeBackVid);
private:
    Mesh& mesh;
    BaseComplex& m_baseComplex;
};

#endif /* LIBCOTRIK_SRC_BASECOMPLEXEDITOR_H_ */
