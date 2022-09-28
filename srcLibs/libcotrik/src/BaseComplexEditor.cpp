/*
 * BaseComplexEditor.cpp
 *
 *  Created on: Jul 10, 2017
 *      Author: cotrik
 */

#include "BaseComplexEditor.h"
#include <algorithm>
#include <fstream>
#include <iostream>


std::vector<size_t> singleNumber(const std::vector<size_t>& nums);

BaseComplexEditor::BaseComplexEditor(Mesh& mesh, BaseComplex& baseComplex)
: mesh(mesh), m_baseComplex(baseComplex)
{


}

BaseComplexEditor::BaseComplexEditor(const BaseComplexEditor& rhs)
: mesh(rhs.GetMesh()), m_baseComplex(rhs.GetBaseComplex())
{


}

BaseComplexEditor::~BaseComplexEditor()
{
    // TODO Auto-generated destructor stub
}

Mesh& BaseComplexEditor::GetMesh() const
{
    return mesh;
}

BaseComplex& BaseComplexEditor::GetBaseComplex() const
{
    return m_baseComplex;
}

void BaseComplexEditor::Run()
{
    //----------------------------------------
//    MoveSingularity(m_baseComplex.SingularityI.E.at(24), 1);
//    WriteBaseComplex_ColorFacesVTK("MoveSingularity.vtk");
    //----------------------------------------
    for (auto & singularEdge : m_baseComplex.SingularityI.E) {
        if (singularEdge.separatedFacePatchIds.size() == 4) {
            EraseSeparatedFaces(singularEdge);
        }
    }

    WriteBaseComplex_ColorFacesVTK("EraseSeparatedFaces.vtk");
    //----------------------------------------
//    RemoveOnePadding();
//    WriteBaseComplex_ColorFacesVTK("RemoveOnePadding.face.vtk");
//    WriteBaseComplex("RemoveOnePadding.hex.vtk");
}

void BaseComplexEditor::CollapseComponent()
{
    for (auto & singularEdge : m_baseComplex.SingularityI.E) {
        if (singularEdge.separatedFacePatchIds.size() == 4) {
            CollapseComponent(singularEdge);
        }
    }
    WriteBaseComplex_ColorFacesVTK("CollapseComponent.vtk");
}

std::vector<size_t> BaseComplexEditor::GetTargetNeighborComponentCidesGroup(const SingularE& singularEdge, const std::vector<size_t>& separatedFacePatchIds) const {
    std::vector<size_t> targetNeighborComponentCidsGroup;
    size_t maxContainComponentFacesCount = 0;
    for (const auto& neighborComponentCidsGroup : singularEdge.neighborComponentCidsGroups) {
        int count = 0;
        for (const auto& componentCid : neighborComponentCidsGroup) {
            const auto& component = m_baseComplex.componentC.at(componentCid);
            for (const auto& componentFaceId : component.Fids) {
                for (const auto separatedFacePatchId : separatedFacePatchIds) {
                    auto& separatedComponentFaceIds = singularEdge.separatedComponentFids.at(separatedFacePatchId);
                    for (auto separatedComponentFaceId : separatedComponentFaceIds) {
                        if (separatedComponentFaceId == componentFaceId) ++count;
                    }
                }
            }
            if (count > maxContainComponentFacesCount) {
                maxContainComponentFacesCount = count;
                targetNeighborComponentCidsGroup = neighborComponentCidsGroup;
            }
        }
    }
    return targetNeighborComponentCidsGroup;
}

bool IsTwoEdgesOnFace(const Mesh& mesh, const Edge& edge1, const Edge& edge2, const Face& face) {
    int count = 0;
    for (auto edgeid : face.Eids)
        if (edgeid == edge1.id || edgeid == edge2.id) ++count;
    return count == 2;
}

void BaseComplexEditor::GetStructureIdSets(const std::vector<size_t>& targetNeighborComponentCidsGroup,
        std::unordered_set<size_t>& vids, std::unordered_set<size_t>& eids, std::unordered_set<size_t>& fids, std::unordered_set<size_t>& cids) const {
    for (auto targetNeighborComponentCid : targetNeighborComponentCidsGroup) {
        const auto& component = m_baseComplex.componentC.at(targetNeighborComponentCid);
        for (auto cellid : component.cids_patch) {
            auto& cell = mesh.C.at(cellid);
            for (auto vid : cell.Vids)
                if (vids.find(vid) == vids.end()) vids.insert(vid);
            for (auto eid : cell.Eids)
                if (eids.find(eid) == eids.end()) eids.insert(eid);
            for (auto fid : cell.Fids)
                if (fids.find(fid) == fids.end()) fids.insert(fid);
            cids.insert(cellid);
            cell.isActive = false;  // make the cell invalid;
        }
    }
}
std::vector<std::vector<std::vector<size_t>>> BaseComplexEditor::GetStructureVids(const SingularE& singularEdge, const std::vector<size_t>& separatedFacePatchIds,
        const std::vector<size_t>& targetNeighborComponentCidsGroup) const {
    std::unordered_set<size_t> vids, eids, fids, cids;
    GetStructureIdSets(targetNeighborComponentCidsGroup, vids, eids, fids, cids);
    return GetStructureVids(singularEdge, vids, eids, fids, cids);
}

size_t GetBeginEdgeIds(const SingularE& singularEdge, const Vertex& begin_v, const std::unordered_set<size_t>& eids) {
    size_t begin_eid = MAXID;
    for (auto eid : begin_v.N_Eids) {
        if (eids.find(eid) == eids.end()) continue;  // cannot find in the set;
        if (eid != singularEdge.es_link.front() && eid != singularEdge.es_link.back()) {
            begin_eid = eid;
            break;
        }
    }
    return begin_eid;
}

size_t GetNextEdgeIds(const Mesh& mesh, const size_t begin_eid, const size_t next_eid, const std::unordered_set<size_t>& fids) {
    size_t res = MAXID;
    const auto& begin_e = mesh.E.at(begin_eid);
    const auto& next_e = mesh.E.at(next_eid);
    for (auto fid : next_e.N_Fids) {
        if (fids.find(fid) == fids.end()) continue;  // cannot find in the set;
        int count = 0;
        const auto& face = mesh.F.at(fid);
        for (auto eid : face.Eids)
            if (eid == begin_eid || eid == next_eid) ++count;
        if (count == 2) {
            for (auto eid : face.Eids) {
                const auto& edge = mesh.E.at(eid);
                if (edge.Vids[0] != begin_e.Vids[0] && edge.Vids[0] != begin_e.Vids[1] &&
                    edge.Vids[1] != begin_e.Vids[0] && edge.Vids[1] != begin_e.Vids[1]) return eid;
            }
        }
    }
    return res;
}

std::vector<size_t> GetOneLineVids(const Mesh& mesh, const Vertex& begin_v, const Edge& begin_e,
        const std::unordered_set<size_t>& vids, const std::unordered_set<size_t>& eids) {
    std::unordered_map<size_t, bool> v_visited, e_visited;
    for (auto i : vids) v_visited[i] = false;
    for (auto i : eids) e_visited[i] = false;
    std::vector<size_t> res;
    res.push_back(begin_v.id);
    v_visited[begin_v.id] = true;
    size_t curr_vid = begin_e.Vids[0] == begin_v.id ? begin_e.Vids[1] : begin_e.Vids[0];
    res.push_back(curr_vid);
    v_visited[curr_vid] = true;
    const auto& curr_v = mesh.V.at(curr_vid);

    e_visited[begin_e.id] = true;
    size_t next_eid = MAXID;
    for (auto next : begin_e.consecutiveEids) {
        if (eids.find(next) == eids.end() || e_visited[next]) continue;
        next_eid = next;
        e_visited[next_eid] = true;
    }

    while (next_eid != MAXID) {
        const auto& next_e = mesh.E.at(next_eid);
        size_t next_vid = next_e.Vids[0] == curr_vid ? next_e.Vids[1] : next_e.Vids[0];
        v_visited[next_vid] = true;
        curr_vid = next_vid;
        res.push_back(curr_vid);

        next_eid = MAXID;
        for (auto next : next_e.consecutiveEids) {
            if (eids.find(next) == eids.end() || e_visited[next]) continue;
            next_eid = next;
            e_visited[next_eid] = true;
        }
    }

    return res;
}

std::vector<std::vector<size_t>> GetNextPlaneVids(const Mesh& mesh, const std::vector<std::vector<size_t>>& planeVids, std::unordered_set<size_t>& vids) {
    std::vector<std::vector<size_t>> res;
    for (auto& row : planeVids)
        for (auto vid : row)
            vids.erase(vid);
    for (auto& row : planeVids) {
        std::vector<size_t> r;
        for (auto& vid : row) {
            const auto& v = mesh.V.at(vid);
            for (auto neighbor_vid : v.N_Vids)
                if (vids.find(neighbor_vid) != vids.end()) r.push_back(neighbor_vid);
        }
        if (r.empty()) break;
        res.push_back(r);
    }
    return res;
}
std::vector<std::vector<std::vector<size_t>>> BaseComplexEditor::GetStructureVids(const SingularE& singularEdge,
        std::unordered_set<size_t>& vids, std::unordered_set<size_t>& eids, std::unordered_set<size_t>& fids, std::unordered_set<size_t>& cids) const {
    std::vector<std::vector<std::vector<size_t>>> res;
    size_t begin_vid = singularEdge.vs_link.front();
    const auto& begin_v = mesh.V.at(begin_vid);
    size_t begin_eid = GetBeginEdgeIds(singularEdge, begin_v, eids);
    const auto& begin_e = mesh.E.at(begin_eid);

    std::vector<std::vector<size_t>> planeVids;
    auto lineVids = GetOneLineVids(mesh, begin_v, begin_e, vids, eids);
    planeVids.push_back(lineVids);
    size_t curr_eid = begin_eid;
    for (size_t i = 1; i < lineVids.size(); ++i) {
        size_t next_vid = singularEdge.vs_link.at(i);
        const auto& next_v = mesh.V.at(next_vid);
        size_t next_eid = GetNextEdgeIds(mesh, curr_eid, singularEdge.es_link.at(i - 1), fids);
        curr_eid = next_eid;
        const auto& next_e = mesh.E.at(next_eid);
        lineVids = GetOneLineVids(mesh, next_v, next_e, vids, eids);
        planeVids.push_back(lineVids);
    }
    res.push_back(planeVids);

    std::unordered_set<size_t> vids_copy = vids;
    planeVids = GetNextPlaneVids(mesh, planeVids, vids_copy);
    while (!planeVids.empty()) {
        res.push_back(planeVids);
        planeVids = GetNextPlaneVids(mesh, planeVids, vids_copy);
    }

    return res;
}


void BaseComplexEditor::CollapseComponent(SingularE& singularEdge) {
    std::vector<size_t> separatedFacePatchIds;
    if (!CanNeighborCompoenntsCollapsed(singularEdge, separatedFacePatchIds)) return;

    std::vector<size_t> targetNeighborComponentCidsGroup = GetTargetNeighborComponentCidesGroup(singularEdge, separatedFacePatchIds);
    if (targetNeighborComponentCidsGroup.empty()) {
        std::cerr << "Error in EraseSeparatedFaces!\n";
        return;
    }

    auto vids = GetStructureVids(singularEdge, separatedFacePatchIds, targetNeighborComponentCidsGroup);
    auto x = vids[0][0].size();
    auto y = vids[0].size();
    auto z = vids.size();
    if (x != z)
        std::cout << "Error in CollapseComponent\n";

    std::unordered_map<size_t, size_t> key_val;
    std::vector<bool> visited(mesh.V.size());

    for (int i = 0; i < 1; ++i)
        for (int j = 0; j < vids[i].size(); ++j)
            for (int k = 0; k < vids[i][j].size(); ++k)
                key_val[vids[i][j][k]] = vids[k][j][k];
    for (int i = 0; i < vids.size(); ++i)
        for (int j = 0; j < vids[i].size(); ++j)
            for (int k = 0; k < 1; ++k)
                key_val[vids[i][j][k]] = vids[i][j][i];

    for (auto& item : key_val) {
        Vertex& key = mesh.V.at(item.first);
        Vertex& val = mesh.V.at(item.second);
        //key.id = val.id;
        key.x = val.x;
        key.y = val.y;
        key.z = val.z;
    }
    for (auto& cell : mesh.C) {
        //if (cell.isActive) continue;
        for (auto& vid : cell.Vids)
            if (key_val.find(vid) != key_val.end()) {
                vid = key_val[vid];
            }
    }
}

bool BaseComplexEditor::GetNext(size_t& vid, int keys[], int& val, std::unordered_set<size_t>& searchVids,
        std::unordered_map<size_t, size_t>& key_val, std::vector<bool>& visited) {

    return true;
}

void BaseComplexEditor::ModifyIdOfSeparatedFacePatch(SingularE& singularEdge, std::vector<size_t>& separatedFacePatchIds, std::unordered_map<size_t, size_t>& key_val) {
    std::unordered_set<size_t> vids;
    for (const auto separatedFacePatchId : separatedFacePatchIds) {
        auto& separatedComponentFaceIds = singularEdge.neighborComponentFidsGroups.at(separatedFacePatchId);
        for (auto separatedComponentFaceId : separatedComponentFaceIds) {
            const auto& componentFace = m_baseComplex.componentF.at(separatedComponentFaceId);
            for (auto faceid : componentFace.fids_patch) {
                auto& face = mesh.F.at(faceid);
                for (auto vid : face.Vids)
                    if (vids.find(vid) != vids.end()) vids.insert(vid);
            }
        }
    }
    for (auto vid : vids) {
        auto& v = mesh.V.at(vid);
        auto& tv = mesh.V.at(key_val[vid]);
        v.id = tv.id;
        v.x = tv.x;
        v.y = tv.y;
        v.z = tv.z;
    }
}
void BaseComplexEditor::MoveSingularity(SingularE& singularEdge, const size_t neigborId)
{
    for (int i = 0; i < singularEdge.separatedFacePatchIds.size(); ++i) {
        if (i == neigborId) {
            m_baseComplex.componentF.at(singularEdge.N_Fids[neigborId]).isActive = false;
            for (auto face_id : m_baseComplex.componentF.at(singularEdge.N_Fids[neigborId]).fids_patch) {
                mesh.F.at(face_id).isActive = false;
            }
        }
        else {
             for (auto & face_id : m_baseComplex.separatedFacePatches.at(singularEdge.separatedFacePatchIds.at(i))) {
                 mesh.F.at(face_id).isActive = false;
                 m_baseComplex.componentF.at(mesh.F.at(face_id).componentFid).isActive = false;
             }
        }
    }
}

std::vector<size_t> singleNumber(const std::vector<size_t>& nums) {
    int aXorb = 0;  // the result of a xor b;
    for (auto item : nums) aXorb ^= item;
    int lastBit = (aXorb & (aXorb - 1)) ^ aXorb;  // the last bit that a diffs b
    int intA = 0, intB = 0;
    for (auto item : nums) {
        // based on the last bit, group the items into groupA(include a) and groupB
        if (item & lastBit) intA = intA ^ item;
        else intB = intB ^ item;
    }
    return std::vector<size_t>{(size_t)intA, (size_t)intB};
}

bool BaseComplexEditor::CanNeighborCompoenntsCollapsed(SingularE& singularEdge, std::vector<size_t>& separatedFacePatchIds) {
    for (auto edgeid : singularEdge.es_link)
        if (mesh.E.at(edgeid).isSharpFeature) return false;

    int count = 0;
    for (int i = 0; i < singularEdge.separatedFacePatchIds.size(); ++i) {
        std::vector<size_t>&  separatedComponentEids = singularEdge.separatedComponentEids.at(i);
        bool foundOtherSingularEdge = false;
        for (auto componentEid : separatedComponentEids) {
            if (!mesh.E.at(m_baseComplex.componentE.at(componentEid).eids_link.front()).isSingularity) continue;
            bool isOnSingularEdge = false;
            for (auto edgeid : singularEdge.es_link) {
                if (componentEid == mesh.E.at(edgeid).componentEid) {
                    isOnSingularEdge = true;
                    break;
                }
            }
            if (!isOnSingularEdge){
                foundOtherSingularEdge = true;
                break;
            }
        }
        if (foundOtherSingularEdge) continue;
        ++count;
    }
    if (count != 2) return false;

    //std::vector<size_t> separatedFacePatchIds;
    for (int i = 0; i < singularEdge.separatedFacePatchIds.size(); ++i) {
        std::vector<size_t>&  separatedComponentEids = singularEdge.separatedComponentEids.at(i);
        bool foundOtherSingularEdge = false;
        for (auto componentEid : separatedComponentEids) {
            if (!mesh.E.at(m_baseComplex.componentE.at(componentEid).eids_link.front()).isSingularity) continue;
            bool isOnSingularEdge = false;
            for (auto edgeid : singularEdge.es_link) {
                if (componentEid == mesh.E.at(edgeid).componentEid) {
                    isOnSingularEdge = true;
                    break;
                }
            }
            if (!isOnSingularEdge){
                foundOtherSingularEdge = true;
                break;
            }
        }
        if (foundOtherSingularEdge) continue;
        for (auto & face_id : m_baseComplex.separatedFacePatches.at(singularEdge.separatedFacePatchIds.at(i))) {
            mesh.F.at(face_id).isActive = false;
            m_baseComplex.componentF.at(mesh.F.at(face_id).componentFid).isActive = false;
        }
        // for connect to the opposite parallel component edge;
        separatedFacePatchIds.push_back(i);
    }

    if (separatedFacePatchIds.empty()) return false;
    return true;
}

void BaseComplexEditor::EraseSeparatedFaces(SingularE& singularEdge)
{
    std::vector<size_t> separatedFacePatchIds;
    if (!CanNeighborCompoenntsCollapsed(singularEdge, separatedFacePatchIds));

    std::vector<size_t> targetNeighborComponentCidsGroup = GetTargetNeighborComponentCidesGroup(singularEdge, separatedFacePatchIds);
    if (targetNeighborComponentCidsGroup.empty()) {
        std::cerr << "Error in EraseSeparatedFaces!\n";
        return;
    }

    std::vector<size_t> singularEdgeSeparatedNeighborFaceEdgeIds;
    for (const auto separatedFacePatchId : separatedFacePatchIds) {
        auto& separatedComponentFaceIds = singularEdge.neighborComponentFidsGroups.at(separatedFacePatchId);
        for (auto separatedComponentFaceId : separatedComponentFaceIds) {
            auto& separatedComponentFace = m_baseComplex.componentF.at(separatedComponentFaceId);
            std::copy(separatedComponentFace.Eids.begin(), separatedComponentFace.Eids.end(), back_inserter(singularEdgeSeparatedNeighborFaceEdgeIds));
        }
    }
    std::sort(singularEdgeSeparatedNeighborFaceEdgeIds.begin(), singularEdgeSeparatedNeighborFaceEdgeIds.end());
    const size_t s = std::distance(singularEdgeSeparatedNeighborFaceEdgeIds.begin(),
            std::unique(singularEdgeSeparatedNeighborFaceEdgeIds.begin(), singularEdgeSeparatedNeighborFaceEdgeIds.end()));
    singularEdgeSeparatedNeighborFaceEdgeIds.resize(s);

    std::vector<size_t> singularEdgeOppositeComponentEdgeIds;
    for (const auto componentCid : targetNeighborComponentCidsGroup) {
        const auto& component = m_baseComplex.componentC.at(componentCid);
        for (auto componentEid : component.Eids) {
            const auto& componentEdge = m_baseComplex.componentE.at(componentEid);
//            const auto vid1 = componentEdge.Vids[0];
//            const auto vid2 = componentEdge.Vids[1];
            const auto vid1 = mesh.V.at(componentEdge.vids_link.front()).component_id;
            const auto vid2 = mesh.V.at(componentEdge.vids_link.back()).component_id;
            bool foundIntersection = false;
            for (auto separatedComponentEid : singularEdgeSeparatedNeighborFaceEdgeIds) {
                const auto& separatedComponentEdge = m_baseComplex.componentE.at(separatedComponentEid);
//                const auto vid_1 = separatedComponentEdge.Vids[0];
//                const auto vid_2 = separatedComponentEdge.Vids[1];
                const auto vid_1 = mesh.V.at(separatedComponentEdge.vids_link.front()).component_id;
                const auto vid_2 = mesh.V.at(separatedComponentEdge.vids_link.back()).component_id;
                if (vid1 != vid_1 && vid1 != vid_2 && vid2 != vid_1 && vid2 != vid_2) {
                    ;
                } else {
                    foundIntersection = true;
                    break;
                }
            }
            if (!foundIntersection) singularEdgeOppositeComponentEdgeIds.push_back(componentEid);
        }
    }
    std::vector<size_t> singularEdgeOppositeComponentEdgeVids;
    for (auto componentEid : singularEdgeOppositeComponentEdgeIds) {
        const auto& componentEdge = m_baseComplex.componentE.at(componentEid);
        const auto vid1 = mesh.V.at(componentEdge.vids_link.front()).component_id;
        const auto vid2 = mesh.V.at(componentEdge.vids_link.back()).component_id;
        singularEdgeOppositeComponentEdgeVids.push_back(vid1);
        singularEdgeOppositeComponentEdgeVids.push_back(vid2);
    }
    std::vector<size_t> oppositeComponentEdgeFrontEndVids = singleNumber(singularEdgeOppositeComponentEdgeVids);
    if (oppositeComponentEdgeFrontEndVids.size() != 2) std::cerr << "Error in EraseSeparatedFaces singleNumber!\n";

    size_t oppositeEdgeFrontVid = oppositeComponentEdgeFrontEndVids.front();
    size_t oppositeEdgeBackVid = oppositeComponentEdgeFrontEndVids.back();
    size_t singularEdgeFrontVid = mesh.V.at(singularEdge.vs_link.front()).component_id;
    size_t singularEdgeBackVid = mesh.V.at(singularEdge.vs_link.back()).component_id;
    if (!isAligned(singularEdge, targetNeighborComponentCidsGroup, oppositeEdgeFrontVid, oppositeEdgeBackVid))
        std::swap(oppositeEdgeFrontVid, oppositeEdgeBackVid);
    Face f(4);
    f.id = mesh.F.size();
    f.Vids = std::vector<size_t> {oppositeEdgeFrontVid, oppositeEdgeBackVid, singularEdgeBackVid, singularEdgeFrontVid};
    mesh.F.push_back(f);
    ComponentFace componentFace;
    componentFace.fids_patch = std::vector<size_t> {f.id};
    componentFace.id = m_baseComplex.componentF.size();
    m_baseComplex.componentF.push_back(componentFace);
}

bool BaseComplexEditor::isAligned(SingularE& singularEdge, const std::vector<size_t>& targetNeighborComponentCidsGroup, size_t oppositeEdgeFrontVid, size_t oppositeEdgeBackVid) {
    size_t singularEdgeFrontVid = mesh.V.at(singularEdge.vs_link.front()).component_id;
    size_t singularEdgeBackVid = mesh.V.at(singularEdge.vs_link.back()).component_id;
    for (const auto componentCid : targetNeighborComponentCidsGroup) {
        const auto& component = m_baseComplex.componentC.at(componentCid);
        for (auto componentFid : component.Fids) {
            const auto& componentFace = m_baseComplex.componentF.at(componentFid);
            int count = 0;
            for (const auto vid : componentFace.Vids) {
                if (vid == oppositeEdgeFrontVid || vid == singularEdgeFrontVid) ++count;
            }
            if (count == 2)
                return true;
        }
    }
    return false;
}
size_t BaseComplexEditor::GetOppositeComponentEdgeId(const ComponentCell& component, size_t componentFaceId1, size_t componentFaceId2) {
    std::vector<size_t> existingComponentEdgeIds = m_baseComplex.componentF.at(componentFaceId1).Eids;
    std::copy(m_baseComplex.componentF.at(componentFaceId1).Eids.begin(), m_baseComplex.componentF.at(componentFaceId1).Eids.end(), back_inserter(existingComponentEdgeIds));
    for (auto componentEdgeId : component.Eids) {
        bool foundConnectionInExistingComponentEdgeIds = false;
        for (auto existingComponentEdgeId : existingComponentEdgeIds) {
            if (existingComponentEdgeId == componentEdgeId) {
                foundConnectionInExistingComponentEdgeIds = true;
                break;
            }
        }
        if (foundConnectionInExistingComponentEdgeIds) continue;

        const auto frontVid = m_baseComplex.componentE.at(componentEdgeId).vids_link.front();
        const auto backVid = m_baseComplex.componentE.at(componentEdgeId).vids_link.back();
        for (auto existingComponentEdgeId : existingComponentEdgeIds) {
            const auto existingFrontVid = m_baseComplex.componentE.at(existingComponentEdgeId).vids_link.front();
            const auto existingBackVid = m_baseComplex.componentE.at(existingComponentEdgeId).vids_link.back();
            if (frontVid == existingFrontVid || backVid == existingFrontVid || frontVid == existingBackVid || backVid == existingBackVid) {
                foundConnectionInExistingComponentEdgeIds = true;
                break;
            }
        }
        if (foundConnectionInExistingComponentEdgeIds) continue;
        return componentEdgeId;
    }
    std::cerr << "Error in GetOppositeComponentEdgeIds!\n";
    return 0;
}

bool BaseComplexEditor::IsOnBoundary(const SingularE& singularEdge) const
{
    return m_baseComplex.mesh.E.at(singularEdge.es_link.front()).isBoundary;
}

bool BaseComplexEditor::IsConnectedToBoundary(const SingularE& singularEdge) const
{
    return m_baseComplex.mesh.V.at(singularEdge.vs_link.front()).isBoundary || m_baseComplex.mesh.V.at(singularEdge.vs_link.back()).isBoundary;
}

void BaseComplexEditor::Remove(SingularE& singularEdge)
{
    for (int i = 0; i < singularEdge.separatedFacePatchIds.size(); ++i) {
        for (auto & face_id : m_baseComplex.separatedFacePatches.at(singularEdge.separatedFacePatchIds.at(i))) {
            mesh.F.at(face_id).isActive = false;
            m_baseComplex.componentF.at(mesh.F.at(face_id).componentFid).isActive = false;
        }
    }
}
void BaseComplexEditor::RemoveOnePadding()
{
    for (auto& singularEdge : m_baseComplex.SingularityI.E) {
        if (!IsOnBoundary(singularEdge) && IsConnectedToBoundary(singularEdge)) {
            Remove(singularEdge);
        }
    }
}

void BaseComplexEditor::WriteBaseComplex_ColorFacesVTK(const char *filename) const
{
    const std::vector<Vertex>& V = mesh.V;
    const std::vector<Edge>& E = mesh.E;
    const std::vector<Face>& F = m_baseComplex.F;
    const std::vector<size_t>& Vids = m_baseComplex.Vids;
    const std::vector<size_t> Eids = m_baseComplex.Eids;
    const std::vector<size_t> Fids = m_baseComplex.Fids;
    const std::vector<size_t> Cids = m_baseComplex.Cids;
    const std::vector<ComponentVertex>& componentV = m_baseComplex.componentV;
    const std::vector<ComponentEdge>& componentE = m_baseComplex.componentE;
    const std::vector<ComponentFace>& componentF = m_baseComplex.componentF;
    const std::vector<ComponentCell>& componentC = m_baseComplex.componentC;

    std::ofstream ofs(filename);
    ofs << "# vtk DataFile Version 2.0" << std::endl
        << filename << std::endl
        << "ASCII" << std::endl << std::endl
        << "DATASET POLYDATA" << std::endl;
    ofs << "POINTS " << V.size() << " double" << std::endl;

    for (size_t i = 0; i < V.size(); i++)
        ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << std::endl;

//    ofs << "VERTICES " << Vids.size() << " " << 2 * Vids.size() << std::endl;
//    for (size_t i = 0; i < Vids.size(); i++)
//        ofs << "1 " << V.at(Vids.at(i)).id << std::endl;
//
//    size_t eid_count = 0;
//    for (auto& componentEdge : componentE)
//        eid_count += 1 + componentEdge.vids_link.size();
//
//    ofs << "LINES " << componentE.size() << " " << eid_count << std::endl;
//    for (auto& componentEdge : componentE) {
//        ofs << componentEdge.vids_link.size();
//        for (auto vid : componentEdge.vids_link)
//            ofs << " " << vid;
//        ofs << "\n";
//    }

    size_t face_count = 0;
    for (const auto& facePatch : componentF) {
        for (const auto& faceid : facePatch.fids_patch) {
            const Face& face = mesh.F.at(faceid);
            if (face.isActive) ++face_count;
        }
    }
    ofs << "POLYGONS " << face_count << " " << 5 * face_count << std::endl;
    for (const auto& facePatch : componentF) {
        for (const auto& faceid : facePatch.fids_patch) {
            const Face& face = mesh.F.at(faceid);
            if (!face.isActive) continue;
            ofs << face.Vids.size();
            for (size_t j = 0; j < face.Vids.size(); j++)
                ofs << " " << face.Vids.at(j);
            ofs << "\n";
        }
    }

//    ofs << "CELL_DATA " << face_count << std::endl
//        << "SCALARS " << " faceid" << " int 1" << std::endl
//        << "LOOKUP_TABLE default" << std::endl;
//    for (size_t i = 0; i < Vids.size(); i++)
//        ofs << (V.at(Vids.at(i)).isSingularity ? componentF.size() : 0) << std::endl;
//    for (auto& componentEdge : componentE)
//        ofs << (E.at(componentEdge.eids_link.front()).isSingularity ? componentF.size() : 0) << std::endl;
//    for (const auto& facePatch : componentF)
//        for (const auto& faceid : facePatch.fids_patch)
//            ofs << facePatch.id << std::endl;


//    ofs << "SCALARS " << " color" << " int 1" << std::endl
//        << "LOOKUP_TABLE default" << std::endl;
//    for (size_t i = 0; i < Vids.size(); i++)
//        ofs << (V.at(Vids.at(i)).isSingularity ? 8 : 0) << std::endl;
//    for (auto& componentEdge : componentE)
//        ofs << (E.at(componentEdge.eids_link.front()).isSingularity ? 8 : 0) << std::endl;
//    for (const auto& facePatch : componentF)
//        for (const auto& faceid : facePatch.fids_patch)
//            ofs << facePatch.id % 8 << std::endl;

}

void BaseComplexEditor::WriteMeshVTK(const char *filename) const
{
    const std::vector<Vertex>& V = mesh.V;
    const std::vector<Edge>& E = mesh.E;
    const std::vector<Face>& F = mesh.F;
    const std::vector<Cell>& C = mesh.C;

    std::ofstream ofs(filename);
    ofs << "# vtk DataFile Version 2.0" << std::endl
        << filename << std::endl
        << "ASCII" << std::endl << std::endl
        << "DATASET UNSTRUCTURED_GRID" << std::endl;
    ofs << "POINTS " << V.size() << " float" << std::endl;

    for (size_t i = 0; i < V.size(); i++)
        ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << std::endl;

    size_t cell_count = 0;
    for (const auto& cell : C)
        if (cell.isActive) ++cell_count;
    ofs << "CELLS " << cell_count << " " << 9 * cell_count << std::endl;
    for (const auto& cell : C) {
        if (!cell.isActive) continue;
        ofs << cell.Vids.size();
        for (size_t j = 0; j < cell.Vids.size(); j++)
            ofs << " " << cell.Vids.at(j);
        ofs << "\n";
    }
    ofs << "CELL_TYPES " << cell_count << "\n";
    for (const auto& cell : C) {
        if (!cell.isActive) continue;
        ofs << "12\n";
    }
}

void BaseComplexEditor::WriteBaseComplex(const char *filename) const
{
    const std::vector<Vertex>& V = mesh.V;
    const std::vector<Edge>& E = mesh.E;
    const std::vector<Face>& F = m_baseComplex.F;
    const std::vector<size_t>& Vids = m_baseComplex.Vids;
    const std::vector<size_t> Eids = m_baseComplex.Eids;
    const std::vector<size_t> Fids = m_baseComplex.Fids;
    const std::vector<size_t> Cids = m_baseComplex.Cids;
    const std::vector<ComponentVertex>& componentV = m_baseComplex.componentV;
    const std::vector<ComponentEdge>& componentE = m_baseComplex.componentE;
    const std::vector<ComponentFace>& componentF = m_baseComplex.componentF;
    const std::vector<ComponentCell>& componentC = m_baseComplex.componentC;

    std::ofstream ofs(filename);
    ofs << "# vtk DataFile Version 2.0" << std::endl
        << filename << std::endl
        << "ASCII" << std::endl << std::endl
        << "DATASET UNSTRUCTURED_GRID" << std::endl;
    ofs << "POINTS " << V.size() << " float" << std::endl;

    for (size_t i = 0; i < V.size(); i++)
        ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << std::endl;
    size_t cell_count = 0;
    //std::cout << "CELLS " << mesh.C.size() << std::endl;
    for (const auto& cell : mesh.C) {
        bool containsBoundary = false;
        for (const auto& faceid : cell.Fids)
            if (mesh.F.at(faceid).isBoundary) {
                containsBoundary = true;
                break;
            }
        if (containsBoundary) continue;
//        bool isAllFacesValid = true;
//        for (const auto& faceid : cell.Fids)
//            if (!mesh.F.at(faceid).isActive) {
//                isAllFacesValid = false;
//                break;
//            }
//        if (isAllFacesValid) ++cell_count;
        ++cell_count;
    }
    std::cout << "cell_count " << cell_count << std::endl;
    ofs << "CELLS " << cell_count << " " << 9 * cell_count << std::endl;
    for (const auto& cell : mesh.C) {
        bool containsBoundary = false;
        for (const auto& faceid : cell.Fids)
            if (mesh.F.at(faceid).isBoundary) {
                containsBoundary = true;
                break;
            }
        if (containsBoundary) continue;
//        bool isAllFacesValid = true;
//        for (const auto& faceid : cell.Fids)
//            if (!mesh.F.at(faceid).isActive) {
//                isAllFacesValid = false;
//                break;
//            }
//        if (!isAllFacesValid) continue;
        ofs << cell.Vids.size();
        for (const auto& vid : cell.Vids)
            ofs << " " << vid;
        ofs << "\n";
    }

    ofs << "CELL_TYPES " << cell_count << "\n";
    for (int i = 0; i < cell_count; ++i)
        ofs << 12 << "\n";
}
