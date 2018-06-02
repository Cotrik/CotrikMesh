/*
 * BaseComplex.cpp
 *
 *  Created on: May 31, 2017
 *      Author: cotrik
 */

#include "BaseComplexQuad.h"
#include "MeshFileWriter.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <stack>

BaseComplexQuad::BaseComplexQuad(Mesh& mesh)
: BaseComplex(mesh)
//: mesh(mesh)
{
    // TODO Auto-generated constructor stub

}

BaseComplexQuad::BaseComplexQuad(const BaseComplexQuad& baseComplex)
: BaseComplex(baseComplex)
//: mesh((Mesh&)baseComplex.GetMesh())
//, components(baseComplex.components)
{

}

BaseComplexQuad::~BaseComplexQuad()
{
    // TODO Auto-generated destructor stub
}

const Mesh& BaseComplexQuad::GetMesh() const
{
    return mesh;
}

size_t BaseComplexQuad::GetNextSingularEdgeId(const size_t vid, const size_t cur_edge_id)
{
    size_t nextEdgeId = MAXID;
    const Edge& currentEdge = mesh.E.at(cur_edge_id);
    int num1 = 0, num2 = 0;
    for (size_t i = 0; i < mesh.V[vid].N_Eids.size(); i++) {
        const size_t next_edge_id = mesh.V[vid].N_Eids[i];
        if (next_edge_id == cur_edge_id) continue;
        const Edge& nextEdge = mesh.E.at(next_edge_id);
        if (nextEdge.N_Cids.size() == currentEdge.N_Cids.size()) {
            num1++;
            nextEdgeId = next_edge_id;
        }
        else if (nextEdge.isSingularity) num2++;
    }
    if (num1 == 1 && num2 == 0)
        return nextEdgeId;

    return MAXID;
}

bool BaseComplexQuad::TraceEdge(const size_t vid, const size_t eid, std::vector<size_t>& Vids, std::vector<size_t>& Eids, std::vector<bool>& is_edge_visited)
{
    size_t current_eid = eid;
    if (!is_edge_visited[current_eid]) {
        Eids.push_back(current_eid);
        is_edge_visited[current_eid] = true;
    }

    size_t current_vid = vid;
    Vids.push_back(current_vid);
    size_t next_eid = GetNextSingularEdgeId(current_vid, current_eid);
    while (next_eid != MAXID) {
        Eids.push_back(next_eid);
        current_eid = next_eid;
        is_edge_visited[current_eid] = true;
        if (mesh.E[current_eid].Vids[0] == current_vid) current_vid = mesh.E[current_eid].Vids[1];
        else
            current_vid = mesh.E[current_eid].Vids[0];

        Vids.push_back(current_vid);
        if (current_eid == eid) {
            //is_circle = true;
            //break;
            return true;
        }
        next_eid = GetNextSingularEdgeId(current_vid, current_eid);
    }

    return false;
}
void BaseComplexQuad::AddSingularV(const size_t current_vid, const size_t singularVid, size_t& ending_singular_vid)
{
    bool existed = false;
    for (size_t j = 0; j < SingularityI.V.size(); j++) {
        if (SingularityI.V[j].id_mesh == current_vid) {
            existed = true;
            ending_singular_vid = j;
        }
    }
    if (!existed) {
        SingularV sv;
        sv.id = singularVid;
        sv.id_mesh = current_vid;
        sv.isBoundary = mesh.V[current_vid].isBoundary;
        SingularityI.V.push_back(sv);
        ending_singular_vid = sv.id;
    }
}

void BaseComplexQuad::AddSingularE(const std::vector<size_t>& leftVids, const std::vector<size_t>& leftEids,
        const std::vector<size_t>& rightVids, const std::vector<size_t>& rightEids,
        const size_t left_singular_vid, const size_t right_singular_vid, const size_t singularEid)
{
    SingularE se;
    se.startend_Vid[0] = left_singular_vid;
    se.startend_Vid[1] = right_singular_vid;

    for (int j = leftEids.size() - 1; j >= 0; j--) {
        se.es_link.push_back(leftEids[j]);
        mesh.E.at(leftEids[j]).singularEid = singularEid;
    }
    for (int j = leftVids.size() - 1; j >= 0; j--)
        se.vs_link.push_back(leftVids[j]);
    for (size_t j = 0; j < rightEids.size(); j++) {
        se.es_link.push_back(rightEids[j]);
        mesh.E.at(rightEids[j]).singularEid = singularEid;
    }
    for (int j = 0; j < rightVids.size(); j++)
        se.vs_link.push_back(rightVids[j]);
    se.isBoundary = mesh.E[leftEids[0]].isBoundary;
    se.id = singularEid;
    se.edge_type = EDGE_TYPE_REGULAR;
    SingularityI.E.push_back(se);
}


void BaseComplexQuad::AddcircularSingularE(const std::vector<size_t>& leftVids, const std::vector<size_t>& leftEids, const size_t singularEid)
{
    SingularE se;
    se.startend_Vid[0] = leftVids.front();
    se.startend_Vid[1] = leftVids.back();
    for (size_t j = 0; j < leftEids.size() - 1; j++) {
        se.es_link.push_back(leftEids[j]);
        mesh.E.at(leftEids[j]).singularEid = singularEid;
    }
    for (size_t j = 0; j < leftVids.size(); j++)
        se.vs_link.push_back(leftVids[j]);
    se.isBoundary = mesh.E[leftEids[0]].isBoundary;
    se.id = singularEid;
    se.edge_type = EDGE_TYPE_CIRCULAR;
    SingularityI.E.push_back(se);
}
void BaseComplexQuad::Build()
{
    ExtractSingularVandE();
    BuildE();
    std::cout << "Finished extracting BaseComplex Edges --> " << Eids.size() << " edges\n";
    BuildV();
    std::cout << "Finished extracting BaseComplex Vertices --> " << Vids.size() << " vertices\n";
    BuildF();
    std::cout << "Finished extracting Component Cells --> " << componentC.size() << " component cells\n";
    std::cout << "Finished extracting Component Faces --> " << componentF.size() << " component faces\n";
    std::cout << "Finished extracting Component Edges --> " << componentE.size() << " component edges\n";
    std::cout << "Finished extracting Component Vertices --> " << Vids.size() << " component vertices\n";

    BuildComponentConnectivities();
    BuildSigularEdge_separatedPatches();
    AlignComponentE();
}

/*
                  3__________________2
                  /|                 /|
                 / |                / |
                /  |               /  |
            0  /___|_____________1/   |
               |   |              |   |
               |   |              |   |
               |   |              |   |
               |   |______________|___|
               |   / 7            |  /6
               |  /               | /
               | /                |/
               |/_________________/
             4                    5
*/
void BaseComplexQuad::AlignComponentE() {
    const size_t parallelEdges[3][4][2] = {
            {{0, 4}, {3, 7}, {2, 6}, {1, 5}},
            {{0, 1}, {3, 2}, {7, 6}, {4, 5}},
            {{0, 3}, {1, 2}, {5, 6}, {4, 7}}
    };
    for (auto& componentCell : componentC) {
        if (componentCell.Vids.empty()) continue; // circular component cell
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 4; ++j){
                auto vid1 = componentCell.Vids.at(parallelEdges[i][j][0]);
                auto vid2 = componentCell.Vids.at(parallelEdges[i][j][1]);
                auto v1 = Vids.at(vid1);
                auto v2 = Vids.at(vid2);
                for (auto componentEdgeId : componentCell.Eids) {
                    auto& componentEdge = componentE.at(componentEdgeId);
                    auto v_1 = componentEdge.eids_link.front();
                    auto v_2 = componentEdge.eids_link.back();
                    if (v1 == v_2 && v2 == v_1) {
                        std::reverse(componentEdge.eids_link.begin(), componentEdge.eids_link.end());
                        std::cout << "Reversed componentEdge " << componentEdgeId << "\n";
                        break;
                    }
                }
            }
    }
}

void BaseComplexQuad::BuildSigularEdge_separatedPatches()
{
    BuildSigularEdge_separatedFacePatches();
//    BuildSigularEdge_separatedEdgePatches();
    BuildSigularEdge_separatedComponentFacePatches();
    BuildSigularEdge_separatedComponentEdgePatches();
}

void BaseComplexQuad::BuildSigularEdge_separatedFacePatches()
{
    for (auto& singularityEdge : SingularityI.E) {
        const size_t edge_id = singularityEdge.es_link[0];
        const Edge& edge = mesh.E.at(edge_id);
        for (auto neighbor_face_id : edge.N_Fids) {
            for (size_t i = 0; i < separatedFacePatches.size(); ++i) {
                auto& patch = separatedFacePatches.at(i);
                bool found = false;
                for (auto face_id : patch) {
                    if (face_id == neighbor_face_id) {
                        singularityEdge.separatedFacePatchIds.push_back(i);
                        singularityEdge.separatedFacePatches.push_back(patch);
                        found = true;
                        break;
                    }
                }
                if (found) break;
            }
        }
    }
}

void BaseComplexQuad::BuildSigularEdge_separatedComponentFacePatches()
{
    for (auto& patch : separatedFacePatches) {
        std::vector<size_t> separatedComponentFaceIds;
        for (auto face_id : patch) {
            const Face& face = mesh.F.at(face_id);
            separatedComponentFaceIds.push_back(face.componentFid);
        }
        std::sort(separatedComponentFaceIds.begin(), separatedComponentFaceIds.end());
        separatedComponentFaceIds.resize(std::distance(separatedComponentFaceIds.begin(),
                std::unique(separatedComponentFaceIds.begin(), separatedComponentFaceIds.end())));
        separatedComponentFacePatches.push_back(separatedComponentFaceIds);
    }

    for (auto& singularityEdge : SingularityI.E) {
        singularityEdge.separatedComponentFids.resize(singularityEdge.separatedFacePatchIds.size());
        for (size_t i = 0; i < singularityEdge.separatedFacePatchIds.size(); ++i) {
            singularityEdge.separatedComponentFids[i] = separatedComponentFacePatches[singularityEdge.separatedFacePatchIds[i]];
        }
    }
}

void BaseComplexQuad::BuildSigularEdge_separatedComponentEdgePatches()
{
    for (auto& patch : separatedComponentFacePatches) {
        std::vector<size_t> separatedComponentEdgeIds;
        for (auto component_face_id : patch) {
            const ComponentFace& component_face = componentF.at(component_face_id);
            std::copy(component_face.Eids.begin(), component_face.Eids.end(), back_inserter(separatedComponentEdgeIds));
        }
        std::sort(separatedComponentEdgeIds.begin(), separatedComponentEdgeIds.end());
        separatedComponentEdgeIds.resize(std::distance(separatedComponentEdgeIds.begin(),
                std::unique(separatedComponentEdgeIds.begin(), separatedComponentEdgeIds.end())));
        separatedComponentEdgePatches.push_back(separatedComponentEdgeIds);
    }

    for (auto& singularityEdge : SingularityI.E) {
        singularityEdge.separatedComponentEids.resize(singularityEdge.separatedFacePatchIds.size());
        for (size_t i = 0; i < singularityEdge.separatedFacePatchIds.size(); ++i) {
            singularityEdge.separatedComponentEids[i] = separatedComponentEdgePatches[singularityEdge.separatedFacePatchIds[i]];
        }
    }
}

void BaseComplexQuad::BuildComponentConnectivities()
{
    for (auto& component : componentF) {
        component.Eids = GetComponentEids(component);
    }
    for (auto& component : componentE) {
        component.Vids = (component.vids_link.front() == component.vids_link.back()) ? std::vector<size_t>() :
                std::vector<size_t>{mesh.V.at(component.vids_link.front()).component_id, mesh.V.at(component.vids_link.back()).component_id};
        component.N_Fids = GetNeighborComponentFids(component);
    }
    for (auto& component : componentF)
        component.Vids = GetComponentVids(component);
    for (auto& component : componentC)
        component.Vids = GetComponentVids(component);

//    for (auto& singularityEdge : SingularityI.E) {
//        singularityEdge.N_Cids = GetNeighborComponentCids(singularityEdge);
//        singularityEdge.N_Fids = GetNeighborComponentFids(singularityEdge);
//        singularityEdge.neighborComponentCidsGroups = GetNeighborComponentCidsGroups(singularityEdge);
//        singularityEdge.neighborComponentFidsGroups = GetNeighborComponentFidsGroups(singularityEdge);
//    }
}

static std::vector<size_t> get_boundary_ids(const std::vector<size_t>& ids)
{
    std::vector<size_t> boundary_ids;
    boundary_ids.reserve(ids.size());
    for (size_t i = 0; i < ids.size() - 1; ++i)
        if (ids[i + 1] != ids[i])
            boundary_ids.push_back(ids[i]);
    if (ids[ids.size() - 2] != ids.back())
        boundary_ids.push_back(ids.back());
    return boundary_ids;
}

std::vector<size_t> BaseComplexQuad::GetComponentFids(const ComponentCell& component) const
{
    std::vector<size_t> fids;
    for (auto& cellid : component.cids_patch) {
        const std::vector<size_t>& cell_fids = mesh.C.at(cellid).Fids;
        std::copy(cell_fids.begin(), cell_fids.end(), back_inserter(fids));
    }
    std::sort(fids.begin(), fids.end());
    std::vector<size_t> boundary_fids = get_boundary_ids(fids);

    std::vector<size_t> boundary_component_fids(boundary_fids.size());
    for (size_t i = 0; i < boundary_fids.size(); ++i)
        boundary_component_fids[i] = mesh.F.at(boundary_fids.at(i)).componentFid;
    std::sort(boundary_component_fids.begin(), boundary_component_fids.end());
    boundary_component_fids.resize(std::distance(boundary_component_fids.begin(),
            std::unique(boundary_component_fids.begin(), boundary_component_fids.end())));
    if (boundary_component_fids.back() == MAXID)
        boundary_component_fids.pop_back();
    return boundary_component_fids;
}

std::vector<size_t> BaseComplexQuad::GetComponentEids(const ComponentCell& component) const
{
    std::vector<size_t> eids;
   for (const auto componentFaceId : component.Fids) {
       std::vector<size_t> ids = GetComponentEids(componentF.at(componentFaceId));
       std::copy(ids.begin(), ids.end(), back_inserter(eids));
   }
   std::sort(eids.begin(), eids.end());
   eids.resize(std::distance(eids.begin(), std::unique(eids.begin(), eids.end())));
   return eids;
}

std::vector<size_t> BaseComplexQuad::GetComponentVids(const ComponentCell& component) const
{
    const ComponentFace& componentFace = componentF.at(component.Fids.front());
    std::vector<size_t> quadVids = componentFace.Vids;  // top face vids;
    std::vector<size_t> vids; // bottom face vids
    for (const auto componentFaceId : component.Fids) {
       std::vector<size_t> faceVids = componentF.at(componentFaceId).Vids;
       std::copy(quadVids.begin(), quadVids.end(), back_inserter(faceVids));
       std::sort(faceVids.begin(), faceVids.end());
       if (std::distance(faceVids.begin(), std::unique(faceVids.begin(), faceVids.end())) == 8) {
           vids = componentF.at(componentFaceId).Vids;
           break;
       }
    }
    if (vids.empty()) {
        std::cerr << "Circular Cell, No Vertices\n";
        return {};
    }
    // if (vids.empty()) std::cerr << "Error in std::vector<size_t> BaseComplexQuad::GetComponentVids(const ComponentCell& component) const\n";
    std::vector<std::pair<size_t, size_t>> vid_pairs;
    for (const auto componentEdgeId : component.Eids) {
        const std::vector<size_t>& ids = componentE.at(componentEdgeId).vids_link;
        int id1 = mesh.V.at(ids.front()).component_id;
        int id2 = mesh.V.at(ids.back()).component_id;
        int count = 0;
        for (auto id : quadVids)
            if (id == id1 || id == id2) ++count;
        if (count == 1) {
            std::pair<size_t, size_t> p;
            for (auto id : quadVids) {
                if (id == id1) {
                    p.first = id1;
                    break;
                }
                if (id == id2) {
                    p.first = id2;
                    break;
                }
            }
            for (auto id : vids) {
                if (id == id1) {
                    p.second = id1;
                    break;
                }
                if (id == id2) {
                    p.second = id2;
                    break;
                }
            }
            vid_pairs.push_back(p);
        }
    }
    std::vector<std::pair<size_t, size_t>> vid_ps;
    for (auto id : quadVids) {
        std::pair<size_t, size_t> p;
        p.first = id;
        for (auto& pp : vid_pairs) {
            if (pp.first == id) {
                p.second = pp.second;
                break;
            }
        }
        vid_ps.push_back(p);
    }

    std::vector<size_t> reorient_vids;
    for (auto& pp : vid_ps)
        reorient_vids.push_back(pp.first);
    for (auto& pp : vid_ps)
        reorient_vids.push_back(pp.second);
    return reorient_vids;
}

std::vector<size_t> BaseComplexQuad::GetComponentVids(const ComponentFace& componentFace) const
{
    std::vector<size_t> vids;
   for (const auto componentEdgeId : componentFace.Eids) {
       const std::vector<size_t>& ids = componentE.at(componentEdgeId).vids_link;
       if (ids.front() != ids.back()) {
           vids.push_back(mesh.V.at(ids.front()).component_id);
           vids.push_back(mesh.V.at(ids.back()).component_id);
       }
   }
   std::sort(vids.begin(), vids.end());
   vids.resize(std::distance(vids.begin(), std::unique(vids.begin(), vids.end())));
   if (vids.size() != 4) return vids;
   // reorder to make a quad orientation.
   size_t v1 = mesh.V.at(componentE.at(componentFace.Eids.front()).vids_link.front()).component_id;
   size_t v2 = mesh.V.at(componentE.at(componentFace.Eids.front()).vids_link.back()).component_id;
   size_t v3 = MAXID;
   size_t v4 = MAXID;
   for (int i = 1; i < componentFace.Eids.size(); ++i) {
       const auto& componentEdgeId = componentFace.Eids.at(i);
       const std::vector<size_t>& ids = componentE.at(componentEdgeId).vids_link;
       int id1 = mesh.V.at(ids.front()).component_id;
       int id2 = mesh.V.at(ids.back()).component_id;
       if (id1 == v2 || id2 == v2) {
           v3 = id1 == v2 ? id2 : id1;
           break;
       }
   }
   for (int i = 1; i < componentFace.Eids.size(); ++i) {
       const auto& componentEdgeId = componentFace.Eids.at(i);
       const std::vector<size_t>& ids = componentE.at(componentEdgeId).vids_link;
       int id1 = mesh.V.at(ids.front()).component_id;
       int id2 = mesh.V.at(ids.back()).component_id;
       if ((id1 == v3 && id2 != v2) || (id2 == v3 && id1 != v2)) {
           v4 = id1 == v3 ? id2 : id1;
           break;
       }
   }
   return {v1, v2, v3, v4};
}

std::vector<size_t> BaseComplexQuad::GetComponentEids(const ComponentFace& componentFace) const
{
    std::vector<size_t> eids;
    for (auto& faceid : componentFace.fids_patch) {
        const std::vector<size_t>& face_eids = mesh.F.at(faceid).Eids;
        std::copy(face_eids.begin(), face_eids.end(), back_inserter(eids));
    }
    std::sort(eids.begin(), eids.end());
    std::vector<size_t> boundary_eids = get_boundary_ids(eids);

    std::vector<size_t> boundary_component_eids(boundary_eids.size());
    for (size_t i = 0; i < boundary_eids.size(); ++i)
        boundary_component_eids[i] = mesh.E.at(boundary_eids.at(i)).componentEid;
    std::sort(boundary_component_eids.begin(), boundary_component_eids.end());
    boundary_component_eids.resize(std::distance(boundary_component_eids.begin(),
            std::unique(boundary_component_eids.begin(), boundary_component_eids.end())));
    if (boundary_component_eids.back() == MAXID) boundary_component_eids.pop_back();
    return boundary_component_eids;
}

std::vector<size_t> BaseComplexQuad::GetNeighborComponentCids(const ComponentFace& componentFace) const
{
    const std::vector<size_t>& cids = mesh.F.at(componentFace.fids_patch.front()).N_Cids;
    std::vector<size_t> componentCids(cids.size());
    for (size_t i = 0; i < cids.size(); ++i)
        componentCids[i] = mesh.C.at(cids[i]).componentCid;
    return componentCids;
}


std::vector<size_t> BaseComplexQuad::GetNeighborComponentCids(const ComponentEdge& componentEdge) const
{
    const std::vector<size_t>& cids = mesh.E.at(componentEdge.eids_link.front()).N_Cids;
    std::vector<size_t> componentCids(cids.size());
    for (size_t i = 0; i < cids.size(); ++i)
        componentCids[i] = mesh.C.at(cids[i]).componentCid;
    return componentCids;
}

std::vector<size_t> BaseComplexQuad::GetNeighborComponentFids(const ComponentEdge& componentEdge) const
{
    const std::vector<size_t>& fids = mesh.E.at(componentEdge.eids_link.front()).N_Fids;
    std::vector<size_t> componentFids(fids.size());
    for (size_t i = 0; i < fids.size(); ++i)
        componentFids[i] = mesh.F.at(fids[i]).componentFid;
    return componentFids;
}

const std::vector<size_t>& BaseComplexQuad::GetComponentEidsLink(SingularE& singularEdge) const
{
    if (!singularEdge.componentEids_link.empty())
        return singularEdge.componentEids_link;

    std::vector<size_t> component_eids_link;
    component_eids_link.push_back(mesh.E.at(singularEdge.es_link.front()).componentEid);
    for (auto edgeid : singularEdge.es_link) {
        const size_t component_eid = mesh.E.at(edgeid).componentEid;
        if (component_eids_link.back() != component_eid)
            component_eids_link.push_back(component_eid);
    }
    singularEdge.componentEids_link = component_eids_link;

    return singularEdge.componentEids_link;
}

static bool TwoComponentCellsHaveCommonFace(const ComponentCell& componentCell1, const ComponentCell componentCell2)
{
    std::vector<size_t> component_face_ids = componentCell1.Fids;
    std::copy(componentCell2.Fids.begin(), componentCell2.Fids.end(), back_inserter(component_face_ids));
    std::sort(component_face_ids.begin(), component_face_ids.end());
    return std::unique(component_face_ids.begin(), component_face_ids.end()) != component_face_ids.end();
}

static bool TwoComponentFacesHaveCommonEdge(const ComponentFace& componentFace1, const ComponentFace componentFace2)
{
    std::vector<size_t> component_edge_ids = componentFace1.Eids;
    std::copy(componentFace2.Eids.begin(), componentFace2.Eids.end(), back_inserter(component_edge_ids));
    std::sort(component_edge_ids.begin(), component_edge_ids.end());
    return std::unique(component_edge_ids.begin(), component_edge_ids.end()) != component_edge_ids.end();
}

const std::vector<std::vector<size_t> >& BaseComplexQuad::GetNeighborComponentCidsGroups(SingularE& singularEdge) const
{
    if (!singularEdge.neighborComponentCidsGroups.empty())
        return singularEdge.neighborComponentCidsGroups;

    const std::vector<size_t> component_cids = GetNeighborComponentCids(singularEdge);
    size_t numOfneighbors = component_cids.size() / singularEdge.componentEids_link.size();
    std::vector<std::vector<size_t> > neighborComponentCids(numOfneighbors);
    for (size_t i = 0; i < numOfneighbors; ++i) {
        size_t start_component_cid = component_cids.at(i);
        neighborComponentCids[i].push_back(start_component_cid);
        for (size_t j = 1; j < singularEdge.componentEids_link.size(); ++j) {
            for (size_t k = j * numOfneighbors; k < (j + 1) * numOfneighbors; ++k) {
                const ComponentCell& start_component = componentC.at(start_component_cid);
                size_t next_component_cid = component_cids.at(k);
                const ComponentCell& next_component = componentC.at(next_component_cid);
                if (TwoComponentCellsHaveCommonFace(start_component, next_component)) {
                    start_component_cid = next_component_cid;
                    neighborComponentCids[i].push_back(start_component_cid);
                    break;
                }
            }
        }
    }

    singularEdge.neighborComponentCidsGroups = neighborComponentCids;
    return singularEdge.neighborComponentCidsGroups;
}

const std::vector<std::vector<size_t> >& BaseComplexQuad::GetNeighborComponentFidsGroups(SingularE& singularEdge) const
{
    if (!singularEdge.neighborComponentFidsGroups.empty())
        return singularEdge.neighborComponentFidsGroups;

    const std::vector<size_t> component_fids = GetNeighborComponentFids(singularEdge);
    size_t numOfneighbors = component_fids.size() / singularEdge.componentEids_link.size();
    std::vector<std::vector<size_t> > neighborComponentFids(numOfneighbors);
    for (size_t i = 0; i < numOfneighbors; ++i) {
        size_t start_component_fid = component_fids.at(i);
        neighborComponentFids[i].push_back(start_component_fid);
        for (size_t j = 1; j < singularEdge.componentEids_link.size(); ++j) {
            for (size_t k = j * numOfneighbors; k < (j + 1) * numOfneighbors; ++k) {
                const ComponentFace& start_component = componentF.at(start_component_fid);
                size_t next_component_fid = component_fids.at(k);
                const ComponentFace& next_component = componentF.at(next_component_fid);
                if (TwoComponentFacesHaveCommonEdge(start_component, next_component)) {
                    start_component_fid = next_component_fid;
                    neighborComponentFids[i].push_back(start_component_fid);
                    break;
                }
            }
        }
    }

    singularEdge.neighborComponentFidsGroups = neighborComponentFids;
    return singularEdge.neighborComponentFidsGroups;
}

std::vector<size_t> BaseComplexQuad::GetNeighborComponentCids(const SingularE& singularEdge) const
{
    const std::vector<size_t>& component_eids_link = GetComponentEidsLink((SingularE&)singularEdge);
    std::vector<size_t> component_cids;
    for (auto component_eid : component_eids_link) {
        const ComponentEdge& component_edge = componentE.at(component_eid);
        std::copy(component_edge.N_Cids.begin(), component_edge.N_Cids.end(), back_inserter(component_cids));
    }
    return component_cids;
}

std::vector<size_t> BaseComplexQuad::GetNeighborComponentFids(const SingularE& singularEdge) const
{
    const std::vector<size_t>& component_eids_link = GetComponentEidsLink((SingularE&)singularEdge);
    std::vector<size_t> component_fids;
    for (auto component_eid : component_eids_link) {
        const ComponentEdge& component_edge = componentE.at(component_eid);
        std::copy(component_edge.N_Fids.begin(), component_edge.N_Fids.end(), back_inserter(component_fids));
    }
    return component_fids;
}

void BaseComplexQuad::ExtractSingularVandE() {
    size_t singularVid = 0;
    for (auto& v : mesh.V) {
        if (v.isSingularity) {
            SingularV sv;
            sv.id = singularVid++;
            sv.id_mesh = v.id;
            sv.isBoundary = v.isBoundary;
            SingularityI.V.push_back(sv);
        }
    }
}

bool BaseComplexQuad::straight_line_test(const size_t h_v1, const size_t h_v2, const size_t h_v3)
{
    std::vector<size_t> nei = mesh.V[h_v1].N_Fids;
    std::vector<size_t> nei2 = mesh.V[h_v2].N_Fids;
    std::vector<size_t> nei3 = mesh.V[h_v3].N_Fids;

    for (int i = 0; i < mesh.V[h_v1].N_Fids.size(); i++)
        for (int j = 0; j < mesh.V[h_v3].N_Fids.size(); j++)
            if (mesh.V[h_v1].N_Fids[i] == mesh.V[h_v3].N_Fids[j]) return false;

    return true;
}

void BaseComplexQuad::AddcircularSingularE(const std::vector<SingularE> &circle_ses)
{
    if (circle_ses.empty()) return;

    std::vector<int> arrayv_test;
    for (int i = 0; i < mesh.V.size(); i++)
        arrayv_test.push_back(-2);
    for (int i = 0; i < SingularityI.E.size(); i++) {
        for (int j = 0; j < SingularityI.E[i].vs_link.size(); j++) {
            int vid = SingularityI.E[i].vs_link[j];
            arrayv_test[vid] = -1;
        }
    }
    for (int i = 0; i < SingularityI.V.size(); i++) {
        int vid = SingularityI.V[i].id_mesh;
        arrayv_test[vid] = -3;
    }
    for (int i = 0; i < circle_ses.size(); i++) {
        for (int j = 0; j < circle_ses[i].vs_link.size(); j++) {
            int vid = circle_ses[i].vs_link[j];
            arrayv_test[vid] = i;
        }
    }
    //test intersection with current singular nodes
    std::vector<SingularE> parallel_circle_ses;
    std::vector<std::vector<int> > Vfinds, Whichones;
    for (int i = 0; i < circle_ses.size(); i++) {
        std::vector<int> vfind;
        std::vector<int> which_ones;
        for (int j = 0; j < circle_ses[i].vs_link.size() - 1; j++) {
            int thisv = circle_ses[i].vs_link[j];
            const Vertex& hv = mesh.V[thisv];
            for (int k = 0; k < mesh.V[thisv].N_Eids.size(); k++) {
                bool found = false;
                int currente = mesh.V[thisv].N_Eids[k];
                if (currente == circle_ses[i].es_link[(1 + j) % circle_ses[i].es_link.size()] || currente == circle_ses[i].es_link[j]) continue;

                int currentv = thisv;
                int nextv = mesh.E[currente].Vids[0];
                if (nextv == thisv) nextv = mesh.E[currente].Vids[1];
                int cout_step = 0;
                std::vector<int> asery;
                asery.push_back(currentv);

                while (true) {
                    cout_step++;
                    asery.push_back(nextv);
                    if (arrayv_test[nextv] == -3) {
                        vfind.push_back(thisv);
                        which_ones.push_back(j);
                        found = true;
                        break;
                    }
                    bool foundnextv = false;
                    for (int m = 0; m < mesh.V[nextv].N_Vids.size(); m++) {
                        int nv = mesh.V[nextv].N_Vids[m];

                        if (nv != currentv) {
                            if (straight_line_test(nv, nextv, currentv) && arrayv_test[nv] != -1) {
                                currentv = nextv;
                                nextv = nv;
                                foundnextv = true;
                                break;
                            }
                        }
                    }
                    if (nextv == thisv || !foundnextv) break;
                }
                if (found) break;
            }
        }

        Vfinds.push_back(vfind);
        Whichones.push_back(which_ones);
    }

    //test intersection with no circle singular edge
    for (int i = 0; i < circle_ses.size(); i++) {
        if (Vfinds[i].size() != 0) continue;

        std::vector<int> vfind;
        std::vector<int> which_ones;
        for (int j = 0; j < circle_ses[i].vs_link.size() - 1; j++) {
            int thisv = circle_ses[i].vs_link[j];

            const Vertex& hv = mesh.V[thisv];
            for (int k = 0; k < mesh.V[thisv].N_Eids.size(); k++) {
                bool found = false;
                int currente = mesh.V[thisv].N_Eids[k];
                if (currente == circle_ses[i].es_link[(1 + j) % circle_ses[i].es_link.size()] || currente == circle_ses[i].es_link[j]) continue;

                int currentv = thisv;
                int nextv = mesh.E[currente].Vids[0];
                if (nextv == thisv) nextv = mesh.E[currente].Vids[1];
                int cout_step = 0;
                std::vector<int> asery;
                asery.push_back(currentv);

                while (true) {
                    cout_step++;
                    asery.push_back(nextv);
                    if (arrayv_test[nextv] == -1) {
                        vfind.push_back(thisv);
                        which_ones.push_back(j);
                        found = true;
                        break;
                    }
                    bool foundnextv = false;
                    for (int m = 0; m < mesh.V[nextv].N_Vids.size(); m++) {
                        int nv = mesh.V[nextv].N_Vids[m];
                        if (nv != currentv) {
                            if (straight_line_test(nv, nextv, currentv)) {
                                currentv = nextv;
                                nextv = nv;
                                foundnextv = true;
                                break;
                            }
                        }
                    }
                    if (nextv == thisv || !foundnextv) break;
                }
                if (found) break;
            }
        }

        Vfinds[i] = vfind;
        Whichones[i] = which_ones;
    }
    //if all the rest circle singular edge are the same length
    int length_circle = -1;
    int first_circle = -1;
    for (int i = 0; i < circle_ses.size(); i++) {
        if (Vfinds[i].size() == 0) {
            length_circle = circle_ses[i].es_link.size();
            first_circle = i;
            break;
        }
    }
    if (length_circle != -1) {
        for (int i = 0; i < circle_ses.size(); i++) {
            if (Vfinds[i].size() == 0) {
                if (length_circle != circle_ses[i].es_link.size()) ; //exit(0);
            }
        }
    }
    //test intersection with circle singular edge
    if (first_circle != -1) {
        std::vector<int> circlevfind, which_ones;
        int onesection = circle_ses[first_circle].es_link.size() / 3; //3 is the minimum length of circle could be!
        for (int j = 0; j < 3; j++) {
            int thisv = circle_ses[first_circle].vs_link[j * onesection];
            circlevfind.push_back(thisv);
            which_ones.push_back(first_circle);

            for (int k = 0; k < mesh.V[thisv].N_Eids.size(); k++) {
                int currente = mesh.V[thisv].N_Eids[k];
                if (currente == circle_ses[first_circle].es_link[(1 + j * onesection) % circle_ses[first_circle].es_link.size()]
                        || currente == circle_ses[first_circle].es_link[j * onesection]) continue;

                int currentv = thisv;
                int nextv = mesh.E[currente].Vids[0];
                if (nextv == thisv) nextv = mesh.E[currente].Vids[1];
                while (true) {
                    if (arrayv_test[nextv] >= 0 && arrayv_test[nextv] != first_circle) {
                        bool added = false;
                        for (int m = 0; m < which_ones.size(); m++)
                            if (which_ones[m] == arrayv_test[nextv]) added = true;
                        if (!added) {
                            circlevfind.push_back(nextv);
                            which_ones.push_back(arrayv_test[nextv]);
                        }
                    }
                    bool foundnextv = false;
                    for (int m = 0; m < mesh.V[nextv].N_Vids.size(); m++) {
                        int nv = mesh.V[nextv].N_Vids[m];
                        if (nv != currentv) {
                            if (straight_line_test(nv, nextv, currentv)) {
                                currentv = nextv;
                                nextv = nv;
                                foundnextv = true;
                                break;
                            }
                        }
                    }
                    if (nextv == thisv || !foundnextv) break;
                }
            }

            for (int k = 0; k < which_ones.size(); k++) {
                Vfinds[which_ones[k]].push_back(circlevfind[k]);
            }
        }
    }

//      Singular_V sv;
//      sv.index_hex=136;
//      sv.is_boundary=-1;
//      sv.index_own=SingularityI.V.size();
//      SingularityI.V.push_back(sv);
//      sv.index_hex=139;
//      sv.is_boundary=-1;
//      sv.index_own=SingularityI.V.size();
//      SingularityI.V.push_back(sv);

//      for(int i=0;i<circle_ses.size();i++)
//          SingularityI.E.push_back(circle_ses[i]);

    for (int i = 0; i < Vfinds.size(); i++) {
        for (int j = 0; j < Vfinds[i].size(); j++) {
            //continue;;

            SingularV sv;
            sv.id_mesh = Vfinds[i][j];
            sv.isBoundary = mesh.V[Vfinds[i][j]].isBoundary;
            sv.fake = true;
            sv.id = SingularityI.V.size();
            SingularityI.V.push_back(sv);

            //* for later use
            SingularE se;
            se.edge_type = i + 1;
            se.startend_Vid[0] = sv.id;
            se.startend_Vid[1] = sv.id + 1;
            if (j == Vfinds[i].size() - 1) se.startend_Vid[1] = sv.id - Vfinds[i].size() + 1;

            int start_p = -1, end_p = -1;
            for (int k = 0; k < circle_ses[i].vs_link.size() - 1; k++) {
                if (circle_ses[i].vs_link[k] == Vfinds[i][j]) start_p = k;
                if (circle_ses[i].vs_link[k] == Vfinds[i][(j + 1) % Vfinds[i].size()]) end_p = k;
            }
            int length_singulare = (end_p + circle_ses[i].vs_link.size() - 1 - start_p + 1) % (circle_ses[i].vs_link.size() - 1);

            for (int k = 0; k < length_singulare; k++)
                se.vs_link.push_back(circle_ses[i].vs_link[(k + start_p) % (circle_ses[i].vs_link.size() - 1)]);

            int start_e = -1, end_e = -1;
            for (int k = 0; k < circle_ses[i].es_link.size(); k++) {
                int e = circle_ses[i].es_link[k];
                if ((circle_ses[i].vs_link[start_p] == mesh.E[e].Vids[0]
                        && circle_ses[i].vs_link[(start_p + 1) % (circle_ses[i].vs_link.size() - 1)] == mesh.E[e].Vids[1])
                        || (circle_ses[i].vs_link[start_p] == mesh.E[e].Vids[1]
                                && circle_ses[i].vs_link[(start_p + 1) % (circle_ses[i].vs_link.size() - 1)] == mesh.E[e].Vids[0])) start_e = k;

                if ((circle_ses[i].vs_link[end_p] == mesh.E[e].Vids[0]
                        && circle_ses[i].vs_link[(end_p - 1 + circle_ses[i].vs_link.size() - 1) % (circle_ses[i].vs_link.size() - 1)] == mesh.E[e].Vids[1])
                        || (circle_ses[i].vs_link[end_p] == mesh.E[e].Vids[1]
                                && circle_ses[i].vs_link[(end_p - 1 + circle_ses[i].vs_link.size() - 1) % (circle_ses[i].vs_link.size() - 1)] == mesh.E[e].Vids[0])) end_e =
                        k;
            }

            length_singulare = (end_e - start_e + circle_ses[i].es_link.size() + 1) % circle_ses[i].es_link.size();
            for (int k = 0; k < length_singulare; k++)
                se.es_link.push_back(circle_ses[i].es_link[(k + start_e) % (circle_ses[i].es_link.size())]);

            SingularityI.E.push_back(se);
            //*/
        }
    }
}

void BaseComplexQuad::BuildSingularityConnectivity()
{
    for (size_t i = 0; i < SingularityI.E.size(); i++) {
        const size_t vid1 = SingularityI.E[i].startend_Vid[0];
        const size_t vid2 = SingularityI.E[i].startend_Vid[1];
        if (vid1 != vid2) {
            SingularityI.V[vid1].N_Vids.push_back(vid2);
            SingularityI.V[vid2].N_Vids.push_back(vid1);
            SingularityI.V[vid1].N_Eids.push_back(i);
            SingularityI.V[vid2].N_Eids.push_back(i);
        }
//        else {
//            // circular
//        }
    }
}

void BaseComplexQuad::BuildV()
{
    if (Fids.empty()) {
        std::cerr << "Fids of Base-Complex is empty()\n";
        return;
    }

    std::vector<bool> is_mesh_vertex_visited(mesh.V.size(), false);
    for (auto faceId : Fids)
        for (auto vertexId : mesh.F.at(faceId).Vids)
            is_mesh_vertex_visited[vertexId] = true;

    Vids.clear();
    Vids.reserve(is_mesh_vertex_visited.size());
    for (size_t i = 0; i < is_mesh_vertex_visited.size(); i++)
        if (is_mesh_vertex_visited[i]) Vids.push_back(i);
    Vids.resize(Vids.size());

    std::vector<bool> is_mesh_edge_in_an_edge_of_base_complex(mesh.E.size(), false);
    for (auto edgeId : Eids)
        is_mesh_edge_in_an_edge_of_base_complex[edgeId] = true;

    std::vector<bool> is_mesh_vertex_in_an_vertex_of_base_complex(mesh.V.size(), false);
    for (auto vertexId : Vids) {
        bool allNeighborEdgesInEids = true;
        for (auto edgeId : mesh.V.at(vertexId).N_Eids) {
            if (!is_mesh_edge_in_an_edge_of_base_complex[edgeId]) {
                allNeighborEdgesInEids = false;
                break;
            }
        }
        if (allNeighborEdgesInEids) is_mesh_vertex_in_an_vertex_of_base_complex[vertexId] = true;
    }

    Vids.clear();
    Vids.reserve(is_mesh_vertex_in_an_vertex_of_base_complex.size());
    for (size_t i = 0; i < is_mesh_vertex_in_an_vertex_of_base_complex.size(); i++)
        if (is_mesh_vertex_in_an_vertex_of_base_complex[i]) Vids.push_back(i);
    Vids.resize(Vids.size());

//    MeshFileWriter writer(mesh, "BaseComplexVertices.vtk");
//    writer.WriteVerticesVtk(Vids);
}

bool BaseComplexQuad::IsFaceMetSingularNode(const Face& face)
{
    for (size_t i = 0; i < face.Vids.size(); i++)
        for (size_t j = 0; j < SingularityI.V.size(); j++)
            if (SingularityI.V[j].id_mesh == face.Vids[i])
                return true;
    return false;
}

bool BaseComplexQuad::IsEdgeMetSingularNode(const Edge& edge)
{
    for (size_t i = 0; i < edge.Vids.size(); i++)
        for (size_t j = 0; j < SingularityI.V.size(); j++)
            if (SingularityI.V[j].id_mesh == edge.Vids[i])
                return true;
    return false;
}

bool BaseComplexQuad::IsEdgeMetSingularEdge(const Edge& edge)
{
    for (size_t j = 0; j < SingularityI.E.size(); j++)
        for (size_t i = 0; i < SingularityI.E[j].es_link.size(); i++)
            if (SingularityI.E[j].es_link[i] == edge.id)
                return true;
    return false;
}

size_t BaseComplexQuad::GetOppositEdgeId(const Face& face, const size_t edgeid)
{
    const Edge& edge = mesh.E.at(edgeid);
    const size_t edge_vid1 = edge.Vids[0];
    const size_t edge_vid2 = edge.Vids[1];
    for (size_t i = 0; i < face.Eids.size(); i++) {
        const size_t eid = face.Eids.at(i);
        const Edge& e = mesh.E.at(eid);
        const size_t e_vid1 = e.Vids[0];
        const size_t e_vid2 = e.Vids[1];
        if (e_vid1 != edge_vid1 && e_vid1 != edge_vid2 && e_vid2 != edge_vid1 && e_vid2 != edge_vid2)
            return eid;
    }

    return MAXID;
}

static bool IsTwoFaceSeperated(const Face& face1, const Face& face2)
{
    std::vector<size_t> vids = face1.Vids;
    std::copy(face2.Vids.begin(), face2.Vids.end(), back_inserter(vids));
    std::sort(vids.begin(), vids.end());
    auto iter = std::unique(vids.begin(), vids.end());
    return iter == vids.end();
}

size_t BaseComplexQuad::GetOppositFaceId(const Cell& cell, const size_t faceid)
{
    const Face& face = mesh.F.at(faceid);
    for (size_t i = 0; i < cell.Fids.size(); i++) {
        const size_t fid = cell.Fids.at(i);
        const Face& f = mesh.F.at(fid);
        if (IsTwoFaceSeperated(face, f))
            return fid;
    }

    return MAXID;
}

bool BaseComplexQuad::IsEdgeMetBaseComplexVertex(const Edge& edge)
{
    for (size_t i = 0; i < edge.Vids.size(); i++)
        for (size_t j = 0; j < Vids.size(); j++)
            if (componentV[j].id == edge.Vids[i])
                return true;
    return false;
}

bool BaseComplexQuad::IsEdgeOnBaseComplexEdge(const Edge& edge)
{
    for (auto edge_id : Eids)
        if (edge_id == edge.id) return true;
    return false;
}

bool BaseComplexQuad::IsFaceOnBaseComplexFace(const Face& face)
{
    for (auto face_id : Fids)
        if (face_id == face.id) return true;
    return false;
}

//bool BaseComplexQuad::IsOnBaseComplex(const Vertex& vertex) const
//{
//    return std::find(Vids.begin(), Vids.end(), vertex.id) != Vids.end();
//}
//
//bool BaseComplexQuad::IsOnBaseComplex(const Edge& edge) const
//{
//    return std::find(Eids.begin(), Eids.end(), edge.id) != Eids.end();
//}
//
//bool BaseComplexQuad::IsOnBaseComplex(const Face& face) const
//{
//    return std::find(Fids.begin(), Fids.end(), face.id) != Fids.end();
//}

bool BaseComplexQuad::IsOnBaseComplex(const Vertex& vertex) const {
    return hashVids.find(vertex.id) != hashVids.end();
}

bool BaseComplexQuad::IsOnBaseComplex(const Edge& edge) const {
    return hashEids.find(edge.id) != hashEids.end();
}

bool BaseComplexQuad::IsOnBaseComplex(const Face& face) const {
    return hashFids.find(face.id) != hashFids.end();
}

void BaseComplexQuad::TraceEdge(const std::vector<size_t>& mesh_edge_ids_on_singular_edge, const Face& start_face, std::vector<bool>& is_mesh_edge_visited)
{
    std::vector<size_t> mesh_edge_ids = mesh_edge_ids_on_singular_edge;
    Face* currentFace = (Face*) &start_face;
    bool is_edge_met_singularity = false;
    while (!is_edge_met_singularity) {
        std::cout << "Edge Id link : ";
        for (auto i : mesh_edge_ids) std::cout << i << " ";
        std::cout << std::endl;
        is_mesh_edge_visited[mesh_edge_ids[0]] = true;
        mesh_edge_ids[0] = GetOppositEdgeId(*currentFace, mesh_edge_ids[0]);
        //is_mesh_face_visited[currentFace->id] = true;

        is_edge_met_singularity = IsEdgeMetSingularNode(mesh.E.at(mesh_edge_ids[0]));
        if (!is_edge_met_singularity) is_edge_met_singularity = IsEdgeMetSingularEdge(mesh.E.at(mesh_edge_ids[0]));

        for (size_t j = 1; j < mesh_edge_ids.size(); j++) {
            const Edge& edge = mesh.E.at(mesh_edge_ids.at(j));
            for (size_t k = 0; k < edge.N_Fids.size(); k++) {
                const Face& next_face = mesh.F.at(edge.N_Fids.at(k));
                if (IsTwoFaceShareCommonEdge(*currentFace, next_face)) {
                    currentFace = (Face*) &next_face;
                    is_mesh_edge_visited[mesh_edge_ids[j]] = true;
                    mesh_edge_ids[j] = GetOppositEdgeId(*currentFace, mesh_edge_ids[j]);
                    //is_mesh_face_visited[currentFace->id] = true;
                    if (!is_edge_met_singularity) is_edge_met_singularity = IsEdgeMetSingularNode(mesh.E.at(mesh_edge_ids[j]));
                    if (!is_edge_met_singularity) is_edge_met_singularity = IsEdgeMetSingularEdge(mesh.E.at(mesh_edge_ids[j]));
                    break;
                }
            }
        }

        std::reverse(mesh_edge_ids.begin(), mesh_edge_ids.end());
        bool found = false;
        const Edge& edge = mesh.E.at(mesh_edge_ids[0]);
        for (size_t k = 0; k < edge.N_Fids.size(); k++) {
            const Face& face = mesh.F.at(edge.N_Fids.at(k));
            const size_t opposite_edge_id = GetOppositEdgeId(*currentFace, mesh_edge_ids[0]);
            //if ((!mesh.E.at(opposite_edge_id).isSingularity || !is_mesh_edge_visited[opposite_edge_id]) && !IsTwoFaceOnOneCell(mesh, face, *currentFace, edge)) {
            //if (!IsTwoFaceOnOneCell(mesh, face, *currentFace, edge)) {
            if (!is_mesh_edge_visited[opposite_edge_id] && !IsTwoFaceOnOneCell(mesh, face, *currentFace, edge)) {
                currentFace = (Face*) &face;
                found = true;
                break;
            }
        }
        if (!found) break;
    }
}

void BaseComplexQuad::TraceFace(const Face& start_face, std::vector<bool>& is_mesh_face_visited)
{
    std::vector<size_t> facePatchIds;
    facePatchIds.reserve(is_mesh_face_visited.size());
    facePatchIds.push_back(start_face.id);
    is_mesh_face_visited[start_face.id] = true;

    bool finished = false;
    while (!finished) {
        finished = true;
        for (auto faceId : facePatchIds) {
            const Face& face = mesh.F.at(faceId);
            for (auto edgeId : face.Eids) {
                const Edge& edge = mesh.E.at(edgeId);
                for (auto edgeNeighborFaceId : edge.N_Fids) {
                    const Face& edgeNeighborFace = mesh.F.at(edgeNeighborFaceId);
                    if (!is_mesh_face_visited[edgeNeighborFaceId] && !IsTwoFaceOnOneCell(mesh, face, edgeNeighborFace, edge)) {
                        is_mesh_face_visited[edgeNeighborFaceId] = true;
                        facePatchIds.push_back(edgeNeighborFaceId);
                        finished = false;
                    }
                }
            }
        }
    }
}

void BaseComplexQuad::TraceFace(const Face& start_face, std::vector<bool>& is_mesh_face_visited, std::vector<size_t>& fids)
{
    std::vector<size_t>& facePatchIds = fids;
    facePatchIds.reserve(is_mesh_face_visited.size());
    facePatchIds.push_back(start_face.id);
    is_mesh_face_visited[start_face.id] = true;

    bool finished = false;
    while (!finished) {
        finished = true;
        for (auto faceId : facePatchIds) {
            const Face& face = mesh.F.at(faceId);
            for (auto edgeId : face.Eids) {
                const Edge& edge = mesh.E.at(edgeId);
                for (auto edgeNeighborFaceId : edge.N_Fids) {
                    const Face& edgeNeighborFace = mesh.F.at(edgeNeighborFaceId);
                    if (!is_mesh_face_visited[edgeNeighborFaceId] && !IsTwoFaceOnOneCell(mesh, face, edgeNeighborFace, edge)
                        && !edge.isSingularity) {
                        is_mesh_face_visited[edgeNeighborFaceId] = true;
                        facePatchIds.push_back(edgeNeighborFaceId);
                        finished = false;
                    }
                }
            }
        }
    }
}

size_t BaseComplexQuad::TraceAlongEdge(const Vertex& start_vertex, const Edge& start_edge,
        std::vector<bool>& is_edge_visited, std::vector<size_t>& vids_link, std::vector<size_t>& eids_link) {
    Vertex* currentVertex = (Vertex*) &start_vertex;
    Edge* currentEdge = (Edge*) &start_edge;
    vids_link.push_back(currentVertex->id);
    std::vector<bool> is_mesh_vertex_visited(mesh.V.size(), false);
    is_mesh_vertex_visited[currentVertex->id] = true;
    while (currentEdge != NULL) {
        eids_link.push_back(currentEdge->id);
        is_edge_visited[currentEdge->id] = true;
        currentVertex = (Vertex*)&mesh.V.at(currentEdge->Vids[0] == currentVertex->id ? currentEdge->Vids[1] : currentEdge->Vids[0]);
        vids_link.push_back(currentVertex->id);
        if (currentVertex->isBoundary || currentVertex->isSingularity || is_mesh_vertex_visited[currentVertex->id]) break;
        is_mesh_vertex_visited[currentVertex->id] = true;
        currentEdge = GetNextContinuousEdge(*currentVertex, *currentEdge);
    }

    return vids_link.back();
}


void BaseComplexQuad::BuildE() {
    std::vector<bool> is_singular_vertex_visited(SingularityI.V.size(), false);
    std::vector<bool> is_mesh_edge_visited(mesh.E.size(), false);
    for (auto& sv : SingularityI.V) {
        is_singular_vertex_visited[sv.id] = true;
        const Vertex& vertex = mesh.V.at(sv.id_mesh);
        for (auto edgeid : vertex.N_Eids) {
            const Edge& edge = mesh.E.at(edgeid);
            if (!is_mesh_edge_visited[edgeid]) {
                std::vector<size_t> vids_link;
                std::vector<size_t> eids_link;
                TraceAlongEdge(vertex, edge, is_mesh_edge_visited, vids_link, eids_link);
                //SingularityI.E[i].separatedComponentFids.push_back(fids);
                separatedVertexIdsLink.push_back(vids_link);
                separatedEdgeIdsLink.push_back(eids_link);
            }
        }
    }

    std::vector<bool> is_mesh_vertex_visited(mesh.V.size(), false);
    Eids.clear();
    Eids.reserve(is_mesh_edge_visited.size());
    for (size_t i = 0; i < is_mesh_edge_visited.size(); i++)
        if (is_mesh_edge_visited[i] || mesh.E.at(i).isBoundary) {
            Eids.push_back(i);
            const Edge& edge = mesh.E.at(i);
            is_mesh_vertex_visited.at(edge.Vids[0]) = true;
            is_mesh_vertex_visited.at(edge.Vids[1]) = true;
        }
    Eids.resize(Eids.size());

    Vids.clear();
    Vids.reserve(mesh.V.size());
    for (size_t i = 0; i < is_mesh_vertex_visited.size(); i++)
        if (is_mesh_vertex_visited[i]) Vids.push_back(i);
    Vids.resize(Vids.size());

    Fids.resize(mesh.F.size());
    for (size_t i = 0; i < Fids.size(); i++)
        Fids[i] = i;
}

void BaseComplexQuad::BuildF() {
    hashVids.insert(Vids.begin(), Vids.end());
    hashEids.insert(Eids.begin(), Eids.end());
    hashFids.insert(Fids.begin(), Fids.end());
    BuildComponentF();
    BuildComponentColor();
    BuildComponentE();
    BuildComponentV();
}

void BaseComplexQuad::BuildComponentV()
{
    for (size_t i = 0; i < Vids.size(); i++) {
        const size_t vertexId = Vids.at(i);
        mesh.V.at(vertexId).component_id = i;
        Vertex v = mesh.V.at(vertexId).xyz();
        v.id = i;
        V.push_back(v);
    }
}

void BaseComplexQuad::BuildComponentE()
{
    std::vector<bool> is_edge_visited(mesh.E.size(), false);
    int eid = 0;
    for (auto& edgeid : Eids) {
        const Edge& edge = mesh.E.at(edgeid);
        if (is_edge_visited[edgeid]) continue;
        std::vector<size_t> eids;
        RegionGrowEdge(edge, is_edge_visited, eids);
        for (auto edge_id : eids)
            mesh.E.at(edge_id).componentEid = eid;
        ComponentEdge componentEdge;
        //componentCell.Vids = c.Vids;
        componentEdge.eids_link = eids;
        componentEdge.id = eid;
        sort(eids.begin(), eids.end());
        BuildComponentEdgeLink(componentEdge);
        componentE.push_back(componentEdge);
        ++eid;
    }
}

void BaseComplexQuad::BuildComponentEdgeLink(ComponentEdge& componentEdge)
{
    std::vector<size_t> vids;
    vids.reserve(componentEdge.eids_link.size() * 2);
    for (auto eid : componentEdge.eids_link) {
        std::copy(mesh.E.at(eid).Vids.begin(), mesh.E.at(eid).Vids.end(), back_inserter(vids));
    }

    size_t start_vid = vids[0];
    std::sort(vids.begin(), vids.end());
    std::vector<size_t> vids_copy = vids;
    vids_copy.resize(std::distance(vids_copy.begin(), std::unique(vids_copy.begin(), vids_copy.end())));
    if(vids_copy.size() != componentEdge.eids_link.size()) {
        // not circular
        std::vector<bool> count(mesh.V.size(), false);
        for (auto vid : vids)
            count[vid] = !count[vid];
        for (auto vid : vids) {
            if (count[vid]) {
                start_vid = vid;
                break;
            }
        }
    }
    componentEdge.vids_link = vids_copy;

    const Vertex& start_vertex = mesh.V.at(start_vid);
    size_t start_eid = 0;
    for (auto neigbor_eid : start_vertex.N_Eids) {
        bool found = false;
        for (auto eid : componentEdge.eids_link)
            if (neigbor_eid == eid) {
                start_eid = eid;
                found = true;
                break;
            }
        if (found) break;
    }
    const Edge& start_edge = mesh.E.at(start_eid);

    std::vector<size_t> vids_link;
    std::vector<size_t> eids_link;
    std::vector<bool> is_edge_visited(mesh.E.size(), false);
    TraceVertex(componentEdge, start_vertex, start_edge, is_edge_visited, vids_link, eids_link);
    componentEdge.vids_link = vids_link;
    componentEdge.eids_link = eids_link;
}

size_t BaseComplexQuad::TraceVertex(ComponentEdge& componentEdge, const Vertex& start_vertex, const Edge& start_edge,
        std::vector<bool>& is_edge_visited, std::vector<size_t>& vids_link, std::vector<size_t> &eids_link) // return end_vertex_id
{
    vids_link.push_back(start_vertex.id);
    eids_link.push_back(start_edge.id);
    is_edge_visited[start_edge.id] = true;
    size_t next_vertex_id = start_edge.Vids[0] == start_vertex.id ? start_edge.Vids[1] : start_edge.Vids[0];
    while (true) {
        const Vertex& next_vertex = mesh.V.at(next_vertex_id);
        bool found_next_edge = false;
        for (auto edgeId : next_vertex.N_Eids) {
            bool is_edgeId_in_eids_link = std::find(componentEdge.eids_link.begin(), componentEdge.eids_link.end(), edgeId) != componentEdge.eids_link.end();
            if (is_edgeId_in_eids_link && !is_edge_visited[edgeId]) {
                const Edge& next_edge = mesh.E.at(edgeId);
                next_vertex_id = next_edge.Vids[0] == next_vertex_id ? next_edge.Vids[1] : next_edge.Vids[0];
                vids_link.push_back(next_vertex.id);
                eids_link.push_back(next_edge.id);
                is_edge_visited[next_edge.id] = true;
                found_next_edge = true;
                break;
            }
        }
        if (!found_next_edge) break;
    }
    vids_link.push_back(next_vertex_id);
    return next_vertex_id;
}

void BaseComplexQuad::BuildComponentF() {
    std::vector<bool> is_face_visited(mesh.F.size(), false);
    int fid = 0;
    for (auto& faceid : Fids) {
        const Face& face = mesh.F.at(faceid);
        if (is_face_visited[face.id]) continue;
        std::vector<size_t> fids;
        RegionGrowFace(face, is_face_visited, fids);
        for (auto face_id : fids) {
            mesh.F.at(face_id).componentFid = fid;
            mesh.F.at(face_id).component_id = fid;
        }
        ComponentFace componentFace;
        //componentCell.Vids = c.Vids;
        componentFace.fids_patch = fids;
        componentFace.id = fid;
        componentF.push_back(componentFace);
        ++fid;
    }
}

void BaseComplexQuad::BuildComponentC()
{
    std::vector<bool> is_cell_visited(mesh.C.size(), false);
    int cid = 0;
    for (auto& cell : mesh.C) {
        if (is_cell_visited[cell.id]) continue;
        std::vector<size_t> cids;
        RegionGrowCell(cell, is_cell_visited, cids);
        std::copy(cids.begin(), cids.end(), back_inserter(Cids));
        for (auto cellid : cids) {
            ((Cell&)mesh.C.at(cellid)).component_id = cid;
        }
//        std::cout << "component " << cid << " : ";
//        for (auto id : cids) std::cout << id << " ";
//        std::cout << std::endl;
        for (auto cellid : cids)
            mesh.C.at(cellid).componentCid = cid;
        ComponentCell componentCell;
        //componentCell.Vids = c.Vids;
        componentCell.cids_patch = cids;
        componentCell.id = cid;
        componentC.push_back(componentCell);
        ++cid;
    }
}
//
//void BaseComplexQuad::BuildComponentE()
//{
//    size_t eid = 0;
//    std::vector<bool> is_edge_visited(mesh.E.size(), false);
//    for (size_t i = 0; i < Vids.size(); i++) {
//        const size_t start_vertex_id_in_Vids = i;
//        const size_t vertexId = Vids.at(i);
//        const Vertex& start_vertex = mesh.V.at(vertexId);
//        for (auto edgeId : start_vertex.N_Eids) {
//            bool is_edgeId_in_Eids = std::find(Eids.begin(), Eids.end(), edgeId) != Eids.end();
//            if (is_edgeId_in_Eids && !is_edge_visited[edgeId]) {
//                const Edge& start_edge = mesh.E.at(edgeId);
//                std::vector<size_t> vids_link;
//                std::vector<size_t> eids_link;
//                const size_t end_vertex_id = TraceVertex(start_vertex, start_edge, is_edge_visited, vids_link, eids_link);
//                const size_t end_vertex_id_in_Vids = std::distance(Vids.begin(), std::find(Vids.begin(), Vids.end(), end_vertex_id));
//                Edge e(2);
//                e.Vids[0] = start_vertex_id_in_Vids;
//                e.Vids[1] = end_vertex_id_in_Vids;
//                e.id = eid;
//                E.push_back(e);
//
//                ComponentEdge componentEdge;
//                componentEdge.Vids = e.Vids;
//                componentEdge.vids_link = vids_link;
//                componentEdge.eids_link = eids_link;
//                componentEdge.id = eid;
//                componentE.push_back(componentEdge);
//                ++eid;
//            }
//        }
//    }
//}
//
//void BaseComplexQuad::BuildComponentF()
//{
//    size_t fid = 0;
//    std::vector<bool> is_face_visited(mesh.F.size(), false);
//    for (size_t i = 0; i < componentE.size(); i++) {
//        for (auto edgeId : componentE[i].eids_link) {
//            const Edge& start_edge = mesh.E.at(edgeId);
//            for (auto faceId : start_edge.N_Fids) {
//                const Face& start_face = mesh.F.at(faceId);
//                bool is_faceId_in_Fids = std::find(Fids.begin(), Fids.end(), faceId) != Fids.end();
//                if (is_faceId_in_Fids && !is_face_visited[start_face.id]) {
//                    std::vector<size_t> fids_patch;
//                    const size_t end_component_eid = TraceEdge(componentE[i].eids_link, start_face, is_face_visited, fids_patch);
//                    Face f(4);
//                    f.Vids[0] = componentE[i].vids_link.front();
//                    f.Vids[1] = componentE[i].vids_link.back();
//                    f.Vids[2] = componentE[end_component_eid].vids_link.front();
//                    f.Vids[3] = componentE[end_component_eid].vids_link.back();
//                    bool foundEdgeWithFVids1and2 = false;
//                    for (const auto& edge : componentE) {
//                        if ((edge.vids_link.front() == f.Vids[1] && edge.vids_link.back() == f.Vids[2]) || (edge.vids_link.back() == f.Vids[1] && edge.vids_link.front() == f.Vids[2])) {
//                            foundEdgeWithFVids1and2 = true;
//                            break;
//                        }
//                    }
//                    if (!foundEdgeWithFVids1and2) std::swap(f.Vids[2], f.Vids[3]);
//                    f.id = fid;
//                    F.push_back(f);
//
//                    ComponentFace componentFace;
//                    componentFace.Vids = f.Vids;
//                    componentFace.fids_patch = fids_patch;
//                    componentFace.id = fid;
//                    componentF.push_back(componentFace);
//                    ++fid;
//                }
//            }
//        }
//    }
//}
//
//void BaseComplexQuad::BuildComponentC()
//{
//    size_t cid = 0;
//    std::vector<bool> is_cell_visited(mesh.C.size(), false);
//    for (size_t i = 0; i < componentF.size(); i++) {
//        for (auto faceId : componentF[i].fids_patch) {
//            const Face& start_face = mesh.F.at(faceId);
//            for (auto cellId : start_face.N_Cids) {
//                const Cell& start_cell = mesh.C.at(cellId);
//                if (!is_cell_visited[start_cell.id]) {
//                    std::vector<size_t> cids_patch;
//                    const size_t end_component_fid = TraceFace(componentF[i].fids_patch, start_cell, is_cell_visited, cids_patch);
////                    std::sort(cids_patch.begin(), cids_patch.end());
////                    cids_patch.resize(std::distance(cids_patch.begin(), std::unique(cids_patch.begin(), cids_patch.end())));
//                    std::copy(cids_patch.begin(), cids_patch.end(), back_inserter(Cids));
////                    std::cout << "component " << cid << " : start_face_id : " << componentF[i].fids_patch[0] << " : cell_id :";
////                    for (auto id : cids_patch) std::cout << id << " ";
////                    std::cout << std::endl;
//                    Cell c(8);
////                    c.Vids[0] = componentE[i].vids_link.front();
////                    c.Vids[1] = componentE[i].vids_link.back();
////                    c.Vids[2] = componentE[end_component_eid].vids_link.front();
////                    c.Vids[3] = componentE[end_component_eid].vids_link.back();
////                    bool foundEdgeWithFVids1and2 = false;
////                    for (const auto& edge : componentE) {
////                        if ((edge.vids_link.front() == f.Vids[1] && edge.vids_link.back() == f.Vids[2]) || (edge.vids_link.back() == f.Vids[1] && edge.vids_link.front() == f.Vids[2])) {
////                            foundEdgeWithFVids1and2 = true;
////                            break;
////                        }
////                    }
////                    if (!foundEdgeWithFVids1and2) std::swap(f.Vids[2], f.Vids[3]);
////                    f.id = cid;
////                    F.push_back(f);
//
//                    ComponentCell componentCell;
//                    componentCell.Vids = c.Vids;
//                    componentCell.cids_patch = cids_patch;
//                    componentCell.id = cid;
//                    componentC.push_back(componentCell);
//                    ++cid;
//                }
//            }
//        }
//    }
//}

size_t BaseComplexQuad::TraceVertex(const Vertex& start_vertex, const Edge& start_edge,
        std::vector<bool>& is_edge_visited, std::vector<size_t>& vids_link, std::vector<size_t> &eids_link) // return end_vertex_id
{
    vids_link.push_back(start_vertex.id);
    eids_link.push_back(start_edge.id);
    is_edge_visited[start_edge.id] = true;
    size_t next_vertex_id = start_edge.Vids[0] == start_vertex.id ? start_edge.Vids[1] : start_edge.Vids[0];
    bool is_next_vertex_id_in_Vids = std::find(Vids.begin(), Vids.end(), next_vertex_id) != Vids.end();
    while (!is_next_vertex_id_in_Vids) {
        const Vertex& next_vertex = mesh.V.at(next_vertex_id);
        for (auto edgeId : next_vertex.N_Eids) {
            bool is_edgeId_in_Eids = std::find(Eids.begin(), Eids.end(), edgeId) != Eids.end();
            if (is_edgeId_in_Eids && !is_edge_visited[edgeId]) {
                const Edge& next_edge = mesh.E.at(edgeId);
                next_vertex_id = next_edge.Vids[0] == next_vertex_id ? next_edge.Vids[1] : next_edge.Vids[0];
                vids_link.push_back(next_vertex.id);
                eids_link.push_back(next_edge.id);
                is_edge_visited[next_edge.id] = true;
                break;
            }
        }
        is_next_vertex_id_in_Vids = std::find(Vids.begin(), Vids.end(), next_vertex_id) != Vids.end();
    }
    vids_link.push_back(next_vertex_id);
    return next_vertex_id;
}

/*
                ---------
                |       |
                |   *<--|----  currentFace
                |       |
currentEdge --> ---------
                |       |
                |   *<--|----  nextFace
                |       |
                ---------
 */
Face* BaseComplexQuad::GetNextFace(const Edge* currentEdge, const Face* currentFace, const std::vector<bool>& is_mesh_edge_visited)
{
    for (size_t k = 0; k < currentEdge->N_Fids.size(); k++) {
        const Face& face = mesh.F.at(currentEdge->N_Fids.at(k));
        const size_t opposite_edge_id = GetOppositEdgeId(*currentFace, currentEdge->id);
        if (!is_mesh_edge_visited[opposite_edge_id] && !IsTwoFaceOnOneCell(mesh, face, *currentFace, *currentEdge)) {
            return (Face*) &face;
        }
    }
    std::cerr << "GetNextFace Error!\n";
    return NULL;
}

Face* BaseComplexQuad::GetNextFace(const Edge& currentEdge, const Face& currentFace)
{
    for (size_t k = 0; k < currentEdge.N_Fids.size(); k++) {
        const Face& face = mesh.F.at(currentEdge.N_Fids.at(k));
        //const size_t opposite_edge_id = GetOppositEdgeId(currentFace, currentEdge.id);
        if (!IsTwoFaceOnOneCell(mesh, face, currentFace, currentEdge)) {
            return (Face*)&face;
        }
    }
    // std::cerr << "GetNextFace Error!\n";
    return NULL;
}


Edge* BaseComplexQuad::GetNextEdge(const Vertex& currentVertex, const Edge& currentEdge)
{
    for (size_t k = 0; k < currentVertex.N_Eids.size(); k++) {
        const Edge& edge = mesh.E.at(currentVertex.N_Eids.at(k));
        if (!IsTwoEdgesOnOneFace(mesh, edge, currentEdge, currentVertex)) {
            return (Edge*)&edge;
        }
    }
    // std::cerr << "GetNextEdge Error!\n";
    return NULL;
}

static bool IsTwoEdgeHasCommonVertex(const Edge& edge1, const Edge& edge2)
{
    return (edge1.Vids[0] == edge2.Vids[0] || edge1.Vids[0] == edge2.Vids[1]
          ||edge1.Vids[1] == edge2.Vids[0] || edge1.Vids[1] == edge2.Vids[1]);
}
/*
                ---------
                |       |
currentEdge --> |   *<--|----  currentFace
                |       |
 commonEdge --> ---------
                |       |
   nextEdge --> |   *<--|----  nextFace
                |       |
                ---------
 */

const Face& BaseComplexQuad::GetNextFace(const Edge& currentEdge, const Face& currentFace, const Edge& nextEdge)
{
    size_t commonEdgeId = currentFace.Eids[0];
    for (auto edgeid : currentFace.Eids) {
        if (edgeid != currentEdge.id && IsTwoEdgeHasCommonVertex(mesh.E.at(edgeid), nextEdge)) {
            commonEdgeId = edgeid;
            break;
        }
    }
    Face* nextFace = GetNextFace(mesh.E.at(commonEdgeId), currentFace);
    return *nextFace;
}

static bool IsTwoFaceHaveCommonEdge(const Face& face1, const Face& face2)
{
    std::vector<size_t> vids = face1.Vids;
    std::copy(face2.Vids.begin(), face2.Vids.end(), back_inserter(vids));
    std::sort(vids.begin(), vids.end());
    auto iter = std::unique(vids.begin(), vids.end());
    return std::distance(vids.begin(), iter) == 6;
}

static bool IsTwoCellsHaveCommonFace(const Cell& cell1, const Cell& cell2)
{
    std::vector<size_t> vids = cell1.Vids;
    std::copy(cell2.Vids.begin(), cell2.Vids.end(), back_inserter(vids));
    std::sort(vids.begin(), vids.end());
    auto iter = std::unique(vids.begin(), vids.end());
    return std::distance(vids.begin(), iter) == 12;
}

/*
                  3__________________2___________________
                  /|                 /|                 /|
currentFace -----/-|----->*         / |       *<-------/-|----- nextFace
                /  |               /  |               /  |
            0  /___|_____________1/___|______________/   |
               |   |              |   |              |   |
               |   |              |   |              |   |
currentCell ---|---|----->*       |   |       *<-----|---|----- nextCell
               |   |______________|___|______________|___|
               |   / 7            |  /               |  /
               |  /               | /                | /
               | /                |/                 |/
               |/_________________/__________________/
             4                    5
*/
const Cell& BaseComplexQuad::GetNextCell(const Face& currentFace, const Cell& currentCell, const Face& nextFace)
{
    size_t commonFaceId = currentCell.Fids[0];
    for (auto faceid : currentCell.Fids) {
        if (faceid != currentFace.id && IsTwoFaceHaveCommonEdge(mesh.F.at(faceid), nextFace)) {
            commonFaceId = faceid;
            break;
        }
    }
    const Face& commonFace = mesh.F.at(commonFaceId);
    size_t nextCellId = commonFace.N_Cids[0] == currentCell.id ? commonFace.N_Cids[1] : commonFace.N_Cids[0];
    const Cell& nextCell = mesh.C.at(nextCellId);
    return nextCell;
}

size_t BaseComplexQuad::TraceEdge(const Edge& start_edge, const Face& start_face, std::vector<bool>& is_face_visited, std::vector<size_t>& fids_patch)
{
    Edge* currentEdge = (Edge*) &start_edge;
    Face* currentFace = (Face*) &start_face;
    bool is_met_basecomplex_edge = false;
    std::vector<bool> is_mesh_edge_visited(mesh.E.size(), false);
    while (!is_met_basecomplex_edge && currentFace != NULL) {
        fids_patch.push_back(currentFace->id);
        is_face_visited[currentFace->id] = true;
        is_mesh_edge_visited[currentEdge->id] = true;
        currentEdge = (Edge*)&mesh.E.at(GetOppositEdgeId(*currentFace, currentEdge->id));
        // currentFace = GetNextFace(currentEdge, currentFace, is_mesh_edge_visited);
        currentFace = GetNextFace(*currentEdge, *currentFace);
        is_met_basecomplex_edge = IsEdgeOnBaseComplexEdge(*currentEdge);
    }

    return currentEdge->id;
}

size_t BaseComplexQuad::TraceEdge(const std::vector<size_t>& start_eids_link, const Face& start_face, std::vector<bool>& is_face_visited, std::vector<size_t>& fids_patch)
{
    Face* currentFace = (Face*) &start_face;
    size_t end_compoent_edge_id = 0;
    for (size_t i = 0; i < start_eids_link.size(); i++) {
        const Edge& currentEdge = mesh.E.at(start_eids_link.at(i));
        const size_t end_edge_id = TraceEdge(currentEdge, *currentFace, is_face_visited, fids_patch);
        if (i == 0) {
            for (const auto& component_edge : componentE) {
                if (component_edge.eids_link.front() == end_edge_id || component_edge.eids_link.back() == end_edge_id) {
                    end_compoent_edge_id = component_edge.id;
                    break;
                }
            }
        }
        if (i < start_eids_link.size() - 1) {
            const Edge& nextEdge = mesh.E.at(start_eids_link.at(i + 1));
            const Face& nextFace = GetNextFace(currentEdge, *currentFace, nextEdge);
            currentFace = (Face*)&nextFace;
        }
    }

    return end_compoent_edge_id;
}

size_t BaseComplexQuad::TraceFace(const std::vector<size_t>& start_fids_patch, const Cell& start_cell, std::vector<bool>& is_cell_visited, std::vector<size_t>& cids_patch)
{
    Cell* currentCell = (Cell*) &start_cell;
    Face* currentFace = (Face*)&mesh.F.at(start_fids_patch.at(0));
    size_t end_compoent_face_id = 0;
    for (size_t i = 0; i < start_fids_patch.size(); i++) {
//        Face* currentFace = (Face*)&mesh.F.at(start_fids_patch.at(i));
        const size_t end_face_id = TraceFace(*currentFace, *currentCell, is_cell_visited, cids_patch);
        if (i == 0) {
            for (const auto& component_face : componentF) {
                bool found = false;
                for (auto fid : component_face.fids_patch) {
                    if (fid == end_face_id || fid == end_face_id) {
                        end_compoent_face_id = component_face.id;
                        found = true;
                        break;
                    }
                }
                if (found) break;
            }
        }
//        if (i < start_fids_patch.size() - 1) {
//            const Face& nextFace = mesh.F.at(start_fids_patch.at(i + 1));
//            const Cell& nextCell = GetNextCell(currentFace, *currentCell, nextFace);
//            currentCell = (Face*)&nextCell;
//        }
        // Not efficient
        if (i < start_fids_patch.size() - 1) {
            const Face& nextFace = mesh.F.at(start_fids_patch.at(i + 1));
            for (auto nextCellId : nextFace.N_Cids) {
                const Cell& nextCell = mesh.C.at(nextCellId);
                bool found = false;
                for (auto cellid : cids_patch) {
                    const Cell& cell = mesh.C.at(cellid);
                    if (IsTwoCellsHaveCommonFace(nextCell, cell)) {
                        currentFace = (Face*)&nextFace;
                        currentCell = (Cell*)&nextCell;
                        found = true;
                        break;
                    }
                }
                if (found) break;
            }
        }
    }

    return end_compoent_face_id;
}

size_t BaseComplexQuad::TraceFace(const Face& start_face, const Cell& start_cell, std::vector<bool>& is_cell_visited, std::vector<size_t>& cids_patch)
{
    Face* currentFace = (Face*) &start_face;
    Cell* currentCell = (Cell*) &start_cell;
    bool is_met_basecomplex_face = false;
    std::vector<bool> is_mesh_face_visited(mesh.F.size(), false);
    while (!is_met_basecomplex_face && currentCell != NULL) {
        cids_patch.push_back(currentCell->id);
        is_cell_visited[currentCell->id] = true;
        is_mesh_face_visited[currentFace->id] = true;
        currentFace = (Face*)&mesh.F.at(GetOppositFaceId(*currentCell, currentFace->id));
        if (currentFace->N_Cids.size() == 1) {
            currentCell = NULL;
            continue;
        }
        const size_t cellId = currentFace->N_Cids[0] == currentCell->id ? currentFace->N_Cids[1] : currentFace->N_Cids[0];
        currentCell = (Cell*)&mesh.C.at(cellId);
        is_met_basecomplex_face = IsFaceOnBaseComplexFace(*currentFace);
    }

    return currentFace->id;
}

void BaseComplexQuad::BuildComponentColor()
{
    if (componentF.size() == 0) return;
    componentF[0].color = 0;

    for (auto faceid : componentF[0].fids_patch) {
        auto& face = mesh.F.at(faceid);
        face.color = 0;
    }
    int last_color = 0;
    for (auto& component : componentF) {
        if (component.color >= 0)
            continue;
        std::vector<char> colors = GetColorsOfNeighborComponents(component);
        int color = 0;
        int j = 0;
        for (int i = 0; i < colors.size(); ++i)
            if (colors[i] == -1) {
                color = i;
                break;
            }
        component.color = color;
        for (auto faceid : component.fids_patch) {
            auto& face = (Face&) mesh.F.at(faceid);
            face.color = color;
        }
    }
}

std::vector<char> BaseComplexQuad::GetColorsOfNeighborComponents(const ComponentCell& component)
{
    std::vector<char> colors(8, -1);
//    int color = 0;
    for (const auto cellid : component.cids_patch) {
        const Cell& cell = mesh.C.at(cellid);
        for (auto faceid : cell.Fids) {
            const Face& face = mesh.F.at(faceid);
            for (auto neighbor_cellid : face.N_Cids) {
                const Cell& neighbor_cell = mesh.C.at(neighbor_cellid);
                if (neighbor_cell.color >= 0) colors[neighbor_cell.color] = 0;
            }
        }
    }
    int last_color = 0;
    for (const auto cellid : component.cids_patch) {
        const Cell& cell = mesh.C.at(cellid);
        for (auto faceid : cell.Fids) {
            const Face& face = mesh.F.at(faceid);
            for (auto neighbor_cellid : face.N_Cids) {
                const Cell& neighbor_cell = mesh.C.at(neighbor_cellid);
                if (neighbor_cell.color >= 0) colors[neighbor_cell.color] = 0;
                else {
                    int color = 0;
                    for (int i = 0; i < colors.size(); ++i)
                        if (colors[i] == -1) {
                            color = i;
                            last_color = i;
                            colors[i] = 0;
                            break;
                        }
                    for (auto cid : componentC.at(neighbor_cell.component_id).cids_patch) {
                        Cell& comCell = (Cell&)mesh.C.at(cid);
                        comCell.color = color;
                    }
                }
            }
        }
    }
    return colors;
}

std::vector<char> BaseComplexQuad::GetColorsOfNeighborComponents(const ComponentFace& component) {
    std::vector<char> colors(4, -1);
    for (const auto faceid : component.fids_patch) {
        const Face& face = mesh.F.at(faceid);
        for (auto edgeid : face.Eids) {
            const Edge& edge = mesh.E.at(edgeid);
            for (auto neighbor_faceid : edge.N_Fids) {
                const Face& neighbor_face = mesh.F.at(neighbor_faceid);
                if (neighbor_face.color >= 0) colors[neighbor_face.color] = 0;
            }
        }
    }
    int last_color = 0;
    for (const auto faceid : component.fids_patch) {
        const Face& face = mesh.F.at(faceid);
        for (auto edgeid : face.Eids) {
            const Edge& edge = mesh.E.at(edgeid);
                for (auto neighbor_faceid : edge.N_Fids) {
                const Face& neighbor_face = mesh.F.at(neighbor_faceid);
                if (neighbor_face.color >= 0) colors[neighbor_face.color] = 0;
                else {
                    int color = 0;
                    for (int i = 0; i < colors.size(); ++i)
                        if (colors[i] == -1) {
                            color = i;
                            last_color = i;
                            colors[i] = 0;
                            break;
                        }
                    for (auto fid : componentF.at(neighbor_face.component_id).fids_patch) {
                        Face& comFace = (Face&)mesh.C.at(fid);
                        comFace.color = color;
                    }
                }
            }
        }
    }
    return colors;
}

size_t BaseComplexQuad::RegionGrowCell(const Cell& start_cell, std::vector<bool>& is_cell_visited, std::vector<size_t>& cids)
{
    is_cell_visited[start_cell.id] = true;
    cids.push_back(start_cell.id);
    for (auto face_id : start_cell.Fids) {
        const Face& face = mesh.F.at(face_id);
        if (!IsOnBaseComplex(face)) {
            const size_t nextCellId = face.N_Cids[0] == start_cell.id ? face.N_Cids[1] : face.N_Cids[0];
            if (!is_cell_visited[nextCellId]) {
                const Cell& nextCell = mesh.C.at(nextCellId);
                RegionGrowCell(nextCell, is_cell_visited, cids);
            }
        }
    }

    return cids.size();
}

size_t BaseComplexQuad::RegionGrowFace(const Face& start_face, std::vector<bool>& is_face_visited, std::vector<size_t>& fids)
{
    is_face_visited[start_face.id] = true;
    fids.push_back(start_face.id);
    for (auto edge_id : start_face.Eids) {
        const Edge& edge = mesh.E.at(edge_id);
        if (!IsOnBaseComplex(edge)) {
            const size_t nextFaceId = edge.N_Fids[0] == start_face.id ? edge.N_Fids[1] : edge.N_Fids[0];
            if (!is_face_visited[nextFaceId]) {
                const Face& nextFace = mesh.F.at(nextFaceId);
                RegionGrowFace(nextFace, is_face_visited, fids);
            }
        }
    }

    return fids.size();
}

size_t BaseComplexQuad::RegionGrowEdge(const Edge& start_edge, std::vector<bool>& is_edge_visited, std::vector<size_t>& eids)
{
    is_edge_visited[start_edge.id] = true;
    eids.push_back(start_edge.id);
    for (auto vertex_id : start_edge.Vids) {
        const Vertex& vertex = mesh.V.at(vertex_id);
        if (!IsOnBaseComplex(vertex)) {
            Edge* nextEdge = GetNextEdge(vertex, start_edge);
            //const size_t nextEdgeId = vertex.N_Eids[0] == start_edge.id ? vertex.N_Eids[1] : vertex.N_Eids[0];
            if (!is_edge_visited[nextEdge->id]) {
                //const Edge& nextEdge = mesh.E.at(nextEdgeId);
                RegionGrowEdge(*nextEdge, is_edge_visited, eids);
            }
        }
    }
    return eids.size();
}


void BaseComplexQuad::WriteSingularV_VTK(const char* filename) const
{
    std::ofstream ofs(filename);

    ofs << "# vtk DataFile Version 2.0" << std::endl << "Singular Vertex" << std::endl;
    ofs << "ASCII" << std::endl;
    ofs << "DATASET POLYDATA" << std::endl;
    ofs << "POINTS " << mesh.V.size() << " double" << std::endl;
    for (size_t i = 0; i < mesh.V.size(); i++) {
        const Vertex& v = mesh.V[i];
        ofs << v.x << " " << v.y << " " << v.z << std::endl;
    }
    ofs << "VERTICES " << SingularityI.V.size() << " " << SingularityI.V.size() * 2 << std::endl;
    for (int i = 0; i < SingularityI.V.size(); i++)
        ofs << "1 " << SingularityI.V[i].id_mesh << std::endl;
}

void BaseComplexQuad::WriteSingularE_VTK(const char *filename) const
{
    std::ofstream ofs(filename);

    ofs << "# vtk DataFile Version 2.0" << std::endl << "Singularities Color" << std::endl;
    ofs << "ASCII" << std::endl;
    ofs << "DATASET POLYDATA" << std::endl;

    ofs << "POINTS " << mesh.V.size() << " double" << std::endl;
    for (size_t i = 0; i < mesh.V.size(); i++) {
        const Vertex& v = mesh.V[i];
        ofs << v.x << " " << v.y << " " << v.z << std::endl;
    }

    int line_num = 0;
    for (int i = 0; i < SingularityI.E.size(); i++)
        line_num += SingularityI.E[i].vs_link.size();

    ofs << "LINES " << SingularityI.E.size() << " " << line_num + SingularityI.E.size() << std::endl;
    for (int i = 0; i < SingularityI.E.size(); i++) {
        ofs << SingularityI.E[i].vs_link.size() << " ";
        for (int j = 0; j < SingularityI.E[i].vs_link.size(); j++)
            ofs << SingularityI.E[i].vs_link[j] << " ";
        ofs << std::endl;
    }

    ofs << "CELL_DATA " << SingularityI.E.size() << std::endl;
    ofs << "SCALARS SingularityEdge int" << std::endl;
    ofs << "LOOKUP_TABLE SingularEdge" << std::endl;
    for (int i = 0; i < SingularityI.E.size(); i++) ofs << i << std::endl;

//    ofs << "CELL_DATA " << SingularityI.E.size() << std::endl;
//    ofs << "SCALARS valence int" << std::endl;
//    ofs << "LOOKUP_TABLE valence" << std::endl;
//    for (int i = 0; i < SingularityI.E.size(); i++) ofs << SingularityI.E[i].es_link.front() << std::endl;
}

void BaseComplexQuad::WriteSingularities_VTK(const char *filename) const
{
    std::ofstream ofs(filename);

    ofs << "# vtk DataFile Version 2.0" << std::endl << "Singularities Color" << std::endl;
    ofs << "ASCII" << std::endl;
    ofs << "DATASET POLYDATA" << std::endl;

    ofs << "POINTS " << mesh.V.size() << " double" << std::endl;
    for (size_t i = 0; i < mesh.V.size(); i++) {
        const Vertex& v = mesh.V[i];
        ofs << v.x << " " << v.y << " " << v.z << std::endl;
    }

    ofs << "VERTICES " << SingularityI.V.size() << " " << SingularityI.V.size() * 2 << std::endl;
    for (int i = 0; i < SingularityI.V.size(); i++)
        ofs << "1 " << SingularityI.V[i].id_mesh << std::endl;

    int line_num = 0;
    for (int i = 0; i < SingularityI.E.size(); i++)
        line_num += SingularityI.E[i].vs_link.size();

    ofs << "LINES " << SingularityI.E.size() << " " << line_num + SingularityI.E.size() << std::endl;
    for (int i = 0; i < SingularityI.E.size(); i++) {
        ofs << SingularityI.E[i].vs_link.size() << " ";
        for (int j = 0; j < SingularityI.E[i].vs_link.size(); j++)
            ofs << SingularityI.E[i].vs_link[j] << " ";
        ofs << std::endl;
    }

    ofs << "CELL_DATA " << SingularityI.V.size() + SingularityI.E.size() << std::endl;
    ofs << "SCALARS SingularityEdge int" << std::endl;
    ofs << "LOOKUP_TABLE SingularEdge" << std::endl;
    for (int i = 0; i < SingularityI.V.size(); i++) ofs << 0 << std::endl;
    for (int i = 0; i < SingularityI.E.size(); i++) ofs << i << std::endl;
}

void BaseComplexQuad::WriteBaseComplex_ColorVerticesVTK(const char *filename) const
{
    const std::vector<Vertex>& V = mesh.V;
    const std::vector<Edge>& E = mesh.E;
    const std::vector<Face>& F = F;

    std::ofstream ofs(filename);
    ofs << "# vtk DataFile Version 2.0" << std::endl
        << filename << std::endl
        << "ASCII" << std::endl << std::endl
        << "DATASET POLYDATA" << std::endl;
    ofs << "POINTS " << V.size() << " double" << std::endl;

    for (size_t i = 0; i < V.size(); i++)
        ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << std::endl;

    ofs << "VERTICES " << Vids.size() << " " << 2 * Vids.size() << std::endl;
    for (size_t i = 0; i < Vids.size(); i++)
        ofs << "1 " << V.at(Vids.at(i)).id << std::endl;

    ofs << "CELL_DATA " << Vids.size() << std::endl
        << "SCALARS " << " vertexid" << " int 1\n"
        << "LOOKUP_TABLE default" << std::endl;
    for (size_t i = 0; i < Vids.size(); i++)
        ofs << i << "\n";

    ofs << "SCALARS " << " singularity" << " int 1\n"
        << "LOOKUP_TABLE default\n";
    for (size_t i = 0; i < Vids.size(); i++)
        ofs << (V.at(Vids.at(i)).isSingularity ? 1 : 0) << std::endl;

}

void BaseComplexQuad::WriteBaseComplex_ColorEdgesVTK(const char *filename) const
{
    const std::vector<Vertex>& V = mesh.V;
    const std::vector<Edge>& E = mesh.E;
    const std::vector<Face>& F = F;

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

    size_t eid_count = 0;
    for (auto& componentEdge : componentE)
        eid_count += 1 + componentEdge.vids_link.size();

    ofs << "LINES " << componentE.size() << " " << eid_count << std::endl;
    for (auto& componentEdge : componentE) {
        ofs << componentEdge.vids_link.size();
        for (auto vid : componentEdge.vids_link)
            ofs << " " << vid;
        ofs << "\n";
    }

    ofs << "CELL_DATA " << /*Vids.size() + */componentE.size() << std::endl
        << "SCALARS " << " edgeid" << " int 1" << std::endl
        << "LOOKUP_TABLE default" << std::endl;
//    for (size_t i = 0; i < Vids.size(); i++)
//        ofs << (V.at(Vids.at(i)).isSingularity ? componentE.size() : 0) << std::endl;
    for (auto& componentEdge : componentE)
        ofs << componentEdge.id << std::endl;


    ofs << "SCALARS " << " color" << " int 1" << std::endl
        << "LOOKUP_TABLE default" << std::endl;
//    for (size_t i = 0; i < Vids.size(); i++)
//        ofs << (V.at(Vids.at(i)).isSingularity ? 8 : 0) << std::endl;
    for (auto& componentEdge : componentE)
        ofs << componentEdge.id % 8 << std::endl;

    ofs << "SCALARS " << " singularity" << " int 1" << std::endl
        << "LOOKUP_TABLE default" << std::endl;
//    for (size_t i = 0; i < Vids.size(); i++)
//        ofs << (V.at(Vids.at(i)).isSingularity ? 1 : 0) << std::endl;
    for (auto& componentEdge : componentE)
        ofs << (E.at(componentEdge.eids_link.front()).isSingularity ? 1 : 0) << std::endl;

}

void BaseComplexQuad::WriteBaseComplex_ColorFacesVTK(const char *filename) const
{
    const std::vector<Vertex>& V = mesh.V;
    const std::vector<Edge>& E = mesh.E;
    const std::vector<Face>& F = F;

    std::ofstream ofs(filename);
    ofs << "# vtk DataFile Version 2.0" << std::endl
        << filename << std::endl
        << "ASCII" << std::endl << std::endl
        << "DATASET POLYDATA" << std::endl;
    ofs << "POINTS " << V.size() << " double" << std::endl;

    for (size_t i = 0; i < V.size(); i++)
        ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << std::endl;

    ofs << "POLYGONS " << Fids.size() << " " << 5 * Fids.size() << std::endl;
    for (const auto& facePatch : componentF) {
        for (const auto& faceid : facePatch.fids_patch) {
            const Face& face = mesh.F.at(faceid);
            ofs << face.Vids.size();
            for (size_t j = 0; j < face.Vids.size(); j++)
                ofs << " " << face.Vids.at(j);
            ofs << "\n";
        }
    }

    ofs << "CELL_DATA " << Fids.size() << std::endl
        << "SCALARS " << " faceid" << " int 1" << std::endl
        << "LOOKUP_TABLE default" << std::endl;

    for (const auto& facePatch : componentF)
        for (const auto& faceid : facePatch.fids_patch)
            ofs << facePatch.id << std::endl;

    ofs << "SCALARS " << " color" << " int 1" << std::endl
        << "LOOKUP_TABLE default" << std::endl;

    for (const auto& facePatch : componentF)
        for (const auto& faceid : facePatch.fids_patch)
            ofs << facePatch.color % 4 << std::endl;

}

void BaseComplexQuad::WriteBaseComplexComponentsVTK(const char *filename) const
{
    const std::vector<Vertex>& V = mesh.V;
    const std::vector<Edge>& E = mesh.E;
    const std::vector<Face>& F = F;

    std::ofstream ofs(filename);
    ofs << "# vtk DataFile Version 2.0" << std::endl
        << filename << std::endl
        << "ASCII" << std::endl << std::endl
        << "DATASET UNSTRUCTURED_GRID" << std::endl;

    ofs << "POINTS " << V.size() << " double" << std::endl;
    for (size_t i = 0; i < V.size(); i++)
        ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << std::endl;

    ofs << "CELLS " << Fids.size() << " " << 5 * Fids.size() << std::endl;
    for (const auto& facePatch : componentF) {
        for (const auto& faceid : facePatch.fids_patch) {
            const Face& face = mesh.F.at(faceid);
            ofs << face.Vids.size();
            for (size_t j = 0; j < face.Vids.size(); j++)
                ofs << " " << face.Vids.at(j);
            ofs << "\n";
        }
    }

    ofs << "CELL_TYPES " << Fids.size() << std::endl;
    for (const auto id : Fids)
        ofs << "9\n";

    ofs << "CELL_DATA " << Fids.size() << std::endl
        << "SCALARS " << " componentid" << " int 1" << std::endl
        << "LOOKUP_TABLE default" << std::endl;
    for (const auto& facePatch : componentF)
        for (const auto& faceid : facePatch.fids_patch)
            ofs << facePatch.id << std::endl;
    ofs << "SCALARS " << " color" << " int 1" << std::endl
        << "LOOKUP_TABLE default" << std::endl;
    for (const auto& facePatch : componentF)
        for (const auto& faceid : facePatch.fids_patch)
            ofs << facePatch.color % 8 << std::endl;
}

void BaseComplexQuad::WriteBaseComplexAllComponentsVTK(const char *filename_prefix) const {
    for (size_t id = 0; id < componentC.size(); ++id)
        WriteBaseComplexComponentsVTK(filename_prefix, id);
}

void BaseComplexQuad::WriteBaseComplexComponentsVTK(const char *filename_prefix, const size_t id) const
{
    const std::vector<Vertex>& V = mesh.V;
    const std::vector<Edge>& E = mesh.E;
    const std::vector<Face>& F = F;

    std::string filename = std::string(filename_prefix) + std::to_string(id) + ".vtk";
    std::ofstream ofs(filename);
    ofs << "# vtk DataFile Version 2.0" << std::endl
        << filename << std::endl
        << "ASCII" << std::endl << std::endl
        << "DATASET UNSTRUCTURED_GRID" << std::endl;

    ofs << "POINTS " << V.size() << " double" << std::endl;
    for (size_t i = 0; i < V.size(); i++)
        ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << std::endl;

    const auto& cellPatch = componentC.at(id);
    ofs << "CELLS " << cellPatch.cids_patch.size() << " " << 9 * cellPatch.cids_patch.size() << std::endl;

    for (const auto& cellid : cellPatch.cids_patch) {
        const Cell& cell = mesh.C.at(cellid);
        ofs << cell.Vids.size();
        for (size_t j = 0; j < cell.Vids.size(); j++)
            ofs << " " << cell.Vids.at(j);
        ofs << "\n";
    }

    ofs << "CELL_TYPES " << cellPatch.cids_patch.size() << std::endl;
    for (const auto cellid : cellPatch.cids_patch)
        ofs << "12\n";

    ofs << "CELL_DATA " << cellPatch.cids_patch.size() << std::endl
        << "SCALARS " << " componentid" << " int 1" << std::endl
        << "LOOKUP_TABLE default" << std::endl;
        for (const auto& cellid : cellPatch.cids_patch)
            ofs << id << std::endl;
    ofs << "SCALARS " << " color" << " int 1" << std::endl
        << "LOOKUP_TABLE default" << std::endl;
        for (const auto& cellid : cellPatch.cids_patch)
            ofs << cellPatch.color % 8 << std::endl;
}

void BaseComplexQuad::WriteBaseComplexHexVTK(const char *filename) const
{
//    const std::vector<Vertex>& V = mesh.V;
//    const std::vector<Edge>& E = mesh.E;
//    const std::vector<Face>& F = F;

    std::ofstream ofs(filename);
    ofs << "# vtk DataFile Version 2.0" << std::endl
        << filename << std::endl
        << "ASCII" << std::endl << std::endl
        << "DATASET UNSTRUCTURED_GRID" << std::endl;

    ofs << "POINTS " << V.size() << " double" << std::endl;
    for (auto& v : V)
        ofs << v.x << " " << v.y << " " << v.z << std::endl;

    ofs << "CELLS " << componentC.size() << " " << 9 * componentC.size() << std::endl;
    for (const auto& component : componentC) {
        ofs << component.Vids.size();
        for (auto vid : component.Vids)
            ofs << " " << vid;
        ofs << "\n";
    }

    ofs << "CELL_TYPES " << componentC.size() << std::endl;
    for (const auto& id : componentC)
        ofs << "12\n";
    ofs << "CELL_DATA " << componentC.size() << std::endl
        << "SCALARS " << " componentid" << " int 1" << std::endl
        << "LOOKUP_TABLE default" << std::endl;
    for (const auto& component : componentC)
            ofs << component.id << std::endl;
    ofs << "SCALARS " << " color" << " int 1" << std::endl
        << "LOOKUP_TABLE default" << std::endl;
    for (const auto& component : componentC)
            ofs << component.color % 8 << std::endl;
}

void BaseComplexQuad::WriteBaseComplexComponentsWithoutSingularitiesVTK(const char *filename) const
{
    const std::vector<Vertex>& V = mesh.V;
    const std::vector<Edge>& E = mesh.E;
    const std::vector<Face>& F = F;

    std::ofstream ofs(filename);
    ofs << "# vtk DataFile Version 2.0" << std::endl
        << filename << std::endl
        << "ASCII" << std::endl << std::endl
        << "DATASET UNSTRUCTURED_GRID" << std::endl;

    ofs << "POINTS " << V.size() << " double" << std::endl;
    for (size_t i = 0; i < V.size(); i++)
        ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << std::endl;

    int count = 0;
    std::vector<size_t> hasSingularities(componentC.size(), false);
    int sid = 0;
    for (const auto& cellPatch : componentC) {
        size_t num = 0;
        for (const auto& cellid : cellPatch.cids_patch) {
            const Cell& cell = mesh.C.at(cellid);
            bool foundSingularity = false;
            for (auto eid : cell.Eids) {
                if (mesh.E.at(eid).isSingularity) {
                    foundSingularity = true;
                }
            }
            if (!foundSingularity) ++num;
            else {
                num = 0;
                break;
            }
        }
        count += num;
        hasSingularities[sid++] = num == 0;
    }

    sid = 0;
    ofs << "CELLS " << count << " " << 9 * count << std::endl;
    for (const auto& cellPatch : componentC) {
        if (hasSingularities[sid++]) continue;
        for (const auto& cellid : cellPatch.cids_patch) {
            const Cell& cell = mesh.C.at(cellid);
            ofs << cell.Vids.size();
            for (size_t j = 0; j < cell.Vids.size(); j++)
                ofs << " " << cell.Vids.at(j);
            ofs << "\n";
        }
    }

    ofs << "CELL_TYPES " << count << std::endl;
    for (size_t i = 0; i < count; ++i)
        ofs << "12\n";

    sid = 0;
    ofs << "CELL_DATA " << count << std::endl
        << "SCALARS " << " componentid" << " int 1" << std::endl
        << "LOOKUP_TABLE default" << std::endl;
    for (const auto& cellPatch : componentC)
        if (!hasSingularities[sid++])
        for (const auto& cellid : cellPatch.cids_patch)
                ofs << cellPatch.id << std::endl;
    sid = 0;
    ofs << "SCALARS " << " color" << " int 1" << std::endl
        << "LOOKUP_TABLE default" << std::endl;
    for (const auto& cellPatch : componentC)
        if (!hasSingularities[sid++])
        for (const auto& cellid : cellPatch.cids_patch)
                ofs << cellPatch.color % 8 << std::endl;
}

void BaseComplexQuad::WriteBaseComplexComponentsWithSingularitiesVTK(const char *filename) const
{
    const std::vector<Vertex>& V = mesh.V;
    const std::vector<Edge>& E = mesh.E;
    const std::vector<Face>& F = F;

    std::ofstream ofs(filename);
    ofs << "# vtk DataFile Version 2.0" << std::endl
        << filename << std::endl
        << "ASCII" << std::endl << std::endl
        << "DATASET UNSTRUCTURED_GRID" << std::endl;

    ofs << "POINTS " << V.size() << " double" << std::endl;
    for (size_t i = 0; i < V.size(); i++)
        ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << std::endl;

    int count = 0;
    std::vector<size_t> hasSingularities(componentC.size(), false);
    int sid = 0;
    for (const auto& cellPatch : componentC) {
        size_t num = 0;
        bool foundSingularity = false;
        for (const auto& cellid : cellPatch.cids_patch) {
            const Cell& cell = mesh.C.at(cellid);
            for (auto eid : cell.Eids)
                if (mesh.E.at(eid).isSingularity)
                    foundSingularity = true;
            if (foundSingularity) break;
        }
        if (foundSingularity) count += cellPatch.cids_patch.size();
        hasSingularities[sid++] = foundSingularity;
    }

    sid = 0;
    ofs << "CELLS " << count << " " << 9 * count << std::endl;
    for (const auto& cellPatch : componentC) {
        if (!hasSingularities[sid++]) continue;
        for (const auto& cellid : cellPatch.cids_patch) {
            const Cell& cell = mesh.C.at(cellid);
            ofs << cell.Vids.size();
            for (size_t j = 0; j < cell.Vids.size(); j++)
                ofs << " " << cell.Vids.at(j);
            ofs << "\n";
        }
    }

    ofs << "CELL_TYPES " << count << std::endl;
    for (size_t i = 0; i < count; ++i)
        ofs << "12\n";

    sid = 0;
    ofs << "CELL_DATA " << count << std::endl
        << "SCALARS " << " componentid" << " int 1" << std::endl
        << "LOOKUP_TABLE default" << std::endl;
    for (const auto& cellPatch : componentC)
        if (hasSingularities[sid++])
        for (const auto& cellid : cellPatch.cids_patch)
                ofs << cellPatch.id << std::endl;
    sid = 0;
    ofs << "SCALARS " << " color" << " int 1" << std::endl
        << "LOOKUP_TABLE default" << std::endl;
    for (const auto& cellPatch : componentC)
        if (hasSingularities[sid++])
        for (const auto& cellid : cellPatch.cids_patch)
                ofs << cellPatch.color % 8 << std::endl;
}

void BaseComplexQuad::WriteBaseComplexSeparatedFacePatchesVTK(const char *filename) const
{
    const std::vector<Vertex>& V = mesh.V;
    const std::vector<Edge>& E = mesh.E;
    const std::vector<Face>& F = F;

    std::ofstream ofs(filename);
    ofs << "# vtk DataFile Version 2.0" << std::endl
        << filename << std::endl
        << "ASCII" << std::endl << std::endl
        << "DATASET POLYDATA" << std::endl;
    ofs << "POINTS " << V.size() << " double" << std::endl;

    for (size_t i = 0; i < V.size(); i++)
        ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << std::endl;

    ofs << "VERTICES " << Vids.size() << " " << 2 * Vids.size() << std::endl;
    for (size_t i = 0; i < Vids.size(); i++)
        ofs << "1 " << V.at(Vids.at(i)).id << std::endl;

    ofs << "LINES " << Eids.size() << " " << 3 * Eids.size() << std::endl;
    for (size_t i = 0; i < Eids.size(); i++)
        ofs << "2 " << E.at(Eids.at(i)).Vids[0] << " " << E.at(Eids.at(i)).Vids[1] << std::endl;

    ofs << "POLYGONS " << Fids.size() << " " << 5 * Fids.size() << std::endl;
    for (const auto& patch : separatedFacePatches) {
        for (const auto& faceid : patch) {
            const Face& face = mesh.F.at(faceid);
            ofs << face.Vids.size();
            for (size_t j = 0; j < face.Vids.size(); j++)
                ofs << " " << face.Vids.at(j);
            ofs << "\n";
        }
    }

    ofs << "CELL_DATA " << Vids.size() + Eids.size() + Fids.size() << std::endl
        << "SCALARS " << " patch" << " int 1" << std::endl
        << "LOOKUP_TABLE default" << std::endl;
    for (size_t i = 0; i < Vids.size(); i++)
        ofs << (V.at(Vids.at(i)).isSingularity ? separatedFacePatches.size() : 0) << std::endl;
    for (size_t i = 0; i < Eids.size(); i++)
        ofs << (E.at(Eids.at(i)).isSingularity ? separatedFacePatches.size() : 0) << std::endl;
    int id = 0;
    for (const auto& patch : separatedFacePatches) {
        for (const auto& faceid : patch)
            ofs << id << "\n";
        ++id;
    }


    ofs << "SCALARS " << " color" << " int 1" << std::endl
        << "LOOKUP_TABLE default" << std::endl;
    for (size_t i = 0; i < Vids.size(); i++)
        ofs << (V.at(Vids.at(i)).isSingularity ? 8 : 0) << std::endl;
    for (size_t i = 0; i < Eids.size(); i++)
        ofs << (E.at(Eids.at(i)).isSingularity ? 8 : 0) << std::endl;
    id = 0;
    for (const auto& patch : separatedFacePatches) {
        for (const auto& faceid : patch)
            ofs << id % 8 << "\n";
        ++id;
    }

}

void BaseComplexQuad::WriteBaseComplex_VTK(const char *filename) const
{
    const std::vector<Vertex>& V = mesh.V;
    const std::vector<Edge>& E = mesh.E;
    const std::vector<Face>& F = mesh.F;

    std::ofstream ofs(filename);
    ofs << "# vtk DataFile Version 2.0\n"
        << filename << "\n"
        << "ASCII\n\n"
        << "DATASET POLYDATA\n";
    ofs << "POINTS " << V.size() << " double\n";

    for (size_t i = 0; i < V.size(); i++)
        ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << "\n";

    ofs << "VERTICES " << Vids.size() << " " << 2 * Vids.size() << "\n";
    for (size_t i = 0; i < Vids.size(); i++)
        ofs << "1 " << V.at(Vids.at(i)).id << "\n";

    ofs << "LINES " << Eids.size() << " " << 3 * Eids.size() << "\n";
    for (size_t i = 0; i < Eids.size(); i++)
        ofs << "2 " << E.at(Eids.at(i)).Vids[0] << " " << E.at(Eids.at(i)).Vids[1] << "\n";

    ofs << "POLYGONS " << Fids.size() << " " << 5 * Fids.size() << "\n";
    for (size_t i = 0; i < Fids.size(); i++) {
        const Face& face = F.at(Fids.at(i));
        ofs << face.Vids.size();
        for (size_t j = 0; j < face.Vids.size(); j++)
            ofs << " " << face.Vids.at(j);
        ofs << "\n";
    }

    ofs << "CELL_DATA " << Vids.size() + Eids.size() + Fids.size() << "\n"
        << "SCALARS " << " Singularity" << " int 1\n"
        << "LOOKUP_TABLE default" << std::endl;
    for (size_t i = 0; i < Vids.size(); i++)
        ofs << (V.at(Vids.at(i)).isSingularity ? 1 : 0) << "\n";
    for (size_t i = 0; i < Eids.size(); i++)
        ofs << (E.at(Eids.at(i)).isSingularity ? 1 : 0) << "\n";
    for (size_t i = 0; i < Fids.size(); i++)
        ofs << 0 << "\n";
}

void BaseComplexQuad::WriteComponentEdge_NeighborComponentFaces_VTK(const char *filename) const
{
    const std::vector<Vertex>& V = mesh.V;
    const std::vector<Edge>& E = mesh.E;
    std::ofstream ofs(filename);
    ofs << "# vtk DataFile Version 2.0" << std::endl
        << filename << std::endl
        << "ASCII" << std::endl << std::endl
        << "DATASET POLYDATA" << std::endl;
    ofs << "POINTS " << V.size() << " double" << std::endl;

    for (size_t i = 0; i < V.size(); i++)
        ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << std::endl;

    ofs << "VERTICES " << Vids.size() << " " << 2 * Vids.size() << std::endl;
    for (size_t i = 0; i < Vids.size(); i++)
        ofs << "1 " << V.at(Vids.at(i)).id << std::endl;

    int edge_num = 0;
    for (const auto & componentEdge : componentE)
        edge_num += componentEdge.vids_link.size() + 1;

    int face_num = 0;
    for (const auto & componentEdge : componentE)
        for (auto & componentFaceId : componentEdge.N_Fids)
            face_num += componentF.at(componentFaceId).fids_patch.size();

//    ofs << "LINES " << componentE.size() << " " << edge_num << std::endl;
//    for (const auto & componentEdge : componentE) {
//        ofs << componentEdge.vids_link.size();
//        for (auto vid : componentEdge.vids_link)
//            ofs << " " << vid;
//        ofs << "\n";
//    }
    ofs << "LINES " << Eids.size() << " " << 3 * Eids.size() << std::endl;
    for (size_t i = 0; i < Eids.size(); i++)
        ofs << "2 " << E.at(Eids.at(i)).Vids[0] << " " << E.at(Eids.at(i)).Vids[1] << std::endl;

    ofs << "POLYGONS " << face_num << " " << 5 * face_num << std::endl;
    for (const auto & componentEdge : componentE) {
        for (auto & componentFaceId : componentEdge.N_Fids) {
            for (const auto& faceid : componentF.at(componentFaceId).fids_patch) {
                const Face& face = mesh.F.at(faceid);
                ofs << face.Vids.size();
                for (size_t j = 0; j < face.Vids.size(); j++)
                    ofs << " " << face.Vids.at(j);
                ofs << "\n";
            }
        }
    }


    ofs << "CELL_DATA " << Vids.size() + Eids.size() + face_num << std::endl
        << "SCALARS " << "NeighborFaces" << " int 1" << std::endl
        << "LOOKUP_TABLE default" << std::endl;
    for (size_t i = 0; i < Vids.size(); i++)
        ofs << (V.at(Vids.at(i)).isSingularity ? separatedFacePatches.size() : 0) << "\n";
//    for (size_t i = 0; i < componentE.size(); i++)
//        ofs << i << "\n";
    for (size_t i = 0; i < Eids.size(); i++)
        ofs << (E.at(Eids.at(i)).componentEid) << std::endl;
    for (size_t i = 0; i < componentE.size(); i++) {
        const auto & componentEdge = componentE[i];
        for (auto & componentFaceId : componentEdge.N_Fids) {
            for (const auto& faceid : componentF.at(componentFaceId).fids_patch) {
                ofs << i << "\n";
            }
        }
    }
}

void BaseComplexQuad::WriteComponentEdge_NeighborComponentCells_VTK(const char *filename) const
{

}

void BaseComplexQuad::WriteComponentFace_NeighborComponentCells_VTK(const char *filename) const
{

}

void BaseComplexQuad::WriteSingularEdge_NeighborSeparatedFacePatches_VTK(const char *filename) const
{

}

void BaseComplexQuad::WriteAllSingularEdge_NeighborSeparatedComponentFacePatches_VTK(const char *filename_prefix) const
{
    int id = 0;
    for (const auto& e : SingularityI.E)
        WriteSingularEdge_NeighborSeparatedComponentFacePatches_VTK((std::string(filename_prefix) + std::to_string(id) + ".vtk").c_str(), id++);
}
void BaseComplexQuad::WriteSingularEdge_NeighborSeparatedComponentFacePatches_VTK(const char *filename) const
{
    const std::vector<Vertex>& V = mesh.V;
    const std::vector<Edge>& E = mesh.E;
    std::ofstream ofs(filename);
    ofs << "# vtk DataFile Version 2.0" << std::endl
        << filename << std::endl
        << "ASCII" << std::endl << std::endl
        << "DATASET POLYDATA" << std::endl;
    ofs << "POINTS " << V.size() << " double" << std::endl;

    for (size_t i = 0; i < V.size(); i++)
        ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << std::endl;

    ofs << "VERTICES " << Vids.size() << " " << 2 * Vids.size() << std::endl;
    for (size_t i = 0; i < Vids.size(); i++)
        ofs << "1 " << V.at(Vids.at(i)).id << std::endl;

    int edge_num = 0;
    for (const auto & singularityEdge : SingularityI.E)
        edge_num += singularityEdge.vs_link.size() + 1;

    int face_num = 0;
//    for (const auto & singularityEdge : SingularityI.E)
//        for (auto & neighborComponentFidsGroup : singularityEdge.neighborComponentFidsGroups)
//            for (auto & componentFid : neighborComponentFidsGroup)
//                face_num += componentF.at(componentFid).fids_patch.size();
    for (const auto & singularityEdge : SingularityI.E)
        for (auto & separatedFacePatchId : singularityEdge.separatedFacePatchIds)
//            for (auto & componentFid : separatedFacePatches.at(separatedFacePatchId))
                face_num += separatedFacePatches.at(separatedFacePatchId).size();

    ofs << "LINES " << SingularityI.E.size() << " " << edge_num << std::endl;
    for (const auto & singularityEdge : SingularityI.E) {
        ofs << singularityEdge.vs_link.size();
        for (const auto vid : singularityEdge.vs_link)
            ofs << " " << vid;
        ofs << "\n";
    }

    ofs << "POLYGONS " << face_num << " " << 5 * face_num << std::endl;
    for (const auto & singularityEdge : SingularityI.E)
        for (auto & separatedFacePatchId : singularityEdge.separatedFacePatchIds)
            for (auto & face_id : separatedFacePatches.at(separatedFacePatchId))
                //for (auto & face_id : componentF.at(componentFid).fids_patch)
                {
                    const Face& face = mesh.F.at(face_id);
                    ofs << face.Vids.size();
                    for (size_t j = 0; j < face.Vids.size(); j++)
                        ofs << " " << face.Vids.at(j);
                    ofs << "\n";
                }

    ofs << "CELL_DATA " << Vids.size() + SingularityI.E.size() + face_num << std::endl
        << "SCALARS " << "SingularEdge_SeparatedFacePatches" << " int 1" << std::endl
        << "LOOKUP_TABLE default" << std::endl;
    for (size_t i = 0; i < Vids.size(); i++)
        ofs << (V.at(Vids.at(i)).isSingularity ? SingularityI.E.size() : 0) << "\n";
    for (size_t i = 0; i < SingularityI.E.size(); i++)
        ofs << i << "\n";
    int count = 0;
    for (const auto & singularityEdge : SingularityI.E) {
        for (auto & separatedFacePatchId : singularityEdge.separatedFacePatchIds)
            for (auto & face_id : separatedFacePatches.at(separatedFacePatchId))
                    ofs << count << "\n";
        ++count;
    }

//    ofs << "SCALARS " << "SingularEdge_SeparatedComponentFacePatches" << " int 1" << std::endl
//        << "LOOKUP_TABLE default\n";
//    for (size_t i = 0; i < Vids.size(); i++)
//        ofs << (V.at(Vids.at(i)).isSingularity ? componentF.size() : 0) << "\n";
//    for (size_t i = 0; i < SingularityI.E.size(); i++)
//        ofs << SingularityI.E.size() << "\n";
//    for (const auto & singularityEdge : SingularityI.E) {
//        for (auto & separatedFacePatchId : singularityEdge.separatedFacePatchIds)
//            for (auto & face_id : separatedFacePatches.at(separatedFacePatchId))
//                    ofs << mesh.F.at(face_id).componentFid << std::endl;
//    }
}

void BaseComplexQuad::WriteSingularEdge_NeighborSeparatedComponentFacePatches_VTK(const char *filename, const int singularEdgeId) const
{
    const std::vector<Vertex>& V = mesh.V;
    const std::vector<Edge>& E = mesh.E;
    std::ofstream ofs(filename);
    ofs << "# vtk DataFile Version 2.0" << std::endl
        << filename << std::endl
        << "ASCII" << std::endl << std::endl
        << "DATASET POLYDATA" << std::endl;
    ofs << "POINTS " << V.size() << " double" << std::endl;

    for (size_t i = 0; i < V.size(); i++)
        ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << std::endl;

    auto & singularityEdge = SingularityI.E.at(singularEdgeId);

    int edge_num = singularityEdge.vs_link.size() + 1;
    int face_num = 0;
    for (auto & separatedFacePatchId : singularityEdge.separatedFacePatchIds)
            face_num += separatedFacePatches.at(separatedFacePatchId).size();

    ofs << "LINES " << 1 << " " << edge_num << std::endl;
    ofs << singularityEdge.vs_link.size();
    for (const auto vid : singularityEdge.vs_link)
        ofs << " " << vid;
    ofs << "\n";

    ofs << "POLYGONS " << face_num << " " << 5 * face_num << std::endl;
    for (auto & separatedFacePatchId : singularityEdge.separatedFacePatchIds)
        for (auto & face_id : separatedFacePatches.at(separatedFacePatchId)) {
            const Face& face = mesh.F.at(face_id);
            ofs << face.Vids.size();
            for (size_t j = 0; j < face.Vids.size(); j++)
                ofs << " " << face.Vids.at(j);
            ofs << "\n";
        }

    ofs << "CELL_DATA " << 1 + face_num << std::endl
        << "SCALARS " << "SingularEdge_SeparatedFacePatches" << " int 1" << std::endl
        << "LOOKUP_TABLE default" << std::endl;
    ofs << singularEdgeId << "\n";
    int count = 0;
    for (auto & separatedFacePatchId : singularityEdge.separatedFacePatchIds)
        for (auto & face_id : separatedFacePatches.at(separatedFacePatchId))
            ofs << count << "\n";
    ++count;
}

/*
                |
                |
     edge       | vertex
    ------------|----------- res
                |
                |
                |
    */
Edge* BaseComplexQuad::GetNextContinuousEdge(const Vertex& vertex, const Edge& edge) const {
    for (auto eid : vertex.N_Eids) {
        auto& e = mesh.E.at(eid);
        if (e.id != edge.id && !IsTwoEdgesOnOneFace(mesh, edge, e, vertex))
            return &e;
    }
    return NULL;
}
