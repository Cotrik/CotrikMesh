/*
 * Mesh.cpp
 *
 *  Created on: Nov 6, 2016
 *      Author: cotrik
 */

#include "Mesh.h"
#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "FeatureLine.h"
//#include "MeshQuality.h"
#include "glm/gtx/intersect.hpp"
#include <algorithm>
#include <map>
#include <iostream>

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkDataSet.h>
#include <vtkMeshQuality.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkTriangleFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkProperty.h>
#include <vtkScalarBarActor.h>
#include <vtkLookupTable.h>
#include <vtkColorTransferFunction.h>
#include <vtkThreshold.h>
#include <vtkExtractEdges.h>

//#include <eigen3/Eigen/Core>
//#include <eigen3/Eigen/Eigen>
//#include <eigen3/Eigen/Dense>
//#include <eigen3/Eigen/Sparse>
//#include <eigen3/Eigen/Cholesky>
//#include <eigen3/Eigen/LU>
//#include <eigen3/Eigen/SVD>
//#include <eigen3/Eigen/SparseCore>
//#include <eigen3/Eigen/SparseLU>
//#include <eigen3/Eigen/SparseQR>
//#include <eigen3/Eigen/SparseCholesky>
//using namespace Eigen;

const size_t MAXID = 0xffffffffffffffff;

Mesh::Mesh()
: m_cellType(HEXAHEDRA)
, avgEdgeLength(0.0)
, numOfSharpEdges(0)
{
    // TODO Auto-generated constructor stub

}

Mesh::~Mesh()
{
    // TODO Auto-generated destructor stub
    V.clear();
    E.clear();
    F.clear();
    C.clear();

    pointScalarFieldNames.clear();
    pointScalarFields.clear();
    cellScalarNameFields.clear();
    cellScalarFields.clear();
}

Mesh::Mesh(const Mesh& r)
: V(r.V)
, E(r.E)
, F(r.F)
, C(r.C)
, m_cellType(r.m_cellType)
, m_cellTypes(r.m_cellTypes)
, pointScalarFieldNames(r.pointScalarFieldNames)
, pointScalarFields(r.pointScalarFields)
, cellScalarNameFields(r.cellScalarNameFields)
, cellScalarFields(r.cellScalarFields)
, avgEdgeLength(r.avgEdgeLength)
, numOfSharpEdges(r.numOfSharpEdges)
{
    m_refIds.resize(V.size());
    for (size_t i = 0; i < V.size(); ++i) m_refIds[i] = i;
}

Mesh::Mesh(const std::vector<Vertex>& V, const std::vector<Cell>& C, ElementType m_cellType)
: V(V)
, C(C)
, m_cellType(m_cellType)
, avgEdgeLength(0.0)
, numOfSharpEdges(0)
{
    m_refIds.resize(V.size());
    for (size_t i = 0; i < V.size(); ++i) m_refIds[i] = i;
}

Mesh::Mesh(const std::vector<Vertex>& V, const std::vector<Face>& F, ElementType m_cellType)
: V(V)
, F(F)
, m_cellType(m_cellType)
, avgEdgeLength(0.0)
, numOfSharpEdges(0)
{
    m_refIds.resize(V.size());
    for (size_t i = 0; i < V.size(); ++i) m_refIds[i] = i;
}

Mesh::Mesh(const Mesh& r, const std::vector<size_t>& cellIds)
: m_cellType(r.m_cellType)
, avgEdgeLength(0.0)
, numOfSharpEdges(0)
{
    V.resize(r.V.size());
    // Read V
    for (vtkIdType i = 0; i < V.size(); i++) {
        Vertex& v = V.at(i);
        const Vertex& rv = r.V.at(i);
        v.x = rv.x;
        v.y = rv.y;
        v.z = rv.z;
        v.id = i;
    }
    // Read C
    C.resize(cellIds.size());
    for (vtkIdType i = 0; i < cellIds.size(); i++)
        C.at(i).Vids = r.C.at(cellIds.at(i)).Vids;

    if (m_cellType == TRIANGLE || r.m_cellType == QUAD) {
        F.resize(C.size());
        for (vtkIdType i = 0; i < cellIds.size(); i++)
            F[i].Vids = C[i].Vids;
    }

    m_refIds.resize(V.size());
    for (size_t i = 0; i < V.size(); ++i) m_refIds[i] = i;
}

void Mesh::BuildE() {
    std::vector<Edge> Es;
    if (m_cellType == HEXAHEDRA) Es.reserve(C.size() * 12);
    else if (m_cellType == TETRAHEDRA) Es.reserve(C.size() * 6);
    else if (m_cellType == QUAD) Es.reserve(C.size() * 4);
    else if (m_cellType == TRIANGLE) Es.reserve(C.size() * 3);
	else if (m_cellType == POLYGON) Es.reserve(C.size() * 4);
	else if (m_cellType == POLYHEDRA) Es.reserve(C.size() * 12);

    for (auto& c : C) {
        Edge e(2);
        if (c.cellType == VTK_HEXAHEDRON) for (size_t j = 0; j < 12; j++) { // there are 12 edges in hex
            e.Vids[0] = c.Vids[HexEdge[j][0]];
            e.Vids[1] = c.Vids[HexEdge[j][1]];
            Es.push_back(e);
        } else if (c.cellType == VTK_QUAD) for (size_t j = 0; j < 4; j++) {// there are 3 edges in tri 
            e.Vids[0] = c.Vids[QuadEdge[j][0]];
            e.Vids[1] = c.Vids[QuadEdge[j][1]];
            Es.push_back(e);
        } else if (c.cellType == VTK_TRIANGLE) for (size_t j = 0; j < 3; j++) { // there are 3 edges in tri
            e.Vids[0] = c.Vids[TriEdge[j][0]];
            e.Vids[1] = c.Vids[TriEdge[j][1]];
            Es.push_back(e);
        } else if (c.cellType == VTK_TETRA) for (size_t j = 0; j < 6; j++) { // there are 6 edges in tet
            e.Vids[0] = c.Vids[TetEdge[j][0]];
            e.Vids[1] = c.Vids[TetEdge[j][1]];
            Es.push_back(e);
        } else if (c.cellType == VTK_WEDGE) for (size_t j = 0; j < 9; j++) { // there are 9 edges in wedge
			e.Vids[0] = c.Vids[WedgeEdge[j][0]];
			e.Vids[1] = c.Vids[WedgeEdge[j][1]];
			Es.push_back(e);
		} else if (c.cellType == VTK_PENTAGONAL_PRISM) for (size_t j = 0; j < 15; j++) { // there are 15 edges in pentahedral
			e.Vids[0] = c.Vids[PentaEdge[j][0]];
			e.Vids[1] = c.Vids[PentaEdge[j][1]];
			Es.push_back(e);
		} else if (c.cellType == VTK_POLYGON) {
            if (c.Vids.size() == 3) {
                for (size_t j = 0; j < 3; j++) {
                    e.Vids[0] = c.Vids[TriEdge[j][0]];
                    e.Vids[1] = c.Vids[TriEdge[j][1]];
                    Es.push_back(e);
                }
            } else if (c.Vids.size() == 4) {
                for (size_t j = 0; j < 4; j++) {
                    e.Vids[0] = c.Vids[QuadEdge[j][0]];
                    e.Vids[1] = c.Vids[QuadEdge[j][1]];
                    Es.push_back(e);
                }
			} else {
				for (size_t j = 0; j < c.Vids.size(); j++) {
					e.Vids[0] = c.Vids[j % c.Vids.size()];
					e.Vids[1] = c.Vids[(j + 1) % c.Vids.size()];
					Es.push_back(e);
				}
			}
        }
    }

    size_t E_N = 0;
    for (auto& e : Es) {
        bool havesame = false;
        const size_t id1 = e.Vids[0];
        const size_t id2 = e.Vids[1];

        for (size_t j = 0; j < V[id1].N_Vids.size(); j++)
            if (V[id1].N_Vids[j] == id2) {
                havesame = true;
                break;
            }
        if (!havesame) {
            e.id = E_N++;
            e.isBoundary = false;

            for (auto cid0 : V[e.Vids[0]].N_Cids)
                for (auto cid1 : V[e.Vids[1]].N_Cids)
                    if (cid0 == cid1) e.N_Cids.push_back(cid0);

            E.push_back(e);
            V[id1].N_Eids.push_back(e.id);
            V[id1].N_Vids.push_back(id2);

            V[id2].N_Eids.push_back(e.id);
            V[id2].N_Vids.push_back(id1);
        }
    }

	BuildE_C();

    for (auto& v : V) {
		std::set<size_t> N_Vids(v.N_Vids.begin(), v.N_Vids.end());
		v.N_Vids.clear();
		for (auto nvid : N_Vids)
			v.N_Vids.push_back(nvid);
    }
}

void Mesh::BuildParallelE() {
    if (m_cellType == HEXAHEDRA) {
        for (size_t i = 0; i < F.size(); i++)
            for (size_t j = 0; j < F[i].N_Ortho_4Eids.size(); j++)
                for (size_t k = 0; k < 4; k++)
                    for (size_t m = 0; m < 4; m++)
                        if (k != m)
                            E[F[i].N_Ortho_4Eids[j][k]].parallelEids.push_back(F[i].N_Ortho_4Eids[j][m]);

        for (size_t i = 0; i < E.size(); i++)
			Util::set_redundent_clearn(E[i].parallelEids);
    } else //if (m_cellType == QUAD) 
	{
        for (auto& e : E) {
            for (auto& fid : e.N_Fids)
				if (F.at(fid).Eids.size() == 4)
                for (auto& edgeid : F.at(fid).Eids) {
                    auto& edge = E.at(edgeid);
                    if (edge.Vids[0] != e.Vids[0] && edge.Vids[1] != e.Vids[1] && edge.Vids[0] != e.Vids[1] && edge.Vids[1] != e.Vids[0]) e.parallelEids.push_back(edgeid);
                }
        }
    }
}

void Mesh::BuildConsecutiveE() {
	if (m_cellType == HEXAHEDRA) {
		for (size_t i = 0; i < E.size(); i++) {
			std::vector<size_t> fes;
			for (size_t j = 0; j < E[i].N_Fids.size(); j++) {
				size_t fid = E[i].N_Fids[j];
				for (size_t k = 0; k < 4; k++)
					fes.push_back(F[fid].Eids[k]);
			}
			Util::set_redundent_clearn(fes);
			for (size_t j = 0; j < 2; j++) {
				int vid = E[i].Vids[j];
				std::vector<size_t> leftes;
				Util::set_exclusion(V[vid].N_Eids, fes, leftes);
				std::copy(leftes.begin(), leftes.end(), back_inserter(E[i].consecutiveEids));
			}
			Util::set_redundent_clearn(E[i].consecutiveEids);
		}
	} else //if (m_cellType == QUAD) 
	{
		for (auto& e : E) {
			std::set<size_t> fids;
			for (auto vid : e.Vids)
				if (!V[vid].isSingularity) fids.insert(V[vid].N_Fids.begin(), V[vid].N_Fids.end());
			for (auto fid : e.N_Fids)
				fids.erase(fid);

			std::set<size_t> consecutive_eids;
			for (auto fid : fids) {
				//if (F.at(fid).Eids.size() != 4) break;
				int count = 0;
				for (auto vid : F.at(fid).Vids)
					if (vid == e.Vids[0] || vid == e.Vids[1]) ++count;
				if (count > 1) continue;
				for (auto edgeid : F.at(fid).Eids) {
					auto& edge = E.at(edgeid);
					//if (edge.N_Fids[0] != e.N_Fids[0] && edge.N_Fids[1] != e.N_Fids[1] && edge.N_Fids[0] != e.Vids[1] && edge.N_Fids[1] != e.N_Fids[0])
					//	consecutive_eids.insert(edgeid);
					if (edge.Vids[0] == e.Vids[0] || edge.Vids[1] == e.Vids[1] || edge.Vids[0] == e.Vids[1] || edge.Vids[1] == e.Vids[0]) {
						bool found = true;
						for (auto nfid : edge.N_Fids) {
							for (auto nnfid : e.N_Fids)
								if (nnfid == nfid) {
									found = false;
									break;
								}
							if (!found)
								break;
						}
						if (!found) continue;
						consecutive_eids.insert(edgeid);
					}
				}
			}
			for (auto eid : consecutive_eids)
				e.consecutiveEids.push_back(eid);
		}
	}
}

void Mesh::BuildOrthogonalE() {
    for (size_t i = 0; i < E.size(); i++) {
        Edge& edge = E.at(i);
        std::vector<size_t> faceEdges;
        for (size_t j = 0; j < E[i].N_Fids.size(); j++) {
            size_t fid = E[i].N_Fids[j];
            for (size_t k = 0; k < 4; k++)
                faceEdges.push_back(F[fid].Eids[k]);
        }
		Util::set_redundent_clearn(faceEdges);
        size_t vid1 = edge.Vids[0];
        size_t vid2 = edge.Vids[1];
        for (size_t j = 0; j < faceEdges.size(); j++) {
            bool hasVid1 = false;
            bool hasVid2 = false;
            const Edge& faceEdge = E.at(faceEdges.at(j));
            for (size_t k = 0; k < faceEdge.Vids.size(); k++) {
                if (faceEdge.Vids.at(k) == vid1) hasVid1 = true;
                if (faceEdge.Vids.at(k) == vid2) hasVid2 = true;
            }
            if (hasVid1 ^ hasVid2) edge.orthogonalEids.push_back(faceEdge.id);
        }
    }
}

void Mesh::unifyOrientation() {
    
    Vertex& v0_r_f = V[F.at(0).Vids[0]];
    Vertex& v1_r_f = V[F.at(0).Vids[1]];
    Vertex& v2_r_f = V[F.at(0).Vids[2]];
    
    const glm::dvec3 v10_r_f = v0_r_f.xyz() - v1_r_f.xyz();
    const glm::dvec3 v12_r_f = v2_r_f.xyz() - v1_r_f.xyz();
    const glm::dvec3 n_r_f = glm::normalize(glm::cross(v12_r_f, v10_r_f));
    
    int fsize = F.size();
    for (int i = 1; i < fsize; i++) {
        Vertex& v0 = V[F.at(i).Vids[0]];
        Vertex& v1 = V[F.at(i).Vids[1]];
        Vertex& v2 = V[F.at(i).Vids[2]];

        const glm::dvec3 v10 = v0.xyz() - v1.xyz();
        const glm::dvec3 v12 = v2.xyz() - v1.xyz();
        const glm::dvec3 n = glm::normalize(glm::cross(v12, v10));
        if (glm::dot(n_r_f, n) < 0.f) {
            std::reverse(F.at(i).Vids.begin(), F.at(i).Vids.end());
        }
    }

    Vertex& v0_r = V[C.at(0).Vids[0]];
    Vertex& v1_r = V[C.at(0).Vids[1]];
    Vertex& v2_r = V[C.at(0).Vids[2]];

    const glm::dvec3 v10_r = v0_r.xyz() - v1_r.xyz();
    const glm::dvec3 v12_r = v2_r.xyz() - v1_r.xyz();
    const glm::dvec3 n_r = glm::normalize(glm::cross(v12_r, v10_r));
    int csize = C.size();
    for (int i = 1; i < csize; i++) {
        Vertex& v0 = V[C.at(i).Vids[0]];
        Vertex& v1 = V[C.at(i).Vids[1]];
        Vertex& v2 = V[C.at(i).Vids[2]];

        const glm::dvec3 v10 = v0.xyz() - v1.xyz();
        const glm::dvec3 v12 = v2.xyz() - v1.xyz();
        const glm::dvec3 n = glm::normalize(glm::cross(v12, v10));
        
        if (glm::dot(n_r, n) < 0.f) {
            std::reverse(C.at(i).Vids.begin(), C.at(i).Vids.end());
        }
    }
}


void Mesh::GetNormalOfSurfaceFaces() {
    for (size_t i = 0; i < F.size(); i++) {
        Face& face = F.at(i);
        const Vertex& v0 = V[face.Vids[0]];
        const Vertex& v1 = V[face.Vids[1]];
        const Vertex& v2 = V[face.Vids[2]];

        const glm::dvec3 v10 = v0.xyz() - v1.xyz();
        const glm::dvec3 v12 = v2.xyz() - v1.xyz();
        const glm::dvec3 n = glm::normalize(glm::cross(v12, v10));
        face.normal = n;
    }
}

void Mesh::GetNormalOfSurfaceVertices() {
    for (size_t i = 0; i < V.size(); i++) {
        Vertex& v = V.at(i);
        v.normal = glm::dvec3(0.0, 0.0, 0.0);
        if (IsVolumetricMesh() && !v.isBoundary) continue;
        glm::dvec3 sumNormal(0.0, 0.0, 0.0);
        size_t faceCount = 0;
        for (size_t j = 0; j < v.N_Fids.size(); j++) {
            const Face& face = F.at(v.N_Fids.at(j));
            if (face.isBoundary) {
                faceCount++;
                v.normal.x += face.normal.x;
                v.normal.y += face.normal.y;
                v.normal.z += face.normal.z;
            }
        }
        v.normal = glm::normalize(v.normal);
    }
}

bool Mesh::IsSurfaceMesh() const {
	return m_cellType == TRIANGLE || m_cellType == QUAD || m_cellType == POLYGON;
}

bool Mesh::IsVolumetricMesh() const {
	return m_cellType == HEXAHEDRA || m_cellType == TETRAHEDRA || m_cellType == POLYHEDRA || m_cellType == PENTAHEDRON;
}

#define CROSSVECTOR3(a,b,c)       {(a)[0]=(b)[1]*(c)[2]-(b)[2]*(c)[1]; \
    (a)[1]=(b)[2]*(c)[0]-(b)[0]*(c)[2]; \
    (a)[2]=(b)[0]*(c)[1]-(b)[1]*(c)[0];}
const double PAI = 3.1415926535898;
void Mesh::ClassifyVertexTypes() {
    std::vector<float> face_angles;
    for (size_t i = 0; i < E.size(); i++) {
        const Edge& edge = E.at(i);
        size_t v1 = edge.Vids[0];
        size_t v2 = edge.Vids[1];

        std::vector<size_t> the_other_vs;
        for (size_t j = 0; j < edge.N_Fids.size(); j++) {
            const Face& face = F.at(edge.N_Fids.at(j));
            if (face.isBoundary) for (size_t k = 0; k < 3; k++) {
                const size_t vid = face.Vids[k];
                if (vid != v1 && vid != v2) the_other_vs.push_back(vid);
            }
        }
        std::vector<size_t> V1s, V2s;
        V1s.push_back(v1);
        V1s.push_back(v2);
        V1s.push_back(the_other_vs[0]);

        V2s.push_back(v1);
        V2s.push_back(v2);
        V2s.push_back(the_other_vs[1]);

        double edge1vector1[3], edge2vector1[3], edge1vector2[3], edge2vector2[3];

        for (size_t j = 0; j < 3; j++) {
            edge1vector1[j] = V[V1s[1]][j] - V[V1s[0]][j];
            edge2vector1[j] = V[V1s[2]][j] - V[V1s[0]][j];

            edge1vector2[j] = V[V2s[1]][j] - V[V2s[0]][j];
            edge2vector2[j] = V[V2s[2]][j] - V[V2s[0]][j];
        }

        glm::dvec3 normal1, normal2;
        CROSSVECTOR3(normal1, edge1vector1, edge2vector1);
        CROSSVECTOR3(normal2, edge2vector2, edge1vector2);
        normal1 = glm::normalize(normal1);
        normal2 = glm::normalize(normal2);

        double cos_angle = normal1[0] * normal2[0] + normal1[1] * normal2[1] + normal1[2] * normal2[2];
        if (cos_angle < -1.0) cos_angle = -1.0;
        else if (cos_angle > 1.0) cos_angle = 1.0;

        double face_angle = acos(cos_angle);

        if (face_angle > PAI) face_angle = PAI;
        else if (face_angle < 0) face_angle = 0;

        face_angles.push_back(PAI - face_angle);

        E[i].face_angle = PAI - face_angle;
    }
    std::vector<std::vector<size_t> > flags(V.size());
    const double angle_threshold = 170.0 / 180 * PAI;
    for (size_t i = 0; i < E.size(); i++) {
        if (E[i].face_angle < angle_threshold) {
            flags[E[i].Vids[0]].push_back(i);
            flags[E[i].Vids[1]].push_back(i);
        }
    }
    for (size_t i = 0; i < flags.size(); i++) {
        if (flags[i].size() == 0) V[i].type = REGULAR;
        else if (flags[i].size() <= 2) {
            V[i].type = FEATURE;
            std::vector<size_t> vs_order;
            for (size_t j = 0; j < flags[i].size(); j++) {
                if (j == 0) {
                    size_t v1 = E[flags[i][j]].Vids[0];
                    size_t v2 = E[flags[i][j]].Vids[1];
                    if (v1 == i) {
                        vs_order.push_back(v2);
                        vs_order.push_back(v1);
                    } else {
                        vs_order.push_back(v1);
                        vs_order.push_back(v2);
                    }
                } else {
                    int v1 = E[flags[i][j]].Vids[0];
                    int v2 = E[flags[i][j]].Vids[1];
                    if (v1 == i) {
                        vs_order.push_back(v1);
                        vs_order.push_back(v2);
                    } else {
                        vs_order.push_back(v2);
                        vs_order.push_back(v1);
                    }
                }
            }
            float length_all = 0;
            for (size_t j = 1; j < vs_order.size(); j++) {
                const glm::dvec3 dir = V[vs_order[j - 1]] - V[vs_order[j]];
                const float dis = glm::length(dir);
                V[i].tangent += glm::dvec3(dis * dir.x, dis * dir.y, dis * dir.z);
            }
            //vector_normalization(tmi_sur.Vs[i].tangent);
        } else if (flags[i].size() > 2) V[i].type = CORNER; //tmi_sur.Vs[i].type=1;//
    }
}
void Mesh::ClearLabelOfSurface() {
    for (size_t i = 0; i < F.size(); i++) {
        Face& face = F.at(i);
        face.label = MAXID;
    }
    for (size_t i = 0; i < E.size(); i++) {
        Edge& edge = E.at(i);
        edge.isSharpFeature = false;
    }
    for (size_t i = 0; i < V.size(); i++) {
        Vertex& v = V.at(i);
        v.type = MAXID;
        v.isCorner = false;
        v.tangent = glm::dvec3(0.0, 0.0, 0.0);
    }
}

void Mesh::LabelSurface() {
    size_t label = 0;
    for (size_t i = 0; i < F.size(); i++) {
        Face& face = F.at(i);
        if (!face.isBoundary || face.label != MAXID) continue;
        LabelFace(face, label);
        label++;
    }
    numberOfPatches = label;

    for (size_t i = 0; i < V.size(); i++) {
        Vertex& v = V.at(i);
        if (v.isBoundary) v.type = REGULAR;
    }
    bool hasBoundary = false;
    for (auto& e : E) {
        if (e.N_Fids.size() == 1) {
            hasBoundary = true;
            break;
        }
    }
    if (hasBoundary) return;
    std::cout << "E size " << E.size() << std::endl;
    for (size_t i = 0; i < E.size(); i++) {
        Edge& edge = E.at(i);
        if (!edge.isBoundary) continue;
        Face* pFace1 = NULL;
        for (size_t j = 0; j < edge.N_Fids.size(); j++) {
            Face& face = F.at(edge.N_Fids.at(j));
            if (face.isBoundary) pFace1 = &face;
        }
        Face* pFace2 = NULL;
        for (size_t j = 0; j < edge.N_Fids.size(); j++) {
            Face& face = F.at(edge.N_Fids.at(j));
            if (face.isBoundary && face.id != pFace1->id) pFace2 = &face;
        }

        const double cos_angle = GetCosAngle(edge, *pFace1, *pFace2);
        // cos(15) = 0.9659 cos(30) = 0.866 cos(45) = 0.707

        if (pFace1->label != pFace2->label || fabs(cos_angle) < cos_angle_threshold) {
            edge.isSharpFeature = true;
            V.at(edge.Vids[0]).type = FEATURE;
            V.at(edge.Vids[1]).type = FEATURE;
        }
    }
    std::cout << "V size: " << V.size() << std::endl;
    for (size_t i = 0; i < V.size(); i++) {
        Vertex& v = V.at(i);
        if (!v.isBoundary) continue;
        size_t numOfSharpEdges = 0;
        for (size_t j = 0; j < v.N_Eids.size(); j++) {
            const Edge& edge = E.at(v.N_Eids.at(j));
            if (edge.isBoundary && edge.isSharpFeature) numOfSharpEdges++;
        }
        if (numOfSharpEdges >= 3) {
            v.isCorner = true;
            v.type = CORNER;
        } else if (numOfSharpEdges < 3) {
            std::vector<size_t> vs_order;
            numOfSharpEdges = 0;
            for (size_t j = 0; j < v.N_Eids.size(); j++) {
                const Edge& edge = E.at(v.N_Eids.at(j));
                if (edge.isBoundary && edge.isSharpFeature) {
                    size_t v1 = edge.Vids[0];
                    size_t v2 = edge.Vids[1];
                    if (numOfSharpEdges == 0) {
                        if (v1 == v.id) {
                            vs_order.push_back(v2);
                            vs_order.push_back(v1);
                        } else {
                            vs_order.push_back(v1);
                            vs_order.push_back(v2);
                        }
                    } else {
                        if (v1 == i) {
                            vs_order.push_back(v1);
                            vs_order.push_back(v2);
                        } else {
                            vs_order.push_back(v2);
                            vs_order.push_back(v1);
                        }
                    }
                    numOfSharpEdges++;
                }
            }
            for (size_t j = 1; j < vs_order.size(); j++) {
                const glm::dvec3 dir = V[vs_order[j - 1]] - V[vs_order[j]];
                const float dis = glm::length(dir);
                v.tangent += glm::dvec3(dis * dir.x, dis * dir.y, dis * dir.z);
            }
        }
    }
    std::cout << "number of sharp edges " << numOfSharpEdges << std::endl;
}

void Mesh::Label2DSurfaceVertices() {
	if (m_cellType != TRIANGLE && m_cellType != QUAD && m_cellType != POLYGON) return;

	auto diff = 180.0 - feature_angle_threshold;
	auto feature_angle_threshold_concave = 180.0 + diff;
	for (auto& v : V) {
		if (v.isBoundary) {
			v.type = FEATURE;
			auto angle = Util::GetAngle(*this, v, v.N_Fids);
			if (angle < feature_angle_threshold || angle > feature_angle_threshold_concave)
				v.type = CORNER;
		}
	}
}

void Mesh::LabelSharpEdges(const bool breakAtConrer/* = false*/) {
    size_t label = 0;
    for (size_t i = 0; i < E.size(); i++) {
        Edge& edge = E.at(i);
        if (((m_cellType == HEXAHEDRA || m_cellType == TETRAHEDRA) && !edge.isBoundary) ||
                !edge.isSharpFeature || edge.label != MAXID) continue;
        LabelEdge(edge, label, breakAtConrer);
        label++;
    }
    this->numOfSharpEdges = label;
}

void Mesh::LabelEdge(Edge& edge, size_t& label, const bool breakAtConrer/* = false*/) {
    edge.label = label;
	auto& v0 = V.at(edge.Vids[0]);
	auto& v1 = V.at(edge.Vids[1]);
	if (!v0.isCorner) v0.label = label;
	v0.labels.insert(label);
	if (!v1.isCorner) v1.label = label;
	v1.labels.insert(label);
    for (size_t i = 0; i < edge.Vids.size(); i++) {
        const Vertex& v = V.at(edge.Vids.at(i));
        if (v.isCorner && breakAtConrer)
            continue;
        std::vector<Edge*> edges;
        for (size_t j = 0; j < v.N_Eids.size(); j++) {
            Edge& edge2 = E.at(v.N_Eids.at(j));
            if (edge2.id != edge.id 
				// && ((m_cellType == HEXAHEDRA || m_cellType == TETRAHEDRA) && edge2.isBoundary) 
				&& edge.isSharpFeature && edge2.label == MAXID) {
                edges.push_back(&edge2);
            }
        }
        for (size_t j = 0; j < edges.size(); j++)
            LabelEdge(*edges.at(j), label, breakAtConrer);
    }
}

void Mesh::ProjectSharpEdgesTo(const std::vector<FeatureLine>& featureLines) {
    std::vector<std::vector<size_t> > featureFaceIds(featureLines.size());
    for (size_t i = 0; i < featureLines.size(); i++) {
        const FeatureLine& fl = featureLines.at(i);
        for (size_t j = 0; j < fl.Vids.size(); j++) {
            const Vertex& target_v = fl.mesh.V.at(fl.Vids.at(j));
            double min_distance = MAXID;
            size_t closestFaceId = 0;
            for (size_t k = 0; k < F.size(); k++) {
                const Face& face = F.at(k);
                if (!face.isBoundary) continue;
                double distance = 0;
                for (size_t l = 0; l < face.Vids.size(); l++) {
                    const Vertex& v = V.at(face.Vids.at(l));
                    distance += glm::length(v.xyz() - target_v.xyz());
                    if (distance > min_distance) break;
                }
                if (distance < min_distance) {
                    min_distance = distance;
                    closestFaceId = face.id;
                }
            }
            featureFaceIds.at(i).push_back(closestFaceId);
        }
        std::sort(featureFaceIds.at(i).begin(), featureFaceIds.at(i).end());
        std::vector<size_t>::iterator iter = std::unique(featureFaceIds.at(i).begin(), featureFaceIds.at(i).end());
        featureFaceIds.at(i).resize(std::distance(featureFaceIds.at(i).begin(), iter));
        std::vector<Face> faces;
        for (size_t j = 0; j < featureFaceIds.at(i).size(); j++) {
            faces.push_back(F.at(featureFaceIds.at(i).at(j)));
        }
        std::string filename = std::string("FeatureFace") + std::to_string(i) + ".vtk";
        MeshFileWriter writer(V, faces, filename.c_str(), QUAD);
        writer.WriteFacesVtk();
    }
}

//double Mesh::GetCosAngle(const Edge& edge, const Face& face1, const Face& face2)
//{
//    /*              v0
//     *             /|\
//     *            / | \
//     *      face1/  |  \face2
//     *          /   |   \
//     *         /    |    \
//     *        -------------
//     *       v2     v1     v3
//     * */
//    const Vertex& v1 = V.at(edge.Vids[0]);
//    const Vertex& v0 = V.at(edge.Vids[1]);
//    size_t v2id = MAXID;
//    for (size_t i = 0; i < face1.Eids.size(); i++) {
//        const Edge& edge1 = E.at(face1.Eids.at(i));
//        if (edge == edge1)
//            continue;
//        if (edge1.Vids[0] == v1.id)
//            v2id = edge1.Vids[1];
//        else if (edge1.Vids[1] == v1.id)
//            v2id = edge1.Vids[0];
//        if (v2id != MAXID)
//            break;
//    }
//    const Vertex& v2 = V.at(v2id);
//
//    size_t v3id = MAXID;
//    for (size_t i = 0; i < face2.Eids.size(); i++) {
//        const Edge& edge2 = E.at(face2.Eids.at(i));
//        if (edge == edge2)
//            continue;
//        if (edge2.Vids[0] == v1.id)
//            v3id = edge2.Vids[1];
//        else if (edge2.Vids[1] == v1.id)
//            v3id = edge2.Vids[0];
//        if (v3id != MAXID)
//            break;
//    }
//    const Vertex& v3 = V.at(v3id);
//    Plane plane1(v2, v1, v0);
//    const Plane plane2(v3, v0, v1);
//
//    return plane1.IntersectionAngle(plane2);
//}

double Mesh::GetCosAngle(const Edge& edge, const Face& face1, const Face& face2)
{
    /*              v0
     *             /|\
     *            / | \
     *      face1/  |  \face2
     *          /   |   \
     *         /    |    \
     *        -------------
     *       v2     v1     v3
     * */
    Vertex* pv0 = &V.at(edge.Vids[0]);
    Vertex* pv1 = &V.at(edge.Vids[1]);
    bool correct_orientation = false;
    const size_t n = face1.Vids.size();
    for (size_t i = 0; i < n; i++) {
        const Vertex& v0 = V.at(face1.Vids[i % n]);
        const Vertex& v1 = V.at(face1.Vids[(i+1) % n]);
        if (v0.id == pv0->id && v1.id == pv1->id) {
            correct_orientation = true;
            break;
        }
    }
    if (!correct_orientation) {
        pv1 = &V.at(edge.Vids[0]);
        pv0 = &V.at(edge.Vids[1]);
    }
    const Vertex& v0 = *pv0;
    const Vertex& v1 = *pv1;

    size_t v2id = MAXID;
    for (size_t i = 0; i < face1.Eids.size(); i++) {
        const Edge& edge1 = E.at(face1.Eids.at(i));
        if (edge == edge1)
            continue;
        if (edge1.Vids[0] == v1.id)
            v2id = edge1.Vids[1];
        else if (edge1.Vids[1] == v1.id)
            v2id = edge1.Vids[0];
        if (v2id != MAXID)
            break;
    }
    const Vertex& v2 = V.at(v2id);

    size_t v3id = MAXID;
    for (size_t i = 0; i < face2.Eids.size(); i++) {
        const Edge& edge2 = E.at(face2.Eids.at(i));
        if (edge == edge2)
            continue;
        if (edge2.Vids[0] == v1.id)
            v3id = edge2.Vids[1];
        else if (edge2.Vids[1] == v1.id)
            v3id = edge2.Vids[0];
        if (v3id != MAXID)
            break;
    }
    const Vertex& v3 = V.at(v3id);
    Plane plane1(v2, v1, v0);
    const Plane plane2(v3, v0, v1);

    return plane1.IntersectionAngle(plane2);
}

double Mesh::Get2DAngle(const Vertex& v, const Edge& e1, const Edge& e2) {
	/*             v
	*             / \
	*            /   \
	*         e1/     \e2
	*          /       \
	*         /         \
	*        /           \
	*       v1           v2
	* */
	auto v1id = e1.Vids[0] == v.id ? e1.Vids[1] : e1.Vids[0];
	auto v2id = e2.Vids[0] == v.id ? e2.Vids[1] : e2.Vids[0];
	auto& v1 = V.at(v1id);
	auto& v2 = V.at(v2id);
	auto vv1 = v1 - v;
	auto vv2 = v2 - v;

	return acos(glm::dot(vv1, vv2));
}

void Mesh::SetCosAngleThreshold(const double cos_angle/* = 0.91*/) {
    cos_angle_threshold = cos_angle;
}

void Mesh::SetFeatureAngleThreshold(const double angle/* = 170.0*/) {
    feature_angle_threshold = angle;
}

void Mesh::LabelFace(Face& face, size_t& label) {
    std::vector<Face*> queue;
    queue.push_back(&face);

    while (!queue.empty()) {
        Face& f = *queue.back();
        queue.pop_back();
        f.label = label;

        for (size_t i = 0; i < f.Eids.size(); i++) {
            auto& edge = E.at(f.Eids.at(i));
            for (size_t j = 0; j < edge.N_Fids.size(); j++) {
                Face& face2 = F.at(edge.N_Fids.at(j));
                if (face2.isBoundary && face2.id != f.id && face2.label == MAXID) {
                    const double cos_angle = GetCosAngle(edge, f, face2);
                    if (cos_angle > cos_angle_threshold) // cos(15) = 0.9659 cos(30) = 0.866
                        queue.push_back(&face2);
                    //else edge.isSharpFeature = true;
                }
            }
        }

    }
    // face.label = label;
    // // std::cout << face.id << std::endl;
    // for (size_t i = 0; i < face.Eids.size(); i++) {
    //     auto& edge = E.at(face.Eids.at(i));
    //     std::vector<Face*> faces;
    //     for (size_t j = 0; j < edge.N_Fids.size(); j++) {
    //         Face& face2 = F.at(edge.N_Fids.at(j));
    //         if (face2.isBoundary && face2.id != face.id && face2.label == MAXID) {
    //             const double cos_angle = GetCosAngle(edge, face, face2);
    //             if (cos_angle > cos_angle_threshold) // cos(15) = 0.9659 cos(30) = 0.866
	// 				faces.push_back(&face2);
	// 			//else edge.isSharpFeature = true;
    //         }
    //     }
    //     for (size_t i = 0; i < faces.size(); i++)
    //         LabelFace(*faces.at(i), label);
    // }
}

void Mesh::RemoveUselessVertices() {
    std::vector<size_t> v_real_index;
    size_t c_size = 0;
    for (size_t i = 0; i < C.size(); i++) {
        const Cell& cell = C.at(i);
        for (size_t j = 0; j < cell.Vids.size(); j++)
            v_real_index.push_back(cell.Vids.at(j));
        c_size++;
    }

    std::sort(v_real_index.begin(), v_real_index.end());
    std::vector<size_t>::iterator iter = std::unique(v_real_index.begin(), v_real_index.end());
    v_real_index.resize(std::distance(v_real_index.begin(), iter));

    std::vector<Vertex> newV(v_real_index.size());
    for (size_t i = 0; i < v_real_index.size(); i++) {
        const Vertex& v = V.at(v_real_index.at(i));
        Vertex& newv = newV.at(i);
        newv.id = i;
        newv.x = v.x;
        newv.y = v.y;
        newv.z = v.z;
    }
    V = newV;
    m_refIds = v_real_index;
    //////////////////////////////////////////////////////
    std::map<size_t, size_t> v_v;
    size_t index = 0;
    for (size_t i = 0; i < v_real_index.size(); i++) {
        v_v[v_real_index.at(i)] = i;
    }

    std::vector<Cell> newC;
    for (size_t i = 0; i < C.size(); i++) {
        Cell c;
        if (m_cellType == HEXAHEDRA) c.Vids.resize(8);
        else if (m_cellType == TETRAHEDRA) c.Vids.resize(4);
        else if (m_cellType == QUAD) c.Vids.resize(4);
        else if (m_cellType == TRIANGLE) c.Vids.resize(3);
        const Cell& cell = C.at(i);
        for (size_t j = 0; j < C.at(i).Vids.size(); j++) {
            c.Vids.at(j) = v_v[cell.Vids.at(j)];
            c.id = i;
        }
		c.cellType = cell.cellType;
        newC.push_back(c);
    }
    C.resize(newC.size());
    C = newC;

    if (m_cellType == QUAD) {
        std::vector<Face> newF;
        for (size_t i = 0; i < F.size(); i++) {
            Face f;
            if (m_cellType == HEXAHEDRA) f.Vids.resize(8);
            else if (m_cellType == TETRAHEDRA) f.Vids.resize(4);
            else if (m_cellType == QUAD) f.Vids.resize(4);
            else if (m_cellType == TRIANGLE) f.Vids.resize(3);
            const Face& cell = F.at(i);
            for (size_t j = 0; j < F.at(i).Vids.size(); j++) {
                f.Vids.at(j) = v_v[cell.Vids.at(j)];
                f.id = i;
            }
            newF.push_back(f);
        }
        F.resize(newC.size());
        F = newF;
    }
}

void Mesh::CompressWithFeaturePreserved() {
    std::set<size_t> vids;
    for (const auto& c : C)
        vids.insert(c.Vids.begin(), c.Vids.end());
    std::vector<size_t> v_real_index;
    std::copy(vids.begin(), vids.end(), std::back_inserter(v_real_index));

    std::vector<Vertex> newV(v_real_index.size());
    for (size_t i = 0; i < v_real_index.size(); i++) {
        const Vertex& v = V.at(v_real_index.at(i));
        Vertex& newv = newV.at(i);
        newv = v;
        newv.id = i;
        newv.type = v.type;
        newv.isCorner = v.isCorner;
		newv.label = v.label;
		newv.patch_id = v.patch_id;
		newv.isSpecial = v.isSpecial;
		newv.isConvex = v.isConvex;
		newv.labels = v.labels;
		newv.patch_ids = v.patch_ids;
		newv.idealValence = v.idealValence;
    }
    V = newV;
    for (size_t i = 0; i < newV.size(); i++) {
        Vertex& v = V.at(i);
        Vertex& newv = newV.at(i);
        v.id = i;
        v.type = newv.type;
        v.isCorner = newv.isCorner;
		v.label = newv.label;
		v.patch_id = newv.patch_id;
		v.isSpecial = newv.isSpecial;
		v.isConvex = newv.isConvex;
		v.labels = newv.labels;
		v.patch_ids = newv.patch_ids;
        v.idealValence = newv.idealValence;
    }
    m_refIds = v_real_index;
    //////////////////////////////////////////////////////
    std::map<size_t, size_t> v_v;
    size_t index = 0;
    for (size_t i = 0; i < v_real_index.size(); i++) {
        v_v[v_real_index.at(i)] = i;
    }

    std::vector<Cell> newC;
    newC.reserve(C.size());
    Cell cell;
    for (auto& c : C) {
        cell.Vids = c.Vids;
        for (auto& vid : cell.Vids)
            vid = v_v[vid];
        cell.id = c.id;
        cell.cellType = c.cellType;
        newC.push_back(cell);
    }
    C = newC;

//    if (m_cellType == QUAD) {
//        std::vector<Face> newF;
//        for (size_t i = 0; i < F.size(); i++) {
//            Face f;
//            if (m_cellType == HEXAHEDRA) f.Vids.resize(8);
//            else if (m_cellType == TETRAHEDRA) f.Vids.resize(4);
//            else if (m_cellType == QUAD) f.Vids.resize(4);
//            else if (m_cellType == TRIANGLE) f.Vids.resize(3);
//            const Face& cell = F.at(i);
//            for (size_t j = 0; j < F.at(i).Vids.size(); j++) {
//                f.Vids.at(j) = v_v[cell.Vids.at(j)];
//                f.id = i;
//            }
//            newF.push_back(f);
//        }
//        F.resize(newC.size());
//        F = newF;
//    }
}

bool Mesh::IsPointInTriangle(const glm::dvec3& P, const glm::dvec3& A, const glm::dvec3& B, const glm::dvec3& C) const {
    // Compute vectors
    const glm::dvec3 v0 = C - A;
    const glm::dvec3 v1 = B - A;
    const glm::dvec3 v2 = P - A;

    // Compute dot products
    const double dot00 = glm::dot(v0, v0);
    const double dot01 = glm::dot(v0, v1);
    const double dot02 = glm::dot(v0, v2);
    const double dot11 = glm::dot(v1, v1);
    const double dot12 = glm::dot(v1, v2);

    // Compute barycentric coordinates
    const double invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);
    const double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    const double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

    // Check if point is in triangle
    //return (u >= 0) && (v >= 0) && (u + v <= 1.0);
    return (u >= -2e-1) && (v >= -2e-1) && (u + v <= 1.3);
}

bool Mesh::IsPointInFace(const glm::dvec3& P, const Face& face) const
{
    const glm::dvec3& A = V.at(face.Vids[0]).xyz();
    const glm::dvec3& B = V.at(face.Vids[1]).xyz();
    const glm::dvec3& C = V.at(face.Vids[2]).xyz();

    // Compute vectors
    const glm::dvec3 v0 = C - A;
    const glm::dvec3 v1 = B - A;
    const glm::dvec3 v2 = P - A;

    // Compute dot products
    const double dot00 = glm::dot(v0, v0);
    const double dot01 = glm::dot(v0, v1);
    const double dot02 = glm::dot(v0, v2);
    const double dot11 = glm::dot(v1, v1);
    const double dot12 = glm::dot(v1, v2);

    // Compute barycentric coordinates
    const double invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);
    const double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    const double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

    // Check if point is in triangle
    return (u >= -2e-1) && (v >= -2e-1) && (u + v <= 1.3);
    //return (u >= -1e-3) && (v >= -1e-3) && (u + v <= 1.001);
}

//bool Mesh::IsPointInTriangle(const Vertex& p, const Vertex& p0, const Vertex& p1, const Vertex& p2) const
//{
//    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(3, 3);
//    Eigen::VectorXd b = Eigen::VectorXd::Zero(3);
//    A(0,0) = p0.x; A(0,1) = p1.x; A(0,2) = p2.x;
//    A(1,0) = p0.y; A(1,1) = p1.y; A(1,2) = p2.y;
//    A(2,0) = p0.z; A(2,1) = p1.z; A(2,2) = p2.z;
//    b(0) = p.x; b(1) = p.y; b(2) = p.z;
//    Eigen::VectorXd x = A.ldlt().solve(b);
//
//    if (x(0) > -1e-3 && x(1) > -1e-3 && x(2) > -1e-3 && (x(0) + x(1) + x(2)) < 1.001)
//    {
//        return true;
//    }
//
//    return false;
//}
//
//double Area(const glm::dvec3& p0, const glm::dvec3& p1, const glm::dvec3& p2)
//{
//    const glm::dvec3 v10 = p0 - p1;
//    const glm::dvec3 v12 = p2 - p1;
//    const glm::dvec3 normal = glm::cross(v12, v10);
//    return 0.5f * glm::length(normal);
//}
//
//bool Mesh::IsPointInTriangle(const glm::dvec3& p, const glm::dvec3& p0, const glm::dvec3& p1, const glm::dvec3& p2) const
//{
////    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(4, 3);
////    Eigen::VectorXd b = Eigen::VectorXd::Zero(4);
////    A(0,0) = p0.x; A(0,1) = p1.x; A(0,2) = p2.x;
////    A(1,0) = p0.y; A(1,1) = p1.y; A(1,2) = p2.y;
////    A(2,0) = p0.z; A(2,1) = p1.z; A(2,2) = p2.z;
////    A(3,0) = 1.0 ; A(3,1) = 1.0 ; A(3,2) = 1.0 ;
////    b(0) = p.x; b(1) = p.y; b(2) = p.z; b(3) = 1.0;
////    Eigen::MatrixXd AT = A.transpose();
////    Eigen::MatrixXd ATA = AT*A;
////    Eigen::VectorXd x = ATA.ldlt().solve(AT*b);
////
////    if (x(0) > -1e-1 && x(1) > -1e-1 && x(2) > -1e-1 && (x(0) + x(1) + x(2)) < 1.1)
////    {
////        return true;
////    }
//    const double a0 = Area(p, p0, p1);
//    const double a1 = Area(p, p1, p2);
//    const double a2 = Area(p, p2, p0);
//    const double a  = Area(p0, p1, p2);
//
//    if (fabs(a - (a0 + a1 + a2))/a < 2e-1)
//        return true;
//    return false;
//}

struct d_i {
    double d;
    glm::dvec3 i;
    d_i(): d(0.0){}
    d_i(const double d, const glm::dvec3 i)
    : d(d), i(i) {}
    d_i(const d_i& r)
    : d(r.d), i(r.i) {}
    bool operator < (const d_i& r) const {
        return d < r.d;
    }
};

glm::dvec3 Mesh::GetProjectLocation(const glm::dvec3& p) const {
    for (size_t i = 0; i < F.size(); i++) {
        const Face& face = F.at(i);
        if (!face.isBoundary) continue;
        Plane plane(V[face.Vids[0]].xyz(), V[face.Vids[1]].xyz(), V[face.Vids[2]].xyz());
        const double d1 = glm::length(V[face.Vids[0]].xyz() - V[face.Vids[1]].xyz());
        const double d2 = glm::length(V[face.Vids[1]].xyz() - V[face.Vids[2]].xyz());
        const double d3 = glm::length(V[face.Vids[2]].xyz() - V[face.Vids[0]].xyz());
        const double d = 0.33333333 * (d1 + d2 + d3);
        glm::dvec3 intersection;
        const double distance = plane.DistanseFromPoint(p, intersection);
        if (distance < 1.0 * avgEdgeLength && distance < 1.0 * d) if (IsPointInFace(intersection, face))
        //if (IsPointInTriangle(intersection, V[face.Vids[0]].xyz(), V[face.Vids[1]].xyz(), V[face.Vids[2]].xyz()))
        return intersection;
    }
    //std::cerr << "\033[1;31mERROR in GetProjectLocation(" << p.x << "," << p.y << "," << p.z << ")\033[0m\n";
    return glm::dvec3(0.0, 0.0, 0.0);
}
glm::dvec3 Mesh::GetProjectLocationOnRefSurface(const glm::dvec3& p, const Vertex& refV) const {
    std::vector<d_i> d_is;
    for (size_t i = 0; i < refV.twoRingNeighborSurfaceFaceIds.size(); i++) {
        const Face& face = F.at(refV.twoRingNeighborSurfaceFaceIds[i]);
        if (!face.isBoundary) continue;
        for (size_t j = 0; j < face.Vids.size(); j++) {
            Plane plane(V[face.Vids[(0 + j) % 4]].xyz(), V[face.Vids[(1 + j) % 4]].xyz(), V[face.Vids[(2 + j) % 4]].xyz());
            const double d1 = glm::length(V[face.Vids[0]].xyz() - V[face.Vids[1]].xyz());
            const double d2 = glm::length(V[face.Vids[1]].xyz() - V[face.Vids[2]].xyz());
            const double d3 = glm::length(V[face.Vids[2]].xyz() - V[face.Vids[0]].xyz());
            const double d = 0.33333333 * (d1 + d2 + d3);
            glm::dvec3 intersection;
            const double distance = plane.DistanseFromPoint(p, intersection);
            if (distance < 1.0 * avgEdgeLength && distance < 1.0 * d) if (IsPointInFace(intersection, face)) d_is.push_back(d_i(distance, intersection));
//            if (face.Vids.size() == 3)
//                break;
        }
    }
    if (d_is.empty()) {
        //std::cerr << "\033[1;31mERROR in GetProjectLocationOfVertex " << p.id << "(" << p.x << "," << p.y << "," << p.z << ")\033[0m\n";
        return glm::dvec3(0.0, 0.0, 0.0);
    }

    std::sort(d_is.begin(), d_is.end());
    //std::vector<d_i>::iterator iter = std::unique(d_is.begin(), d_is.end());
    return d_is[0].i;
}

glm::dvec3 Mesh::GetProjectLocation(const Vertex& p) const {
    std::vector<d_i> d_is;
    for (size_t i = 0; i < F.size(); i++) {
        const Face& face = F.at(i);
        if (!face.isBoundary) continue;
        for (size_t j = 0; j < face.Vids.size(); j++) {
            Plane plane(V[face.Vids[(0 + j) % 4]].xyz(), V[face.Vids[(1 + j) % 4]].xyz(), V[face.Vids[(2 + j) % 4]].xyz());
            const double d1 = glm::length(V[face.Vids[0]].xyz() - V[face.Vids[1]].xyz());
            const double d2 = glm::length(V[face.Vids[1]].xyz() - V[face.Vids[2]].xyz());
            const double d3 = glm::length(V[face.Vids[2]].xyz() - V[face.Vids[0]].xyz());
            const double d = 0.33333333 * (d1 + d2 + d3);
            glm::dvec3 intersection;
            const double distance = plane.DistanseFromPoint(p, intersection);
            if (distance < 1.0 * avgEdgeLength && distance < 1.0 * d) if (IsPointInFace(intersection, face))
            //if (IsPointInTriangle(intersection, V[face.Vids[(0 + j) % 4]].xyz(), V[face.Vids[(1 + j) % 4]].xyz(), V[face.Vids[(2 + j) % 4]].xyz()))
            //return intersection;
            //if (glm::dot(p.normal, face.normal) > 0)
            d_is.push_back(d_i(distance, intersection));

            if (face.Vids.size() == 3) break;
        }
    }
    if (d_is.empty()) {
        //std::cerr << "\033[1;31mERROR in GetProjectLocationOfVertex " << p.id << "(" << p.x << "," << p.y << "," << p.z << ")\033[0m\n";
        return glm::dvec3(0.0, 0.0, 0.0);
    }

    std::sort(d_is.begin(), d_is.end());
    //std::vector<d_i>::iterator iter = std::unique(d_is.begin(), d_is.end());
    return d_is[0].i;
}
glm::dvec3 Mesh::GetProjectLocationFast(const Vertex& p) const {
    std::vector<d_i> d_is;
    for (size_t i = 0; i < p.twoRingNeighborSurfaceFaceIds.size(); i++) {
        const Face& face = F.at(p.twoRingNeighborSurfaceFaceIds[i]);
        if (!face.isBoundary) continue;
        for (size_t j = 0; j < face.Vids.size(); j++) {
            Plane plane(V[face.Vids[(0 + j) % 4]].xyz(), V[face.Vids[(1 + j) % 4]].xyz(), V[face.Vids[(2 + j) % 4]].xyz());
            const double d1 = glm::length(V[face.Vids[0]].xyz() - V[face.Vids[1]].xyz());
            const double d2 = glm::length(V[face.Vids[1]].xyz() - V[face.Vids[2]].xyz());
            const double d3 = glm::length(V[face.Vids[2]].xyz() - V[face.Vids[0]].xyz());
            const double d = 0.33333333 * (d1 + d2 + d3);
            glm::dvec3 intersection;
            const double distance = plane.DistanseFromPoint(p, intersection);
            if (distance < 1.0 * avgEdgeLength && distance < 1.0 * d) if (IsPointInFace(intersection, face))
            //if (IsPointInTriangle(intersection, V[face.Vids[(0 + j) % 4]].xyz(), V[face.Vids[(1 + j) % 4]].xyz(), V[face.Vids[(2 + j) % 4]].xyz()))
            //return intersection;
            //if (glm::dot(p.normal, face.normal) > 0)
            d_is.push_back(d_i(distance, intersection));

            if (face.Vids.size() == 3) break;
        }
    }
    if (d_is.empty()) {
        //std::cerr << "\033[1;31mERROR in GetProjectLocationOfVertex " << p.id << "(" << p.x << "," << p.y << "," << p.z << ")\033[0m\n";
        return glm::dvec3(0.0, 0.0, 0.0);
    }

    std::sort(d_is.begin(), d_is.end());
    //std::vector<d_i>::iterator iter = std::unique(d_is.begin(), d_is.end());
    return d_is[0].i;
}
glm::dvec3 Mesh::GetProjectLocationOnTargetSurface(const glm::dvec3& p, const Vertex& refV, const Mesh& targetSurfaceMesh) const {
    std::vector<d_i> d_is;
    for (size_t i = 0; i < refV.twoRingNeighborSurfaceFaceIds.size(); i++) {
        const Face& face = targetSurfaceMesh.F.at(refV.twoRingNeighborSurfaceFaceIds[i]);
        if (!face.isBoundary) continue;
        for (size_t j = 0; j < face.Vids.size(); j++) {
            Plane plane(targetSurfaceMesh.V[face.Vids[(0 + j) % 4]].xyz(), targetSurfaceMesh.V[face.Vids[(1 + j) % 4]].xyz(), targetSurfaceMesh.V[face.Vids[(2 + j) % 4]].xyz());
            const double d1 = glm::length(targetSurfaceMesh.V[face.Vids[0]].xyz() - targetSurfaceMesh.V[face.Vids[1]].xyz());
            const double d2 = glm::length(targetSurfaceMesh.V[face.Vids[1]].xyz() - targetSurfaceMesh.V[face.Vids[2]].xyz());
            const double d3 = glm::length(targetSurfaceMesh.V[face.Vids[2]].xyz() - targetSurfaceMesh.V[face.Vids[0]].xyz());
            const double d = 0.33333333 * (d1 + d2 + d3);
            glm::dvec3 intersection;
            const double distance = plane.DistanseFromPoint(p, intersection);
            if (distance < 1.0 * avgEdgeLength && distance < 1.0 * d) if (targetSurfaceMesh.IsPointInFace(intersection, face)) d_is.push_back(d_i(distance, intersection));
//            if (face.Vids.size() == 3)
//                break;
        }
    }
    if (d_is.empty()) {
        //std::cerr << "\033[1;31mERROR in GetProjectLocationOfVertex " << p.id << "(" << p.x << "," << p.y << "," << p.z << ")\033[0m\n";
        return glm::dvec3(0.0, 0.0, 0.0);
    }

    std::sort(d_is.begin(), d_is.end());
    //std::vector<d_i>::iterator iter = std::unique(d_is.begin(), d_is.end());
    return d_is[0].i;
}

void Mesh::ProjectTo(const Mesh& mesh) {
    glm::dvec3 vm(0.0, 0.0, 0.0);
    GetAvgEdgeLength();
#pragma omp parallel for
    for (size_t i = 0; i < V.size(); i++) {
        Vertex& v = V.at(i);
        if (!v.isBoundary) continue;
        const glm::dvec3 newv = mesh.GetProjectLocation(v);
        if (vm != newv) v = newv;
    }
}

void Mesh::FastProjectTo(const Mesh& mesh) {
    glm::dvec3 vm(0.0, 0.0, 0.0);
    GetAvgEdgeLength();
#pragma omp parallel for
    for (size_t i = 0; i < V.size(); i++) {
        Vertex& v = V.at(i);
        if (!v.isBoundary) continue;
        const glm::dvec3 newv = mesh.GetProjectLocationFast(v);
        if (vm != newv) v = newv;
    }
}

void Mesh::ProjectToRefMesh(const Mesh& refMesh) {
    glm::dvec3 vm(0.0, 0.0, 0.0);
    GetAvgEdgeLength();
#pragma omp parallel for
    for (size_t i = 0; i < V.size(); i++) {
        Vertex& v = V.at(i);
        if (!v.isBoundary) continue;
        const Vertex& refV = refMesh.V.at(m_refIds[i]);
        if (!refV.isBoundary) continue;
        const glm::dvec3 newv = refMesh.GetProjectLocation(refV);
        if (vm != newv) v = newv;
    }
}

void Mesh::ProjectToTargetSurface(const Mesh& refMesh, const Mesh& targetSurfaceMesh) {
    glm::dvec3 vm(0.0, 0.0, 0.0);
    GetAvgEdgeLength();
#pragma omp parallel for
    for (size_t i = 0; i < V.size(); i++) {
        Vertex& v = V.at(i);
        if (!v.isBoundary) continue;
        const Vertex& refV = refMesh.V.at(m_refIds[i]);
        if (!refV.isBoundary) continue;
        const glm::dvec3 newv = refMesh.GetProjectLocationOnTargetSurface(v.xyz(), refV, targetSurfaceMesh);
        if (vm != newv) v = newv;
    }
}

double Mesh::GetAvgEdgeLength() {
    if (avgEdgeLength != 0) return avgEdgeLength;
    double sum = 0;
    for (auto& e : E) {
        // if (!e.isBoundary) continue;
        e.length = glm::length(V.at(e.Vids[0]).xyz() - V.at(e.Vids[1]).xyz());
		sum += e.length;
    }
    avgEdgeLength = sum / E.size();
    return avgEdgeLength;
}

glm::dvec3 Mesh::GetFaceCenter(const Face& f) {
    glm::dvec3 center(0.0, 0.0, 0.0);
    for (size_t k = 0; k < f.Vids.size(); k++) {
        const Vertex& f_v = V.at(f.Vids[k]);
        center += f_v.xyz();
    }
    center.x /= f.Vids.size();
    center.y /= f.Vids.size();
    center.z /= f.Vids.size();
    return center;
}

double Mesh::SmoothVolume(const SmoothMethod smoothMethod/* = LAPLACE_EDGE*/) {
    std::vector<Vertex> oldV = V;
    std::vector<Vertex> newV = V;
    //while (iters-- != 0)
    {
#pragma omp parallel for
        for (size_t i = 0; i < V.size(); i++) {
            const Vertex& v = V.at(i);
            if (v.isBoundary) continue;
            //if (v.type == REGULAR || (v.type == FEATURE && treatSharpFeatureAsRegular) || (v.type == CORNER && treatCornerAsRegular))
            {
                if (smoothMethod == LAPLACE_EDGE) {
                    newV.at(i) = LapLace(v);
                } else if (smoothMethod == LAPLACE_FACE_CENTER) {
                    glm::dvec3 sum(0.0, 0.0, 0.0);
                    int count = 0;
                    for (size_t j = 0; j < v.N_Fids.size(); j++) {
                        const Face& f = F.at(v.N_Fids[j]);
                        if (!f.isBoundary) continue;
                        sum += GetFaceCenter(f);
                        count++;
                    }
                    newV.at(i).x = sum.x / count;
                    newV.at(i).y = sum.y / count;
                    newV.at(i).z = sum.z / count;
                }
            }
        }

        for (size_t i = 0; i < V.size(); i++) {
            Vertex& v = V.at(i);
            if (v.isBoundary) continue;
            const Vertex& newv = newV.at(i);
            v = newv.xyz();
        }
    }

    double energy = 0;
    for (size_t i = 0; i < V.size(); i++) {
        const Vertex& v = V.at(i);
        if (v.isBoundary) continue;
        const Vertex& oldv = oldV.at(i);
        const double distance = glm::length(v.xyz() - oldv.xyz());
        energy += distance * distance;
    }

    std::cout << "Volume Energy = " << energy << std::endl;
    return energy;
}

double Mesh::SmoothSurface(size_t iters/* = 1*/, const SmoothMethod smoothMethod/* = LAPLACE_EDGE*/,
        const bool preserveSharpFeature/* = false*/, const bool treatSharpFeatureAsRegular/* = false*/, const bool treatCornerAsRegular/* = false*/)
{
    std::vector<Vertex> oldV = V;
    std::vector<Vertex> newV = V;
    while (iters-- != 0)
    {
#pragma omp parallel for
        for (size_t i = 0; i < V.size(); i++) {
            const Vertex& v = V.at(i);
            if (!v.isBoundary)
                continue;
            if (v.type == REGULAR || (v.type == FEATURE && treatSharpFeatureAsRegular) || (v.type == CORNER && treatCornerAsRegular))
            {
                if (smoothMethod == LAPLACE_EDGE) {
                    newV.at(i) = LapLace(v, treatSharpFeatureAsRegular, treatCornerAsRegular);
                }
                else if (smoothMethod == LAPLACE_FACE_CENTER) {
                    glm::dvec3 sum(0.0, 0.0, 0.0);
                    int count = 0;
                    for (size_t j = 0; j < v.N_Fids.size(); j++) {
                        const Face& f = F.at(v.N_Fids[j]);
                        if (!f.isBoundary)
                            continue;
                        sum += GetFaceCenter(f);
                        count++;
                    }
                    newV.at(i).x = sum.x/count;
                    newV.at(i).y = sum.y/count;
                    newV.at(i).z = sum.z/count;
                }
            }
            else if (v.type == FEATURE) {
                if (preserveSharpFeature)
                    continue;
                if (smoothMethod == LAPLACE_EDGE) {
                    newV.at(i) = LapLace(v);
                }
                else if (smoothMethod == LAPLACE_FACE_CENTER) {
                    glm::dvec3 sum(0.0, 0.0, 0.0);
                    int count = 0;
                    for (size_t j = 0; j < v.N_Fids.size(); j++) {
                        const Face& f = F.at(v.N_Fids[j]);
                        if (!f.isBoundary)
                            continue;
                        sum += GetFaceCenter(f);
                        count++;
                    }
                    newV.at(i).x = sum.x/count;
                    newV.at(i).y = sum.y/count;
                    newV.at(i).z = sum.z/count;
                }
            }
        }

        for (size_t i = 0; i < V.size(); i++) {
            Vertex& v = V.at(i);
            if (!v.isBoundary)
                continue;
            const Vertex& newv = newV.at(i);
            v = newv.xyz();
        }
    }

    double energy = 0;
    for (size_t i = 0; i < V.size(); i++) {
        const Vertex& v = V.at(i);
        if (!v.isBoundary)
            continue;
        const Vertex& oldv = oldV.at(i);
        const double distance = glm::length(v.xyz() - oldv.xyz());
        energy += distance * distance;
    }

    std::cout << "Energy = " << energy << std::endl;
    return energy;
}

double Mesh::GetMinScaledJacobian(double& avgSJ) const {
    double minScaledJacobian = 1;
    double sum = 0;
    for (int i = 0; i < C.size(); i++) {
        const Cell& cell = C.at(i);
        float scaledJacobian = GetScaledJacobian(cell);
        minScaledJacobian = minScaledJacobian < scaledJacobian ? minScaledJacobian : scaledJacobian;
        sum += scaledJacobian;
    }
    avgSJ = sum / C.size();
    return minScaledJacobian;
}

size_t Mesh::GetQuality(const char* filename, double& minValue, double& avgValue, const double minSJ/* = 0.0*/) {
    vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName(filename);
    reader->Update();
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGridMesh = reader->GetOutput();
    vtkSmartPointer<vtkMeshQuality> qualityFilter = vtkSmartPointer<vtkMeshQuality>::New();
#if VTK_MAJOR_VERSION <= 5
    qualityFilter->SetInputConnection(mesh->GetProducerPort());
#else
    ///qualityFilter->SetInputData(mesh);
    qualityFilter->SetInputConnection(reader->GetOutputPort());
#endif
    minValue = 0;
    double maxValue = 0;
    avgValue = 0;
    std::vector<double> metrics;
    //qualityFilter->SetTriangleQualityMeasureToArea();
    //qualityFilter->SetHexQualityMeasureToEdgeRatio()
    //qualityFilter->SetHexQualityMeasureToMedAspectFrobenius()
    //qualityFilter->SetHexQualityMeasureToMaxAspectFrobenius()
    //qualityFilter->SetHexQualityMeasureToMaxEdgeRatios()
    //qualityFilter->SetHexQualityMeasureToSkew()
    //qualityFilter->SetHexQualityMeasureToTaper()
    //qualityFilter->SetHexQualityMeasureToVolume()
    //qualityFilter->SetHexQualityMeasureToStretch()
    //qualityFilter->SetHexQualityMeasureToDiagonal()
    //qualityFilter->SetHexQualityMeasureToDimension()
    //qualityFilter->SetHexQualityMeasureToOddy()
    //qualityFilter->SetHexQualityMeasureToCondition()
    //qualityFilter->SetHexQualityMeasureToJacobian()
    qualityFilter->SetHexQualityMeasureToScaledJacobian();
    //qualityFilter->SetHexQualityMeasureToShear()
    //qualityFilter->SetHexQualityMeasureToShape()
    //qualityFilter->SetHexQualityMeasureToRelativeSizeSquared()
    //qualityFilter->SetHexQualityMeasureToShapeAndSize()
    //qualityFilter->SetHexQualityMeasureToShearAndSize()
    //qualityFilter->SetHexQualityMeasureToDistortion()
    qualityFilter->Update();

    vtkSmartPointer<vtkDoubleArray> qualityArray = vtkDoubleArray::SafeDownCast(qualityFilter->GetOutput()->GetCellData()->GetArray("Quality"));
    minValue = qualityArray->GetValue(0);
    maxValue = qualityArray->GetValue(0);
    avgValue = 0;
    size_t numOfInvertedElements = 0;
    for (vtkIdType i = 0; i < qualityArray->GetNumberOfTuples(); i++) {
        double val = qualityArray->GetValue(i);
        if (minValue > val) minValue = val;
        if (maxValue < val) maxValue = val;
        avgValue += val;

        if (val < minSJ) {
            numOfInvertedElements++;
            //badCellIds.push_back(i);
        } else if (val > 1) numOfInvertedElements++;
        //std::cout << "value " << i << " : " << val << std::endl;
    }
    avgValue /= qualityArray->GetNumberOfTuples();
    return numOfInvertedElements;
}

#include "verdict.h"
size_t Mesh::GetQualityVerdict(double& minValue, double& avgValue, const double minSJ/* = 0.0*/) {
    size_t numOfInvertedElements = 0;
    avgValue = 0.0;
    //badCellIds.clear();
    double minScaledJacobian = 1;
    //const std::vector<Vertex>& V = mesh.V;
    //const std::vector<Cell>& C = mesh.C;
    double coordinates[8][3];
    for (size_t i = 0; i < C.size(); i++) {
        const Cell& c = C[i];
        for (size_t j = 0; j < 8; j++) {
            const Vertex& v = V[c.Vids[j]];
            coordinates[j][0] = v.x;
            coordinates[j][1] = v.y;
            coordinates[j][2] = v.z;
        }
        double scaledJacobian = v_hex_scaled_jacobian(8, coordinates);
        if (scaledJacobian < minScaledJacobian) minScaledJacobian = scaledJacobian;
        //minScaledJacobian = minScaledJacobian < scaledJacobian ? minScaledJacobian : scaledJacobian;
        if (scaledJacobian < minSJ) {
            numOfInvertedElements++;
            //badCellIds.push_back(i);
        }
        avgValue += scaledJacobian;
    }
    minValue = minScaledJacobian;
    avgValue /= C.size();

    return numOfInvertedElements;
}
// preseverQuality, so we cannot be parallel
double Mesh::SmoothAndProjectSurface(const Mesh& mesh, size_t iters/* = 1*/, const SmoothMethod smoothMethod/* = LAPLACE_EDGE*/,
        const bool preserveSharpFeature/* = false*/, const bool treatSharpFeatureAsRegular/* = false*/, const bool treatCornerAsRegular/* = false*/,
        const bool preserveQuality/* = false*/)
{
    double targetMinSJ = 0;
    double targetAvgSJ = 0;
    //double targetMinSJ = GetMinScaledJacobian(targetAvgSJ);
    //MeshFileWriter writer(*this, "temp.vtk");
    //writer.WriteFile();
    //GetQuality("temp.vtk", targetMinSJ, targetAvgSJ);
    GetQualityVerdict(targetMinSJ, targetAvgSJ);
    //GetQuality(targetMinSJ, targetAvgSJ);
    std::cout << "targetMinSJ = " << targetMinSJ << " targetAvgSJ = " << targetAvgSJ << std::endl;
    double energy = 0;
    std::vector<Vertex> oldV = V;
    std::vector<Vertex> newV = V;
    while (iters-- != 0)
    {
//#pragma omp parallel for
        for (size_t i = 0; i < V.size(); i++) {
            const Vertex& v = V.at(i);
            if (!v.isBoundary)
                continue;
            if (v.type == REGULAR || (v.type == FEATURE && treatSharpFeatureAsRegular) || (v.type == CORNER && treatCornerAsRegular))
            {
                if (smoothMethod == LAPLACE_EDGE) {
                    const glm::dvec3 smoothedv = LapLace(v, treatSharpFeatureAsRegular, treatCornerAsRegular);
                    const glm::dvec3 newv = mesh.GetProjectLocation(smoothedv);
                    if (!preserveQuality) {
                        V.at(i) = newv;
                    }
                    else {
                        const glm::dvec3 oldv = V.at(i).xyz();
                        double oldMinSJ = 0;
                        double oldAvgSJ = 0;
                        GetQuality(V.at(i), oldMinSJ, oldAvgSJ);
                        double newMinSJ = 0;
                        double newAvgSJ = 0;
                        V.at(i) = newv;
                        GetQuality(V.at(i), newMinSJ, newAvgSJ);
                        if (newMinSJ < targetMinSJ || newAvgSJ < oldAvgSJ)
                            V.at(i) = oldv;
                        else{
                            const double distance = glm::length(newv- oldv);
                            energy += distance * distance;
                        }
                    }
                }
                else if (smoothMethod == LAPLACE_FACE_CENTER) {
                    glm::dvec3 sum(0.0, 0.0, 0.0);
                    int count = 0;
                    for (size_t j = 0; j < v.N_Fids.size(); j++) {
                        const Face& f = F.at(v.N_Fids[j]);
                        if (!f.isBoundary)
                            continue;
                        sum += GetFaceCenter(f);
                        count++;
                    }
                    newV.at(i).x = sum.x/count;
                    newV.at(i).y = sum.y/count;
                    newV.at(i).z = sum.z/count;
                }
            }
            else if (v.type == FEATURE) {
                if (preserveSharpFeature)
                    continue;
                if (smoothMethod == LAPLACE_EDGE) {
                    newV.at(i) = LapLace(v);
                }
                else if (smoothMethod == LAPLACE_FACE_CENTER) {
                    glm::dvec3 sum(0.0, 0.0, 0.0);
                    int count = 0;
                    for (size_t j = 0; j < v.N_Fids.size(); j++) {
                        const Face& f = F.at(v.N_Fids[j]);
                        if (!f.isBoundary)
                            continue;
                        sum += GetFaceCenter(f);
                        count++;
                    }
                    newV.at(i).x = sum.x/count;
                    newV.at(i).y = sum.y/count;
                    newV.at(i).z = sum.z/count;
                }
            }
        }
    }

    std::cout << "Energy = " << energy << std::endl;
    return energy;
}

double Mesh::ProjectSurface(const Mesh& mesh, size_t iters/* = 1*/, const SmoothMethod smoothMethod/* = LAPLACE_EDGE*/,
        const bool preserveSharpFeature/* = false*/, const bool treatSharpFeatureAsRegular/* = false*/, const bool treatCornerAsRegular/* = false*/,
        const bool preserveQuality/* = false*/)
{
    double targetAvgSJ = 0;
    double targetMinSJ = GetMinScaledJacobian(targetAvgSJ);
    std::cout << "targetMinSJ = " << targetMinSJ << " targetAvgSJ = " << targetAvgSJ << std::endl;
    double energy = 0;
    std::vector<Vertex> oldV = V;
    std::vector<Vertex> newV = V;
    while (iters-- != 0)
    {
//#pragma omp parallel for
        for (size_t i = 0; i < V.size(); i++) {
            const Vertex& v = V.at(i);
            if (!v.isBoundary)
                continue;
            if (v.type == REGULAR || (v.type == FEATURE && treatSharpFeatureAsRegular) || (v.type == CORNER && treatCornerAsRegular))
            {
                if (smoothMethod == LAPLACE_EDGE) {
                    //const glm::dvec3 smoothedv = LapLace(v, treatSharpFeatureAsRegular, treatCornerAsRegular);
                    const glm::dvec3 smoothedv = v.xyz();
                    const glm::dvec3 newv = mesh.GetProjectLocation(smoothedv);
                    if (!preserveQuality) {
                        V.at(i) = newv;
                    }
                    else {
                        const glm::dvec3 oldv = V.at(i).xyz();
                        double oldMinSJ = 0;
                        double oldAvgSJ = 0;
                        GetQuality(V.at(i), oldMinSJ, oldAvgSJ);
                        double newMinSJ = 0;
                        double newAvgSJ = 0;
                        V.at(i) = newv;
                        GetQuality(V.at(i), newMinSJ, newAvgSJ);
                        if (newMinSJ < targetMinSJ || newAvgSJ < oldAvgSJ)
                            V.at(i) = oldv;
                        else{
                            const double distance = glm::length(newv- oldv);
                            energy += distance * distance;
                        }
                    }
                }
                else if (smoothMethod == LAPLACE_FACE_CENTER) {
                    glm::dvec3 sum(0.0, 0.0, 0.0);
                    int count = 0;
                    for (size_t j = 0; j < v.N_Fids.size(); j++) {
                        const Face& f = F.at(v.N_Fids[j]);
                        if (!f.isBoundary)
                            continue;
                        sum += GetFaceCenter(f);
                        count++;
                    }
                    newV.at(i).x = sum.x/count;
                    newV.at(i).y = sum.y/count;
                    newV.at(i).z = sum.z/count;
                }
            }
            else if (v.type == FEATURE) {
                if (preserveSharpFeature)
                    continue;
                if (smoothMethod == LAPLACE_EDGE) {
                    newV.at(i) = LapLace(v);
                }
                else if (smoothMethod == LAPLACE_FACE_CENTER) {
                    glm::dvec3 sum(0.0, 0.0, 0.0);
                    int count = 0;
                    for (size_t j = 0; j < v.N_Fids.size(); j++) {
                        const Face& f = F.at(v.N_Fids[j]);
                        if (!f.isBoundary)
                            continue;
                        sum += GetFaceCenter(f);
                        count++;
                    }
                    newV.at(i).x = sum.x/count;
                    newV.at(i).y = sum.y/count;
                    newV.at(i).z = sum.z/count;
                }
            }
        }
    }

    std::cout << "Energy = " << energy << std::endl;
    return energy;
}

double Mesh::SmoothVolume(const Mesh& mesh, size_t iters/* = 1*/, const SmoothMethod smoothMethod/* = LAPLACE_EDGE*/,
        const bool preserveQuality/* = false*/)
{
    double targetAvgSJ = 0;
    double targetMinSJ = GetMinScaledJacobian(targetAvgSJ);
    std::cout << "targetMinSJ = " << targetMinSJ << " targetAvgSJ = " << targetAvgSJ << std::endl;
    double energy = 0;
    std::vector<Vertex> oldV = V;
    std::vector<Vertex> newV = V;
    while (iters-- != 0)
    {
//#pragma omp parallel for
        for (size_t i = 0; i < V.size(); i++) {
            const Vertex& v = V.at(i);
            if (v.isBoundary)
                continue;
            //if (v.type == REGULAR || (v.type == FEATURE && treatSharpFeatureAsRegular) || (v.type == CORNER && treatCornerAsRegular))
            {
                if (smoothMethod == LAPLACE_EDGE) {
                    const glm::dvec3 smoothedv = LapLace(v);
                    //const glm::dvec3 newv = mesh.GetProjectLocation(smoothedv);
                    const glm::dvec3 newv = smoothedv;
                    if (!preserveQuality) {
                        V.at(i) = newv;
                    }
                    else {
                        const glm::dvec3 oldv = V.at(i).xyz();
                        double oldMinSJ = 0;
                        double oldAvgSJ = 0;
                        GetQuality(V.at(i), oldMinSJ, oldAvgSJ);
                        double newMinSJ = 0;
                        double newAvgSJ = 0;
                        V.at(i) = newv;
                        GetQuality(V.at(i), newMinSJ, newAvgSJ);
                        if (newMinSJ < targetMinSJ || newAvgSJ < oldAvgSJ)
                            V.at(i) = oldv;
                        else{
                            const double distance = glm::length(newv- oldv);
                            energy += distance * distance;
                        }
                    }
                }
                else if (smoothMethod == LAPLACE_FACE_CENTER) {
                    glm::dvec3 sum(0.0, 0.0, 0.0);
                    int count = 0;
                    for (size_t j = 0; j < v.N_Fids.size(); j++) {
                        const Face& f = F.at(v.N_Fids[j]);
                        if (!f.isBoundary)
                            continue;
                        sum += GetFaceCenter(f);
                        count++;
                    }
                    newV.at(i).x = sum.x/count;
                    newV.at(i).y = sum.y/count;
                    newV.at(i).z = sum.z/count;
                }
            }
        }
    }

    std::cout << "Volum Energy = " << energy << std::endl;
    return energy;
}


const int V_T[8][4] =
{
    {0, 3, 4, 1},
    {1, 0, 5, 2},
    {2, 1, 6, 3},
    {3, 2, 7, 0},
    {4, 7, 5, 0},
    {5, 4, 6, 1},
    {6, 5, 7, 2},
    {7, 6, 4, 3}
};

static float cal_volume_Tet_real(float v0[3],float v1[3],float v2[3],float v3[3])
{
    float v1v0[3], v2v0[3], v3v0[3];
    for (int i = 0; i < 3; i++) {
        v1v0[i] = v1[i] - v0[i];
        v2v0[i] = v2[i] - v0[i];
        v3v0[i] = v3[i] - v0[i];
    }

    float norm1 = sqrt(v1v0[0] * v1v0[0] + v1v0[1] * v1v0[1] + v1v0[2] * v1v0[2]);
    float norm2 = sqrt(v2v0[0] * v2v0[0] + v2v0[1] * v2v0[1] + v2v0[2] * v2v0[2]);
    float norm3 = sqrt(v3v0[0] * v3v0[0] + v3v0[1] * v3v0[1] + v3v0[2] * v3v0[2]);

    float volume = v1v0[0] * (v2v0[1] * v3v0[2] - v2v0[2] * v3v0[1]) - v1v0[1] * (v2v0[0] * v3v0[2] - v2v0[2] * v3v0[0]) + v1v0[2] * (v2v0[0] * v3v0[1] - v2v0[1] * v3v0[0]);
    return volume;
}

bool JudgeDirection(const Mesh& mesh, const Cell& c) {
    float v[8][3];
    for (int i = 0; i < c.Vids.size(); i++) {
        v[i][0] = mesh.V.at(c.Vids.at(i)).x;
        v[i][1] = mesh.V.at(c.Vids.at(i)).y;
        v[i][2] = mesh.V.at(c.Vids.at(i)).z;
    }

    float VL[8];
    for (int i = 0; i < 8; i++) {
        const int* p = V_T[i];
        VL[i] = cal_volume_Tet_real(v[p[0]], v[p[1]], v[p[2]], v[p[3]]);
    }

    if (VL[0] + VL[1] + VL[2] + VL[3] + VL[4] + VL[5] + VL[6] + VL[7] < 0) return false;
    else return true;
}

static float _GetScaledJacobian(const glm::dvec3& i, const glm::dvec3& j, const glm::dvec3& k) {
    const glm::mat3x3 m(i, j, k);
    return glm::determinant(m);
}

const float Mesh::GetScaledJacobian(const Cell& c) const {
    Cell c1(c);
    if (!JudgeDirection(*this, c)) {
        std::swap(c1.Vids[0], c1.Vids[3]);
        std::swap(c1.Vids[1], c1.Vids[2]);
        std::swap(c1.Vids[4], c1.Vids[7]);
        std::swap(c1.Vids[5], c1.Vids[6]);
    }
    float minScaledJacobian = 1;
    for (int n = 0; n < 7; n++) {
        const Vertex& o = V.at(c1.Vids.at(n));

        const Vertex& i = V.at(c1.Vids.at(HexPoint_Points[n][0]));
        const Vertex& j = V.at(c1.Vids.at(HexPoint_Points[n][1]));
        const Vertex& k = V.at(c1.Vids.at(HexPoint_Points[n][2]));

        const glm::dvec3 ei(i.x - o.x, i.y - o.y, i.z - o.z);
        const glm::dvec3 ej(j.x - o.x, j.y - o.y, j.z - o.z);
        const glm::dvec3 ek(k.x - o.x, k.y - o.y, k.z - o.z);

        const float length_i = glm::length(ei);
        const float length_j = glm::length(ej);
        const float length_k = glm::length(ek);

        const glm::dvec3 ni(ei.x / length_i, ei.y / length_i, ei.z / length_i);
        const glm::dvec3 nj(ej.x / length_j, ej.y / length_j, ej.z / length_j);
        const glm::dvec3 nk(ek.x / length_k, ek.y / length_k, ek.z / length_k);

        float scaledJacobian = _GetScaledJacobian(ni, nj, nk);
        minScaledJacobian = minScaledJacobian < scaledJacobian ? minScaledJacobian : scaledJacobian;
    }

    return minScaledJacobian;
}
//void Mesh::GetQuality(const Vertex& v, double& minSJ, double& avgSJ)
//{
//    double min = 1.0;
//    double sum = 0;
//    for (size_t i = 0; i < v.N_Cids.size(); i++) {
//        const Cell& cell = C.at(v.N_Cids.at(i));
//        const float scaledJacobian = GetScaledJacobian(cell);
//        if (scaledJacobian < min)
//            min = scaledJacobian;
//        sum += scaledJacobian;
//    }
//    minSJ = min;
//    avgSJ = sum / v.N_Cids.size();
//}

void Mesh::GetQuality(const Vertex& v, double& minSJ, double& avgSJ) {
    double minScaledJacobian = 1;
    double coordinates[8][3];
    for (size_t i = 0; i < v.N_Cids.size(); i++) {
        const Cell& c = C.at(v.N_Cids.at(i));
        for (size_t j = 0; j < 8; j++) {
            const Vertex& v = V[c.Vids[j]];
            coordinates[j][0] = v.x;
            coordinates[j][1] = v.y;
            coordinates[j][2] = v.z;
        }
        double scaledJacobian = v_hex_scaled_jacobian(8, coordinates);
        if (scaledJacobian < minScaledJacobian) minScaledJacobian = scaledJacobian;
        avgSJ += scaledJacobian;
    }
    minSJ = minScaledJacobian;
    avgSJ /= v.N_Cids.size();

//    std::vector<size_t> cellIds(v.N_Cids.size());
//    for (size_t i = 0; i < v.N_Cids.size(); i++) {
//        const Cell& cell = C.at(v.N_Cids.at(i));
//        cellIds[i] = cell.id;
//    }
    //OutputBadCells(cellIds, "temp.vtk");
    //GetQuality("temp.vtk", minSJ, avgSJ);
}

void Mesh::OutputBadCells(const std::vector<size_t>& badCellIds, const char* filename) {
    std::vector<Cell> cells(badCellIds.size());
    for (size_t i = 0; i < badCellIds.size(); i++)
        cells.at(i) = C.at(badCellIds.at(i));
    MeshFileWriter writer(V, cells, filename);
    writer.FixMesh();
    writer.WriteFile();
}

glm::dvec3 Mesh::LapLace(const Vertex& v, const bool treatSharpFeatureAsRegular/* = false*/, const bool treatCornerAsRegular/* = false*/) {
    if (v.type == REGULAR || (v.type == FEATURE && treatSharpFeatureAsRegular) || (v.type == CORNER && treatCornerAsRegular)) {
        glm::dvec3 sum(0.0, 0.0, 0.0);
        int count = 0;
        for (size_t j = 0; j < v.N_Vids.size(); j++) {
            const Vertex& n_v = V.at(v.N_Vids[j]);
            if (!n_v.isBoundary) continue;
            sum += n_v.xyz();
            count++;
        }
        return glm::dvec3(sum.x / count, sum.y / count, sum.z / count);
    } else if (v.type == FEATURE) {
        glm::dvec3 sum(0.0, 0.0, 0.0);
        int count = 0;
        for (size_t j = 0; j < v.N_Vids.size(); j++) {
            const Vertex& n_v = V.at(v.N_Vids[j]);
            if (!n_v.isBoundary) continue;
            if (n_v.type != FEATURE && n_v.type != CORNER) continue;
            sum += n_v.xyz();
            count++;
        }
        if (count == 1) {
            int c = 0;
            sum = glm::dvec3(0.0, 0.0, 0.0);
            for (size_t j = 0; j < v.N_Vids.size(); j++) {
                const Vertex& n_v = V.at(v.N_Vids[j]);
                if (!n_v.isBoundary) continue;
                sum += n_v.xyz();
                c++;
            }
            return glm::dvec3(sum.x / c, sum.y / c, sum.z / c);
        } else if (count == 2) {
            return glm::dvec3(sum.x / count, sum.y / count, sum.z / count);
        } else {
            //std::cerr << "\033[1;31mERROR\033[0m in smoothing Vertex of Sharp Feature!!" << std::endl;
            return v.xyz();
        }
    } else if (v.type == CORNER) return v.xyz();
    else {
        glm::dvec3 sum(0.0, 0.0, 0.0);
        int count = 0;
        for (size_t j = 0; j < v.N_Vids.size(); j++) {
            const Vertex& n_v = V.at(v.N_Vids[j]);
            //if (!n_v.isBoundary)
            //    continue;
            sum += n_v.xyz();
            count++;
        }
        return glm::dvec3(sum.x / count, sum.y / count, sum.z / count);
    }
}

//glm::dvec3 Mesh::LapLace(const Vertex& v)
//{
//    if (v.type == REGULAR || v.type == FEATURE) {
//        glm::dvec3 sum(0.0, 0.0, 0.0);
//        int count = 0;
//        for (size_t j = 0; j < v.N_Vids.size(); j++) {
//            const Vertex& n_v = V.at(v.N_Vids[j]);
//            if (!n_v.isBoundary)
//                continue;
//            sum += n_v.xyz();
//            count++;
//        }
//        return glm::dvec3(sum.x/count, sum.y/count, sum.z/count);
//    }
//    else if (v.type == FEATURE) {
//        glm::dvec3 sum(0.0, 0.0, 0.0);
//        int count = 0;
//        for (size_t j = 0; j < v.N_Vids.size(); j++) {
//            const Vertex& n_v = V.at(v.N_Vids[j]);
//            if (!n_v.isBoundary)
//                continue;
//            if (n_v.type != FEATURE && n_v.type != CORNER)
//                continue;
//            sum += n_v.xyz();
//            count++;
//        }
//        if (count == 1) {
//            int c = 0;
//            sum = glm::dvec3 (0.0, 0.0, 0.0);
//            for (size_t j = 0; j < v.N_Vids.size(); j++) {
//                const Vertex& n_v = V.at(v.N_Vids[j]);
//                if (!n_v.isBoundary)
//                    continue;
//                sum += n_v.xyz();
//                c++;
//            }
//            return glm::dvec3(sum.x/c, sum.y/c, sum.z/c);
//        }
//        else if (count == 2) {
//            return glm::dvec3(sum.x/count, sum.y/count, sum.z/count);
//        }
//        else {
//            std::cout << "\033[1;31mERROR\033[0m in smoothing Vertex of Sharp Feature!!" << std::endl;
//        }
//    }
//    else
//        return v.xyz();
//}

static void set_cross(const std::vector<size_t>& set1, const std::vector<size_t>& set2, std::vector<size_t> &result_set) {
    result_set.clear();
    for (size_t i = 0; i < set1.size(); i++) {
        bool inside = false;
        for (size_t j = 0; j < set2.size(); j++)
            if (set2[j] == set1[i]) {
                inside = true;
                break;
            }
        if (inside) result_set.push_back(set1[i]);
    }
}

std::vector<size_t> unique(const std::vector<size_t>& s) {
	std::vector<size_t> res;
	std::set<size_t> ss(s.begin(), s.end());
	std::copy(ss.begin(), ss.end(), std::back_inserter(res));
	return res;
}

void Mesh::BuildF()
{
    if (m_cellType == TRIANGLE || m_cellType == QUAD || m_cellType == POLYGON) {
        F.resize(C.size());
        for (size_t i = 0; i < C.size(); i++) {
            F[i].Vids = C[i].Vids;
            F[i].id = C[i].id;
            //F[i].isBoundary = true;
        }
        if (m_cellType == TRIANGLE) {
            for (size_t i = 0; i < F.size(); i++) {
                F[i].Eids.resize(3);
                for (size_t j = 0; j < 3; j++) {
                    bool found = false;
                    for (size_t k = 0; k < V[F[i].Vids[j]].N_Eids.size(); k++) {
                        size_t ne = V[F[i].Vids[j]].N_Eids[k];
                        for (size_t m = 0; m < V[F[i].Vids[(j + 1) % 3]].N_Eids.size(); m++) {
                            if (ne == V[F[i].Vids[(j + 1) % 3]].N_Eids[m]) {
                                E[ne].N_Fids.push_back(i);
                                F[i].Eids[j] = ne;
                                found = true;
                                break;
                            }
                        }
                        if (found) break;
                    }
                }
            }
            ////////////////////////////////////////////////
            std::vector<bool> Vs_flags(V.size(), false);
            for (int i = 0; i < F.size(); i++) {
                for (int j = 0; j < 3; j++)
                    Vs_flags[F[i].Vids[j]] = true;
                for (int j = 0; j < F[i].N_Cids.size(); j++) {
                    int nhid = F[i].N_Cids[j];
                    for (int k = 0; k < 4; k++) {
                        bool have_true = false;
                        for (int m = 0; m < 3; m++)
                            if (Vs_flags[F[C[nhid].Fids[k]].Vids[m]]) {
                                have_true = true;
                                break;
                            }
                        if (!have_true) {
                            F[i].N_Fids.push_back(C[nhid].Fids[k]);
                            break;
                        }
                    }
                }
            }
        }
        else if (m_cellType == QUAD) {
            for (size_t i = 0; i < F.size(); i++) {
                F[i].Eids.resize(4);
                for (size_t j = 0; j < 4; j++) {
                    bool found = false;
                    for (size_t k = 0; k < V[F[i].Vids[j]].N_Eids.size(); k++) {
                        size_t ne = V[F[i].Vids[j]].N_Eids[k];
                        for (size_t m = 0; m < V[F[i].Vids[(j + 1) % 4]].N_Eids.size(); m++) {
                            if (ne == V[F[i].Vids[(j + 1) % 4]].N_Eids[m]) {
                                E[ne].N_Fids.push_back(i);
                                F[i].Eids[j] = ne;
                                found = true;
                                break;
                            }
                        }
                        if (found) break;
                    }
                }
            }
            ////////////////////////////////////////////////
            std::vector<bool> Vs_flags(V.size(), false);
            for (int i = 0; i < F.size(); i++) {
                for (int j = 0; j < 4; j++)
                    Vs_flags[F[i].Vids[j]] = true;
                for (int j = 0; j < F[i].N_Cids.size(); j++) {
                    int nhid = F[i].N_Cids[j];
                    for (int k = 0; k < 4; k++) {
                        bool have_true = false;
                        for (int m = 0; m < 4; m++)
                            if (Vs_flags[F[C[nhid].Fids[k]].Vids[m]]) {
                                have_true = true;
                                break;
                            }
                        if (!have_true) {
                            F[i].N_Fids.push_back(C[nhid].Fids[k]);
                            break;
                        }
                    }
                }
            }
        }
        else if (m_cellType == POLYGON) {
            for (size_t i = 0; i < F.size(); i++) {
                auto& f = F[i];
                f.Eids.resize(f.Vids.size());
                if (f.Vids.size() == 4) {
                    for (size_t j = 0; j < 4; j++) {
                        bool found = false;
                        for (size_t k = 0; k < V[F[i].Vids[j]].N_Eids.size(); k++) {
                            size_t ne = V[F[i].Vids[j]].N_Eids[k];
                            for (size_t m = 0; m < V[F[i].Vids[(j + 1) % 4]].N_Eids.size(); m++) {
                                if (ne == V[F[i].Vids[(j + 1) % 4]].N_Eids[m]) {
                                    E[ne].N_Fids.push_back(i);
                                    F[i].Eids[j] = ne;
                                    found = true;
                                    break;
                                }
                            }
                            if (found) break;
                        }
                    }
                } else if (f.Vids.size() == 3) {
                    for (size_t j = 0; j < 3; j++) {
                        bool found = false;
                        for (size_t k = 0; k < V[F[i].Vids[j]].N_Eids.size(); k++) {
                            size_t ne = V[F[i].Vids[j]].N_Eids[k];
                            for (size_t m = 0; m < V[F[i].Vids[(j + 1) % 3]].N_Eids.size(); m++) {
                                if (ne == V[F[i].Vids[(j + 1) % 3]].N_Eids[m]) {
                                    E[ne].N_Fids.push_back(i);
                                    F[i].Eids[j] = ne;
                                    found = true;
                                    break;
                                }
                            }
                            if (found) break;
                        }
                    }
                } else {
					for (size_t j = 0; j < f.Vids.size(); j++) {
						bool found = false;
						for (size_t k = 0; k < V[f.Vids[j]].N_Eids.size(); k++) {
							size_t ne = V[f.Vids[j]].N_Eids[k];
							for (size_t m = 0; m < V[f.Vids[(j + 1) % f.Vids.size()]].N_Eids.size(); m++) {
								if (ne == V[f.Vids[(j + 1) % f.Vids.size()]].N_Eids[m]) {
									E[ne].N_Fids.push_back(i);
									f.Eids[j] = ne;
									found = true;
									break;
								}
							}
							if (found) break;
						}
					}
				}
            }
            ////////////////////////////////////////////////
            std::vector<bool> Vs_flags(V.size(), false);
            for (int i = 0; i < F.size(); i++) {
                const auto& f = F[i];
                if (f.Vids.size() == 4) {
                    for (int j = 0; j < 4; j++)
                        Vs_flags[F[i].Vids[j]] = true;
                    for (int j = 0; j < F[i].N_Cids.size(); j++) {
                        int nhid = F[i].N_Cids[j];
                        for (int k = 0; k < 4; k++) {
                            bool have_true = false;
                            for (int m = 0; m < 4; m++)
                                if (Vs_flags[F[C[nhid].Fids[k]].Vids[m]]) {
                                    have_true = true;
                                    break;
                                }
                            if (!have_true) {
                                F[i].N_Fids.push_back(C[nhid].Fids[k]);
                                break;
                            }
                        }
                    }
                } else if (f.Vids.size() == 3) {
                    for (int j = 0; j < 3; j++)
                        Vs_flags[F[i].Vids[j]] = true;
                    for (int j = 0; j < F[i].N_Cids.size(); j++) {
                        int nhid = F[i].N_Cids[j];
                        for (int k = 0; k < 4; k++) {
                            bool have_true = false;
                            for (int m = 0; m < 3; m++)
                                if (Vs_flags[F[C[nhid].Fids[k]].Vids[m]]) {
                                    have_true = true;
                                    break;
                                }
                            if (!have_true) {
                                F[i].N_Fids.push_back(C[nhid].Fids[k]);
                                break;
                            }
                        }
                    }
                }
            }
        }
    } else if (m_cellType == HEXAHEDRA) {
        F.reserve(C.size() * 6);
        int F_N = 0;
        for (size_t i = 0; i < C.size(); i++) {
            std::vector<Face> hf(6, Face(4, 4));
            for (size_t j = 0; j < 6; j++)
                for (size_t k = 0; k < 4; k++)
                    hf[j].Vids[k] = C[i].Vids[HexFaces[j][k]];

            for (size_t j = 0; j < 6; j++) {
                bool have = false;
                for (size_t m = 0; m < 4; m++) {
                    for (size_t n = 0; n < V[hf[j].Vids[m]].N_Fids.size(); n++) {
                        const size_t F_id = V[hf[j].Vids[m]].N_Fids[n];
                        bool all_have = true;
                        for (size_t p = 0; p < 4; p++) {
                            bool exist_v = false;
                            for (size_t q = 0; q < 4; q++)
                                if (hf[j].Vids[p] == F[F_id].Vids[q])
                                    exist_v = true;
                            if (!exist_v)
                                all_have = false;
                        }
                        if (all_have) {
                            have = true;
                            F[F_id].N_Cids.push_back(i);
                        }
                    }
                }
                if (!have) {
                    hf[j].id = F_N++;
                    for (size_t k = 0; k < 4; k++) {
                        size_t id1 = hf[j].Vids[k];
                        size_t id2 = hf[j].Vids[(k + 1) % 4];
                        bool found = false;
                        for (size_t m = 0; m < V[id1].N_Eids.size(); m++) {
                            int edge1 = V[id1].N_Eids[m];
                            for (size_t n = 0; n < V[id2].N_Eids.size(); n++) {
                                size_t edge2 = V[id2].N_Eids[n];
                                if (edge1 == edge2) {
                                    hf[j].Eids[k] = edge1;
                                    found = true;
                                }
                                if (found)
                                    break;
                            }
                            if (found)
                                break;
                        }
                    }
                    F.push_back(hf[j]);
                    V[hf[j].Vids[0]].N_Fids.push_back(hf[j].id);
                    V[hf[j].Vids[1]].N_Fids.push_back(hf[j].id);
                    V[hf[j].Vids[2]].N_Fids.push_back(hf[j].id);
                    V[hf[j].Vids[3]].N_Fids.push_back(hf[j].id);

                    F[F.size() - 1].N_Cids.push_back(i);
                }
            }
        }

        for (size_t i = 0; i < F.size(); i++) {
            std::vector<size_t> N_Cids = F[i].N_Cids;
            F[i].N_Cids.clear();
            for (size_t j = 0; j < N_Cids.size(); j++) {
                bool already = false;
                for (size_t k = 0; k < F[i].N_Cids.size(); k++) {
                    if (N_Cids[j] == F[i].N_Cids[k]) {
                        already = true;
                        break;
                    }
                }
                if (!already)
                    F[i].N_Cids.push_back(N_Cids[j]);
            }
        }

        for (size_t i = 0; i < F.size(); i++) {
            for (size_t j = 0; j < 4; j++) {
                bool found = false;
                for (size_t k = 0; k < V[F[i].Vids[j]].N_Eids.size(); k++) {
                    size_t ne = V[F[i].Vids[j]].N_Eids[k];
                    for (size_t m = 0; m < V[F[i].Vids[(j + 1) % 4]].N_Eids.size(); m++) {
                        if (ne == V[F[i].Vids[(j + 1) % 4]].N_Eids[m]) {
                            E[ne].N_Fids.push_back(i);
                            F[i].Eids[j] = ne;
                            found = true;
                            break;
                        }
                    }
                    if (found)
                        break;
                }
            }
        }

        for (size_t i = 0; i < F.size(); i++)
            for (size_t j = 0; j < F[i].N_Cids.size(); j++)
                C[F[i].N_Cids[j]].N_Fids.push_back(i);

        for (size_t i = 0; i < C.size(); i++) {
            C[i].Fids.resize(6);
            std::vector<size_t> f_ids = C[i].N_Fids;
            C[i].N_Fids.clear();
            for (size_t j = 0; j < f_ids.size(); j++) {
                bool havesame = false;
                for (size_t k = j + 1; k < f_ids.size(); k++)
                    if (f_ids[j] == f_ids[k])
                        havesame = true;
                if (!havesame) {
                    C[i].N_Fids.push_back(f_ids[j]);
                    C[i].Fids[C[i].N_Fids.size() - 1] = F[f_ids[j]].id;
                }
            }
        }
        ////////////////////////////////////////////////
        std::vector<bool> Vs_flags(V.size(), false);
        for (int i = 0; i < F.size(); i++) {
            for (int j = 0; j < 4; j++)
                Vs_flags[F[i].Vids[j]] = true;
            for (int j = 0; j < F[i].N_Cids.size(); j++) {
                int nhid = F[i].N_Cids[j];
                for (int k = 0; k < 6; k++) {
                    bool have_true = false;
                    for (int m = 0; m < 4; m++)
                        if (Vs_flags[F[C[nhid].Fids[k]].Vids[m]]) {
                            have_true = true;
                            break;
                        }
                    if (!have_true) {
                        F[i].N_Fids.push_back(C[nhid].Fids[k]);
                        break;
                    }
                }
            }

            std::vector<size_t> N_Ortho_4Vs1(4), N_Ortho_4Vs2(4);
            for (size_t k = 0; k < 4; k++) {
                size_t fvid = F[F[i].N_Fids[0]].Vids[k];
                N_Ortho_4Vs1[k] = fvid;
            }
            F[i].N_Ortho_4Vids.push_back(N_Ortho_4Vs1);
            for (size_t j = 0; j < F[i].N_Fids.size(); j++) {
                if (j == 1)
                    for (size_t k = 0; k < 4; k++)
                        Vs_flags[F[F[i].N_Fids[1]].Vids[k]] = true;
                for (size_t k = 0; k < 4; k++) {
                    int fvid = N_Ortho_4Vs1[k];
                    for (size_t m = 0; m < V[fvid].N_Vids.size(); m++) {
                        if (Vs_flags[V[fvid].N_Vids[m]]) {
                            N_Ortho_4Vs2[k] = V[fvid].N_Vids[m];
                            break;
                        }
                    }
                }
                F[i].N_Ortho_4Vids.push_back(N_Ortho_4Vs2);
                N_Ortho_4Vs1 = N_Ortho_4Vs2;
                for (size_t k = 0; k < 4; k++)
                    Vs_flags[N_Ortho_4Vs1[k]] = false;
            }

            std::vector<size_t> N_4Eids(4);
            for (size_t j = 1; j < F[i].N_Ortho_4Vids.size(); j++) {
                for (size_t k = 0; k < 4; k++) {
                    std::vector<size_t> sharedEids;
                    set_cross(V[F[i].N_Ortho_4Vids[j - 1][k]].N_Eids, V[F[i].N_Ortho_4Vids[j][k]].N_Eids, sharedEids);
                    N_4Eids[k] = sharedEids[0];
                }
                F[i].N_Ortho_4Eids.push_back(N_4Eids);
            }
        }
    } else if (m_cellType == POLYHEDRA) {
		F.reserve(C.size() * 6);
		int F_N = 0;
		for (auto& c : C) {
			std::vector<Face> faces;
			Face temp_face;
			auto iter = cell_faces.find(c.cellType);
			auto iterEdge = cell_edges.find(c.cellType);
			for (auto& face : iter->second) {
				temp_face.Vids.resize(face.size());
				for (auto i = 0; i < face.size(); ++i)
					temp_face.Vids[i] = c.Vids[face[i]];
				temp_face.Eids.resize(face.size());
				faces.push_back(temp_face);
			}

			for (auto& f : faces) {
				bool have = false;
				for (auto fvid : f.Vids) {
					for (auto F_id : V[fvid].N_Fids) {
						bool all_have = true;
						for (size_t p : f.Vids) {
							bool exist_v = false;
							for (size_t q : F[F_id].Vids)
								if (p == q) exist_v = true;
							if (!exist_v) all_have = false;
						}
						if (all_have) {
							have = true;
							F[F_id].N_Cids.push_back(c.id);
						}
					}
				}
				if (!have) {
					f.id = F_N++;
					for (size_t k = 0; k < f.Vids.size(); k++) {
						size_t id1 = f.Vids[k];
						size_t id2 = f.Vids[(k + 1) % f.Vids.size()];
						bool found = false;
						for (auto edge1 : V[id1].N_Eids) {
							for (auto edge2 : V[id2].N_Eids) {
								if (edge1 == edge2) {
									f.Eids[k] = edge1;
									found = true;
									break;
								}
							}
							if (found)	break;
						}
					}
					F.push_back(f);
					for (auto fvid : f.Vids)
						V[fvid].N_Fids.push_back(f.id);
					F[F.size() - 1].N_Cids.push_back(c.id);
				}
			}
		}

		for (auto& f : F)
			f.N_Cids = unique(f.N_Cids);

		for (auto& f : F) {
			for (size_t j = 0; j < f.Vids.size(); j++) {
				bool found = false;
				for (auto neid : V[f.Vids[j]].N_Eids) {
					for (auto neid1 : V[f.Vids[(j + 1) % f.Vids.size()]].N_Eids) {
						if (neid == neid1) {
							E[neid].N_Fids.push_back(f.id);
							f.Eids[j] = neid;
							found = true;
							break;
						}
					}
					if (found) break;
				}
			}
		}

		for (auto& f : F)
			for (auto cid : f.N_Cids)
				C[cid].N_Fids.push_back(f.id);

		for (auto& c : C) {
			c.Fids.reserve(6);
			std::vector<size_t> f_ids = c.N_Fids;
			c.N_Fids.clear();
			for (size_t j = 0; j < f_ids.size(); j++) {
				bool havesame = false;
				for (size_t k = j + 1; k < f_ids.size(); k++)
					if (f_ids[j] == f_ids[k])
						havesame = true;
				if (!havesame) {
					c.N_Fids.push_back(f_ids[j]);
					c.Fids.push_back(F[f_ids[j]].id);
				}
			}
		}
		////////////////////////////////////////////////
		std::vector<bool> Vs_flags(V.size(), false);
		for (auto& f : F) {
			for (auto fvid : f.Vids)
				Vs_flags[fvid] = true;
			for (auto ncid : f.N_Cids) {
				auto& nc = C.at(ncid);
				for (auto ncfid : nc.Fids) {
					auto& ncf = F.at(ncfid);
					bool have_true = false;
					for (auto ncfvid : ncf.Vids)
						if (Vs_flags[ncfvid]) {
							have_true = true;
							break;
						}
					if (!have_true) {
						f.N_Fids.push_back(ncfid);
						break;
					}
				}
			}
		}
	} else if (m_cellType == TETRAHEDRA) {
        F.reserve(C.size() * 4);
        int F_N = 0;
        for (size_t i = 0; i < C.size(); i++) {
            std::vector<Face> hf(4, Face(3, 3));
            for (size_t j = 0; j < 4; j++)
                for (size_t k = 0; k < 3; k++)
                    hf[j].Vids[k] = C[i].Vids[TetFaces[j][k]];

            for (size_t j = 0; j < 4; j++) {
                bool have = false;
                for (size_t m = 0; m < 3; m++) {
                    for (size_t n = 0; n < V[hf[j].Vids[m]].N_Fids.size(); n++) {
                        const size_t F_id = V[hf[j].Vids[m]].N_Fids[n];
                        bool all_have = true;
                        for (size_t p = 0; p < 3; p++) {
                            bool exist_v = false;
                            for (size_t q = 0; q < 3; q++)
                                if (hf[j].Vids[p] == F[F_id].Vids[q])
                                    exist_v = true;
                            if (!exist_v)
                                all_have = false;
                        }
                        if (all_have) {
                            have = true;
                            F[F_id].N_Cids.push_back(i);
                        }
                    }
                }
                if (!have) {
                    hf[j].id = F_N++;
                    for (size_t k = 0; k < 3; k++) {
                        size_t id1 = hf[j].Vids[k];
                        size_t id2 = hf[j].Vids[(k + 1) % 3];
                        bool found = false;
                        for (size_t m = 0; m < V[id1].N_Eids.size(); m++) {
                            int edge1 = V[id1].N_Eids[m];
                            for (size_t n = 0; n < V[id2].N_Eids.size(); n++) {
                                size_t edge2 = V[id2].N_Eids[n];
                                if (edge1 == edge2) {
                                    hf[j].Eids[k] = edge1;
                                    found = true;
                                }
                                if (found)
                                    break;
                            }
                            if (found)
                                break;
                        }
                    }
                    F.push_back(hf[j]);
                    V[hf[j].Vids[0]].N_Fids.push_back(hf[j].id);
                    V[hf[j].Vids[1]].N_Fids.push_back(hf[j].id);
                    V[hf[j].Vids[2]].N_Fids.push_back(hf[j].id);

                    F[F.size() - 1].N_Cids.push_back(i);
                }
            }
        }

        for (size_t i = 0; i < F.size(); i++) {
            std::vector<size_t> N_Cids = F[i].N_Cids;
            F[i].N_Cids.clear();
            for (size_t j = 0; j < N_Cids.size(); j++) {
                bool already = false;
                for (size_t k = 0; k < F[i].N_Cids.size(); k++) {
                    if (N_Cids[j] == F[i].N_Cids[k]) {
                        already = true;
                        break;
                    }
                }
                if (!already)
                    F[i].N_Cids.push_back(N_Cids[j]);
            }
        }

        for (size_t i = 0; i < F.size(); i++) {
            for (size_t j = 0; j < 3; j++) {
                bool found = false;
                for (size_t k = 0; k < V[F[i].Vids[j]].N_Eids.size(); k++) {
                    size_t ne = V[F[i].Vids[j]].N_Eids[k];
                    for (size_t m = 0; m < V[F[i].Vids[(j + 1) % 3]].N_Eids.size(); m++) {
                        if (ne == V[F[i].Vids[(j + 1) % 3]].N_Eids[m]) {
                            E[ne].N_Fids.push_back(i);
                            F[i].Eids[j] = ne;
                            found = true;
                            break;
                        }
                    }
                    if (found)
                        break;
                }
            }
        }

        for (size_t i = 0; i < F.size(); i++)
            for (size_t j = 0; j < F[i].N_Cids.size(); j++)
                C[F[i].N_Cids[j]].N_Fids.push_back(i);

        for (size_t i = 0; i < C.size(); i++) {
            C[i].Fids.resize(4);
            std::vector<size_t> f_ids = C[i].N_Fids;
            C[i].N_Fids.clear();
            for (size_t j = 0; j < f_ids.size(); j++) {
                bool havesame = false;
                for (size_t k = j + 1; k < f_ids.size(); k++)
                    if (f_ids[j] == f_ids[k])
                        havesame = true;
                if (!havesame) {
                    C[i].N_Fids.push_back(f_ids[j]);
                    C[i].Fids[C[i].N_Fids.size() - 1] = F[f_ids[j]].id;
                }
            }
        }
        ////////////////////////////////////////////////
        std::vector<bool> Vs_flags(V.size(), false);
        for (int i = 0; i < F.size(); i++) {
            for (int j = 0; j < 3; j++)
                Vs_flags[F[i].Vids[j]] = true;
            for (int j = 0; j < F[i].N_Cids.size(); j++) {
                int nhid = F[i].N_Cids[j];
                for (int k = 0; k < 4; k++) {
                    bool have_true = false;
                    for (int m = 0; m < 3; m++)
                        if (Vs_flags[F[C[nhid].Fids[k]].Vids[m]]) {
                            have_true = true;
                            break;
                        }
                    if (!have_true) {
                        F[i].N_Fids.push_back(C[nhid].Fids[k]);
                        break;
                    }
                }
            }
        }
    } else if (m_cellType == TETRAHEDRA) {
        F.reserve(C.size() * 4);
        int F_N = 0;
        for (size_t i = 0; i < C.size(); i++) {
            std::vector<Face> hf(4, Face(3, 3));
            for (size_t j = 0; j < 4; j++)
                for (size_t k = 0; k < 3; k++)
                    hf[j].Vids[k] = C[i].Vids[TetFaces[j][k]];

            for (size_t j = 0; j < 4; j++) {
                bool have = false;
                for (size_t m = 0; m < 3; m++) {
                    for (size_t n = 0; n < V[hf[j].Vids[m]].N_Fids.size(); n++) {
                        const size_t F_id = V[hf[j].Vids[m]].N_Fids[n];
                        bool all_have = true;
                        for (size_t p = 0; p < 3; p++) {
                            bool exist_v = false;
                            for (size_t q = 0; q < 3; q++)
                                if (hf[j].Vids[p] == F[F_id].Vids[q])
                                    exist_v = true;
                            if (!exist_v)
                                all_have = false;
                        }
                        if (all_have) {
                            have = true;
                            F[F_id].N_Cids.push_back(i);
                        }
                    }
                }
                if (!have) {
                    hf[j].id = F_N++;
                    for (size_t k = 0; k < 3; k++) {
                        size_t id1 = hf[j].Vids[k];
                        size_t id2 = hf[j].Vids[(k + 1) % 3];
                        bool found = false;
                        for (size_t m = 0; m < V[id1].N_Eids.size(); m++) {
                            int edge1 = V[id1].N_Eids[m];
                            for (size_t n = 0; n < V[id2].N_Eids.size(); n++) {
                                size_t edge2 = V[id2].N_Eids[n];
                                if (edge1 == edge2) {
                                    hf[j].Eids[k] = edge1;
                                    found = true;
                                }
                                if (found)
                                    break;
                            }
                            if (found)
                                break;
                        }
                    }
                    F.push_back(hf[j]);
                    V[hf[j].Vids[0]].N_Fids.push_back(hf[j].id);
                    V[hf[j].Vids[1]].N_Fids.push_back(hf[j].id);
                    V[hf[j].Vids[2]].N_Fids.push_back(hf[j].id);

                    F[F.size() - 1].N_Cids.push_back(i);
                }
            }
        }

        for (size_t i = 0; i < F.size(); i++) {
            std::vector<size_t> N_Cids = F[i].N_Cids;
            F[i].N_Cids.clear();
            for (size_t j = 0; j < N_Cids.size(); j++) {
                bool already = false;
                for (size_t k = 0; k < F[i].N_Cids.size(); k++) {
                    if (N_Cids[j] == F[i].N_Cids[k]) {
                        already = true;
                        break;
                    }
                }
                if (!already)
                    F[i].N_Cids.push_back(N_Cids[j]);
            }
        }

        for (size_t i = 0; i < F.size(); i++) {
            for (size_t j = 0; j < 3; j++) {
                bool found = false;
                for (size_t k = 0; k < V[F[i].Vids[j]].N_Eids.size(); k++) {
                    size_t ne = V[F[i].Vids[j]].N_Eids[k];
                    for (size_t m = 0; m < V[F[i].Vids[(j + 1) % 3]].N_Eids.size(); m++) {
                        if (ne == V[F[i].Vids[(j + 1) % 3]].N_Eids[m]) {
                            E[ne].N_Fids.push_back(i);
                            F[i].Eids[j] = ne;
                            found = true;
                            break;
                        }
                    }
                    if (found)
                        break;
                }
            }
        }

        for (size_t i = 0; i < F.size(); i++)
            for (size_t j = 0; j < F[i].N_Cids.size(); j++)
                C[F[i].N_Cids[j]].N_Fids.push_back(i);

        for (size_t i = 0; i < C.size(); i++) {
            C[i].Fids.resize(4);
            std::vector<size_t> f_ids = C[i].N_Fids;
            C[i].N_Fids.clear();
            for (size_t j = 0; j < f_ids.size(); j++) {
                bool havesame = false;
                for (size_t k = j + 1; k < f_ids.size(); k++)
                    if (f_ids[j] == f_ids[k])
                        havesame = true;
                if (!havesame) {
                    C[i].N_Fids.push_back(f_ids[j]);
                    C[i].Fids[C[i].N_Fids.size() - 1] = F[f_ids[j]].id;
                }
            }
        }
        ////////////////////////////////////////////////
        std::vector<bool> Vs_flags(V.size(), false);
        for (int i = 0; i < F.size(); i++) {
            for (int j = 0; j < 3; j++)
                Vs_flags[F[i].Vids[j]] = true;
            for (int j = 0; j < F[i].N_Cids.size(); j++) {
                int nhid = F[i].N_Cids[j];
                for (int k = 0; k < 4; k++) {
                    bool have_true = false;
                    for (int m = 0; m < 3; m++)
                        if (Vs_flags[F[C[nhid].Fids[k]].Vids[m]]) {
                            have_true = true;
                            break;
                        }
                    if (!have_true) {
                        F[i].N_Fids.push_back(C[nhid].Fids[k]);
                        break;
                    }
                }
            }
        }
    }
}

void Mesh::BuildV_V()
{

}
void Mesh::BuildV_E()
{

}
void Mesh::BuildV_F()
{

}
void Mesh::BuildV_C() {
//    if (m_cellType == HEXAHEDRA || m_cellType == TETRAHEDRA || m_cellType == WEDGE || m_cellType == PENTAHEDRON || m_cellType == POLYHEDRA)
        for (auto& c : C)
            for (auto vid : c.Vids)
                V[vid].N_Cids.push_back(c.id);
    /*else */if (m_cellType == TRIANGLE || m_cellType == QUAD || m_cellType == POLYGON)
        for (auto& c : C)
			for (auto vid : c.Vids)
                V[vid].N_Fids.push_back(c.id);
}
void Mesh::BuildE_V() {
    for (auto& e: E) {
        for (auto vid : e.Vids) {
            for (auto nvid : V[vid].N_Vids) {
                e.N_Vids.push_back(nvid);
            }
        }
		for (auto vid : e.Vids) {
            for (auto iter = e.N_Vids.begin(); iter != e.N_Vids.end();) {
                if (*iter == vid) iter = e.N_Vids.erase(iter);
                else ++iter;
            }
        }
    }
}
void Mesh::BuildE_E() {
    for (size_t i = 0; i < E.size(); i++) {
        for (size_t j = 0; j < E[i].Vids.size(); j++) {
            for (size_t k = 0; k < V[E[i].Vids[j]].N_Eids.size(); k++) {
                E[i].N_Eids.push_back(V[E[i].Vids[j]].N_Eids[k]);
            }
        }
        for (size_t j = 0; j < E[i].Vids.size(); j++) {
            for (std::vector<size_t>::iterator iter = E[i].N_Eids.begin(); iter != E[i].N_Eids.end();) {
                if (*iter == E[i].id) iter = E[i].N_Eids.erase(iter);
                else ++iter;
            }
        }
    }
}
void Mesh::BuildE_F() {
    for (size_t i = 0; i < E.size(); i++) {
        for (size_t j = 0; j < E[i].Vids.size(); j++) {
            for (size_t k = 0; k < V[E[i].Vids[j]].N_Fids.size(); k++) {
                const size_t N_Fid = V[E[i].Vids[j]].N_Fids[k];
                if (std::find(E[i].N_Fids.begin(), E[i].N_Fids.end(), N_Fid) == E[i].N_Fids.end()) E[i].N_Fids.push_back(N_Fid);
            }
        }
    }
}
void Mesh::BuildE_C() {
	for (auto& e : E)
		for (auto cid0 : V[e.Vids[0]].N_Cids)
			for (auto cid1 : V[e.Vids[1]].N_Cids)
				if (cid0 == cid1) e.N_Cids.push_back(cid0);

	for (auto& e : E) {
		std::set<size_t> N_Cids(e.N_Cids.begin(), e.N_Cids.end());
		e.N_Cids.clear();
		for (auto ncid : N_Cids)
			e.N_Cids.push_back(ncid);
	}
}
void Mesh::BuildF_V()
{

}
void Mesh::BuildF_E()
{

}
void Mesh::BuildF_F()
{
    if (!F[0].N_Fids.empty())
        return;

    for (size_t i = 0; i < F.size(); i++) {
        Face& face = F.at(i);
        for (size_t j = 0; j < face.Vids.size(); j++) {
            const Vertex& vertex = V.at(face.Vids.at(j));
            for (size_t k = 0; k < vertex.N_Fids.size(); k++) {
                Face& neighborFace = F.at(vertex.N_Fids.at(k));
                if (neighborFace.id != face.id)
                    face.N_Fids.push_back(neighborFace.id);
            }
        }

        std::sort(face.N_Fids.begin(), face.N_Fids.end());
        std::vector<size_t>::iterator iter = std::unique(face.N_Fids.begin(), face.N_Fids.end());
        face.N_Fids.resize(std::distance(face.N_Fids.begin(), iter));
    }
}
void Mesh::BuildF_C()
{

}
void Mesh::BuildC_V()
{

}
void Mesh::BuildC_E()
{
    if (m_cellType == HEXAHEDRA)
        for (size_t i = 0; i < C.size(); i++) {
            Cell& cell = C.at(i);
            cell.Eids.resize(12);
            for (size_t k = 0; k < 12; k++) {
                Edge e1(2);
                e1.Vids[0] = cell.Vids[HexEdge[k][0]];
                e1.Vids[1] = cell.Vids[HexEdge[k][1]];
                for (size_t j = 0; j < cell.Fids.size(); j++) {
                    bool found = false;
                    const Face& face = F.at(cell.Fids.at(j));
                    for (size_t n = 0; n < face.Eids.size(); n++) {
                        const Edge& e2 = E[face.Eids[n]];
                        if (e1 == e2) {
                            cell.Eids[k] = e2.id;
                            found = true;
                            break;
                        }
                    }
                    if (found)
                        break;
                }
            }
        }
    else if (m_cellType == TETRAHEDRA)
        for (size_t i = 0; i < C.size(); i++) {
            Cell& cell = C.at(i);
            cell.Eids.resize(6);
            for (size_t k = 0; k < 6; k++) {
                Edge e1(2);
                e1.Vids[0] = cell.Vids[TetEdge[k][0]];
                e1.Vids[1] = cell.Vids[TetEdge[k][1]];
                for (size_t j = 0; j < cell.Fids.size(); j++) {
                    bool found = false;
                    const Face& face = F.at(cell.Fids.at(j));
                    for (size_t n = 0; n < face.Eids.size(); n++) {
                        const Edge& e2 = E[face.Eids[n]];
                        if (e1 == e2) {
                            cell.Eids[k] = e2.id;
                            found = true;
                            break;
                        }
                    }
                    if (found)
                        break;
                }
            }
        }
	else if (m_cellType == POLYHEDRA)
		for (auto & c : C) {
			c.Eids.resize(6);
			auto iter = cell_edges.find(c.cellType);
			c.Eids.resize(iter->second.size());
			for (size_t k = 0; k < iter->second.size(); k++) {
				auto& vids = iter->second[k];
				Edge e1(vids);
				for (size_t j = 0; j < c.Fids.size(); j++) {
					bool found = false;
					const Face& face = F.at(c.Fids.at(j));
					for (auto eid : face.Eids) {
						const Edge& e2 = E[eid];
						if (e1 == e2) {
							c.Eids[k] = e2.id;
							found = true;
							break;
						}
					}
					if (found) break;
				}
			}
		}
}
void Mesh::BuildC_F()
{

}
void Mesh::BuildC_C()
{

}

void Mesh::BuildAllConnectivities()
{
    BuildV_C();
    BuildE(); // BuildV_V(); BuildV_E(); BuildV_F();
    BuildF(); // BuildF_C();
    BuildC_E();
    BuildF_F();
//    if (m_cellType != HEXAHEDRA){
//        BuildE_V();
//        BuildE_E();
//        BuildE_F();
//    }
}

//static void set_cross(const std::vector<size_t> set1, const std::vector<size_t> set2, std::vector<size_t> &result_set)
//{
//    result_set.clear();
//    for (int i = 0; i < set1.size(); i++)
//    {
//        bool inside = false;
//        for (int j = 0; j < set2.size(); j++)
//        {
//            if (set2[j] == set1[i])
//            {
//                inside = true;
//                break;
//            }
//        }
//        if (inside)
//            result_set.push_back(set1[i]);
//    }
//}

void Mesh::FixOrientation() {
	if (m_cellType == POLYGON || m_cellType == TRIANGLE || m_cellType == QUAD)
	for (auto& e : E) {
		auto& f0 = F[e.N_Fids[0]];
		auto& f1 = F[e.N_Fids[1]];
		auto feid = MAXID;
		for (auto eid : f1.Eids)
			if (eid == e.id) {
				feid = eid;
				break;
			}
		auto fvid0 = MAXID;
		auto fvid1 = MAXID;
		for (auto i = 0; i < f1.Vids.size(); ++i) {
			auto fvid0_ = f1.Vids.at(i);
			auto fvid1_ = f1.Vids.at((i + 1) % f1.Vids.size());
			if ((fvid0_ == e.Vids[0] && fvid1_ == e.Vids[1]) || (fvid0_ == e.Vids[1] && fvid1_ == e.Vids[0])) {
				fvid0 = fvid0_;
				fvid1 = fvid1_;
				break;
			}
		}
		auto& v0 = V[fvid0];
		auto& v1 = V[fvid1];
		auto f0center = 0.25 * (V[f0.Vids[0]].xyz() + V[f0.Vids[1]].xyz() + V[f0.Vids[2]].xyz() + V[f0.Vids[3]].xyz());
		auto f1center = 0.25 * (V[f1.Vids[0]].xyz() + V[f1.Vids[1]].xyz() + V[f1.Vids[2]].xyz() + V[f1.Vids[3]].xyz());
		auto vec_f = glm::cross(f1center - f0center, v1.xyz() - v0.xyz());
		auto vec_o = glm::cross(V[f0.Vids[2]].xyz() - V[f0.Vids[1]].xyz(), V[f0.Vids[0]].xyz() - V[f0.Vids[1]].xyz());
		auto sign = glm::dot(vec_f, vec_o);
		if (sign < 0) std::swap(e.N_Fids[0], e.N_Fids[1]);
	}
}

void Mesh::ExtractBoundary() {
    if (m_cellType == POLYGON || m_cellType == TRIANGLE || m_cellType == QUAD) {
        for (auto& e : E)
            if (e.N_Fids.size() == 1) {
                hasBoundary = true;
                break;
            }
        for (auto& e : E) {
            if (e.N_Fids.size() == 1) {
                V.at(e.Vids[0]).isBoundary = true;
                V.at(e.Vids[1]).isBoundary = true;
                e.isBoundary = true;
                V.at(e.Vids[0]).type = FEATURE;
                V.at(e.Vids[1]).type = FEATURE;
            }
        }
//        for (auto& v : V) {
//            if (v.N_Fids.size() == 1) {
//                v.type = CORNER;
//                v.isCorner = true;
//            }
//        }
        for (size_t i = 0; i < F.size(); i++)
            F[i].isBoundary = true;
        return;
    }

    for (size_t i = 0; i < E.size(); i++) {
        std::vector<size_t> neighborhs = E[i].N_Cids;
        E[i].N_Cids.clear();
        size_t pre_ = neighborhs[0];
        size_t aft_;
        E[i].N_Cids.push_back(pre_);
        std::vector<bool> tags(neighborhs.size(), false);
        tags[0] = true;
        for (size_t m = 1; m < neighborhs.size(); m++) {
            for (size_t j = 0; j < neighborhs.size(); j++) {
                if (tags[j]) continue;

                std::vector<size_t> fs_pre, fs_cur;
                if (m_cellType == HEXAHEDRA) for (int k = 0; k < 6; k++) {
                    fs_pre.push_back(C[pre_].Fids[k]);
                    fs_cur.push_back(C[neighborhs[j]].Fids[k]);
                }
                else if (m_cellType == TETRAHEDRA) for (int k = 0; k < 4; k++) {
                    fs_pre.push_back(C[pre_].Fids[k]);
                    fs_cur.push_back(C[neighborhs[j]].Fids[k]);
                }
                std::vector<size_t> sharedf;
                set_cross(fs_pre, fs_cur, sharedf);
                if (sharedf.size()) {
                    E[i].N_Cids.push_back(neighborhs[j]);
                    tags[j] = true;
                    pre_ = neighborhs[j];
                }
            }
        }
        if (E[i].N_Cids.size() != neighborhs.size()) {
            //cout<<"ERROR here"<<endl;//boundary edge have not been handled yet
            E[i].N_Cids = neighborhs;
        }
    }
    for (auto& f : F)
        if (f.N_Cids.size() == 1) {
            hasBoundary = true;
            break;
        }
    for (size_t i = 0; i < F.size(); i++) {
        F[i].isBoundary = false;
        if (F[i].N_Cids.size() == 1) {
            F[i].isBoundary = true;
            if (m_cellType == HEXAHEDRA)
                for (size_t j = 0; j < 4; j++) {
                    E[F[i].Eids[j]].isBoundary = true;
                    V[F[i].Vids[j]].isBoundary = true;
                }
            else if (m_cellType == TETRAHEDRA)
                for (size_t j = 0; j < 3; j++) {
                    E[F[i].Eids[j]].isBoundary = true;
                    V[F[i].Vids[j]].isBoundary = true;
                }
        }
    }
    for (size_t i = 0; i < C.size(); i++) {
        C[i].isBoundary = false;
        for (size_t j = 0; j < C[i].N_Fids.size(); j++) {
            size_t fid = C[i].N_Fids[j];
            if (F[fid].isBoundary)
                C[i].isBoundary = true;
        }
    }
}

size_t Mesh::ExtractLayers() {
    std::vector<bool> savedCellBoundary(C.size(), false);
    for (size_t i = 0; i < C.size(); i++)
        savedCellBoundary.at(i) = C[i].isBoundary;

    std::vector<Cell> surfaceCells;
    std::vector<Cell> innerCells;

    int layerCount = 1;
    for (size_t i = 0; i < C.size(); i++) {
        if (C[i].isBoundary) surfaceCells.push_back(C[i]);
        else innerCells.push_back(C[i]);
    }
    while (!innerCells.empty()) {
        // -------------------------------------------------------
        // build outLayer
        {
            Layer layer;
            for (size_t i = 0; i < surfaceCells.size(); i++)
                layer.Cids.push_back(surfaceCells[i].id);

            for (size_t i = 0; i < layer.Cids.size(); i++) {
                const Cell& cell = C[layer.Cids[i]];
                for (size_t j = 0; j < cell.Fids.size(); j++)
                    layer.Fids.push_back(cell.Fids[j]);
                for (size_t j = 0; j < cell.Eids.size(); j++)
                    layer.Eids.push_back(cell.Eids[j]);
                for (size_t j = 0; j < cell.Vids.size(); j++)
                    layer.Vids.push_back(cell.Vids[j]);
            }
            std::sort(layer.Fids.begin(), layer.Fids.end());
            std::vector<size_t>::iterator iterF = std::unique(layer.Fids.begin(), layer.Fids.end());
            layer.Fids.resize(std::distance(layer.Fids.begin(), iterF));

            std::sort(layer.Eids.begin(), layer.Eids.end());
            std::vector<size_t>::iterator iterE = std::unique(layer.Eids.begin(), layer.Eids.end());
            layer.Eids.resize(std::distance(layer.Eids.begin(), iterE));

            std::sort(layer.Vids.begin(), layer.Vids.end());
            std::vector<size_t>::iterator iterV = std::unique(layer.Vids.begin(), layer.Vids.end());
            layer.Vids.resize(std::distance(layer.Vids.begin(), iterV));

            layer.fixed.resize(layer.Vids.size(), false);
            for (size_t i = 0; i < layer.Vids.size(); i++)
                if (V[layer.Vids[i]].isBoundary) layer.fixed.at(i) = true;
            layers.push_back(layer);
        }
        // build innerLayer
        {
            Layer layer;
            for (size_t i = 0; i < innerCells.size(); i++)
                layer.Cids.push_back(innerCells[i].id);

            for (size_t i = 0; i < layer.Cids.size(); i++){
                const Cell& cell = C[layer.Cids[i]];
                for (size_t j = 0; j < cell.Fids.size(); j++)
                    layer.Fids.push_back(cell.Fids[j]);
                for (size_t j = 0; j < cell.Eids.size(); j++)
                    layer.Eids.push_back(cell.Eids[j]);
                for (size_t j = 0; j < cell.Vids.size(); j++)
                    layer.Vids.push_back(cell.Vids[j]);
            }
            std::sort(layer.Fids.begin(), layer.Fids.end());
            std::vector<size_t>::iterator iterF = std::unique(layer.Fids.begin(), layer.Fids.end());
            layer.Fids.resize(std::distance(layer.Fids.begin(), iterF));

            std::sort(layer.Eids.begin(), layer.Eids.end());
            std::vector<size_t>::iterator iterE = std::unique(layer.Eids.begin(), layer.Eids.end());
            layer.Eids.resize(std::distance(layer.Eids.begin(), iterE));

            std::sort(layer.Vids.begin(), layer.Vids.end());
            std::vector<size_t>::iterator iterV = std::unique(layer.Vids.begin(), layer.Vids.end());
            layer.Vids.resize(std::distance(layer.Vids.begin(), iterV));

            layer.fixed.resize(layer.Vids.size(), false);
            for (size_t i = 0; i < layer.Vids.size(); i++)
                if (V[layer.Vids[i]].isBoundary) layer.fixed.at(i) = true;
            innerLayers.push_back(layer);
        }
        // -------------------------------------------------------
        std::string surfaceCellsFileName = std::string("OutLayer") + std::to_string(layerCount) + ".vtk";
        MeshFileWriter surfaceCellsWriter(V, surfaceCells, surfaceCellsFileName.c_str(), HEXAHEDRA);
        surfaceCellsWriter.SetFixFlag(true);
        surfaceCellsWriter.WriteFile();

        std::string innerCellsFileName = std::string("InnerLayer") + std::to_string(layerCount) + ".vtk";
        MeshFileWriter innerCellsWriter(V, innerCells, innerCellsFileName.c_str(), HEXAHEDRA);
        innerCellsWriter.SetFixFlag(true);
        innerCellsWriter.WriteFile();

        MeshFileReader innerCellsReader(innerCellsFileName.c_str());
        Mesh& mesh = (Mesh&)innerCellsReader.GetMesh();
        mesh.BuildAllConnectivities();
        mesh.ExtractBoundary();

        std::vector<Cell> newC = innerCells;
        surfaceCells.clear();
        innerCells.clear();
        for (size_t i = 0; i < newC.size(); i++) {
            if (mesh.C[i].isBoundary) surfaceCells.push_back(newC[i]);
            else innerCells.push_back(newC[i]);
        }

        layerCount++;
    }

    for (size_t i = 0; i < C.size(); i++)
        C[i].isBoundary = savedCellBoundary.at(i);

    return layerCount;
}

//size_t Mesh::ExtractLayers()
//{
//    std::vector<Cell> surfaceCells;
//    std::vector<Cell> innerCells;
//
//    int layerCount = 1;
//    for (size_t i = 0; i < C.size(); i++)
//    {
//        if (C[i].isBoundary)
//            surfaceCells.push_back(Cell(C[i].Vids));
//        else
//            innerCells.push_back(Cell(C[i].Vids));
//    }
//    while (!innerCells.empty()) {
//        std::string surfaceCellsFileName = std::string("OutLayer") + std::to_string(layerCount) + ".vtk";
//        MeshFileWriter surfaceCellsWriter(V, surfaceCells, surfaceCellsFileName.c_str(), HEXAHEDRA);
//        surfaceCellsWriter.SetFixFlag(true);
//        surfaceCellsWriter.WriteFile();
//
//        std::string innerCellsFileName = std::string("InnerLayer") + std::to_string(layerCount) + ".vtk";
//        MeshFileWriter innerCellsWriter(V, innerCells, innerCellsFileName.c_str(), HEXAHEDRA);
//        innerCellsWriter.WriteFile();
//
//        MeshFileReader innerCellsReader(innerCellsFileName.c_str());
//        Mesh& mesh = (Mesh&)innerCellsReader.GetMesh();
//        mesh.BuildAllConnectivities();
//        mesh.ExtractBoundary();
//
//        surfaceCells.clear();
//        innerCells.clear();
//        for (size_t i = 0; i < mesh.C.size(); i++)
//        {
//            if (mesh.C[i].isBoundary)
//                surfaceCells.push_back(Cell(mesh.C[i].Vids));
//            else
//                innerCells.push_back(Cell(mesh.C[i].Vids));
//        }
//
//        layerCount++;
//    }
//
//    {
//        Layer layer;
//        for (size_t i = 0; i < C.size(); i++)
//            if (C[i].isBoundary)
//                layer.Cids.push_back(C[i].id);
//        for (size_t i = 0; i < layer.Cids.size(); i++){
//            const Cell& cell = C[layer.Cids[i]];
//            for (size_t j = 0; j < cell.Fids.size(); j++)
//                layer.Fids.push_back(cell.Fids[j]);
//            for (size_t j = 0; j < cell.Eids.size(); j++)
//                layer.Eids.push_back(cell.Eids[j]);
//            for (size_t j = 0; j < cell.Vids.size(); j++)
//                layer.Vids.push_back(cell.Vids[j]);
//        }
//        std::sort(layer.Fids.begin(), layer.Fids.end());
//        std::vector<size_t>::iterator iterF = std::unique(layer.Fids.begin(), layer.Fids.end());
//        layer.Fids.resize(std::distance(layer.Fids.begin(), iterF));
//
//        std::sort(layer.Eids.begin(), layer.Eids.end());
//        std::vector<size_t>::iterator iterE = std::unique(layer.Eids.begin(), layer.Eids.end());
//        layer.Eids.resize(std::distance(layer.Eids.begin(), iterE));
//
//        std::sort(layer.Vids.begin(), layer.Vids.end());
//        std::vector<size_t>::iterator iterV = std::unique(layer.Vids.begin(), layer.Vids.end());
//        layer.Vids.resize(std::distance(layer.Vids.begin(), iterV));
//
//        layer.fixed.resize(layer.Vids.size(), false);
//        for (size_t i = 0; i < layer.Vids.size(); i++)
//            if (V[layer.Vids[i]].isBoundary)
//                layer.fixed.at(i) = true;
//        layers.push_back(layer);
//    }
//
//    return layerCount;
//}

void Mesh::ExtractSingularities() {
    if (m_cellType == HEXAHEDRA) {
        for (size_t i = 0; i < E.size(); i++) {
            Edge& edge = E.at(i);
            if (edge.isBoundary) {
                if (edge.N_Cids.size() != 2) edge.isSingularity = true;
            }
            else {
                if (edge.N_Cids.size() != 4) edge.isSingularity = true;
            }

            if (edge.isSingularity) {
                Vertex& v1 = V.at(edge.Vids.at(0));
                Vertex& v2 = V.at(edge.Vids.at(1));
                v1.isSingularity = true;
                v2.isSingularity = true;
            }
        }
    }
    else if (m_cellType == QUAD) {
        for (size_t i = 0; i < V.size(); i++) {
            Vertex& v = V.at(i);
            if (v.isBoundary) {
                if (v.N_Fids.size() != 2) v.isSingularity = true;
            } else {
                if (v.N_Fids.size() != 4) v.isSingularity = true;
            }
        }
    }
}

void Mesh::ExtractTwoRingNeighborSurfaceFaceIdsForEachVertex(int N/* = 2*/)
{
    for (size_t i = 0; i < V.size(); i++) {
        Vertex& v = V[i];
        if (!v.isBoundary) continue;
        v.twoRingNeighborSurfaceFaceIds.clear();
        int n = N;
        std::vector<size_t> surfaceFaceVids(1, v.id);
        while (n-- != 0) {
            for (size_t j = 0; j < surfaceFaceVids.size(); j++) {
                const Vertex& vv = V[surfaceFaceVids[j]];
                for (size_t k = 0; k < vv.N_Fids.size(); k++) {
                    const Face& neighboringF = F[vv.N_Fids[k]];
                    if (neighboringF.isBoundary) v.twoRingNeighborSurfaceFaceIds.push_back(neighboringF.id);
                }
            }
            for (size_t j = 0; j < v.twoRingNeighborSurfaceFaceIds.size(); j++) {
                const Face& f = F[v.twoRingNeighborSurfaceFaceIds[j]];
                for (size_t k = 0; k < f.Vids.size(); k++)
                    surfaceFaceVids.push_back(f.Vids[k]);
            }
            std::sort(surfaceFaceVids.begin(), surfaceFaceVids.end());
            std::vector<size_t>::iterator iter = std::unique(surfaceFaceVids.begin(), surfaceFaceVids.end());
            surfaceFaceVids.resize(std::distance(surfaceFaceVids.begin(), iter));
        }
        std::sort(v.twoRingNeighborSurfaceFaceIds.begin(), v.twoRingNeighborSurfaceFaceIds.end());
        std::vector<size_t>::iterator iter = std::unique(v.twoRingNeighborSurfaceFaceIds.begin(), v.twoRingNeighborSurfaceFaceIds.end());
        v.twoRingNeighborSurfaceFaceIds.resize(std::distance(v.twoRingNeighborSurfaceFaceIds.begin(), iter));
    }
//    for (size_t i = 0; i < V.size(); i++) {
//        Vertex& v = V[i];
//        if (!v.isBoundary) continue;
//        v.twoRingNeighborSurfaceFaceIds.clear();
//        for (size_t j = 0; j < v.N_Vids.size(); j++) {
//            const Vertex& neighboringV = V[v.N_Vids[j]];
//            for (size_t k = 0; k < neighboringV.N_Fids.size(); k++) {
//                const Face& neighboringF = F[neighboringV.N_Fids[k]];
//                if (neighboringF.isBoundary) v.twoRingNeighborSurfaceFaceIds.push_back(neighboringF.id);
//            }
//        }
//        std::sort(v.twoRingNeighborSurfaceFaceIds.begin(), v.twoRingNeighborSurfaceFaceIds.end());
//        std::vector<size_t>::iterator iter = std::unique(v.twoRingNeighborSurfaceFaceIds.begin(), v.twoRingNeighborSurfaceFaceIds.end());
//        v.twoRingNeighborSurfaceFaceIds.resize(std::distance(v.twoRingNeighborSurfaceFaceIds.begin(), iter));
//    }
}

void Mesh::Zoom(const glm::dvec3& ref, const double scale/* = 1*/) {
    if (scale == 1.0)  return;
//    if (glm::length(m_center) == 0)
//        m_center = GetCenter(V);
    for (size_t i = 0; i < V.size(); i++) {
        Vertex& v = V.at(i);
        const glm::dvec3 dir = v - ref;
        const glm::dvec3 d(dir.x, dir.y, dir.z);
        v.x = scale * d.x + ref.x;
        v.y = scale * d.y + ref.y;
        v.z = scale * d.z + ref.z;
    }
}

bool Mesh::HasBoundary() const {
    return hasBoundary;
}

const unsigned int HexEdges[12][2] =
{
    { 0, 1 },
    { 1, 2 },
    { 2, 3 },
    { 3, 0 },
    { 4, 5 },
    { 5, 6 },
    { 6, 7 },
    { 7, 4 },
    { 0, 4 },
    { 1, 5 },
    { 2, 6 },
    { 3, 7 },
};

const std::vector<std::vector<size_t>> HexEdge = {
	{ 0, 1 },
	{ 1, 2 },
	{ 2, 3 },
	{ 3, 0 },
	{ 4, 5 },
	{ 5, 6 },
	{ 6, 7 },
	{ 7, 4 },
	{ 0, 4 },
	{ 1, 5 },
	{ 2, 6 },
	{ 3, 7 },
};


const unsigned int TetEdges[6][2] =
{
    { 0, 1 },
    { 0, 2 },
    { 0, 3 },
    { 1, 2 },
    { 1, 3 },
    { 2, 3 }
};

const std::vector<std::vector<size_t>> TetEdge = {
	{ 0, 1 },
	{ 0, 2 },
	{ 0, 3 },
	{ 1, 2 },
	{ 1, 3 },
	{ 2, 3 }
};

//        0
//        /\
//       /| \
//   	/ |  \
//    1/__|___\2
//     | 3|   |
//     |  /\  |
//     | /  \ |
//     |/    \|
//    4/______\5
const unsigned int WedgeEdges[9][2] =
{
	{ 0, 1 },
	{ 1, 2 },
	{ 2, 1 },
	{ 3, 4 },
	{ 4, 5 },
	{ 5, 3 },
	{ 0, 3 },
	{ 1, 4 },
	{ 2, 5 }
};

const std::vector<std::vector<size_t>> WedgeEdge = {
	{ 0, 1 },
	{ 1, 2 },
	{ 2, 1 },
	{ 3, 4 },
	{ 4, 5 },
	{ 5, 3 },
	{ 0, 3 },
	{ 1, 4 },
	{ 2, 5 }
};

const std::vector<std::vector<size_t>> WedgeFace =
{
	{ 0, 1, 2 },
	{ 3, 4, 5 },
	{ 0, 3, 5, 2 },
	{ 2, 5, 4, 1 },
	{ 1, 4, 3, 0 },
};

const unsigned int PentaEdges[15][2] =
{
	{ 0, 1 },
	{ 1, 2 },
	{ 2, 3 },
	{ 3, 4 },
    { 4, 0 },
	{ 5, 6 },
	{ 6, 7 },
	{ 7, 8 },
    { 8, 9 },
    { 9, 5 },
	{ 0, 5 },
	{ 1, 6 },
	{ 2, 7 },
	{ 3, 8 },
    { 4, 9 }
};

const std::vector<std::vector<size_t>> PentaEdge = {
	{ 0, 1 },
	{ 1, 2 },
	{ 2, 3 },
	{ 3, 4 },
	{ 4, 0 },
	{ 5, 6 },
	{ 6, 7 },
	{ 7, 8 },
	{ 8, 9 },
	{ 9, 5 },
	{ 0, 5 },
	{ 1, 6 },
	{ 2, 7 },
	{ 3, 8 },
	{ 4, 9 }
};

const std::vector<std::vector<size_t>> PentaFace =
{
	{ 0, 1, 2, 3, 4 },
	{ 5, 6, 7, 8, 9 },
	{ 0, 5, 6, 1 },
	{ 1, 6, 7, 2 },
	{ 2, 7, 8, 3 },
	{ 3, 8, 9, 4 },
	{ 4, 9, 5, 0 },
};

const unsigned int QuadEdge[4][2] =
{
    { 0, 1 },
    { 1, 2 },
    { 2, 3 },
    { 3, 0 }
};

const unsigned int TriEdge[3][2] =
{
    { 0, 1 },
    { 1, 2 },
    { 2, 0 },
};

const unsigned int HexPoint_Points[8][3] =
{
    { 1, 3, 4 },
    { 2, 0, 5 },
    { 3, 1, 6 },
    { 0, 2, 7 },
    { 7, 5, 0 },
    { 4, 6, 1 },                                  //const unsigned int HexEdge[12][2] =
    { 5, 7, 2 },                                  //{
    { 6, 4, 3 }                                   //    { 0, 1 }, // 0
};                                                //    { 1, 2 }, // 1
                                                  //    { 2, 3 }, // 2
const unsigned int HexPoint_Edges[8][3] =         //    { 3, 0 }, // 3
{                                                 //    { 4, 5 }, // 4
    { 0, 3, 8 },  //{ (0,1), (0,3), (0,4) },      //    { 5, 6 }, // 5
    { 1, 0, 9 },  //{ (1,2), (1,0), (1,5) },      //    { 6, 7 }, // 6
    { 2, 1, 10},  //{ (2,3), (2,1), (2,6) },      //    { 7, 4 }, // 7
    { 3, 2, 11},  //{ (3,0), (3,2), (3,7) },      //    { 0, 4 }, // 8
    { 7, 4, 8 },  //{ (4,7), (4,5), (4,0) },      //    { 1, 5 }, // 9
    { 4, 5, 9 },  //{ (5,4), (5,6), (5,1) },      //    { 2, 6 }, // 10
    { 5, 6, 10},  //{ (6,5), (6,7), (6,2) },      //    { 3, 7 }, // 11
    { 6, 7, 11}   //{ (7,6), (7,4), (7,3) },      //};
};

const unsigned int TetPoint_Points[4][3] =
{
    { 1, 2, 3 },
    { 0, 2, 3 },
    { 0, 1, 3 },
    { 0, 1, 2 }
};

const unsigned int QuadPoint_Points[4][2] =
{
    { 1, 3 },
    { 2, 0 },
    { 3, 1 },
    { 0, 2 }
};

const unsigned int TriPoint_Points[3][2] =
{
    { 1, 2 },
    { 2, 0 },
    { 0, 1 }
};

const unsigned int HexFaces1[6][4] =
{
    {0, 3, 2, 1},
    {0, 1, 5, 4},
    {0, 4, 7, 3},
    {1, 2, 6, 5},
    {3, 7, 6, 2},
    {4, 5, 6, 7}
};

const unsigned int HexFaces[6][4] =
{
    {0, 3, 2, 1},
    {4, 5, 6, 7},
    {0, 4, 7, 3},
    {1, 2, 6, 5},
    {0, 1, 5, 4},
    {3, 7, 6, 2}

//    {0, 3, 2, 1},
//    {0, 1, 5, 4},
//    {0, 4, 7, 3},
//    {1, 2, 6, 5},
//    {2, 3, 7, 6},
//    {4, 5, 6, 7}
};

const std::vector<std::vector<size_t>> HexFace = {
	{0, 3, 2, 1},
	{4, 5, 6, 7},
	{0, 4, 7, 3},
	{1, 2, 6, 5},
	{0, 1, 5, 4},
	{3, 7, 6, 2}
};

const unsigned int TetFaces[4][3] =
{
    {0, 2, 1},
    {0, 3, 2},
    {0, 1, 3},
    {1, 2, 3}
};

const std::vector<std::vector<size_t>> TetFace = {
    {0, 2, 1},
    {0, 3, 2},
    {0, 1, 3},
    {1, 2, 3}
};

const std::map<size_t, std::vector<std::vector<size_t>>> cell_faces = {
	{ VTK_TETRA, TetFace },
    { VTK_HEXAHEDRON, HexFace },
    { VTK_WEDGE, WedgeFace },
    { VTK_PENTAGONAL_PRISM, PentaFace }
};

const std::map<size_t, std::vector<std::vector<size_t>>> cell_edges = {
	{ VTK_TETRA, TetEdge },
	{ VTK_HEXAHEDRON, HexEdge },
	{ VTK_WEDGE, WedgeEdge },
	{ VTK_PENTAGONAL_PRISM, PentaEdge }
};

const unsigned int HexPoint_Faces[8][3][4] =
{
    {
    {0, 3, 2, 1},
    {0, 4, 7, 3},
    {0, 1, 5, 4}},

    {
    {0, 3, 2, 1},
    {0, 1, 5, 4},
    {1, 2, 6, 5}},

    {
    {0, 3, 2, 1},
    {1, 2, 6, 5},
    {3, 7, 6, 2}},

    {
    {0, 3, 2, 1},
    {3, 7, 6, 2},
    {0, 4, 7, 3}},
    //
    {
    {4, 5, 6, 7},
    {0, 4, 7, 3},
    {0, 1, 5, 4}},

    {
    {4, 5, 6, 7},
    {0, 1, 5, 4},
    {1, 2, 6, 5}},

    {
    {4, 5, 6, 7},
    {1, 2, 6, 5},
    {3, 7, 6, 2}},

    {
    {4, 5, 6, 7},
    {3, 7, 6, 2},
    {0, 4, 7, 3}}
};

const unsigned int TetPoint_Faces[4][3][3] =
{
    {{0, 1, 2},
    {0, 2, 3},
    {0, 3, 1}},

    {{0, 1, 2},
    {1, 3, 2},
    {0, 3, 1}},

    {{0, 1, 2},
    {0, 2, 3},
    {1, 3, 2}},

    {{1, 3, 2},
    {0, 2, 3},
    {0, 3, 1}}
};

const unsigned long hexTet[5][4] =
{
    {0, 4, 5, 7},
    {2, 5, 6, 7},
    {0, 2, 3, 7},
    {0, 1, 2, 5},
    {0, 2, 5, 7}
};

///////////////////////////////////////////////
const unsigned int HexPoint_Points_CW[8][3] =
{
    { 1, 3, 4 },
    { 0, 2, 5 },
    { 1, 3, 6 },
    { 0, 2, 7 },
    { 0, 5, 7 },
    { 1, 4, 6 },
    { 2, 5, 7 },
    { 3, 4, 6 }
};

const unsigned int TetPoint_Points_CW[4][3] =
{
    { 1, 2, 3 },
    { 0, 2, 3 },
    { 0, 1, 3 },
    { 0, 1, 2 }
};

const unsigned int HexFaces_CW[6][4] =
{
    {0, 1, 2, 3},
    {0, 4, 5, 1},
    {0, 3, 7, 4},
    {1, 5, 6, 2},
    {2, 6, 7, 3},
    {4, 5, 6, 7}
};

const unsigned int TetFaces_CW[4][3] =
{
    {0, 1, 2},
    {0, 2, 3},
    {0, 3, 1},
    {1, 3, 2}
};

const unsigned int HexPoint_Faces_CW[8][3][4] =
{
    {{0, 1, 2, 3},
    {0, 3, 7, 4},
    {0, 4, 5, 1}},

    {{0, 1, 2, 3},
    {0, 4, 5, 1},
    {1, 5, 6, 2}},

    {{0, 1, 2, 3},
    {1, 5, 6, 2},
    {2, 6, 7, 3}},

    {{0, 1, 2, 3},
    {2, 6, 7, 3},
    {0, 3, 7, 4}},
    //
    {{4, 5, 6, 7},
    {0, 3, 7, 4},
    {0, 4, 5, 1}},

    {{4, 5, 6, 7},
    {0, 4, 5, 1},
    {1, 5, 6, 2}},

    {{4, 5, 6, 7},
    {1, 5, 6, 2},
    {2, 6, 7, 3}},

    {{4, 5, 6, 7},
    {2, 6, 7, 3},
    {0, 3, 7, 4}}
};

const unsigned int TetPoint_Faces_CW[4][3][3] =
{
    {{0, 1, 2},
    {0, 2, 3},
    {0, 3, 1}},

    {{0, 1, 2},
    {1, 3, 2},
    {0, 3, 1}},

    {{0, 1, 2},
    {0, 2, 3},
    {1, 3, 2}},

    {{1, 3, 2},
    {0, 2, 3},
    {0, 3, 1}}
};

const unsigned long INVALID_NUM = 3277;

const unsigned int HexTrippleEdge[12][6] =
{
    { 0, 1, 4, 3, 5, 2},
    { 1, 2, 5, 0, 6, 3},
    { 2, 3, 6, 1, 7, 0},
    { 3, 0, 7, 2, 4, 1},
    { 0, 4, 3, 1, 7, 5},
    { 1, 5, 0, 2, 4, 6},
    { 2, 6, 1, 3, 5, 7},
    { 3, 7, 2, 0, 6, 4},
    { 4, 5, 7, 0, 6, 1},
    { 5, 6, 4, 1, 7, 2},
    { 6, 7, 5, 2, 4, 3},
    { 7, 4, 6, 3, 5, 0},
};

bool Mesh::IsPointInside(const glm::dvec3& orig, const glm::dvec3 dir) const {
    bool bInside = false;
    for (size_t i = 0; i < F.size(); i++) {
        const Face& tri = F.at(i);
        if (!tri.isBoundary) continue;
        const glm::dvec3& v0 = V.at(tri.Vids.at(0));
        const glm::dvec3& v1 = V.at(tri.Vids.at(1));
        const glm::dvec3& v2 = V.at(tri.Vids.at(2));
        glm::dvec3 position;
        if (glm::intersectRayTriangle(orig, dir, v0, v1, v2, position) || glm::intersectRayTriangle(orig, dir, v0, v2, v1, position)) {
            bInside = !bInside;
        }
    }
    return bInside;
}
