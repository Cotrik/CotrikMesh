/*
 * Dual.h
 *
 *  Created on: Oct 30, 2017
 *      Author: cotrik
 */

#ifndef LIBCOTRIK_SRC_DUAL_H_
#define LIBCOTRIK_SRC_DUAL_H_

#include "Mesh.h"
#include "MeshFileReader.h"
#include "MeshFileWriter.h"

class DualVertex: public Vertex {
public:
    DualVertex()
    : Vertex()
    {}
    DualVertex(const DualVertex& r)
    : Vertex(r)
    {}
    DualVertex(const glm::vec3& v)
    : Vertex(v)
    {}
    virtual ~DualVertex()
    {}
    DualVertex& operator =(const DualVertex& r) {
        if (r == *this) return *this;
        Vertex::operator =(r);
        return *this;
    }
    DualVertex& operator =(const glm::vec3& r) {
        if (r == *this) return *this;
        Vertex::operator =(r);
        return *this;
    }
    void BuildFrom(const Mesh& mesh, const size_t cellId, bool isBoundary, size_t surfaceId = MAXID) {
        this->id = cellId;
        this->isBoundary = isBoundary;
        if (!isBoundary) {
            const Cell& cell = mesh.C.at(cellId);
            const size_t facesSize = cell.N_Fids.size();
            // Compute center x,y,z coordinates
            for (size_t i = 0; i < facesSize; i++) {
                const Face& face = mesh.F.at(cell.N_Fids.at(i));
                glm::vec3 faceCenter(0.0, 0.0, 0.0);
                const size_t faceVidsSize = face.Vids.size();
                for (size_t j = 0; j < faceVidsSize; j++) {
                    const Vertex& v = mesh.V.at(face.Vids.at(j));
                    faceCenter.x += v.x;
                    faceCenter.y += v.y;
                    faceCenter.z += v.z;
                }
                x += faceCenter.x / faceVidsSize;
                y += faceCenter.y / faceVidsSize;
                z += faceCenter.z / faceVidsSize;
            }
            x /= facesSize;
            y /= facesSize;
            z /= facesSize;
            // Compute N_Eids and N_Vids
            for (size_t i = 0; i < facesSize; i++) {
                const Face& face = mesh.F.at(cell.N_Fids.at(i));
                if (!face.isBoundary) {
                    for (size_t j = 0; j < face.N_Cids.size(); j++)
                        if (face.N_Cids.at(j) != cellId) {
                            N_Vids.push_back(face.N_Cids.at(j));
                            break;
                        }
                } else {
                    size_t surfaceNodeId = mesh.C.size();
                    for (size_t j = 0; j < mesh.F.size(); j++)
                        if (mesh.F.at(j).isBoundary) {
                            if (face.id == mesh.F.at(j).id) {
                                N_Vids.push_back(surfaceNodeId);
                                break;
                            }
                            surfaceNodeId++;
                        }
                }
                N_Eids.push_back(face.id);
            }
        } else {
            const Face& face = mesh.F.at(surfaceId);
            glm::vec3 faceCenter(0.0, 0.0, 0.0);
            const size_t faceVidsSize = face.Vids.size();
            for (size_t j = 0; j < faceVidsSize; j++) {
                const Vertex& v = mesh.V.at(face.Vids.at(j));
                faceCenter.x += v.x;
                faceCenter.y += v.y;
                faceCenter.z += v.z;
            }
            x += faceCenter.x / faceVidsSize;
            y += faceCenter.y / faceVidsSize;
            z += faceCenter.z / faceVidsSize;

            N_Vids.push_back(face.N_Cids.at(0));
            N_Eids.push_back(surfaceId);
        }
    }
};

class DualEdge: public Edge {
public:
    DualEdge()
    : Edge()
    {}
    DualEdge(const DualEdge& r)
    : Edge(r)
    , vids_link(r.vids_link)
    , eids_link(r.eids_link)
    {}
    DualEdge(size_t vnum)
    : Edge(vnum)
    {}
    DualEdge(const std::vector<size_t>& Vids)
    : Edge(Vids)
    {}
    virtual ~DualEdge()
    {}
public:
    bool operator ==(const DualEdge& e) const {
        return ((Vids[0] == e.Vids[0] && Vids[1] == e.Vids[1]) || (Vids[0] == e.Vids[1] && Vids[1] == e.Vids[0]));
    }

    std::vector<size_t> vids_link;
    std::vector<size_t> eids_link;
};

class DualFace : public Face
{
public:
    DualFace()
    : Face()
    {}
    DualFace(const DualFace& r)
    : Face(r)
    {}
    DualFace(size_t vnum)
    : Face(vnum)
    {}
    DualFace(size_t vnum, size_t eNum)
    : Face(vnum, eNum)
    {}
    DualFace(const std::vector<size_t>& Vids)
    : Face(Vids)
    {}
    DualFace(const std::vector<size_t>& Vids, const std::vector<size_t> Eids)
    : Face(Vids, Eids)
    {}
    virtual ~DualFace()
    {}
};

class DualCell : public Cell
{
public:
    DualCell()
    : Cell()
    {}
    DualCell(const DualCell& r)
    : Cell(r)
    {}
    DualCell(size_t vnum)
    : Cell(vnum)
    {}
    DualCell(size_t vNum, size_t eNum, size_t fNum)
    : Cell(vNum, eNum, fNum)
    {}
    DualCell(const std::vector<size_t>& Vids)
    : Cell(Vids)
    {}
    DualCell(const std::vector<size_t>& Vids, const std::vector<size_t> Eids, const std::vector<size_t> Fids)
    : Cell(Vids, Eids, Fids)
    {}
    virtual ~DualCell()
    {}
};

class Dual
{
public:
    Dual(Mesh& mesh);
    virtual ~Dual();
    void WriteFaces() const;
    void WriteEdges() const;
private:
    Dual(Dual&);
    Dual();
public:
    void Build();
private:
    void BuildV();
    void BuildE();
    void BuildF();
    void BuildC();
public:
    std::vector<DualVertex> V;
    std::vector<DualEdge> E;
    std::vector<DualFace> F;
    std::vector<DualFace> C;
private:
    Mesh& mesh;
};

std::vector<std::vector<size_t>> getParallelEdgeIds(const Mesh& mesh, const Cell& cell);
#endif /* LIBCOTRIK_SRC_DUAL_H_ */
