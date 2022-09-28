/*
 * Component.h
 *
 *  Created on: May 31, 2017
 *      Author: cotrik
 */

#ifndef LIBCOTRIK_SRC_COMPONENT_H_
#define LIBCOTRIK_SRC_COMPONENT_H_

#include "Mesh.h"

class Component
{
public:
    Component(const Mesh& mesh);
    Component(const Component& component);
    virtual ~Component();
private:
    Component();

public:
    const Mesh& GetMesh() const;
    void Extract();
    void ExtractSingularNodes();
    void ExtractSingularEdges();


public:
    std::vector<size_t> Vids;  // vertex ids of mesh
    std::vector<size_t> Eids;  // edge ids of mesh
    std::vector<size_t> Fids;  // face ids of mesh
    std::vector<size_t> Cids;  // cell ids of mesh

    std::vector<size_t> cornerVids;  // 8 vertex ids of this component, each id is a vertex id of the mesh
    std::vector<std::vector<size_t> > boundaryEids;  // 12 edge ids of this component, each edge contains a link set of edge ids of mesh
    std::vector<std::vector<size_t> > boundaryFids;  // 6 face ids of this component, each face contains a set of face ids of mesh
private:
    const Mesh& mesh;
};

#endif /* LIBCOTRIK_SRC_COMPONENT_H_ */
