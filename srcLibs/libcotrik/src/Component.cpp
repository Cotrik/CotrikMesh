/*
 * Component.cpp
 *
 *  Created on: May 31, 2017
 *      Author: cotrik
 */

#include "Component.h"

Component::Component(const Mesh& mesh)
: mesh(mesh)
{
    // TODO Auto-generated constructor stub

}

Component::Component(const Component& component)
: Vids(component.Vids)
, Eids(component.Eids)
, Fids(component.Fids)
, Cids(component.Cids)
, cornerVids(component.cornerVids)
, boundaryEids(component.boundaryEids)
, boundaryFids(component.boundaryFids)
, mesh(component.GetMesh())
{

}

Component::~Component()
{
    // TODO Auto-generated destructor stub
}

const Mesh& Component::GetMesh() const
{
    return mesh;
}

void Component::Extract()
{

}
