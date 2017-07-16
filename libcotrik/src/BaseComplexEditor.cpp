/*
 * BaseComplexEditor.cpp
 *
 *  Created on: Jul 10, 2017
 *      Author: cotrik
 */

#include "BaseComplexEditor.h"

BaseComplexEditor::BaseComplexEditor(const Mesh& mesh, BaseComplex& baseComplex)
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

const Mesh& BaseComplexEditor::GetMesh() const
{
    return mesh;
}

BaseComplex& BaseComplexEditor::GetBaseComplex() const
{
    return m_baseComplex;
}
