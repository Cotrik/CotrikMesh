/*
 * BaseComplexEditor.h
 *
 *  Created on: Jul 10, 2017
 *      Author: cotrik
 */

#ifndef LIBCOTRIK_SRC_BASECOMPLEXEDITOR_H_
#define LIBCOTRIK_SRC_BASECOMPLEXEDITOR_H_

#include "BaseComplex.h"

class BaseComplexEditor
{
public:
    BaseComplexEditor(const Mesh& mesh, BaseComplex& baseComplex);
    BaseComplexEditor(const BaseComplexEditor& rhs);
    virtual ~BaseComplexEditor();
private:
    BaseComplexEditor();
public:
    const Mesh& GetMesh() const;
    BaseComplex& GetBaseComplex() const;
private:
    const Mesh& mesh;
    BaseComplex& m_baseComplex;
};

#endif /* LIBCOTRIK_SRC_BASECOMPLEXEDITOR_H_ */
