/*
 * CollapseSimplifier.h
 *
 *  Created on: Jan 20, 2019
 *      Author: cotrik
 */

#ifndef LIBCOTRIK_SRC_COLLAPSESIMPLIFIER_H_
#define LIBCOTRIK_SRC_COLLAPSESIMPLIFIER_H_

#include "Simplifier.h"
class CollapseSimplifier : public Simplifier {
public:
    CollapseSimplifier(Mesh& mesh);
    virtual ~CollapseSimplifier();
private:
    CollapseSimplifier();
    CollapseSimplifier(const CollapseSimplifier&);
    CollapseSimplifier& operator = (const CollapseSimplifier&);
public:
    void Run();
};

#endif /* LIBCOTRIK_SRC_COLLAPSESIMPLIFIER_H_ */
