/*
 * DoubletSimplifier.h
 *
 *  Created on: Jan 20, 2019
 *      Author: cotrik
 */

#ifndef LIBCOTRIK_SRC_DOUBLET_SIMPLIFIER_H_
#define LIBCOTRIK_SRC_DOUBLET_SIMPLIFIER_H_

#include "Simplifier.h"
class DoubletSimplifier : public Simplifier {
public:
    DoubletSimplifier(Mesh& mesh);
    virtual ~DoubletSimplifier();
private:
    DoubletSimplifier();
    DoubletSimplifier(const DoubletSimplifier&);
    DoubletSimplifier& operator = (const DoubletSimplifier&);
public:
    void Run(std::set<size_t>& canceledFids);
    void RunCollective(std::set<size_t>& canceledFids);
};

#endif /* LIBCOTRIK_SRC_DOUBLET_SIMPLIFIER_H_ */
