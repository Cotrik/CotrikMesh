/*
 * SplitSimplifier.h
 *
 *  Created on: Jan 20, 2019
 *      Author: cotrik
 */

#ifndef LIBCOTRIK_SRC_SPLITSIMPLIFIER_H_
#define LIBCOTRIK_SRC_SPLITSIMPLIFIER_H_

#include "Simplifier.h"
class SplitSimplifier : public Simplifier {
public:
    SplitSimplifier(Mesh& mesh);
    virtual ~SplitSimplifier();
private:
    SplitSimplifier();
    SplitSimplifier(const SplitSimplifier&);
    SplitSimplifier& operator = (const SplitSimplifier&);
public:
    void Run();
};

#endif /* LIBCOTRIK_SRC_SPLITSIMPLIFIER_H_ */
