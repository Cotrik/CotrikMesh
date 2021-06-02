/*
 * DiagnalCollapseSimplifier.h
 *
 *  Created on: Jan 20, 2019
 *      Author: cotrik
 */

#ifndef LIBCOTRIK_SRC_DIAGNAL_COLLAPSE_SIMPLIFIER_H_
#define LIBCOTRIK_SRC_DIAGNAL_COLLAPSE_SIMPLIFIER_H_

#include "Simplifier.h"
class DiagnalCollapseSimplifier : public Simplifier {
public:
    DiagnalCollapseSimplifier(Mesh& mesh);
    virtual ~DiagnalCollapseSimplifier();
private:
    DiagnalCollapseSimplifier();
    DiagnalCollapseSimplifier(const DiagnalCollapseSimplifier&);
    DiagnalCollapseSimplifier& operator = (const DiagnalCollapseSimplifier&);
public:
    void Run1(std::set<size_t>& canceledFids);
    void Run2(std::set<size_t>& canceledFids);
    void Run3(std::set<size_t>& canceledFids);
    void Run(std::set<size_t>& canceledFids);
    void RunCollective(std::set<size_t>& canceledFids);
    void CollapseSinglets(std::set<size_t>& canceledFids);
private:
    void CollapseDiagnal(const Vertex& v0, const Vertex& v2);
    bool CanCollapseDiagnal(const Vertex& v0, const Vertex& v2);
};

#endif /* LIBCOTRIK_SRC_DIAGNAL_COLLAPSE_SIMPLIFIER_H_ */
