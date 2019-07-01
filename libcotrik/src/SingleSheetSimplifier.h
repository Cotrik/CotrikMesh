/*
 * SingleSheetSimplifier.h
 *
 *  Created on: Jan 20, 2019
 *      Author: cotrik
 */

#ifndef LIBCOTRIK_SRC_SINGLE_SHEET_SIMPLIFIER_H_
#define LIBCOTRIK_SRC_SINGLE_SHEET_SIMPLIFIER_H_

#include "SheetSimplifier.h"
class SingleSheetSimplifier : public SheetSimplifier {
public:
    SingleSheetSimplifier(Mesh& mesh);
    virtual ~SingleSheetSimplifier();
private:
    SingleSheetSimplifier();
    SingleSheetSimplifier(const SingleSheetSimplifier&);
    SingleSheetSimplifier& operator = (const SingleSheetSimplifier&);
public:
    void Run(std::set<size_t>& canceledFids);
    std::set<size_t> GetCanceledEdgeIds(const std::vector<size_t>& linkEids, std::map<size_t, size_t>& canceledFaceIds);
    bool CanCollapseWithFeaturePreserved(const BaseComplexSheetQuad& baseComplexSheets, std::map<size_t, size_t>& canceledFaceIds);
};

#endif /* LIBCOTRIK_SRC_SINGLE_SHEET_SIMPLIFIER_H_ */
