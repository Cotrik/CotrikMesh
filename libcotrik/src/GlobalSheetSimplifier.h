/*
 * GlobalSheetSimplifier.h
 *
 *  Created on: Jan 20, 2019
 *      Author: cotrik
 */

#ifndef LIBCOTRIK_SRC_GLOBAL_SHEET_SIMPLIFIER_H_
#define LIBCOTRIK_SRC_GLOBAL_SHEET_SIMPLIFIER_H_

#include "SheetSimplifier.h"
class GlobalSheetSimplifier : public SheetSimplifier {
public:
    GlobalSheetSimplifier(Mesh& mesh);
    virtual ~GlobalSheetSimplifier();
private:
    GlobalSheetSimplifier();
    GlobalSheetSimplifier(const GlobalSheetSimplifier&);
    GlobalSheetSimplifier& operator = (const GlobalSheetSimplifier&);
public:
    void Run(std::set<size_t>& canceledFids);
    std::set<size_t> GetCanceledEdgeIds(const std::vector<size_t>& linkEids, std::map<size_t, size_t>& canceledFaceIds);
    bool CanCollapseWithFeaturePreserved(const BaseComplexSheetQuad& baseComplexSheets, std::map<size_t, size_t>& canceledFaceIds);
};

#endif /* LIBCOTRIK_SRC_GLOBAL_SHEET_SIMPLIFIER_H_ */
