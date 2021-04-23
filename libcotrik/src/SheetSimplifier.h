/*
 * SheetSimplifier.h
 *
 *  Created on: Jan 20, 2019
 *      Author: cotrik
 */

#ifndef LIBCOTRIK_SRC_SHEETSIMPLIFIER_H_
#define LIBCOTRIK_SRC_SHEETSIMPLIFIER_H_

#include "Simplifier.h"
class SheetSimplifier : public Simplifier {
public:
    SheetSimplifier(Mesh& mesh);
    virtual ~SheetSimplifier();
private:
    SheetSimplifier();
    SheetSimplifier(const SheetSimplifier&);
    SheetSimplifier& operator = (const SheetSimplifier&);
public:
    void Run(std::set<size_t>& canceledFids);
    std::set<size_t> GetAllParallelEdgeIds(const size_t eid);
    std::set<size_t> GetCanceledEdgeIds(const BaseComplexSheetQuad& baseComplexSheets, std::map<size_t, size_t>& canceledFaceIds,
            size_t sheetId);
    bool CanCollapseWithFeaturePreserved(const BaseComplexSheetQuad& baseComplexSheets, std::map<size_t, size_t>& canceledFaceIds,
            size_t sheetId);
    void CollapseWithFeaturePreserved(std::unordered_map<size_t, size_t>& key_edgeId, std::unordered_map<std::string, size_t>& key_faceId,
        std::map<size_t, size_t>& canceledFaceIds, std::set<size_t>& canceledEdgeIds);
    void ExtractAndCollapse(std::set<size_t>& canceledFids);
    void GetSheetsProvisions(std::vector<std::set<size_t>>& canceledEdgesIds, std::vector<std::map<size_t, size_t>>& canceledFacesIds);
    void CollapseSelectedSheets(std::set<size_t>& canceledFids, std::vector<std::set<size_t>>& canceledEdgesIds, std::vector<std::map<size_t, size_t>>& canceledFacesIds);
};

#endif /* LIBCOTRIK_SRC_SHEETSIMPLIFIER_H_ */
