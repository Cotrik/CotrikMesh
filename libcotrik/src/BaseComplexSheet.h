/*
 * BaseComplexSheet.h
 *
 *  Created on: Oct 5, 2017
 *      Author: cotrik
 */

#ifndef LIBCOTRIK_SRC_BASECOMPLEXSHEET_H_
#define LIBCOTRIK_SRC_BASECOMPLEXSHEET_H_

#include "BaseComplex.h"
#include <unordered_set>

class BaseComplexSheet
{
public:
    BaseComplexSheet(BaseComplex& baseComplex);
    virtual ~BaseComplexSheet();
private:
    BaseComplexSheet(const BaseComplexSheet&);
    BaseComplexSheet();
    BaseComplexSheet& operator = (const BaseComplexSheet&);
public:
    void Extract();
//    std::unordered_set<std::vector<std::vector<size_t>>> ExtractSheetsComponentEdgeIdsSets();
//    std::unordered_set<std::vector<std::vector<size_t>>> ExtractSheetsComponentFaceIdsSets();
//    std::unordered_set<std::vector<std::vector<size_t>>> ExtractSheetsComponentCellIdsSets();
//    void ExtractSets();
    void WriteSheetsEdgesVTK(const char *filename) const;
    void WriteSheetsFacesVTK(const char *filename) const;
    void WriteSheetsCellsVTK(const char *filename) const;
    void WriteSheetCellsVTK(const char *filename, const size_t sheet_id) const;
    void WriteSheetFacesVTK(const char *filename, const size_t sheet_id) const;
    void WriteSheetEdgesVTK(const char *filename, const size_t sheet_id) const;
    void WriteSheetFacesAndEdgesVTK(const char *filename, const size_t sheet_id) const;
    void WriteAllSheetsCellsVTK(const char *filename_prefix) const;
    void WriteAllSheetsFacesVTK(const char *filename_prefix) const;
    void WriteAllSheetsEdgesVTK(const char *filename_prefix) const;
    void WriteAllSheetsFacesAndEdgesVTK(const char *filename_prefix) const;
    void WriteSheetCellsDualVTK(const char *filename, const size_t sheet_id) const;
    void WriteAllSheetsCellsDualVTK(const char *filename_prefix) const;
    void GetParallelComponents(const ComponentEdge & componentEdge,
            std::vector<size_t>& sheetComponentEdgeIds, std::vector<size_t>& sheetComponentFaceIds, std::vector<size_t>& sheetComponentCellIds);
    std::vector<size_t> GetParallelComponentEdgeIds(const ComponentEdge & componentEdge);
    std::unordered_set<size_t> GetParallelEdgeIds(const size_t sheet_id) const;
private:
    BaseComplex& baseComplex;
    std::vector<std::vector<size_t>> sheets_componentEdgeIds;  // each sheet consists a list of base-complex componentEdge ids;
    std::vector<std::vector<size_t>> sheets_componentFaceIds;  // each sheet consists a list of base-complex componentFace ids;
    std::vector<std::vector<size_t>> sheets_componentCellIds;  // each sheet consists a list of base-complex componentCell ids;
};

#endif /* LIBCOTRIK_SRC_BASECOMPLEXSHEET_H_ */
