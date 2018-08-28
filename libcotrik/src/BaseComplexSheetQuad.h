/*
 * BaseComplexSheet.h
 *
 *  Created on: Oct 5, 2017
 *      Author: cotrik
 */

#ifndef LIBCOTRIK_SRC_BASECOMPLEXSHEETQUAD_H_
#define LIBCOTRIK_SRC_BASECOMPLEXSHEETQUAD_H_

#include "BaseComplexQuad.h"
#include "RefinedDualQuad.h"
#include <unordered_set>
#include <unordered_map>

enum SheetType {
    SIMPLE = 0,
    NON_SIMPLE
};

class BaseComplexSheetQuad
{
public:
    BaseComplexSheetQuad(BaseComplexQuad& baseComplex);
    virtual ~BaseComplexSheetQuad();
private:
    BaseComplexSheetQuad(const BaseComplexSheetQuad&);
    BaseComplexSheetQuad();
    BaseComplexSheetQuad& operator = (const BaseComplexSheetQuad&);
public:
    void Extract();
    std::unordered_set<std::vector<std::vector<size_t>>> ExtractSheetsComponentEdgeIdsSets();
    std::unordered_set<std::vector<std::vector<size_t>>> ExtractSheetsComponentFaceIdsSets();
    std::unordered_set<std::vector<std::vector<size_t>>> ExtractSheetsComponentCellIdsSets();
    void ExtractSets();
    void ExtractSheetDecompositions(const bool bfs = false);
    void WriteSheetsEdgesVTK(const char *filename) const;
    void WriteSheetsFacesVTK(const char *filename) const;

    void ExtractSheetConnectivities();
//    void ExtractMainSheetConnectivities();
    void ExtractMainSheetConnectivities(int main_sheets_id = 0);
    void ComputeComplexity();
    void ComputeImportance();
    void WriteSheetsConnectivitiesMatrixVTK(const char *filename) const;
    void WriteSheetsConnectivitiesMatrixMat(const char *filename) const;

    //void WriteSheetsCellsVTK(const char *filename) const;
    //void WriteSheetCellsVTK(const char *filename, const size_t sheet_id) const;
    void WriteSheetFacesVTK(const char *filename, const size_t sheet_id) const;
    void WriteSheetEdgesVTK(const char *filename, const size_t sheet_id) const;
    void WriteSheetFacesAndEdgesVTK(const char *filename, const size_t sheet_id) const;
    //void WriteAllSheetsCellsVTK(const char *filename_prefix) const;
    void WriteAllSheetsFacesVTK(const char *filename_prefix) const;
    void WriteAllSheetsEdgesVTK(const char *filename_prefix) const;
    void WriteAllSheetsFacesAndEdgesVTK(const char *filename_prefix) const;
    void WriteSelfIntersectingEdges(const char *filename) const;
    //void WriteSheetFacesDualVTK(const char *filename, const size_t sheet_id) const;
    void WriteSheetFacesDualVTK(const char *filename, const size_t sheet_id, const RefinedDualQuad& dual) const;
    void WriteSheetFacesDualVTK(const char *filename, const size_t sheet_id, const RefinedDualQuad& dual, const std::unordered_set<size_t>& dualEdgeIds) const;
    void WriteSheetFacesDualVTK(const char *filename, const size_t sheet_id, const RefinedDualQuad& dual, const std::unordered_set<size_t>& dualEdgeIds,
            const std::unordered_set<size_t>& singularVertexIds) const;
    void WriteAllSheetsFacesDualVTK(const char *filename_prefix) const;
    void WriteSheetDecompositionsFile(const char *filename) const;
    void WriteSheetDecompositionsDuaVTK(const char *filename, const int scalar = 0, const int main_sheets_id = 0) const;
    void GetParallelComponents(const ComponentEdge & componentEdge,
            std::vector<size_t>& sheetComponentEdgeIds, std::vector<size_t>& sheetComponentFaceIds, std::vector<size_t>& sheetComponentCellIds);
    std::vector<size_t> GetParallelComponentEdgeIds(const ComponentEdge & componentEdge);
    std::unordered_set<size_t> GetParallelEdgeIds(const size_t sheet_id) const;
    std::unordered_set<size_t> GetParallelSingularEdgeIds(const size_t sheet_id) const;
    std::unordered_set<size_t> GetSheetBoundaryFaceComponentIds(size_t sheetId) const;
    std::unordered_set<size_t> GetSheetBoundaryEdgeComponentIds(size_t sheetId) const;
    bool HasCoveredSheetComponents(size_t sheetId, const std::vector<bool>& componentCovered) const;
    bool IsSheetRedundant(size_t sheetId, const std::vector<bool>& componentCovered) const;
    std::vector<size_t> GetCoverSheetIds(size_t beginSheetId) const;
    std::vector<size_t> GetCoverSheetIdsBFS(size_t beginSheetId) const;
    std::unordered_set<size_t> GetNeighborSheetIds(size_t sheetId) const;
    size_t GetMinComponentIntersectionNeighborSheetId(size_t sheetId, const std::unordered_set<size_t>& neighborSheetIds, const std::unordered_set<size_t>& resSet) const;
    std::vector<size_t> GetMinComponentIntersectionNeighborSheetIds(size_t sheetId, const std::unordered_set<size_t>& neighborSheetIds, const std::unordered_set<size_t>& resSet) const;
    void WriteSheetNeighborsSheetIdsJS(const char* filename) const;
    std::unordered_set<size_t> GetDualEdgeIds(const size_t sheet_id) const;
    void RemoveSheetSetsRedundancy();
    void RemoveSheetSetsRedundancy(std::vector<size_t>& coverSheetIds);
    const std::vector<std::vector<size_t>>& Get_sheets_coverSheetIds() const;
    const size_t GetNumOfSheets() const;
    bool IsAdjacent(const int component_id1, const int component_id2) const;
    void ComputeComplexityDrChen(int sheetid = 0);
    std::unordered_set<size_t> GetCommonComponentFaceIds(const std::unordered_set<size_t>& common_component_edge_ids) const;
    std::unordered_set<size_t> GetCommonComponentEdgeIds(size_t sheetid1, size_t sheetid2) const;
    std::unordered_set<size_t> GetCommonComponentFaceIds(size_t sheetid1, size_t sheetid2) const;
    size_t GetNumOfIntersections(const std::unordered_set<size_t>& common_component_cell_ids) const;
    void WriteComplexityMat(const std::vector<std::vector<float>>& M, const char* filename) const;
    void WriteDiagonalMat(const std::vector<std::vector<float>>& M, const char* filename) const;
    void WriteAdjacentMat(const std::vector<std::vector<float>>& M, const char* filename) const;
    void WriteIntersectingMat(const std::vector<std::vector<float>>& M, const char* filename) const;
    void WriteHybridMat(const std::vector<std::vector<float>>& M, const char* filename) const;
    int GetNumberOfSingularities(int sheet_id) const;
    float ComputeComplexityUnbalancedMatrix(int mainsheetid = 0);
    void ExtractSheetDecompositionsAll();
    SheetType GetSheetType();
    void VerifySheetDecompositions();
private:
    BaseComplexQuad& baseComplex;
    std::vector<std::vector<size_t>> sheets_componentEdgeIds;  // each sheet consists a list of base-complex componentEdge ids;
    std::vector<std::vector<size_t>> sheets_componentFaceIds;  // each sheet consists a list of base-complex componentFace ids;
    std::vector<std::vector<size_t>> sheets_componentCellIds;  // each sheet consists a list of base-complex componentCell ids;
    std::vector<std::unordered_set<size_t>> sheets_hashComponentFaceIds;  // each sheet consists a list of base-complex componentCell ids;
    //std::vector<SheetType> sheets_types;
    SheetType sheetType = SIMPLE;

    std::vector<std::unordered_set<size_t>> faceComponent_sheetIds;
    std::unordered_map<size_t, std::unordered_set<size_t>> edgeComponent_neighborSheetIds;
    std::vector<std::unordered_set<size_t>> sheets_BoundaryEdgeComponentIds;
    std::vector<std::vector<size_t>> sheets_coverSheetIds; // begin with sheetId, find a list of neighbor sheets that cover all components;

    std::vector<std::vector<size_t>> sheets_connectivities;
    std::vector<std::vector<float>> sheets_connectivities_float;
    std::vector<std::vector<size_t>> all_sheets_connectivities;
    std::vector<float> sheets_importance;
    bool m_flag = false;
};

#endif /* LIBCOTRIK_SRC_BASECOMPLEXSHEETQUAD_H_ */
