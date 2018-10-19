/*
 * BaseComplexSheet.h
 *
 *  Created on: Oct 5, 2017
 *      Author: cotrik
 */

#ifndef LIBCOTRIK_SRC_BASECOMPLEXSHEET_H_
#define LIBCOTRIK_SRC_BASECOMPLEXSHEET_H_

#include "BaseComplex.h"
#include "RefinedDual.h"
#include <unordered_set>
#include <unordered_map>

enum SheetsConnectivity {
    SheetsConnectivity_UNKNOWN = 0,
    SheetsConnectivity_NEIGHBOR,
    SheetsConnectivity_INTERSECT,
    SheetsConnectivity_NEIGHBOR_AND_INTERSECT
};

enum Relation {
    Relation_UNKNOWN = 0,
    ADJACENT,
    INTERSECTING,
    HYBRID
};


// const int SheetsConnectivity_HYBRID = SheetsConnectivity_NEIGHBOR_AND_INTERSECT;
// #define SheetsConnectivity_HYBRID SheetsConnectivity_NEIGHBOR_AND_INTERSECT

struct sheetIds_overlaps_complexity {
	std::set<size_t> sheetIds;
	std::vector<size_t> sorted_sheetIds;
	size_t overlaps;
	float complexity;
};

class BaseComplexSheet {
public:
    BaseComplexSheet(BaseComplex& baseComplex);
    virtual ~BaseComplexSheet();
private:
    BaseComplexSheet(const BaseComplexSheet&);
    BaseComplexSheet();
    BaseComplexSheet& operator = (const BaseComplexSheet&);
public:
    void Extract();
    std::unordered_set<std::vector<std::vector<size_t>>> ExtractSheetsComponentEdgeIdsSets();
    std::unordered_set<std::vector<std::vector<size_t>>> ExtractSheetsComponentFaceIdsSets();
    std::unordered_set<std::vector<std::vector<size_t>>> ExtractSheetsComponentCellIdsSets();
    void ExtractSets();
    void ExtractSheetDecompositionsAll();
    void ExtractSheetDecompositions(const bool bfs = false);
    void ExtractSheetConnectivities();
    void VerifySheetDecompositions();
    /////////////////////////////////////////////
    void ExtractMainSheets();
    std::vector<std::vector<std::vector<std::unordered_set<size_t>>>> Get_sheet_intersecting_component_ids_groups() const;
    /////////////////////////////////////////////
    void ExtractMainSheetConnectivities(int main_sheets_id = 0);
    void ComputeComplexity();
    void ComputeComplexityDrChen(int sheetid = 0);
    void ComputeComplexityUnbalancedMatrix(int mainsheetid = 0, bool write_mat = false);
    void ComputeImportance();
    void WriteSheetDecompositionsFile(const char *filename) const;
    void WriteSheetDecompositionsVTK(const char *filename) const;
    void WriteSheetsEdgesVTK(const char *filename) const;
    void WriteSheetsFacesVTK(const char *filename) const;
    void WriteSheetsCellsVTK(const char *filename) const;
    void WriteSheetCellsVTK(const char *filename, const size_t sheet_id) const;
    void WriteSheetFacesVTK(const char *filename, const size_t sheet_id) const;
    void WriteSheetEdgesVTK(const char *filename, const size_t sheet_id) const;
    void WriteSheetFacesAndEdgesVTK(const char *filename, const size_t sheet_id) const;
    void WriteSheetCellsWithConnectedSingularitiesVTK(const char *filename, const size_t sheet_id) const;
    void WriteSheetCellsConnectedSingularitiesVTK(const char *filename, const size_t sheet_id) const;
    void WriteAllSheetsCellsVTK(const char *filename_prefix) const;
    void WriteAllSheetsFacesVTK(const char *filename_prefix) const;
    void WriteAllSheetsEdgesVTK(const char *filename_prefix) const;
    void WriteAllSheetsFacesAndEdgesVTK(const char *filename_prefix) const;
    void WriteAllSheetsCellsWithConnectedSingularitiesVTK(const char *filename_prefix) const;
    void WriteAllSheetsCellsConnectedSingularitiesVTK(const char *filename_prefix) const;
    void WriteSheetCellsDualVTK(const char *filename, const size_t sheet_id) const;
    void WriteSheetCellsDualVTK(const char *filename, const size_t sheet_id, const RefinedDual& dual) const;
    void WriteSheetCellsDualVTK(const char *filename, const size_t sheet_id, const RefinedDual& dual, const std::unordered_set<size_t>& dualFaceIds) const;
    void WriteSheetCellsDualVTK(const char *filename, const size_t sheet_id, const RefinedDual& dual, const std::unordered_set<size_t>& dualFaceIds,
            const std::unordered_set<size_t>& parallelSingularEdgeIds) const;
    void WriteAllSheetsCellsDualVTK(const char *filename_prefix) const;
    void WriteSheetsConnectivitiesMatrixVTK(const char *filename) const;
    void WriteSheetsConnectivitiesMatrixMat(const char* filename) const;
    void WriteAllDominantSheetsConnectivitiesMatrixMat(const char* filename_prefix);
    void WriteDominantSheetsConnectivitiesMatrixMat(const char* filename) const;
    void GetParallelComponents(const ComponentEdge & componentEdge,
            std::vector<size_t>& sheetComponentEdgeIds, std::vector<size_t>& sheetComponentFaceIds, std::vector<size_t>& sheetComponentCellIds);
    std::vector<size_t> GetParallelComponentEdgeIds(const ComponentEdge & componentEdge);
    std::unordered_set<size_t> GetParallelEdgeIds(const size_t sheet_id) const;
    std::unordered_set<size_t> GetParallelSingularEdgeIds(const size_t sheet_id) const;
    std::unordered_set<size_t> GetSheetBoundaryFaceComponentIds(size_t sheetId) const;
    bool HasCoveredSheetComponents(size_t sheetId, const std::vector<bool>& componentCovered) const;
    std::vector<size_t> GetCoverSheetIds(size_t beginSheetId) const;
    std::vector<size_t> GetCoverSheetIdsBFS(size_t beginSheetId) const;
    std::unordered_set<size_t> GetNeighborSheetIds(size_t sheetId) const;
    size_t GetMinComponentIntersectionNeighborSheetId(size_t sheetId, const std::unordered_set<size_t>& neighborSheetIds, const std::unordered_set<size_t>& resSet) const;
    std::vector<size_t> GetMinComponentIntersectionNeighborSheetIds(size_t sheetId, const std::unordered_set<size_t>& neighborSheetIds, const std::unordered_set<size_t>& resSet) const;
    void WriteSheetNeighborsSheetIdsJS(const char* filename);
    std::unordered_set<size_t> GetDualFaceIds(const size_t sheet_id) const;
    void RemoveSheetSetsRedundancy();
    void RemoveSheetSetsRedundancy(std::vector<size_t>& coverSheetIds);
    const std::vector<std::vector<size_t>>& Get_sheets_coverSheetIds() const;
    const size_t GetNumOfSheets() const;
    bool IsAdjacent(const int component_id1, const int component_id2) const;
    std::unordered_set<size_t> GetCommonComponentCellIds(const std::unordered_set<size_t>& common_component_face_ids) const;
    std::unordered_set<size_t> GetCommonComponentFaceIds(size_t sheetid1, size_t sheetid2) const;
    std::unordered_set<size_t> GetCommonComponentCellIds(size_t sheetid1, size_t sheetid2) const;
    size_t GetNumOfIntersections(const std::unordered_set<size_t>& common_component_cell_ids) const;
    std::vector<std::unordered_set<size_t>> GetIntersectionGroups(const std::unordered_set<size_t>& common_component_cell_ids) const;
    void WriteComplexityMat(const std::vector<std::vector<float>>& M, const char* filename) const;
    void WriteDiagonalMat(const std::vector<std::vector<float>>& M, const char* filename) const;
    void WriteAdjacentMat(const std::vector<std::vector<float>>& M, const char* filename) const;
    void WriteIntersectingMat(const std::vector<std::vector<float>>& M, const char* filename) const;
    void WriteHybridMat(const std::vector<std::vector<float>>& M, const char* filename) const;
    //void WriteHybridMat(const std::vector<std::vector<float>>& M, const char* filename) const;
    int GetNumberOfSingularities(int sheet_id) const;
    bool IsSheetRedundant(size_t sheetId, const std::vector<bool>& componentCovered) const;
    bool Verify();

//private:
    BaseComplex& baseComplex;
    std::vector<std::vector<size_t>> sheets_componentEdgeIds;  // each sheet consists a list of base-complex componentEdge ids;
    std::vector<std::vector<size_t>> sheets_componentFaceIds;  // each sheet consists a list of base-complex componentFace ids;
    std::vector<std::vector<size_t>> sheets_componentCellIds;  // each sheet consists a list of base-complex componentCell ids;
    //std::vector<std::unordered_set<size_t>> sheets_hashComponentCellIds;  // each sheet consists a list of base-complex componentCell ids;

    std::vector<std::unordered_set<size_t>> cellComponent_sheetIds;
    std::unordered_map<size_t, std::unordered_set<size_t>> faceComponent_neighborSheetIds;
    std::vector<std::unordered_set<size_t>> sheets_BoundaryFaceComponentIds;
    std::vector<std::vector<size_t>> sheets_coverSheetIds; // begin with sheetId, find a list of neighbor sheets that cover all components;

    std::vector<std::vector<size_t>> sheets_connectivities;
    std::vector<std::vector<float>> sheets_connectivities_float;
    std::vector<std::vector<size_t>> all_sheets_connectivities;
    std::vector<float> sheets_importance;
    bool m_flag = false;

    std::vector<sheetIds_overlaps_complexity> socs;
};

#endif /* LIBCOTRIK_SRC_BASECOMPLEXSHEET_H_ */
