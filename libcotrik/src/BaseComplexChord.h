/*
 * BaseComplexChord.h
 *
 *  Created on: Oct 5, 2017
 *      Author: cotrik
 */

#ifndef LIBCOTRIK_SRC_BASECOMPLEXCHORD_H_
#define LIBCOTRIK_SRC_BASECOMPLEXCHORD_H_

#include "BaseComplex.h"

class BaseComplexChord
{
public:
    BaseComplexChord(BaseComplex& baseComplex);
    virtual ~BaseComplexChord();
private:
    BaseComplexChord(const BaseComplexChord&);
    BaseComplexChord();
    BaseComplexChord& operator = (const BaseComplexChord&);
public:
    void Extract();
    void WriteChordsEdgesVTK(const char *filename) const;
    void WriteChordsFacesVTK(const char *filename) const;
    void WriteChordsCellsVTK(const char *filename) const;
    void WriteChordCellsVTK(const char *filename, const size_t sheet_id) const;
    void WriteChordFacesVTK(const char *filename, const size_t sheet_id) const;
    void WriteChordEdgesVTK(const char *filename, const size_t sheet_id) const;
    void WriteChordFacesAndEdgesVTK(const char *filename, const size_t sheet_id) const;
    void WriteChordCurvesVTK(const char *filename, const size_t sheet_id) const;
    void WriteChordFramesVTK(const char *filename, const size_t sheet_id) const;
    void WriteAllChordsCellsVTK(const char *filename_prefix) const;
    void WriteAllChordsFacesVTK(const char *filename_prefix) const;
    void WriteAllChordsEdgesVTK(const char *filename_prefix) const;
    void WriteAllChordsFacesAndEdgesVTK(const char *filename_prefix) const;
    void WriteAllChordsCurvesVTK(const char *filename_prefix) const;
    void WriteAllChordsFramesVTK(const char *filename_prefix) const;
    void GetParallelComponents(const ComponentFace & componentFace,
            std::vector<size_t>& chordComponentEdgeIds, std::vector<size_t>& chordComponentFaceIds, std::vector<size_t>& chordComponentCellIds);
    std::vector<size_t> GetParallelComponentFaceIds(const ComponentFace & componentFace) const;
    std::vector<size_t> GetLinkedChordComponentFaceIds(const std::vector<size_t>& chordComponentFaceIds) const;
    std::vector<size_t> GetLinkedChordComponentCellIds(const std::vector<size_t>& chordComponentCellIds, const std::vector<size_t>& linkedChordComponentFaceIds) const;
    std::vector<glm::vec3> GetLinkedChordCurveVertices(const std::vector<size_t>& linkedChordComponentFaceIds, const std::vector<size_t>& linkedChordComponentCellIds) const;
    glm::vec3 GetCenter(const ComponentFace& c) const;
    glm::vec3 GetCenter(const ComponentCell& c) const;
private:
    BaseComplex& baseComplex;
    std::vector<std::vector<size_t>> chords_componentEdgeIds;  // each sheet consists a list of base-complex componentEdge ids;
    std::vector<std::vector<size_t>> chords_componentFaceIds;  // each sheet consists a list of base-complex componentFace ids;
    std::vector<std::vector<size_t>> chords_componentCellIds;  // each sheet consists a list of base-complex componentCell ids;
};

#endif /* LIBCOTRIK_SRC_BASECOMPLEXCHORD_H_ */
