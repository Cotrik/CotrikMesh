#ifndef Mesh_File_Writer_H_
#define Mesh_File_Writer_H_

#include <vector>
#include <set>
#include <string>
#include <fstream>

#include "Mesh.h"
#include "FeatureLine.h"
class MeshFileWriter {
public:
	MeshFileWriter(const Mesh& mesh, const char* pFileName = "tempfile");
	MeshFileWriter(const std::vector<Vertex>& v, const std::vector<Cell>& c, const char* pFileName, const ElementType cellType = HEXAHEDRA);
	MeshFileWriter(const std::vector<Vertex>& v, const std::vector<Face>& f, const char* pFileName, const ElementType cellType = QUAD);
	MeshFileWriter();
	~MeshFileWriter();

	void WriteMeshBoundaryInfo(const char* pFileName = NULL);
	void WriteCellData(const std::vector<int>& cellData, const char* dataName = "cell_scalar", bool first = true);
	void WriteCellData(const std::vector<float>& cellData, const char* dataName = "cell_scalar", bool first = true);
	void WriteCellData(const std::vector<double>& cellData, const char* dataName = "cell_scalar", bool first = true);
	void WritePointData(const std::vector<int>& pointData, const char* dataName = "point_scalar", bool first = true);
	void WritePointData(const std::vector<float>& pointData, const char* dataName = "point_scalar", bool first = true);
	void WritePointData(const std::vector<double>& pointData, const char* dataName = "point_scalar", bool first = true);
	void WriteEdgesVtk();
	void WriteFramesVtk();
	void WriteSharpEdgesVtk();
	void WriteCornersVtk();
	void WriteVertexFeatureVtk();
	void WriteEdgesVtk(const std::vector<size_t>& edgeIds);
    void WriteEdgesVtk(const std::set<size_t>& edgeIds);
	void WriteEdgesVtk(const std::vector<std::set<size_t>>& edgeIds);
	void WriteEdgesVtk(const std::vector<std::vector<size_t>>& edgeIds);
	void WriteFacesVtk();
	void WriteFacesVtk(const std::vector<size_t>& faceIds);
    void WriteFacesVtk(const std::set<size_t>& faceIds);
    void WriteFacesVtk(const std::vector<std::set<size_t>>& faceIds);
	void WriteFacesVtk(const std::vector<std::vector<size_t>>& faceIds);
	void WriteFaceSegmentsVtk();	
	void WriteCellsVtk(const std::vector<size_t>& cellIds);
	void WriteCellsVtk(const std::set<size_t>& cellIds);
	void WriteCellsVtu(const std::vector<size_t>& cellIds);
	void WriteCellsVtu(const std::vector<size_t>& cellIds, const std::vector<int>& cellData, const char* dataName);
	void WriteVerticesVtk();
	void WriteVerticesVtk(const std::vector<size_t>& vertexIds);
    void WriteVerticesVtk(const std::set<size_t>& vertexIds);
	void WriteVerticesVtk(const std::vector<std::vector<size_t>>& vertexIds);
	void WriteLinksVtk(const std::vector<std::vector<size_t>>& linkVids);
	void WriteSurfaceOff();
	void WriteMatrixVTK(const std::vector<std::vector<size_t>>& m) const;
	void WriteFeatureLinesVtk(const std::vector<FeatureLine>& featureLines);
	void WriteCellsWithFeatureLinesVtk(const std::vector<FeatureLine>& featureLines);
    void WriteCellsWithFeatureVtk(const std::vector<FeatureLine>& featureLines, const std::vector<size_t>& featureVids);

public:
	void WriteFile();
	void WriteMeshFile();
	void WriteVtkFile();
	void WriteVtuFile();
	void WriteVtkPolyDataFile();
	void WriteOffFile();
	void WriteObjFile();
	void WriteStlFile();
	void WritePlyFile();
	void WriteSingularitiesVtkFile();
    void SetFixFlag(bool bFixed = true);
	void FixMesh();

protected:
	std::string m_strFileName;
	Mesh m_mesh;
	bool m_bFixed;
};

#endif // Mesh_File_Writer_H_
