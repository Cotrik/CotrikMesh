#ifndef __Mesh_File_Reader_H__
#define __Mesh_File_Reader_H__

#include <vector>
#include <string>
#include <fstream>
#include "Mesh.h"

class MeshFileReader
{
public:
//	MeshFileReader();
	MeshFileReader(const char* pFileName);
	~MeshFileReader();

public:
	const Mesh& GetMesh() const;
	void GetScalarFields();
	void GetPointsScalarFields();

private:
	void ReadOffFile();
	void ReadMeshFile();
	void ReadVtkFile();
	void ReadObjFile();
	void ReadStlFile();
	bool HasCellType(const VTKCellType cellType) const;

private:
	std::string m_strFileName;
	Mesh m_mesh;
};

#endif // __Mesh_File_Reader_H__
