/*
 * UnstructuredVTKWriter.cpp
 *
 *  Created on: Sep 2, 2018
 *      Author: davim
 */

#include "UnstructuredVTKWriter.h"

void WriteIntScalarField(std::ofstream& ofs, const std::string& fieldName, const std::vector<size_t>& field);
void WriteDoubleScalarField(std::ofstream& ofs, const std::string& fieldName, const std::vector<double>& field);
void WriteIntVectorField(std::ofstream& ofs, const std::string& fieldName, const std::vector<std::vector<size_t>>& field);
void WriteDoubleVectorField(std::ofstream& ofs, const std::string& fieldName, const std::vector<std::vector<double>>& field);

UnstructuredVTKWriter::UnstructuredVTKWriter(const std::vector<Vertex>& V, const std::vector<Cell>& C, const char* pFileName)
: MeshFileWriter(V, C, pFileName, HEXAHEDRA)
{

}

UnstructuredVTKWriter::UnstructuredVTKWriter(const Mesh& mesh, const char* pFileName/* = "tempfile"*/)
: MeshFileWriter(mesh, pFileName)
{

}

UnstructuredVTKWriter::~UnstructuredVTKWriter() {
	// TODO Auto-generated destructor stub
}

void UnstructuredVTKWriter::WriteFile() const {
	const std::vector<Vertex>& V = m_mesh.V;
	const std::vector<Cell>& C = m_mesh.C;
	const size_t vnum = V.size();
	const size_t cnum = C.size();
	// ------ Head
	std::ofstream ofs(m_strFileName.c_str());
	ofs << "# vtk DataFile Version 3.0\n"
		<< m_strFileName.c_str() << "\n"
		<< "ASCII\n\n"
		<< "DATASET UNSTRUCTURED_GRID\n";
	// ------ Points ------
	ofs << "POINTS " << vnum << " double\n";
	for (auto & v : V)
		ofs << v.x << " " << v.y << " " << v.z << "\n";
	// ------ Cells ------
	size_t total_number_of_vertices = 0;
	for (auto& cell : C)
		total_number_of_vertices += 1 + cell.Vids.size();
	ofs << "CELLS " << cnum << " " << total_number_of_vertices << "\n";
	for (auto& c : C) {
		ofs << c.Vids.size();
		for (auto vid : c.Vids)
			ofs << " " << vid;
		ofs << "\n";
	}
	// ------ CellTypes ------
	ofs << "CELL_TYPES " << cnum << "\n";
	for (auto& c : C) {
		ofs << c.cellType << "\n";
	}
	// ------ PointData ------
	if (!point_intScalarFields.empty() || !point_doubleScalarFields.empty()
		|| !point_intVectorFields.empty() || !point_doubleVectorFields.empty()) {
	    ofs << "POINT_DATA " << vnum << "\n";
	    for (auto & p : point_intScalarFields) WriteIntScalarField(ofs, p.first, p.second);
	    for (auto & p : point_doubleScalarFields) WriteDoubleScalarField(ofs, p.first, p.second);
	    for (auto & p : point_intVectorFields) WriteIntVectorField(ofs, p.first, p.second);
	    for (auto & p : point_doubleVectorFields) WriteDoubleVectorField(ofs, p.first, p.second);
	}
	// ------ CellData ------
	if (!cell_intScalarFields.empty() || !cell_doubleScalarFields.empty()
		|| !cell_intVectorFields.empty() || !cell_doubleVectorFields.empty()) {
	    ofs << "CELL_DATA " << cnum << "\n";
	    for (auto & p : cell_intScalarFields) WriteIntScalarField(ofs, p.first, p.second);
	    for (auto & p : cell_doubleScalarFields) WriteDoubleScalarField(ofs, p.first, p.second);
	    for (auto & p : cell_intVectorFields) WriteIntVectorField(ofs, p.first, p.second);
	    for (auto & p : cell_doubleVectorFields) WriteDoubleVectorField(ofs, p.first, p.second);
	}
}

void WriteIntScalarField(std::ofstream& ofs, const std::string& fieldName, const std::vector<size_t>& field) {
	ofs << "SCALARS " << fieldName << " int 1\n"
		<< "LOOKUP_TABLE default\n";
	for (auto& scalar : field)
		ofs << scalar << "\n";
}

void WriteDoubleScalarField(std::ofstream& ofs, const std::string& fieldName, const std::vector<double>& field) {
	ofs << "SCALARS " << fieldName << " double 1\n"
		<< "LOOKUP_TABLE default\n";
	for (auto& scalar : field)
		ofs << scalar << "\n";
}

// Only for int vector3
void WriteIntVectorField(std::ofstream& ofs, const std::string& fieldName, const std::vector<std::vector<size_t>>& field) {
	ofs << "VECTORS " << fieldName << " int\n";
	for (auto& vec : field) {
		for (auto v : vec)
			ofs << v << " ";
		ofs << "\n";
	}
}

// Only for double vector3
void WriteDoubleVectorField(std::ofstream& ofs, const std::string& fieldName, const std::vector<std::vector<double>>& field) {
	ofs << "VECTORS " << fieldName << " double\n";
	for (auto& vec : field) {
		for (auto v : vec)
			ofs << v << " ";
		ofs << "\n";
	}
}

void UnstructuredVTKWriter::AddPointScalarFields(const std::string fieldName, const std::vector<size_t>& field) {
    point_intScalarFields.push_back(std::pair<std::string, std::vector<size_t>>(fieldName, field));
}

void UnstructuredVTKWriter::AddPointScalarFields(const std::string fieldName, const std::vector<double>& field) {
    point_doubleScalarFields.push_back(std::pair<std::string, std::vector<double>>(fieldName, field));
}

void UnstructuredVTKWriter::AddCellScalarFields(const std::string fieldName, const std::vector<size_t>& field) {
    cell_intScalarFields.push_back(std::pair<std::string, std::vector<size_t>>(fieldName, field));
}

void UnstructuredVTKWriter::AddCellScalarFields(const std::string fieldName, const std::vector<double>& field) {
    cell_doubleScalarFields.push_back(std::pair<std::string, std::vector<double>>(fieldName, field));
}

void UnstructuredVTKWriter::AddPointVectorFields(const std::string fieldName, const std::vector<std::vector<size_t>>& field) {
    point_intVectorFields.push_back(std::pair<std::string, std::vector<std::vector<size_t>>>(fieldName, field));
}

void UnstructuredVTKWriter::AddPointVectorFields(const std::string fieldName, const std::vector<std::vector<double>>& field) {
    point_doubleVectorFields.push_back(std::pair<std::string, std::vector<std::vector<double>>>(fieldName, field));
}

void UnstructuredVTKWriter::AddCellVectorFields(const std::string fieldName, const std::vector<std::vector<size_t>>& field) {
    cell_intVectorFields.push_back(std::pair<std::string, std::vector<std::vector<size_t>>>(fieldName, field));
}

void UnstructuredVTKWriter::AddCellVectorFields(const std::string fieldName, const std::vector<std::vector<double>>& field) {
    cell_doubleVectorFields.push_back(std::pair<std::string, std::vector<std::vector<double>>>(fieldName, field));
}
