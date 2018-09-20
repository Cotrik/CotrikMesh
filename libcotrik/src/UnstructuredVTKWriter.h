/*
 * UnstructuredVTKWriter.h
 *
 *  Created on: Sep 2, 2018
 *      Author: davim
 */

#ifndef LIBCOTRIK_SRC_UNSTRUCTUREDVTKWRITER_H_
#define LIBCOTRIK_SRC_UNSTRUCTUREDVTKWRITER_H_

#include "MeshFileWriter.h"
#include <map>
#include <set>

class UnstructuredVTKWriter : public MeshFileWriter {
public:
	UnstructuredVTKWriter(const std::vector<Vertex>& V, const std::vector<Cell>& C, const char* pFileName);
	UnstructuredVTKWriter(const Mesh& mesh, const char* pFileName = "tempfile");
	virtual ~UnstructuredVTKWriter();
private:
	UnstructuredVTKWriter();
	UnstructuredVTKWriter(const UnstructuredVTKWriter&);
	UnstructuredVTKWriter& operator = (const UnstructuredVTKWriter&);
public:
	void AddPointScalarFields(const std::string fieldName, const std::vector<size_t>& field);
	void AddPointScalarFields(const std::string fieldName, const std::vector<double>& field);
	void AddCellScalarFields(const std::string fieldName, const std::vector<size_t>& field);
	void AddCellScalarFields(const std::string fieldName, const std::vector<double>& field);
	void AddPointVectorFields(const std::string fieldName, const std::vector<std::vector<size_t>>& field);
	void AddPointVectorFields(const std::string fieldName, const std::vector<std::vector<double>>& field);
	void AddCellVectorFields(const std::string fieldName, const std::vector<std::vector<size_t>>& field);
	void AddCellVectorFields(const std::string fieldName, const std::vector<std::vector<double>>& field);

	void WriteFile() const;
private:
	std::vector<std::pair<std::string, std::vector<size_t>>> point_intScalarFields;
	std::vector<std::pair<std::string, std::vector<size_t>>> cell_intScalarFields;
	std::vector<std::pair<std::string, std::vector<std::vector<size_t>>>> point_intVectorFields;
	std::vector<std::pair<std::string, std::vector<std::vector<size_t>>>> cell_intVectorFields;

	std::vector<std::pair<std::string, std::vector<double>>> point_doubleScalarFields;
	std::vector<std::pair<std::string, std::vector<double>>> cell_doubleScalarFields;
	std::vector<std::pair<std::string, std::vector<std::vector<double>>>> point_doubleVectorFields;
	std::vector<std::pair<std::string, std::vector<std::vector<double>>>> cell_doubleVectorFields;
};

#endif /* LIBCOTRIK_SRC_UNSTRUCTUREDVTKWRITER_H_ */
