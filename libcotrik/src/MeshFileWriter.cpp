#include "MeshFileWriter.h"
#include <map>
#include <iostream>
#include <iomanip>
#include <memory>
#include <algorithm>
#include "stdio.h"
using namespace std;
//#include <boost/smart_ptr.hpp>

#include <vtkGenericDataObjectReader.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkOBJReader.h>
#include <vtkSTLReader.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkDataSetAttributes.h>
#include <vtkSmartPointer.h>
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkTetra.h>
#include <vtkCellArray.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkDataSetMapper.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkPointData.h>
#include <vtkVertexGlyphFilter.h>
MeshFileWriter::MeshFileWriter()
{

}
MeshFileWriter::~MeshFileWriter()
{

}

MeshFileWriter::MeshFileWriter(const Mesh& mesh, const char* pFileName)
: m_strFileName(pFileName)
, m_mesh(mesh)
, m_bFixed(false)
{

}

MeshFileWriter::MeshFileWriter(const std::vector<Vertex>& v, const std::vector<Cell>& c,
    const char* pFileName, const ElementType cellType/* = HEXAHEDRA*/)
: m_strFileName(pFileName)
, m_mesh(v, c, cellType)
, m_bFixed(false)
{

}

MeshFileWriter::MeshFileWriter(const std::vector<Vertex>& v, const std::vector<Face>& f,
    const char* pFileName, const ElementType cellType/* = HEXAHEDRA*/)
: m_strFileName(pFileName)
, m_mesh(v, f, cellType)
, m_bFixed(false)
{

}

void MeshFileWriter::WriteFile()
{
    if (m_bFixed) FixMesh();
    if (m_strFileName.find(".vtk") != m_strFileName.npos)       WriteVtkFile();
	else if (m_strFileName.find(".vtu") != m_strFileName.npos)  WriteVtuFile();
    else if (m_strFileName.find(".off") != m_strFileName.npos)  WriteOffFile();
    else if (m_strFileName.find(".mesh") != m_strFileName.npos) WriteMeshFile();
    else if (m_strFileName.find(".obj") != m_strFileName.npos)  WriteObjFile();
    else if (m_strFileName.find(".stl") != m_strFileName.npos)  WriteStlFile();
}

void MeshFileWriter::WriteMeshFile()
{
    const std::vector<Vertex>& V = m_mesh.V;
    const std::vector<Cell>& C = m_mesh.C;
    const size_t vnum = V.size();
    const size_t cnum = C.size();

    std::ofstream ofs(m_strFileName.c_str());
    ofs << "MeshVersionFormatted 2" << endl;
    ofs << "Dimension 3" << endl;
    ofs << "Vertices " << vnum << endl;

    for (size_t i = 0; i < vnum; i++)
        ofs << std::fixed << setprecision(7) << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << " 0" << endl;

    if (m_mesh.m_cellType == TRIANGLE) ofs << "Triangles ";
    else if (m_mesh.m_cellType == QUAD) ofs << "Quadrilaterals ";
    else if (m_mesh.m_cellType == TETRAHEDRA) ofs << "Tetrahedra ";
    else if (m_mesh.m_cellType == HEXAHEDRA) ofs << "Hexahedra ";
    ofs << cnum << std::endl;

    for (size_t i = 0; i < cnum; i++){
        for (size_t j = 0; j < C.at(i).Vids.size(); j++)
            ofs << C.at(i).Vids.at(j) + 1 << " ";
        ofs << "0" << std::endl;
    }
    ofs << "End" << std::endl;
}

void MeshFileWriter::WriteVtkFile() {
    if (m_mesh.m_cellType == POLYGON)
        return WriteVtkPolyDataFile();
    const std::vector<Vertex>& V = m_mesh.V;
    const std::vector<Cell>& C = m_mesh.C;
    const size_t vnum = V.size();
    const size_t cnum = C.size();

    std::ofstream ofs(m_strFileName.c_str());
    ofs << "# vtk DataFile Version 3.0\n"
        << m_strFileName.c_str() << "\n"
        << "ASCII\n\n"
        << "DATASET UNSTRUCTURED_GRID\n";
    ofs << "POINTS " << vnum << " double\n";
    for (size_t i = 0; i < vnum; i++)
        ofs << std::fixed << setprecision(7) << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << "\n";
    ofs << "CELLS " << cnum << " ";

    vtkIdType idType = VTK_TRIANGLE;
    if (m_mesh.m_cellType == TRIANGLE) ofs << 4*cnum << "\n";
    else if (m_mesh.m_cellType == QUAD) {idType = VTK_QUAD;  ofs << 5*cnum << "\n";}
    else if (m_mesh.m_cellType == TETRAHEDRA) {idType = VTK_TETRA; ofs << 5*cnum << "\n";}
    else if (m_mesh.m_cellType == HEXAHEDRA) {idType = VTK_HEXAHEDRON; ofs << 9*cnum << "\n";}
    else if (m_mesh.m_cellType == POLYHEDRA) {
        size_t sum = 0;
        for (auto& cell : C)
            sum += 1 + cell.Vids.size();
        ofs << sum << "\n";
    }

    for (size_t i = 0; i < cnum; i++){
        ofs << C.at(i).Vids.size();
        for (size_t j = 0; j < C.at(i).Vids.size(); j++)
            ofs << " " << C.at(i).Vids.at(j);
        ofs << "\n";
    }
    ofs << "CELL_TYPES " << cnum << "\n";
    if (m_mesh.m_cellType != POLYHEDRA)
        for (size_t i = 0; i < cnum; i++)
            ofs << idType << "\n";
    else
        for (auto cellType : m_mesh.m_cellTypes)
            ofs << cellType << "\n";
}

void MeshFileWriter::WriteVtuFile() {
	vtkSmartPointer<vtkPoints> points =	vtkSmartPointer<vtkPoints>::New();
	for (auto& v : m_mesh.V) 
		points->InsertNextPoint(v.x, v.y, v.z);
	vtkSmartPointer<vtkUnstructuredGrid> ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	ugrid->SetPoints(points);
	for (auto& c : m_mesh.C) {
		vtkSmartPointer<vtkCellArray> faces = vtkSmartPointer<vtkCellArray>::New();
		for (auto fid : c.Fids) {
			auto& f = m_mesh.F.at(fid);
			faces->InsertNextCell(f.Vids.size(), (vtkIdType*)f.Vids.data());
		}
		ugrid->InsertNextCell(VTK_POLYHEDRON, c.Vids.size(), (vtkIdType*)c.Vids.data(), c.Fids.size(), faces->GetPointer());
	}
	vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
#if VTK_MAJOR_VERSION <= 5
	writer->SetInput(ugrid);
#else
	writer->SetInputData(ugrid);
#endif
	writer->SetFileName(m_strFileName.c_str());
	writer->SetDataModeToAscii();
	writer->Update();
}

void MeshFileWriter::WriteCellsVtu(const std::vector<size_t>& cellIds) {
	vtkSmartPointer<vtkPoints> points =	vtkSmartPointer<vtkPoints>::New();
	for (auto& v : m_mesh.V) 
		points->InsertNextPoint(v.x, v.y, v.z);
	vtkSmartPointer<vtkUnstructuredGrid> ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	ugrid->SetPoints(points);
	for (auto cellId : cellIds) {
		auto& c = m_mesh.C.at(cellId);
		vtkSmartPointer<vtkCellArray> faces = vtkSmartPointer<vtkCellArray>::New();
		for (auto fid : c.Fids) {
			auto& f = m_mesh.F.at(fid);
			faces->InsertNextCell(f.Vids.size(), (vtkIdType*)f.Vids.data());
		}
		ugrid->InsertNextCell(VTK_POLYHEDRON, c.Vids.size(), (vtkIdType*)c.Vids.data(), c.Fids.size(), faces->GetPointer());
	}
	vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
#if VTK_MAJOR_VERSION <= 5
	writer->SetInput(ugrid);
#else
	writer->SetInputData(ugrid);
#endif
	writer->SetFileName(m_strFileName.c_str());
	writer->SetDataModeToAscii();
	writer->Update();
}

void MeshFileWriter::WriteCellsVtu(const std::vector<size_t>& cellIds, const std::vector<int>& cellDataField, const char* dataName) {
	vtkSmartPointer<vtkPoints> points =	vtkSmartPointer<vtkPoints>::New();
	for (auto& v : m_mesh.V) 
		points->InsertNextPoint(v.x, v.y, v.z);
	vtkSmartPointer<vtkUnstructuredGrid> ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	ugrid->SetPoints(points);
	for (auto cellId : cellIds) {
		auto& c = m_mesh.C.at(cellId);
		vtkSmartPointer<vtkCellArray> faces = vtkSmartPointer<vtkCellArray>::New();
		for (auto fid : c.Fids) {
			auto& f = m_mesh.F.at(fid);
			faces->InsertNextCell(f.Vids.size(), (vtkIdType*)f.Vids.data());
		}
		ugrid->InsertNextCell(VTK_POLYHEDRON, c.Vids.size(), (vtkIdType*)c.Vids.data(), c.Fids.size(), faces->GetPointer());
	}
	vtkSmartPointer<vtkIntArray> cellData = vtkSmartPointer<vtkIntArray>::New();
	cellData->SetName(dataName);
	cellData->SetNumberOfComponents(1);
	cellData->SetNumberOfTuples(cellIds.size());
	for (int i = 0; i < cellIds.size(); i++) {
		float rgb[1] = { cellDataField[i] };
		cellData->InsertTuple(i, rgb);
	}

	ugrid->GetCellData()->SetScalars(cellData);

	vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
#if VTK_MAJOR_VERSION <= 5
	writer->SetInput(ugrid);
#else
	writer->SetInputData(ugrid);
#endif
	writer->SetFileName(m_strFileName.c_str());
	writer->SetDataModeToAscii();
	writer->Update();
}
void MeshFileWriter::WriteVtkPolyDataFile()
{
    const std::vector<Vertex>& V = m_mesh.V;
    const std::vector<Cell>& C = m_mesh.C;
    const size_t vnum = V.size();
    const size_t cnum = C.size();

    std::ofstream ofs(m_strFileName.c_str());
    ofs << "# vtk DataFile Version 2.0" << endl
        << m_strFileName.c_str() << endl
        << "ASCII" << endl << endl
        << "DATASET POLYDATA" << endl;
    ofs << "POINTS " << vnum << " double" << endl;
    for (size_t i = 0; i < vnum; i++)
        ofs << std::fixed << setprecision(7) << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << endl;
    size_t numOfV = 0;
    for (auto& c : C)
        numOfV += 1 + c.Vids.size();
    ofs << "POLYGONS " << cnum << " " << numOfV << "\n";

    // vtkIdType idType = VTK_TRIANGLE;
    // if (m_mesh.m_cellType == TRIANGLE) ofs << 4*cnum << std::endl;
    // else if (m_mesh.m_cellType == QUAD) {idType = VTK_QUAD;  ofs << 5*cnum << std::endl;}
    // else if (m_mesh.m_cellType == TETRAHEDRA) {idType = VTK_TETRA; ofs << 5*cnum << std::endl;}
    // else if (m_mesh.m_cellType == HEXAHEDRA) {idType = VTK_HEXAHEDRON; ofs << 9*cnum << std::endl;}

    for (auto& c : C) {
        ofs << c.Vids.size();
        for (auto vid : c.Vids)
            ofs << " " << vid;
        ofs << std::endl;
    }
}
void MeshFileWriter::WriteOffFile()
{
    const std::vector<Vertex>& V = m_mesh.V;
    const std::vector<Cell>& C = m_mesh.C;
    const size_t vnum = V.size();
    const size_t cnum = C.size();

    std::ofstream ofs(m_strFileName.c_str());
    ofs << "OFF\n";
    ofs << vnum << " " << cnum << " 0\n";
    for (size_t i = 0; i < vnum; i++)
        ofs << std::fixed << setprecision(7) << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << endl;

    for (size_t i = 0; i < cnum; i++) {
        if (m_mesh.m_cellType == POLYGON) ofs << C.at(i).Vids.size();
        else if (m_mesh.m_cellType == TRIANGLE) ofs << 3;
        else if (m_mesh.m_cellType == QUAD) ofs << 4;
        else if (m_mesh.m_cellType == TETRAHEDRA) ofs << 5;
        else if (m_mesh.m_cellType == HEXAHEDRA) ofs << 10;

        for (size_t j = 0; j < C.at(i).Vids.size(); j++)
            ofs << " " << C.at(i).Vids.at(j);
        if (m_mesh.m_cellType == TETRAHEDRA) ofs << " 0\n";
        else if (m_mesh.m_cellType == HEXAHEDRA) ofs << " 0 0\n";
        else ofs << "\n";
    }
}

void MeshFileWriter::WriteObjFile()
{

}

void MeshFileWriter::WriteStlFile()
{

}

void MeshFileWriter::WriteCellData(const std::vector<int>& cellData, const char* dataName, bool first/* = true*/) {
    std::ofstream ofs(m_strFileName.c_str(), std::ios_base::app);
	if (first) { ofs << "CELL_DATA " << cellData.size() << std::endl; }
    ofs << "SCALARS " << dataName << " int 1" << std::endl
        << "LOOKUP_TABLE default" << std::endl;
    for (auto data : cellData)
        ofs << data << std::endl;
}

void MeshFileWriter::WriteCellData(const std::vector<float>& cellData, const char* dataName, bool first/* = true*/) {
	std::ofstream ofs(m_strFileName.c_str(), std::ios_base::app);
	if (first) { ofs << "CELL_DATA " << cellData.size() << std::endl; }
	ofs << "SCALARS " << dataName << " float 1" << std::endl
		<< "LOOKUP_TABLE default" << std::endl;
	for (auto data : cellData)
		ofs << data << std::endl;
}

void MeshFileWriter::WriteCellData(const std::vector<double>& cellData, const char* dataName, bool first/* = true*/) {
	std::ofstream ofs(m_strFileName.c_str(), std::ios_base::app);
	if (first) { ofs << "CELL_DATA " << cellData.size() << std::endl; }
	ofs << "SCALARS " << dataName << " double 1" << std::endl
		<< "LOOKUP_TABLE default" << std::endl;
	for (auto data : cellData)
		ofs << data << std::endl;
}

void MeshFileWriter::WritePointData(const std::vector<int>& pointData, const char* dataName, bool first/* = true*/) {
    std::ofstream ofs(m_strFileName.c_str(), std::ios_base::app);
	if (first) { ofs << "POINT_DATA " << pointData.size() << std::endl; }
    ofs << "SCALARS " << dataName << " int 1" << std::endl
        << "LOOKUP_TABLE default" << std::endl;
	for (auto data : pointData)
		ofs << data << std::endl;
}

void MeshFileWriter::WritePointData(const std::vector<float>& pointData, const char* dataName, bool first/* = true*/) {
	std::ofstream ofs(m_strFileName.c_str(), std::ios_base::app);
	if (first) { ofs << "POINT_DATA " << pointData.size() << std::endl; }
	ofs << "SCALARS " << dataName << " float 1" << std::endl
		<< "LOOKUP_TABLE default" << std::endl;
	for (auto data : pointData)
		ofs << data << std::endl;
}

void MeshFileWriter::WritePointData(const std::vector<double>& pointData, const char* dataName, bool first/* = true*/) {
	std::ofstream ofs(m_strFileName.c_str(), std::ios_base::app);
	if (first) { ofs << "POINT_DATA " << pointData.size() << std::endl; }
	ofs << "SCALARS " << dataName << " double 1" << std::endl
		<< "LOOKUP_TABLE default" << std::endl;
	for (auto data : pointData)
		ofs << data << std::endl;
}

void MeshFileWriter::SetFixFlag(bool bFixed)
{
	m_bFixed = bFixed;
}
void MeshFileWriter::FixMesh()
{
	std::vector<unsigned long> v_real_index;
	size_t c_size = 0;
    for (auto& c : m_mesh.C) {
        if (c.Vids.size() > 8) continue;
        else {
            for (auto vid : c.Vids)
                v_real_index.push_back(vid);
            ++c_size;
        }
    }

	std::sort(v_real_index.begin(), v_real_index.end());
	std::vector<unsigned long>::iterator iter = std::unique(v_real_index.begin(), v_real_index.end());
	v_real_index.resize(std::distance(v_real_index.begin(), iter));

	std::vector<Vertex> V(v_real_index.size());
	for (unsigned long i = 0; i < v_real_index.size(); i++)	{
		const Vertex& v = m_mesh.V.at(v_real_index.at(i));
		V.at(i) = v;
	}
	m_mesh.V = V;
	//////////////////////////////////////////////////////
	std::map<size_t, size_t> v_v;
	for (size_t i = 0; i < v_real_index.size(); i++)
		v_v[v_real_index.at(i)] = i;

    for (auto& c : m_mesh.C) {
        const size_t cellSize = c.Vids.size();
        if (cellSize == 8 || cellSize == 3 || cellSize == 4)
            for (auto& vid : c.Vids)
                vid = v_v[vid];
    }
}

void MeshFileWriter::WriteEdgesVtk()
{
    const std::vector<Vertex>& V = m_mesh.V;
    const std::vector<Edge>& E = m_mesh.E;

    std::ofstream ofs(m_strFileName.c_str());
    ofs << "# vtk DataFile Version 2.0" << endl
        << m_strFileName.c_str() << endl
        << "ASCII" << endl << endl
        << "DATASET POLYDATA" << endl;
    ofs << "POINTS " << V.size() << " float" << endl;
    ofs << std::fixed << setprecision(7);
    for (size_t i = 0; i < V.size(); i++)
        ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << endl;
    ofs << "LINES " << E.size() << " " << 3 * E.size() << std::endl;
    for (size_t i = 0; i < E.size(); i++)
        ofs << "2 " << E.at(i).Vids[0] << " " << E.at(i).Vids[1] << std::endl;
    //ofs << "CELL_DATA " << E.size() << std::endl
    //    << "SCALARS " << " Singularity" << " int 1" << std::endl
    //    << "LOOKUP_TABLE default" << std::endl;
    //for (size_t i = 0; i < E.size(); i++)
    //    ofs << (E.at(i).isSingularity ? 1 : 0) << std::endl;

    //ofs << "SCALARS " << " ESingularity" << " float" << std::endl
    //    << "LOOKUP_TABLE default" << std::endl;
    //for (size_t i = 0; i < E.size(); i++)
    //    ofs << E.at(i).energySingularity << std::endl;

    //ofs << "SCALARS " << " EOrthogonality" << " float" << std::endl
    //    << "LOOKUP_TABLE default" << std::endl;
    //for (size_t i = 0; i < E.size(); i++)
    //    ofs << E.at(i).energyOrthogonality << std::endl;

    //ofs << "SCALARS " << " EStraightness" << " float" << std::endl
    //    << "LOOKUP_TABLE default" << std::endl;
    //for (size_t i = 0; i < E.size(); i++)
    //    ofs << E.at(i).energyStraightness << std::endl;
}

void MeshFileWriter::WriteFramesVtk()
{
    const std::vector<Vertex>& V = m_mesh.V;
    const std::vector<Edge>& E = m_mesh.E;

    std::ofstream ofs(m_strFileName.c_str());
    ofs << "# vtk DataFile Version 2.0" << endl
        << m_strFileName.c_str() << endl
        << "ASCII" << endl << endl
        << "DATASET POLYDATA" << endl;
    ofs << "POINTS " << V.size() << " float" << endl;
    ofs << std::fixed << setprecision(7);
    for (size_t i = 0; i < V.size(); i++)
        ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << endl;
    ofs << "LINES " << E.size() << " " << 3 * E.size() << std::endl;
    for (size_t i = 0; i < E.size(); i++)
        ofs << "2 " << E.at(i).Vids[0] << " " << E.at(i).Vids[1] << std::endl;
    ofs << "CELL_DATA " << E.size() << std::endl
        << "SCALARS " << " Frames" << " int 1" << std::endl
        << "LOOKUP_TABLE default" << std::endl;
    for (size_t i = 0; i < E.size(); i++)
        ofs << i%3 << std::endl;
}

void MeshFileWriter::WriteSharpEdgesVtk()
{
    const std::vector<Vertex>& V = m_mesh.V;
    const std::vector<Edge>& E = m_mesh.E;

    std::ofstream ofs(m_strFileName.c_str());
    ofs << "# vtk DataFile Version 2.0" << endl
        << m_strFileName.c_str() << endl
        << "ASCII" << endl << endl
        << "DATASET POLYDATA" << endl;
    ofs << "POINTS " << V.size() << " float" << endl;
    ofs << std::fixed << setprecision(7);
    for (size_t i = 0; i < V.size(); i++)
        ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << endl;
    size_t numOfSharpEdges = 0;
    for (size_t i = 0; i < E.size(); i++) {
        const Edge& edge = E.at(i);
        if (edge.isSharpFeature)
            ++numOfSharpEdges;
    }

    ofs << "LINES " << numOfSharpEdges << " " << 3 * numOfSharpEdges << std::endl;
    for (size_t i = 0; i < E.size(); i++) {
        const Edge& edge = E.at(i);
        if (edge.isSharpFeature)
            ofs << "2 " << edge.Vids[0] << " " << edge.Vids[1] << std::endl;
    }

    ofs << "CELL_DATA " << numOfSharpEdges << std::endl
        << "SCALARS " << " Feature" << " int 1" << std::endl
        << "LOOKUP_TABLE default" << std::endl;
    for (size_t i = 0; i < E.size(); i++) {
        const Edge& edge = E.at(i);
        if (edge.isSharpFeature)
            ofs << edge.label << std::endl;
    }
}

void MeshFileWriter::WriteCornersVtk()
{
    const std::vector<Vertex>& V = m_mesh.V;
    const std::vector<Edge>& E = m_mesh.E;

    std::ofstream ofs(m_strFileName.c_str());
    ofs << "# vtk DataFile Version 2.0" << endl
        << m_strFileName.c_str() << endl
        << "ASCII" << endl << endl
        << "DATASET POLYDATA" << endl;
    ofs << "POINTS " << V.size() << " float" << endl;
    ofs << std::fixed << setprecision(7);
    for (size_t i = 0; i < V.size(); i++)
        ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << endl;
    size_t numOfCorners = 0;
    for (size_t i = 0; i < V.size(); i++) {
        const Vertex& v = V.at(i);
        if (v.isCorner)
            ++numOfCorners;
    }

    ofs << "VERTICES " << numOfCorners << " " << 2 * numOfCorners << std::endl;
    for (size_t i = 0; i < V.size(); i++) {
        const Vertex& v = V.at(i);
        if (v.isCorner)
            ofs << "1 " << v.id << std::endl;
    }
}

void MeshFileWriter::WriteVertexFeatureVtk() {
    const std::vector<Vertex>& V = m_mesh.V;
    const std::vector<Edge>& E = m_mesh.E;

    std::ofstream ofs(m_strFileName.c_str());
    ofs << "# vtk DataFile Version 2.0" << endl
        << m_strFileName.c_str() << endl
        << "ASCII" << endl << endl
        << "DATASET POLYDATA" << endl;
    ofs << "POINTS " << V.size() << " double" << endl;
    for (auto& v : V)
        ofs << v.x << " " << v.y << " " << v.z << endl;
    ofs << "VERTICES " << V.size() << " " << 2 * V.size() << std::endl;
    for (auto& v : V)
        ofs << "1 " << v.id << std::endl;
	ofs << "CELL_DATA " << V.size() << std::endl
		<< "SCALARS " << " feature" << " int 1" << std::endl
		<< "LOOKUP_TABLE default" << std::endl;
	for (auto& v : V)
		ofs << v.type << std::endl;
	ofs << "SCALARS " << " label" << " int 1" << std::endl
		<< "LOOKUP_TABLE default" << std::endl;
	int numOfSharpEdges = 0;
	for (auto& v : V)
		if (v.label != MAXID && v.label > numOfSharpEdges) numOfSharpEdges = v.label;
	numOfSharpEdges++;
	for (auto& v : V) {
		if (v.label == MAXID) ofs << numOfSharpEdges << std::endl;
		else ofs << v.label << std::endl;
	}
	ofs << "SCALARS " << " patch_id" << " int 1" << std::endl
		<< "LOOKUP_TABLE default" << std::endl;
	int patch_id = 0;
	for (auto& v : V)
		if (v.patch_id != MAXID && v.patch_id > numOfSharpEdges) patch_id = v.label;
	patch_id++;
	for (auto& v : V) {
		if (v.patch_id == MAXID) ofs << patch_id << std::endl;
		else ofs << v.patch_id << std::endl;
	}
}

void MeshFileWriter::WriteEdgesVtk(const std::vector<size_t>& edgeIds) {
    const std::vector<Vertex>& V = m_mesh.V;
    const std::vector<Edge>& E = m_mesh.E;

    std::ofstream ofs(m_strFileName.c_str());
    ofs << "# vtk DataFile Version 2.0" << endl
        << m_strFileName.c_str() << endl
        << "ASCII" << endl << endl
        << "DATASET POLYDATA" << endl;
    ofs << "POINTS " << V.size() << " float" << endl;
    ofs << std::fixed << setprecision(7);
    for (size_t i = 0; i < V.size(); i++)
        ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << endl;
    ofs << "LINES " << edgeIds.size() << " " << 3 * edgeIds.size() << std::endl;
    for (size_t i = 0; i < edgeIds.size(); i++)
        ofs << "2 " << E.at(edgeIds.at(i)).Vids[0] << " " << E.at(edgeIds.at(i)).Vids[1] << std::endl;
    ofs << "CELL_DATA " << edgeIds.size() << std::endl
        << "SCALARS " << " Singularity" << " int 1" << std::endl
        << "LOOKUP_TABLE default" << std::endl;
    for (size_t i = 0; i < edgeIds.size(); i++)
        ofs << (E.at(edgeIds.at(i)).isSingularity ? 1 : 0) << std::endl;

//    ofs << "SCALARS " << " ESingularity" << " float" << std::endl
//        << "LOOKUP_TABLE default" << std::endl;
//    for (size_t i = 0; i < edgeIds.size(); i++)
//        ofs << E.at(edgeIds.at(i)).energySingularity << std::endl;
//
//    ofs << "SCALARS " << " EOrthogonality" << " float" << std::endl
//        << "LOOKUP_TABLE default" << std::endl;
//    for (size_t i = 0; i < edgeIds.size(); i++)
//        ofs << E.at(edgeIds.at(i)).energyOrthogonality << std::endl;
//
//    ofs << "SCALARS " << " EStraightness" << " float" << std::endl
//        << "LOOKUP_TABLE default" << std::endl;
//    for (size_t i = 0; i < edgeIds.size(); i++)
//        ofs << E.at(edgeIds.at(i)).energyStraightness << std::endl;
}


void MeshFileWriter::WriteEdgesVtk(const std::set<size_t>& edgeIds) {
    const std::vector<Vertex>& V = m_mesh.V;
    const std::vector<Edge>& E = m_mesh.E;

    std::ofstream ofs(m_strFileName.c_str());
    ofs << "# vtk DataFile Version 2.0" << endl
        << m_strFileName.c_str() << endl
        << "ASCII" << endl << endl
        << "DATASET POLYDATA" << endl;
    ofs << "POINTS " << V.size() << " float" << endl;
    ofs << std::fixed << setprecision(7);
    for (size_t i = 0; i < V.size(); i++)
        ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << endl;
    ofs << "LINES " << edgeIds.size() << " " << 3 * edgeIds.size() << std::endl;
    for (auto edgeId : edgeIds)
        ofs << "2 " << E.at(edgeId).Vids[0] << " " << E.at(edgeId).Vids[1] << std::endl;
    ofs << "CELL_DATA " << edgeIds.size() << std::endl
        << "SCALARS " << " Singularity" << " int 1" << std::endl
        << "LOOKUP_TABLE default" << std::endl;
    for (auto edgeId : edgeIds)
        ofs << (E.at(edgeId).isSingularity ? 1 : 0) << std::endl;
}

void MeshFileWriter::WriteEdgesVtk(const std::vector<std::set<size_t>>& edgeIds) {
    const std::vector<Vertex>& V = m_mesh.V;
    const std::vector<Edge>& E = m_mesh.E;

    std::ofstream ofs(m_strFileName.c_str());
    ofs << "# vtk DataFile Version 2.0" << endl
        << m_strFileName.c_str() << endl
        << "ASCII" << endl << endl
        << "DATASET POLYDATA" << endl;
    ofs << "POINTS " << V.size() << " float" << endl;
    ofs << std::fixed << setprecision(7);
    for (size_t i = 0; i < V.size(); i++)
        ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << endl;
	auto numberOfLines = 0;
	for (auto& eids : edgeIds) numberOfLines += eids.size();
    ofs << "LINES " << numberOfLines << " " << 3 * numberOfLines << std::endl;
	for (auto& eids : edgeIds)
    for (auto edgeId : eids)
        ofs << "2 " << E.at(edgeId).Vids[0] << " " << E.at(edgeId).Vids[1] << std::endl;
    ofs << "CELL_DATA " << numberOfLines << std::endl
        << "SCALARS " << "id" << " int 1" << std::endl
        << "LOOKUP_TABLE default" << std::endl;
	int id = 0;
	for (auto& eids : edgeIds) {
		for (auto edgeId : eids)
			ofs << id << std::endl;
		++id;
	}
	ofs << "SCALARS " << "size" << " int 1" << std::endl
		<< "LOOKUP_TABLE default" << std::endl;
	for (auto& eids : edgeIds) {
		for (auto edgeId : eids)
			ofs << eids.size() << std::endl;
	}
}

void MeshFileWriter::WriteEdgesVtk(const std::vector<std::vector<size_t>>& edgeIds) {
    const std::vector<Vertex>& V = m_mesh.V;
    const std::vector<Edge>& E = m_mesh.E;

    std::ofstream ofs(m_strFileName.c_str());
    ofs << "# vtk DataFile Version 2.0" << endl
        << m_strFileName.c_str() << endl
        << "ASCII" << endl << endl
        << "DATASET POLYDATA" << endl;
    ofs << "POINTS " << V.size() << " float" << endl;
    ofs << std::fixed << setprecision(7);
    for (size_t i = 0; i < V.size(); i++)
        ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << endl;
	auto numberOfLines = 0;
	for (auto& eids : edgeIds) numberOfLines += eids.size();
    ofs << "LINES " << numberOfLines << " " << 3 * numberOfLines << std::endl;
	for (auto& eids : edgeIds)
    for (auto edgeId : eids)
        ofs << "2 " << E.at(edgeId).Vids[0] << " " << E.at(edgeId).Vids[1] << std::endl;
    ofs << "CELL_DATA " << numberOfLines << std::endl
        << "SCALARS " << "id" << " int 1" << std::endl
        << "LOOKUP_TABLE default" << std::endl;
	int id = 0;
	for (auto& eids : edgeIds) {
		for (auto edgeId : eids)
			ofs << id << std::endl;
		++id;
	}
	ofs << "SCALARS " << "size" << " int 1" << std::endl
		<< "LOOKUP_TABLE default" << std::endl;
	for (auto& eids : edgeIds) {
		for (auto edgeId : eids)
			ofs << eids.size() << std::endl;
	}
}

void MeshFileWriter::WriteFacesVtk()
{
    const std::vector<Vertex>& V = m_mesh.V;
    const std::vector<Face>& F = m_mesh.F;

    std::ofstream ofs(m_strFileName.c_str());
    ofs << "# vtk DataFile Version 2.0" << endl
        << m_strFileName.c_str() << endl
        << "ASCII" << endl << endl
        << "DATASET POLYDATA" << endl;
    ofs << "POINTS " << V.size() << " float" << endl;
    ofs << std::fixed << setprecision(7);
    for (size_t i = 0; i < V.size(); i++)
        ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << endl;
    if (m_mesh.m_cellType == HEXAHEDRA || m_mesh.m_cellType == QUAD) {
        ofs << "POLYGONS " << F.size() << " " << 5 * F.size() << std::endl;
        for (size_t i = 0; i < F.size(); i++) {
            ofs << "4";
            const Face& face = F.at(i);
            for (size_t j = 0; j < face.Vids.size(); j++)
                ofs << " " << face.Vids.at(j);
            ofs << "\n";
        }
    }
    else if (m_mesh.m_cellType == TETRAHEDRA || m_mesh.m_cellType == TRIANGLE) {
        ofs << "POLYGONS " << F.size() << " " << 4 * F.size() << std::endl;
        for (size_t i = 0; i < F.size(); i++) {
            ofs << "3";
            const Face& face = F.at(i);
            for (size_t j = 0; j < face.Vids.size(); j++)
                ofs << " " << face.Vids.at(j);
            ofs << "\n";
        }
    } else if (m_mesh.m_cellType == POLYGON) {
		size_t fvnum = 0;
		for (auto& f : F) fvnum += 1 + f.Vids.size();
		ofs << "POLYGONS " << F.size() << " " << fvnum << std::endl;
		for (auto& f : F) {
			ofs << f.Vids.size();
			for (auto vid : f.Vids)
				ofs << " " << vid;
			ofs << "\n";
		}
	}
    //if (m_mesh.m_cellType == TRIANGLE || m_mesh.m_cellType == QUAD)
    {
//        ofs << "POINT_DATA " << V.size() << std::endl
//            << "SCALARS " << " Type" << " int 1" << std::endl
//            << "LOOKUP_TABLE default" << std::endl;
//        for (size_t i = 0; i < V.size(); i++)
//            ofs << (V.at(i).type == MAXID ? 0 : V.at(i).type) << std::endl;

        ofs << "CELL_DATA " << F.size() << std::endl
            << "SCALARS " << " Label" << " int 1" << std::endl
            << "LOOKUP_TABLE default" << std::endl;
        for (size_t i = 0; i < F.size(); i++)
            ofs << (F.at(i).label == MAXID ? 0 : F.at(i).label)  << std::endl;
    }
}
void MeshFileWriter::WriteFacesVtk(const std::vector<size_t>& faceIds)
{
    const std::vector<Vertex>& V = m_mesh.V;
    const std::vector<Face>& F = m_mesh.F;

    std::ofstream ofs(m_strFileName.c_str());
    ofs << "# vtk DataFile Version 2.0" << endl
        << m_strFileName.c_str() << endl
        << "ASCII" << endl << endl
        << "DATASET POLYDATA" << endl;
    ofs << "POINTS " << V.size() << " float" << endl;
    ofs << std::fixed << setprecision(7);
    for (size_t i = 0; i < V.size(); i++)
        ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << endl;
    if (m_mesh.m_cellType == HEXAHEDRA || m_mesh.m_cellType == QUAD) {
        ofs << "POLYGONS " << faceIds.size() << " " << 5 * faceIds.size() << std::endl;
        for (size_t i = 0; i < faceIds.size(); i++) {
            ofs << "4";
            const Face& face = F.at(faceIds.at(i));
            for (size_t j = 0; j < face.Vids.size(); j++)
                ofs << " " << face.Vids.at(j);
            ofs << "\n";
        }
    } else if (m_mesh.m_cellType == TETRAHEDRA || m_mesh.m_cellType == TRIANGLE) {
        ofs << "POLYGONS " << faceIds.size() << " " << 4 * faceIds.size() << std::endl;
        for (size_t i = 0; i < faceIds.size(); i++) {
            ofs << "3";
            const Face& face = F.at(faceIds.at(i));
            for (size_t j = 0; j < face.Vids.size(); j++)
                ofs << " " << face.Vids.at(j);
            ofs << "\n";
        }
    } else if (m_mesh.m_cellType == POLYGON) {
		size_t fvnum = 0; 
		for (auto& fid : faceIds) fvnum += 1 + F[fid].Vids.size();
		ofs << "POLYGONS " << faceIds.size() << " " << fvnum << std::endl;
		for (auto& fid : faceIds) {
			const Face& face = F[fid].Vids;
			ofs << face.Vids.size();
			for (auto vid : face.Vids)
				ofs << " " << vid;
			ofs << "\n";
		}
	}
    //if (m_mesh.m_cellType == TRIANGLE || m_mesh.m_cellType == QUAD)
    {
//        ofs << "POINT_DATA " << V.size() << std::endl
//            << "SCALARS " << " Type" << " int 1" << std::endl
//            << "LOOKUP_TABLE default" << std::endl;
//        for (size_t i = 0; i < V.size(); i++)
//            ofs << (V.at(i).type == MAXID ? 0 : V.at(i).type) << std::endl;

        ofs << "CELL_DATA " << faceIds.size() << std::endl
            << "SCALARS " << " Label" << " int 1" << std::endl
            << "LOOKUP_TABLE default" << std::endl;
        for (size_t i = 0; i < faceIds.size(); i++)
            ofs << (F.at(faceIds.at(i)).label == MAXID ? 0 : F.at(faceIds.at(i)).label)  << std::endl;
        ofs << "SCALARS " << " Color" << " int 1" << std::endl
            << "LOOKUP_TABLE default" << std::endl;
        for (size_t i = 0; i < faceIds.size(); i++)
            ofs << (F.at(faceIds.at(i)).label == MAXID ? 0 : F.at(faceIds.at(i)).label % 8)  << std::endl;
    }
}

void MeshFileWriter::WriteFacesVtk(const std::set<size_t>& faceIds) {
    const std::vector<Vertex>& V = m_mesh.V;
    const std::vector<Face>& F = m_mesh.F;

    std::ofstream ofs(m_strFileName.c_str());
    ofs << "# vtk DataFile Version 2.0" << endl
        << m_strFileName.c_str() << endl
        << "ASCII" << endl << endl
        << "DATASET POLYDATA" << endl;
    ofs << "POINTS " << V.size() << " float" << endl;
    ofs << std::fixed << setprecision(7);
    for (size_t i = 0; i < V.size(); i++)
        ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << endl;
    if (m_mesh.m_cellType == HEXAHEDRA || m_mesh.m_cellType == QUAD) {
        ofs << "POLYGONS " << faceIds.size() << " " << 5 * faceIds.size() << std::endl;
        for (auto faceId : faceIds) {
            ofs << "4";
            const Face& face = F.at(faceId);
            for (size_t j = 0; j < face.Vids.size(); j++)
                ofs << " " << face.Vids.at(j);
            ofs << "\n";
        }
    } else if (m_mesh.m_cellType == TETRAHEDRA || m_mesh.m_cellType == TRIANGLE) {
        ofs << "POLYGONS " << faceIds.size() << " " << 4 * faceIds.size() << std::endl;
        for (auto faceId : faceIds) {
            ofs << "3";
            const Face& face = F.at(faceId);
            for (size_t j = 0; j < face.Vids.size(); j++)
                ofs << " " << face.Vids.at(j);
            ofs << "\n";
        }
    }
    //if (m_mesh.m_cellType == TRIANGLE || m_mesh.m_cellType == QUAD)
    {
        //        ofs << "POINT_DATA " << V.size() << std::endl
        //            << "SCALARS " << " Type" << " int 1" << std::endl
        //            << "LOOKUP_TABLE default" << std::endl;
        //        for (size_t i = 0; i < V.size(); i++)
        //            ofs << (V.at(i).type == MAXID ? 0 : V.at(i).type) << std::endl;

        ofs << "CELL_DATA " << faceIds.size() << std::endl
            << "SCALARS " << " Label" << " int 1" << std::endl
            << "LOOKUP_TABLE default" << std::endl;
        for (auto faceId : faceIds)
            ofs << (F.at(faceId).label == MAXID ? 0 : F.at(faceId).label) << std::endl;
        ofs << "SCALARS " << " Color" << " int 1" << std::endl
            << "LOOKUP_TABLE default" << std::endl;
        for (auto faceId : faceIds)
            ofs << (F.at(faceId).label == MAXID ? 0 : F.at(faceId).label % 8) << std::endl;
    }
}

void MeshFileWriter::WriteFacesVtk(const std::vector<std::set<size_t>>& faceIds) {
    const std::vector<Vertex>& V = m_mesh.V;
    const std::vector<Face>& F = m_mesh.F;

    std::ofstream ofs(m_strFileName.c_str());
    ofs << "# vtk DataFile Version 2.0" << endl
        << m_strFileName.c_str() << endl
        << "ASCII" << endl << endl
        << "DATASET POLYDATA" << endl;
    ofs << "POINTS " << V.size() << " float" << endl;
    ofs << std::fixed << setprecision(7);
    for (size_t i = 0; i < V.size(); i++)
        ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << endl;
    size_t numOfFaces = 0;
    for (auto& s : faceIds) numOfFaces += s.size();
    if (m_mesh.m_cellType == HEXAHEDRA || m_mesh.m_cellType == QUAD) {
        ofs << "POLYGONS " << numOfFaces << " " << 5 * numOfFaces << std::endl;
        for (auto& s : faceIds)
        for (auto faceId : s) {
            ofs << "4";
            const Face& face = F.at(faceId);
            for (size_t j = 0; j < face.Vids.size(); j++)
                ofs << " " << face.Vids.at(j);
            ofs << "\n";
        }
    } else if (m_mesh.m_cellType == TETRAHEDRA || m_mesh.m_cellType == TRIANGLE) {
        ofs << "POLYGONS " << numOfFaces << " " << 4 * numOfFaces << std::endl;
        for (auto& s : faceIds)
        for (auto faceId : s) {
            ofs << "3";
            const Face& face = F.at(faceId);
            for (size_t j = 0; j < face.Vids.size(); j++)
                ofs << " " << face.Vids.at(j);
            ofs << "\n";
        }
    }
    //if (m_mesh.m_cellType == TRIANGLE || m_mesh.m_cellType == QUAD)
    {
        //        ofs << "POINT_DATA " << V.size() << std::endl
        //            << "SCALARS " << " Type" << " int 1" << std::endl
        //            << "LOOKUP_TABLE default" << std::endl;
        //        for (size_t i = 0; i < V.size(); i++)
        //            ofs << (V.at(i).type == MAXID ? 0 : V.at(i).type) << std::endl;

        ofs << "CELL_DATA " << numOfFaces << std::endl
            << "SCALARS " << " Id" << " int 1" << std::endl
            << "LOOKUP_TABLE default" << std::endl;
        int id = 0;
        for (auto& s : faceIds) {
            for (auto faceId : s)
                ofs << id << std::endl;
            ++id;
        }
        //ofs << "SCALARS " << " Label" << " int 1" << std::endl
        //    << "LOOKUP_TABLE default" << std::endl;
        //for (auto& s : faceIds)
        //for (auto faceId : s)
        //    ofs << (F.at(faceId).label == MAXID ? 0 : F.at(faceId).label) << std::endl;
    }
}

void MeshFileWriter::WriteFacesVtk(const std::vector<std::vector<size_t>>& faceIds) {
	const std::vector<Vertex>& V = m_mesh.V;
	const std::vector<Face>& F = m_mesh.F;

	std::ofstream ofs(m_strFileName.c_str());
	ofs << "# vtk DataFile Version 2.0" << endl
		<< m_strFileName.c_str() << endl
		<< "ASCII" << endl << endl
		<< "DATASET POLYDATA" << endl;
	ofs << "POINTS " << V.size() << " float" << endl;
	ofs << std::fixed << setprecision(7);
	for (size_t i = 0; i < V.size(); i++)
		ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << endl;
	size_t numOfFaces = 0;
	for (auto& s : faceIds) numOfFaces += s.size();
	if (m_mesh.m_cellType == HEXAHEDRA || m_mesh.m_cellType == QUAD) {
		ofs << "POLYGONS " << numOfFaces << " " << 5 * numOfFaces << std::endl;
		for (auto& s : faceIds)
			for (auto faceId : s) {
				ofs << "4";
				const Face& face = F.at(faceId);
				for (size_t j = 0; j < face.Vids.size(); j++)
					ofs << " " << face.Vids.at(j);
				ofs << "\n";
			}
	} else if (m_mesh.m_cellType == TETRAHEDRA || m_mesh.m_cellType == TRIANGLE) {
		ofs << "POLYGONS " << numOfFaces << " " << 4 * numOfFaces << std::endl;
		for (auto& s : faceIds)
			for (auto faceId : s) {
				ofs << "3";
				const Face& face = F.at(faceId);
				for (size_t j = 0; j < face.Vids.size(); j++)
					ofs << " " << face.Vids.at(j);
				ofs << "\n";
			}
	} else if (m_mesh.m_cellType == POLYGON) {
		size_t numOfFaceVids = 0;
		for (auto& s : faceIds)
			for (auto fid : s) numOfFaceVids += 1 + F.at(fid).Vids.size();
		ofs << "POLYGONS " << numOfFaces << " " << numOfFaceVids << std::endl;
		for (auto& s : faceIds)
			for (auto faceId : s) {
				ofs << F.at(faceId).Vids.size();
				const Face& face = F.at(faceId);
				for (size_t j = 0; j < face.Vids.size(); j++)
					ofs << " " << face.Vids.at(j);
				ofs << "\n";
			}
	}
	
	ofs << "CELL_DATA " << numOfFaces << std::endl
		<< "SCALARS " << " Id" << " int 1" << std::endl
		<< "LOOKUP_TABLE default" << std::endl;
	int id = 0;
	for (auto& s : faceIds) {
		for (auto faceId : s)
			ofs << id << std::endl;
		++id;
	}
}

void MeshFileWriter::WriteFaceSegmentsVtk()
{
    const std::vector<Vertex>& V = m_mesh.V;
    const std::vector<Face>& F = m_mesh.F;

    std::ofstream ofs(m_strFileName.c_str());
    ofs << "# vtk DataFile Version 2.0" << endl
        << m_strFileName.c_str() << endl
        << "ASCII" << endl << endl
        << "DATASET POLYDATA" << endl;
    ofs << "POINTS " << V.size() << " float" << endl;
    ofs << std::fixed << setprecision(7);
    for (size_t i = 0; i < V.size(); i++)
        ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << endl;

    size_t count = 0;
    if (m_mesh.m_cellType == HEXAHEDRA || m_mesh.m_cellType == QUAD) {
        ofs << "POLYGONS " << F.size() << " " << 5 * F.size() << std::endl;
        for (size_t i = 0; i < F.size(); i++) {
            ofs << "4";
            const Face& face = F.at(i);
            for (size_t j = 0; j < face.Vids.size(); j++)
                ofs << " " << face.Vids.at(j);
            ofs << "\n";
        }
    }
    else if (m_mesh.m_cellType == TETRAHEDRA || m_mesh.m_cellType == TRIANGLE) {
        for (size_t i = 0; i < F.size(); i++) {
            const Face& face = F.at(i);
            int scalar0 = m_mesh.pointScalarFields[0][face.Vids[0]];
            int scalar1 = m_mesh.pointScalarFields[0][face.Vids[1]];
            int scalar2 = m_mesh.pointScalarFields[0][face.Vids[2]];
            if (scalar0 == scalar1 && scalar1 == scalar2)
                continue;
            count++;
        }

        ofs << "POLYGONS " << count << " " << 4 * count << std::endl;
        for (size_t i = 0; i < F.size(); i++) {
            const Face& face = F.at(i);
            int scalar0 = m_mesh.pointScalarFields[0][face.Vids[0]];
            int scalar1 = m_mesh.pointScalarFields[0][face.Vids[1]];
            int scalar2 = m_mesh.pointScalarFields[0][face.Vids[2]];
            if (scalar0 == scalar1 && scalar1 == scalar2) continue;
            ofs << "3";
            for (size_t j = 0; j < face.Vids.size(); j++)
                ofs << " " << face.Vids.at(j);
            ofs << "\n";
        }
    }
    //if (m_mesh.m_cellType == TRIANGLE || m_mesh.m_cellType == QUAD)
//    {
//        ofs << "POINT_DATA " << V.size() << std::endl
//            << "SCALARS " << " Type" << " int 1" << std::endl
//            << "LOOKUP_TABLE default" << std::endl;
//        for (size_t i = 0; i < V.size(); i++)
//            ofs << (V.at(i).type == MAXID ? 0 : V.at(i).type) << std::endl;
//
//        ofs << "CELL_DATA " << F.size() << std::endl
//            << "SCALARS " << " Label" << " int 1" << std::endl
//            << "LOOKUP_TABLE default" << std::endl;
//        for (size_t i = 0; i < F.size(); i++)
//            ofs << (F.at(i).label == MAXID ? 0 : F.at(i).label)  << std::endl;
//    }
}

void MeshFileWriter::WriteCellsVtk(const std::vector<size_t>& cellIds) {
	const std::vector<Vertex>& V = m_mesh.V;
	const std::vector<Cell>& C = m_mesh.C;

	std::ofstream ofs(m_strFileName.c_str());
	ofs << "# vtk DataFile Version 4.0" << endl
		<< m_strFileName.c_str() << endl
		<< "ASCII" << endl << endl
		<< "DATASET UNSTRUCTURED_GRID" << endl;
	ofs << "POINTS " << V.size() << " double" << endl;
	ofs << std::fixed << setprecision(7);
	for (size_t i = 0; i < V.size(); i++)
		ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << endl;
	size_t numOfCellVids = 0;
	for (auto cellId : cellIds) {
		const auto& cell = C.at(cellId);
		numOfCellVids += cell.Vids.size() + 1;
	}
	ofs << "CELLS " << cellIds.size() << " " << numOfCellVids << std::endl;
	for (auto cellId : cellIds) {
		const auto& cell = C.at(cellId);
		ofs << cell.Vids.size();
		for (auto vid : cell.Vids)
			ofs << " " << vid;
		ofs << "\n";
	}
	ofs << "CELL_TYPES " << cellIds.size() << std::endl;
	for (auto cellId : cellIds) {
		const auto& cell = C.at(cellId);
		ofs << cell.cellType << std::endl;
	}
}

void MeshFileWriter::WriteCellsVtk(const std::set<size_t>& cellIds) {
	const std::vector<Vertex>& V = m_mesh.V;
	const std::vector<Cell>& C = m_mesh.C;

	std::ofstream ofs(m_strFileName.c_str());
	ofs << "# vtk DataFile Version 4.0" << endl
		<< m_strFileName.c_str() << endl
		<< "ASCII" << endl << endl
		<< "DATASET UNSTRUCTURED_GRID" << endl;
	ofs << "POINTS " << V.size() << " double" << endl;
	ofs << std::fixed << setprecision(7);
	for (size_t i = 0; i < V.size(); i++)
		ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << endl;
	size_t numOfCellVids = 0;
	for (auto cellId : cellIds) {
		const auto& cell = C.at(cellId);
		numOfCellVids += cell.Vids.size() + 1;
	}
	ofs << "CELLS " << cellIds.size() << " " << numOfCellVids << std::endl;
	for (auto cellId : cellIds) {
		const auto& cell = C.at(cellId);
		ofs << cell.Vids.size();
		for (auto vid : cell.Vids)
			ofs << " " << vid;
		ofs << "\n";
	}	
	ofs << "CELL_TYPES " << cellIds.size() << std::endl;
	for (auto cellId : cellIds) {
		const auto& cell = C.at(cellId);
		ofs << cell.cellType << std::endl;
	}
}

void MeshFileWriter::WriteVerticesVtk()
{
    const std::vector<Vertex>& V = m_mesh.V;

    std::ofstream ofs(m_strFileName.c_str());
    ofs << "# vtk DataFile Version 2.0" << endl
        << m_strFileName.c_str() << endl
        << "ASCII" << endl << endl
        << "DATASET POLYDATA" << endl;
    ofs << "POINTS " << V.size() << " float" << endl;
    ofs << std::fixed << setprecision(7);
    for (size_t i = 0; i < V.size(); i++)
        ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << endl;
    ofs << "VERTICES " << V.size() << " " << 2 * V.size() << std::endl;
    for (size_t i = 0; i < V.size(); i++)
        ofs << "1 " << V.at(i).id << std::endl;
}

void MeshFileWriter::WriteVerticesVtk(const std::vector<size_t>& vertexIds)
{
    const std::vector<Vertex>& V = m_mesh.V;

    std::ofstream ofs(m_strFileName.c_str());
    ofs << "# vtk DataFile Version 2.0" << endl
        << m_strFileName.c_str() << endl
        << "ASCII" << endl << endl
        << "DATASET POLYDATA" << endl;
    ofs << "POINTS " << V.size() << " float" << endl;
    ofs << std::fixed << setprecision(7);
    for (size_t i = 0; i < V.size(); i++)
        ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << endl;
    ofs << "VERTICES " << vertexIds.size() << " " << 2 * vertexIds.size() << std::endl;
    for (size_t i = 0; i < vertexIds.size(); i++)
        ofs << "1 " << V.at(vertexIds.at(i)).id << std::endl;
    ofs << "CELL_DATA " << vertexIds.size() << std::endl
        << "SCALARS " << " Singularity" << " int 1" << std::endl
        << "LOOKUP_TABLE default" << std::endl;
    for (size_t i = 0; i < vertexIds.size(); i++)
        ofs << (V.at(vertexIds.at(i)).isSingularity ? 1 : 0) << std::endl;
    ofs <<  "SCALARS " << " Valence" << " int 1" << std::endl
        << "LOOKUP_TABLE default" << std::endl;
    for (size_t i = 0; i < vertexIds.size(); i++)
        ofs << V.at(vertexIds.at(i)).N_Fids.size() << std::endl;
}
void MeshFileWriter::WriteVerticesVtk(const std::set<size_t>& vertexIds) {
    const std::vector<Vertex>& V = m_mesh.V;

    std::ofstream ofs(m_strFileName.c_str());
    ofs << "# vtk DataFile Version 2.0" << endl
        << m_strFileName.c_str() << endl
        << "ASCII" << endl << endl
        << "DATASET POLYDATA" << endl;
    ofs << "POINTS " << V.size() << " float" << endl;
    ofs << std::fixed << setprecision(7);
    for (size_t i = 0; i < V.size(); i++)
        ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << endl;
    ofs << "VERTICES " << vertexIds.size() << " " << 2 * vertexIds.size() << std::endl;
    for (auto vertexId : vertexIds)
        ofs << "1 " << vertexId << std::endl;
    ofs << "CELL_DATA " << vertexIds.size() << std::endl
        << "SCALARS " << " Singularity" << " int 1" << std::endl
        << "LOOKUP_TABLE default" << std::endl;
    for (auto vertexId : vertexIds)
        ofs << (V.at(vertexId).isSingularity ? 1 : 0) << std::endl;
}

void MeshFileWriter::WriteVerticesVtk(const std::vector<std::vector<size_t>>& vertexIds)
{
    const std::vector<Vertex>& V = m_mesh.V;

    std::ofstream ofs(m_strFileName.c_str());
    ofs << "# vtk DataFile Version 2.0" << endl
        << m_strFileName.c_str() << endl
        << "ASCII" << endl << endl
        << "DATASET POLYDATA" << endl;
    ofs << "POINTS " << V.size() << " float" << endl;
    ofs << std::fixed << setprecision(7);
    for (size_t i = 0; i < V.size(); i++)
        ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << endl;
	size_t total_num_of_vids = 0;
	for (auto& vids : vertexIds) total_num_of_vids += vids.size();
    ofs << "VERTICES " << total_num_of_vids << " " << 2 * total_num_of_vids << std::endl;
	for (auto& vids : vertexIds)
		for (auto vid : vids)
			ofs << "1 " << vid << std::endl;
    ofs << "CELL_DATA " << total_num_of_vids << std::endl
        << "SCALARS " << " id" << " int 1" << std::endl
        << "LOOKUP_TABLE default" << std::endl;
	int id = 0;
	for (auto& vids : vertexIds) {
		for (auto vid : vids)
			ofs << id << std::endl;
		++id;
	}
}
void MeshFileWriter::WriteLinksVtk(const std::vector<std::vector<size_t>>& linkVids) {
    const std::vector<Vertex>& V = m_mesh.V;
    std::ofstream ofs(m_strFileName.c_str());
    ofs << "# vtk DataFile Version 2.0" << endl
        << m_strFileName.c_str() << endl
        << "ASCII" << endl << endl
        << "DATASET POLYDATA" << endl;
    ofs << "POINTS " << V.size() << " float" << endl;
    ofs << std::fixed << setprecision(7);
    for (size_t i = 0; i < V.size(); i++)
        ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << endl;
	size_t numOfLinkVids = 0;
	for (auto& link : linkVids) numOfLinkVids += 1 + link.size();
    ofs << "LINES " << linkVids.size() << " " << numOfLinkVids << std::endl;
	for (auto& link : linkVids) {
		ofs << link.size();
		for (auto vid : link) ofs << " " << vid;
		ofs << std::endl;
	}
    ofs << "CELL_DATA " << linkVids.size() << std::endl
        << "SCALARS " << " id" << " int 1" << std::endl
        << "LOOKUP_TABLE default" << std::endl;
    for (size_t i = 0; i < linkVids.size(); i++)
        ofs << i << std::endl;
}
void MeshFileWriter::WriteSurfaceOff()
{
    const std::vector<Vertex>& V = m_mesh.V;
    const std::vector<Face>& F = m_mesh.F;

    size_t surfaceNum = 0;
    for (size_t i = 0; i < F.size(); i++)
        if (F.at(i).isBoundary)
            surfaceNum++;

    std::ofstream ofs(m_strFileName.c_str());
    ofs << "OFF\n";
    ofs << V.size() << " " << surfaceNum << " 0\n";
    ofs << std::fixed << setprecision(7);
    for (size_t i = 0; i < V.size(); i++)
        ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << endl;
    for (size_t i = 0; i < F.size(); i++){
        const Face& face = F.at(i);
        if (face.isBoundary){
            ofs << face.Vids.size();
            for (size_t j = 0; j < face.Vids.size(); j++)
                ofs << " " << face.Vids.at(j);
            ofs << "\n";
        }
    }
}

static void GetQuadVertexIds(size_t i, size_t j, int n, size_t quadVids[]){
    quadVids[0] = i * (n + 1) + j;
    quadVids[1] = i * (n + 1) + j + 1;
    quadVids[2] = (i + 1) * (n + 1) + j + 1;
    quadVids[3] = (i + 1) * (n + 1) + j;
}

void MeshFileWriter::WriteMatrixVTK(const std::vector<std::vector<size_t>>& m) const {
    const size_t n = m.size();
    const size_t numOfVertices = (n + 1) * (n + 1);
    const size_t numOfCells = n * n;
    std::vector<glm::vec2> V(numOfVertices);
    for (size_t i = 0; i <= n; ++i)
        for (size_t j = 0; j <= n; ++j)
            V[i * (n + 1) + j] = glm::vec2(i, j);

    std::ofstream ofs(m_strFileName);
    ofs << "# vtk DataFile Version 3.0\n"
        << m_strFileName << "\n"
        << "ASCII\n\n"
        << "DATASET POLYDATA\n";

    ofs << "POINTS " << numOfVertices << " double\n";
    for (size_t i = 0; i < V.size(); i++)
        ofs << V.at(i).x << " " << V.at(i).y << " 0\n";

    ofs << "POLYGONS " << numOfCells << " " << 5 * numOfCells << "\n";
    size_t quadVids[4] = {0};
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j) {
            GetQuadVertexIds(i, j, n, quadVids);
            ofs << "4 " << quadVids[0] << " " << quadVids[1] << " " << quadVids[2] << " " << quadVids[3] << "\n";
        }

    ofs << "CELL_DATA " << numOfCells << "\n"
        << "SCALARS " << "label" << " int 1\n"
        << "LOOKUP_TABLE default\n";
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            ofs << m[i][j] << "\n";
}

void MeshFileWriter::WriteFeatureLinesVtk(const std::vector<FeatureLine>& featureLines) {
    const std::vector<Vertex>& V = m_mesh.V;
    const std::vector<Edge>& E = m_mesh.E;

    std::ofstream ofs(m_strFileName);
    ofs << "# vtk DataFile Version 2.0" << std::endl
        << m_strFileName << std::endl
        << "ASCII" << std::endl << std::endl
        << "DATASET UNSTRUCTURED_GRID" << std::endl;
    ofs << "POINTS " << V.size() << " double" << std::endl;
    for (auto& v : V)
        ofs << v.x << " " << v.y << " " << v.z << std::endl;
    size_t numOfSharpVertices = 0;
    for (const FeatureLine& fl : featureLines)
        numOfSharpVertices += fl.Vids.size();

    ofs << "CELLS " << featureLines.size() << " " << numOfSharpVertices + featureLines.size() << std::endl;
    for (const auto& fl : featureLines) {
        ofs << fl.Vids.size();
        for (const auto vid : fl.Vids)
            ofs << " " << vid;
        ofs << std::endl;
    }

    ofs << "CELL_TYPES " << featureLines.size() << std::endl;
    for (const FeatureLine& fl : featureLines)
        ofs << 4 << std::endl;

    ofs << "CELL_DATA " << featureLines.size() << std::endl
        << "SCALARS " << " Feature" << " int 1" << std::endl
        << "LOOKUP_TABLE default" << std::endl;
    for (size_t i = 0; i < featureLines.size(); i++) {
        ofs << i << std::endl;
    }
}

void MeshFileWriter::WriteCellsWithFeatureLinesVtk(const std::vector<FeatureLine>& featureLines) {
    const std::vector<Vertex>& V = m_mesh.V;
    const std::vector<Edge>& E = m_mesh.E;
    const std::vector<Cell>& C = m_mesh.C;

    vtkIdType idType = VTK_TRIANGLE;
    if (m_mesh.m_cellType == QUAD) idType = VTK_QUAD;
    else if (m_mesh.m_cellType == TETRAHEDRA) idType = VTK_TETRA;
    else if (m_mesh.m_cellType == HEXAHEDRA) idType = VTK_HEXAHEDRON;

    std::ofstream ofs(m_strFileName);
    ofs << "# vtk DataFile Version 2.0" << std::endl
        << m_strFileName << std::endl
        << "ASCII" << std::endl << std::endl
        << "DATASET UNSTRUCTURED_GRID" << std::endl;
    ofs << "POINTS " << V.size() << " double" << std::endl;
    for (auto& v : V)
        ofs << v.x << " " << v.y << " " << v.z << std::endl;
    size_t numOfSharpVertices = 0;
    for (const FeatureLine& fl : featureLines)
        numOfSharpVertices += fl.Vids.size();

    size_t numOfCorners = 0;
    for (const Vertex& v : V)
        if (v.isCorner) ++numOfCorners;

    ofs << "CELLS " << C.size() + featureLines.size() + numOfCorners << " " << C.size() * (C.front().Vids.size() + 1) +
            numOfSharpVertices + featureLines.size() + numOfCorners * 2 << std::endl;
    for (const auto& c : C) {
        ofs << c.Vids.size();
        for (const auto vid : c.Vids)
            ofs << " " << vid;
        ofs << std::endl;
    }
    for (const auto& fl : featureLines) {
        ofs << fl.Vids.size();
        for (const auto vid : fl.Vids)
            ofs << " " << vid;
        ofs << std::endl;
    }
    for (const Vertex& v : V)
        if (v.isCorner) ofs << "1 " << v.id << std::endl;

    ofs << "CELL_TYPES " << C.size() + featureLines.size() + numOfCorners << std::endl;
    if (m_mesh.m_cellType != POLYHEDRA)
        for (size_t i = 0; i < C.size(); i++)
            ofs << idType << "\n";
    else
        for (auto cellType : m_mesh.m_cellTypes)
            ofs << cellType << "\n";
    for (const FeatureLine& fl : featureLines)
        ofs << VTK_POLY_LINE << std::endl;
    for (size_t i = 0; i < numOfCorners; i++)
        ofs << VTK_VERTEX << std::endl;

    ofs << "CELL_DATA " << C.size() + featureLines.size() + numOfCorners << std::endl
        << "SCALARS " << " Feature" << " int 1" << std::endl
        << "LOOKUP_TABLE default" << std::endl;
    for (size_t i = 0; i < C.size(); i++)
        ofs << 0 << "\n";
    for (size_t i = 0; i < featureLines.size(); i++)
        ofs << i + 1 << std::endl;
    for (size_t i = 0; i < numOfCorners; i++)
        ofs << featureLines.size() + 1 << std::endl;
}

void MeshFileWriter::WriteCellsWithFeatureVtk(const std::vector<FeatureLine>& featureLines, const std::vector<size_t>& featureVids) {
    const std::vector<Vertex>& V = m_mesh.V;
    const std::vector<Edge>& E = m_mesh.E;
    const std::vector<Cell>& C = m_mesh.C;

    vtkIdType idType = VTK_TRIANGLE;
    if (m_mesh.m_cellType == QUAD) idType = VTK_QUAD;
    else if (m_mesh.m_cellType == TETRAHEDRA) idType = VTK_TETRA;
    else if (m_mesh.m_cellType == HEXAHEDRA) idType = VTK_HEXAHEDRON;

    std::ofstream ofs(m_strFileName);
    ofs << "# vtk DataFile Version 2.0" << std::endl
        << m_strFileName << std::endl
        << "ASCII" << std::endl << std::endl
        << "DATASET UNSTRUCTURED_GRID" << std::endl;
    ofs << "POINTS " << V.size() << " double" << std::endl;
    for (auto& v : V)
        ofs << v.x << " " << v.y << " " << v.z << std::endl;
    size_t numOfSharpVertices = 0;
    for (const FeatureLine& fl : featureLines)
        numOfSharpVertices += fl.Vids.size();

    ofs << "CELLS " << C.size() + featureLines.size() + featureVids.size() << " " << C.size() * (C.front().Vids.size() + 1) +
            numOfSharpVertices + featureLines.size() + featureVids.size() * 2 << std::endl;
    for (const auto& c : C) {
        ofs << c.Vids.size();
        for (const auto vid : c.Vids)
            ofs << " " << vid;
        ofs << std::endl;
    }
    for (const auto& fl : featureLines) {
        ofs << fl.Vids.size();
        for (const auto vid : fl.Vids)
            ofs << " " << vid;
        ofs << std::endl;
    }
    for (const auto& vid : featureVids)
        ofs << "1 " << vid << std::endl;

    ofs << "CELL_TYPES " << C.size() + featureLines.size() + featureVids.size() << std::endl;
    if (m_mesh.m_cellType != POLYHEDRA)
        for (size_t i = 0; i < C.size(); i++)
            ofs << idType << "\n";
    else
        for (auto cellType : m_mesh.m_cellTypes)
            ofs << cellType << "\n";
    for (const FeatureLine& fl : featureLines)
        ofs << VTK_POLY_LINE << std::endl;
    for (size_t i = 0; i < featureVids.size(); i++)
        ofs << VTK_VERTEX << std::endl;

    ofs << "CELL_DATA " << C.size() + featureLines.size() + featureVids.size() << std::endl
        << "SCALARS " << " Feature" << " int 1" << std::endl
        << "LOOKUP_TABLE default" << std::endl;
    for (size_t i = 0; i < C.size(); i++)
        ofs << 0 << "\n";
    for (size_t i = 0; i < featureLines.size(); i++)
        ofs << i + 1 << std::endl;
    for (size_t i = 0; i < featureVids.size(); i++)
        ofs << featureLines.size() + 1 << std::endl;
}
