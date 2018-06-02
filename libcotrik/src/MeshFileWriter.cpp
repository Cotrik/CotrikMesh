#include "MeshFileWriter.h"
#include <map>
#include <iostream>
#include <iomanip>
#include <memory>
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

void MeshFileWriter::WriteVtkFile()
{
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
    ofs << "POINTS " << vnum << " float" << endl;
    for (size_t i = 0; i < vnum; i++)
        ofs << std::fixed << setprecision(7) << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << endl;
    ofs << "POLYGONS " << cnum << " ";

    vtkIdType idType = VTK_TRIANGLE;
    if (m_mesh.m_cellType == TRIANGLE) ofs << 4*cnum << std::endl;
    else if (m_mesh.m_cellType == QUAD) {idType = VTK_QUAD;  ofs << 5*cnum << std::endl;}
    else if (m_mesh.m_cellType == TETRAHEDRA) {idType = VTK_TETRA; ofs << 5*cnum << std::endl;}
    else if (m_mesh.m_cellType == HEXAHEDRA) {idType = VTK_HEXAHEDRON; ofs << 9*cnum << std::endl;}

    for (size_t i = 0; i < cnum; i++){
        ofs << C.at(i).Vids.size();
        for (size_t j = 0; j < C.at(i).Vids.size(); j++)
            ofs << " " << C.at(i).Vids.at(j);
        ofs << std::endl;
    }
//    ofs << "CELL_TYPES " << cnum << endl;
//    for (size_t i = 0; i < cnum; i++)
//        ofs << idType << std::endl;
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

void MeshFileWriter::WriteCellData(const std::vector<int>& cellData, const char* dataName)
{
    std::ofstream ofs(m_strFileName.c_str(), std::ios_base::app);
    ofs << "CELL_DATA " << cellData.size() << std::endl
            << "SCALARS " << dataName << " int 1" << std::endl
            << "LOOKUP_TABLE default" << std::endl;
    for (size_t i = 0; i < cellData.size(); i++)
    {
        ofs << cellData[i] << std::endl;
    }
    ofs.close();
}

void MeshFileWriter::WritePointData(const std::vector<int>& pointData, const char* dataName)
{
    std::ofstream ofs(m_strFileName.c_str(), std::ios_base::app);
    ofs << "POINT_DATA " << pointData.size() << std::endl
            << "SCALARS " << dataName << " int 1" << std::endl
            << "LOOKUP_TABLE default" << std::endl;
    for (size_t i = 0; i < pointData.size(); i++)
    {
        ofs << pointData[i] << std::endl;
    }
    ofs.close();
}

void MeshFileWriter::WritePointData(const std::vector<float>& pointData, const char* dataName)
{
    std::ofstream ofs(m_strFileName.c_str(), std::ios_base::app);
    ofs << "POINT_DATA " << pointData.size() << std::endl
            << "SCALARS " << dataName << " float 1" << std::endl
            << "LOOKUP_TABLE default" << std::endl;
    for (size_t i = 0; i < pointData.size(); i++)
    {
        ofs << pointData[i] << std::endl;
    }
    ofs.close();
}

void MeshFileWriter::SetFixFlag(bool bFixed)
{
	m_bFixed = bFixed;
}
void MeshFileWriter::FixMesh()
{
	std::vector<unsigned long> v_real_index;
	size_t c_size = 0;
	for (size_t i = 0; i < m_mesh.C.size(); i++)
	{
		if (m_mesh.C.at(i).Vids.size() > 8)
			continue;
		else
		{
			for (size_t j = 0; j < m_mesh.C.at(i).Vids.size(); j++)
			{
				v_real_index.push_back(m_mesh.C.at(i).Vids.at(j));
			}
			c_size++;
		}
	}

	std::sort(v_real_index.begin(), v_real_index.end());
	std::vector<unsigned long>::iterator iter = std::unique(v_real_index.begin(), v_real_index.end());
	v_real_index.resize(std::distance(v_real_index.begin(), iter));

	std::vector<Vertex> V(v_real_index.size());
	for (unsigned long i = 0; i < v_real_index.size(); i++)
	{
		const Vertex& v = m_mesh.V.at(v_real_index.at(i));
		V.at(i) = v;
	}
	m_mesh.V = V;
	//////////////////////////////////////////////////////
	std::map<size_t, size_t> v_v;
	size_t index = 0;
	for (size_t i = 0; i < v_real_index.size(); i++)
	{
		v_v[v_real_index.at(i)] = i;
	}

	std::vector<Cell> C;
	for (size_t i = 0; i < m_mesh.C.size(); i++)
	{
		Cell c;
		const size_t cellSize = m_mesh.C.at(i).Vids.size();
		if (cellSize == 8 || cellSize == 3 || cellSize == 4)
		{
			for (size_t j = 0; j < m_mesh.C.at(i).Vids.size(); j++)
			{
				const Cell& cell = m_mesh.C.at(i);
				c.Vids.push_back(v_v[cell.Vids.at(j)]);
			}
			C.push_back(c);
		}
	}
	C.resize(C.size());
	m_mesh.C = C;
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
    ofs << "CELL_DATA " << E.size() << std::endl
        << "SCALARS " << " Singularity" << " int 1" << std::endl
        << "LOOKUP_TABLE default" << std::endl;
    for (size_t i = 0; i < E.size(); i++)
        ofs << (E.at(i).isSingularity ? 1 : 0) << std::endl;

    ofs << "SCALARS " << " ESingularity" << " float" << std::endl
        << "LOOKUP_TABLE default" << std::endl;
    for (size_t i = 0; i < E.size(); i++)
        ofs << E.at(i).energySingularity << std::endl;

    ofs << "SCALARS " << " EOrthogonality" << " float" << std::endl
        << "LOOKUP_TABLE default" << std::endl;
    for (size_t i = 0; i < E.size(); i++)
        ofs << E.at(i).energyOrthogonality << std::endl;

    ofs << "SCALARS " << " EStraightness" << " float" << std::endl
        << "LOOKUP_TABLE default" << std::endl;
    for (size_t i = 0; i < E.size(); i++)
        ofs << E.at(i).energyStraightness << std::endl;
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

void MeshFileWriter::WriteEdgesVtk(const std::vector<size_t>& edgeIds)
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
    }
    else if (m_mesh.m_cellType == TETRAHEDRA || m_mesh.m_cellType == TRIANGLE) {
        ofs << "POLYGONS " << faceIds.size() << " " << 4 * faceIds.size() << std::endl;
        for (size_t i = 0; i < faceIds.size(); i++) {
            ofs << "3";
            const Face& face = F.at(faceIds.at(i));
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
        for (size_t i = 0; i < faceIds.size(); i++)
            ofs << (F.at(faceIds.at(i)).label == MAXID ? 0 : F.at(faceIds.at(i)).label)  << std::endl;
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
