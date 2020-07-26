#include <fstream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;
#include <string.h>

#include "MeshFileReader.h"

#include <vtkGenericDataObjectReader.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkOBJReader.h>
#include <vtkSTLReader.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkDataSetAttributes.h>
#include <vtkSmartPointer.h>

unsigned int GetObjectNumberFromStream(ifstream& file, const char* strObjectKeyWord)
{
	if (strObjectKeyWord == NULL)
	{
		unsigned int num = 0;
		file >> num;
		return num;
	}

	string str;
	// read VerticesNumber
	unsigned int objectNum = 0;
	while (getline(file, str))
	{
		if (str.find(strObjectKeyWord) != str.npos)
		{
			bool bBeginFlag = false;
			//bool bEndFlag = false;
			for (unsigned int i = strlen(strObjectKeyWord); i < str.size(); i++)
			{
				char c = str.at(i);
				if (c >= '0' && c <= '9')
				{
					bBeginFlag = true;
					objectNum = objectNum*10 + (c - '0');
				}
				else
				{
					if (bBeginFlag)
					{
						//bEndFlag = true;
						break;
					}
				}
			}
			break;
		}
	}

	if (objectNum == 0)
	{
		file >> objectNum;
		getline(file, str);
	}
	return objectNum;
}


MeshFileReader::~MeshFileReader()
{
}

MeshFileReader::MeshFileReader(const char* pFileName)
: m_strFileName(pFileName)
{
    if (m_strFileName.find(".vtk") != m_strFileName.npos)
        ReadVtkFile();
    else if (m_strFileName.find(".off") != m_strFileName.npos)
        ReadOffFile();
    else if (m_strFileName.find(".mesh") != m_strFileName.npos)
        ReadMeshFile();
    else if (m_strFileName.find(".obj") != m_strFileName.npos)
        ReadObjFile();
    else if (m_strFileName.find(".stl") != m_strFileName.npos)
        ReadStlFile();
//    if (m_strFileName.find(".ply"))
//        ReadPlyFile();
}
void MeshFileReader::ReadVtkFile()
{
    // Get all data from the file
    vtkSmartPointer<vtkGenericDataObjectReader> reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
    reader->SetFileName(m_strFileName.c_str());
    reader->Update();

    // All of the standard data types can be checked and obtained like this:
    if (reader->IsFilePolyData())
    {
        vtkPolyData* output = reader->GetPolyDataOutput();
        const vtkIdType vnum = output->GetNumberOfPoints();
        const vtkIdType cnum = output->GetNumberOfPolys();
        std::cout << m_strFileName << " PolyData: " << vnum << " points " << cnum << " polys" << std::endl;
        std::vector<Vertex>& V = m_mesh.V;
        V.resize(vnum);
        double p[3];
        for (vtkIdType i = 0; i < vnum; i++) {
            output->GetPoint(i, p);
            V.at(i).x = p[0];
            V.at(i).y = p[1];
            V.at(i).z = p[2];
            V.at(i).id = i;
        }
//        vtkSmartPointer<vtkCellTypes> cellTypes = vtkSmartPointer<vtkCellTypes>::New();
//        output->GetCellTypes(cellTypes.GetPointer());
        m_mesh.m_cellType = POLYGON;
        std::vector<Cell>& C = m_mesh.C;
        for (vtkIdType i = 0; i < cnum; i++) {
            vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
            output->GetPolys()->GetNextCell(idList);
            const vtkIdType csize = idList->GetNumberOfIds();
            Cell c(csize);
            for (vtkIdType j = 0; j < csize; j++)
                c.Vids.at(j) = idList->GetId(j);
            c.id = i;
            if (csize == 3) c.cellType = VTK_TRIANGLE;
            else if (csize == 4) c.cellType = VTK_QUAD;
            else c.cellType = VTK_POLYGON;
            C.push_back(c);
        }
        bool hasTriangle = HasCellType(VTK_TRIANGLE);
        bool hasQuad = HasCellType(VTK_QUAD);
        if (hasTriangle && !hasQuad) m_mesh.m_cellType = TRIANGLE;
        if (!hasTriangle && hasQuad) m_mesh.m_cellType = QUAD;
        C.resize(C.size());
    }
    else if (reader->IsFileUnstructuredGrid())
    {
        vtkUnstructuredGrid* output = reader->GetUnstructuredGridOutput();
        const vtkIdType vnum = output->GetNumberOfPoints();
        const vtkIdType cnum = output->GetNumberOfCells();
        std::cout << m_strFileName << " UnstructuredGrid: " << vnum << " points " << cnum << " cells" << std::endl;
        std::vector<Vertex>& V = m_mesh.V;
        V.resize(vnum);
        m_mesh.m_cellTypes.resize(cnum);
        // Read V
        double p[3];
        for (vtkIdType i = 0; i < vnum; i++) {
            output->GetPoint(i, p);
            V.at(i).x = p[0];
            V.at(i).y = p[1];
            V.at(i).z = p[2];
            V.at(i).id = i;
            V.at(i).cellType = VTK_VERTEX;
        }
        // Read C
        std::vector<Cell>& C = m_mesh.C;
        for (vtkIdType i = 0; i < cnum; i++) {
            vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
            output->GetCellPoints(i, idList);
            const vtkIdType csize = idList->GetNumberOfIds();
            Cell c(csize);
            c.id = i;
            if (csize == 4) c.cellType = VTK_TETRA;
			else if (csize == 6) c.cellType = VTK_WEDGE;
            else if (csize == 8) c.cellType = VTK_HEXAHEDRON;
			else if (csize == 10) c.cellType = VTK_PENTAGONAL_PRISM;
			else if (csize == 12) c.cellType = VTK_HEXAGONAL_PRISM;
            else c.cellType = VTK_POLYHEDRON;
            for (vtkIdType j = 0; j < csize; j++)
                c.Vids.at(j) = idList->GetId(j);
            C.push_back(c);
        }
        C.resize(C.size());
        // Read CellType
        const vtkIdType cellType = output->GetCellType(0);
        if (cellType == VTK_TRIANGLE) m_mesh.m_cellType = TRIANGLE;
        else if (cellType == VTK_QUAD) m_mesh.m_cellType = QUAD;
        else if (cellType == VTK_TETRA) m_mesh.m_cellType = TETRAHEDRA;
        else if (cellType == VTK_HEXAHEDRON) m_mesh.m_cellType = HEXAHEDRA;
		for (vtkIdType i = 0; i < cnum; i++) {
			m_mesh.m_cellTypes[i] = output->GetCellType(i);
			m_mesh.C[i].cellType = (VTKCellType)output->GetCellType(i);
		}
        for (vtkIdType i = 0; i < cnum; i++)
            if (output->GetCellType(i) != cellType) {
                m_mesh.m_cellType = POLYHEDRA;
                break;
            }
        if (m_mesh.m_cellType == TRIANGLE || m_mesh.m_cellType == QUAD) {
            m_mesh.F.resize(C.size());
            for (vtkIdType i = 0; i < cnum; i++)
                m_mesh.F[i].Vids = m_mesh.C[i].Vids;
        }
    }
}

void MeshFileReader::ReadObjFile()
{
    // Get all data from the file
    vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
    reader->SetFileName(m_strFileName.c_str());
    reader->Update();

    vtkPolyData* output = reader->GetOutput();
    const vtkIdType vnum = output->GetNumberOfPoints();
    const vtkIdType cnum = output->GetNumberOfPolys();
    std::cout << m_strFileName << " PolyData: " << vnum << " points " << cnum << " polys" << std::endl;
    std::vector<Vertex>& V = m_mesh.V;
    V.resize(vnum);
    double p[3];
    for (vtkIdType i = 0; i < vnum; i++)
    {
        output->GetPoint(i, p);
        V.at(i).x = p[0];
        V.at(i).y = p[1];
        V.at(i).z = p[2];
        V.at(i).id = i;
    }
    // Read CellType
    const vtkIdType cellType = output->GetCellType(0);
    if (cellType == VTK_TRIANGLE) m_mesh.m_cellType = TRIANGLE;
    else if (cellType == VTK_QUAD) m_mesh.m_cellType = QUAD;
    else if (cellType == VTK_POLYGON) m_mesh.m_cellType = POLYGON;

    std::vector<Cell>& C = m_mesh.C;
    for (vtkIdType i = 0; i < cnum; i++)
    {
        vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
        output->GetCellPoints(i, idList);
        const vtkIdType csize = idList->GetNumberOfIds();
        Cell c(csize);
        for (vtkIdType j = 0; j < csize; j++)
            c.Vids.at(j) = idList->GetId(j);
        C.push_back(c);
    }
    C.resize(C.size());

    if (m_mesh.m_cellType == TRIANGLE || m_mesh.m_cellType == QUAD) {
        m_mesh.F.resize(C.size());
        for (vtkIdType i = 0; i < cnum; i++)
            m_mesh.F[i].Vids = m_mesh.C[i].Vids;
    }
}

void MeshFileReader::ReadStlFile()
{
    // Get all data from the file
    vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName(m_strFileName.c_str());
    reader->Update();

    vtkPolyData* output = reader->GetOutput();
    const vtkIdType vnum = output->GetNumberOfPoints();
    const vtkIdType cnum = output->GetNumberOfPolys();
    std::cout << m_strFileName << " PolyData: " << vnum << " points " << cnum << " polys" << std::endl;
    std::vector<Vertex>& V = m_mesh.V;
    V.resize(vnum);
    double p[3];
    for (vtkIdType i = 0; i < vnum; i++)
    {
        output->GetPoint(i, p);
        V.at(i).x = p[0];
        V.at(i).y = p[1];
        V.at(i).z = p[2];
        V.at(i).id = i;
    }
    // Read CellType
    const vtkIdType cellType = output->GetCellType(0);
    if (cellType == VTK_TRIANGLE) m_mesh.m_cellType = TRIANGLE;
    else if (cellType == VTK_QUAD) m_mesh.m_cellType = QUAD;

    std::vector<Cell>& C = m_mesh.C;
    for (vtkIdType i = 0; i < cnum; i++)
    {
        vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
        output->GetCellPoints(i, idList);
        const vtkIdType csize = idList->GetNumberOfIds();
        Cell c(csize);
        for (vtkIdType j = 0; j < csize; j++)
            c.Vids.at(j) = idList->GetId(j);
        C.push_back(c);
    }
    C.resize(C.size());

    if (m_mesh.m_cellType == TRIANGLE || m_mesh.m_cellType == QUAD) {
        m_mesh.F.resize(C.size());
        for (vtkIdType i = 0; i < cnum; i++)
            m_mesh.F[i].Vids = m_mesh.C[i].Vids;
    }
}

void MeshFileReader::ReadOffFile()
{
    std::ifstream ifs(m_strFileName.c_str());
    std::string str;
    ifs >> str;

    long long vnum;
    long long cnum;
    long long t;
    ifs >> vnum >> cnum >> t;
    std::vector<Vertex>& V = m_mesh.V;
    std::vector<Cell>& C = m_mesh.C;
    V.resize(vnum);
    getline(ifs, str);
    for (long long i = 0; i < vnum; i++) {
        getline(ifs, str);
        stringstream ss(str.c_str());
        ss >> V.at(i).x >> V.at(i).y >> V.at(i).z;
        V.at(i).id = i;
    }
//    for (long long i = 0; i < vnum; i++)
//        ifs >> V.at(i).x >> V.at(i).y >> V.at(i).z;
    int a, b;
    for (long long i = 0; i < cnum; i++) {
        int csize;
        ifs >> csize;
        // Read CellType
        Cell c(csize);
        if (csize == 3) c.cellType = VTK_TRIANGLE;
        else if (csize == 4) c.cellType = VTK_QUAD;
        else if (csize == 5) {c.cellType = VTK_TETRA; csize = 4;}
        else if (csize == 10) {c.cellType = VTK_HEXAHEDRON; csize = 8;}

        for (auto& vid : c.Vids)
            ifs >> vid;
        if (c.cellType == VTK_TETRA) ifs >> a;
        else if (c.cellType == VTK_HEXAHEDRON) ifs >> a >> b;
        C.push_back(c);
    }
    C.resize(C.size());

    m_mesh.m_cellType = POLYGON;
    bool hasTriangle = HasCellType(VTK_TRIANGLE);
    bool hasQuad = HasCellType(VTK_QUAD);
    if (hasTriangle && !hasQuad) m_mesh.m_cellType = TRIANGLE;
    else if (!hasTriangle && hasQuad) m_mesh.m_cellType = QUAD;
    else if (!hasTriangle && !hasQuad) m_mesh.m_cellType = POLYHEDRA;

    if (m_mesh.m_cellType == TRIANGLE || m_mesh.m_cellType == QUAD) {
        m_mesh.F.resize(C.size());
        for (vtkIdType i = 0; i < cnum; i++) {
            m_mesh.F[i].Vids = m_mesh.C[i].Vids;
            m_mesh.F[i].id = i;
        }
    }
}

void MeshFileReader::ReadMeshFile()
{
    std::ifstream ifs(m_strFileName.c_str());
    std::string str;
    while (ifs >> str)
        if (str == "Vertices")
            break;

    int vnum;
    ifs >> vnum;
    std::vector<Vertex>& V = m_mesh.V;
    V.resize(vnum);
    for (long long i = 0; i < vnum; i++){
        ifs >> V.at(i).x >> V.at(i).y >> V.at(i).z >> str;
        V.at(i).id = i;
    }
    // read Cells

    int cnum;
    int csize = 0;
    ifs >> str >> cnum;
    std::vector<Cell>& C = m_mesh.C;

    if (str == "Hexahedra") {
        C.resize(cnum, Cell(8));
        m_mesh.m_cellType = HEXAHEDRA;
        csize = 8;
    } else if (str == "Tetrahedra") {
        cnum = cnum / 4;
        C.resize(cnum, Cell(4));
        m_mesh.m_cellType = TETRAHEDRA;
        csize = 4;
    } else if (str == "Quadrilaterals") {
        C.resize(cnum, Cell(4));
        m_mesh.m_cellType = QUAD;
        csize = 4;
    } else if (str == "Triangles") {
        C.resize(cnum, Cell(3));
        m_mesh.m_cellType = TRIANGLE;
        csize = 3;
    }
    else if (str == "Edges") {
        while (ifs >> str)
            if (str == "Hexahedra")
                break;
        ifs >> cnum;
        C.resize(cnum, Cell(8)); m_mesh.m_cellType = HEXAHEDRA; csize = 8;}

    for (long long i = 0; i < cnum; i++)
    {
        for (int j = 0; j < csize; j++){
            ifs >> C.at(i).Vids.at(j);
            C.at(i).Vids.at(j)--;
        }
        ifs >> str;
    }

    if (m_mesh.m_cellType == TRIANGLE || m_mesh.m_cellType == QUAD) {
        m_mesh.F.resize(C.size());
        for (vtkIdType i = 0; i < cnum; i++)
            m_mesh.F[i].Vids = m_mesh.C[i].Vids;
    }
}

const Mesh& MeshFileReader::GetMesh() const
{
	return m_mesh;
}

void MeshFileReader::GetScalarFields()
{
    vtkSmartPointer<vtkUnstructuredGridReader> pReader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    pReader->SetFileName(m_strFileName.c_str());
    pReader->Update();
	vtkCellData* cellData = pReader->GetOutput()->GetCellData();
	vtkDataSetAttributes* attribute = vtkDataSetAttributes::SafeDownCast(cellData);
	vtkIntArray* scalarDataInt = vtkIntArray::SafeDownCast(attribute->GetScalars(pReader->GetScalarsNameInFile(0)));
	if (scalarDataInt)
	{
		int nc = scalarDataInt->GetNumberOfTuples();
		std::cout << "There are " << nc << " components in " << pReader->GetScalarsNameInFile(0) << std::endl;
		std::vector<double> scalarField(nc);
		for (int i = 0; i < nc; i++)
			scalarField.at(i) = scalarDataInt->GetValue(i);
		m_mesh.cellScalarFields.push_back(scalarField);
	}
	pReader->SetScalarsName(pReader->GetScalarsNameInFile(1));
	pReader->Update();
	cellData->Update();
	attribute = vtkDataSetAttributes::SafeDownCast(cellData);
	vtkFloatArray* scalarDataFloat = vtkFloatArray::SafeDownCast(attribute->GetScalars(pReader->GetScalarsNameInFile(1)));
	if (scalarDataFloat)
	{
		int nc = scalarDataFloat->GetNumberOfTuples();
		std::cout << "There are " << nc << " components in scalarDataFloat"	<< std::endl;
		std::vector<double> scalarField(nc);
		for (int i = 0; i < nc; i++)
			scalarField.at(i) = scalarDataFloat->GetValue(i);
		m_mesh.cellScalarFields.push_back(scalarField);
	}
}

void MeshFileReader::GetPointsScalarFields()
{
    vtkSmartPointer<vtkUnstructuredGridReader> pReader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    pReader->SetFileName(m_strFileName.c_str());
    pReader->Update();
	vtkPointData* pointData = pReader->GetOutput()->GetPointData();
	vtkDataSetAttributes* attribute = vtkDataSetAttributes::SafeDownCast(pointData);
	vtkIntArray* scalarDataInt = vtkIntArray::SafeDownCast(attribute->GetScalars(pReader->GetScalarsNameInFile(0)));
	if (scalarDataInt)
	{
		int nc = scalarDataInt->GetNumberOfTuples();
		std::cout << "There are " << nc << " components in " << pReader->GetScalarsNameInFile(0) << std::endl;
		std::vector<double> scalarField(nc);
		for (int i = 0; i < nc; i++)
			scalarField.at(i) = scalarDataInt->GetValue(i);
		m_mesh.pointScalarFields.push_back(scalarField);
	}
//	pReader->SetScalarsName(pReader->GetScalarsNameInFile(1));
//	pReader->Update();
//	pointData->Update();
//	attribute = vtkDataSetAttributes::SafeDownCast(pointData);
//	vtkFloatArray* scalarDataFloat = vtkFloatArray::SafeDownCast(attribute->GetScalars(pReader->GetScalarsNameInFile(1)));
//	if (scalarDataFloat)
//	{
//		int nc = scalarDataFloat->GetNumberOfTuples();
//		std::cout << "There are " << nc << " components in scalarDataFloat"	<< std::endl;
//		std::vector<double> scalarField(nc);
//		for (int i = 0; i < nc; i++)
//			scalarField.at(i) = scalarDataFloat->GetValue(i);
//		m_mesh.cellScalarFields.push_back(scalarField);
//	}
}

bool MeshFileReader::HasCellType(const VTKCellType cellType) const {
    for (auto& c : m_mesh.C)
        if (c.cellType == cellType) return true;
    return false;
}
