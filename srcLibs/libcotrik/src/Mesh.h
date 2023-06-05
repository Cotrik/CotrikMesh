/*
 * Mesh.h
 *
 *  Created on: Nov 6, 2016
 *      Author: cotrik
 */

#ifndef MESH_H_
#define MESH_H_

#include "Util.h"
#include <glm/glm.hpp>
#include <vtkCellType.h>
class FeatureLine;
enum ElementType {
    POLYGON,
    TRIANGLE,
    QUAD,
    TETRAHEDRA,
    HEXAHEDRA,
    POLYHEDRA,
	WEDGE,
	PENTAHEDRON,
};

extern const size_t MAXID;
extern const unsigned int HexEdges[12][2];
extern const unsigned int HexFaces[6][4];

extern const unsigned int TetEdges[6][2];
extern const unsigned int TetFaces[4][3];

extern const unsigned int WedgeEdges[9][2];
extern const unsigned int PentaEdges[15][2];

extern const std::vector<std::vector<size_t>> HexFace;
extern const std::vector<std::vector<size_t>> TetFace;
extern const std::vector<std::vector<size_t>> WedgeFace;
extern const std::vector<std::vector<size_t>> PentaFace;

extern const std::vector<std::vector<size_t>> HexEdge;
extern const std::vector<std::vector<size_t>> TetEdge;
extern const std::vector<std::vector<size_t>> WedgeEdge;
extern const std::vector<std::vector<size_t>> PentaEdge;

extern const unsigned int TriEdge[3][2];
extern const unsigned int QuadEdge[4][2];

extern const unsigned int HexPoint_Points[8][3];
extern const unsigned int HexPoint_Edges[8][3];

extern const std::map<size_t, std::vector<std::vector<size_t>>> cell_faces;
extern const std::map<size_t, std::vector<std::vector<size_t>>> cell_edges;

class NeighborInfo
{
public:
    NeighborInfo(){}
    NeighborInfo(const NeighborInfo& rhs)
    : N_Vids(rhs.N_Vids)
    , N_Eids(rhs.N_Eids)
    , N_Fids(rhs.N_Fids)
    , N_Cids(rhs.N_Cids)
    {}
    NeighborInfo& operator = (const NeighborInfo& rhs)
    {
        N_Vids = rhs.N_Vids;
        N_Eids = rhs.N_Eids;
        N_Fids = rhs.N_Fids;
        N_Cids = rhs.N_Cids;
        return *this;
    }
    ~NeighborInfo(){}
    std::vector<size_t> N_Vids;  // neighboring vertices ids
    std::vector<size_t> N_Eids;  // neighboring edges ids
    std::vector<size_t> N_Fids;  // neighboring faces ids
    std::vector<size_t> N_Cids;  // neighboring cells ids
};

struct GeoInfo
{
public:
    GeoInfo()
    : id(MAXID)
    , component_id(MAXID)
    , color(-1)
    , isBoundary(false)
    , isSingularity(false)
    , isActive(true)
    , cellType(VTK_HEXAHEDRON)
    {}
    GeoInfo(const GeoInfo& r)
    : id(r.id)
    , component_id(r.component_id)
    , color(r.color)
    , isBoundary(r.isBoundary)
    , isSingularity(r.isSingularity)
    , isActive(r.isActive)
    , cellType(r.cellType)
    {}
    GeoInfo& operator = (const GeoInfo& rhs)
    {
        id = rhs.id;
        component_id = rhs.component_id;
        color = rhs.color;
        isBoundary = rhs.isBoundary;
        isSingularity = rhs.isSingularity;
        isActive = rhs.isActive;
        cellType = rhs.cellType;
        return *this;
    }
    virtual ~GeoInfo()
    {}
public:
    size_t id;
    size_t component_id;
    char color;
    bool isBoundary;
    bool isSingularity;
    bool isActive;
    VTKCellType cellType;
};

enum SurfaceVertexType
{
    REGULAR = 0,
    FEATURE,
    CORNER
};

enum SmoothMethod
{
    LAPLACE_EDGE = 0,
    LAPLACE_FACE_CENTER,
    LAPLACE_FACE_CENTER_ANGLE,
    LAPLACE_BELTRAMI
};

class Vertex : public glm::dvec3, public GeoInfo, public NeighborInfo
{
public:
    Vertex()
    : hvid(MAXID)
    , triVid(MAXID)
    , type(REGULAR)
	, label(MAXID)
	, patch_id(MAXID)
    , isCorner(false)
	, isSpecial(false)
    , isConvex(false)
    , idealValence(0)
    {}
    Vertex(const Vertex& r)
    : glm::dvec3(r)
    , GeoInfo(r)
    , NeighborInfo(r)
    , normal(r.normal)
    , tangent(r.tangent)
    , hvid(r.hvid)
    , triVid(r.triVid)
    , type(r.type)
	, label(r.label)
	, patch_id(r.patch_id)
	, labels(r.labels)
	, patch_ids(r.patch_ids)
    , isCorner(r.isCorner)
	, isSpecial(r.isSpecial)
    , isConvex(r.isConvex)
    , idealValence(r.idealValence)
    {}
    Vertex(const glm::dvec3& v)
    : glm::dvec3(v)
    , hvid(MAXID)
    , triVid(MAXID)
    , type(REGULAR)
	, label(MAXID)
	, patch_id(MAXID)
    , isCorner(false)
	, isSpecial(false)
    , isConvex(false)
    {}
    virtual ~Vertex()
    {}
	Vertex& operator = (const Vertex& r) {
		if (r == *this)
			return *this;
		x = r.x;
		y = r.y;
		z = r.z;

		return *this;
	}
	Vertex& operator = (const glm::dvec3& r) {
		if (r == *this)
			return *this;
		x = r.x;
		y = r.y;
		z = r.z;

		return *this;
	}
	glm::dvec3 xyz() const {
		return glm::dvec3(x, y, z);
	}
    void xyz(glm::dvec3 r) {
		x = r.x;
        y = r.y;
        z = r.z;
	}
    glm::dvec3 normal;

    // for edge cone triangle surface
    glm::dvec3 tangent;
    size_t hvid;
    size_t triVid;
    size_t type; // 0 regular, 1 feature, 2 corner.
	size_t label; // the same as the edge line label;
	size_t patch_id; // face patch id;
	std::set<size_t> labels; // for corner
	std::set<size_t> patch_ids;
    bool isCorner;
	bool isSpecial;
	bool isConvex;
    bool isVisited;
    bool isMovable = true;
	int idealValence = 0;
    std::vector<size_t> twoRingNeighborSurfaceFaceIds; // for surface projecting;
    std::vector<size_t> oneRingNeighborVertices; // for 2D surface smoothing

    double prescribed_length;
    bool smoothLocal;
};

class Edge : public GeoInfo, public NeighborInfo
{
public:
    Edge()
    : length(0.0)
    , energySingularity(0.0)
    , energyOrthogonality(0.0)
    , energyStraightness(0.0)
    , face_angle(0.0)
    , isSharpFeature(false)
    , label(MAXID)
    {}
    Edge(const Edge& r)
    : GeoInfo(r)
    , NeighborInfo(r)
    , Vids(r.Vids)
    , parallelEids(r.parallelEids)
    , consecutiveEids(r.consecutiveEids)
    , orthogonalEids(r.orthogonalEids)
    , length(r.length)
    , energySingularity(r.energySingularity)
    , energyOrthogonality(r.energyOrthogonality)
    , energyStraightness(r.energyStraightness)
    , face_angle(r.face_angle)
    , isSharpFeature(r.isSharpFeature)
    , label(r.label)
    {}
    Edge(size_t vnum)
    : length(0.0)
    , energySingularity(0.0)
    , energyOrthogonality(0.0)
    , energyStraightness(0.0)
    , face_angle(0.0)
    , isSharpFeature(false)
    , label(MAXID)
    {
        Vids.resize(vnum);
    }
    Edge(const std::vector<size_t>& Vids)
    : Vids(Vids)
    , length(0.0)
    , energySingularity(0.0)
    , energyOrthogonality(0.0)
    , energyStraightness(0.0)
    , face_angle(0.0)
    , isSharpFeature(false)
    , label(MAXID)
    {
    }
    virtual ~Edge()
    {}
public:
    bool operator == (const Edge& e) const {
        return ((Vids[0] == e.Vids[0] && Vids[1] == e.Vids[1]) || (Vids[0] == e.Vids[1] && Vids[1] == e.Vids[0]) );
    }
	bool operator < (const Edge& e) const {
		if (*this == e) return false;
		return Vids[0] < e.Vids[0];
	}
public:
    std::vector<size_t> Vids;
    std::vector<size_t> parallelEids;
    std::vector<size_t> consecutiveEids;
    std::vector<size_t> orthogonalEids;
    double length;                             // target length that is used in frameOpt
    double energySingularity;
    double energyOrthogonality;
    double energyStraightness;

    // for edge cone triangle surface
    double face_angle;

    bool isSharpFeature;
    size_t label;     // patch number starting from 0, MAXID is invalid
    size_t componentEid = MAXID;
    size_t singularEid = MAXID;
    bool isVisited = false;
};

class Face : public GeoInfo, public NeighborInfo
{
public:
    Face()
    : label(MAXID)
    {}
    Face(const Face& r)
    : GeoInfo(GeoInfo(r))
    , NeighborInfo(NeighborInfo(r))
    , Vids(r.Vids)
    , Eids(r.Eids)
    , normal(r.normal)
    , label(r.label)
    {}
    Face(size_t vnum)
    : label(MAXID)
    {
        Vids.resize(vnum);
    }
    Face(size_t vnum, size_t eNum)
    : label(MAXID)
    {
        Vids.resize(vnum);
        Eids.resize(eNum);
    }
    Face(const std::vector<size_t>& Vids)
    : Vids(Vids)
    , label(MAXID)
    {
    }

    Face(const std::vector<size_t>& Vids, const std::vector<size_t> Eids)
    : Vids(Vids)
    , Eids(Eids)
    , label(MAXID)
    {
    }
    virtual ~Face()
    {}

public:
    std::vector<size_t> Vids;
    std::vector<size_t> Eids;
    std::vector<std::vector<size_t> > N_Ortho_4Vids;   // Neighboring Orthogonal 4 Vids of Edges
    std::vector<std::vector<size_t> > N_Ortho_4Eids;   // Neighboring Orthogonal 4 Eids of Edges

    glm::dvec3 normal;
    size_t label;     // patch number starting from 0, MAXID is invalid
    size_t componentFid = MAXID;
    bool isNegative = false;
    bool isVisited = false;
};

class Cell : public GeoInfo, public NeighborInfo
{
public:
    Cell()
    {}
    Cell(const Cell& r)
    : GeoInfo(GeoInfo(r))
    , NeighborInfo(NeighborInfo(r))
    , Vids(r.Vids)
    , Eids(r.Eids)
    , Fids(r.Fids)
    {}
    Cell(size_t vnum)
    {
        Vids.resize(vnum);
    }
    Cell(size_t vnum, size_t eNum, size_t fNum)
    {
        Vids.resize(vnum);
        Eids.resize(eNum);
        Fids.resize(fNum);
    }
    Cell(const std::vector<size_t>& Vids)
    : Vids(Vids)
    {
    }

    Cell(const std::vector<size_t>& Vids, const std::vector<size_t> Eids, const std::vector<size_t> Fids)
    : Vids(Vids)
    , Eids(Eids)
    , Fids(Fids)
    {
    }

    virtual ~Cell()
    {}

public:
    std::vector<size_t> Vids;
    std::vector<size_t> Eids;
    std::vector<size_t> Fids;
    //glm::dvec3 cv;                 // center vertex's x,y,z coordinate of the Cell;  For Node of Frame
    size_t componentCid = MAXID;
    double qualityValue = 0;
};

class Layer
{
public:
    Layer(){}
    Layer(const Layer& layer)
    : Vids(layer.Vids)
    , Eids(layer.Eids)
    , Fids(layer.Fids)
    , Cids(layer.Cids)
    , fixed(layer.fixed)
    {}
    ~Layer(){}
    std::vector<size_t> Vids;
    std::vector<size_t> Eids;
    std::vector<size_t> Fids;
    std::vector<size_t> Cids;
    std::vector<bool>  fixed;  // refer to this->Vids;
};

class Plane
{
public:
    Plane(const double a, const double b, const double c, const double d)
        : a(a), b(b), c(c), d(d) {
    }
    Plane(const Vertex& p1, const Vertex& p2, const Vertex& p3) {
        a = ((p2.y - p1.y) * (p3.z - p1.z) - (p2.z - p1.z) * (p3.y - p1.y));
        b = ((p2.z - p1.z) * (p3.x - p1.x) - (p2.x - p1.x) * (p3.z - p1.z));
        c = ((p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x));
        d = (0 - (a * p1.x + b * p1.y + c * p1.z));
    }
    Plane(const Plane& plane)
    : a(plane.a), b(plane.b), c(plane.c), d(plane.d) {
    }
    ~Plane(){}
    double IntersectionAngle(const Plane& p) const  {
        double cosangle = (a*p.a + b*p.b + c*p.c) / (sqrt(a*a + b*b + c*c)* sqrt(p.a*p.a + p.b*p.b + p.c*p.c));
        return cosangle;
    }
    double DistanseFromPoint(const glm::dvec3& p, glm::dvec3& intersection) const {
        const double distance = fabs(a*p.x + b*p.y + c*p.z + d) / sqrt(a*a + b*b + c*c);
        const double t = (a*p.x + b*p.y + c*p.z + d) / (a*a + b*b + c*c);
        intersection = glm::dvec3(p.x - a*t, p.y - b*t, p.z - c*t);
        return distance;
    }
    bool operator < (const Plane& right) const {
        return d < right.d;
    }

public:
    double a;
    double b;
    double c;
    double d;
};

class Line
{
public:
    Line(const glm::dvec3& o, const glm::dvec3& d)
        : o(o), dir(d) {
    }
    Line(const Vertex& p1, const Vertex& p2)
        : o(p1.xyz()), dir(p2.xyz() - p1.xyz()) {
    }
    Line(const Line& r)
    : o(r.o), dir(r.dir)    {
    }
    double Perpendicular(const glm::dvec3& a, glm::dvec3& intersection) const {
        const double t = ((a.x - o.x) * dir.x + (a.y - o.y) * dir.y + (a.z - o.z) * dir.z) / (dir.x * dir.x + dir.y * dir.y + dir.z * dir.z);
        intersection.x = o.x + dir.x * t;
        intersection.y = o.y + dir.y * t;
        intersection.z = o.z + dir.z * t;
        return glm::length(intersection - a);
    }
    ~Line(){}
//    bool operator < (const Plane& right) const
//    {
//        return d < right.d;
//    }
public:
    glm::dvec3 o;
    glm::dvec3 dir;
};

class Mesh {
public:
    Mesh();
    Mesh(const Mesh& r);
    Mesh(const std::vector<Vertex>& V, const std::vector<Cell>& C, ElementType m_cellType);
    Mesh(const std::vector<Vertex>& V, const std::vector<Face>& F, ElementType m_cellType = QUAD);
    Mesh(const Mesh& r, const std::vector<size_t>& cellIds);
    virtual ~Mesh();

public:
    void BuildAllConnectivities(); // Get Neighboring Info, including, E, F, C, V_V, V_E, V_F, V_C, E_V, E_F, E_C, F_V, F_E, F_F, F_C, C_V, C_E, C_F, C_C
    void ExtractBoundary();
    inline bool HasBoundary() const;
	inline bool IsSurfaceMesh() const;
	inline bool IsVolumetricMesh() const;
    size_t ExtractLayers();
    void ExtractSingularities();
    void ExtractTwoRingNeighborSurfaceFaceIdsForEachVertex(int N = 2);
    void BuildParallelE();
    void BuildConsecutiveE();
    void BuildOrthogonalE();
    void ReOrientSurfaceFaces();
    double getConvexVerdict(std::vector<size_t> Vids);
    void unifyOrientation();
    void GetNormalOfSurfaceFaces();             // must ExtractBoundary(); first
    void GetNormalOfSurfaceVertices();          // must ExtractBoundary(); first
    void RemoveUselessVertices();
    void CompressWithFeaturePreserved();
    void ClassifyVertexTypes();
    void LabelSurface();
	void Label2DSurfaceVertices();
    void ClearLabelOfSurface();
    void LabelSharpEdges(const bool breakAtConrer = false);
    void ProjectSharpEdgesTo(const std::vector<FeatureLine>& featureLines);
    void SetCosAngleThreshold(const double cos_angle = 0.91);
	void SetFeatureAngleThreshold(const double angle = 170.0);
    //bool IsPointInTriangle(const glm::dvec3& P, const glm::dvec3& A, const glm::dvec3& B, const glm::dvec3& C) const;
    bool IsPointInTriangle(const Vertex& p, const Vertex& p0, const Vertex& p1, const Vertex& p2) const;
    bool IsPointInTriangle(const glm::dvec3& p, const glm::dvec3& p0, const glm::dvec3& p1, const glm::dvec3& p2) const;
    bool IsPointInFace(const glm::dvec3& P, const Face& face) const;
    glm::dvec3 GetProjectLocation(const glm::dvec3& p) const;
    glm::dvec3 GetProjectLocation(const Vertex& p) const;
    glm::dvec3 GetProjectLocationFast(const Vertex& p) const;
    glm::dvec3 GetProjectLocationOnRefSurface(const glm::dvec3& p, const Vertex& refV) const;
    glm::dvec3 GetProjectLocationOnTargetSurface(const glm::dvec3& p, const Vertex& refV, const Mesh& targetSurfaceMesh) const;
    void ProjectTo(const Mesh& mesh);
    void FastProjectTo(const Mesh& mesh);
    void ProjectToRefMesh(const Mesh& refMesh);
    void ProjectToTargetSurface(const Mesh& refMesh, const Mesh& targetSurfaceMesh);

    double GetAvgEdgeLength();
    double SmoothVolume(const SmoothMethod smoothMethod = LAPLACE_EDGE);
    double SmoothSurface(size_t iters = 1, const SmoothMethod smoothMethod = LAPLACE_EDGE,
            const bool preserveSharpFeature = false, const bool treatSharpFeatureAsRegular = false, const bool treatCornerAsRegular = false);
    double SmoothAndProjectSurface(const Mesh& mesh, size_t iters = 1, const SmoothMethod smoothMethod = LAPLACE_EDGE,
            const bool preserveSharpFeature = false, const bool treatSharpFeatureAsRegular = false, const bool treatCornerAsRegular = false,
            const bool preserveQuality = false);
    double ProjectSurface(const Mesh& mesh, size_t iters = 1, const SmoothMethod smoothMethod = LAPLACE_EDGE,
            const bool preserveSharpFeature = false, const bool treatSharpFeatureAsRegular = false, const bool treatCornerAsRegular = false,
            const bool preserveQuality = false);
    double SmoothVolume(const Mesh& mesh, size_t iters = 1, const SmoothMethod smoothMethod = LAPLACE_EDGE,
            const bool preserveQuality = false);
    double SmoothVolume(size_t iters = 1, const SmoothMethod smoothMethod = LAPLACE_EDGE, const bool preserveQuality = false);
    glm::dvec3 GetFaceCenter(const Face& f);
    glm::dvec3 LapLace(const Vertex& v, const bool treatSharpFeatureAsRegular = false, const bool treatCornerAsRegular = false);
    void ConvertSurfaceToTriangleMesh();

    void Zoom(const glm::dvec3& ref, const double scale = 1);
    const float GetScaledJacobian(const Cell& c) const;
    double GetMinScaledJacobian(double& avgSJ) const;
    bool IsPointInside(const glm::dvec3& orig, const glm::dvec3 dir = glm::dvec3(0, 0, 1)) const;
protected:
    virtual void BuildE();
    virtual void BuildF();
    // -------------
    virtual void BuildV_V();
    virtual void BuildV_E();
    virtual void BuildV_F();
    virtual void BuildV_C();
    // -------------
    virtual void BuildE_V();
    virtual void BuildE_E();
    virtual void BuildE_F();
    virtual void BuildE_C();
    // -------------
    virtual void BuildF_V();
    virtual void BuildF_E();
    virtual void BuildF_F();
    virtual void BuildF_C();
    // -------------
    virtual void BuildC_V();
    virtual void BuildC_E();
    virtual void BuildC_F();
    virtual void BuildC_C();
	// -------------
    void LabelFace(Face& face, size_t& label);
    void LabelEdge(Edge& edge, size_t& label, const bool breakAtConrer = false);
public:
	virtual void FixOrientation();
    double GetCosAngle(const Edge& edge, const Face& face1, const Face& face2);
	double Get2DAngle(const Vertex& v, const Edge& e1, const Edge& e2);
public:
    void GetQuality(const Vertex& v, double& minSJ, double& avgSJ);
    size_t GetQuality(const char* filename, double& minValue, double& avgValue, const double minSJ = 0.0);
    size_t GetQualityVerdict(double& minValue, double& avgValue, const double minSJ = 0.0);
    void OutputBadCells(const std::vector<size_t>& badCellIds, const char* filename);
public:
    void SetOneRingNeighborhood(); // 2D surface smoothing
    void ArrangeFaceVerticesAntiClockwise();
    void ExtractOneRingNeighbors(Vertex& source);
    double GetQuadFaceArea(std::vector<size_t>& Vids);
    double GetQuadMeshArea();
    void SetIdealValence(size_t vid);
    double GetAngle(size_t vid, size_t vid1, size_t vid2);
public:
    std::vector<Vertex> V;
    std::vector<Edge> E;
    std::vector<Face> F;
    std::vector<Cell> C;
    ElementType m_cellType;
    std::vector<size_t> m_cellTypes;

    std::vector<std::string> pointScalarFieldNames;
    std::vector<std::string> cellScalarNameFields;
    std::vector<std::vector<double> > pointScalarFields;
    std::vector<std::vector<double> > cellScalarFields;

    std::vector<Layer> layers;
    std::vector<Layer> innerLayers;

    double avgEdgeLength;
    size_t numOfSharpEdges;

    glm::dvec3 m_center;
    std::vector<size_t> m_refIds;
    double cos_angle_threshold = 0.984807753;
	double feature_angle_threshold = 170.0;
    size_t numberOfPatches = 0;

    bool hasBoundary = false;
    bool isManifold = true;
    bool smoothGlobal = true;
    double totalArea = 0.0;
    double prescribed_length = 0.0;
    bool isPlanar = false;
};

#endif /* MESH_H_ */
