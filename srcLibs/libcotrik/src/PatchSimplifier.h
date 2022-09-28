#ifndef PATCH_SIMPLIFIER_H
#define PATCH_SIMPLIFIER_H

#include "Simplifier.h"
#include "PatchSimplifierOperations.h"
// #include "AngleBasedSmoothQuadMesh.h"
class PatchSimplifier : public Simplifier {
public:
	PatchSimplifier(Mesh& mesh);
	virtual ~PatchSimplifier();
private:
	PatchSimplifier();
	PatchSimplifier(const PatchSimplifier&);
	PatchSimplifier& operator = (const PatchSimplifier&);
public:
	void Run();
	bool Simplify(int& iter);
	bool SimplifyMesh(int& iter);
	bool CheckCorners();
	void AngleBasedSmoothing(std::vector<glm::dvec3>& delta_coords);
	void ResampleBoundaryVertices(std::vector<glm::dvec3>& delta_coords);
	void RemapBoundaryVertices(std::vector<glm::dvec3>& delta_coords);
	double GetMeshEnergy();
	double GetVertexEnergy(int vid);
	void SmoothBoundary();
	void SmoothMesh();
	void RefineMesh();
	void SetOriginalRefinedMesh();
	void SmoothMesh(bool smoothGlobal_);
	void SetPosition(Vertex& v);
	void SetPositionBoundary(Vertex& v);

	void VertexRotate(std::multiset<SimplificationOperationStruct, bool(*)(SimplificationOperationStruct, SimplificationOperationStruct)>& SimplificationOps);
	void EdgeRotate(std::multiset<SimplificationOperationStruct, bool(*)(SimplificationOperationStruct, SimplificationOperationStruct)>& SimplificationOps);
	void EdgeCollapse(std::multiset<SimplificationOperationStruct, bool(*)(SimplificationOperationStruct, SimplificationOperationStruct)>& SimplificationOps);
	void DiagonalCollapse(std::multiset<SimplificationOperationStruct, bool(*)(SimplificationOperationStruct, SimplificationOperationStruct)>& SimplificationOps);
	void GetOperations(std::multiset<SimplificationOperationStruct, bool(*)(SimplificationOperationStruct, SimplificationOperationStruct)>& SimplificationOps, std::string OpType);
	void PerformOperations(std::multiset<SimplificationOperationStruct, bool(*)(SimplificationOperationStruct, SimplificationOperationStruct)>& SimplificationOps, std::set<size_t>& canceledFids);

	std::vector<Vertex> refinedV;
	double refinementFactor = 1;
	std::vector<size_t> smoothVids;
	std::vector<size_t> origBoundaryVids;
	bool smoothGlobal = false;

	int originalFaces = 0;
	// SmoothAlgorithm* smoothing_algorithm;

	std::map<size_t, std::set<size_t>> origLabel_vids;
	std::map<size_t, std::set<size_t>> origPatch_vids;
	std::map<size_t, std::set<size_t>> origSharpEdgeVid_NVids;

	std::map<size_t, std::set<size_t>> label_vids;
	std::map<size_t, std::set<size_t>> sharpEdgeVid_NVids;
	std::vector<Vertex> centerVertices;
	
	
};

#endif // !PATCH_SIMPLIFIER_H
