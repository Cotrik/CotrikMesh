#ifndef PATCH_SIMPLIFIER_H
#define PATCH_SIMPLIFIER_H

#include "Simplifier.h"
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
	bool CheckCorners();
	void AngleBasedSmoothing(std::vector<glm::dvec3>& delta_coords);
	void ResampleBoundaryVertices(std::vector<glm::dvec3>& delta_coords);
	void RemapBoundaryVertices(std::vector<glm::dvec3>& delta_coords);
	double GetMeshEnergy();
	double GetVertexEnergy(int vid);
	bool IsFaceNegative(int fid, int vid, glm::dvec3 false_coord);
	void SmoothBoundary();
	void SmoothMesh();

	// std::vector<Vertex> refinedV;
	double refinementFactor = 0.5;
	std::vector<size_t> smoothVids;
	std::vector<size_t> origBoundaryVids;
	bool smoothGlobal = false;

	int originalFaces = 0;
	// SmoothAlgorithm* smoothing_algorithm;
};

#endif // !PATCH_SIMPLIFIER_H
