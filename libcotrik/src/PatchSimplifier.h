#ifndef PATCH_SIMPLIFIER_H
#define PATCH_SIMPLIFIER_H

#include "Simplifier.h"
//#include "AngleBasedSmoothQuadMesh.h"
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
	bool SimplifyCollective(int& iter);
	bool ChordCollapseSimplify(int& iter);
	bool SeparatrixSplitSimplify(int& iter);
	bool SeparatrixCollapseSimplify(int& iter);
	bool HalfSeparatrixCollapseSimplify(int& iter);
	bool EdgeRotateSimplify(int& iter);
	bool CheckCorners();
	void smoothMesh(int iters_, bool global);
//	SmoothAlgorithm* smoothing_algorithm;
};

#endif // !PATCH_SIMPLIFIER_H
