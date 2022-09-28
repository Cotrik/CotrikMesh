#ifndef DUAL_MESH_H
#define DUAL_MESH_H

#include "Mesh.h"
class DualMesh : public Mesh {
public:
	DualMesh();
	virtual ~DualMesh();
public:
	void Build(const Mesh& mesh);
	void BuildV(const Mesh& mesh);
	void BuildBoundaryV(const Mesh& mesh);
	void BuildF(const Mesh& mesh, std::map<size_t, size_t>& boundaryFid_dualVid, std::map<size_t, size_t>& vid_dualFid);
	void BuildC(const Mesh& mesh);
	void BuildC(const Mesh& mesh, std::map<size_t, size_t>& boundaryFid_dualVid, std::map<size_t, size_t>& vid_dualFid);
	void Build_skip_singularities(const Mesh& mesh);
	//virtual void FixOrientation(const Mesh& mesh);
	virtual void FixOrientation(Mesh& mesh);
	virtual void FixOrientation();
	virtual void BuildConnection();
	virtual void BuildConnection(const Mesh& mesh);
	virtual void BuildE(const Mesh& mesh);
	// -------------   
	virtual void BuildV_V(const Mesh& mesh);
	virtual void BuildV_E(const Mesh& mesh);
	virtual void BuildV_F(const Mesh& mesh);
	virtual void BuildV_C(const Mesh& mesh);
	// -------------
	virtual void BuildE_V(const Mesh& mesh);
	virtual void BuildE_E(const Mesh& mesh);
	virtual void BuildE_F(const Mesh& mesh);
	virtual void BuildE_C(const Mesh& mesh);
	// -------------
	virtual void BuildF_V(const Mesh& mesh);
	virtual void BuildF_E(const Mesh& mesh);
	virtual void BuildF_F(const Mesh& mesh);
	virtual void BuildF_C(const Mesh& mesh);
	// -------------
	virtual void BuildC_V(const Mesh& mesh);
	virtual void BuildC_E(const Mesh& mesh);
	virtual void BuildC_F(const Mesh& mesh);
	virtual void BuildC_C(const Mesh& mesh);
};

std::vector<size_t> get_link_fids(const Mesh& mesh, const Vertex& v);
std::vector<size_t> get_link_fids(const Mesh& mesh, const Vertex& v, std::map<size_t, size_t>& boundaryEid_dualVid);

#endif // !DUAL_MESH_H