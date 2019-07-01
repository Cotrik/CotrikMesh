/*
* Simplifier.h
*
*  Created on: Dec 31, 2018
*      Author: cotrik
*/

#ifndef SIMPLIFIER_H
#define SIMPLIFIER_H

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "BaseComplexQuad.h"
#include "BaseComplexSheetQuad.h"
#include "Patches.h"

#include <unordered_map>
#include <unordered_set>
#include <set>
#include <queue>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <math.h>

class Simplifier{
public:
	Simplifier(Mesh& mesh);
	virtual ~Simplifier();
private:
	Simplifier();
	Simplifier(const Simplifier&);
	Simplifier& operator = (const Simplifier&);

public:
	std::set<size_t> get_canceledEdgeIds(const BaseComplexSheetQuad& baseComplexSheets,
		std::map<size_t, size_t>& canceledFaceIds, size_t sheetId);
	std::set<size_t> get_allParallelEdgeIds(const size_t eid);
	std::map<size_t, size_t> get_canceledFaceIds(const std::set<size_t>& canceledEdgeIds);
	bool can_collapse(const BaseComplexSheetQuad& baseComplexSheets, std::map<size_t, size_t>& canceledFaceIds, size_t sheetId);
	bool can_collapse_with_feature_preserved(const BaseComplexSheetQuad& baseComplexSheets,
		std::map<size_t, size_t>& canceledFaceIds, size_t sheetId);
	void collapse(std::unordered_map<size_t, size_t>& key_edgeId, std::unordered_map<std::string, size_t>& key_faceId,
		std::map<size_t, size_t>& canceledFaceIds, std::set<size_t>& canceledEdgeIds);
	void collapse_with_feature_preserved(std::unordered_map<size_t, size_t>& key_edgeId,
		std::unordered_map<std::string, size_t>& key_faceId,
		std::map<size_t, size_t>& canceledFaceIds, std::set<size_t>& canceledEdgeIds);
	bool can_collapse_vids_with_feature_preserved(const std::vector<size_t>& vids);
	bool can_collapse_vids_with_feature_preserved(const std::set<size_t>& eids);
	void global_simplify(std::set<size_t>& canceledFids);
	void global_simplify1(std::set<size_t>& canceledFids);

	bool collapse_sheet(std::set<size_t>& canceledFids, size_t parallel_edge_id);

	std::unordered_map<size_t, size_t> get_key_edgeId(const Mesh& mesh);
	std::unordered_map<std::string, size_t> get_key_faceId(const Mesh& mesh);
	std::string get_facekey(const Face& f);
	Mesh Refine(const Mesh& hex_mesh, int clockwise);

	void get_parallel_edgeids(size_t start_edge_id, size_t start_face_id,
		std::set<size_t>& parallel_edgeids, std::set<size_t>& parallel_faceids);
	size_t get_faceid(size_t vid, size_t exclude_vid);
	size_t get_diagnal_vid(size_t vid, size_t fid);
	size_t get_diagnal_vid(size_t vid, const std::vector<size_t>& fids);

	std::vector<size_t> get_collapse_vids(size_t vid, size_t eid);
	std::vector<size_t> get_split_vids(size_t vid, size_t eid);

	bool can_collapse_vids(const std::vector<size_t>& vids, size_t target_vid);
	void collapse_vids(const std::vector<size_t>& vids, size_t target_vid);
	bool can_collapse(const std::vector<size_t>& linkVids, const std::vector<size_t>& linkEids);
	void collapse(const std::vector<size_t>& linkVids, const std::vector<size_t>& linkEids);

	bool can_collapse_vids_with_feature_preserved(const std::vector<size_t>& vids, size_t target_vid);
	void collapse_vids_with_feature_preserved(std::vector<size_t>& vids, size_t target_vid);
	std::set<size_t> get_regionVids(const std::vector<size_t>& linkVids);
	bool can_collapse_with_feature_preserved(const std::vector<size_t>& linkVids, const std::vector<size_t>& linkEids);
	bool can_collapse_with_feature_preserved(const std::vector<size_t>& linkVids, const std::vector<size_t>& linkEids,
		size_t v_front_fvid, size_t v_back_fvid);
    bool can_collapse_with_feature_preserved(const std::vector<size_t>& linkVids, const std::vector<size_t>& linkEids,
        size_t v_front_fvid);
	void collapse_with_feature_preserved(const std::vector<size_t>& linkVids, const std::vector<size_t>& linkEids);
	void WriteSharpEdgesVtk(const char* filename, const Mesh& m_mesh, const std::vector<FeatureLine>& featureLines);
	void get_feature();

	void update(std::set<size_t>& canceledFids);
	void update(const BaseComplexQuad& baseComplexQuad);
	void init();
	void align_feature();

	size_t get_diagnal_vid(size_t vid, const std::vector<size_t>& fids, size_t end_vid);
	size_t get_id(size_t singular_vid, size_t target_vid1, size_t target_vid2);
	std::vector<size_t> get_insert_vids(const std::vector<size_t>& linkVids, const std::vector<size_t>& linkVids1);
	std::vector<size_t> get_insert_fids(const std::vector<size_t>& linkVids1, const std::vector<size_t>& linkVids2);
	bool is_neighbor(const Vertex& v, size_t vid);
	bool can_split_with_feature_preserved(size_t vid, size_t eid);
	bool can_split_with_feature_preserved(const std::vector<size_t>& linkVids, const std::vector<size_t>& linkEids,
		size_t v_front_fvid, size_t v_back_fvid);
	bool split_with_feature_preserved(const std::vector<size_t>& linkVids, const std::vector<size_t>& linkEids,
		size_t v_front_fvid, size_t v_back_fvid);

	void strict_simplify(BaseComplexQuad& baseComplex, std::set<size_t>& canceledFids);
	void strict_simplify_reverse(BaseComplexQuad& baseComplex, std::set<size_t>& canceledFids);
	void loose_simplify(BaseComplexQuad& baseComplex, std::set<size_t>& canceledFids);
	void loose_simplify_reverse(BaseComplexQuad& baseComplex, std::set<size_t>& canceledFids);
	void loose_simplify_random(BaseComplexQuad& baseComplex, std::set<size_t>& canceledFids);
	int get_id(std::unordered_set<int>& ids, int n);
	void half_simplify(BaseComplexQuad& baseComplex, std::set<size_t>& canceledFids);
	void sheet_simplify(BaseComplexQuad& baseComplex, std::set<size_t>& canceledFids);
	void simplify(BaseComplexQuad& baseComplex, std::set<size_t>& canceledFids);
	void smooth_project();
	void smooth_project(int resolution);
    void smooth_project1(int resolution);
	double laplacian_positive_cotan_weight(const Vertex& vi, const Edge& e);

	Mesh RefineWithFeaturePreserved(const Mesh& hex_mesh, int clockwise);
	Mesh RefineWithFeaturePreserved2(const Mesh& hex_mesh, int clockwise);
	std::vector<size_t> get_ids(const std::string str);
	std::map<size_t, std::set<size_t>> get_patchid_fids();
	std::map<size_t, std::set<size_t>> get_patchid_vids(const std::map<size_t, std::set<size_t>>& patchid_fids);
	std::set<size_t> get_rotate_fids();
	std::vector<size_t> get_neighbor_fids(const Vertex& v, const std::set<size_t>& patch_fids);
	std::set<size_t> get_rotate_eids();
	std::vector<size_t> GetLinkVids(const std::vector<size_t>& vids);
	std::vector<size_t> GetLinkEVids(const std::vector<size_t>& eids);
	std::vector<size_t> GetLinkVidsFromEids(const std::vector<size_t>& eids);
	std::vector<size_t> get_boundary_eids(const Face& f0, const Face& f1, const Edge& exclude_e);
	std::vector<size_t> get_rotate_vids(const std::vector<size_t>& boundary_eids,
		size_t diag_vid0, size_t diag_vid1);
	void insert_rotate_faces(const std::vector<size_t>& rotate_vids);
	bool can_rotate1(const Edge& e);
	bool can_rotate(const Edge& e);
	void rotate(const Face& f0, const Face& f1, const Edge& e, std::set<size_t>& canceledFids);
	bool rotate(const Face& f0, const Face& f1, std::set<size_t>& canceledFids);
	bool rotate1(const std::vector<size_t>& fids, std::set<size_t>& canceledFids);
	bool rotate(const std::vector<size_t>& fids, std::set<size_t>& canceledFids);
	void insert_rotate_fids(std::set<size_t>& canceledFids);
	void remove_doublet(std::set<size_t>& canceledFids);
	void collapse_diagnal1(std::set<size_t>& canceledFids);
	void collapse_diagnal(std::set<size_t>& canceledFids);
	virtual bool simplify(int& iter);
	double GetAngle(const Vertex& v, const Face& c);
	std::set<size_t> get_convex_corners();

	std::set<size_t> get_rotate_eids(const Vertex& v, const std::vector<size_t>& fids);
	std::set<size_t> get_rotate_eids_();
	double get_angle(const Vertex& v, const Vertex& v0, const Vertex& v1);
	std::vector<size_t> get_neighbor_vids(const Vertex& v, size_t fid);
	double get_angle(const Vertex& v, const std::vector<size_t>& fids);
	bool is_convex(const Vertex& v, const std::vector<size_t>& fids);
	bool is_concave(const Vertex& v, const std::vector<size_t>& fids);
	size_t get_ideal_valence(const Vertex& v, const std::vector<size_t>& fids);
	void rotate(const Edge& e, const Vertex& v, std::set<size_t>& canceledFids);
	void rotate(std::set<size_t>& canceledFids);

	bool hasSingularities() const;

	void Collapse(size_t vid, size_t target_vid);
public:
	Mesh & mesh;
	Mesh origMesh;
	static std::vector<size_t> userCorners;
	static std::vector<size_t> canceledCorners;
	static double angle/* = 1605*/;
	static int maxValence/* = 5*/;
	static int minValence/* = 3*/;
	static int iters/* = 10000*/;
	static int smoothIters/* = 20*/;
	static int resolution/* = 3*/;
	static bool featurePreserved/* = true*/;
	static bool COLLAPSE/* = true*/;
	static bool SPLIT/* = false*/;
	static bool CONFORMAL/* = true*/;
	static bool GLOBAL/* = true*/;
	static bool ROTATE/* = true*/;
	static bool COLLAPSE_DIAGNAL/* = true*/;
	static bool REMOVE_DOUBLET/* = true*/;
	static bool SHEET_SPLIT/* = true*/;
	static bool HALF/* = false*/;
	static bool TRIP/* = true*/;
	static bool checkCorner/* = true*/;
	static bool writeFile/* = true*/;
};

const extern double PI/* = 3.1415926535*/;
//extern int maxValence/* = 5*/;
//extern int minValence/* = 35*/;
//extern int iters/* = 100005*/;
//extern int smoothIters/* = 205*/;
//extern bool featurePreserved/* = true5*/;
//extern double angle/* = 1605*/;
//extern bool COLLAPSE/* = true5*/;
//extern bool SPLIT/* = false5*/;
//extern bool conformal/* = true5*/;
//extern bool global/* = true5*/;
//extern bool ROTATE/* = true5*/;
//extern bool COLLAPSE_DIAGNAL/* = true5*/;
//extern bool REMOVE_DOUBLET/* = true5*/;

#endif // !SIMPLIFIER_H
