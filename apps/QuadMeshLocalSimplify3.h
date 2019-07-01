/*
 * QuadMeshLocalSimplify.cpp
 *
 *  Created on: Nov 6, 2018
 *      Author: cotrik
 */
#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "BaseComplexQuad.h"
#include "BaseComplexSheetQuad.h"
#include "Patches.h"
#include "ArgumentManager.h"

#include <unordered_map>
#include <unordered_set>
#include <set>
#include <queue>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <math.h>

extern int maxValence/* = 5*/;
extern int minValence/* = 3*/;
extern int smoothIters/* = 20*/;
extern bool featurePreserved/* = true*/;
extern double angle/* = 160*/;
extern const double PI/* = 3.1415926535*/;
extern Mesh origMesh;
extern std::vector<size_t> userCorners;
extern std::vector<size_t> canceledCorners;
// std::vector<size_t> userCorners = { 296, 312, 441, 426 }; // fandisk
// std::vector<size_t> canceledCorners = { 310, 314, 443 };  // fandisk
extern bool COLLAPSE/* = false*/;
extern bool SPLIT/* = true*/;
extern bool conformal/* = true*/;
extern bool global/* = true*/;
extern bool ROTATE/* = true*/;
extern bool COLLAPSE_DIAGNAL/* = true*/;
extern bool REMOVE_DOUBLET/* = true*/;

std::set<size_t> get_canceledEdgeIds(const BaseComplexSheetQuad& baseComplexSheets,
	std::map<size_t, size_t>& canceledFaceIds, size_t sheetId);

bool can_collapse(const BaseComplexSheetQuad& baseComplexSheets, std::map<size_t, size_t>& canceledFaceIds, size_t sheetId);
bool can_collapse_with_feature_preserved(const BaseComplexSheetQuad& baseComplexSheets,
	std::map<size_t, size_t>& canceledFaceIds, size_t sheetId);
void collapse(Mesh& mesh, std::unordered_map<size_t, size_t>& key_edgeId, std::unordered_map<std::string, size_t>& key_faceId,
	std::map<size_t, size_t>& canceledFaceIds, std::set<size_t>& canceledEdgeIds);
void collapse_with_feature_preserved(Mesh& mesh, std::unordered_map<size_t, size_t>& key_edgeId,
	std::unordered_map<std::string, size_t>& key_faceId,
	std::map<size_t, size_t>& canceledFaceIds, std::set<size_t>& canceledEdgeIds);
bool can_collapse_vids_with_feature_preserved(Mesh& mesh, const std::vector<size_t>& vids);
bool can_collapse_vids_with_feature_preserved(Mesh& mesh, const std::set<size_t>& eids);
void global_simplify(Mesh& mesh, std::set<size_t>& canceledFids);

std::unordered_map<size_t, size_t> get_key_edgeId(const Mesh& mesh);
std::unordered_map<std::string, size_t> get_key_faceId(const Mesh& mesh);
std::string get_facekey(const Face& f);
Mesh Refine(const Mesh& hex_mesh, int clockwise);

void get_parallel_edgeids(const Mesh& mesh, size_t start_edge_id, size_t start_face_id,
	std::set<size_t>& parallel_edgeids, std::set<size_t>& parallel_faceids);
size_t get_faceid(const Mesh& mesh, size_t vid, size_t exclude_vid);
size_t get_diagnal_vid(const Mesh& mesh, size_t vid, size_t fid);
size_t get_diagnal_vid(const Mesh& mesh, size_t vid, const std::vector<size_t>& fids);

std::vector<size_t> get_collapse_vids(Mesh& mesh, size_t vid, size_t eid);
std::vector<size_t> get_split_vids(Mesh& mesh, size_t vid, size_t eid);

bool can_collapse_vids(Mesh& mesh, const std::vector<size_t>& vids, size_t target_vid);
void collapse_vids(Mesh& mesh, const std::vector<size_t>& vids, size_t target_vid);
bool can_collapse(Mesh& mesh, const std::vector<size_t>& linkVids, const std::vector<size_t>& linkEids);
void collapse(Mesh& mesh, const std::vector<size_t>& linkVids, const std::vector<size_t>& linkEids);

bool can_collapse_vids_with_feature_preserved(Mesh& mesh, const std::vector<size_t>& vids, size_t target_vid);
void collapse_vids_with_feature_preserved(Mesh& mesh, const std::vector<size_t>& vids, size_t target_vid);
std::set<size_t> get_regionVids(const Mesh& mesh, const std::vector<size_t>& linkVids);
bool can_collapse_with_feature_preserved(Mesh& mesh, const std::vector<size_t>& linkVids, const std::vector<size_t>& linkEids,
	size_t v_front_fvid, size_t v_back_fvid);
void collapse_with_feature_preserved(Mesh& mesh, const std::vector<size_t>& linkVids, const std::vector<size_t>& linkEids);
void WriteSharpEdgesVtk(const char* filename, const Mesh& m_mesh, const std::vector<FeatureLine>& featureLines);
void get_feature(Mesh& mesh);

void update(Mesh& mesh, std::set<size_t>& canceledFids);

size_t get_diagnal_vid(const Mesh& mesh, size_t vid, const std::vector<size_t>& fids, size_t end_vid);
size_t get_id(const Mesh& mesh, size_t singular_vid, size_t target_vid1, size_t target_vid2);
std::vector<size_t> get_insert_vids(Mesh& mesh, const std::vector<size_t>& linkVids, const std::vector<size_t>& linkVids1);
std::vector<size_t> get_insert_fids(Mesh& mesh, const std::vector<size_t>& linkVids1, const std::vector<size_t>& linkVids2);
bool is_neighbor(const Vertex& v, size_t vid);
bool can_split_with_feature_preserved(Mesh& mesh, size_t vid, size_t eid);
bool can_split_with_feature_preserved(Mesh& mesh, const std::vector<size_t>& linkVids, const std::vector<size_t>& linkEids,
	size_t v_front_fvid, size_t v_back_fvid);
bool split_with_feature_preserved(Mesh& mesh, const std::vector<size_t>& linkVids, const std::vector<size_t>& linkEids,
	size_t v_front_fvid, size_t v_back_fvid);

void strict_simplify(Mesh& mesh, BaseComplexQuad& baseComplex, std::set<size_t>& canceledFids);
void loose_simplify(Mesh& mesh, BaseComplexQuad& baseComplex, std::set<size_t>& canceledFids);
void simplify(Mesh& mesh, BaseComplexQuad& baseComplex, std::set<size_t>& canceledFids);
void collapse_(Mesh& mesh, BaseComplexQuad& baseComplex, std::set<size_t>& canceledFids);
void smooth_project(Mesh& origMesh, Mesh& mesh);
Mesh RefineWithFeaturePreserved(const Mesh& hex_mesh, int clockwise);
Mesh RefineWithFeaturePreserved2(const Mesh& hex_mesh, int clockwise);
std::vector<size_t> get_ids(const std::string str);
std::map<size_t, std::set<size_t>> get_patchid_fids(const Mesh& mesh);
std::map<size_t, std::set<size_t>> get_patchid_vids(const Mesh& mesh, const std::map<size_t, std::set<size_t>>& patchid_fids);
std::set<size_t> get_rotate_fids(const Mesh& mesh);
std::vector<size_t> get_neighbor_fids(const Vertex& v, const std::set<size_t>& patch_fids);
std::set<size_t> get_rotate_eids(const Mesh& mesh);
std::vector<size_t> GetLinkVids(const Mesh& mesh, const std::vector<size_t>& vids);
std::vector<size_t> GetLinkEVids(const Mesh& mesh, const std::vector<size_t>& eids);
std::vector<size_t> GetLinkVidsFromEids(const Mesh& mesh, const std::vector<size_t>& eids);
std::vector<size_t> get_boundary_eids(const Face& f0, const Face& f1, const Edge& exclude_e);
std::vector<size_t> get_rotate_vids(const Mesh& mesh, const std::vector<size_t>& boundary_eids,
	size_t diag_vid0, size_t diag_vid1);
void insert_rotate_faces(Mesh& mesh, const std::vector<size_t>& rotate_vids);
bool can_rotate1(const Mesh& mesh, const Edge& e);
bool can_rotate(const Mesh& mesh, const Edge& e);
void rotate(Mesh& mesh, const Face& f0, const Face& f1, const Edge& e, std::set<size_t>& canceledFids);
bool rotate(Mesh& mesh, const Face& f0, const Face& f1, std::set<size_t>& canceledFids);
bool rotate1(Mesh& mesh, const std::vector<size_t>& fids, std::set<size_t>& canceledFids);
bool rotate(Mesh& mesh, const std::vector<size_t>& fids, std::set<size_t>& canceledFids);
void insert_rotate_fids(Mesh& mesh, std::set<size_t>& canceledFids);
void remove_doublet(Mesh& mesh, std::set<size_t>& canceledFids);
void collapse_diagnal1(Mesh& mesh, std::set<size_t>& canceledFids);
void collapse_diagnal(Mesh& mesh, std::set<size_t>& canceledFids);
bool simplify(Mesh& mesh, int& iter);
bool simplify2(Mesh& mesh, int& iter);
double GetAngle(const Mesh& mesh, const Vertex& v, const Face& c);
std::set<size_t> get_convex_corners(const Mesh& mesh);

std::set<size_t> get_rotate_eids(const Mesh& mesh, const Vertex& v, const std::vector<size_t>& fids);
std::set<size_t> get_rotate_eids_(const Mesh& mesh);
double get_angle(const Vertex& v, const Vertex& v0, const Vertex& v1);
std::vector<size_t> get_neighbor_vids(const Mesh& mesh, const Vertex& v, size_t fid);
bool is_convex(const Mesh& mesh, const Vertex& v, const std::vector<size_t>& fids);
std::set<size_t> get_rotate_eids(const Mesh& mesh);
void rotate(Mesh& mesh, const Edge& e, const Vertex& v, std::set<size_t>& canceledFids);
void rotate(Mesh& mesh, std::set<size_t>& canceledFids);
