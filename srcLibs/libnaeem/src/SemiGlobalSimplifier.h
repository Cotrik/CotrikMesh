/*
* SemiGlobalSimplifier.h
*
*  Created on: October 25, 2021
*      Author: https://github.com/naeem014
*/

#ifndef SEMI_GLOBAL_SIMPLIFIER_H_
#define SEMI_GLOBAL_SIMPLIFIER_H_

#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <queue>
#include <memory>
#include <mutex>
#include <map>
#include <functional>
#include <stdlib.h>
#include <time.h>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include "Mesh.h"
#include "Memento.h"
#include "BaseComplexQuad.h"
#include "ChordExtractor.h"
#include "ChordCollapse.h"
// #include "SimplificationOperation.h"
#include "DiagonalCollapse.h"
#include "DiagonalThreeFivePair.h"
#include "DirectSeparatrixCollapse.h"
#include "SeparatrixCollapse.h"
#include "EdgeRotation.h"
#include "VertexRotation.h"
#include "EdgeCollapse.h"
#include "VertexSplit.h"
#include "QuadSplit.h"
#include "ThreeFivePair.h"
#include "MeshUtil.h"
#include "SingularityPair.h"
#include "Smooth.h"
#include "PQueue.h"
#include "PQueue.cpp"
#include "Renderer.h"
#include "KDTree.h"

const double PI = 3.1415926535;

struct SingularityLink {
    std::vector<size_t> linkVids;
    std::vector<size_t> linkEids;
    int frontId;
    int backId;
    int id;
    int a = 0;
    int b = 0;
    double rank = 1.0;
    int rots = 0;
    int volt = 0;
    bool diagonal = false;
    double* delta;
};

struct Separatrix {
    size_t threeId;
    size_t fiveId;
    size_t frontId;
    size_t backId;
    int b1 = 0;
    int b2 = 0;
    std::vector<size_t> vidsA;
    std::vector<size_t> vidsB;
    bool empty = true;
    int minSize = -1;
    int rots = 0;
    std::map<size_t, double> dirMap;
    bool diagonal = false;
    bool crossb1 = false;
    bool crossb2 = false;
    size_t moveDir;
    double score = 0;
    /*std::vector<int> frontdirMap;
    std::vector<int> backdirMap;
    std::vector<int> score;
    std::vector<int> frontScore;
    std::vector<int> backScore;*/
};

struct Path {
    std::vector<size_t> vids;
    std::unordered_set<size_t> vset;
    double singularity_score = 0;
    double qem_score = 0;
    double score = 0;
    int b1 = 0;
    int b2 = 0;
    int rots = 0;
    bool is_diagonal = false;
    bool intersect_b1 = false;
    bool intersect_b2 = false;
    int boundary_distance = 0;
    size_t move_dir;
};

struct p_comp {
    public: 
    
    bool operator() (Path& p1, Path& p2) {
        return p1.score < p2.score;
    }
};


struct SeparatrixComparator {
    public:

    bool operator()(Separatrix& l, Separatrix& r) {
        return l.score < r.score;
    }
};

struct LinkComparator {
    public:

    bool operator()(SingularityLink& l, SingularityLink& r) {
        // if (l.rank == r.rank) {
        //     return (l.a+l.b) > (r.a+r.b);
        // }
        // return fabs(*l.delta + l.rank) > fabs(*l.delta + r.rank);
        return fabs(l.rank) > fabs(r.rank);
        // return (l.a+l.b) >= (r.a+r.b) && l.rank >= r.rank;
    }
};

struct SingularityGroup {
    SingularityLink l1;
    SingularityLink l2;
    double rank = 0.0;
    int l1l2Rots = 0;
    int l2l1Rots = 0;
    int l1l2Volt = 0;
    int l2l1Volt = 0;
};

struct GroupComparator {
    public:

    bool operator()(SingularityGroup& l, SingularityGroup& r) {
        return l.rank > r.rank;
    }
};

struct vQEM {
    double qem = 0;
    double lower_bound = 0;
    double upper_bound = 0;
};

struct vMesh {
    // std::unordered_map<size_t, Vertex> vmap;
    // std::unordered_map<size_t, Edge> emap;
    // std::unordered_map<size_t, Face> fmap;
    std::map<size_t, Vertex> vmap;
    std::map<size_t, Edge> emap;
    std::map<size_t, Face> fmap;
    
    int max_vid, max_fid, maxvid, maxfid;
    Mesh* mesh;

    vMesh() {}

    vMesh(Mesh* mesh_) {
        mesh = mesh_;
        max_vid = mesh->V.size();
        max_fid = mesh->F.size();
        maxvid = mesh->V.size()-1;
        maxfid = mesh->F.size()-1;
    }
    vMesh& operator=(const vMesh& other) {
        if (this == &other) return *this;
        vmap = other.vmap;
        fmap = other.fmap;
        mesh = other.mesh;
        max_vid = other.max_vid;
        max_fid = other.max_fid;
        maxvid = other.maxvid;
        maxfid = other.maxfid;
        return *this;
    }

    void SetMesh(Mesh* mesh_) {
        mesh = mesh_;
        max_vid = mesh->V.size();
        max_fid = mesh->F.size();
        maxvid = mesh->V.size()-1;
        maxfid = mesh->F.size()-1;
    }

    Vertex& AddVertex(glm::dvec3 coords, size_t id = 0) {
        max_vid += 4;
        // if (id == 0) {
        //     id = max_vid;
        // }
        // size_t hash = std::hash<size_t>{}(std::chrono::high_resolution_clock::now().time_since_epoch().count() + max_vid);
        vmap[max_vid] = Vertex(coords);
        vmap[max_vid].id = max_vid;
        return vmap[max_vid];
    }

    Face& AddFace(std::vector<size_t> vids, double threshold_shape = 0.0, double threshold_size = 0.0, double avg_area = 0.0) {
        max_fid += 4;
        // size_t hash = 0;
        // for (auto vid: vids) {
        //     hash ^= std::hash<size_t>{}(vid);
        // }
        // std::cout << "hash: " << hash << std::endl;
        fmap[max_fid] = Face(vids);
        fmap[max_fid].id = max_fid;
        fmap[max_fid].threshold_shape = threshold_shape;
        fmap[max_fid].threshold_size = threshold_size;
        fmap[max_fid].avg_area = avg_area;
        return fmap[max_fid];
    }

    void setVertex(size_t vid) {
        if (vmap.find(vid) != vmap.end()) return;
        auto& v = mesh->V.at(vid);
        vmap[vid] = Vertex(v.xyz());
        vmap[vid].id = v.id; vmap[vid].isBoundary = v.isBoundary; vmap[vid].type = v.type;
        vmap[vid].N_Fids = v.N_Fids;
    }

    void setFace(size_t fid) {
        if (fmap.find(fid) != fmap.end()) return;
        auto& f = mesh->F.at(fid);
        fmap[fid] = Face(f.Vids);
        fmap[fid].id = f.id;
        fmap[fid].threshold_shape = f.threshold_shape;
        fmap[fid].threshold_size = f.threshold_size;
        fmap[fid].avg_area = f.avg_area;
    }

    Vertex& getVertex(size_t vid, bool useVM = true) {
        if (useVM && vmap.find(vid) != vmap.end()) return vmap[vid];
        return mesh->V.at(vid);
    }

    Edge& getEdge(size_t eid) {
        if (emap.find(eid) != emap.end()) return emap[eid];
        return mesh->E.at(eid);
    }

    Face& getFace(size_t fid, bool useVM = true) {
        // std::cout << "getting face with fid: " << fid << std::endl;
        if (useVM && fmap.find(fid) != fmap.end()) return fmap[fid];
        return mesh->F.at(fid);
    }

    double getAngle(size_t vid1, size_t vid2, size_t vid3, bool useVM = true) {
        auto& v1 = getVertex(vid1, useVM);
        auto& v2 = getVertex(vid2, useVM);
        auto& v3 = getVertex(vid3, useVM);
        
        auto A = v2.xyz() - v1.xyz();
        auto B = v3.xyz() - v1.xyz();
        double angle = atan2(glm::length(glm::cross(A, B)), glm::dot(A, B));
        if (angle < 0) angle += 2*PI;
        return angle * 180.0 / PI; 
    }
    
    int getIdealValence(size_t vid, bool useVM = true, bool log = false) {
        auto& v = getVertex(vid, useVM);
        size_t start = v.N_Fids.at(0);
        for (int i = 0; i < v.N_Fids.size(); i++) {
            auto& f = getFace(v.N_Fids.at(i), useVM);
            int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
            auto& fv = getVertex(f.Vids.at((idx+1)%4), useVM);
            if (fv.isBoundary || fv.type == FEATURE) {
                start = v.N_Fids.at(i);
                break;
            }
        }
        std::vector<size_t> N_Fids;
        for (int i = 0; i < v.N_Fids.size(); i++) {
            auto& f = getFace(start, useVM);
            int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
            N_Fids.push_back(f.id);
            for (auto fid: v.N_Fids) {
                auto& tf = getFace(fid, useVM);              
                int idx2 = std::distance(tf.Vids.begin(), std::find(tf.Vids.begin(), tf.Vids.end(), v.id));
                if (f.Vids.at((idx+3)%4) == tf.Vids.at((idx2+1)%4)) {
                    start = fid;
                    break;
                }
            }
        }
        int idealValence = 0;
        double total_angle = 0.0;
        // std::cout << "calculating ideal valence for: " << v.id << "(" << v.N_Fids.size() << ")" << " N_Fids: " << N_Fids.size() << std::endl;
        for (auto fid: N_Fids) {
            auto& f = getFace(fid, useVM);
            // std::cout << "f(" << f.id << "): ";
            // for (auto fvid: f.Vids) std::cout << fvid << " "; std::cout << std::endl;
            int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
            // auto& v1 = getVertex(f.Vids.at((idx+1)%4), useVM);
            // auto& v2 = getVertex(f.Vids.at((idx+3)%4), useVM);
            // total_angle += getAngle(v.id, v1.id, v2.id, useVM);
            total_angle += getAngle(v.id, f.Vids.at((idx+1)%4), f.Vids.at((idx+2)%4), useVM);
            total_angle += getAngle(v.id, f.Vids.at((idx+2)%4), f.Vids.at((idx+3)%4), useVM);
            if (getVertex(f.Vids.at((idx+3)%4), useVM).isBoundary || getVertex(f.Vids.at((idx+3)%4), useVM).type == FEATURE) {
                idealValence += (floor(total_angle / 90.0) + floor(std::fmod(total_angle, 90.0) / 45));
                // std::cout << "total_angle: " << total_angle << " idealValence: " << idealValence << std::endl;
                total_angle = 0.0;
            }
        }
        // std::cout << "************************************" << std::endl;
        // if (log) std::cout << "v nfids: " << v.N_Fids.size() << std::endl;
        // v.idealValence = 0;
        /*int idealValence = 0;
        double total_angle = 0.0;
        std::vector<size_t> features;
        // size_t next = v.N_Fids.at(0);
        for (int i = 0; i < v.N_Fids.size(); i++) {
            auto& f = getFace(v.N_Fids.at(i), useVM);
            int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), vid));
            if (log) {
                // std::cout << "useVM: " << useVM << std::endl;
                // std::cout << "vid: " << v.id << " idx: " << idx << " f vids: " << f.Vids.size() << std::endl;
                // std::cout << "f Vids: "; 
                // for (auto fvid: f.Vids) std::cout << fvid << " "; 
                // std::cout << std::endl;
                // std::cout << "mesh->V: " << mesh->V.size() << " vmap: " << vmap.size() << std::endl;
                // for (auto it = vmap.begin(); it != vmap.end(); it++) {
                    // std::cout << it->first << " " << it->second.id << std::endl;
                // }
            }
            auto& b1 = getVertex(f.Vids.at((idx+1)%4), useVM);
            // if (log) {
            //     std::cout << "b1: " << b1.id << std::endl;
            // }
            auto& b2 = getVertex(f.Vids.at((idx+3)%4), useVM);
            // if (log) {
            //     std::cout << "b2: " << b2.id << std::endl;
            // }
            if (b1.isBoundary || b1.type == FEATURE) features.push_back(f.Vids.at((idx+1)%4));
            if (b2.isBoundary || b2.type == FEATURE) features.push_back(f.Vids.at((idx+3)%4));
            // for (auto fvid: f.Vids) {
            //     auto& fv = getVertex(fvid, useVM);
            //     if (fv.id != vid && fv.isBoundary || fv.type == FEATURE) features.push_back(fvid);
            // }
            // auto& fv = getVertex(f.Vids.at((idx+1)%4), useVM);
            // if (fv.isBoundary || fv.type == FEATURE) features.push_back(fv.id);
            // for (auto fid: v.N_Fids) {
            //     auto& tf = getFace(fid, useVM);
            //     int _idx = std::distance(tf.Vids.begin(), std::find(tf.Vids.begin(), tf.Vids.end(), vid));
            //     if (f.Vids.at((idx+3)%4) == tf.Vids.at((_idx+1)%4)) {
            //         next = fid;
            //         break;
            //     }
            // }
        }
        std::sort(features.begin(), features.end());
        auto it = std::unique(features.begin(), features.end());
        features.resize(std::distance(features.begin(), it));
        // if (log) {
        //     std::cout << "features: " << features.size() << std::endl;
        //     for (auto id: features) std::cout << id << " "; std::cout << std::endl;
        // }
        if (features.size() == 2) {
            double angle = getAngle(vid, features.at(0), features.at(1), useVM);
            // if (log) std::cout << "angle: " << angle << std::endl;
            idealValence += (floor(angle / 90.0) + floor(std::fmod(angle, 90.0) / 45));
            // if (log) std::cout << "ideal valence: " << idealValence << std::endl;
        } else {
            for (int i = 0; i < features.size(); i++) {
                double angle = getAngle(vid, features.at(i), features.at((i+1)%features.size()), useVM);
                // if (log) std::cout << "angle: " << angle << std::endl;
                idealValence += (floor(angle / 90.0) + floor(std::fmod(angle, 90.0) / 45));
                // if (log) std::cout << "ideal valence: " << idealValence << std::endl;

            }
        }
        // for (int i = 0; i < v.N_Fids.size(); i++) {
        //     auto& f = getFace(v.N_Fids.at(i), useVM);
        //     int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), vid));
        //     auto& fv = getVertex(f.Vids.at((idx+1)%4), useVM);
        //     if (fv.isBoundary || fv.type == FEATURE) {
        //         // double total_angle = 0.0;
        //         features.push_back(fv.id);
        //         // size_t next = v.N_Fids.at(i);
        //         // for (int j = 0; j < v.N_Fids.size(); j++) {
        //             // std::cout << "i " << i << " j " << j  << " i+v.N_Fids.size() " << i+v.N_Fids.size() << std::endl;
        //             // auto& _f = getFace(next, useVM);
        //             // int _idx = std::distance(_f.Vids.begin(), std::find(_f.Vids.begin(), _f.Vids.end(), vid));
        //             // total_angle += getAngle(v.id, _f.Vids.at((_idx+1)%_f.Vids.size()), _f.Vids.at((_idx+3)%_f.Vids.size()));
        //             // auto& _fv = getVertex(_f.Vids.at((_idx+3)%_f.Vids.size()), useVM);
        //             for (auto fid: v.N_Fids) {
        //                 auto& tf = getFace(fid, useVM);
        //                 int _idx2 = std::distance(tf.Vids.begin(), std::find(tf.Vids.begin(), tf.Vids.end(), vid));
        //                 if (_f.Vids.at((idx+3)%_f.Vids.size()) == tf.Vids.at((_idx2+1)%tf.Vids.size())) {
        //                     next = fid;
        //                     break;
        //                 }
        //             }
        //             if (_fv.isBoundary || _fv.type == FEATURE) {
        //                 break;
        //             }
        //         }
        //         // std::cout << "total angle: " << total_angle << std::endl;
        //         idealValence += (floor(total_angle / 90.0) + floor(std::fmod(total_angle, 90.0) / 45));
        //     }
        // }
        // if (log) std::cout << "boundary vertex: " << idealValence << " " << v.N_Fids.size() << std::endl;*/
        return idealValence;
    }

    void Update() {
        for (auto it = fmap.begin(); it != fmap.end(); it++) {
            auto& f = it->second;
            if (f.id > maxfid) {
                int nId = mesh->F.size();
                // std::cout << "face id: " << it->first << " new id: " << nId << std::endl;
                for (int i = 0; i < f.Vids.size(); i++) {
                    auto& v = getVertex(f.Vids.at(i));
                    // std::cout << "new Vertex N_Fids before: "; for (auto nfid: v.N_Fids) std::cout << nfid << " "; std::cout << std::endl;
                    int idx = std::distance(v.N_Fids.begin(), std::find(v.N_Fids.begin(), v.N_Fids.end(), it->first));
                    v.N_Fids.at(idx) = nId;
                    // std::cout << "new Vertex N_Fids after: "; for (auto nfid: v.N_Fids) std::cout << nfid << " "; std::cout << std::endl;
                }
                f.id = nId;
                mesh->F.push_back(f);
                // mesh->F.resize(mesh->F.size()+1);
                // mesh->F.at(mesh->F.size()-1) = f;
            } else {
                // std::cout << "vmap f: " << it->first << " " << f.id << std::endl;
                // for (auto fvid: f.Vids) std::cout << fvid << " "; std::cout << std::endl;
                mesh->F.at(it->first).Vids = f.Vids;
                // std::cout << "mesh f:" << std::endl; 
                // for (auto fvid: mesh->F.at(it->first).Vids) std::cout << fvid << " "; std::cout << std::endl;

            }
        }
        for (auto it = vmap.begin(); it != vmap.end(); it++) {
            if (it->second.id > maxvid) {
                auto& v = it->second;
                for (int i = 0; i < v.N_Fids.size(); i++) {
                    auto& f = mesh->F.at(v.N_Fids.at(i));
                    // std::cout << "Face id: " << f.id << " it->first: " << it->first << " v.id: " << v.id << std::endl;
                    // std::cout << "face Vids before: "; for (auto fvid: f.Vids) std::cout << fvid << " "; std::cout << std::endl;
                    int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
                    f.Vids.at(idx) = mesh->V.size();
                    // int idx2 = (idx+1)%4;
                    // if (idx2 > idx && f.Vids.at(idx2) == v.id) {
                    //     int rot = 4-idx2;
                    //     std::rotate(f.Vids.begin(), f.Vids.begin()+f.Vids.size()-rot, f.Vids.end());
                    // }
                    // std::cout << "face Vids after: "; for (auto fvid: f.Vids) std::cout << fvid << " "; std::cout << std::endl;
                }                
                v.id = mesh->V.size();
                mesh->V.push_back(v);
                // mesh->V.resize(mesh->V.size()+1);
                // mesh->V.at(mesh->V.size()-1) = v;
            } else {
                // std::cout << "Inside Update" << std::endl;
                // std::cout << "v id: " << it->first << "(" << it->second.N_Fids.size() << ")" << std::endl; 
                // std::cout << "regular: " << (it->second.type == REGULAR) <<  " feature: " << (it->second.type == FEATURE) << " isBoundary: " << it->second.isBoundary << std::endl;
                mesh->V.at(it->first).xyz(it->second.xyz());
                mesh->V.at(it->first).N_Fids = it->second.N_Fids;
                mesh->V.at(it->first).isBoundary = it->second.isBoundary;
                mesh->V.at(it->first).type = it->second.type;
            }
        }
    }

    std::vector<std::vector<size_t>> planes(const Vertex& v, bool useVM = true) {
        std::vector<std::vector<size_t>> p;
        if (v.N_Fids.empty()) return p;
        size_t start = v.N_Fids.at(0);
        // std::cout << "start: " << start << std::endl;
        int countFeatures = [&]() {
            int count = 0;
            for (auto fid: v.N_Fids) {
                // std::cout << "useVM: " << useVM << std::endl;
                auto& f = getFace(fid, useVM);
                int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
                // std::cout << "idx: " << idx << std::endl;
                // std::cout << "f Vids: " << f.Vids.size() << std::endl;
                if (getVertex(f.Vids[(idx + 1) % f.Vids.size()], useVM).isBoundary || getVertex(f.Vids[(idx + 1) % f.Vids.size()], useVM).type == FEATURE) {
                    count++;
                    start = fid;
                }
            }
            return count;
        }();
        // std::cout << "countFeatures: " << countFeatures << std::endl;
        if (countFeatures < 2) {
            p.push_back(v.N_Fids);
        } else {
            p.resize(countFeatures);
            // std::cout << "p size: " << p.size() << std::endl;
            int k = 0;
            for (int i = 0; i < v.N_Fids.size(); i++) {
                auto& f = getFace(start, useVM);
                // std::cout << "setting face at k: " << k << std::endl;
                p[k].push_back(f.id);
                int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
                auto& next = getVertex(f.Vids[(idx + 3) % f.Vids.size()], useVM);
                if (next.isBoundary || next.type == FEATURE) {
                    k++;
                    // std::cout << "k: " << k << std::endl;
                }
                for (auto fid: v.N_Fids) {
                    if (fid == start) continue;
                    auto& tf = getFace(fid, useVM);
                    int tidx = std::distance(tf.Vids.begin(), std::find(tf.Vids.begin(), tf.Vids.end(), v.id));
                    if (tf.Vids[(tidx + 1) % tf.Vids.size()] == next.id) {
                        start = fid;
                        break;
                    }
                }
            }
        }
        // std::cout << "planes " << p.size() << std::endl;
        return p;
    }

    int idealValence(const Vertex& v, bool useVM = true, std::vector<size_t> qids = {}) {
        const double PI = 3.1415926535;
        auto angle = [&] (Vertex& prev, Vertex& next) {
            glm::dvec3 p = prev - v;
            glm::dvec3 n = next - v;
            return std::acos(glm::dot(p, n) / (glm::length(p) * glm::length(n)));
        };
        if (v.N_Fids.empty()) return 0;
        /*size_t start = v.N_Fids.at(0);
        int countFeatures = [&]() {
            int count = 0;
            for (auto fid: v.N_Fids) {
                auto& f = getFace(fid, useVM);
                int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
                if (getVertex(f.Vids[(idx + 1) % f.Vids.size()], useVM).isBoundary || getVertex(f.Vids[(idx + 1) % f.Vids.size()], useVM).type == FEATURE) {
                    count++;
                    start = fid;
                }
            }
            return count;
        }();
        if (countFeatures < 2) {
            double sum = 0.0;
            for (auto fid: v.N_Fids) {
                auto& f = getFace(fid, useVM);
                int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
                sum += angle(getVertex(f.Vids[(idx + 1) % f.Vids.size()], useVM), getVertex(f.Vids[(idx + 2) % f.Vids.size()], useVM));
                sum += angle(getVertex(f.Vids[(idx + 2) % f.Vids.size()], useVM), getVertex(f.Vids[(idx + 3) % f.Vids.size()], useVM));
            }
            int turns = (int) std::round(sum / (PI / 2.0));
            return std::max(turns, 1);
        } else {
            double sum = 0.0;
            std::vector<std::pair<std::vector<size_t>, int>> ideals(countFeatures);
            int k = 0;
            for (int i = 0; i < v.N_Fids.size(); i++) {
                auto& f = getFace(start, useVM);
                ideals[k].first.push_back(f.id);
                int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
                auto& prev = getVertex(f.Vids[(idx + 1) % f.Vids.size()], useVM);
                auto& diag = getVertex(f.Vids[(idx + 2) % f.Vids.size()], useVM);
                auto& next = getVertex(f.Vids[(idx + 3) % f.Vids.size()], useVM);
                sum += angle(prev, diag);
                sum += angle(diag, next);
                if (next.isBoundary || next.type == FEATURE) {
                    int turns = (int) std::round(sum / (PI / 2.0));
                    ideals[k].second = std::max(turns, 1);
                    sum = 0.0;
                    k++;
                }
                for (auto fid: v.N_Fids) {
                    if (fid == start) continue;
                    auto& tf = getFace(fid, useVM);
                    int tidx = std::distance(tf.Vids.begin(), std::find(tf.Vids.begin(), tf.Vids.end(), v.id));
                    if (tf.Vids[(tidx + 1) % tf.Vids.size()] == next.id) {
                        start = fid;
                        break;
                    }
                }
            }*/
            auto p = planes(v, useVM);
            // std::cout << "p size: " << p.size() << std::endl;
            int total = 0;
            for (auto plane: p) {
                double sum = 0.0;
                for (auto id: plane) {
                    auto& f = getFace(id, useVM);
                    int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
                    auto& prev = getVertex(f.Vids[(idx + 1) % f.Vids.size()], useVM);
                    auto& diag = getVertex(f.Vids[(idx + 2) % f.Vids.size()], useVM);
                    auto& next = getVertex(f.Vids[(idx + 3) % f.Vids.size()], useVM);
                    sum += angle(prev, diag);
                    sum += angle(diag, next);
                }
                int turns = (int) std::round(sum / (PI / 2.0));
                int ideal = std::max(turns, 1);
                total += ideal;
                for (auto qid: qids) {
                    if (std::find(plane.begin(), plane.end(), qid) != plane.end()) {
                        return ideal;        
                    }
                }
            }
            // int total = 0;
            // for (auto ideal: ideals) {
            //     for (auto qid: qids) {
            //         if (std::find(ideal.first.begin(), ideal.first.end(), qid) != ideal.first.end()) {
            //             return ideal.second;
            //         }
            //     }
            //     total += ideal.second;
            // }
            return total;
        // }
    }

    int valence(const Vertex& v, bool useVM = true, std::vector<size_t> qids = {}) {
        if (v.N_Fids.empty()) return 0;
        
        // std::cout << "getting planes" << std::endl;
        auto p = planes(v, useVM);
        // std::cout << "planes size: " << p.size() << std::endl;
        int total = 0;
        for (auto plane: p) {
            for (auto qid: qids) {
                if (std::find(plane.begin(), plane.end(), qid) != plane.end()) {
                    return plane.size();
                }
            }
            total += plane.size();
        }
        // std::cout << "total: " << total << std::endl;
        return total;
    }

    int virtualValence(const Vertex& v, bool useVM = true, std::vector<size_t> qids = {}) {
        // if (v.isBoundary || v.type == FEATURE) {
            return 4 + valence(v, useVM, qids) - idealValence(v, useVM, qids);
        // } else {
            // return valence(v);
        // }
    }

    bool corner(const Vertex& v, bool useVM = true) {
        if (!v.isBoundary && v.type != FEATURE) return false;
        return (v.isBoundary && idealValence(v, useVM) != 2) || (v.type == FEATURE && idealValence(v, useVM) != 4);
    }

    void Clear() {
        vmap.clear();
        emap.clear();
        fmap.clear();
    }
};

struct vInfo {
    std::vector<size_t> N_Vids;
    std::vector<size_t> N_Fids;
    vInfo(Mesh* mesh, size_t vid, vMesh* m = nullptr, bool useVM = true) {
        const auto& getVertex = [&] (size_t vid) {
            if (m == nullptr) return mesh->V.at(vid);
            return m->getVertex(vid, useVM);
        };
        const auto& getFace = [&] (size_t fid) {
            if (m == nullptr) return mesh->F.at(fid);
            return m->getFace(fid, useVM);
        };
        const auto& v = getVertex(vid);
        size_t start = v.N_Fids.at(0);
        if (v.isBoundary) {
            for (int i = 0; i < v.N_Fids.size(); i++) {
                const auto& f = getFace(v.N_Fids.at(i));
                int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
                const auto& fv = getVertex(f.Vids.at((idx+1)%4));
                if (fv.isBoundary) {
                    start = v.N_Fids.at(i);
                    break;
                }
            }
        }
        
        for (int i = 0; i < v.N_Fids.size(); i++) {
            const auto& f = getFace(start);
            int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
            N_Vids.push_back(f.Vids.at((idx+1)%4));
            N_Fids.push_back(f.id);
            for (auto fid: v.N_Fids) {
                const auto& tf = getFace(fid);              
                int idx2 = std::distance(tf.Vids.begin(), std::find(tf.Vids.begin(), tf.Vids.end(), v.id));
                if (f.Vids.at((idx+3)%4) == tf.Vids.at((idx2+1)%4)) {
                    start = fid;
                    break;
                }
            }
        }
        if (v.isBoundary) {
            const auto& f = getFace(start);
            int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
            N_Vids.push_back(f.Vids.at((idx+3)%4));
        }
    }

    std::vector<size_t> vids(int startId = -1) {
        if (startId == -1 || N_Vids.at(0) == startId) return N_Vids;
        int idx = std::distance(N_Vids.begin(), std::find(N_Vids.begin(), N_Vids.end(), startId));
        std::rotate(N_Vids.begin(), N_Vids.begin()+idx, N_Vids.end());
        std::rotate(N_Fids.begin(), N_Fids.begin()+idx, N_Fids.end());
        return N_Vids;
    }

    std::vector<size_t> fids(int startId = -1) {
        if (startId == -1 || N_Fids.at(0) == startId) return N_Fids;
        int idx = std::distance(N_Fids.begin(), std::find(N_Fids.begin(), N_Fids.end(), startId));
        std::rotate(N_Vids.begin(), N_Vids.begin()+idx, N_Vids.end());
        std::rotate(N_Fids.begin(), N_Fids.begin()+idx, N_Fids.end());
        return N_Fids;
    }

    int nvids() {
        return N_Vids.size();
    }

    int nfids() {
        return N_Fids.size();
    }

};


struct tfPair {
    size_t tId, fId;
    bool diag = false;
    tfPair(size_t tId_, size_t fId_, bool diag_) {
        tId = tId_;
        fId = fId_;
        diag = diag_;
    }

    void setIds(size_t tId_, size_t fId_) {
        tId = tId_;
        fId = fId_;
    }

    bool isValid() {
        return (tId != std::numeric_limits<size_t>::max() && fId != std::numeric_limits<size_t>::max());
    }
};

struct Operation {
    std::string name;
    size_t vid;
    std::vector<size_t> vids;
    bool clockwise;
    glm::dvec3 coords = {0, 0, 0};

    Operation(std::string name_, size_t vid_, std::vector<size_t> vids_, bool clockwise_, glm::dvec3 coords_ = {0, 0, 0}) {
        name = name_;
        vid = vid_;
        vids = vids_;
        clockwise = clockwise_;
        coords = coords_;
    }

    bool isValid() {
        if (name == "Flip" || name == "Collapse" || name == "Split") return vids.size() == 2;
        return true;
    };
};

class SemiGlobalSimplifier {
    public:
        // Constructors and Destructor
        SemiGlobalSimplifier();
        SemiGlobalSimplifier(Mesh& mesh_, MeshUtil& mu_, Smoother& smoother_, KDTree& kdTree_);
        ~SemiGlobalSimplifier();

        // MeshUtil setters and getters
        void SetMembers(Mesh& mesh_, MeshUtil& mu_, Smoother& Smoother_);
        void SetIters(int iters_);
        void SetFaceMetrics();

        // Simplification Operations
        bool FixBoundary();
        bool RemoveDoublets();
        bool FixValences();
        void SetSimplificationOperations();
        void SetDiagonalCollapseOperations();
        bool SetBoundaryDirectSeparatrixOperations(bool looseCollapse);
        bool SetDirectSeparatrixOperations(bool looseCollapse);
        void SetSeparatrixOperations();
        void SetBoundarySeparatrixOperations();
        void SetHalfSeparatrixOperations();
        void SetChordCollapseOperations();
        void SetEdgeRotationOperations();
        void SetVertexRotationOperations();
        void SetEdgeCollapseOperations();
        void SetVertexSplitOperations();
        void SetQuadSplitOperations();
        void GetSingularityPairs();
        void SetSingularityLinks(std::vector<SingularityLink>& SingularityLinks, std::vector<std::vector<size_t>>& SingularityMap, BaseComplexQuad& bc);
        void SelectSingularityGroups(std::vector<SingularityGroup>& Groups, std::vector<SingularityLink>& SingularityLinks, std::vector<std::vector<size_t>>& SingularityMap);
        void ResolveSingularityGroups(std::vector<SingularityGroup>& Groups, BaseComplexQuad& bc);
        void ResolveSingularityPairs();
        void CheckAndResolveThreeFivePair(size_t vid);
        void CheckAndResolveFiveThreePair(size_t vid);
        std::vector<SingularityLink> TraceSingularityLinks(Vertex& v, BaseComplexQuad& bc);
        void TraceSingularityLinks(size_t vid, std::vector<SingularityLink>& links, bool traceDiagonals = false);
        void TraceLink(const Vertex& v, const Edge& edge, std::vector<SingularityLink>& links, std::vector<size_t> vids = {}, std::vector<size_t> eids = {});
        void TraceLink(const Vertex& v, const Face& face, std::vector<SingularityLink>& links, std::vector<size_t> vids = {});
        void TraceAlongEdge(const Vertex& start_vertex, const Edge& start_edge, std::vector<size_t>& vids_link, bool checkBoundary = true);
        void TraceAlongEdge(const Vertex& start_vertex, const Edge& start_edge, std::vector<size_t>& vids_link, std::vector<size_t>& eids_link);
        void TraceAlongDiagonal(const Vertex& start_vertex, const Face& start_face, std::vector<size_t>& vids_link);
        void GetDiagonalPath(SingularityLink& l);
        int GetEdgeRots(size_t eid1, size_t eid2);
        int GetFaceRots(size_t fid1, size_t fid2, size_t vid);
        void SelectLinks(std::vector<SingularityLink>& links);
        void SelectDirectPairLink(std::vector<size_t> threeFiveIds, std::vector<SingularityLink>& links, std::vector<size_t> verticesToAvoid = {});
        void SelectDiagonalPairLink(std::vector<size_t> threeFiveIds, std::vector<SingularityLink>& links, std::vector<size_t> verticesToAvoid = {});
        void PerformGlobalOperations();
        void Smooth(vMesh* m = nullptr);

        bool ResolveHighValences();
        void AlignSingularities();
        void AlignSingularities(size_t vid, std::queue<size_t>& Singularities, std::vector<bool>& isAvailable, BaseComplexQuad& bc, bool checkValence = true);
        void AlignAndResolveSingularities(bool checkValence = true);
        bool IsPair(size_t vid);
        std::vector<size_t> GetPair(size_t vid);
        int MovePair(std::vector<size_t> threeFiveIds, std::vector<size_t>& secondaryPath, bool checkValence = false);
        bool ResolveIsolatedSingularities(BaseComplexQuad& bc);
        void GenerateSingularityPair(SingularityLink& l1, SingularityLink& l2);
        void ResolveSingularities();
        void GetSingularityGroups(std::vector<size_t> Singularities, BaseComplexQuad& bc);
        int PullSingularity(SingularityLink& l1, SingularityLink& l2);
        std::vector<int> GetTraverseInfo(SingularityLink& l1, SingularityLink& l2);
        std::vector<size_t> TraversePath(size_t prev, size_t current, std::vector<int> Rots);
        std::vector<size_t> GetSecondaryPath(int offset, std::vector<size_t>& mainPath, BaseComplexQuad& bc);
        std::vector<SingularityLink> GetLinks(size_t sid, BaseComplexQuad& bc, bool checkValence = true);
        std::vector<SingularityLink> SelectLinks(std::vector<SingularityLink> links, int valence, bool checkValence = true);
        std::vector<SingularityLink> GetCrossLinks(SingularityLink& l, BaseComplexQuad& bc, bool checkValence = true);
        bool ValidateLink(SingularityLink& l);
        bool ValidatePath(std::vector<size_t> p, bool checkValences = true);
        int ResolveSingularity(size_t sid, BaseComplexQuad& bc);
        int MoveSingularity(SingularityGroup& sg);
        int MoveSingularity(int offset, std::vector<size_t>& mainPath, std::vector<size_t>& secondaryPath);
        bool IsExclusive(SingularityLink& l1, SingularityLink& l2);
        bool IsExclusive(size_t vid, std::vector<size_t> a, std::vector<size_t> b);

        bool doesCrossBoundary(std::vector<size_t> in, bool isVertex);

        void PrototypeBoundary(bool checkValence = true);
        void PrototypeA(size_t vid, BaseComplexQuad& bc, bool checkValence = true);
        void PrototypeB();
        void PrototypeC();
        void PrototypeD();
        void PrototypeE();
        void PrototypeF(int idxOffset = 2);
        void PrototypeG(int vid, BaseComplexQuad& bc);
        bool PrototypeH(int idxOffset = 2);
        bool PrototypeI();
        void PrototypeJ();
        void PrototypeMoveSingularity(SingularityLink& l);

        SingularityLink PrototypeGetLink(size_t vid);
        SingularityLink PrototypePairLink(std::vector<size_t> threeFiveIds, std::vector<size_t>& verticesToAvoid, bool isDiagonal = false);
        bool PrototypePairIds(std::vector<size_t>& threeFiveIds, std::vector<size_t>& vids, size_t toMoveIdx, size_t destIdx);
        void PrototypeResolvePairIds(std::vector<size_t>& threeFiveIds, std::vector<size_t>& vids, size_t toMoveIdx, size_t destIdx);
        void PrototypeMovePair(std::vector<size_t>& threeFiveIds, SingularityLink& l, bool isDiagonal = false);
        bool PrototypeIsLinkValid(SingularityLink& l, std::vector<size_t> verticesToAvoid = {});

        int PrototypeGetRotations(size_t vid, size_t start, size_t end);
        void PrototypeSaveMesh(SingularityLink& l1, SingularityLink& l2, std::string in);
        void PrototypeSaveMesh(const std::vector<SingularityLink>& links, std::string in);
        void PrototypeResolveGroup(SingularityLink& l1, SingularityLink& l2);
        int PrototypeCancelThreeFivePair(SingularityLink& l1, SingularityLink& l2);
        int PrototypeCancelSingularity(size_t vid, BaseComplexQuad& bc);
        bool PrototypeCheckBoundarySingularity(size_t vid);
        int PrototypeCancelSingularityPair(SingularityLink& l, BaseComplexQuad& bc);
        SingularityLink PrototypeGetLink(size_t vid, BaseComplexQuad& bc, size_t vertexToSkip = 0, std::vector<size_t> edgesToCheck = {}, bool checkValence = true, bool boundary = true);
        bool PrototypeIsLinkValid(SingularityGroup& s);
        int PrototypeGetElementPrediction(SingularityGroup& s);


        void PrototypeK();
        void PrototypeTraceSeparatrix(size_t vid, size_t eid, Separatrix& s, bool checkValence = false);
        void PrototypeTransportSeparatrix(Separatrix& s, bool includeIters, std::vector<std::vector<size_t>>& vec);
        void PrototypeTransportDiagonalPair(size_t threeId, size_t fiveId, bool includeIters, std::vector<size_t> verticesToAvoid = {});
        void PrototypeTraceDiagSeparatrix(size_t vid, size_t eid, Separatrix& s, bool checkBoundary, std::vector<size_t> verticesToAvoid);
        void PrototypeTransportDirectPair(size_t threeId, size_t fiveId, bool includeIters, std::vector<size_t> verticesToAvoid = {});
        void PrototypeTraceDirectSeparatrix(size_t vid, size_t eid, Separatrix& s, std::vector<size_t> verticesToAvoid, int rotOffset, bool isRotEdge);
        void PrototypeTransportDirectLink(Separatrix& primary, std::vector<std::vector<size_t>>& vec);
        void PrototypeTraceDirectLinkSeparatrix(size_t vid, size_t eid, Separatrix& s, size_t vertexToAvoid, int rotOffset, bool isRotEdge);
        void PrototypeMoveLinkPair(Separatrix& p, Separatrix& s, size_t singularity, size_t moveDir, size_t sourceDir);
        void PrototypeTransportLink(Separatrix& primary);
        void PrototypeMoveLinkPair(size_t singularityId, size_t moveDir, std::vector<size_t> verticesToAvoid);
        size_t GetEdgeId(size_t vid, size_t vid2);
        size_t GetNextCCedge(size_t vid, size_t startEid);
        size_t GetCCedgeAt(size_t vid, size_t eid, int counter);
        void TraceAlongEdge(const Vertex& start_vertex, const Edge& start_edge, std::vector<size_t>& vids_link, int b1, int b2, int rotOffset, bool checkBoundary = true);
        void PrototypeSaveSeparatrices(std::vector<std::vector<size_t>> vec, std::string in = "test");

        void SaveEdgeMesh();


        void PrototypeExecute();
        std::vector<bool> PrototypeSetSingularities();
        void PrototypeExtractSingularityLinks(size_t vid, std::vector<Separatrix>& separatrices);
        void PrototypeTraceLinks(size_t vid, std::vector<Separatrix>& seps);
        bool PrototypeGetPairIds(size_t singularityId, size_t moveDir, std::vector<size_t>& threeFiveIds, int secVid = -1);
        void PrototypeSetDirMap(Separatrix& s);
        void PrototypeSet33DirMap(Separatrix& s);
        void PrototypeSet55DirMap(Separatrix& s);
        void PrototypeSet35DirMap(Separatrix& s);
        std::pair<size_t,double> PrototypeGetDir(size_t vid, std::vector<Separatrix>& separatrices);
        Separatrix* PrototypeGetSeparatrix(size_t vid, std::vector<Separatrix>& separatrices, std::pair<size_t,double>& dirMap);
        void PrototypeMove(Separatrix& s);
        
        void TraceAlongEdge(const Vertex& start_vertex, const Edge& start_edge, bool& crossBoundary, int& boundary_distance, std::vector<size_t>& vids_link, bool checkBoundary = false);
        void TraceAlongDiagonal(const Vertex& start_vertex, const Face& start_face, int& boundary_distance, std::vector<size_t>& vids_link);
        size_t GettCCFaceAt(size_t vid, size_t eid, int counter);


        void func1();
        std::vector<Path> PrototypeGetPaths();
        void PrototypeGetPathsAt(size_t vid, std::vector<Path>& paths);
        void PrototypeMoveAlongPath(Path& p);
        

        bool PerformOperation(const Operation& op, vMesh* m);
        void MovePair(tfPair& p, size_t dest, vMesh& m);
        bool TestFlips();

        bool CheckMeshValidity();

        void RenderMesh() {
            renderer.Render();
            std::this_thread::sleep_for(std::chrono::milliseconds(1000));
        }

        MeshUtil* mu;
        int iters = 0;
        
    private:
        Mesh* mesh;
        Smoother* smoother;
        KDTree* kd;
        std::unique_ptr<SurfaceProjector> sp;
        // SurfaceProjector sp;
        MeshCaretaker caretaker;
        Renderer renderer;

        void CheckValidity();
        size_t GetFaceID(size_t vid, size_t exclude_vid);
        size_t GetFaceId(size_t vid1, size_t vid2);
        size_t GetDiagonalV(size_t vid, size_t fid);
        size_t GetFaceV(size_t vid, size_t fid, int offset);
        bool Contains(std::vector<size_t> v, size_t val);
        bool Contains(std::vector<size_t> v, std::vector<size_t> v2);
        std::vector<size_t> GetThreeFivePairIds(size_t vid, size_t mainId, size_t secondaryId);
        void MoveSingularities(size_t& toMoveId, size_t& sourceId, size_t& secondaryId, size_t& sourceDir, size_t& secondaryDir, std::vector<size_t> secondaryPath);
        void SetSecondaryPath(size_t& secondaryId, size_t& toMoveId, size_t& sourceId, std::vector<size_t>& mainPath, std::vector<size_t>& secondaryPath, BaseComplexQuad& bc);
        std::vector<std::shared_ptr<SimplificationOperation>> Ops;
        PQueue<std::shared_ptr<SimplificationOperation>> Op_Q;
        std::recursive_mutex mtx;
        std::mutex mx;
        double delta = 0.0;

};

#endif