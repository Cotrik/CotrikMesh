/*
 * ExtractQuadMesh.cpp
 *
 *  Created on: Oct 20, 2018
 *      Author: cotrik
 */

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "FrameField.h"
#include "PolyLine.h"
#include "FrameOpt.h"
#include "ArgumentManager.h"
#include <iostream>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <stack>
#include <algorithm>

size_t sqsize = 100;
Mesh triMesh;
Mesh quadMesh;
std::string space = "xy";
struct VertexUV : public Vertex {
    size_t father;
    glm::dvec2 uv;
};
void read(const char* filename, std::vector<VertexUV>& V, std::vector<Face>& F) {
    std::ifstream ifs(filename);
    std::string line;
    while (getline(ifs, line)) {
        std::stringstream ss(line);
        std::string keyword;
        VertexUV v;
        ss >> keyword;
        if (keyword == "Vertex") {
            ss >> v.id >> v.x >> v.y >> v.z;
            --v.id;
            std::string father, uv, vv;
            ss >> father >> uv >> vv;

            for (auto& c : father)
                if (c == '(' || c == ')') c = ' ';
            std::istringstream ss_father(father);
            Vertex u_v;
            ss_father >> father >> v.father;
            --v.father;

            for (auto& c : uv)
                if (c == '(' || c == ')') c = ' ';
            std::istringstream ss_uv(uv);
            ss_uv >> uv >> v.uv.x;
//            std::string x,y;
//            ss_uv >> uv >> x >> y;
//            v.uv.x = std::stod(x);
//            v.uv.y = std::stod(y);
            //ss_uv >> uv >> v.uv.x >> v.uv.y;
            for (auto& c : vv)
                if (c == '(' || c == ')') c = ' ';
            std::istringstream ss_vv(vv);
            ss_vv >> v.uv.y;

            V.push_back(v);
        } else if (keyword == "Face") {
            Face f;
            ss >> f.id;
            --f.id;
            size_t vid;
            while (ss >> vid)
                f.Vids.push_back(vid - 1);

            F.push_back(f);
        } else if (keyword == "Edge") {
            break;
        }
    }
}

void getLines(const std::vector<VertexUV>& V, std::vector<std::vector<size_t>>& L) {

}

std::set<size_t> get_boundary_vids(const Mesh& uvtriMesh);
std::vector<size_t> get_singular_vids(const Mesh& uvtriMesh, const std::set<size_t>& boundary_vids);
void write(const char* filename, const std::vector<VertexUV>& V, const std::vector<Face>& F) {
    std::vector<Vertex> triV;
    Vertex x;
    for (auto& v : V) {
        x.id = v.id;
        x.x = v.x;
        x.y = v.y;
        x.z = v.z;
        triV.push_back(x);
    }
    std::vector<Cell> C;
    Cell c;
    for (auto& f : F) {
        c.id = f.id;
        c.Vids = f.Vids;
        C.push_back(c);
    }
    MeshFileWriter writer(triV, C, filename, F.front().Vids.size() == 3 ? TRIANGLE : QUAD);
    writer.WriteFile();
}

void generate_tri_mesh(const std::vector<VertexUV>& V, const std::vector<Face>& F, Mesh& triMesh) {
    std::vector<Vertex> triV;
    Vertex x;
    for (auto& v : V) {
        x.id = v.id;
        x.x = v.x;
        x.y = v.y;
        x.z = v.z;
        triV.push_back(x);
    }
    std::vector<Cell> C;
    Cell c;
    for (auto& f : F) {
        c.id = f.id;
        c.Vids = f.Vids;
        C.push_back(c);
    }

    triMesh.V = triV;
    triMesh.C = C;
    triMesh.m_cellType = TRIANGLE;
    triMesh.BuildAllConnectivities();
    triMesh.ExtractBoundary();
}

void generate_tri_mesh(const std::vector<Vertex>& V, const std::vector<Face>& F, Mesh& triMesh) {
    std::vector<Vertex> triV;
    Vertex x;
    for (auto& v : V) {
        x.id = v.id;
        x.x = v.x;
        x.y = v.y;
        x.z = v.z;
        triV.push_back(x);
    }
    std::vector<Cell> C;
    Cell c;
    for (auto& f : F) {
        c.id = f.id;
        c.Vids = f.Vids;
        C.push_back(c);
    }

    triMesh.V = triV;
    triMesh.C = C;
    triMesh.m_cellType = TRIANGLE;
    triMesh.BuildAllConnectivities();
    triMesh.ExtractBoundary();
}

void generate_quad_mesh(const std::vector<Vertex>& V, const std::vector<Face>& F, Mesh& quadMesh) {
    std::vector<Cell> C;
    Cell c;
    for (auto& f : F) {
        c.id = f.id;
        c.Vids = f.Vids;
        C.push_back(c);
    }

    quadMesh.V = V;
    quadMesh.C = C;
    quadMesh.m_cellType = QUAD;
    quadMesh.BuildAllConnectivities();
    quadMesh.ExtractBoundary();
}

void write(const char* filename, const std::vector<Vertex>& V, const std::vector<Face>& F) {
    std::vector<Vertex> triV;
    Vertex x;
    for (auto& v : V) {
        x.id = v.id;
        x.x = v.x;
        x.y = v.y;
        x.z = v.z;
        triV.push_back(x);
    }
    std::vector<Cell> C;
    Cell c;
    for (auto& f : F) {
        c.id = f.id;
        c.Vids = f.Vids;
        C.push_back(c);
    }
    MeshFileWriter writer(triV, C, filename, F.front().Vids.size() == 3 ? TRIANGLE : QUAD);
    //writer.FixMesh();
    writer.WriteFile();
}

#include "glm/gtc/matrix_access.hpp"
bool IsVertextInTriangle(const Vertex& p, const Vertex& p0, const Vertex& p1, const Vertex& p2, glm::dvec2& uv) {
    glm::dvec2 v02((p0.x - p2.x), (p0.y - p2.y));
    glm::dvec2 v12((p1.x - p2.x), (p1.y - p2.y));

    glm::dvec2 r(p.x, p.y);
    glm::dvec2 r3(p2.x, p2.y);

    glm::dvec2 rr3 = r - r3;

    glm::dmat2x2 T(v02, v12);
    glm::dmat2x2 T_inverse = glm::inverse(T);
    glm::dvec2 lambda = T_inverse * (r - r3);

    if (lambda.x > -1e-5 && lambda.y > -1e-5 && (lambda.x + lambda.y) < 1.00001) {
        uv.x = lambda.x;
        uv.y = lambda.y;
        return true;
    }

    return false;
}


bool IsVertextInTriangle(const Vertex& p, const VertexUV& p0, const VertexUV& p1, const VertexUV& p2, glm::dvec2& uv) {
    glm::dvec2 v02((p0.uv.x - p2.uv.x), (p0.uv.y - p2.uv.y));
    glm::dvec2 v12((p1.uv.x - p2.uv.x), (p1.uv.y - p2.uv.y));

    glm::dvec2 r(p.x, p.y);
    glm::dvec2 r3(p2.uv.x, p2.uv.y);

    glm::dvec2 rr3 = r - r3;

    glm::dmat2x2 T(v02, v12);
    glm::dmat2x2 T_inverse = glm::inverse(T);
    glm::dvec2 lambda = T_inverse * (r - r3);

    if (lambda.x > -1e-5 && lambda.y > -1e-5 && (lambda.x + lambda.y) < 1.00001) {
        uv.x = lambda.x;
        uv.y = lambda.y;
        return true;
    }

    return false;
}

void generate_quad_mesh(std::vector<Vertex>& V, std::vector<Face>& F) {
    V.resize((sqsize + 1) * (sqsize + 1));
    F.reserve(sqsize * sqsize);
    Face f;
    size_t fid = 0;
    for (int i = 0; i <= sqsize; ++i) {
        double x = double(i) / sqsize;
        for (int j = 0; j <= sqsize; ++j) {
            double y = double(j) / sqsize;
            size_t id = i * (sqsize + 1) + j;
            V[id] = glm::dvec3(x, y, 0.0);
            V[id].id = id;
            if (i != 0 && j != 0) {
                f.Vids = {id - sqsize - 2, id - sqsize - 1, id, id - 1};
                f.id = fid++;
                F.push_back(f);
            }
        }
    }
}

void generate_quad_mesh2(std::vector<Vertex>& V, std::vector<Face>& F) {
    V.resize((sqsize + 1) * (sqsize + 1));
    F.reserve(sqsize * sqsize);
    Face f;
    size_t fid = 0;
    for (int i = 0; i <= sqsize; ++i) {
        double x = 2 * double(i) / sqsize - 1.0;
        for (int j = 0; j <= sqsize; ++j) {
            double y = 2 * double(j) / sqsize - 1.0;
            size_t id = i * (sqsize + 1) + j;
            V[id] = glm::dvec3(x, y, 0.0);
            V[id].id = id;
            if (i != 0 && j != 0) {
                f.Vids = {id - sqsize - 2, id - sqsize - 1, id, id - 1};
                f.id = fid++;
                F.push_back(f);
            }
        }
    }
}

void generate_quad_mesh(std::vector<Vertex>& V, std::vector<Face>& F, double minx, double miny, double maxx, double maxy) {
    double lenx = maxx - minx;
    double leny = maxy - miny;
    double len = (lenx > leny) ? leny : lenx;
    double d = len / sqsize;
    double offset = d / 100;
    minx += offset;
    miny += offset;
    maxx -= offset;
    maxy -= offset;
    int nx = lenx / d;
    int ny = leny / d;
    V.resize((nx + 1) * (ny + 1));
    F.reserve(nx * ny);
    Face f;
    size_t fid = 0;
    for (int i = 0; i <= nx; ++i) {
        double x = lenx * double(i) * d;
        for (int j = 0; j <= ny; ++j) {
            double y = double(j) * d;
            size_t id = i * (ny + 1) + j;
            V[id] = glm::dvec3(x, y, 0.0);
            V[id].id = id;
            if (i != 0 && j != 0) {
                f.Vids = {id - ny - 2, id - ny - 1, id, id - 1};
                f.id = fid++;
                F.push_back(f);
            }
        }
    }
}

struct Parameters {
    Parameters() {};
    Parameters(glm::dvec2 uv, size_t triangleid) : uv(uv), triangleid(triangleid) {}
    glm::dvec2 uv;
    size_t triangleid;
};

struct MapPoint {
    MapPoint() {};
    MapPoint(const glm::dvec2 uv, const glm::dvec3 xyz, size_t tri_id) : uv(uv), xyz(xyz), tri_id(tri_id) {}
    //MapPoint(const glm::dvec2 uv, const glm::dvec3 xyz, size_t tri_id) : uv(uv), xyz(xyz), tri_id(tri_id) {}
    glm::dvec2 uv;
    glm::dvec3 xyz;
    size_t tri_id;
};

void MapbackToOrigTri(const std::vector<VertexUV>& triV, const std::vector<Face>& triF,
        const std::vector<Vertex>& quadV, const std::vector<Face>& quadF,
        const std::vector<std::vector<Parameters>>& p, std::vector<Vertex>& new_quadV) {
    //parameterization of new OrigTet to new Cube Mesh
    std::cout << "MapbackToOrigTri" << std::endl;
    new_quadV.resize(quadV.size());
    for (int i = 0; i < quadV.size(); i++) {
        auto& params = p[i];
        if (params.empty()) {
            continue;
        }
        for (auto& param: params) {
            const auto triangleid = param.triangleid;
            const auto& triVids = triF.at(triangleid).Vids;

            const Vertex& origT0 = triV[triVids.at(0)];
            const Vertex& origT1 = triV[triVids.at(1)];
            const Vertex& origT2 = triV[triVids.at(2)];

    //        glm::dvec3 p02((origT0.x - origT2.x), (origT0.y - origT2.y), (origT0.z - origT2.z));
    //        glm::dvec3 p12((origT1.x - origT2.x), (origT1.y - origT2.y), (origT0.z - origT2.z));
    //        glm::dmat2x3 OrigT(p02, p12);
    //        glm::dvec3 new_r3(origT2.x, origT2.y, origT2.z);

    //        const glm::dvec2& rambda = param.uv;
    //        glm::dvec3 new_r = OrigT * rambda;
    //        new_r += new_r3;
    //        Vertex baryCenter(new_r.x, new_r.y, new_r.z);
            const double rambda0 = param.uv.x;
            const double rambda1 = param.uv.y;
            const double rambda2 = 1.0 - rambda0 - rambda1;

            Vertex& v = new_quadV.at(i);
    //        v.x = baryCenter.x;
    //        v.y = baryCenter.y;
    //        v.z = baryCenter.z;
            v.x = rambda0 * origT0.x + rambda1 * origT1.x + rambda2 * origT2.x;
            v.y = rambda0 * origT0.y + rambda1 * origT1.y + rambda2 * origT2.y;
            v.z = rambda0 * origT0.z + rambda1 * origT1.z + rambda2 * origT2.z;
        }
    }
}

void MapbackToOrigTri(const std::vector<VertexUV>& triV, const std::vector<Face>& triF,
        std::vector<Vertex>& quadV, std::vector<Face>& quadF,
        const std::vector<std::vector<Parameters>>& p, const Face& overlapFace,
        std::vector<Vertex>& new_quadV) {
    std::vector<size_t> triangleids;
    for (auto& param: p[overlapFace.Vids[0]]) {
        triangleids.push_back(param.triangleid);
    }
    int count = 0;
    for (auto triangleid : triangleids) {
        for (auto vid : overlapFace.Vids) {
            auto& params = p[vid];
            for (auto& param: params) {
                const auto triangleid_ = param.triangleid;
                if (triangleid_ != triangleid) continue;
                const auto& triVids = triF.at(triangleid).Vids;

                const Vertex& origT0 = triV[triVids.at(0)];
                const Vertex& origT1 = triV[triVids.at(1)];
                const Vertex& origT2 = triV[triVids.at(2)];
                const double rambda0 = param.uv.x;
                const double rambda1 = param.uv.y;
                const double rambda2 = 1.0 - rambda0 - rambda1;

                Vertex v;
                v.x = rambda0 * origT0.x + rambda1 * origT1.x + rambda2 * origT2.x;
                v.y = rambda0 * origT0.y + rambda1 * origT1.y + rambda2 * origT2.y;
                v.z = rambda0 * origT0.z + rambda1 * origT1.z + rambda2 * origT2.z;
                v.id = new_quadV.size();
                new_quadV.push_back(v);
            }
        }
        Face f;
        f.Vids = {new_quadV.size() - 4, new_quadV.size() - 3, new_quadV.size() - 2, new_quadV.size() - 1};
        quadF.push_back(f);
        ++count;
    }
}

std::set<size_t> get_triangle_neighbor_triangleids(const Mesh& triMesh, size_t trangleid) {
    std::set<size_t> res;
    for (auto vid : triMesh.F[trangleid].Vids) {
        auto& triangleids = triMesh.V[vid].N_Fids;
        res.insert(triangleids.begin(), triangleids.end());
    }
//    auto x = res;
//    for (auto fid : x) {
//        for (auto vid : triMesh.F[fid].Vids){
//            auto& triangleids = triMesh.V[vid].N_Fids;
//            res.insert(triangleids.begin(), triangleids.end());
//        }
//    }
    return res;
}
void MapbackToOrigTri2(const std::vector<VertexUV>& triV, const std::vector<Face>& triF,
        std::vector<Vertex>& quadV, std::vector<Face>& quadF,
        const std::vector<std::vector<Parameters>>& p, const Face& overlapFace,
        std::vector<Vertex>& new_quadV) {

    std::vector<size_t> triangleids;
    for (auto& param: p[overlapFace.Vids[0]]) {
        triangleids.push_back(param.triangleid);
    }
    size_t target_triangleid = triangleids.front();
    for (auto triangleid : triangleids) {
        auto triangle_neighbor_triangleids = get_triangle_neighbor_triangleids(triMesh, triangleid);
        bool all_quadvids_in_triangle_neighbor_triangleids = true;
        for (auto vid : overlapFace.Vids) {
            auto& params = p[vid];
            bool found = false;
            for (auto& param: params) {
                if (triangle_neighbor_triangleids.find(param.triangleid) != triangle_neighbor_triangleids.end()) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                all_quadvids_in_triangle_neighbor_triangleids = false;
                break;
            }
        }
        if (all_quadvids_in_triangle_neighbor_triangleids) {
            target_triangleid = triangleid;
            break;
        }
    }

    int count = 0;
    //size_t triangleid = target_triangleid;
    auto triangle_neighbor_triangleids = get_triangle_neighbor_triangleids(triMesh, target_triangleid);

    {
        for (auto vid : overlapFace.Vids) {
            auto& params = p[vid];
            for (auto& param: params) {
                for (auto triangleid : triangle_neighbor_triangleids) {
                    const auto triangleid_ = param.triangleid;
                    if (triangleid_ != triangleid) continue;
                    const auto& triVids = triF.at(triangleid).Vids;

                    const Vertex& origT0 = triV[triVids.at(0)];
                    const Vertex& origT1 = triV[triVids.at(1)];
                    const Vertex& origT2 = triV[triVids.at(2)];
                    const double rambda0 = param.uv.x;
                    const double rambda1 = param.uv.y;
                    const double rambda2 = 1.0 - rambda0 - rambda1;

                    Vertex v;
                    v.x = rambda0 * origT0.x + rambda1 * origT1.x + rambda2 * origT2.x;
                    v.y = rambda0 * origT0.y + rambda1 * origT1.y + rambda2 * origT2.y;
                    v.z = rambda0 * origT0.z + rambda1 * origT1.z + rambda2 * origT2.z;
                    v.id = new_quadV.size();
                    new_quadV.push_back(v);
                }
            }
        }
        Face f;
        f.Vids = {new_quadV.size() - 4, new_quadV.size() - 3, new_quadV.size() - 2, new_quadV.size() - 1};
        quadF.push_back(f);
        ++count;
    }
}

void MapbackToOrigTri(const std::vector<VertexUV>& triV, const std::vector<Face>& triF,
        std::vector<Vertex>& quadV, std::vector<Face>& quadF,
        const std::vector<std::vector<Parameters>>& p, const Face& overlapFace, const size_t target_triangleid,
        std::vector<Vertex>& new_quadV) {
    auto triangleids = get_triangle_neighbor_triangleids(triMesh, target_triangleid);
    for (auto vid : overlapFace.Vids) {
        auto& params = p[vid];
        for (auto& param: params) {
            for (auto triangleid : triangleids) {
                const auto triangleid_ = param.triangleid;
                if (triangleid_ != triangleid) continue;
                const auto& triVids = triF.at(triangleid).Vids;

                const Vertex& origT0 = triV[triVids.at(0)];
                const Vertex& origT1 = triV[triVids.at(1)];
                const Vertex& origT2 = triV[triVids.at(2)];
                const double rambda0 = param.uv.x;
                const double rambda1 = param.uv.y;
                const double rambda2 = 1.0 - rambda0 - rambda1;

                Vertex v;
                v.x = rambda0 * origT0.x + rambda1 * origT1.x + rambda2 * origT2.x;
                v.y = rambda0 * origT0.y + rambda1 * origT1.y + rambda2 * origT2.y;
                v.z = rambda0 * origT0.z + rambda1 * origT1.z + rambda2 * origT2.z;
                v.id = new_quadV.size();
                new_quadV.push_back(v);
            }
        }
    }
    Face f;
    f.Vids = {new_quadV.size() - 4, new_quadV.size() - 3, new_quadV.size() - 2, new_quadV.size() - 1};
    quadF.push_back(f);
}

void MapbackToOrigTri(const std::vector<VertexUV>& triV, const std::vector<Face>& triF,
        const std::map<size_t, Parameters>& p, std::vector<Vertex>& new_quadV) {
    for (auto& item : p) {
        auto& param = item.second;
        const auto triangleid = param.triangleid;
        const auto& triVids = triF.at(triangleid).Vids;

        const Vertex& origT0 = triV[triVids.at(0)];
        const Vertex& origT1 = triV[triVids.at(1)];
        const Vertex& origT2 = triV[triVids.at(2)];
        const double rambda0 = param.uv.x;
        const double rambda1 = param.uv.y;
        const double rambda2 = 1.0 - rambda0 - rambda1;

        Vertex& v = new_quadV.at(item.first);
        v.x = rambda0 * origT0.x + rambda1 * origT1.x + rambda2 * origT2.x;
        v.y = rambda0 * origT0.y + rambda1 * origT1.y + rambda2 * origT2.y;
        v.z = rambda0 * origT0.z + rambda1 * origT1.z + rambda2 * origT2.z;
    }
}

void get_extend(const std::vector<Face>& quadF, const std::vector<Face>& nonoverlap_quadF, const std::set<size_t>& inside_quadF_set,
        std::set<size_t>& nonoverlapquad_extendquadids, std::set<size_t>& nonoverlapquad_extendquad_vids) {
    for (auto& f : nonoverlap_quadF)
        for (auto vid : f.Vids) {
            const auto& v = quadMesh.V.at(vid);
            const auto& fids = v.N_Fids;
            nonoverlapquad_extendquadids.insert(fids.begin(), fids.end());
        }
    for (auto& f : nonoverlap_quadF)
        nonoverlapquad_extendquadids.erase(f.id);
    std::vector<size_t> needtoerase;
    for (auto& fid : nonoverlapquad_extendquadids)
        if (inside_quadF_set.find(fid) == inside_quadF_set.end())
            needtoerase.push_back(fid);
    for (auto fid : needtoerase)
        nonoverlapquad_extendquadids.erase(fid);

    for (auto fid : nonoverlapquad_extendquadids) {
        auto& f = quadF.at(fid);
        nonoverlapquad_extendquad_vids.insert(f.Vids.begin(), f.Vids.end());
    }
}

void get_inside_quads(const std::vector<VertexUV>& triV, const std::vector<Face>& triF, const std::vector<Vertex>& quadV,
        std::vector<std::vector<Parameters>>& vertex_triangles,
        std::map<size_t, size_t>& overlap_quad_vids, std::vector<bool>& vertex_inside) {
    for (auto& v : quadV) {
        bool inside = false;
        int count = 0;
        for (auto& tri : triF) {
            glm::dvec2 uv;
            if (IsVertextInTriangle(v, triV[tri.Vids[0]], triV[tri.Vids[1]], triV[tri.Vids[2]], uv)) {
                vertex_triangles[v.id].push_back(Parameters(uv, tri.id));
                inside = true;
                ++count;
                //break;
            }
        }
        if (inside) vertex_inside[v.id] = true;
        if (count > 1) {
            overlap_quad_vids[v.id] = count;
            // std::cout << v.id << " --> " << count << "\n";
        }
    }
}

void get_inside_quads(const std::vector<Face>& quadF, const std::vector<bool>& vertex_inside,
        std::vector<Face>& inside_quadF, std::vector<bool>& quad_inside) {
    for (auto& f : quadF) {
        bool all_inside = true;
        for (auto vid : f.Vids) {
            if (!vertex_inside[vid]) {
                all_inside = false;
                break;
            }
        }
        if (all_inside) quad_inside[f.id] = true;
        if (all_inside) inside_quadF.push_back(f);
    }
}

void get_next(const std::map<size_t, Parameters>& nonoverlap_Vids, std::map<size_t, Parameters>& nonoverlapquad_Vids,
        const std::map<size_t, Parameters>& remaining_nonoverlap_Vids, const std::vector<std::vector<Parameters>>& vertex_triangles,
        const std::set<size_t>& inside_quadF_set,
        std::map<size_t, Parameters>& all_nextvids, std::set<size_t>& nextvids, std::set<size_t>& nextquads) {
    for (auto item : remaining_nonoverlap_Vids) {
        // std:: cout << "neighbor quads = " << quadMesh.V[vid].N_Fids.size() << "\n";
        const auto vid = item.first;
        for (auto fid : quadMesh.V[vid].N_Fids) {
            bool is_quad_contains_vid_in_nonoverlapquad_Vids = false;

            auto& n_f = quadMesh.F.at(fid);
            size_t nonoverlapquad_V_triangleid = 0;
            for (auto n_vid : n_f.Vids)
                if (nonoverlapquad_Vids.find(n_vid) != nonoverlapquad_Vids.end()) {
                    is_quad_contains_vid_in_nonoverlapquad_Vids = true;
                    nonoverlapquad_V_triangleid = nonoverlapquad_Vids[n_vid].triangleid;
                    break;
                }
            if (is_quad_contains_vid_in_nonoverlapquad_Vids && inside_quadF_set.find(fid) != inside_quadF_set.end() && nonoverlapquad_V_triangleid < triMesh.F.size()) {
                nextquads.insert(fid);
                auto triangleids = get_triangle_neighbor_triangleids(triMesh, nonoverlapquad_V_triangleid);
//                    std::cout << "----------\n";
//                    std::cout << "---trianglceids\n";
//                    for (auto id : triangleids) {
//                        std::cout << id << " ";
//                    }
//                    std::cout << "\n";
//                    std::cout << "---quad vids\n";
//                    for (auto n_vid : n_f.Vids) {
//                        auto& params = vertex_triangles[n_vid];
//                        for (auto& param : params)
//                            std::cout << param.triangleid << " ";
//                    }
//                    std::cout << "\n\n";
                for (auto n_vid : n_f.Vids) {
                    if (nonoverlap_Vids.find(n_vid) == nonoverlap_Vids.end()) {
                        auto& params = vertex_triangles[n_vid];
                        if (!params.empty()) nextvids.insert(n_vid);
                        //std::cout << "params.size() = " << params.size() << "\n";
                        for (auto& param : params) {
                            auto it = triangleids.find(param.triangleid);
                            if (it != triangleids.end()) {
                                all_nextvids.insert(std::make_pair(n_vid, param));
                            }
                        }
                    }
                }
            }
        }
    }
}

std::set<size_t> get_nonoverlap_quadEids(const Mesh& quadMesh, std::vector<Face>& nonoverlap_quadF) {
    std::set<size_t> nonoverlap_quadEids;
    for (auto& f : nonoverlap_quadF)
        nonoverlap_quadEids.insert(quadMesh.F[f.id].Eids.begin(), quadMesh.F[f.id].Eids.end());
    return nonoverlap_quadEids;
}

double get_avg_edge_length(const std::vector<Vertex>& new_quadV, const std::set<size_t>& nonoverlap_quadEids) {
    double avg_edge_length = 0;
    for (auto eid : nonoverlap_quadEids) {
        auto& e = quadMesh.E.at(eid);
        avg_edge_length += glm::length(new_quadV[e.Vids[0]].xyz() - new_quadV[e.Vids[1]].xyz());
    }
    avg_edge_length /= nonoverlap_quadEids.size();
    return avg_edge_length;
}

std::vector<size_t> get_bad_quadids(const std::vector<Vertex>& new_quadV, const std::vector<Face>& next_quadF, double avg_edge_length) {
    std::vector<size_t> bad_quadids;
    for (auto& f : next_quadF) {
        for (auto eid : quadMesh.F[f.id].Eids) {
            auto& e = quadMesh.E.at(eid);
            auto length = glm::length(new_quadV[e.Vids[0]].xyz() - new_quadV[e.Vids[1]].xyz());
            if (length > avg_edge_length * 3) {
                bad_quadids.push_back(f.id);
                break;
            }
        }
    }
    return bad_quadids;
}

std::vector<Face> get_nonoverlap_quadF(const std::vector<Face>& inside_quadF, const std::map<size_t, size_t>& overlap_quad_vids) {
    std::vector<Face> nonoverlap_quadF;
    for (auto& f : inside_quadF) {
        bool all_nonoverlap = true;
        for (auto vid : f.Vids) {
            if (overlap_quad_vids.find(vid) != overlap_quad_vids.end()) {
                all_nonoverlap = false;
                break;
            }
        }
        if (all_nonoverlap) nonoverlap_quadF.push_back(f);
    }
    return nonoverlap_quadF;
}

std::vector<Face> get_next_quadF_good(const std::vector<Face>& next_quadF, const std::vector<size_t>& bad_quadids) {
    std::vector<Face> next_quadF_good;
    for (auto& f : next_quadF) {
        bool isbad = false;
        for (auto fid : bad_quadids)
            if (f.id == fid) {
                isbad = true;
                break;
            }
        if (!isbad) next_quadF_good.push_back(f);
    }
    return next_quadF_good;
}

std::map<size_t, Parameters> get_nonoverlap_Vids(const std::vector<Vertex>& quadV, const std::vector<std::vector<Parameters>>& vertex_triangles) {
    std::map<size_t, Parameters> nonoverlap_Vids;
    for (auto& v : quadV)
        if (vertex_triangles[v.id].size() == 1) nonoverlap_Vids.insert(std::make_pair(v.id, vertex_triangles[v.id].front()));
    return nonoverlap_Vids;
}

std::map<size_t, Parameters> get_overlap_Vids(const std::vector<Vertex>& quadV, const std::vector<std::vector<Parameters>>& vertex_triangles,
        std::map<size_t, Parameters>& nonoverlap_Vids) {
    std::map<size_t, Parameters> overlap_Vids;
    for (auto& v : quadV)
        if (vertex_triangles[v.id].size() == 2) {
            auto& params = vertex_triangles[v.id];
            overlap_Vids.insert(std::make_pair(v.id, params[0].triangleid == nonoverlap_Vids[v.id].triangleid ? params[1] : params[0]));
        }
    return overlap_Vids;
}

std::map<size_t, Parameters> get_nonoverlapquad_Vids(const std::vector<Face>& nonoverlap_quadF, const std::vector<std::vector<Parameters>>& vertex_triangles) {
    std::map<size_t, Parameters> nonoverlapquad_Vids;
    for (auto& f : nonoverlap_quadF)
        for (auto vid : f.Vids)
            nonoverlapquad_Vids.insert(std::make_pair(vid, vertex_triangles[vid].front()));
    return nonoverlapquad_Vids;
}

std::map<size_t, Parameters> get_remaining_nonoverlap_Vids(const std::set<size_t>& nonoverlapquad_extendquad_vids,
        const std::map<size_t, Parameters>& nonoverlapquad_Vids, const std::vector<std::vector<Parameters>>& vertex_triangles) {
    std::map<size_t, Parameters> remaining_nonoverlap_Vids;
//for (auto item : nonoverlap_Vids) {
//    const auto vid = item.first;
//    if (nonoverlapquad_Vids.find(vid) == nonoverlapquad_Vids.end())
//        remaining_nonoverlap_Vids.insert(item);
//}
    for (auto vid : nonoverlapquad_extendquad_vids)
        if (nonoverlapquad_Vids.find(vid) == nonoverlapquad_Vids.end() && !vertex_triangles[vid].empty()) remaining_nonoverlap_Vids.insert(
                std::make_pair(vid, vertex_triangles[vid].front()));
    return remaining_nonoverlap_Vids;
}

std::map<size_t, Parameters> get_remaining_overlap_Vids(const std::set<size_t>& nonoverlapquad_extendquad_vids,
        const std::map<size_t, Parameters>& nonoverlapquad_Vids, const std::vector<std::vector<Parameters>>& vertex_triangles) {
    std::map<size_t, Parameters> remaining_nonoverlap_Vids;
    for (auto vid : nonoverlapquad_extendquad_vids)
        if (nonoverlapquad_Vids.find(vid) == nonoverlapquad_Vids.end()) remaining_nonoverlap_Vids.insert(
                std::make_pair(vid, Parameters()));
    return remaining_nonoverlap_Vids;
}

void update(const std::vector<Face>& next_quadF_good, std::map<size_t, Parameters>& all_nextvids,
        std::map<size_t, Parameters>& nonoverlap_Vids, std::map<size_t, Parameters>& nonoverlapquad_Vids, std::vector<Face>& nonoverlap_quadF) {
    for (auto& f : next_quadF_good) {
        for (auto& vid : f.Vids) {
            if (nonoverlap_Vids.find(vid) == nonoverlap_Vids.end()) {
                nonoverlap_Vids[vid] = all_nextvids[vid];
                nonoverlapquad_Vids[vid] = all_nextvids[vid];
            }
        }
        nonoverlap_quadF.push_back(f);
    }
}

void bfs (std::vector<VertexUV>& V, std::vector<Face>& F,
        std::vector<Vertex>& quadV, std::vector<Vertex>& new_quadV, const std::vector<Face>& quadF,
        std::vector<Face>& inside_quadF, std::vector<Face>& nonoverlap_quadF,
        std::vector<std::vector<Parameters>>& vertex_triangles, double avg_edge_length,
        std::set<size_t>& inside_quadF_set, std::set<size_t>& nonoverlapquad_extendquadids, std::set<size_t>& nonoverlapquad_extendquad_vids,
        std::map<size_t, Parameters>& nonoverlap_Vids, std::map<size_t, Parameters>& nonoverlapquad_Vids, std::map<size_t, Parameters>& remaining_nonoverlap_Vids,
        int max_iters = 20) {
//    std::set<size_t> nonoverlap_quadEids = get_nonoverlap_quadEids(quadMesh, nonoverlap_quadF);
//    double avg_edge_length = get_avg_edge_length(new_quadV, nonoverlap_quadEids);
//    std::cout << "avg_edge_length = " << avg_edge_length << "\n\n";
//    //----------------------------------------------------
//    std::set<size_t> inside_quadF_set;
//    for (auto& f : inside_quadF)
//        inside_quadF_set.insert(f.id);
//    std::set<size_t> nonoverlapquad_extendquadids;
//    std::set<size_t> nonoverlapquad_extendquad_vids;
//    get_extend(quadF, nonoverlap_quadF, inside_quadF_set, nonoverlapquad_extendquadids, nonoverlapquad_extendquad_vids);
//
//    std::map<size_t, Parameters> nonoverlap_Vids = get_nonoverlap_Vids(quadV, vertex_triangles);
//    std::map<size_t, Parameters> nonoverlapquad_Vids = get_nonoverlapquad_Vids(nonoverlap_quadF, vertex_triangles);
//    std::map<size_t, Parameters> remaining_nonoverlap_Vids = get_remaining_nonoverlap_Vids(nonoverlapquad_extendquad_vids, nonoverlapquad_Vids, vertex_triangles);
    int iter = 0;
    while (remaining_nonoverlap_Vids.size() && iter < max_iters) {
        std::cout << "iter = " << iter << "\n";
        std::cout << "nonoverlap_Vids.size() = " << nonoverlap_Vids.size() << "\n";
        std::cout << "nonoverlapquad_Vids.size() = " << nonoverlapquad_Vids.size() << "\n";
        std::cout << "remaining_nonoverlap_Vids.size() = " << remaining_nonoverlap_Vids.size() << "\n";
        std::string fname = std::string("uv.nextvids") + std::to_string(iter) + ".vtk";
        MeshFileWriter nextwriter(quadMesh, fname.c_str());
        std::vector<size_t> xx;
        for (auto& item : remaining_nonoverlap_Vids)
            xx.push_back(item.first);
        //nextwriter.WriteVerticesVtk(xx);

        std::map<size_t, Parameters> all_nextvids;
        std::set<size_t> nextvids;
        std::set<size_t> nextquads;
        get_next(nonoverlap_Vids, nonoverlapquad_Vids, remaining_nonoverlap_Vids, vertex_triangles, inside_quadF_set, all_nextvids, nextvids, nextquads);
        std::cout << "all_nextvids.size() = " << all_nextvids.size() << "\n";
        std::cout << "nextvids.size() = " << nextvids.size() << "\n";
        std::cout << "nextquads.size() = " << nextquads.size() << "\n\n";

        std::vector<Face> next_quadF;
        for (auto nextquad : nextquads)
            next_quadF.push_back(quadF.at(nextquad));
        fname = std::string("uv.nextquads") + std::to_string(iter) + ".vtk";
        //write(fname.c_str(), quadV, next_quadF);

        // nonoverlap_Vids.insert(all_nextvids.begin(), all_nextvids.end());
        MapbackToOrigTri(V, F, all_nextvids, new_quadV);
        std::vector<size_t> bad_quadids = get_bad_quadids(new_quadV, next_quadF, avg_edge_length);
        std::cout << "bad_quadids.size() = " << bad_quadids.size() << "\n\n";

        std::vector<Face> next_quadF_good = get_next_quadF_good(next_quadF, bad_quadids);
        if (iter < 10) fname = std::string("next0") + std::to_string(iter) + ".quad.vtk";
        else fname = std::string("next") + std::to_string(iter) + ".quad.vtk";
        write(fname.c_str(), new_quadV, next_quadF_good);
        update(next_quadF_good, all_nextvids, nonoverlap_Vids, nonoverlapquad_Vids, nonoverlap_quadF);
        get_extend(quadF, nonoverlap_quadF, inside_quadF_set, nonoverlapquad_extendquadids, nonoverlapquad_extendquad_vids);
        remaining_nonoverlap_Vids = get_remaining_nonoverlap_Vids(nonoverlapquad_extendquad_vids, nonoverlapquad_Vids, vertex_triangles);
        ++iter;
    }
}
std::vector<Face> get_alloverlap_quadF(std::vector<VertexUV>& V, std::vector<Face>& F,
        std::vector<Vertex>& quadV, std::vector<Vertex>& new_quadV, std::vector<Face>& quadF,
        const std::vector<Face>& inside_quadF, std::map<size_t, size_t>& overlap_quad_vids,
        std::vector<std::vector<Parameters>>& vertex_triangles) {
    std::vector<Face> alloverlap_quadF;
    for (auto& f : inside_quadF) {
        bool all_overlap = true;
        std::set<size_t> triangleids;
        //std::map<size_t, glm::dvec2> vId_triangleid_uv;
        for (auto vid : f.Vids) {
            if (overlap_quad_vids.find(vid) == overlap_quad_vids.end()) {
                all_overlap = false;
            } else {
                for (auto& param : vertex_triangles[vid]) {
                    triangleids.insert(param.triangleid);
    //                    size_t key = (vid << 16) + param.triangleid;
    //                    vId_triangleid_uv[key] = param.uv;
                }
            }
        }
        if (all_overlap && triangleids.size() == 2) {
            alloverlap_quadF.push_back(f);
            //MapbackToOrigTri(V, F, quadV, quadF, vertex_triangles, f, new_quadV);
        }
    //        if (!all_overlap) {
    //            alloverlap_quadF.push_back(f);
    //            MapbackToOrigTri2(V, F, quadV, quadF, vertex_triangles, f, new_quadV);
    //        }
    }
    return alloverlap_quadF;
}
int main1(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: ExtractQuadMesh vertex_uv.m output.vtk\n";
        return -1;
    }
    const char* filename = argv[1];
    std::vector<VertexUV> V;
    std::vector<Face> F;
    read(filename, V, F);
    write(argv[2], V, F);

    generate_tri_mesh(V, F, triMesh);

    std::vector<Vertex> triV(V.size());
    for (auto& v : V)
        triV[v.id] = glm::dvec3(v.uv.x, v.uv.y, 0.0f);
    write("uv.tri.vtk", triV, F);

    std::vector<Vertex> quadV;
    std::vector<Face> quadF;
    generate_quad_mesh(quadV, quadF);
    generate_quad_mesh(quadV, quadF, quadMesh);

    std::vector<bool> vertex_inside(quadV.size(), false);
    std::vector<std::vector<Parameters>> vertex_triangles(quadV.size());
    std::map<size_t, size_t> overlap_quad_vids;
    get_inside_quads(V, F, quadV, vertex_triangles, overlap_quad_vids, vertex_inside);

    std::vector<Face> inside_quadF;
    std::vector<bool> quad_inside(quadF.size(), false);
    get_inside_quads(quadF, vertex_inside, inside_quadF, quad_inside);
    write("uv.quad.vtk", quadV, inside_quadF);

    std::set<size_t> inside_quadF_set;
    for (auto& f : inside_quadF)
        inside_quadF_set.insert(f.id);

    auto new_quadV = quadV;
    MapbackToOrigTri(V, F, quadV, quadF, vertex_triangles, new_quadV);
    write("inside.quad.vtk", new_quadV, inside_quadF);

    std::vector<Face> nonoverlap_quadF = get_nonoverlap_quadF(inside_quadF, overlap_quad_vids);
    write("nonoverlap.quad.vtk", new_quadV, nonoverlap_quadF);
    write("uv.nonoverlap.quad.vtk", quadV, nonoverlap_quadF);
    double avg_edge_length = 0.0;

        //----------------------------------------------------
        std::set<size_t> nonoverlap_quadEids = get_nonoverlap_quadEids(quadMesh, nonoverlap_quadF);
        avg_edge_length = get_avg_edge_length(new_quadV, nonoverlap_quadEids);
        std::cout << "avg_edge_length = " << avg_edge_length << "\n\n";
        std::set<size_t> nonoverlapquad_extendquadids;
        std::set<size_t> nonoverlapquad_extendquad_vids;
        get_extend(quadF, nonoverlap_quadF, inside_quadF_set, nonoverlapquad_extendquadids, nonoverlapquad_extendquad_vids);

        std::map<size_t, Parameters> nonoverlap_Vids = get_nonoverlap_Vids(quadV, vertex_triangles);
        std::map<size_t, Parameters> nonoverlapquad_Vids = get_nonoverlapquad_Vids(nonoverlap_quadF, vertex_triangles);
        std::map<size_t, Parameters> remaining_nonoverlap_Vids = get_remaining_nonoverlap_Vids(nonoverlapquad_extendquad_vids, nonoverlapquad_Vids, vertex_triangles);
        std::cout << "-->nonoverlap_Vids.size() = " << nonoverlap_Vids.size() << "\n";
        std::cout << "-->nonoverlapquad_Vids.size() = " << nonoverlapquad_Vids.size() << "\n";
        bfs(V, F, quadV, new_quadV, quadF, inside_quadF, nonoverlap_quadF, vertex_triangles,
                avg_edge_length, inside_quadF_set, nonoverlapquad_extendquadids, nonoverlapquad_extendquad_vids,
                nonoverlap_Vids, nonoverlapquad_Vids, remaining_nonoverlap_Vids);
        //----------------------------------------------------

    std::vector<Face> overlap_quadF = get_alloverlap_quadF(V, F, quadV, new_quadV, quadF, inside_quadF, overlap_quad_vids, vertex_triangles);
    std::cout << "overlap_quadF.size() = " << overlap_quadF.size() << "\n";
    //std::vector<Face> overlap_quadF;
    std::cout << "new_quadV.size() = " << new_quadV.size() << "\n";
    std::cout << "quadF.size() = " << quadF.size() << "\n";
//    for (int i = sqsize * sqsize; i < quadF.size(); ++i) {
//        overlap_quadF.push_back(quadF.at(i));
//    }
    write("overlap.quad.vtk", new_quadV, overlap_quadF);
    write("uv.overlap.quad.vtk", quadV, overlap_quadF);

    auto fff = overlap_quadF.at(89);
    overlap_quadF.clear();
    overlap_quadF.push_back(fff);
        //----------------------------------------------------
        std::set<size_t> overlap_quadEids = get_nonoverlap_quadEids(quadMesh, overlap_quadF);
        std::set<size_t> overlapquad_extendquadids;
        std::set<size_t> overlapquad_extendquad_vids;
        get_extend(quadF, overlap_quadF, inside_quadF_set, overlapquad_extendquadids, overlapquad_extendquad_vids);
        std::cout << "overlapquad_extendquadids.size() = " << overlapquad_extendquadids.size() << "\n";
        {
            std::vector<Face> extend_overlap_quadF;
            for (auto id : overlapquad_extendquadids)
                extend_overlap_quadF.push_back(quadF.at(id));
            write("uv.extend.quad.vtk", quadV, extend_overlap_quadF);
        }

        std::map<size_t, Parameters> overlap_Vids = get_overlap_Vids(quadV, vertex_triangles, nonoverlap_Vids);
        std::map<size_t, Parameters> overlapquad_Vids;
        std::set<size_t> overlap_quadF_Vids_set;
        for (auto& f : overlap_quadF)
            overlap_quadF_Vids_set.insert(f.Vids.begin(), f.Vids.end());
        for (auto& vid : overlap_quadF_Vids_set) {
            auto& params = vertex_triangles[vid];
            overlapquad_Vids.insert(std::make_pair(vid, params[0].triangleid == nonoverlap_Vids[vid].triangleid ? params[1] : params[0]));
        }
        // overlap_Vids = overlapquad_Vids;
        std::map<size_t, Parameters> remaining_overlap_Vids = get_remaining_overlap_Vids(overlapquad_extendquad_vids, overlapquad_Vids, vertex_triangles);
        bfs(V, F, quadV, new_quadV, quadF, inside_quadF, overlap_quadF, vertex_triangles,
                avg_edge_length, inside_quadF_set, overlapquad_extendquadids, overlapquad_extendquad_vids,
                overlap_Vids, overlapquad_Vids, remaining_overlap_Vids,100);
        //----------------------------------------------------

    return 0;
}

std::vector<size_t> get_triangle_ids(const std::vector<VertexUV>& triV, const std::vector<Face>& triF,
        const std::vector<Vertex>& quadV, const std::vector<Face>& quadF, const Face& f) {
    std::vector<size_t> res;

    const Vertex& v0 = quadV.at(f.Vids[0]);
    const Vertex& v2 = quadV.at(f.Vids[2]);
    Vertex center = 0.5 * (v0.xyz() + v2.xyz());
    for (auto& tri : triF) {
        glm::dvec2 uv;
        if (IsVertextInTriangle(center, triV[tri.Vids[0]], triV[tri.Vids[1]], triV[tri.Vids[2]], uv))
            res.push_back(tri.id);
    }

    return res;
}

std::vector<size_t> get_inside_quadids(const std::vector<VertexUV>& triV, const std::vector<Face>& triF,
        const std::vector<Vertex>& quadV, const std::vector<Face>& quadF) {
    std::vector<size_t> res;
    for (auto& f : quadF) {
        const Vertex& v0 = quadV.at(f.Vids[0]);
        const Vertex& v2 = quadV.at(f.Vids[2]);
        Vertex center = 0.5 * (v0.xyz() + v2.xyz());
        for (auto& tri : triF) {
            glm::dvec2 uv;
            if (IsVertextInTriangle(center, triV[tri.Vids[0]], triV[tri.Vids[1]], triV[tri.Vids[2]], uv)) {
                res.push_back(f.id);
                break;
            }
        }
    }
    return res;
}

void MapbackToOrigTri(const std::vector<VertexUV>& triV, const std::vector<Face>& triF,
        const std::vector<Vertex>& quadV, const std::vector<Face>& quadF,
        std::vector<Vertex>& new_quadV, std::vector<Face>& new_quadF) {

    std::vector<bool> vertex_inside(quadV.size(), false);
    std::vector<std::vector<Parameters>> vertex_triangles(quadV.size());
    std::map<size_t, size_t> overlap_quad_vids;
    get_inside_quads(triV, triF, quadV, vertex_triangles, overlap_quad_vids, vertex_inside);

    std::vector<Face> inside_quadF;
    std::vector<bool> quad_inside(quadF.size(), false);
    get_inside_quads(quadF, vertex_inside, inside_quadF, quad_inside);
    write("uv.quad.vtk", quadV, inside_quadF);

    std::set<size_t> inside_quadF_set;
    for (auto& f : inside_quadF)
        inside_quadF_set.insert(f.id);

    auto& p = vertex_triangles;
    std::cout << "MapbackToOrigTri" << std::endl;
    size_t new_vid = 0;
    size_t new_fid = 0;
    for (const auto& f : inside_quadF) {
        auto triangle_ids = get_triangle_ids(triV, triF, quadV, quadF, f);
        for (auto triangle_id : triangle_ids) {
            bool is_good_quad = true;
            for (auto vid : f.Vids) {
                auto& params = p[vid];
                if (params.empty()) {
                    std::cerr << "ERROR!\n";
                }
                glm::dvec2 uv;
                size_t triangleid;
                auto tri_face_ids = get_triangle_neighbor_triangleids(triMesh, triangle_id);
                bool is_good_quadvid = false;
                for (auto tri_face_id : tri_face_ids) {
                    bool found = false;
                    for (auto& param: params) {
                        if (param.triangleid == tri_face_id) {
                            uv = param.uv;
                            triangleid = param.triangleid;
                            found = true;
                            is_good_quadvid = true;
                            break;
                        }
                    }
                    if (found) break;
                }
                if (!is_good_quadvid) {
                    is_good_quad = false;
                    break;
                }
                const auto& triVids = triF.at(triangleid).Vids;

                const Vertex& origT0 = triV[triVids.at(0)];
                const Vertex& origT1 = triV[triVids.at(1)];
                const Vertex& origT2 = triV[triVids.at(2)];
                const double rambda0 = uv.x;
                const double rambda1 = uv.y;
                const double rambda2 = 1.0 - rambda0 - rambda1;

                Vertex v;
                v.x = rambda0 * origT0.x + rambda1 * origT1.x + rambda2 * origT2.x;
                v.y = rambda0 * origT0.y + rambda1 * origT1.y + rambda2 * origT2.y;
                v.z = rambda0 * origT0.z + rambda1 * origT1.z + rambda2 * origT2.z;
                v.id = new_vid++;
                new_quadV.push_back(v);
            }
            if (!is_good_quad) continue;
            Face new_face;
            new_face.Vids = {new_vid - 4, new_vid - 3, new_vid - 2, new_vid - 1};
            new_face.id = new_fid++;
            new_quadF.push_back(new_face);
        }
    }
}

void compress(std::vector<Vertex>& V, std::vector<Face>& F, const double eps = 1e-3) {
    std::unordered_map<size_t, size_t> vid_vid;
#pragma omp parallel for
    for (int i = V.size() - 1; i >= 0; --i) {
        for (int j = 0; j <= i; ++j) {
            float length = glm::length(V[i].xyz() - V[j].xyz());
            if (length < eps) {
#pragma omp critical
                vid_vid[i] = j;
                break;
            }
        }
    }

    std::set<size_t> vid_set;
    for (auto& f : F)
        for (auto& vid : f.Vids) {
            vid = vid_vid[vid];
            if (vid_set.find(vid) == vid_set.end()) vid_set.insert(vid);
        }

    std::unordered_map<std::string, std::vector<size_t>> face_vids;
#pragma omp parallel for
    for (size_t i = 0; i < F.size(); i++) {
        auto& f = F[i];
        std::set<size_t> vids(f.Vids.begin(), f.Vids.end());
        std::string key;
        for (auto vid : vids)
            key += std::to_string(vid) + "@";
#pragma omp critical
        face_vids[key] = f.Vids;
    }
    std::vector<Face> newF(face_vids.size());
    size_t fid = 0;
    for (auto& item : face_vids)
        newF[fid++] = item.second;
    F = newF;

    std::vector<size_t> v_real_index(vid_set.begin(), vid_set.end());
    std::vector<Vertex> newV(v_real_index.size());
//#pragma omp parallel for
    for (unsigned long i = 0; i < v_real_index.size(); i++) {
        const auto& v = V.at(v_real_index.at(i));
        newV.at(i) = v;
    }
    V = newV;
    //////////////////////////////////////////////////////
    std::unordered_map<size_t, size_t> v_v;
    for (size_t i = 0; i < v_real_index.size(); i++)
        v_v[v_real_index.at(i)] = i;
    for (auto& f : F)
        for (auto& vid : f.Vids)
            vid = v_v[vid];
}

void clean(std::vector<Vertex>& new_quadV, std::vector<Face>& new_quadF) {
    std::vector<Face> newnew_quadF;
    for (auto& f : new_quadF) {
        bool isbad = false;
        double minl = 1000000;
        double maxl = -1000000;
        for (auto i = 0; i < 4; ++i) {
            const auto& v0 = new_quadV[f.Vids[i]];
            const auto& v1 = new_quadV[f.Vids[(i + 1) % 4]];
            auto length = glm::length(v0.xyz() - v1.xyz());
            minl = std::min(length, minl);
            maxl = std::max(length, maxl);
        }
        if (maxl < 2* minl) newnew_quadF.push_back(f);
    }
    new_quadF = newnew_quadF;

    std::set<size_t> nonoverlap_quadEids = get_nonoverlap_quadEids(quadMesh, new_quadF);
    auto avg_edge_length = get_avg_edge_length(new_quadV, nonoverlap_quadEids);
    std::cout << "avg_edge_length = " << avg_edge_length << "\n\n";

    compress(new_quadV, new_quadF, avg_edge_length/1000);
}


int main2(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: ExtractQuadMesh vertex_uv.m output.vtk\n";
        return -1;
    }
    if (argc > 3)
    sqsize = std::stoi(argv[3]);
    const char* filename = argv[1];
    std::vector<VertexUV> V;
    std::vector<Face> F;
    read(filename, V, F);
    write(argv[2], V, F);

    generate_tri_mesh(V, F, triMesh);
    {
        std::vector<size_t> edgeids;
        //std::vector<size_t> edgeids_set;
        std::set<size_t> boundary_vids;
        for (auto& e : triMesh.E)
            if (e.N_Fids.size() == 1) {
                edgeids.push_back(e.id);
                boundary_vids.insert(e.Vids[0]);
                boundary_vids.insert(e.Vids[1]);
            }
        {
            MeshFileWriter writer(triMesh, "tri.boundary.vtk");
            writer.WriteEdgesVtk(edgeids);
//            std::ofstream ofs("boundary_vids");
//            for (auto vid : boundary_vids)
//                ofs << vid << " ";
//            ofs << "\n";
        }
        std::vector<size_t> singular_vids;
        for (auto vid : boundary_vids) {
            auto& v = triMesh.V.at(vid);
            std::vector<size_t> nvids;
            for (auto nvid : v.N_Vids)
                if (boundary_vids.find(nvid) != boundary_vids.end()) nvids.push_back(nvid);
            auto& v1 = triMesh.V.at(nvids[0]);
            auto& v2 = triMesh.V.at(nvids[1]);
            if (nvids.size() == 2) {
                auto length = glm::length(v1.xyz() - v2.xyz());
                if (length < 1e-3) singular_vids.push_back(vid);
            } else if (nvids.size() == 3) {
                auto& v3 = triMesh.V.at(nvids[2]);
                auto length1 = glm::length(v1.xyz() - v2.xyz());
                auto length2 = glm::length(v1.xyz() - v3.xyz());
                auto length3 = glm::length(v2.xyz() - v3.xyz());
                if (length1 < 1e-3 || length2 < 1e-3 || length3 < 1e-3) singular_vids.push_back(vid);
            }

//            if (vid == 542) {
//                std::cerr << length << "\n";
//                for (auto nvid : v.N_Vids)
//                    std::cerr << nvid << " ";
//                std::cerr << "\n---------------\n";
//                for (auto nvid : nvids)
//                    std::cerr << nvid << " ";
//                std::cerr << "\n";
//            }
        }
        {
            MeshFileWriter writer(triMesh, "singularities.quad.vtk");
            writer.WriteVerticesVtk(singular_vids);
        }
    }

    std::vector<Vertex> triV(V.size());
    for (auto& v : V)
        triV[v.id] = glm::dvec3(v.uv.x, v.uv.y, 0.0f);
    write("uv.tri.vtk", triV, F);
    write("tri.vtk", V, F);

    double minx = 1000000, miny = 1000000, maxx = -1000000, maxy = -1000000;
    for (auto& v : V) {
        if (v.uv.x < minx) minx = v.uv.x;
        if (v.uv.y < miny) miny = v.uv.y;
        if (v.uv.x > maxx) maxx = v.uv.x;
        if (v.uv.y > maxy) maxy = v.uv.y;
    }
    std::vector<Vertex> quadV;
    std::vector<Face> quadF;
    if (argc > 4) generate_quad_mesh2(quadV, quadF);
    else generate_quad_mesh(quadV, quadF);
    // generate_quad_mesh(quadV, quadF, minx, miny, maxx, maxy);
    generate_quad_mesh(quadV, quadF, quadMesh);


    std::vector<Vertex> new_quadV;
    std::vector<Face> new_quadF;
    MapbackToOrigTri(V, F, quadV, quadF, new_quadV, new_quadF);
    clean(new_quadV, new_quadF);
    write(argv[2], new_quadV, new_quadF);

    {
        MeshFileReader reader(argv[2]);
        Mesh& mesh = (Mesh&)reader.GetMesh();
        mesh.RemoveUselessVertices();
        mesh.BuildAllConnectivities();
        mesh.ExtractBoundary();

        MeshFileWriter writer(mesh, "quad.boundary.vtk");
        std::vector<size_t> edgeids;
        for (auto& e : mesh.E)
            if (e.N_Fids.size() == 1) edgeids.push_back(e.id);
        writer.WriteEdgesVtk(edgeids);
    }
    return 0;
}


std::set<size_t> uvTriMeshSingularCutEdgeIds;
std::map<size_t, size_t> uvEid_origEid;
std::set<size_t> cutEdgeIds;
std::set<size_t> cutVertexIds;
std::vector<size_t> grid_vids;
std::vector<size_t> get_uvTriMeshCutEdgeIds(const Mesh& origTriMesh, const Mesh& uvTriMesh, const std::vector<VertexUV>& V) {
    std::vector<size_t> res;
    for (auto& e : uvTriMesh.E) {
        auto& v0 = uvTriMesh.V.at(e.Vids[0]);
        auto& v1 = uvTriMesh.V.at(e.Vids[1]);
        if (V[v0.id].father >= origTriMesh.V.size() || V[v1.id].father >= origTriMesh.V.size()) continue;
        auto& orig_v0 = origTriMesh.V.at(V[v0.id].father);
        auto& orig_v1 = origTriMesh.V.at(V[v1.id].father);

        std::vector<size_t> orig_v0_neids = orig_v0.N_Eids;
        std::vector<size_t> orig_v1_neids = orig_v1.N_Eids;
        std::sort(orig_v0_neids.begin(), orig_v0_neids.end());
        std::sort(orig_v1_neids.begin(), orig_v1_neids.end());

        std::vector<size_t> eids;
        std::set_intersection(orig_v0_neids.begin(), orig_v0_neids.end(), orig_v1_neids.begin(), orig_v1_neids.end(), std::back_inserter(eids));
        if (eids.size() != 1) {
            std::cerr << "Error in get_uvTriMeshSingularCutEdgeIds()\n";
        } else {
            const auto orig_eid = eids.front();
            const auto& orig_e = origTriMesh.E.at(orig_eid);
            uvEid_origEid[e.id] = orig_eid;
            if (e.isBoundary && !orig_e.isBoundary) {
                res.push_back(e.id);
                cutVertexIds.insert(e.Vids.begin(), e.Vids.end());
            }
        }
    }
    cutEdgeIds.insert(res.begin(), res.end());
    return res;
}

std::vector<size_t> get_uvTriMeshSingularCutEdgeIds(const Mesh& origTriMesh, const Mesh& uvTriMesh, const std::vector<VertexUV>& V) {
    std::vector<size_t> res;
    auto remain = cutEdgeIds;
    std::vector<std::set<size_t>> cut_lines;
    while (!remain.empty()) {
        std::set<size_t> cut_line_eids;
        std::queue<size_t> q;
        q.push(*remain.begin());
        while (!q.empty()) {
            auto n = q.size();
            for (auto i = 0; i < n; ++i) {
                auto eid = q.front();
                q.pop();
                auto& e = uvTriMesh.E.at(eid);
                for (auto vid : e.Vids) {
                    auto& v = uvTriMesh.V.at(vid);
                    for (auto neid : v.N_Eids) {
                        auto& ne = uvTriMesh.E.at(neid);
                        if (!ne.isBoundary) continue;
                        if (cutEdgeIds.find(neid) != cutEdgeIds.end() && cut_line_eids.find(neid) == cut_line_eids.end()) {
                            cut_line_eids.insert(neid);
                            q.push(neid);
                        }
                    }
                }
            }
        }
        for (auto eid : cut_line_eids) remain.erase(eid);
        cut_lines.push_back(cut_line_eids);
    }

    std::cout << "cut_lines.size() = " << cut_lines.size() << "\n";
    auto uv_boundary_vids = get_boundary_vids(uvTriMesh);
    auto singular_vids = get_singular_vids(uvTriMesh, uv_boundary_vids);
    std::set<size_t> singular_vids_set(singular_vids.begin(), singular_vids.end());
    for (auto& cut_line : cut_lines) {
        bool foundSingularV = false;
        for (auto eid : cut_line) {
            auto& e = uvTriMesh.E.at(eid);
            if (singular_vids_set.find(e.Vids[0]) != singular_vids_set.end() || singular_vids_set.find(e.Vids[1]) != singular_vids_set.end()) {
                foundSingularV = true;
                break;
            }
        }
        if (foundSingularV) std::copy(cut_line.begin(), cut_line.end(), std::back_inserter(res));
    }
    uvTriMeshSingularCutEdgeIds.insert(res.begin(), res.end());
    return res;
}


std::map<size_t, Parameters> endingTriangleid_neighborParams;

std::set<size_t> get_boundary_vids(const Mesh& uvtriMesh) {
    std::set<size_t> boundary_vids;
    for (auto& e : uvtriMesh.E)
        if (e.N_Fids.size() == 1) {
            boundary_vids.insert(e.Vids[0]);
            boundary_vids.insert(e.Vids[1]);
        }
    return boundary_vids;
}

std::set<size_t> get_boundary_eids(const Mesh& uvtriMesh) {
    std::set<size_t> edgeids;
    for (auto& e : uvtriMesh.E)
        if (e.N_Fids.size() == 1)
            edgeids.insert(e.id);
    return edgeids;
}

std::vector<size_t> get_singular_vids(const Mesh& uvtriMesh, const std::set<size_t>& boundary_vids) {
    std::vector<size_t> singular_vids;
    for (auto vid : boundary_vids) {
        auto& v = uvtriMesh.V.at(vid);
        std::vector<size_t> nvids;
        for (auto nvid : v.N_Vids)
            if (boundary_vids.find(nvid) != boundary_vids.end()) nvids.push_back(nvid);
        auto& v1 = uvtriMesh.V.at(nvids[0]);
        auto& v2 = uvtriMesh.V.at(nvids[1]);
        if (nvids.size() == 2) {
            auto length = glm::length(v1.xyz() - v2.xyz());
            if (length < 1e-3) singular_vids.push_back(vid);
        } else if (nvids.size() == 3) {
            auto& v3 = uvtriMesh.V.at(nvids[2]);
            auto length1 = glm::length(v1.xyz() - v2.xyz());
            auto length2 = glm::length(v1.xyz() - v3.xyz());
            auto length3 = glm::length(v2.xyz() - v3.xyz());
            if (length1 < 1e-3 || length2 < 1e-3 || length3 < 1e-3) singular_vids.push_back(vid);
        }
    }

    return singular_vids;
}

double get_step_size(const std::vector<VertexUV>& V) {
    double minx = 1000000, miny = 1000000, maxx = -1000000, maxy = -1000000;
    for (auto& v : V) {
        if (v.uv.x < minx) minx = v.uv.x;
        if (v.uv.y < miny) miny = v.uv.y;
        if (v.uv.x > maxx) maxx = v.uv.x;
        if (v.uv.y > maxy) maxy = v.uv.y;
    }
    double lenx = maxx - minx;
    double leny = maxy - miny;
    double len = (lenx > leny) ? leny : lenx;
    double d = len / sqsize;
    return d;
}

bool get_parameter(const Mesh& uvTriMesh, const glm::dvec2& p, const Face& tri, Parameters& res) {
    bool success = true;
    Vertex v = glm::dvec3(p.x, p.y, 0.0);
    glm::dvec2 uv;
    if (IsVertextInTriangle(v, uvTriMesh.V[tri.Vids[0]], uvTriMesh.V[tri.Vids[1]], uvTriMesh.V[tri.Vids[2]], uv)) res = Parameters(uv, tri.id);
    else {
        std::cerr << "Error in get_parameters triangle_id " << tri.id << "\n";
        success = false;
    }
    return success;
}

glm::dvec3 get_origin_coordinate(const Mesh& origTriMesh, const glm::dvec2& p, const Parameters& param) {
    const auto& triF = origTriMesh.F.at(param.triangleid);
    const auto& triVids = triF.Vids;

    const Vertex& origT0 = origTriMesh.V[triVids.at(0)];
    const Vertex& origT1 = origTriMesh.V[triVids.at(1)];
    const Vertex& origT2 = origTriMesh.V[triVids.at(2)];
    const double rambda0 = param.uv.x;
    const double rambda1 = param.uv.y;
    const double rambda2 = 1.0 - rambda0 - rambda1;

    glm::dvec3 v;
    v.x = rambda0 * origT0.x + rambda1 * origT1.x + rambda2 * origT2.x;
    v.y = rambda0 * origT0.y + rambda1 * origT1.y + rambda2 * origT2.y;
    v.z = rambda0 * origT0.z + rambda1 * origT1.z + rambda2 * origT2.z;
    return v;
}

std::vector<Parameters> get_parameters(const Mesh& uvTriMesh, const std::vector<VertexUV>& V, const glm::dvec2& p,
        const std::set<size_t>& constraint_triangle_ids) {
    std::vector<Parameters> res;
    Vertex v = glm::dvec3(p.x, p.y, 0.0);
    for (auto& tri_id : constraint_triangle_ids) {
        auto& tri = uvTriMesh.F.at(tri_id);
        glm::dvec2 uv;
        if (IsVertextInTriangle(v, V[tri.Vids[0]], V[tri.Vids[1]], V[tri.Vids[2]], uv))
            res.push_back(Parameters(uv, tri.id));
    }

    return res;
}

bool get_geodesic_vertex(const Mesh& origTriMesh, const Mesh& uvTriMesh, const std::vector<VertexUV>& V,
        const glm::dvec2& p, const std::set<size_t>& constraint_triangle_ids, Parameters& out_param, Vertex& out_v) {
    bool res = false;
    std::set<size_t> triangle_ids = {out_param.triangleid};
    {
        auto params = get_parameters(uvTriMesh, V, p, triangle_ids);
        for (auto& param : params) {
            if (triangle_ids.find(param.triangleid) != triangle_ids.end()) {
                out_param = param;
                res = true;
                out_v = get_origin_coordinate(uvTriMesh, p, param);
                break;
            }
        }
        if (res)
            return res;
    }
    auto params = get_parameters(uvTriMesh, V, p, constraint_triangle_ids);
    for (auto& param : params) {
        if (constraint_triangle_ids.find(param.triangleid) != constraint_triangle_ids.end()) {
            out_param = param;
            res = true;
            out_v = get_origin_coordinate(uvTriMesh, p, param);
            break;
        }
    }
    return res;
}

bool get_geodesic_vertex(const Mesh& origTriMesh, const Mesh& uvTriMesh, const std::vector<VertexUV>& V,
    const glm::dvec2& p, const std::set<size_t>& constraint_triangle_ids, Parameters& out_param, MapPoint& out_v) {
    bool res = false;
    std::set<size_t> triangle_ids = { out_param.triangleid };
    {
        auto params = get_parameters(uvTriMesh, V, p, triangle_ids);
        for (auto& param : params) {
            if (triangle_ids.find(param.triangleid) != triangle_ids.end()) {
                out_param = param;
                res = true;
                out_v.xyz = get_origin_coordinate(uvTriMesh, p, param);
                out_v.uv = out_param.uv;
                out_v.tri_id = out_param.triangleid;
                break;
            }
        }
        if (res)
            return res;
    }

    auto params = get_parameters(uvTriMesh, V, p, constraint_triangle_ids);
    for (auto& param : params) {
        const auto& ending_tri = uvTriMesh.F.at(out_param.triangleid);
        const auto& v0 = V[ending_tri.Vids[0]];
        const auto& v1 = V[ending_tri.Vids[1]];
        const auto& v2 = V[ending_tri.Vids[2]];
        const auto& orig_v0 = origTriMesh.V[v0.father];
        const auto& orig_v1 = origTriMesh.V[v1.father];
        const auto& orig_v2 = origTriMesh.V[v2.father];
        //const vector<size_t>& orig_ending_tri_vids = {origTriMesh.V[v0.father]};

        //if (!orig_v0.isBoundary && !orig_v1.isBoundary && !orig_v2.isBoundary) 
        bool found_e0 = cutEdgeIds.find(ending_tri.Eids[0]) != cutEdgeIds.end();
        bool found_e1 = cutEdgeIds.find(ending_tri.Eids[1]) != cutEdgeIds.end();
        bool found_e2 = cutEdgeIds.find(ending_tri.Eids[2]) != cutEdgeIds.end();
        if (found_e0 || found_e1 || found_e2) {
            std::vector<size_t> pre_vids = uvTriMesh.F.at(out_param.triangleid).Vids;
            std::vector<size_t> cur_vids = uvTriMesh.F.at(param.triangleid).Vids;
            std::sort(pre_vids.begin(), pre_vids.end());
            std::sort(cur_vids.begin(), cur_vids.end());

            std::vector<size_t> vids;
            std::set_intersection(pre_vids.begin(), pre_vids.end(), cur_vids.begin(), cur_vids.end(), std::back_inserter(vids));
            if (vids.size() <= 1) {
                continue;
            }
        }
        if (constraint_triangle_ids.find(param.triangleid) != constraint_triangle_ids.end()) {
            out_param = param;
            res = true;
            out_v.xyz = get_origin_coordinate(uvTriMesh, p, param);
            out_v.uv = out_param.uv;
            out_v.tri_id = out_param.triangleid;
            break;
        }
    }
    return res;
}


size_t get_opposite_edge_id(const Mesh& origTriMesh, const Mesh& uvTriMesh, const std::vector<VertexUV>& V, const VertexUV& v0, const VertexUV& v1) {
    if (uvTriMesh.V[v0.id].isSingularity) {
        for (auto& e: uvTriMesh.E) {
            const auto& v_0 = V[e.Vids[0]];
            const auto& v_1 = V[e.Vids[1]];
            if (e.Vids[0] == v0.id && e.Vids[1] != v1.id && glm::length(v_1.xyz() - v1.xyz()) < 1e-5) return e.id;
            if (e.Vids[1] == v0.id && e.Vids[0] != v1.id && glm::length(v_0.xyz() - v1.xyz()) < 1e-5) return e.id;
        }
    } else if (uvTriMesh.V[v1.id].isSingularity) {
        for (auto& e: uvTriMesh.E) {
            const auto& v_0 = V[e.Vids[0]];
            const auto& v_1 = V[e.Vids[1]];
            if (e.Vids[0] == v1.id && e.Vids[1] != v0.id && glm::length(v_1.xyz() - v0.xyz()) < 1e-5) return e.id;
            if (e.Vids[1] == v1.id && e.Vids[0] != v0.id && glm::length(v_0.xyz() - v0.xyz()) < 1e-5) return e.id;
        }
    }

//    for (auto& e: uvTriMesh.E) {
//        if ((V[e.Vids[0]].xyz() == v1.xyz() && V[e.Vids[1]].xyz() == v2.xyz()) ||
//                (V[e.Vids[0]].xyz() == v2.xyz() && V[e.Vids[1]].xyz() == v1.xyz())) {
//            if (e.Vids[0] != v1.id && e.Vids[0] != v2.id && e.Vids[1] != v1.id && e.Vids[1] != v2.id) {
//                return e.id;
//            }
//        }
//    }

    for (auto& e: uvTriMesh.E) {
        const auto& v_0 = V[e.Vids[0]];
        const auto& v_1 = V[e.Vids[1]];
        if ((glm::length(v_0.xyz() - v0.xyz()) < 1e-5 && glm::length(v_1.xyz() - v1.xyz()) < 1e-5) ||
                (glm::length(v_0.xyz() - v1.xyz()) < 1e-5 && glm::length(v_1.xyz() - v0.xyz()) < 1e-5)) {
            if (e.Vids[0] != v0.id && e.Vids[0] != v1.id && e.Vids[1] != v0.id && e.Vids[1] != v1.id) {
                return e.id;
            }
        }
    }
    return MAXID;
}

Parameters g_param;

bool is_on_cut_triangle(const Mesh& origTriMesh, const Mesh& uvTriMesh, const std::vector<VertexUV>& V, const Parameters& param,
        const glm::dvec2& p) {
    if (param.triangleid > uvTriMesh.F.size())
        return false;
    //std::cout << "ending triangle = " <<  param.triangleid << "\n";
    bool res = false;
    const auto& ending_tri = uvTriMesh.F.at(param.triangleid);
    const auto& v0 = V[ending_tri.Vids[0]];
    const auto& v1 = V[ending_tri.Vids[1]];
    const auto& v2 = V[ending_tri.Vids[2]];
    const auto& orig_v0 = origTriMesh.V[v0.father];
    const auto& orig_v1 = origTriMesh.V[v1.father];
    const auto& orig_v2 = origTriMesh.V[v2.father];
    //const vector<size_t>& orig_ending_tri_vids = {origTriMesh.V[v0.father]};

    //if (!orig_v0.isBoundary && !orig_v1.isBoundary && !orig_v2.isBoundary) 
    bool found_e0 = cutEdgeIds.find(ending_tri.Eids[0]) != cutEdgeIds.end();
    bool found_e1 = cutEdgeIds.find(ending_tri.Eids[1]) != cutEdgeIds.end();
    bool found_e2 = cutEdgeIds.find(ending_tri.Eids[2]) != cutEdgeIds.end();
    if (found_e0 || found_e1 || found_e2)
    {
        // std::cout << "triangle " << param.triangleid << " is on the cut. ";
        res = true;
        auto opposite_edge_id = MAXID;
        if (uvTriMesh.V[v0.id].isSingularity || uvTriMesh.V[v1.id].isSingularity || uvTriMesh.V[v2.id].isSingularity) {
            if (uvTriMesh.V[v0.id].isSingularity) {
                if (uvTriMesh.V[v1.id].isBoundary) opposite_edge_id = get_opposite_edge_id(origTriMesh, uvTriMesh, V, v0, v1);
                else  opposite_edge_id = get_opposite_edge_id(origTriMesh, uvTriMesh, V, v0, v2);
            } else if (uvTriMesh.V[v1.id].isSingularity) {
                if (uvTriMesh.V[v0.id].isBoundary) opposite_edge_id = get_opposite_edge_id(origTriMesh, uvTriMesh, V, v0, v1);
                else  opposite_edge_id = get_opposite_edge_id(origTriMesh, uvTriMesh, V, v1, v2);
            } else if (uvTriMesh.V[v2.id].isSingularity) {
                if (uvTriMesh.V[v0.id].isBoundary) opposite_edge_id = get_opposite_edge_id(origTriMesh, uvTriMesh, V, v0, v2);
                else  opposite_edge_id = get_opposite_edge_id(origTriMesh, uvTriMesh, V, v1, v2);
            }
        } else if ((uvTriMesh.V[v0.id].isBoundary && uvTriMesh.V[v1.id].isBoundary)
                || (uvTriMesh.V[v0.id].isBoundary && uvTriMesh.V[v2.id].isBoundary)
                || (uvTriMesh.V[v1.id].isBoundary && uvTriMesh.V[v2.id].isBoundary)) {
            //if (uvTriMesh.V[v0.id].isBoundary && uvTriMesh.V[v1.id].isBoundary && (!orig_v0.isBoundary || !orig_v1.isBoundary))
            //    opposite_edge_id = get_opposite_edge_id(origTriMesh, uvTriMesh, V, v0, v1);
            //else if (uvTriMesh.V[v0.id].isBoundary && uvTriMesh.V[v2.id].isBoundary && (!orig_v0.isBoundary || !orig_v2.isBoundary))
            //    opposite_edge_id = get_opposite_edge_id(origTriMesh, uvTriMesh, V, v0, v2);
            //else if (uvTriMesh.V[v1.id].isBoundary && uvTriMesh.V[v2.id].isBoundary && (!orig_v1.isBoundary || !orig_v2.isBoundary))
            //    opposite_edge_id = get_opposite_edge_id(origTriMesh, uvTriMesh, V, v1, v2);
            //if (uvTriMesh.V[v0.id].isBoundary && uvTriMesh.V[v1.id].isBoundary)
            //    opposite_edge_id = get_opposite_edge_id(origTriMesh, uvTriMesh, V, v0, v1);
            //else if (uvTriMesh.V[v0.id].isBoundary && uvTriMesh.V[v2.id].isBoundary)
            //    opposite_edge_id = get_opposite_edge_id(origTriMesh, uvTriMesh, V, v0, v2);
            //else if (uvTriMesh.V[v1.id].isBoundary && uvTriMesh.V[v2.id].isBoundary)
            //    opposite_edge_id = get_opposite_edge_id(origTriMesh, uvTriMesh, V, v1, v2);
            if (found_e0) {
                const auto& e = uvTriMesh.E[ending_tri.Eids[0]];
                opposite_edge_id = get_opposite_edge_id(origTriMesh, uvTriMesh, V, V[e.Vids[0]], V[e.Vids[1]]);
            } else if (found_e1) {
                const auto& e = uvTriMesh.E[ending_tri.Eids[1]];
                opposite_edge_id = get_opposite_edge_id(origTriMesh, uvTriMesh, V, V[e.Vids[0]], V[e.Vids[1]]);
            } else if (found_e2) {
                const auto& e = uvTriMesh.E[ending_tri.Eids[2]];
                opposite_edge_id = get_opposite_edge_id(origTriMesh, uvTriMesh, V, V[e.Vids[0]], V[e.Vids[1]]);
            }
        } else if (!orig_v0.isBoundary && !orig_v1.isBoundary) {
            opposite_edge_id = get_opposite_edge_id(origTriMesh, uvTriMesh, V, v0, v1);
        } else if (!orig_v0.isBoundary && !orig_v2.isBoundary) {
            opposite_edge_id = get_opposite_edge_id(origTriMesh, uvTriMesh, V, v0, v2);
        } else if (!orig_v1.isBoundary && !orig_v2.isBoundary) {
            opposite_edge_id = get_opposite_edge_id(origTriMesh, uvTriMesh, V, v1, v2);
        } else {
            ;
        }
        if (opposite_edge_id != MAXID) {
            auto opposite_triangle_id = uvTriMesh.E.at(opposite_edge_id).N_Fids.front();
            Parameters param__;
            auto& opposite_triangle = uvTriMesh.F.at(opposite_triangle_id);
            //if (!get_parameter(uvTriMesh, p, opposite_triangle, param__))
            {
                const auto& v_0 = V[opposite_triangle.Vids[0]];
                const auto& v_1 = V[opposite_triangle.Vids[1]];
                const auto& v_2 = V[opposite_triangle.Vids[2]];
//                param__.uv = 0.33333 * (v_0.uv + v_1.uv + v_2.uv);
                param__.uv = param.uv.x *v_0.uv + param.uv.y * v_1.uv + (1.0 - param.uv.x - param.uv.y) * v_2.uv;
            }
            // std::cout << " param__.triangleid = " << param__.triangleid << "\n";
            param__.triangleid = opposite_triangle_id;
            endingTriangleid_neighborParams[param.triangleid] = param__;
            g_param = Parameters(param__.uv, param__.triangleid);
            // std::cout << "opposite triangle = " << uvTriMesh.E.at(opposite_edge_id).N_Fids.front() <<"\n";
        } else {
            res = false;
            std::cerr << "error in get_opposite_edge_id\n";
        }
    }
    return res;
}

/*
[ v0.x v1.x v2.x ]   [ baryPosition.x ]
[ v0.y v1.y v2.y ] * [ baryPosition.y ]
[ v0.z v1.z v2.z ]   [ (1.0 - baryPosition.x - baryPosition.y) ]
*/
bool IsPointInTriangle(const glm::dvec3& v, const glm::dvec3& v0, const glm::dvec3& v1, const glm::dvec3& v2, glm::dvec3& lambda) {
    glm::dmat3x3 m = {
            v0.x, v1.x, v2.x,
            v0.y, v1.y, v2.y,
            v0.z, v1.z, v2.z
    };

    lambda = glm::inverse(m) * glm::dvec3(v.x, v.y, v.z);

    if (lambda.x > -1e-5 && lambda.y > -1e-5 && (lambda.x + lambda.y) < 1.00001) {
        return true;
    }
    std::cerr << "Error in IsPointInTriangle()\n";
    return false;
}

bool IsPointInTriangle_xy(const glm::dvec3& p, const glm::dvec3& p0, const glm::dvec3& p1, const glm::dvec3& p2, glm::dvec3& uv) {
    glm::dvec2 v02((p0.x - p2.x), (p0.y - p2.y));
    glm::dvec2 v12((p1.x - p2.x), (p1.y - p2.y));

    glm::dvec2 r(p.x, p.y);
    glm::dvec2 r3(p2.x, p2.y);

    glm::dvec2 rr3 = r - r3;

    glm::dmat2x2 T(v02, v12);
    glm::dmat2x2 T_inverse = glm::inverse(T);
    glm::dvec2 lambda = T_inverse * (r - r3);

    if (lambda.x > -1e-5 && lambda.y > -1e-5 && (lambda.x + lambda.y) < 1.00001) {
        uv.x = lambda.x;
        uv.y = lambda.y;
        uv.z = 1.0 - lambda.x - lambda.y;
        return true;
    }

    return false;
}

bool IsPointInTriangle_yz(const glm::dvec3& p, const glm::dvec3& p0, const glm::dvec3& p1, const glm::dvec3& p2, glm::dvec3& uv) {
    glm::dvec2 v02((p0.y - p2.y), (p0.z - p2.z));
    glm::dvec2 v12((p1.y - p2.y), (p1.z - p2.z));

    glm::dvec2 r(p.y, p.z);
    glm::dvec2 r3(p2.y, p2.z);

    glm::dvec2 rr3 = r - r3;

    glm::dmat2x2 T(v02, v12);
    glm::dmat2x2 T_inverse = glm::inverse(T);
    glm::dvec2 lambda = T_inverse * (r - r3);

    if (lambda.x > -1e-5 && lambda.y > -1e-5 && (lambda.x + lambda.y) < 1.00001) {
        uv.x = lambda.x;
        uv.y = lambda.y;
        uv.z = 1.0 - lambda.x - lambda.y;
        return true;
    }

    return false;
}

bool IsPointInTriangle_xz(const glm::dvec3& p, const glm::dvec3& p0, const glm::dvec3& p1, const glm::dvec3& p2, glm::dvec3& uv) {
    glm::dvec2 v02((p0.x - p2.x), (p0.z - p2.z));
    glm::dvec2 v12((p1.x - p2.x), (p1.z - p2.z));

    glm::dvec2 r(p.x, p.z);
    glm::dvec2 r3(p2.x, p2.z);

    glm::dvec2 rr3 = r - r3;

    glm::dmat2x2 T(v02, v12);
    glm::dmat2x2 T_inverse = glm::inverse(T);
    glm::dvec2 lambda = T_inverse * (r - r3);

    if (lambda.x > -1e-5 && lambda.y > -1e-5 && (lambda.x + lambda.y) < 1.00001) {
        uv.x = lambda.x;
        uv.y = lambda.y;
        uv.z = 1.0 - lambda.x - lambda.y;
        return true;
    }

    return false;
}

bool IsPointInTriangle(const glm::dvec3& p, const glm::dvec3& p0, const glm::dvec3& p1, const glm::dvec3& p2, glm::dvec3& uv,
    std::string s = "xy") {
    if (s == "xy") return IsPointInTriangle_xy(p, p0, p1, p2, uv);
    else if (s == "yz") return IsPointInTriangle_yz(p, p0, p1, p2, uv);
    else if (s == "xz") return IsPointInTriangle_xz(p, p0, p1, p2, uv);
    else {
        std::cerr << "Error in IsPointInTriangle\n";
        return false;
    }
}

//bool IsPointInTriangle(const glm::dvec3& p, const glm::dvec3& p0, const glm::dvec3& p1, const glm::dvec3& p2, glm::dvec3& uv,
//    std::string s = "xy") {
//    if (s == "xy") return IsPointInTriangle_xy(glm::dvec3(p), p0, p1, p2, uv);
//    else if (s == "yz") return IsPointInTriangle_yz(glm::dvec3(p), p0, p1, p2, uv);
//    else if (s == "xz") return IsPointInTriangle_xz(glm::dvec3(p), p0, p1, p2, uv);
//    else {
//        std::cerr << "Error in IsPointInTriangle\n";
//        return false;
//    }
//}

#include "Eigen/Core"
#include "Eigen/Eigen"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/Cholesky"
#include "Eigen/LU"
#include "Eigen/SVD"
#include "Eigen/SparseCore"
bool IsPointInTriangle_Robust(const glm::dvec3& p, const glm::dvec3& p0, const glm::dvec3& p1, const glm::dvec3& p2, glm::dvec3& uvw) {
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(3, 3);
    A(0,0) = p0.x; A(0,1) = p1.x; A(0,2) = p2.x;
    A(1,0) = p0.y; A(1,1) = p1.y; A(1,2) = p2.y;
    A(2,0) = p0.z; A(2,1) = p1.z; A(2,2) = p2.z;

    Eigen::VectorXd b = Eigen::VectorXd::Zero(3);
    b(0) = p.x; b(1) = p.y; b(2) = p.z;
    Eigen::VectorXd lambda = A.llt().solve(b);
    if (lambda(0) > 0 && lambda(1) > 0 && lambda(2)  > 0
            && lambda(0) + lambda(1) + lambda(2) <= 1.000001
            && lambda(0) + lambda(1) + lambda(2) >= 0.999999)
    {
        uvw.x = lambda(0);
        uvw.y = lambda(1);
        uvw.z = lambda(2);
        return true;
    }
    std::cerr << "Error in IsPointInTriangle_Robust()\n";
    return false;
}
glm::dvec3 get_next_point_on_orig_mesh(const Vertex& v1, const Vertex& v2) {
    const auto dir = v2 - v1;
    return v2 + dir;
}

glm::dvec3 get_next_point_on_orig_mesh(const glm::dvec3& v1, const glm::dvec3& v2) {
    const auto dir = v2 - v1;
    return v2 + dir;
}

std::vector<MapPoint> get_geodesic_line(const Mesh& origTriMesh, const Mesh& uvTriMesh, const std::vector<VertexUV>& V,
        const glm::dvec2& p, const glm::dvec2& dir, const Parameters& param) {
    std::vector<MapPoint> res;
    std::set<size_t> cut_triangles;
    auto triangle_ids = get_triangle_neighbor_triangleids(uvTriMesh, param.triangleid);
    cut_triangles.insert(triangle_ids.begin(), triangle_ids.end());
    auto origin_coordinate = get_origin_coordinate(uvTriMesh, p, param);
    MapPoint v(param.uv, origin_coordinate, param.triangleid);
    res.push_back(v);
    Parameters next_param;
    auto next_p = p + dir;
    next_param.triangleid = param.triangleid;
    bool isDebugPoint = false;
    while (get_geodesic_vertex(origTriMesh, uvTriMesh, V, next_p, triangle_ids, next_param, v)) {
        triangle_ids = get_triangle_neighbor_triangleids(uvTriMesh, next_param.triangleid);
        res.push_back(v);
        next_p += dir;

        if (glm::length(v.xyz - glm::dvec3(0.79193, 0.803864, 0)) < 1e-4) {
            isDebugPoint = true;
        }
    }
    if (glm::length(v.xyz - glm::dvec3(0.79193, 0.803864, 0)) < 1e-4) {
        isDebugPoint = true;
    }
    if (next_param.triangleid == 874) {
        isDebugPoint = true;
    }
    while (is_on_cut_triangle(origTriMesh, uvTriMesh, V, next_param, next_p)) {
        cut_triangles.insert(next_param.triangleid);
        auto dir1 = res.back().xyz - res.at(res.size() - 2).xyz;
        const auto next_point_on_orig_mesh = get_next_point_on_orig_mesh(res.at(res.size() - 2).xyz, res.back().xyz);
        auto opposite_tri = uvTriMesh.F.at(g_param.triangleid);
        const auto& v0 = uvTriMesh.V[opposite_tri.Vids[0]];
        const auto& v1 = uvTriMesh.V[opposite_tri.Vids[1]];
        const auto& v2 = uvTriMesh.V[opposite_tri.Vids[2]];
        auto cut_triangles_bak = cut_triangles;
        cut_triangles.insert(opposite_tri.id);
        if (cut_triangles_bak == cut_triangles) break;
        glm::dvec3 lambda(0.33333, 0.333333, 0.333333);
        if (!IsPointInTriangle(next_point_on_orig_mesh, v0.xyz(), v1.xyz(), v2.xyz(), lambda, space))
            lambda = glm::dvec3(0.33333, 0.333333, 0.333333);
        auto uv_tri = uvTriMesh.F.at(next_param.triangleid);
        const auto& v_0 = V[opposite_tri.Vids[0]];
        const auto& v_1 = V[opposite_tri.Vids[1]];
        const auto& v_2 = V[opposite_tri.Vids[2]];
        next_p.x = lambda.x * v_0.uv.x + lambda.y * v_1.uv.x + lambda.z * v_2.uv.x;
        next_p.y = lambda.x * v_0.uv.y + lambda.y * v_1.uv.y + lambda.z * v_2.uv.y;
        v.xyz = next_point_on_orig_mesh;
        v.uv = next_p;
        v.tri_id = opposite_tri.id;
        res.push_back(v);
        //next_p += dir;
        triangle_ids = get_triangle_neighbor_triangleids(uvTriMesh, opposite_tri.id);
        next_param.triangleid = opposite_tri.id;
        const auto d = fabs(dir.x + dir.y);
        const std::vector<glm::dvec2> dirs = {{d, 0}, {-d, 0}, {0, d}, {0, -d}};
        auto next_param_bak = next_param;
        auto triangle_ids_bak = triangle_ids;
        auto next_p_bak = next_p;
        bool correct_dir = false;
        auto best_dir = dir;
        auto best_angle = -4.0;
        for (auto& dd : dirs) {
            if (correct_dir) break;
            next_param = next_param_bak;
            triangle_ids = triangle_ids_bak;
            next_p = next_p_bak;
            int count = 0;
            while (get_geodesic_vertex(origTriMesh, uvTriMesh, V, next_p, triangle_ids, next_param, v)) {
                triangle_ids = get_triangle_neighbor_triangleids(uvTriMesh, next_param.triangleid);
                res.push_back(v);
                next_p += dd;
                ++count;
                if (count == 3) {
                    auto dir2 = res.back().xyz - res.at(res.size() - 2).xyz;
                    auto angle = glm::dot(glm::normalize(dir1), glm::normalize(dir2));
                    if (angle < 1.0) {
                        if (angle > best_angle) {
                            best_angle = angle;
                            best_dir = dd;
                        }
                        res.pop_back();
                        res.pop_back();
                        res.pop_back();
                        break;
                    } else {correct_dir = true;};
                }
            }
        }
//        if (!correct_dir) std::cerr << "Error correct_dir!\n";

        next_param = next_param_bak;
        triangle_ids = triangle_ids_bak;
        next_p = next_p_bak;
        int count = 0;
        auto prev_parm = next_param;
        while (get_geodesic_vertex(origTriMesh, uvTriMesh, V, next_p, triangle_ids, next_param, v)) {
            //std::vector<size_t> pre_vids = uvTriMesh.F.at(prev_parm.triangleid).Vids;
            //std::vector<size_t> cur_vids = uvTriMesh.F.at(next_param.triangleid).Vids;
            //std::sort(pre_vids.begin(), pre_vids.end());
            //std::sort(cur_vids.begin(), cur_vids.end());

            //std::vector<size_t> vids;
            //if (glm::length(v.xyz - glm::dvec3(0.798628, 0.800812, 0)) < 1e-4) {
            //    isDebugPoint = true;
            //}
            //std::set_intersection(pre_vids.begin(), pre_vids.end(), cur_vids.begin(), cur_vids.end(), std::back_inserter(vids));
            //if (vids.size() <= 1) {
            //    next_param = prev_parm;
            //    break;
            //}
            //prev_parm = next_param;
            triangle_ids = get_triangle_neighbor_triangleids(uvTriMesh, next_param.triangleid);
            res.push_back(v);
            next_p += best_dir;
        }
        if (glm::length(v.xyz - glm::dvec3(0.798628, 0.800812, 0)) < 1e-4) {
            isDebugPoint = true;
        }
    }
    return res;
}


glm::dvec2 next_dir(const glm::dvec2& dir) {
    glm::dvec2 next = dir;
    if (dir.x > 0) {
        next.x = 0;
        next.y = -dir.x;
    } else if (dir.y < 0) {
        next.x = dir.y;
        next.y = 0;
    } else if (dir.x < 0) {
        next.x = 0;
        next.y = -dir.x;
    } else if (dir.y > 0) {
        next.x = dir.y;
        next.y = 0;
    }
    return next;
}

glm::dvec2 next_dir2(const glm::dvec2& dir) {
    glm::dvec2 next = dir;
    if (dir.x > 0) {
        next.x = 0;
        next.y = dir.x;
    } else if (dir.y < 0) {
        next.x = -dir.y;
        next.y = 0;
    } else if (dir.x < 0) {
        next.x = 0;
        next.y = dir.x;
    } else if (dir.y > 0) {
        next.x = -dir.y;
        next.y = 0;
    }
    return next;
}

glm::dvec2 next_dir3(const glm::dvec2& dir) {
    glm::dvec2 next = dir;
    if (dir.x > 0) {
        next.x = 0;
        next.y = -dir.x;
    } else if (dir.y < 0) {
        next.x = dir.y;
        next.y = 0;
    } else if (dir.x < 0) {
        next.x = 0;
        next.y = -dir.x;
    } else if (dir.y > 0) {
        next.x = dir.y;
        next.y = 0;
    }
    return next;
}

std::vector<MapPoint> get_geodesic_line2(const Mesh& origTriMesh, const Mesh& uvTriMesh, const std::vector<VertexUV>& V,
    const glm::dvec2& p, const glm::dvec2& dir, const Parameters& param) {
    std::vector<MapPoint> res;
    auto triangle_ids = get_triangle_neighbor_triangleids(uvTriMesh, param.triangleid);
    auto origin_coordinate = get_origin_coordinate(uvTriMesh, p, param);
    MapPoint v(param.uv, origin_coordinate, param.triangleid);
    res.push_back(v);
    Parameters next_param;
    auto next_p = p + dir;
    next_param.triangleid = param.triangleid;
    while (get_geodesic_vertex(origTriMesh, uvTriMesh, V, next_p, triangle_ids, next_param, v)) {
        triangle_ids = get_triangle_neighbor_triangleids(uvTriMesh, next_param.triangleid);
        res.push_back(v);
        next_p += dir;
    }
    while (is_on_cut_triangle(origTriMesh, uvTriMesh, V, next_param, next_p)) {
        auto dir1 = res.back().xyz - res.at(res.size() - 2).xyz;
        const auto next_point_on_orig_mesh = get_next_point_on_orig_mesh(res.at(res.size() - 2).xyz, res.back().xyz);
        auto opposite_tri = uvTriMesh.F.at(g_param.triangleid);
        const auto& v0 = uvTriMesh.V[opposite_tri.Vids[0]];
        const auto& v1 = uvTriMesh.V[opposite_tri.Vids[1]];
        const auto& v2 = uvTriMesh.V[opposite_tri.Vids[2]];
        glm::dvec3 lambda(0.33333, 0.333333, 0.333333);
        if (!IsPointInTriangle(next_point_on_orig_mesh, v0.xyz(), v1.xyz(), v2.xyz(), lambda, space))
            lambda = glm::dvec3(0.33333, 0.333333, 0.333333);
        auto uv_tri = uvTriMesh.F.at(next_param.triangleid);
        const auto& v_0 = V[opposite_tri.Vids[0]];
        const auto& v_1 = V[opposite_tri.Vids[1]];
        const auto& v_2 = V[opposite_tri.Vids[2]];
        next_p.x = lambda.x * v_0.uv.x + lambda.y * v_1.uv.x + lambda.z * v_2.uv.x;
        next_p.y = lambda.x * v_0.uv.y + lambda.y * v_1.uv.y + lambda.z * v_2.uv.y;
        v.xyz = next_point_on_orig_mesh;
        v.uv = next_p;
        v.tri_id = opposite_tri.id;
        res.push_back(v);
        //next_p += dir;
        triangle_ids = get_triangle_neighbor_triangleids(uvTriMesh, opposite_tri.id);
        next_param.triangleid = opposite_tri.id;
        

        bool is_on_singular_cut = false;
        for (auto eid : opposite_tri.Eids)
            if (uvTriMeshSingularCutEdgeIds.find(eid) != uvTriMeshSingularCutEdgeIds.end()) is_on_singular_cut = true;
        glm::dvec2 next_dir = dir;
        if (is_on_singular_cut) {
            glm::dvec2 next = dir;
            if (dir.x > 0) {
                next.x = 0;
                next.y = -dir.x;
            } else if (dir.y < 0) {
                next.x = dir.y;
                next.y = 0;
            } else if (dir.x < 0) {
                next.x = 0;
                next.y = -dir.x;
            } else if (dir.y > 0) {
                next.x = dir.y;
                next.y = 0;
            }
            next_dir = next;
        } else next_dir = dir;
        while (get_geodesic_vertex(origTriMesh, uvTriMesh, V, next_p, triangle_ids, next_param, v)) {
            triangle_ids = get_triangle_neighbor_triangleids(uvTriMesh, next_param.triangleid);
            res.push_back(v);
            next_p += next_dir;
        }
    }
    return res;
}


std::vector<std::vector<MapPoint>> get_geodesic_lines(const Mesh& origTriMesh, const Mesh& uvTriMesh, const std::vector<VertexUV>& V,
        const glm::dvec2& current_p, const glm::dvec2& dir, const std::set<size_t>& constraint_triangle_ids) {
    const Vertex v_p = glm::dvec3(current_p.x, current_p.y, 0.0);
    const glm::dvec2 next_p = current_p + dir;
    auto params = get_parameters(uvTriMesh, V, next_p, constraint_triangle_ids);
    std::vector<std::vector<MapPoint>> lines;
    for (auto& param : params) {
        std::vector<MapPoint> line = get_geodesic_line(origTriMesh, uvTriMesh, V, next_p, dir, param);
        //line.insert(line.begin(), v_p);
        lines.push_back(line);
    }
    return lines;
}

std::vector<std::vector<std::vector<MapPoint>>> get_geodesic_lines(const Mesh& origTriMesh, const Mesh& uvTriMesh, const std::vector<VertexUV>& V) {
    std::vector<std::vector<std::vector<MapPoint>>> singularity_line_vs; //singularity_line_vs[i][j][k] // i: line id; j : line id; k: vid;
    auto uv_boundary_vids = get_boundary_vids(uvTriMesh);
    auto singular_vids = get_singular_vids(uvTriMesh, uv_boundary_vids);
    auto d = get_step_size(V);

    for (auto singular_vid : singular_vids) {
        auto& singular_v = V[singular_vid];
        std::set<size_t> triangle_ids(uvTriMesh.V[singular_vid].N_Fids.begin(), uvTriMesh.V[singular_vid].N_Fids.end());
        const std::vector<glm::dvec2> dirs = {{d, 0}, {-d, 0}, {0, d}, {0, -d}};
        for (const auto& dir : dirs) {
            auto line_vs = get_geodesic_lines(origTriMesh, uvTriMesh, V, singular_v.uv, dir, triangle_ids);
            singularity_line_vs.push_back(line_vs);
        }
    }

    std::cout << "----------------------------------\n\n";
    return singularity_line_vs;
}

std::vector<std::vector<std::vector<MapPoint>>> get_grid_lines(const Mesh& origTriMesh, const Mesh& uvTriMesh, const std::vector<VertexUV>& V) {
    std::vector<std::vector<std::vector<MapPoint>>> singularity_line_vs; //singularity_line_vs[i][j][k] // i: line id; j : line id; k: vid;
    auto uv_boundary_vids = get_boundary_vids(uvTriMesh);
    auto singular_vids = get_singular_vids(uvTriMesh, uv_boundary_vids);
    auto d = get_step_size(V);

    for (auto singular_vid : grid_vids) {
        auto& singular_v = V[singular_vid];
        std::set<size_t> triangle_ids(uvTriMesh.V[singular_vid].N_Fids.begin(), uvTriMesh.V[singular_vid].N_Fids.end());
        const std::vector<glm::dvec2> dirs = { { d, 0 },{ -d, 0 },{ 0, d },{ 0, -d } };
        for (const auto& dir : dirs) {
            auto line_vs = get_geodesic_lines(origTriMesh, uvTriMesh, V, singular_v.uv, dir, triangle_ids);
            singularity_line_vs.push_back(line_vs);
        }
    }

    std::cout << "----------------------------------\n\n";
    return singularity_line_vs;
}

void write_lines(const char* filename, const std::vector<std::vector<std::vector<Vertex>>>& lines) {
    size_t vnum = 0, lnum = 0, lvnum = 0;
    for (auto& singularity_lines : lines)
        for (auto& singularity_line : singularity_lines) {
            vnum += singularity_line.size();
            lvnum += singularity_line.size() + 1;
            ++lnum;
        }
    std::ofstream ofs(filename);
    ofs << "# vtk DataFile Version 3.0\n"
        << filename << "\n"
        << "ASCII\n\n"
        << "DATASET POLYDATA\n";
    ofs << "POINTS " << vnum << " double\n";

    for (auto& singularity_lines : lines)
        for (auto& singularity_line : singularity_lines)
            for (auto& v : singularity_line)
                ofs << v.x << " " << v.y << " " << v.z << "\n";
    size_t vid = 0;
    ofs << "LINES " << lnum << " " << lvnum << "\n";
    for (auto& singularity_lines : lines)
        for (auto& singularity_line : singularity_lines) {
            ofs << singularity_line.size() << " ";
            for (auto& v : singularity_line) ofs << vid++ << " ";
        }
}

void write_lines(const char* filename, const std::vector<std::vector<std::vector<MapPoint>>>& lines) {
    size_t vnum = 0, lnum = 0, lvnum = 0;
    for (auto& singularity_lines : lines)
        for (auto& singularity_line : singularity_lines) {
            vnum += singularity_line.size();
            lvnum += singularity_line.size() + 1;
            ++lnum;
        }
    std::ofstream ofs(filename);
    ofs << "# vtk DataFile Version 3.0\n"
        << filename << "\n"
        << "ASCII\n\n"
        << "DATASET POLYDATA\n";
    ofs << "POINTS " << vnum << " double\n";

    for (auto& singularity_lines : lines)
        for (auto& singularity_line : singularity_lines)
            for (auto& v : singularity_line)
                ofs << v.xyz.x << " " << v.xyz.y << " " << v.xyz.z << "\n";
    size_t vid = 0;
    ofs << "LINES " << lnum << " " << lvnum << "\n";
    for (auto& singularity_lines : lines)
        for (auto& singularity_line : singularity_lines) {
            ofs << singularity_line.size() << " ";
            for (auto& v : singularity_line) ofs << vid++ << " ";
        }
}

void write_lines(const char* filename, std::vector<std::vector<MapPoint>>& lines) {
    size_t vnum = 0, lnum = 0, lvnum = 0;
    for (auto& singularity_line : lines) {
        vnum += singularity_line.size();
        lvnum += singularity_line.size() + 1;
        ++lnum;
    }
    std::ofstream ofs(filename);
    ofs << "# vtk DataFile Version 3.0\n"
        << filename << "\n"
        << "ASCII\n\n"
        << "DATASET POLYDATA\n";
    ofs << "POINTS " << vnum << " double\n";

    for (auto& singularity_line : lines)
        for (auto& v : singularity_line)
            ofs << v.xyz.x << " " << v.xyz.y << " " << v.xyz.z << "\n";
    size_t vid = 0;
    ofs << "LINES " << lnum << " " << lvnum << "\n";
    for (auto& singularity_line : lines) {
        ofs << singularity_line.size() << " ";
        for (auto& v : singularity_line) ofs << vid++ << " ";
    }
}

void write_lines_uv(const char* filename, const std::vector<std::vector<std::vector<MapPoint>>>& lines) {
    size_t vnum = 0, lnum = 0, lvnum = 0;
    for (auto& singularity_lines : lines)
        for (auto& singularity_line : singularity_lines) {
            vnum += singularity_line.size();
            lvnum += singularity_line.size() + 1;
            ++lnum;
        }
    std::ofstream ofs(filename);
    ofs << "# vtk DataFile Version 3.0\n"
        << filename << "\n"
        << "ASCII\n\n"
        << "DATASET POLYDATA\n";
    ofs << "POINTS " << vnum << " double\n";

    for (auto& singularity_lines : lines)
        for (auto& singularity_line : singularity_lines)
            for (auto& v : singularity_line)
                ofs << v.uv.x << " " << v.uv.y << " " << 0 << "\n";
    size_t vid = 0;
    ofs << "LINES " << lnum << " " << lvnum << "\n";
    for (auto& singularity_lines : lines)
        for (auto& singularity_line : singularity_lines) {
            ofs << singularity_line.size() << " ";
            for (auto& v : singularity_line) ofs << vid++ << " ";
        }
}

std::set<size_t> get_corner_vids(const Mesh& origTriMesh, const Mesh& uvTriMesh, const std::vector<VertexUV>& V) {
    std::set<size_t> corner_vids;
    for (auto& v : uvTriMesh.V) {
        if (v.isBoundary && cutVertexIds.find(v.id) == cutVertexIds.end()) {
            bool is_on_cut = false;
            for (auto& nv : v.N_Vids)
                if (cutVertexIds.find(nv) != cutVertexIds.end()) {
                    is_on_cut = true;
                    break;
                }
            if (is_on_cut) continue;
            std::vector<size_t> N_Vids;
            for (auto& nv : v.N_Vids)
                if (uvTriMesh.V.at(nv).isBoundary)
                    N_Vids.push_back(nv);
            auto v1 = uvTriMesh.V.at(N_Vids[0]).xyz();
            auto v2 = uvTriMesh.V.at(N_Vids[1]).xyz();
            auto d1 = glm::normalize(v1 - v);
            auto d2 = -glm::normalize(v2 - v);
            auto cosangle = glm::dot(d1, d2);
            if (cosangle < 0.5) corner_vids.insert(v.id);
        }
    }
    std::vector<size_t> res(corner_vids.begin(), corner_vids.end());
    {
        MeshFileWriter writer(uvTriMesh, "corners.vtk");
        writer.WriteVerticesVtk(res);
    }
    return corner_vids;
}

std::vector<MapPoint> get_boundary_line(const Mesh& origTriMesh, const Mesh& uvTriMesh, const std::vector<VertexUV>& V, 
    std::set<size_t>& ancor_vids, size_t start_vid, size_t next_vid, size_t& end_vid) {
    std::vector<MapPoint> res;
    MapPoint mp;
    mp.xyz = origTriMesh.V[start_vid].xyz();
    res.push_back(mp);
    while (ancor_vids.find(next_vid) == ancor_vids.end()) {
        mp.xyz = origTriMesh.V[next_vid].xyz();
        res.push_back(mp);
        auto& v = origTriMesh.V[next_vid];
        for (auto nv : v.N_Vids)
            if (nv != start_vid && origTriMesh.V.at(nv).isBoundary) {
                start_vid = next_vid;
                next_vid = nv;
                break;
            }
    }
    mp.xyz = origTriMesh.V[next_vid].xyz();
    res.push_back(mp);
    end_vid = next_vid;
    return res;
}

std::vector<std::vector<MapPoint>> get_boundary_lines(const Mesh& origTriMesh, const Mesh& uvTriMesh,
    const std::vector<VertexUV>& V, std::set<size_t>& ancor_vids, std::map<size_t, std::vector<MapPoint>>& nodepair_line_) {
    std::vector<std::vector<MapPoint>> res;
    std::map<size_t, std::vector<MapPoint>> nodepair_line;
    std::set<size_t> orig_ancor_vids;
    std::map<size_t, size_t> orig_uv_ancor_vids;
    for (auto vid : ancor_vids) {
        orig_ancor_vids.insert(V.at(vid).father);
        orig_uv_ancor_vids[V.at(vid).father] = vid;
    }
    for (auto vid : orig_ancor_vids) {
        auto& v = origTriMesh.V.at(vid);
        std::vector<size_t> N_Vids;
        for (auto& nv : v.N_Vids)
            if (origTriMesh.V.at(nv).isBoundary)
                N_Vids.push_back(nv);
        auto v1 = origTriMesh.V.at(N_Vids[0]);
        auto v2 = origTriMesh.V.at(N_Vids[1]);
        {
            size_t end_vid = MAXID;
            auto line1 = get_boundary_line(origTriMesh, uvTriMesh, V, orig_ancor_vids, v.id, v1.id, end_vid);
            size_t p1 = (vid << 32) + end_vid;
            size_t p2 = (end_vid << 32) + vid;
            if (nodepair_line.find(p1) == nodepair_line.end() && nodepair_line.find(p2) == nodepair_line.end()) {
                nodepair_line[p1] = line1;
            }
        }
        {
            size_t end_vid = MAXID;
            auto line = get_boundary_line(origTriMesh, uvTriMesh, V, orig_ancor_vids, v.id, v2.id, end_vid);
            size_t p1 = (vid << 32) + end_vid;
            size_t p2 = (end_vid << 32) + vid;
            if (nodepair_line.find(p1) == nodepair_line.end() && nodepair_line.find(p2) == nodepair_line.end()) {
                nodepair_line[p1] = line;
            }
        }
    }
    for (auto& item : nodepair_line)
        res.push_back(item.second);

    for (auto& item : nodepair_line) {
        size_t p1 = orig_uv_ancor_vids[item.first >> 32];
        size_t p2 = orig_uv_ancor_vids[item.first & 0xffffffff];
        nodepair_line_[(p1 << 32) + p2] = item.second;
    }
    return res;
}

std::vector<glm::vec2> get_parameters(std::vector<MapPoint>& l1, std::vector<MapPoint>& l2, std::vector<MapPoint>& l3, std::vector<MapPoint>& l4) {
    std::vector<glm::vec2> res;
    return res;
}
void generate_quad_mesh(const Mesh& origTriMesh, const Mesh& uvTriMesh, const std::vector<VertexUV>& V,
    std::vector<std::vector<MapPoint>>& lines, const std::vector<std::vector<size_t>>& combs,
    std::vector<Vertex>& quadV, std::vector<Face>& quadF) {
    //std::vector<glm::dvec3> newV;
    //std::set<glm::dvec3> coarseV;
    float d = 0;
    for (auto& e : uvTriMesh.E)
        d += glm::length(uvTriMesh.V.at(e.Vids[0]).xyz() - uvTriMesh.V.at(e.Vids[1]).xyz());
    d /= uvTriMesh.E.size();

    Vertex temp_v;
    for (auto& line : lines) {
        auto f = glm::dvec3(line.front().xyz);
        auto b = glm::dvec3(line.back().xyz);
        int id = -1;
        for (auto& v : quadV) {
            if (glm::length(v - f) < 1e-5) {
                id = v.id;
                break;
            }
        }
        if (id == -1) {
            temp_v = f;
            temp_v.id = quadV.size();
            quadV.push_back(temp_v);
        }
        id = -1;
        for (auto& v : quadV) {
            if (glm::length(v - b) < 1e-5) {
                id = v.id;
                break;
            }
        }
        if (id == -1) {
            temp_v = b;
            temp_v.id = quadV.size();
            quadV.push_back(temp_v);
        }
    }
    std::cout << "quadV.size() = " << quadV.size() << "\n";

    std::vector<std::vector<MapPoint>> new_lines;
    for (auto& line : lines) {
        std::vector<MapPoint> new_line;
        new_line.push_back(line.front());
        auto gap = round(double (line.size()) / 10);
        for (int i = gap; i < line.size() && new_line.size() < 11; i += gap) {
            new_line.push_back(line.at(i));
            //double min_distance = 10000000;
            //double closest_vid;
            //for (auto&v : uvTriMesh.V) {
            //    auto distance = glm::length(glm::dvec3(v.xyz()) - line.at(i).xyz);
            //    if (distance < min_distance) {
            //        min_distance = distance;
            //        closest_vid = v.id;
            //    }
            //}
            //grid_vids.push_back(closest_vid);
        }
        new_line.push_back(line.back());
        new_lines.push_back(new_line);
    }
    lines = new_lines;
}

std::vector<std::vector<MapPoint>> get_nonoverlap_geodesic_lines(Mesh& origTriMesh, Mesh& uvTriMesh, std::vector<VertexUV>& V, 
    std::vector<std::vector<std::vector<MapPoint>>>& lines) {
    std::set<size_t> corner_vids = get_corner_vids(origTriMesh, uvTriMesh, V);
    auto ancor_vids = corner_vids;
    auto eps = 10 * glm::length(lines.front().front().front().xyz - (lines.front().front().begin() + 1)->xyz);
    std::map<size_t, std::vector<MapPoint>> nodepair_line;
    std::map<size_t, std::pair<size_t, std::vector<MapPoint>>> nodepair_midTriId_line;
    auto uv_boundary_vids = get_boundary_vids(uvTriMesh);
    auto singular_vids = get_singular_vids(uvTriMesh, uv_boundary_vids);
    auto d = get_step_size(V);
    for (auto& singularity_lines : lines) {
        for (auto& singularity_line : singularity_lines) {
            auto begin_v = singularity_line.front().xyz;
            glm::dvec3 begin_singular_v;
            size_t begin_singular_vid;
            size_t begin_vid = 0;
            for (auto singular_vid : singular_vids) {
                auto singular_v = glm::dvec3(uvTriMesh.V.at(singular_vid).xyz());
                auto dis = glm::length(begin_v - singular_v);
                if (dis < eps) {
                    begin_singular_v = singular_v;
                    begin_singular_vid = singular_vid;
                    singularity_line.front().xyz = singular_v;
                    break;
                }
            }
            for (size_t vid = 1; vid < singularity_line.size(); ++vid) {
                auto& v = singularity_line[vid];
                for (auto singular_vid : singular_vids) {
                    auto singular_v = glm::dvec3(uvTriMesh.V.at(singular_vid).xyz());
                    auto dis = glm::length(v.xyz - singular_v);
                    if (dis < eps && singular_vid != begin_singular_vid) {
                        size_t closest_id = vid;
                        auto min_dis = dis;
                        for (auto id = vid; (id < singularity_line.size() - 1) && id < vid + 20; ++id) {
                            auto& v_ = singularity_line[id];
                            auto dis_ = glm::length(v_.xyz - singular_v);
                            if (dis_ < min_dis) {
                                min_dis = dis_;
                                closest_id = id;
                            }
                        }
                        vid = closest_id;
                        singularity_line.at(vid).xyz = singular_v;
                        size_t p1 = (begin_singular_vid << 32) + singular_vid;
                        size_t p2 = (singular_vid << 32) + begin_singular_vid;
                        //size_t mid_tri_id = (singularity_line.begin() + begin_vid + (vid - begin_vid) / 2)->tri_id;
                        if (nodepair_line.find(p1) == nodepair_line.end() && nodepair_line.find(p2) == nodepair_line.end()) {
                            nodepair_line[p1] = std::vector<MapPoint>(singularity_line.begin() + begin_vid, 
                                singularity_line.begin() + vid + 1);
                            //nodepair_midTriId_line[p1] = std::make_pair(mid_tri_id, std::vector<MapPoint>(singularity_line.begin() + begin_vid,
                            //    singularity_line.begin() + vid + 1));
                        } 
                        //else if (nodepair_midTriId_line.find(p1) != nodepair_midTriId_line.end()) {
                        //    auto it = nodepair_midTriId_line.find(p1);
                        //    auto ntris = get_triangle_neighbor_triangleids(uvTriMesh, mid_tri_id);
                        //    if (ntris.find(it->second.first) == ntris.end()) {
                        //        nodepair_line[p1] = std::vector<MapPoint>(singularity_line.begin() + begin_vid,
                        //            singularity_line.begin() + vid + 1);
                        //        nodepair_midTriId_line[p1] = std::make_pair(mid_tri_id, std::vector<MapPoint>(singularity_line.begin() + begin_vid,
                        //            singularity_line.begin() + vid + 1));
                        //    }
                        //} else if (nodepair_midTriId_line.find(p2) != nodepair_midTriId_line.end()) {
                        //    auto it = nodepair_midTriId_line.find(p2);
                        //    auto ntris = get_triangle_neighbor_triangleids(uvTriMesh, mid_tri_id);
                        //    if (ntris.find(it->second.first) == ntris.end()) {
                        //        nodepair_line[p1] = std::vector<MapPoint>(singularity_line.begin() + begin_vid,
                        //            singularity_line.begin() + vid + 1);
                        //        nodepair_midTriId_line[p1] = std::make_pair(mid_tri_id, std::vector<MapPoint>(singularity_line.begin() + begin_vid,
                        //            singularity_line.begin() + vid + 1));
                        //    }
                        //}
                        begin_singular_v = singular_v;
                        begin_singular_vid = singular_vid; 
                        begin_vid = vid;
                        begin_v = singularity_line[vid].xyz;
                        --vid;
                        break;
                    }
                }
            }
            auto end_v = singularity_line.back().xyz;
            glm::dvec3 end_singular_v;
            size_t end_singular_vid = -1;
            for (auto singular_vid : singular_vids) {
                auto singular_v = glm::dvec3(uvTriMesh.V.at(singular_vid).xyz());
                if (glm::length(end_v - singular_v) < eps && singular_vid != begin_singular_vid) {
                    end_singular_v = singular_v;
                    end_singular_vid = singular_vid;
                    break;
                }
            }
            if (end_singular_vid == -1) {
                bool is_ending_in_singularities = false;
                auto& ending_tri = uvTriMesh.F.at(singularity_line.back().tri_id);
                for (auto singular_vid : singular_vids) {
                    auto& singular_v = uvTriMesh.V.at(singular_vid);
                    std::set<size_t> tri_ids(singular_v.N_Fids.begin(), singular_v.N_Fids.end());
                    if (tri_ids.find(ending_tri.id) != tri_ids.end()) {
                        is_ending_in_singularities = true;
                        break;
                    }
                }
                //if (!uvTriMesh.E.at(ending_tri.Eids[0]).isBoundary &&
                //    !uvTriMesh.E.at(ending_tri.Eids[1]).isBoundary &&
                //    !uvTriMesh.E.at(ending_tri.Eids[2]).isBoundary) {
                //    is_ending_in_singularities = true;
                //}
                if (uvTriMeshSingularCutEdgeIds.find(ending_tri.Eids[0]) == uvTriMeshSingularCutEdgeIds.end() &&
                    uvTriMeshSingularCutEdgeIds.find(ending_tri.Eids[1]) == uvTriMeshSingularCutEdgeIds.end() &&
                    uvTriMeshSingularCutEdgeIds.find(ending_tri.Eids[2]) == uvTriMeshSingularCutEdgeIds.end() && (
                    cutEdgeIds.find(ending_tri.Eids[0]) != cutEdgeIds.end() ||
                    cutEdgeIds.find(ending_tri.Eids[1]) != cutEdgeIds.end() ||
                    cutEdgeIds.find(ending_tri.Eids[2]) != cutEdgeIds.end() )) {
                    is_ending_in_singularities = true;
                }
                if (is_ending_in_singularities) continue;
                auto closest_vid = ending_tri.Vids.front();
                double min_dis = 1000000.0;
                //Vertex new_v = glm::dvec3(singularity_line.back().xyz);
                //new_v.id = uvTriMesh.V.size();
                //new_v.isBoundary = true;
                for (auto id : ending_tri.Vids) {
                    if (!uvTriMesh.V.at(id).isBoundary) continue;
                    auto v_ = glm::dvec3(uvTriMesh.V.at(id).xyz());
                    auto dis_ = glm::length(v_ - singularity_line.back().xyz);
                    if (dis_ < min_dis) {
                        closest_vid = id;
                        min_dis = dis_;
                    }
                    //new_v.N_Vids.push_back(id);
                    //uvTriMesh.V.at(id).N_Vids.push_back(new_v.id);
                }
                //uvTriMesh.V.push_back(new_v);
                //uvTriMesh.V.at(closest_vid) = glm::dvec3(singularity_line.back().xyz);
                //singularity_line.back().xyz = glm::dvec3(uvTriMesh.V.at(closest_vid).xyz());

                //size_t p1 = (begin_singular_vid << 32) + closest_vid;
                //size_t p2 = (closest_vid << 32) + begin_singular_vid;
                //size_t pt1 = (begin_singular_vid << 32) + ending_tri.id;
                //size_t pt2 = (ending_tri.id << 32) + begin_singular_vid;
                //if (nodepair_line.find(p1) == nodepair_line.end() && nodepair_line.find(p2) == nodepair_line.end() &&
                //    nodepair_line.find(pt1) == nodepair_line.end() && nodepair_line.find(pt2) == nodepair_line.end()) {
                //    nodepair_line[p1] = std::vector<MapPoint>(singularity_line.begin() + begin_vid, singularity_line.end());
                //    nodepair_line[pt1] = std::vector<MapPoint>();
                //    ancor_vids.insert(closest_vid);
                //}

                size_t p1 = (begin_singular_vid << 32) + closest_vid;
                size_t p2 = (closest_vid << 32) + begin_singular_vid;
                size_t pt1 = (begin_singular_vid << 32) + ending_tri.id;
                size_t pt2 = (ending_tri.id << 32) + begin_singular_vid;
                if (nodepair_line.find(p1) == nodepair_line.end() && nodepair_line.find(p2) == nodepair_line.end() &&
                    nodepair_line.find(pt1) == nodepair_line.end() && nodepair_line.find(pt2) == nodepair_line.end()) {                 
                    Vertex new_orig_v = glm::dvec3(singularity_line.back().xyz);
                    new_orig_v.id = origTriMesh.V.size();
                    new_orig_v.isBoundary = true;

                    Vertex new_v = glm::dvec3(singularity_line.back().xyz);
                    new_v.id = uvTriMesh.V.size();
                    new_v.isBoundary = true;
                    for (auto id : ending_tri.Vids) {
                        if (!uvTriMesh.V.at(id).isBoundary) continue;
                        new_v.N_Vids.push_back(id);
                        new_orig_v.N_Vids.push_back(V[id].father);
                        //uvTriMesh.V.at(id).N_Vids.push_back(new_v.id);
                        //origTriMesh.V.at(V[id].father).N_Vids.push_back(new_orig_v.id);
                        uvTriMesh.V.at(id).N_Vids.insert(uvTriMesh.V.at(id).N_Vids.begin(), new_v.id);
                        origTriMesh.V.at(V[id].father).N_Vids.insert(origTriMesh.V.at(V[id].father).N_Vids.begin(), new_orig_v.id);
                    }
                    origTriMesh.V.push_back(new_orig_v);
                    uvTriMesh.V.push_back(new_v);

                    VertexUV newv;
                    newv.x = singularity_line.back().xyz.x;
                    newv.y = singularity_line.back().xyz.y;
                    newv.z = singularity_line.back().xyz.z;
                    newv.id = V.size();
                    newv.uv = singularity_line.back().uv;
                    newv.triVid = singularity_line.back().tri_id;
                    newv.father = new_orig_v.id;
                    V.push_back(newv);
                    p1 = (begin_singular_vid << 32) + new_v.id;
                    nodepair_line[p1] = std::vector<MapPoint>(singularity_line.begin() + begin_vid, singularity_line.end());
                    nodepair_line[pt1] = std::vector<MapPoint>();
                    ancor_vids.insert(new_v.id);
                }
            }
        }
    }
    {
        MeshFileWriter writer(uvTriMesh, "ancor.vids.vtk");
        writer.WriteVerticesVtk(std::vector<size_t>(ancor_vids.begin(), ancor_vids.end()));
    }
    for (auto& item : nodepair_line) {
        ;
    }
    std::vector<std::vector<MapPoint>> res;
    std::map<size_t, std::vector<MapPoint>> nodepair_line_;
    auto boundary_lines = get_boundary_lines(origTriMesh, uvTriMesh, V, ancor_vids, nodepair_line_);
    for (auto& item : nodepair_line)
        if (!item.second.empty()) {
            res.push_back(item.second);
            //nodepair_line_.insert(item);
        }
    {
        std::vector<std::vector<MapPoint>> new_lines;
        for (auto& line : res) {
            std::vector<MapPoint> new_line;
            new_line.push_back(line.front());
            auto gap = round(double(line.size()) / 4);
            for (int i = gap; i < line.size() && new_line.size() < 5; i += gap) {
                new_line.push_back(line.at(i));
                double min_distance = 10000000;
                double closest_vid;
                for (auto&v : uvTriMesh.V) {
                    auto distance = glm::length(glm::dvec3(v.xyz()) - line.at(i).xyz);
                    if (distance < min_distance) {
                        min_distance = distance;
                        closest_vid = v.id;
                    }
                }
                grid_vids.push_back(closest_vid);
            }
            new_line.push_back(line.back());
            new_lines.push_back(new_line);
        }
    }
    std::copy(boundary_lines.begin(), boundary_lines.end(), std::back_inserter(res));

    std::vector<size_t> ps;
    for (auto& item : nodepair_line)
        if (!item.second.empty()) ps.push_back(item.first);
    for (auto& item : nodepair_line_)
        if (!item.second.empty()) ps.push_back(item.first);

    std::cout << "############################\n";
    auto combs = Util::combine(ps.size(), 4);
    for (auto& comb : combs) {
        std::set<size_t> s;
        for (auto n : comb) {
            auto node = ps[n];
            auto n1 = node >> 32;
            auto n2 = node & 0xffffffff;
            s.insert(n1);
            s.insert(n2);
        }
        if (s.size() == 4) {
            size_t singularity_count = 0;
            for (auto n : s)
                for (auto svid : singular_vids)
                    if (svid == n) {
                        ++singularity_count;
                        break;
                    }
            if (singularity_count == 4) continue;

            size_t boundary_count = 0;
            for (auto n : s)
                if (origTriMesh.V.at(V.at(n).father).isBoundary) {
                    ++boundary_count;
                }
            if (boundary_count == 4) continue;
            for (auto n : comb) std::cout << n << " ";
            std::cout << "\n";
        }
    }
    std::cout << "############################\n";
    std::vector<Vertex> quadV;
    std::vector<Face> quadF;
    generate_quad_mesh(origTriMesh, uvTriMesh, V, res, combs, quadV, quadF);
    return res;
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cout << "Usage: ExtractQuadMesh origTri.vtk tri.uv.m output.quad.vtk grid=<100>\n";
        return -1;
    }
    if (argc >= 5) sqsize = std::stoi(argv[4]);
    if (argc >= 6) space = argv[5];
    MeshFileReader origTriReader(argv[1]);
    Mesh& origTriMesh = (Mesh&)origTriReader.GetMesh();
    origTriMesh.RemoveUselessVertices();
    origTriMesh.BuildAllConnectivities();
    origTriMesh.ExtractBoundary();

    std::vector<VertexUV> V;
    std::vector<Face> F;
    read(argv[2], V, F);
    Mesh uvTriMesh;
    generate_tri_mesh(V, F, uvTriMesh);

    std::vector<Vertex> triV(V.size());
    for (auto& v : V) {
        triV[v.id] = glm::dvec3(v.uv.x, v.uv.y, 0.0);
        triV[v.id].id = v.id;
    }
    Mesh uvMesh;
    generate_tri_mesh(triV, F, uvMesh);
    write("uv.tri.vtk", triV, F);
    write("tri.vtk", V, F);

    auto uv_boundary_eids = get_boundary_eids(uvTriMesh);
    auto uv_boundary_vids = get_boundary_vids(uvTriMesh);
    auto singular_vids = get_singular_vids(uvTriMesh, uv_boundary_vids);
    for (auto svid : singular_vids)
        uvTriMesh.V[svid].isSingularity = true;
    {
        std::cout << "------------------------------------\n";
        std::cout << "Writing tri.boundary.vtk...\n";
        MeshFileWriter writer(uvTriMesh, "tri.boundary.vtk");
        writer.WriteEdgesVtk(uv_boundary_eids);
        std::cout << "Finished writing tri.boundary.vtk...\n";
    }
    {
        std::cout << "------------------------------------\n";
        std::cout << "Writing uv.tri.boundary.vtk...\n";
        MeshFileWriter writer(uvMesh, "uv.tri.boundary.vtk");
        writer.WriteEdgesVtk(uv_boundary_eids);
        std::cout << "Finished writing uv.tri.boundary.vtk...\n";
    }
    {
        std::cout << "------------------------------------\n";
        std::cout << "Writing singularities.vtk...\n";
        MeshFileWriter writer(uvTriMesh, "singularities.vtk");
        writer.WriteVerticesVtk(singular_vids);
        std::cout << "Finished writing singularities.vtk...\n";
    }
    {
        std::cout << "------------------------------------\n";
        std::cout << "Writing uv.singularities.vtk...\n";
        MeshFileWriter writer(uvMesh, "uv.singularities.vtk");
        writer.WriteVerticesVtk(singular_vids);
        std::cout << "Finished writing uv.singularities.vtk...\n";
    }

    auto cut_eids = get_uvTriMeshCutEdgeIds(origTriMesh, uvTriMesh, V);
    std::cout << "cut_eids.size() = " << cut_eids.size() << "\n";
    {
        MeshFileWriter writer(uvTriMesh, "cut.eids.vtk");
        writer.WriteEdgesVtk(cut_eids);
    }
    auto singular_cut_eids = get_uvTriMeshSingularCutEdgeIds(origTriMesh, uvTriMesh, V);

    {
        MeshFileWriter writer(uvTriMesh, "cut.singular.eids.vtk");
        writer.WriteEdgesVtk(singular_cut_eids);
    }
    std::cout << "uvTriMeshSingularCutEdgeIds.size() = " << uvTriMeshSingularCutEdgeIds.size() << "\n";
    auto singularity_line_vs = get_geodesic_lines(origTriMesh, uvTriMesh, V);
    write_lines("lines.vtk", singularity_line_vs);
    write_lines_uv("lines_uv.vtk", singularity_line_vs);

    auto lines = get_nonoverlap_geodesic_lines(origTriMesh, uvTriMesh, V, singularity_line_vs);
    write_lines("non_overlap_lines.vtk", lines);
    //std::set<size_t> corner_vids = get_corner_vids(origTriMesh, uvTriMesh, V);


    //auto grid_line_vs = get_grid_lines(origTriMesh, uvTriMesh, V);
    //write_lines("grid_lines.vtk", grid_line_vs);
    //write_lines_uv("lines_uv.vtk", singularity_line_vs);
    return 0;
}
