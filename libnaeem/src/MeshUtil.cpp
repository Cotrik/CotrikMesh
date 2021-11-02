#include <algorithm>
#include <math.h>
#include "MeshUtil.h"

#define PI 3.14159265

MeshUtil::MeshUtil() {}

MeshUtil::MeshUtil(Mesh& mesh_) : mesh(mesh_) {}

MeshUtil::~MeshUtil() {}

void MeshUtil::CheckValidity() {
    if (mesh.V.size() == 0 || mesh.F.size() == 0 || mesh.C.size() == 0) {
        std::cout << "No mesh to use for MeshUtils." << std::endl;
        exit(0);
    }
}

void MeshUtil::SetMesh(Mesh& mesh_) {
    mesh = mesh_;
    SetMeshArea();
}

void MeshUtil::SetMeshArea() {
    CheckValidity();

    double area = 0.0;
    for (auto& f: mesh.F) {
        area += GetFaceArea(f.id);
    }
    mesh.totalArea = area;
}

double MeshUtil::GetMeshArea() {
    CheckValidity();

    if (mesh.totalArea == 0.0) {
        SetMeshArea();
    }
    return mesh.totalArea;
}

double MeshUtil::GetFaceArea(int fid) {
    CheckValidity();

    Face& f = mesh.F.at(fid);
    glm::dvec3 a = mesh.V.at(f.Vids.at(0)).xyz() - mesh.V.at(f.Vids.at(2)).xyz();
    glm::dvec3 b = mesh.V.at(f.Vids.at(1)).xyz() - mesh.V.at(f.Vids.at(3)).xyz();

    return 0.5 * glm::length(glm::cross(a,b));
}

double MeshUtil::GetVertexEnergy(int vid) {
    CheckValidity();

    auto& v = mesh.V.at(vid);
    int n = v.N_Vids.size();
    double ideal_angle = 180 - (360 / n);
    double e = 0.0;
    for (auto eid: v.N_Eids) {
        e += ideal_angle / GetInteriorAngleAtEdge(vid, eid);
    }

    return std::max((double) n, e);
}

double MeshUtil::GetInteriorAngleAtEdge(int vid, int eid) {
    CheckValidity();

    Vertex& v = mesh.V.at(vid);
    Edge& e = mesh.E.at(eid);
    int vid2 = e.Vids.at(0) == vid ? e.Vids.at(1) : e.Vids.at(0);
    std::vector<int> nvids;
    for (auto fid: e.N_Fids) {
        Face& f = mesh.F.at(fid);
        for (auto feid: f.Eids) {
            if (feid == eid) continue;
            auto& fe = mesh.E.at(feid);
            if (std::find(fe.Vids.begin(), fe.Vids.end(), vid) != fe.Vids.end()) {
                int nvid = fe.Vids.at(0) == vid ? fe.Vids.at(1) : fe.Vids.at(0);
                nvids.push_back(nvid);
            }
        }
    }

    Vertex& v_n = mesh.V.at(vid2);
    glm::dvec3 a = glm::normalize(v.xyz() - v_n.xyz());
    glm::dvec3 b = glm::normalize(mesh.V.at(nvids.at(0)).xyz() - v_n.xyz());
    glm::dvec3 c = glm::normalize(mesh.V.at(nvids.at(1)).xyz() - v_n.xyz());

    double alpha1 = acos(std::max(-1.0, std::min(glm::dot(a,b), 1.0))) * 180.0 / PI;
    double alpha2 = acos(std::max(-1.0, std::min(glm::dot(a,c), 1.0))) * 180.0 / PI;
    // std::cout << "alpha1: " << alpha1 << " alpha2: " << alpha2 << " total: " << alpha1 + alpha2 << std::endl;

    return alpha1 + alpha2;
}

