#include <algorithm>
#include <math.h>
#include <vtkSmartPointer.h>
#include <vtkPolygon.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>

#include "ParallelFor.h"
#include "MeshUtil.h"

#define PI 3.14159265

MeshUtil::MeshUtil() {}

MeshUtil::MeshUtil(const MeshUtil& r) {
    mesh = r.mesh;
}

MeshUtil::MeshUtil(Mesh& mesh_) : mesh(&mesh_) {}

MeshUtil::~MeshUtil() {}

void MeshUtil::CheckValidity() {
    if (mesh == NULL) {
        std::cout << "No mesh to use for MeshUtils." << std::endl;
        exit(0);
    }
    if (mesh->V.size() == 0 || mesh->F.size() == 0 || mesh->C.size() == 0) {
        std::cout << "No mesh to use for MeshUtils." << std::endl;
        exit(0);
    }
}

void MeshUtil::SetMembers(Mesh& mesh_) {
    mesh = &mesh_;
    SetMeshArea();
}

vtkSmartPointer<vtkPolyData> MeshUtil::GetPolyData() {
    CheckValidity();

    vtkNew<vtkPoints> points;
    vtkNew<vtkCellArray> cells;
    auto polyData = vtkSmartPointer<vtkPolyData>::New();

    for (auto& v: mesh->V) {
        vtkIdType pId = v.id;
        points->InsertPoint(pId, v.x, v.y, v.z);
    }
    for (auto& c: mesh->C) {
        vtkNew<vtkPolygon> p;
        p->GetPointIds()->SetNumberOfIds(c.Vids.size());
        for (int i = 0; i < c.Vids.size(); i++) {
            p->GetPointIds()->SetId(i, c.Vids.at(i));
        }
        cells->InsertNextCell(p);
    }
    
    polyData->SetPoints(points);
    polyData->SetPolys(cells);

    return polyData;
}

vtkSmartPointer<vtkPolyData> MeshUtil::GetPolyData(Mesh& mesh_, size_t vid) {
    CheckValidity();

    // vtkNew<vtkPoints> points;
    // vtkNew<vtkCellArray> cells;
    auto polyData = vtkSmartPointer<vtkPolyData>::New();
    return polyData;
    // std::vector<size_t> vertices;
    // std::vector<size_t> faces = mesh_.V.at(vid).N_Fids;
    // for (auto fid: faces) {
    //     auto& f = mesh_.F.at(fid);
    //     if (f.N_Fids.empty() || f.Vids.empty()) continue;
    //     AddContents(vertices, f.Vids);
    //     vtkNew<vtkPolygon> p;
    //     p->GetPointIds()->SetNumberOfIds(f.Vids.size());
    //     for (int i = 0; i < f.Vids.size(); i++) {
    //         p->GetPointIds()->SetId(i, f.Vids.at(i));
    //     }
    //     cells->InsertNextCell(p);
    // }
    // for (auto vid: vertices) {
    //     auto& v = mesh_.V.at(vid);
    //     vtkIdType pId = v.id;
    //     points->InsertPoint(pId, v.x, v.y, v.z);
    // }

    // polyData->SetPoints(points);
    // polyData->SetPolys(cells);
    // return polyData;
}

vtkSmartPointer<vtkPolyData> MeshUtil::GetPolyData(Mesh& mesh_, std::vector<size_t> V) {
    CheckValidity();

    vtkNew<vtkPoints> points;
    vtkNew<vtkCellArray> cells;
    auto polyData = vtkSmartPointer<vtkPolyData>::New();
    
    std::vector<size_t> vertices;
    std::vector<size_t> faces;
    for (auto vid: V) {
        AddContents(faces, mesh_.V.at(vid).N_Fids);
    }
    for (auto fid: faces) {
        auto& f = mesh_.F.at(fid);
        if (f.N_Fids.empty() || f.Vids.empty()) continue;
        AddContents(vertices, f.Vids);
        vtkNew<vtkPolygon> p;
        p->GetPointIds()->SetNumberOfIds(f.Vids.size());
        for (int i = 0; i < f.Vids.size(); i++) {
            p->GetPointIds()->SetId(i, f.Vids.at(i));
        }
        cells->InsertNextCell(p);
    }
    for (auto vid: vertices) {
        auto& v = mesh_.V.at(vid);
        vtkIdType pId = v.id;
        points->InsertPoint(pId, v.x, v.y, v.z);
    }

    polyData->SetPoints(points);
    polyData->SetPolys(cells);
    return polyData;
}

void MeshUtil::SetMeshArea() {
    CheckValidity();

    double area = 0.0;
    for (auto& f: mesh->F) {
        area += GetFaceArea(f.id);
    }
    mesh->totalArea = area;
}

double MeshUtil::GetMeshArea() {
    CheckValidity();

    if (mesh->totalArea == 0.0) {
        SetMeshArea();
    }
    return mesh->totalArea;
}

double MeshUtil::GetFaceArea(int fid) {
    CheckValidity();

    Face& f = mesh->F.at(fid);
    glm::dvec3 a = mesh->V.at(f.Vids.at(0)).xyz() - mesh->V.at(f.Vids.at(2)).xyz();
    glm::dvec3 b = mesh->V.at(f.Vids.at(1)).xyz() - mesh->V.at(f.Vids.at(3)).xyz();

    return 0.5 * glm::length(glm::cross(a,b));
}

double MeshUtil::GetVertexEnergy(int vid) {
    CheckValidity();

    auto& v = mesh->V.at(vid);
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

    Vertex& v = mesh->V.at(vid);
    Edge& e = mesh->E.at(eid);
    int vid2 = e.Vids.at(0) == vid ? e.Vids.at(1) : e.Vids.at(0);
    std::vector<int> nvids;
    for (auto fid: e.N_Fids) {
        Face& f = mesh->F.at(fid);
        for (auto feid: f.Eids) {
            if (feid == eid) continue;
            auto& fe = mesh->E.at(feid);
            if (std::find(fe.Vids.begin(), fe.Vids.end(), vid) != fe.Vids.end()) {
                int nvid = fe.Vids.at(0) == vid ? fe.Vids.at(1) : fe.Vids.at(0);
                nvids.push_back(nvid);
            }
        }
    }

    Vertex& v_n = mesh->V.at(vid2);
    glm::dvec3 a = glm::normalize(v.xyz() - v_n.xyz());
    glm::dvec3 b = glm::normalize(mesh->V.at(nvids.at(0)).xyz() - v_n.xyz());
    glm::dvec3 c = glm::normalize(mesh->V.at(nvids.at(1)).xyz() - v_n.xyz());

    double alpha1 = acos(std::max(-1.0, std::min(glm::dot(a,b), 1.0))) * 180.0 / PI;
    double alpha2 = acos(std::max(-1.0, std::min(glm::dot(a,c), 1.0))) * 180.0 / PI;
    // std::cout << "alpha1: " << alpha1 << " alpha2: " << alpha2 << " total: " << alpha1 + alpha2 << std::endl;

    return alpha1 + alpha2;
}

bool MeshUtil::IsSharpFeature(size_t vid) {
    auto& v = mesh->V.at(vid);
    if (v.type == FEATURE || v.isBoundary) {
        int count = 0;
        for (auto nvid: v.N_Vids) {
            if (mesh->V.at(nvid).type == FEATURE || mesh->V.at(nvid).isBoundary) count += 1; 
        }
        if (count != 2) {
            mesh->SetIdealValence(vid);
            if (v.idealValence != 4) return true;
        }
    }
    return false;
}


std::vector<size_t> MeshUtil::GetDifference(std::vector<size_t> a, std::vector<size_t> b) {
    std::vector<size_t> diff;
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());
    std::set_difference(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(diff));
    return diff;
}

std::vector<size_t> MeshUtil::GetUnion(std::vector<size_t> a, std::vector<size_t> b) {
    std::vector<size_t> uni;
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());
    std::set_union(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(uni));
    return uni;
}

std::vector<size_t> MeshUtil::GetIntersection(std::vector<size_t> a, std::vector<size_t> b) {
    std::vector<size_t> itn;
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());
    std::set_intersection(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(itn));
    return itn;
}

std::vector<size_t> MeshUtil::GetIntersectionParallel(std::vector<size_t> a, std::vector<size_t> b) {
    std::vector<size_t> itn;
    for (auto el: a) {
        if (std::find(b.begin(), b.end(), el) != b.end()) itn.push_back(el);
    }
    return itn;
}

void MeshUtil::AddContents(std::vector<size_t>& a, std::vector<size_t> b) {
    std::vector<size_t> temp = a;
    temp = GetUnion(temp, b);
    a.clear();
    a.insert(a.begin(), temp.begin(), temp.end());
}

void MeshUtil::UpdateContents(std::vector<size_t>& a, std::vector<size_t> b) {
    std::vector<size_t> temp = a;
    temp = GetDifference(temp, b);
    a.clear();
    a.insert(a.begin(), temp.begin(), temp.end());
}