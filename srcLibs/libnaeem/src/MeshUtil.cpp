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

MeshUtil::MeshUtil(Mesh& mesh_) {
    mesh = &mesh_;
}

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
    // points->Allocate(mesh->V.size());
    vtkNew<vtkCellArray> cells;
    // cells->Allocate(cells->EstimateSize(mesh->F.size(), 4));
    auto polyData = vtkSmartPointer<vtkPolyData>::New();

    PARALLEL_FOR_BEGIN(0, mesh->V.size()) {
        auto& v = mesh->V.at(i);
        // vtkIdType pId = v.id;
        // points->SetPoint(pId, mesh->V.at(i).x, mesh->V.at(i).y, mesh->V.at(i).z);
        std::lock_guard<std::mutex> lock(Mutex);
        points->InsertPoint(v.id, v.x, v.y, v.z);
    } PARALLEL_FOR_END();
    // for (auto& v: mesh->V) {
    //     vtkIdType pId = v.id;
    //     points->InsertPoint(pId, v.x, v.y, v.z);
    // }
    PARALLEL_FOR_BEGIN(0, mesh->F.size()) {
        vtkNew<vtkPolygon> p;
        p->GetPointIds()->SetNumberOfIds(mesh->F.at(i).Vids.size());
        for (int j = 0; j < mesh->F.at(i).Vids.size(); j++) {
            p->GetPointIds()->SetId(j, mesh->F.at(i).Vids.at(j));
        }
        std::lock_guard<std::mutex> lock(Mutex);
        cells->InsertNextCell(p);
    } PARALLEL_FOR_END();
    // for (auto& c: mesh->C) {
    //     vtkNew<vtkPolygon> p;
    //     p->GetPointIds()->SetNumberOfIds(c.Vids.size());
    //     for (int i = 0; i < c.Vids.size(); i++) {
    //         p->GetPointIds()->SetId(i, c.Vids.at(i));
    //     }
    //     cells->InsertNextCell(p);
    // }
    
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
        if (v.isBoundary || count != 2) {
            mesh->SetIdealValence(vid);
            if ((v.isBoundary && v.idealValence != 2) || (v.type == FEATURE && v.idealValence != 4)) return true;
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

// function to set v_scores

void MeshUtil::SetV_Scores() {
    v_scores.clear();
    v_scores.resize(mesh->V.size());
}

// function to set qems of all vertices

void MeshUtil::SetQEMs() {
    SetV_Scores();
    for (auto& v: mesh->V) {
        CalculateQEM(v.id);
    }
}

// function to calculate QEM of a vertex

double MeshUtil::CalculateQEM(size_t vid) {
    // if (v_scores.size() < mesh->V.size()) v_scores.resize(mesh->V.size());
    auto& v = mesh->V.at(vid);
    // v_score->Q = glm::dmat4(0.0f);
    // v_scores.at(vid).Q = glm::dmat4(0.0f);
    auto& p = v.xyz();
    glm::dmat4 Q = glm::dmat4(0.0f);
    std::unordered_set<size_t> faces;
    for (auto fid: v.N_Fids) {
        faces.insert(fid);
        for (auto nfid: mesh->F.at(fid).N_Fids) {
            faces.insert(nfid);
        }
    }

    for (auto fid: faces) {
        auto& f = mesh->F.at(fid);
        if (f.Vids.empty()) continue;
        // auto& v1 = mesh->V.at(f.Vids.at(0)).xyz();
        // auto& v2 = mesh->V.at(f.Vids.at(1)).xyz();
        // auto& v3 = mesh->V.at(f.Vids.at(2)).xyz();
        // auto& v4 = mesh->V.at(f.Vids.at(3)).xyz();
        // auto& e1 = mesh->V.at(f.Vids.at(1)).xyz() - mesh->V.at(f.Vids.at(0)).xyz();
        // auto& e2 = mesh->V.at(f.Vids.at(2)).xyz() - mesh->V.at(f.Vids.at(1)).xyz();
        // auto& e3 = v1 - v3;
        // auto& e4 = v1 - v4;
        auto& n = glm::normalize(glm::cross((mesh->V.at(f.Vids.at(1)).xyz() - mesh->V.at(f.Vids.at(0)).xyz()), (mesh->V.at(f.Vids.at(2)).xyz() - mesh->V.at(f.Vids.at(1)).xyz())));

        // calculate QEM
        // glm::dmat4 q = glm::dmat4(0.0f);

        // Calculate the QEM matrix elements
        double a = glm::dot(n, n);
        double b = -2.0 * glm::dot(n, p);
        double c = glm::dot(p, p);

        // Add the matrix elements to the QEM matrix
        Q[0][0] += a*a;
        Q[0][1] += a*b;
        Q[0][2] += a*n[0]*c - b*n[0];
        Q[0][3] += a*n[1]*c - b*n[1];
        Q[1][1] += b*b;
        Q[1][2] += b*c;
        Q[2][2] += c*c;
        Q[3][0] += n[0]*c;
        Q[3][1] += n[1]*c;
        Q[3][2] += n[2]*c;

        // q[0] = glm::dvec4(n.x * n.x, n.x * n.y, n.x * n.z, -n.x * glm::dot(v1, n));
        // q[1] = glm::dvec4(n.y * n.x, n.y * n.y, n.y * n.z, -n.y * glm::dot(v1, n));
        // q[2] = glm::dvec4(n.z * n.x, n.z * n.y, n.z * n.z, -n.z * glm::dot(v1, n));
        // q[3] = glm::dvec4(-n.x * glm::dot(v1, n), -n.y * glm::dot(v1, n), -n.z * glm::dot(v1, n), glm::dot(v1, n) * glm::dot(v1, n));
        // q += glm::dmat4(glm::outerProduct(e1, e1))
        //   +  glm::dmat4(glm::outerProduct(e2, e2))
        //   +  glm::dmat4(glm::outerProduct(e3, e3))
        //   +  glm::dmat4(glm::outerProduct(e4, e4));

        // Q += q;
        // v_scores.at(vid).Q += q;
    }
    // v_scores.at(vid).qem_score = glm::dot(glm::dvec4(p.x, p.y, p.z, 1.0), v_scores.at(vid).Q * glm::dvec4(p.x, p.y, p.z, 1.0));
    // return glm::dot(glm::dvec4(p.x, p.y, p.z, 1.0), v_scores.at(vid).Q * glm::dvec4(p.x, p.y, p.z, 1.0));
    return glm::dot(glm::dvec4(p.x, p.y, p.z, 1.0), Q * glm::dvec4(p.x, p.y, p.z, 1.0));
}

// function to calculate QEM cost of a vertex

double MeshUtil::GetQEMcost(size_t vid) {
    /*auto& v = mesh->V.at(vid);
    auto& q = v_scores.at(vid).Q;
    auto& p = v.xyz();
    return glm::dot(glm::dvec4(p.x, p.y, p.z, 1.0), q * glm::dvec4(p.x, p.y, p.z, 1.0));*/
    return v_scores.at(vid).qem_score;
}

// function to get a pointer to VertexScore object for vid

VertexScore* MeshUtil::GetVertexScore(size_t vid, bool calculate_qem) {
    // if (v_scores.size() < mesh->V.size()) v_scores.resize(mesh->V.size());
    VertexScore v_score;
    if (calculate_qem) v_score.qem_score = CalculateQEM(vid);
    return &v_score;
    // return &v_scores.at(vid);
}

// function to get EdgeId given two vertex ids

size_t MeshUtil::GetEdgeId(size_t vid, size_t vid2) {
    size_t res;
    auto& v = mesh->V.at(vid);
    for (auto eid: v.N_Eids) {
        auto& e = mesh->E.at(eid);
        if (e.Vids.at(0) == vid2 || e.Vids.at(1) == vid2) return eid;
    }
    return res;
}

// function to get vertex index in face given an offset
int MeshUtil::GetVIdx(size_t vid, size_t fid) {
    auto& f = mesh->F.at(fid);
    return std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), vid));
}

// function to get vertex index in face

size_t MeshUtil::GetFaceV(size_t fid, int idx, int offset) {
	auto& f = mesh->F.at(fid);
	return f.Vids.at((idx+offset)%f.Vids.size());
}



// function to get Edge in order starting from a counter

size_t MeshUtil::GetCCedgeAt(size_t vid, size_t eid, int counter) {
    auto& v = mesh->V.at(vid);
    for (int i = 0; i < counter; i++) {
        auto& e = mesh->E.at(eid);
        size_t prev = e.Vids.at(0) == v.id ? e.Vids.at(1) : e.Vids.at(0);
        for (auto fid: e.N_Fids) {
            auto& f = mesh->F.at(fid);
            if (GetFaceV(f.id, GetVIdx(v.id, f.id), 1) == prev) {
                eid = GetDifference(GetIntersection(v.N_Eids, f.Eids), std::vector<size_t>{e.id}).at(0);
                break;
            }
        }
    }
    return eid;
}

// function to get face in order starting from a counter

size_t MeshUtil::GettCCFaceAt(size_t vid, size_t eid, int counter) {
    size_t res;
    auto& v = mesh->V.at(vid);
    for (int i = 0; i < counter; i++) {
        auto& e = mesh->E.at(eid);
        size_t prev = e.Vids.at(0) == v.id ? e.Vids.at(1) : e.Vids.at(0);
        for (auto fid: e.N_Fids) {
            auto& f = mesh->F.at(fid);
            if (GetFaceV(v.id, f.id, 1) == prev) {
                res = fid;
                eid = GetEdgeId(v.id, GetFaceV(f.id, GetVIdx(v.id, f.id), 3));
                break;
            }
        }
    }
    return res;
}

// function to get directions of a vertex successively

std::vector<size_t> MeshUtil::GetVertexDirections(size_t vid, std::vector<size_t> exceptions) {
    auto& v = mesh->V.at(vid);
    std::vector<size_t> neighbors;
    if (v.N_Vids.size() == 0) return neighbors;
    size_t eid = GetEdgeId(v.id, v.N_Vids.at(0));
    for (int i = 1; i <= v.N_Vids.size(); i++) {
        auto& f = mesh->F.at(GettCCFaceAt(v.id, eid, 1));
        int idx = GetVIdx(v.id, f.id);
        neighbors.push_back(GetFaceV(f.id, idx, 1));
        neighbors.push_back(GetFaceV(f.id, idx, 2));
        eid = GetCCedgeAt(v.id, eid, 1);
    }
    return GetDifference(neighbors, exceptions);
}