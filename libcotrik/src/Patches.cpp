/*
 * Patches.cpp
 *
 *  Created on: May 17, 2017
 *      Author: cotrik
 */

#include "Patches.h"
#include "algorithm"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <memory>
#include "stdio.h"
using namespace std;

Patches::Patches(Mesh& mesh)
: mesh(mesh)
{
    // TODO Auto-generated constructor stub

}

Patches::~Patches()
{
    // TODO Auto-generated destructor stub
}

void Patches::Extract()
{
    mesh.LabelSurface();
    ExtractPatches();
}

void Patches::Extract2()
{
    this->LabelSurface();
    ExtractPatches();
}

void Patches::ExtractPatches()
{
    const Patch temp(mesh);
    patches.resize(mesh.numberOfPatches, temp);
    for (size_t i = 0; i < mesh.F.size(); i++) {
        const Face& face = mesh.F.at(i);
        if (!face.isBoundary) continue;
        patches.at(face.label).faceIds.push_back(face.id);
    }

    for (size_t p = 0; p < patches.size(); p++) {
        Patch& patch = patches.at(p);
        std::vector<size_t>& faceIds = patch.faceIds;
        std::vector<size_t>& edgeIds = patch.edgeIds;
        std::vector<size_t>& vertexIds = patch.vertexIds;

        faceIds.resize(faceIds.size());
        std::sort(faceIds.begin(), faceIds.end());

        for (size_t i = 0; i < faceIds.size(); i++) {
            const Face& face = mesh.F.at(faceIds.at(i));
            for (size_t j = 0; j < face.Eids.size(); j++) {
                const Edge& edge = mesh.E.at(face.Eids.at(j));
                int count = 0;
                for (size_t k = 0; k < edge.N_Fids.size(); k++) {
                    const size_t n_face_id = edge.N_Fids.at(k);
                    const Face& n_face = mesh.F.at(n_face_id);
                    if (!n_face.isBoundary) continue;
                    const bool found = std::find(faceIds.begin(), faceIds.end(), n_face_id) != faceIds.end();
                    if (found) ++count;
                }
                if (count == 1)
                    edgeIds.push_back(edge.id);
            }

            for (size_t j = 0; j < face.Vids.size(); j++) {
                const Vertex& vertex = mesh.V.at(face.Vids.at(j));
                vertexIds.push_back(vertex.id);
            }
        }
        std::sort(edgeIds.begin(), edgeIds.end());

        std::sort(vertexIds.begin(), vertexIds.end());
        const std::vector<size_t>::iterator iter = std::unique(vertexIds.begin(), vertexIds.end());
        vertexIds.resize(std::distance(vertexIds.begin(), iter));
    }
}

void Patch::LabelFace(Face& initialFace, Face& face, size_t& label)
{
    face.label = label;
    for (size_t i = 0; i < face.Eids.size(); i++) {
        const Edge& edge = mesh.E.at(face.Eids.at(i));
        std::vector<Face*> faces;
        for (size_t j = 0; j < edge.N_Fids.size(); j++) {
            Face& face2 = mesh.F.at(edge.N_Fids.at(j));
            if (std::find(faceIds.begin(), faceIds.end(), face2.id) != faceIds.end())
            if (face2.isBoundary && face2.id != face.id && face2.label == MAXID) {
                const double cos_angle = glm::dot(initialFace.normal, face2.normal);//mesh.GetCosAngle(edge, initialFace, face2);
                //std::cout << "cos_angle = " << cos_angle << std::endl;
                if (cos_angle > cosangle/*mesh.cos_angle_threshold*/) // cos(15) = 0.9659 cos(30) = 0.866
                faces.push_back(&face2);
            }
        }
        for (size_t i = 0; i < faces.size(); i++)
            LabelFace(initialFace, *faces.at(i), label);
    }
}

void Patches::LabelFace(Face& initialFace, Face& face, size_t& label)
{
    face.label = label;
    for (size_t i = 0; i < face.Eids.size(); i++) {
        const Edge& edge = mesh.E.at(face.Eids.at(i));
        std::vector<Face*> faces;
        for (size_t j = 0; j < edge.N_Fids.size(); j++) {
            Face& face2 = mesh.F.at(edge.N_Fids.at(j));
            if (face2.isBoundary && face2.id != face.id && face2.label == MAXID) {
                const double cos_angle = glm::dot(initialFace.normal, face2.normal);//mesh.GetCosAngle(edge, initialFace, face2);
                //std::cout << "cos_angle = " << cos_angle << std::endl;
                if (cos_angle > cosangle/*mesh.cos_angle_threshold*/) // cos(15) = 0.9659 cos(30) = 0.866
                faces.push_back(&face2);
            }
        }
        for (size_t i = 0; i < faces.size(); i++)
            LabelFace(initialFace, *faces.at(i), label);
    }
}

void Patches::LabelSurface()
{
    mesh.LabelSurface();
    // Sort the face according to the normal
    struct normal_faceid{
        normal_faceid(){};
        normal_faceid(const normal_faceid& rhs):n(rhs.n), id(rhs.id){};
        glm::vec3 n;
        size_t id = 0;
        bool operator < (const normal_faceid& rhs) const {
            double max = fabs(n.x);
            if (fabs(n.y) > max) max = fabs(n.y);
            if (fabs(n.z) > max) max = fabs(n.z);

            double max_rhs = fabs(rhs.n.x);
            if (fabs(rhs.n.y) > max_rhs) max_rhs = fabs(rhs.n.y);
            if (fabs(rhs.n.z) > max_rhs) max_rhs = fabs(rhs.n.z);

            return (max > max_rhs);
        }
    };
//    normal_faceid temp;
//    std::vector<normal_faceid> n_fids(mesh.F.size(), temp);
//    for (size_t i = 0; i < mesh.F.size(); i++) {
//        n_fids[i].n = mesh.F[i].normal;
//        n_fids[i].id = mesh.F[i].id;
//    }
//    std::sort(n_fids.begin(), n_fids.end());
//    size_t label = 0;
//    for (size_t i = 0; i < n_fids.size(); i++) {
//        Face& face = mesh.F.at(n_fids.at(i).id);
//        if (!face.isBoundary || face.label != MAXID)
//            continue;
//        LabelFace(face, face, label);
//        label++;
//    }
//    mesh.numberOfPatches = label;

    const Patch temp(mesh);
    patches.resize(mesh.numberOfPatches, temp);
    for (size_t i = 0; i < mesh.F.size(); i++) {
        const Face& face = mesh.F.at(i);
        if (!face.isBoundary) continue;
        patches.at(face.label).faceIds.push_back(face.id);
    }

    size_t label = 0;
    for (size_t p = 0; p < patches.size(); p++) {
        Patch& patch = patches.at(p);
        std::vector<size_t>& faceIds = patch.faceIds;
        std::vector<size_t>& edgeIds = patch.edgeIds;
        std::vector<size_t>& vertexIds = patch.vertexIds;

        patch.SetGlobalCosAngle(cosangle);

        normal_faceid temp;
        std::vector<normal_faceid> n_fids;
        normal_faceid n_fid;
        for (size_t i = 0; i < faceIds.size(); i++) {
            n_fid.n = mesh.F[faceIds[i]].normal;
            n_fid.id = mesh.F[faceIds[i]].id;
            n_fids.push_back(n_fid);
        }
        std::sort(n_fids.begin(), n_fids.end());

        for (size_t i = 0; i < n_fids.size(); i++) {
            Face& face = mesh.F.at(n_fids.at(i).id);
            face.label = MAXID;
        }

        for (size_t i = 0; i < n_fids.size(); i++) {
            Face& face = mesh.F.at(n_fids.at(i).id);
            if (!face.isBoundary || face.label != MAXID)
                continue;
            patch.LabelFace(face, face, label);
            label++;
        }
    }
    mesh.numberOfPatches = label;
    patches.clear();
}

void Patch::WriteMeshFile(const char* filename) const {
    const std::vector<Vertex>& V = mesh.V;
    const std::vector<Face>& F = mesh.F;

    std::ofstream ofs(filename);
    ofs << "# vtk DataFile Version 2.0\n"
        << filename << "\n"
        << "ASCII\n\n"
        << "DATASET POLYDATA\n" ;
    ofs << "POINTS " << V.size() << " float\n";
    ofs << std::fixed << setprecision(7);
    for (size_t i = 0; i < V.size(); i++)
        ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << "\n";
    if (mesh.m_cellType == HEXAHEDRA || mesh.m_cellType == QUAD) {
        ofs << "POLYGONS " << faceIds.size() << " " << 5 * faceIds.size() << std::endl;
        for (auto i : faceIds) {
            ofs << "4";
            const Face& face = F.at(i);
            for (size_t j = 0; j < face.Vids.size(); j++)
                ofs << " " << face.Vids.at(j);
            ofs << "\n";
        }
    }
    else if (mesh.m_cellType == TETRAHEDRA || mesh.m_cellType == TRIANGLE) {
        ofs << "POLYGONS " << faceIds.size() << " " << 4 * faceIds.size() << std::endl;
        for (auto i : faceIds) {
            ofs << "3";
            const Face& face = F.at(i);
            for (size_t j = 0; j < face.Vids.size(); j++)
                ofs << " " << face.Vids.at(j);
            ofs << "\n";
        }
    }
    {
        ofs << "CELL_DATA " << faceIds.size() << std::endl
            << "SCALARS " << " Label" << " int 1" << std::endl
            << "LOOKUP_TABLE default" << std::endl;
        for (auto i : faceIds)
            ofs << (F.at(i).label == MAXID ? 0 : F.at(i).label)  << std::endl;
    }
}

void Patches::WriteMeshFile(const char* filename_prefix) const {
    size_t id = 0;
    for (const auto& p : patches)
        p.WriteMeshFile((string(filename_prefix) + to_string(id++) + ".vtk").c_str());
}
