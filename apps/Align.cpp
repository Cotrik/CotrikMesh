/*
 * Align.cpp
 *
 *  Created on: Nov 19, 2017
 *      Author: cotrik
 */


#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "FrameField.h"
#include "PolyLine.h"
#include "FrameOpt.h"
#include "FeatureLine.h"
#include "MeshQuality.h"
#include "ArgumentManager.h"
#include <iostream>
#include <iomanip>
#include <set>
#include <unordered_set>
#include <algorithm>

void align(const Mesh& orig_mesh, const Mesh& polycube_mesh, Mesh& aligned_polycube_mesh);
void align(const Mesh& orig_mesh, Mesh& polycube_mesh, const int smoothIters = 100);

void WriteVtkFile(const Mesh& m_mesh, const string&  m_strFileName);

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: Align orig=<orig.tet.vtk> polycube=<polycube.tet.vtk> output=<polycube.aligned.tet.vtk> smoothIters=<100> angle=<5|10|15|20|25|30>\n";
        return -1;
    }
    ArgumentManager argumentManager(argc, argv);
    string orig = "orig.tet.vtk";
    string polycube = "polycube.tet.vtk";
    string output = "polycube.aligned.tet.vtk";
    int smoothIters = 100;
    double angle = 10.0;
    {
        const std::string strorig = argumentManager.get("orig");
        if (!strorig.empty()) orig = strorig;

        const std::string strpolycube = argumentManager.get("polycube");
        if (!strpolycube.empty()) polycube = strpolycube;

        const std::string stroutput = argumentManager.get("output");
        if (!stroutput.empty()) output = strpolycube;

        const std::string strsmoothIters = argumentManager.get("smoothIters");
        if (!strsmoothIters.empty()) smoothIters = std::stoi(strsmoothIters);

        const std::string strAngle = argumentManager.get("angle");
        if (!strAngle.empty()) angle = std::stod(strAngle);
    }

    MeshFileReader orig_reader(orig.c_str());
    Mesh& orig_mesh = (Mesh&)orig_reader.GetMesh();
    orig_mesh.BuildAllConnectivities();
    orig_mesh.ExtractBoundary();
    orig_mesh.SetCosAngleThreshold(cos(angle));
    orig_mesh.LabelSurface();
    orig_mesh.LabelSharpEdges(true);
    orig_mesh.ExtractSingularities();

    MeshFileReader polycube_reader(polycube.c_str());
    Mesh& polycube_mesh = (Mesh&)polycube_reader.GetMesh();
    polycube_mesh.BuildAllConnectivities();
    polycube_mesh.ExtractBoundary();
    polycube_mesh.SetCosAngleThreshold(cos(angle));
    polycube_mesh.LabelSurface();
    polycube_mesh.LabelSharpEdges(true);
    polycube_mesh.ExtractSingularities();

    //Mesh aligned_polycube_mesh;
    //align(orig_mesh, polycube_mesh, aligned_polycube_mesh);
    //MeshFileWriter writer(aligned_polycube_mesh.V, polycube_mesh.C, output.c_str(), TETRAHEDRA);
    //MeshFileWriter writer(aligned_polycube_mesh, output.c_str());
    //writer.WriteFile();
    //WriteVtkFile(aligned_polycube_mesh, output);

    align(orig_mesh, polycube_mesh, smoothIters);
    MeshFileWriter writer(polycube_mesh, output.c_str());
    writer.WriteFile();

    return 0;
}
std::unordered_set<size_t> getSharpEdgeVids(const Mesh& mesh) {
    std::unordered_set<size_t> res;
    for (const auto& e: mesh.E) {
        if (!e.isSharpFeature) continue;
        if (res.find(e.Vids[0]) == res.end()) res.insert(e.Vids[0]);
        if (res.find(e.Vids[1]) == res.end()) res.insert(e.Vids[1]);
    }
    return res;
}
std::unordered_set<size_t> getNeighborVIds(const Mesh& mesh, const Vertex& v, int neighbors, std::unordered_set<size_t>& orig_sharp_vids) {
    std::unordered_set<size_t> res;
    res.insert(v.id);
    while (neighbors--) {
        auto n = res.size();
        std::unordered_set<size_t> rescopy = res;
        for (auto i : rescopy) {
            for (auto neiborVid : mesh.V.at(i).N_Vids) {
                const auto& vi = mesh.V.at(neiborVid);
                if (vi.isCorner) return {vi.id};
                if (!vi.isBoundary) continue;
                //if (orig_sharp_vids.find(neiborVid) == orig_sharp_vids.end()) continue;

                if (res.find(vi.id) == res.end()) res.insert(vi.id);
            }
        }
    }
    std::vector<size_t> rem;
    for (auto vid : res)
        if (orig_sharp_vids.find(vid) == orig_sharp_vids.end()) rem.push_back(vid);
    for (auto vid : rem) res.erase(vid);
    return res;
}

void align(const Mesh& orig_mesh, const Mesh& polycube_mesh, Mesh& aligned_polycube_mesh) {
    auto orig_sharp_vids = getSharpEdgeVids(orig_mesh);
    aligned_polycube_mesh = polycube_mesh;
    // Align corners;
    for (auto& v : aligned_polycube_mesh.V) {
        if (!v.isCorner) continue;
        if (orig_mesh.V.at(v.id).isCorner) continue;  // already align
        auto neighborVids = getNeighborVIds(orig_mesh, v, 2, orig_sharp_vids);
        {
            double min_distance = INT_MAX;
            size_t newCornerId = *neighborVids.begin();
            for (auto neighborVid : neighborVids) {
                auto& neighborV = aligned_polycube_mesh.V.at(neighborVid);
                double distance = glm::length(neighborV.xyz() - v.xyz());
                if (distance < min_distance) {
                    min_distance = distance;
                    newCornerId = neighborVid;
                }
            }
            v.isCorner = false;
            auto& newCorner = aligned_polycube_mesh.V.at(newCornerId);
            newCorner.isCorner = true;
            newCorner = v.xyz();
        }
    }
    // Align sharp edges;
    int smooth_iters = 100;
    while (smooth_iters--) {
        for (auto vid : orig_sharp_vids) {
            auto& v = aligned_polycube_mesh.V.at(vid);
            if (v.isCorner) continue;
            //if (orig_sharp_vids.find(vid) == orig_sharp_vids.end()) continue;
            std::vector<size_t> neighborVids(2, 0);
            int count = 0;
            for (auto neighborVid : v.N_Vids)
                if (orig_sharp_vids.find(neighborVid) != orig_sharp_vids.end()) neighborVids[count++] = neighborVid;
            auto& v1 = aligned_polycube_mesh.V.at(neighborVids[0]);
            auto& v2 = aligned_polycube_mesh.V.at(neighborVids[1]);
            v = (v1.xyz() + v2.xyz()) * 0.5f;
        }
    }
    // Align Surface;
}

void align(const Mesh& orig_mesh, Mesh& polycube_mesh, const int smoothIters) {
    auto orig_sharp_vids = getSharpEdgeVids(orig_mesh);
    auto& aligned_polycube_mesh = polycube_mesh;
    // Align corners;
    for (auto& v : aligned_polycube_mesh.V) {
        if (!v.isCorner) continue;
        if (orig_mesh.V.at(v.id).isCorner) continue;  // already align
        auto neighborVids = getNeighborVIds(orig_mesh, v, 2, orig_sharp_vids);
        {
            double min_distance = INT_MAX;
            size_t newCornerId = *neighborVids.begin();
            for (auto neighborVid : neighborVids) {
                auto& neighborV = aligned_polycube_mesh.V.at(neighborVid);
                double distance = glm::length(neighborV.xyz() - v.xyz());
                if (distance < min_distance) {
                    min_distance = distance;
                    newCornerId = neighborVid;
                }
            }
            v.isCorner = false;
            auto& newCorner = aligned_polycube_mesh.V.at(newCornerId);
            newCorner.isCorner = true;
            newCorner = v.xyz();
        }
    }
    // Align sharp edges;
    std::vector<FeatureLine> featureLines(orig_mesh.numOfSharpEdges, FeatureLine(orig_mesh));
    for (size_t i = 0; i < orig_mesh.numOfSharpEdges; i++)
        featureLines.at(i).Extract(i);
    int smooth_iters = smoothIters;
    while (smooth_iters--) {
        for (auto featureLine : featureLines) {
            if (featureLine.Vids.front() != featureLine.Vids.back())
            for (auto i = 1; i < featureLine.Vids.size() - 1; ++i) {
                auto vid = featureLine.Vids.at(i);
                auto& v = aligned_polycube_mesh.V.at(vid);
                if (v.isCorner) continue;
                auto& v1 = aligned_polycube_mesh.V.at(featureLine.Vids.at(i - 1));
                auto& v2 = aligned_polycube_mesh.V.at(featureLine.Vids.at(i + 1));
                v = (v1.xyz() + v2.xyz()) * 0.5f;
            }
            else
                for (auto i = 0; i < featureLine.Vids.size(); ++i) {
                    auto vid = featureLine.Vids.at(i);
                    auto& v = aligned_polycube_mesh.V.at(vid);
                    if (v.isCorner) continue;
                    auto& v1 = i == 0 ? aligned_polycube_mesh.V.at(featureLine.Vids.back()) : aligned_polycube_mesh.V.at(featureLine.Vids.at(i - 1));
                    auto& v2 = aligned_polycube_mesh.V.at(featureLine.Vids.at((i + 1) % featureLine.Vids.size()));
                    v = (v1.xyz() + v2.xyz()) * 0.5f;
                }
        }
    }
//    int smooth_iters = 10;
//    while (smooth_iters--) {
//        for (auto vid : orig_sharp_vids) {
//            auto& v = aligned_polycube_mesh.V.at(vid);
//            if (v.isCorner) continue;
//            if (orig_sharp_vids.find(vid) == orig_sharp_vids.end()) continue;
//            std::vector<size_t> neighborVids(2, 0);
//            int count = 0;
//            for (auto neighborVid : v.N_Vids)
//                if (orig_sharp_vids.find(neighborVid) != orig_sharp_vids.end()) neighborVids[count++] = neighborVid;
//            auto& v1 = aligned_polycube_mesh.V.at(neighborVids[0]);
//            auto& v2 = aligned_polycube_mesh.V.at(neighborVids[1]);
//            v = (v1.xyz() + v2.xyz()) * 0.5f;
//        }
//    }

    // Align Surface;
}

void WriteVtkFile(const Mesh& m_mesh, const string&  m_strFileName)
{
    const std::vector<Vertex>& V = m_mesh.V;
    const std::vector<Cell>& C = m_mesh.C;
    const size_t vnum = V.size();
    const size_t cnum = C.size();

    std::ofstream ofs(m_strFileName.c_str());
    ofs << "# vtk DataFile Version 3.0" << endl
        << m_strFileName.c_str() << endl
        << "ASCII" << endl << endl
        << "DATASET UNSTRUCTURED_GRID" << endl;
    ofs << "POINTS " << vnum << " double" << endl;
    for (size_t i = 0; i < vnum; i++)
        ofs << std::fixed << setprecision(7) << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << endl;
    ofs << "CELLS " << cnum << " " << 5*cnum << std::endl;

    for (size_t i = 0; i < cnum; i++){
        ofs << C.at(i).Vids.size();
        for (size_t j = 0; j < C.at(i).Vids.size(); j++)
            ofs << " " << C.at(i).Vids.at(j);
        ofs << std::endl;
    }
    ofs << "CELL_TYPES " << cnum << endl;
    for (size_t i = 0; i < cnum; i++)
        ofs << 10 << std::endl;
}
