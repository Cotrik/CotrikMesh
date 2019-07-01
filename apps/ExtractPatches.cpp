/*
 * ExtractPatches.cpp
 *
 *  Created on: May 17, 2017
 *      Author: cotrik
 */
#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "FrameField.h"
#include "PolyLine.h"
#include "FrameOpt.h"
#include "FeatureLine.h"
#include "MeshQuality.h"
#include "Patches.h"
#include "Util.h"
#include "ArgumentManager.h"

const double PI = 3.1415926535;

std::vector<FACE_TYPE> GetFaceTypes(const Mesh& mesh) {
    std::vector<FACE_TYPE> faceTypes(mesh.F.size());
    for (unsigned long i = 0; i < mesh.F.size(); i++) {
        const Face& face = mesh.F.at(i);
        faceTypes.at(i) = Util::GetFaceType(face.normal);
    }
    return faceTypes;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cout << "Usage: ExtractPatches <input_tri_file> <output_tri_vtk_file> min_numOfFaces_per_patch=<10> localangle=<170> globalangle=<90>" << std::endl;
        return -1;
    }
    ArgumentManager argumentManager(argc, argv);
    std::string filename = argv[1];
    std::string output_filename = argv[2];
    int min_numOfFaces_per_patch = 10;
    double localangle = 170.0;
    double globalangle = 90.0;

    double coslocalangle = cos((180.0 - localangle) * PI / 180.0);
    const std::string strLocalAngle = argumentManager.get("localangle");
    if (!strLocalAngle.empty()) coslocalangle = cos((180.0 - std::stod(strLocalAngle)) * PI / 180.0);

    double cosglobalangle = cos((180.0 - globalangle) * PI / 180.0);
    const std::string strGlobalaAngle = argumentManager.get("globalangle");
    if (!strGlobalaAngle.empty()) cosglobalangle = cos((180.0 - std::stod(strGlobalaAngle)) * PI / 180.0);

    const std::string strmin_numOfFaces_per_patch = argumentManager.get("min_numOfFaces_per_patch");
    if (!strmin_numOfFaces_per_patch.empty()) min_numOfFaces_per_patch = std::stoi(strmin_numOfFaces_per_patch);

    std::cout << "---------------------------------------" << std::endl;
    std::cout << "filename = " << filename << std::endl;
    std::cout << "output_filename = " << output_filename << std::endl;
    std::cout << "min_numOfFaces_per_patch = " << min_numOfFaces_per_patch << std::endl;
    std::cout << "localangle = " << strLocalAngle << std::endl;
    std::cout << "globalangle = " << strGlobalaAngle << std::endl;
    std::cout << "---------------------------------------" << std::endl;

    MeshFileReader reader(filename.c_str());
    Mesh& mesh = (Mesh&)reader.GetMesh();
    mesh.RemoveUselessVertices();
    mesh.BuildAllConnectivities();
    mesh.ExtractBoundary();
//    mesh.ExtractLayers();
    mesh.ExtractSingularities();
    mesh.SetCosAngleThreshold(coslocalangle);
//    mesh.LabelSurface();
    mesh.GetNormalOfSurfaceFaces();
    mesh.GetNormalOfSurfaceVertices();
//    mesh.ExtractTwoRingNeighborSurfaceFaceIdsForEachVertex(3);
//    size_t genus =
//    std::cout << "genus = " <<  1 - (mesh.V.size() - mesh.E.size() + mesh.F.size() - mesh.C.size()) << std::endl;
//    std::cout << "#V:" << mesh.V.size() << " - #E:" << mesh.E.size() << " + #F:" << mesh.F.size() << " - #C:" << mesh.C.size() << " = " << "1 - genus" << std::endl;


    Patches patches(mesh);
    patches.SetGlobalCosAngle(cosglobalangle);
    if (!strGlobalaAngle.empty()) patches.Extract2();
    else patches.Extract();
    std::vector<size_t> edgeIds;
    for (size_t i = 0; i < patches.patches.size(); i++) {
        const Patch& patch = patches.patches.at(i);
        std::copy(patch.edgeIds.begin(), patch.edgeIds.end(), back_inserter(edgeIds));
    }
    patches.WriteMeshFile(output_filename.substr(0, output_filename.size() - 4).c_str());
    std::vector<size_t> faceids;
    for (auto& f : mesh.F)
        if (f.isBoundary) faceids.push_back(f.id);
    MeshFileWriter facesFileWriter(mesh, output_filename.c_str());
    facesFileWriter.WriteFacesVtk(faceids);
//    MeshFileWriter edgesFileWriter(mesh, "PatchesEdges.vtk");
//    edgesFileWriter.WriteEdgesVtk(edgeIds);

    std::vector<size_t> roundVertexIds;
    std::map<size_t, size_t> roundVertexIds_count;
    for (const Patch& patch : patches.patches) {
        for (auto eid : patch.edgeIds) {
            auto& e = mesh.E.at(eid);
            ++roundVertexIds_count[e.Vids[0]];
            ++roundVertexIds_count[e.Vids[1]];
        }
    }
    std::vector<std::vector<size_t>> sharpEdgeVertexIds;

    for (auto& v : mesh.V)
        if (v.type == CORNER || v.type == FEATURE) roundVertexIds.push_back(v.id);
    for (auto& p : roundVertexIds_count) {
        //std::cout << p.second << "\n";
        if (p.second >=6) roundVertexIds.push_back(p.first);
    }
    for (auto& eid : edgeIds) {
        auto& e = mesh.E.at(eid);
//        for (auto& fid : e.N_Fids) {
//            sharpEdgeVertexIds.push_back({fid, e.Vids[0]});
//            sharpEdgeVertexIds.push_back({fid, e.Vids[1]});
//        }
          sharpEdgeVertexIds.push_back({e.Vids[0], e.Vids[1]});
    }
    {
//        std::ofstream ofs("corners.txt");
//        for (auto c : roundVertexIds)
//            ofs << c << "\n";
        MeshFileWriter writer(mesh, "Corners.vtk");
        writer.WriteVerticesVtk(roundVertexIds);
    }
    {
//        std::ofstream ofs("sharpEdges.txt");
//        for (auto& e : sharpEdgeVertexIds)
//            ofs << e.front() << " " << e.back() << "\n";
    }
    return 0;
}




