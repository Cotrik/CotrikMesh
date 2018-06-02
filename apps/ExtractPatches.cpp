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
#include "ArgumentManager.h"
#include <iostream>

const double PI = 3.1415926535;

enum FACE_TYPE
{
    FACE_UNKNOWN = 0,
    FACE_X = 1,
    FACE_Y = 2,
    FACE_Z = 3,
};

FACE_TYPE GetFaceType(const glm::vec3& normal)
{
    const float l = glm::length(normal);
    if (fabs(normal.x/l) > 0.8)
        return FACE_X;
    else if (fabs(normal.y/l) > 0.8)
        return FACE_Y;
    else if (fabs(normal.z/l) > 0.8)
        return FACE_Z;

    FACE_TYPE faceType = FACE_X;
    double max = fabs(normal.x / l);
    if (fabs(normal.y / l) > max) {
        max = fabs(normal.y / l);
        faceType = FACE_Y;

    }
    if (fabs(normal.z / l) > max) {
        max = fabs(normal.z / l);
        faceType = FACE_Z;
    }
    return faceType;
}

void GetFaceTypes(const Mesh& mesh, std::vector<FACE_TYPE>& faceTypes)
{
    faceTypes.resize(mesh.F.size());
    for (unsigned long i = 0; i < mesh.F.size(); i++) {
        const Face& face = mesh.F.at(i);
        faceTypes.at(i) = GetFaceType(face.normal);
    }
}

int main(int argc, char* argv[])
{
    if (argc < 3)
    {
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
    MeshFileWriter facesFileWriter(mesh, output_filename.c_str());
    facesFileWriter.WriteFacesVtk();
    MeshFileWriter edgesFileWriter(mesh, "PatchesEdges.vtk");
    edgesFileWriter.WriteEdgesVtk(edgeIds);
}




