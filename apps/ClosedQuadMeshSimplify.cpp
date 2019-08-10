/*
 * ClosedQuadMeshSimplify.cpp
 *
 *  Created on: June 19, 2019
 *      Author: cotrik
 */

#include "Simplifier.h"
#include "PatchSimplifier.h"
#include "ArgumentManager.h"

void setup(ArgumentManager& argumentManager, Simplifier& s);

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cout << "Usage: ClosedQuadMeshSimplify quad.vtk simplified.vtk iters=<1> maxValence=<3> maxValence=<5> "
            << "featurePreserved=true smoothIters=20 angle=<160> userCorners=\"\" canceledCorners=\"\" checkCorner=true" <<
            " collapse=false split=true conformal=true global=true rotate=true remove_doublet=true collapse_diagnal=true" <<
            " sheet_split=true half=false trip=true writeFile=true min_numOfFaces_per_patch=<10> localangle=<160> globalangle=<90>" << std::endl;
        return -1;
    }
    ArgumentManager argumentManager(argc, argv);
    std::string input = argv[1];
    std::string output = argv[2];
    int min_numOfFaces_per_patch = 10;
    double localangle = 160.0;
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
    std::cout << "input = " << input << std::endl;
    std::cout << "output = " << output << std::endl;
    std::cout << "min_numOfFaces_per_patch = " << min_numOfFaces_per_patch << std::endl;
    std::cout << "localangle = " << strLocalAngle << std::endl;
    std::cout << "globalangle = " << strGlobalaAngle << std::endl;
    std::cout << "---------------------------------------" << std::endl;

    MeshFileReader reader(input.c_str());
    Mesh& mesh = (Mesh&)reader.GetMesh();
    mesh.RemoveUselessVertices();
    mesh.BuildAllConnectivities();
    mesh.ExtractBoundary();
    mesh.ExtractSingularities();
    mesh.SetCosAngleThreshold(coslocalangle);
    mesh.GetNormalOfSurfaceFaces();
    mesh.GetNormalOfSurfaceVertices();

    Patches patches(mesh);
    patches.SetGlobalCosAngle(cosglobalangle);
    if (!strGlobalaAngle.empty()) patches.Extract2();
    else patches.Extract();
    std::vector<size_t> edgeIds;
    for (const Patch& patch : patches.patches)
        std::copy(patch.edgeIds.begin(), patch.edgeIds.end(), back_inserter(edgeIds));
    patches.WriteMeshFile("Patch");

    for (size_t j = 0; j < patches.patches.size(); j++) {
        auto patch_input = std::string("Patch") + std::to_string(j) + ".vtk";
        MeshFileReader reader(patch_input.c_str());
        Mesh& patch_mesh = (Mesh&) reader.GetMesh();
        patch_mesh.RemoveUselessVertices();
        PatchSimplifier simplifier(patch_mesh);
        setup(argumentManager, simplifier);
        simplifier.Run();
        int i = 2;
        while (i--/*mesh.V.size() < simplifier.origMesh.V.size()*/) {
            simplifier.init();
            patch_mesh = simplifier.RefineWithFeaturePreserved(patch_mesh, 0);
            size_t id = 0;
            for (auto& c : simplifier.mesh.C) {
                c.cellType = VTK_QUAD;
                c.id = id++;
            }
        }
        simplifier.init();
        simplifier.smooth_project(Simplifier::resolution);
        {
            auto patch_output = std::string("PatchOut") + std::to_string(j) + ".vtk";
            MeshFileWriter writer(patch_mesh, patch_output.c_str());
            writer.WriteFile();
            std::cout << "V = " << patch_mesh.V.size() << std::endl;
            std::cout << "E = " << patch_mesh.E.size() << std::endl;
            std::cout << "F = " << patch_mesh.F.size() << std::endl;
        }
    }

    return 0;
}

void setup(ArgumentManager& argumentManager, Simplifier& s) {
    argumentManager.get("iters", Simplifier::iters);
    argumentManager.get("maxValence", Simplifier::maxValence);
    argumentManager.get("minValence", Simplifier::minValence);
    argumentManager.get("smoothIters", Simplifier::smoothIters);
    argumentManager.get("resolution", Simplifier::resolution);

    argumentManager.get("angle", Simplifier::angle);

    argumentManager.get("featurePreserved", Simplifier::featurePreserved);
    argumentManager.get("collapse", Simplifier::COLLAPSE);
    argumentManager.get("split", Simplifier::SPLIT);
    argumentManager.get("half", Simplifier::HALF);
    argumentManager.get("trip", Simplifier::TRIP);
    argumentManager.get("rotate", Simplifier::ROTATE);
    argumentManager.get("remove_doublet", Simplifier::REMOVE_DOUBLET);
    argumentManager.get("collapse_diagnal", Simplifier::COLLAPSE_DIAGNAL);
    argumentManager.get("sheet_split", Simplifier::SHEET_SPLIT);
    argumentManager.get("conformal", Simplifier::CONFORMAL);
    argumentManager.get("global", Simplifier::GLOBAL);
    argumentManager.get("checkCorner", Simplifier::checkCorner);
    argumentManager.get("writeFile", Simplifier::writeFile);

    argumentManager.get("userCorners", Simplifier::userCorners);
    argumentManager.get("canceledCorners", Simplifier::canceledCorners);

    std::cout << "iters = " << Simplifier::iters << std::endl;
    std::cout << "maxValence = " << Simplifier::maxValence << std::endl;
    std::cout << "minValence = " << Simplifier::minValence << std::endl;
    std::cout << "smoothIters = " << Simplifier::smoothIters << std::endl;
    std::cout << "resolution = " << Simplifier::resolution << std::endl;

    std::cout << "angle = " << Simplifier::angle << std::endl;

    std::cout << "featurePreserved = " << Simplifier::featurePreserved << std::endl;
    std::cout << "collapse = " << Simplifier::COLLAPSE << std::endl;
    std::cout << "split = " << Simplifier::SPLIT << std::endl;
    std::cout << "rotate = " << Simplifier::ROTATE << std::endl;
    std::cout << "remove_doublet = " << Simplifier::REMOVE_DOUBLET << std::endl;
    std::cout << "collapse_diagnal = " << Simplifier::COLLAPSE_DIAGNAL << std::endl;
    std::cout << "sheet_split = " << Simplifier::SHEET_SPLIT << std::endl;
    std::cout << "conformal = " << Simplifier::CONFORMAL << std::endl;
    std::cout << "global = " << Simplifier::GLOBAL << std::endl;
    std::cout << "half = " << Simplifier::HALF << std::endl;
    std::cout << "trip = " << Simplifier::TRIP << std::endl;
    std::cout << "checkCorner = " << Simplifier::checkCorner << std::endl;
    std::cout << "writeFile = " << Simplifier::writeFile << std::endl;
    std::cout << "---------------------------------------" << std::endl;
}
