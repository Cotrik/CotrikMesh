/*
 * QuadPatchSimplify.cpp
 *
 *  Created on: Dec 27, 2018
 *      Author: cotrik
 */

#include "Simplifier.h"
#include "PatchSimplifier.h"
#include "ArgumentManager.h"
#include "AngleBasedSmoothQuadMesh.h"
#include <ctime>

void setup(ArgumentManager& argumentManager, Simplifier& s);

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: QuadPatchSimplify quad.vtk simplified.vtk iters=<10000> maxValence=<3> maxValence=<6> " 
			<< "featurePreserved=true smoothIters=20 angle=<160> userCorners=\"\" canceledCorners=\"\" checkCorner=true" <<
			" collapse=true split=false conformal=true global=true rotate=true remove_doublet=true collapse_diagnal=true" <<
			" sheet_split=true half=false trip=false writeFile=false" << std::endl;
        return -1;
    }
	
    ArgumentManager argumentManager(argc, argv);
    std::string input = argv[1];
    std::string output = argv[2];

    std::cout << "---------------------------------------" << std::endl;
    std::cout << "input  = " << input << std::endl;
    std::cout << "output = " << output << std::endl;

    MeshFileReader reader(input.c_str());
    Mesh& mesh = (Mesh&) reader.GetMesh();
    mesh.RemoveUselessVertices();
    PatchSimplifier simplifier(mesh);
    setup(argumentManager, simplifier);

    std::clock_t start;
    double duration;
    start = std::clock();

    simplifier.Run();
	{
		MeshFileWriter writer(mesh, "VertexFeature.vtk");
		writer.WriteVertexFeatureVtk();
	}
	simplifier.init();

    duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
    std::cout << "Simplification time: " << duration << " seconds" << std::endl;
    // mesh.smoothGlobal = true;
    // while (simplifier.findCrossQuads()) {
    //     simplifier.smooth_project(2);
    // }
	simplifier.smooth_project();
	{
        // simplifier.RefineMesh();
        // SmoothAlgorithm smoothAlgo(mesh, simplifier.origMesh, 1000, 1, true, true);
        // smoothAlgo.smoothMesh();

		MeshFileWriter writer(mesh, output.c_str());
		writer.WriteFile();
        for (auto& f: mesh.F) {
            if (f.Vids.size() < 4) {
                std::cout << "Face " << f.id << " has less than 4 vertices" << std::endl;
            }
        }
		std::cout << "V = " << mesh.V.size() << std::endl;
		std::cout << "E = " << mesh.E.size() << std::endl;
		std::cout << "F = " << mesh.F.size() << std::endl;
	}
    return 0;
}

void setup(ArgumentManager& argumentManager, Simplifier& s) {
    argumentManager.get("iters", Simplifier::iters);
    argumentManager.get("maxValence", Simplifier::maxValence);
    argumentManager.get("minValence", Simplifier::minValence);
    argumentManager.get("smoothIters", Simplifier::smoothIters);

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


// TO DO:
// Fix 5 5 split
// 3 5 chord collapsing
