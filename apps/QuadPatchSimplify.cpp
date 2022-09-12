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
    simplifier.smoothGlobal = true;
    simplifier.SmoothMesh();
    // simplifier.smoothGlobal = true;
    // simplifier.SmoothMesh(true);
    simplifier.RefineMesh();
    simplifier.init();
    // SurfaceMapper sm(mesh, simplifier.origMesh);
    // for (auto& v: mesh.V) {
    //     v = sm.MapPoint(v.xyz());
    // }
    simplifier.smoothGlobal = true;
    simplifier.SmoothMesh();
    // simplifier.smooth_project(2);
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
// reproduce the results
// 5 5 split
// collapse causes overlaps
// write the paper
// make the illustrations
// start from siggraph asia
// methodology part can be extended
// results and discussion

// Introduction: Decide what style we want for introduction. Mathematical notations move to later section in overview.
// what is quad mesh and why is it important, problem why we need to simplify it, talk about issues of existing methods
// because of those limitations we propose semi-global simplification, integrate with existing local operators and then we integrate
// a unified operation and decide which operations to perform and achieve boundary configuration. describe briefly.
// optional paragraph: the rest of the paper is structured as follows.

// Related work: Quad mesh generation, cite a lot of papers. Talk about techniques for mesh generation. emphasize the issues of some of
// quad mesh generation. Quad layout in some papers. Structure simplification of quad mesh. Add additional papers if we can find more. 

// Background: how separatrices relate to poincare index, talk about valence, discrete index of each singularity. 

// Our pipeline contains two types of operation; semi-global and local and integrate them into once complete piepline. 
// 5-5 operator is optional which user can turn on or off. In our experiments we turn it off. 
// local smoothing optional every other step.
// Evaluation with certain things turned on or off
// Impact of ranking the operations.