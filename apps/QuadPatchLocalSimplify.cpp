/*
 * QuadPatchLocalSimplify.cpp
 *
 *  Created on: October 7, 2020
 *      Author: https://github.com/naeem014
 */

#include "ArgumentManager.h"
#include "LocalSimplifier.h"


void setup(ArgumentManager& manager, Simplifier& simplifier);

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: QuadPatchSimplify quad.vtk simplified.vtk iters=<10000> " 
			<< "featurePreserved=true angle=<160>" << std::endl;
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
    LocalSimplifier simplifier(mesh);
    setup(argumentManager, simplifier);
    simplifier.Simplify();
    MeshFileWriter writer(mesh, output.c_str());
    writer.WriteFile();
    std::cout << "V = " << mesh.V.size() << std::endl;
    std::cout << "E = " << mesh.E.size() << std::endl;
    std::cout << "F = " << mesh.F.size() << std::endl;
    return 0;
}

void setup(ArgumentManager& manager, Simplifier& simplifier) {
    manager.get("iters", Simplifier::iters);
    manager.get("angle", Simplifier::angle);
    manager.get("featurePreserved", Simplifier::featurePreserved);

    std::cout << "iters = " << Simplifier::iters << std::endl;
    std::cout << "angle = " << Simplifier::angle << std::endl;
    std::cout << "featurePreserved = " << Simplifier::featurePreserved << std::endl;
}