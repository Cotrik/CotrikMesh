#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "BaseComplex.h"
#include "BaseComplexSheet.h"
#include "BaseComplexChord.h"
#include "BaseComplexEditor.h"
#include "SingularityGraph.h"
#include "ArgumentManager.h"

#include <iostream>

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cout << "Usage: CheckManifold quad_mesh" << std::endl;
        return -1;
    }
    ArgumentManager argumentManager(argc, argv);
    std::string filename = argv[1];
    std::cout << "---------------------------------------" << std::endl;
    std::cout << "filename = " << filename << std::endl;
    std::cout << "---------------------------------------" << std::endl;

    MeshFileReader reader(filename.c_str());
    Mesh& mesh = (Mesh&)reader.GetMesh();
    mesh.BuildAllConnectivities();
    mesh.ExtractBoundary();

    for (auto& e : mesh.E)
        if (e.N_Fids.size() > 2) {
            std::cerr << "None-manifold\n";
            return 0;
        }
    std::cerr << "Manifold\n";

    return 0;
}
