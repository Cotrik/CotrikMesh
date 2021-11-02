#include "Mesh.h"
#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "QuadSurfaceMapper.h"
#include "MeshUtil.h"
#include "SemiGlobalSimplifier.h"

int main(int argc, char* argv[]) {
    std::string source_f = argv[1];
    // std::string target_f = argv[2];
    std::string output_f = argv[2];

    MeshFileReader source_reader(source_f.c_str());
    Mesh& source = (Mesh&) source_reader.GetMesh();
    source.RemoveUselessVertices();
    source.BuildAllConnectivities();
    // for (auto& el: source.F) {
    //     std:cout << el.N_Fids.size() << std::endl;
    // }


    // for (auto& e: source.E) {
    //     std::cout << e.N_Fids.size() << std::endl;
    // }
    SemiGlobalSimplifier sg(source);
    sg.SetSimplificationOperations();
    // for (auto& v: source.V) {
    //     // vu.GetVertexEnergy(v.id);
    //     std::cout << mu.GetVertexEnergy(v.id) << std::endl;
    // }

    // MeshFileReader target_reader(target_f.c_str());
    // Mesh& target = (Mesh&) target_reader.GetMesh();
    // target.RemoveUselessVertices();
    // target.BuildAllConnectivities();

    // SurfaceMapper sm(target);
    // sm.Map();

    // source.RemoveUselessVertices();
    // std::cout << source.C.size() << std::endl;
    // std::cout << source.F.size() << std::endl;
    std::vector<Cell> newC;
    for (auto& f: source.F) {
        if (f.N_Fids.size() == 0) continue;
        Cell c;
        c.id = newC.size();
        c.Vids = f.Vids;
        newC.push_back(c);
    }
    source.E.clear();
    source.F.clear();
    source.C.clear();
    source.C.insert(source.C.begin(), newC.begin(), newC.end());
    // source.BuildAllConnectivities();
    // std::cout << newC.size() << std::endl;
    std::cout << source.C.size() << std::endl;

    MeshFileWriter writer(source, output_f.c_str());
    writer.WriteFile();

    return 0;
}