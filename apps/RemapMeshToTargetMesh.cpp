#include "Mesh.h"
#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "QuadSurfaceMapper.h"
#include "MeshUtil.h"
#include "SemiGlobalSimplifier.h"
#include "FeatureExtractor.h"

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

    FeatureExtractor fe(source, 20.0);
    fe.Extract();

    SemiGlobalSimplifier sg(source);
    sg.SetDirectSeparatrixOperations(); 
    // sg.SetSimplificationOperations();
    // sg.SetDiagonalCollapseOperations();
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
    std::cout << "# F in input mesh: " << source.C.size() << std::endl;
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
    std::cout <<  "# F in output mesh: " << source.C.size() << std::endl;

    // std::cout << "Writing output file" << std::endl;
    // std::ofstream ofs(output_f.c_str());
    // ofs << "# vtk DataFile Version 3.0\n"
    //     << output_f.c_str() << ".vtk\n"
    //     << "ASCII\n\n"
    //     << "DATASET UNSTRUCTURED_GRID\n";
    // ofs << "POINTS " << source.V.size() << " double\n";
    // std::vector<size_t> c_indices;
    // for (auto& v: source.V) {
    //     if (v.type == FEATURE) c_indices.push_back(v.id);
    // }
    // // std::vector<size_t> c_indices = {12, 296};
    // // std::cout << c_indices.size() << std::endl;
    // for (size_t i = 0; i < source.V.size(); i++) {
    //     ofs << std::fixed << std::setprecision(7) <<  source.V.at(i).x << " " <<  source.V.at(i).y << " " <<  source.V.at(i).z << "\n";
    // }
    // ofs << "CELLS " << c_indices.size() << " " << 2 * c_indices.size() << std::endl;
    // for (size_t i = 0; i < c_indices.size(); i++) {
    //     ofs << "1 " << c_indices.at(i) << std::endl;
    // }
    // ofs << "CELL_TYPES " << c_indices.size() << "\n";
    // for (size_t i = 0; i < c_indices.size(); i++) {
    //     ofs << "1" << std::endl;
    // }

    MeshFileWriter writer(source, output_f.c_str());
    writer.WriteFile();

    return 0;
}