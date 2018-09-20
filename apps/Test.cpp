/*
 * Test.cpp
 *
 *  Created on: Jan 8, 2018
 *      Author: cotrik
 */

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "BaseComplex.h"
#include <iostream>
#include "ArgumentManager.h"

int TotalVertexValenceOfSurface(const BaseComplex& bc) {
    int sum = 0;
    for (auto& sv : bc.SingularityI.V) {
        if (!sv.isBoundary) continue;
        int count = 0;
        for (auto vid : bc.mesh.V.at(sv.id_mesh).N_Vids)
            if (bc.mesh.V.at(vid).isBoundary) ++count;
        sum += 4 - count;
    }
    return sum;
}

int TotalEdgeValenceOfSurface(const BaseComplex& bc) {
    int sum = 0;
    for (auto& se : bc.SingularityI.E) {
        if (!se.isBoundary) continue;
        sum += 2 - bc.mesh.E.at(se.es_link.front()).N_Cids.size();
    }
    return sum;
}

int TotalEdgeValenceOfInterior(const BaseComplex& bc) {
    int sum = 0;
    for (auto& se : bc.SingularityI.E) {
        if (se.isBoundary) continue;
        sum += 4 - bc.mesh.E.at(se.es_link.front()).N_Cids.size();
    }
    return sum;
}

int Genus(const Mesh& mesh) {
    return  1 - (mesh.V.size() - mesh.E.size() + mesh.F.size() - mesh.C.size());
}

int main(int argc, char* argv[]){
    if (argc < 2) {
        std::cout << "Usage: Test hex.vtk" << std::endl;
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
    mesh.ExtractSingularities();
    //"Info: cos10 = 0.984807753; cos15 = 0.965925826; cos20 = 0.939692621; cos25 = 0.906307787; cos30 = 0.866025404\n\n";
    const double cosangle = 0.866025404;
    mesh.SetCosAngleThreshold(cosangle);
    mesh.LabelSurface();
    mesh.LabelSharpEdges(true);
    // For extracting singularity Graph
    mesh.BuildParallelE();
    mesh.BuildConsecutiveE();
    mesh.BuildOrthogonalE();

    BaseComplex baseComplex(mesh);
    baseComplex.Build();

    std::cerr << "Singularities " << Genus(mesh) << "\t" << TotalVertexValenceOfSurface(baseComplex) << "\t" << TotalEdgeValenceOfSurface(baseComplex)
              << "\t" << TotalEdgeValenceOfInterior(baseComplex) << "\n";
    std::cerr << "Singularities " << baseComplex.SingularityI.V.size() << "\t" << baseComplex.SingularityI.E.size() << "\n";
    return 0;
}
//int main(int argc, char* argv[]) {
//    if (argc < 2) {
//        std::cout << "Usage: Test input_file output_file\n";
//        return -1;
//    }
//    ArgumentManager argumentManager(argc, argv);
//    MeshFileReader reader(argv[1]);
//    Mesh& mesh = (Mesh&) reader.GetMesh();
//    for (auto& cell : mesh.C)
//        if (cell.Vids.size() == 8)
//            std::swap(cell.Vids[5], cell.Vids[7]);
//        else if (cell.Vids.size() == 6)
//            std::swap(cell.Vids[4], cell.Vids[5]);
//    MeshFileWriter writer(mesh, argv[2]);
//    writer.WriteFile();
//    return 0;
//}
