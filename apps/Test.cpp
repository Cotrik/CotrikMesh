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

int TotalValence(const BaseComplex& bc) {
    int sum = 0;
    for (auto& sv : bc.SingularityI.V) {
        auto valence = bc.mesh.V.at(sv.id_mesh).N_Cids.size();
        if (sv.isBoundary) sum += 4 - valence;
        else sum += 8 - valence;
    }
    for (auto& se : bc.SingularityI.E) {
        auto valence = bc.mesh.E.at(se.es_link.front()).N_Cids.size();
        if (se.isBoundary) sum -= 2 - valence;
        else sum -= 4 - valence;
    }
    return sum;
}

int TotalValence1(const BaseComplex& bc) {
    int sum = 0;
    for (auto& sv : bc.SingularityI.V) {
        auto valence = bc.mesh.V.at(sv.id_mesh).N_Cids.size();
        if (sv.isBoundary) sum += 4 - valence;
        else sum += 8 - valence;
    }
    for (auto& se : bc.SingularityI.E) {
        for (auto x : se.es_link) {
            auto valence = bc.mesh.E.at(se.es_link.front()).N_Cids.size();
            if (se.isBoundary) sum -= 2 - valence;
            else sum -= 4 - valence;
        }
        for (auto x = 1; x < se.vs_link.size() - 1; ++x) {
            auto valence = bc.mesh.V.at(se.vs_link.at(x)).N_Cids.size();
            if (se.isBoundary) sum -= 4 - valence;
            else sum -= 8 - valence;
        }
    }
    return sum;
}

double TotalIndex(const BaseComplex& bc) {
    double sum = 0;
    for (auto& sv : bc.SingularityI.V) {
        auto valence = bc.mesh.V.at(sv.id_mesh).N_Cids.size();
        if (sv.isBoundary) {
            sum += 0.5 - valence * 0.125;
            //std::cout << "1 " << valence << " " << 0.5 - valence * 0.125 << "\n";
        }
        else {
            sum += 1.0 - valence * 0.125;
            //std::cout << "0 " << valence << " " << 0.5 - valence * 0.125 << "\n";
        }
    }
    std::cout << "\n";
    for (auto& se : bc.SingularityI.E) {
        if (se.vs_link.empty() || se.vs_link.front() == se.vs_link.back()) continue;
        auto valence = bc.mesh.E.at(se.es_link.front()).N_Cids.size();
        if (se.isBoundary) {
            sum -= 0.5 - valence * 0.25;
            //std::cout << "1 " << valence << " " << -(0.5 - valence * 0.25) << "\n";
        }
        else {
            sum -= 1.0 - valence * 0.25;
            //std::cout << "0 " << valence << " " << -(1.0 - valence * 0.25) << "\n";
        }
    }
    return sum;
}

double TotalIndex1(const BaseComplex& bc) {
    double sum = 0;
    for (auto& sv : bc.SingularityI.V) {
        auto valence = bc.mesh.V.at(sv.id_mesh).N_Cids.size();
        if (sv.isBoundary) sum += 0.5 - valence * 0.125;
        else sum += 1.0 - valence * 0.125;
    }
    for (auto& se : bc.SingularityI.E) {
        for (auto x : se.es_link) {
            auto valence = bc.mesh.E.at(se.es_link.front()).N_Cids.size();
            if (se.isBoundary) sum -= 0.5 - valence * 0.25;
            else sum -= 1.0 - valence * 0.25;
        }
        for (auto x = 1; x < se.vs_link.size() - 1; ++x) {
            auto valence = bc.mesh.V.at(se.vs_link.at(x)).N_Cids.size();
            if (se.isBoundary) sum -= 0.5 - valence * 0.125;
            else sum -= 1.0 - valence * 0.125;
        }
    }
    return sum;
}

int Genus(const Mesh& mesh) {
    return  1 - (int(mesh.V.size()) - mesh.E.size() + mesh.F.size() - mesh.C.size());
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
    mesh.RemoveUselessVertices();
    mesh.BuildAllConnectivities();
    mesh.ExtractBoundary();
    mesh.ExtractSingularities();
//    std::cout << "-----------------\n";
    std::cout << "#V = " << mesh.V.size() << "\t";
    std::cout << "#E = " << mesh.E.size() << "\t";
    std::cout << "#F = " << mesh.F.size() << "\t";
    std::cout << "#C = " << mesh.C.size() << "\n";
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
//    std::cout << "-----------------\n";
    std::cout << "#BV=" << baseComplex.Vids.size() << "\t";
    std::cout << "#BE=" << baseComplex.Eids.size() << "\t";
    std::cout << "#BF=" << baseComplex.Fids.size() << "\t";
    std::cout << "#BC=" << baseComplex.Cids.size() << "\n";
//    std::cout << "-----------------\n";
    std::cout << "#VB=" << baseComplex.componentV.size() << "\t";
    std::cout << "#EB=" << baseComplex.componentE.size() << "\t";
    std::cout << "#FB=" << baseComplex.componentF.size() << "\t";
    std::cout << "#CB=" << baseComplex.componentC.size() << "\n";

    std::cout << "#VS=" << baseComplex.SingularityI.V.size() << "\t";
    std::cout << "#ES=" << baseComplex.SingularityI.E.size() << "\n";

    std::cerr << "Genus = " << Genus(mesh) << "\n";
    std::cerr << "TotalVertexValenceOfSurface = " << TotalVertexValenceOfSurface(baseComplex) << "\n";
    std::cerr << "TotalIndex = " << TotalIndex(baseComplex) << "\n";
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
