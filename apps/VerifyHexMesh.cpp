/*
 * VerifyHexMesh.cpp
 *
 *  Created on: Nov 28, 2017
 *      Author: cotrik
 */

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include <iostream>
#include <unordered_set>
#include "ArgumentManager.h"
int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cerr << "Usage: VerifyHexMesh input.vtk\n";
        return -1;
    }

    MeshFileReader reader(argv[1]);
    reader.GetPointsScalarFields();
    Mesh& mesh = (Mesh&)reader.GetMesh();
    mesh.BuildAllConnectivities();
    mesh.ExtractBoundary();
    mesh.ExtractSingularities();
    mesh.BuildParallelE();
    mesh.BuildConsecutiveE();
    mesh.BuildOrthogonalE();

    bool flag = true;
    // Veriy foreach vi 3h=2f
    for (auto& v : mesh.V) {
        if (v.isBoundary) continue;
        auto f = v.N_Fids.size();
        auto h = v.N_Cids.size();
        if (3 * h != 2 * f) {
            cout << "1.\t3h != 2f for vertex" << v.id << "\n";
            flag = false;
            break;
        }
    }
    if (flag) cout << "1.\t3h = 2f\n";

    flag = true;
    // Veriy foreach vi 3e - f = 6
    for (auto& v : mesh.V) {
        if (v.isBoundary) continue;
        auto f = v.N_Fids.size();
        auto e = v.N_Eids.size();
        if (3 * e - f != 6) {
            cout << "2.\t3e - f != 6 for vertex" << v.id << "\n";
            flag = false;
            break;
        }
    }
    if (flag) cout << "2.\t3e - f = 6\n";
    auto E = mesh.E.size();
    auto F = mesh.F.size();
    auto H = mesh.C.size();
    for (auto& e : mesh.E) if (e.isBoundary) --E;
    for (auto& f : mesh.F) if (f.isBoundary) --F;

    if (3 * H != F)  cout << "3.\t3H != F\n";
    else cout << "3.\t 3H = F\n";

    if (3 * E - 2 * F != 3)  cout << "4.\t3E - 2F != 3\n";
    else cout << "4.\t3E - 2F = 3\n";

    cout << "5.\tE = " << E << ", F = " << F << ", H = " << H << "\n";
    cout << "6.\t3E = " << 3 * E << ", 2F = " << 2 * F << ", 3H = " << 3 * H << "\n";
    return 0;
}


