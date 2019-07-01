/*
 * ExtractSurface.cpp
 *
 *  Created on: Oct 16, 2018
 *      Author: cotrik
 */

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "ArgumentManager.h"

#include <iostream>

int main(int argc, char* argv[]){
    if (argc < 2) {
        std::cout << "Usage: ExtractSurface input output compressed" << std::endl;
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

    std::vector<Cell> C;
    int i = 0;
    Cell cell;
    for (auto& f : mesh.F) {
        if (f.isBoundary) {
            cell.Vids = f.Vids;
            cell.id = i++;
            C.push_back(cell);
        }
    }
    auto cellType = TRIANGLE;
    if (mesh.m_cellType == HEXAHEDRA) cellType = QUAD;
    else if (mesh.m_cellType == TETRAHEDRA) cellType = TRIANGLE;
    else if (mesh.m_cellType == QUAD) cellType = QUAD;
    else if (mesh.m_cellType == TRIANGLE) cellType = TRIANGLE;

    Mesh surface(mesh.V, C, cellType);
    if (argc >= 4) surface.RemoveUselessVertices();
    MeshFileWriter writer(surface, argv[2]);
    writer.WriteFile();

    if (mesh.m_cellType == QUAD || mesh.m_cellType == TRIANGLE) {
        std::vector<size_t> edgeids;
        //std::set<size_t> boundary_vids;
        for (auto& e : mesh.E)
            if (e.N_Fids.size() == 1) {
                edgeids.push_back(e.id);
                //boundary_vids.insert(e.Vids[0]);
                //boundary_vids.insert(e.Vids[1]);
            }
        {
            MeshFileWriter writer(mesh, "boundary.vtk");
            writer.WriteEdgesVtk(edgeids);
        }
    }
    return 0;
}
