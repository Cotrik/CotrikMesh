/*
 * ExtractFaceSegment.cpp
 *
 *  Created on: Dec 5, 2016
 *      Author: cotrik
 */

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include <iostream>
#include "ArgumentManager.h"
int main(int argc, char* argv[])
{
    if (argc < 3)
    {
        std::cerr << "Usage: ExtractFaceSegment input.vtk output.vtk\n";
        return -1;
    }

    MeshFileReader reader(argv[1]);
    reader.GetPointsScalarFields();
    Mesh& mesh = (Mesh&)reader.GetMesh();
    //mesh.BuildAllConnectivities();

    std::vector<Cell> segs;
    for (size_t i = 0; i < mesh.C.size(); i++) {
        const Cell& cell = mesh.C.at(i);
        int scalar0 = mesh.pointScalarFields[0][cell.Vids[0]];
        int scalar1 = mesh.pointScalarFields[0][cell.Vids[1]];
        int scalar2 = mesh.pointScalarFields[0][cell.Vids[2]];
        int scalar3 = mesh.pointScalarFields[0][cell.Vids[3]];
        if (scalar0 == scalar1 && scalar1 == scalar2 && scalar2 == scalar3)
            continue;
        segs.push_back(cell);
    }
    MeshFileWriter writer(mesh.V, segs, argv[2], TETRAHEDRA);
    writer.WriteFile();

    return 0;
}
