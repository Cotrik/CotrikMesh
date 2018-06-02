/*
 * Test.cpp
 *
 *  Created on: Jan 8, 2018
 *      Author: cotrik
 */

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "SmoothAlgorithm.h"
#include <iostream>
#include "ArgumentManager.h"

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: Test input_file output_file\n";
        return -1;
    }
    ArgumentManager argumentManager(argc, argv);
    MeshFileReader reader(argv[1]);
    Mesh& mesh = (Mesh&) reader.GetMesh();
    for (auto& cell : mesh.C)
        if (cell.Vids.size() == 8)
            std::swap(cell.Vids[5], cell.Vids[7]);
        else if (cell.Vids.size() == 6)
            std::swap(cell.Vids[4], cell.Vids[5]);
    MeshFileWriter writer(mesh, argv[2]);
    writer.WriteFile();
    return 0;
}
