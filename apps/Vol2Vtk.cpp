/*
 * Vol2Vtk.cpp
 *
 *  Created on: Nov 11, 2018
 *      Author: cotrik
 */

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include <iostream>
#include <sstream>
using namespace std;

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cout << "Usage: vol2vtk vol_file polydata.vtk" << std::endl;
        return -1;
    }

    ifstream file(argv[1]);

    std::vector <Vertex> V;
    std::vector <Cell> C;

    std::string str;
    // read cells
    unsigned long cellsNum = 0;
    while (getline(file, str)) {
        if (str.find("surfaceelementsgi") != str.npos) {
            getline(file, str);
            stringstream ss(str.c_str());
            ss >> cellsNum;
            std::cout << "#faces : " << cellsNum << std::endl;
            for (unsigned long i = 0; i < cellsNum; i++) {
                getline(file, str);
                stringstream ss(str.c_str());
                stringstream cell_strstream(str.c_str());
                unsigned long numOfVids;
                cell_strstream >> numOfVids >> numOfVids >> numOfVids >> numOfVids >> numOfVids;
                Cell cell(numOfVids);
                for (auto & vid : cell.Vids) {
                    cell_strstream >> vid;
                    --vid;
                }
                cell.cellType = numOfVids == 3 ? VTK_TRIANGLE : VTK_QUAD;
                C.push_back(cell);
            }
            C.resize(C.size());
            break;
        }
    }

    // read vertices
    unsigned long verticesNum = 0;
    while (getline(file, str)) {
        if (str.find("points") != str.npos) {
            getline(file, str);
            stringstream strstream(str.c_str());
            strstream >> verticesNum;
            std::cout << "points : " << verticesNum << std::endl;
            for (unsigned long i = 0; i < verticesNum; i++) {
                getline(file, str);
                stringstream vertex_strstream(str.c_str());
                Vertex v;
                vertex_strstream >> v.x >> v.y >> v.z;
                V.push_back(v);
            }
            V.resize(V.size());
            break;
        }
    }

    MeshFileWriter writer(V, C, argv[2], POLYGON);
    writer.WriteFile();
    file.close();
}



