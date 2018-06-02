/*
 * Rotate.cpp
 *
 *  Created on: Jan 10, 2018
 *      Author: cotrik
 */

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include <iostream>
#include <fstream>
#include <string.h>

void SetRotMatrix(double rx[3][3], double ry[3][3], double rz[3][3], double theta_x, double theta_y, double theta_z)
{
  rx[0][0] = 1;
  rx[1][1] = cos(theta_x);
  rx[1][2] = sin(theta_x);
  rx[2][2] = cos(theta_x);
  rx[2][1] = -sin(theta_x);

  ry[1][1] = 1;
  ry[0][0] = cos(theta_y);
  ry[0][2] = -sin(theta_y);
  ry[2][2] = cos(theta_y);
  ry[2][0] = sin(theta_y);

  rz[2][2] = 1;
  rz[0][0] = cos(theta_z);
  rz[0][1] = sin(theta_z);
  rz[1][1] = cos(theta_z);
  rz[1][0] = -sin(theta_z);
}
const double PI = 3.1415926536;

int main(int argc, char *argv[]) {
    int num = 10;
    if (argc == 3) num = atoi(argv[2]);
    else if (argc != 2) {
        std::cerr << "Usage: Rotate input numOfRotationsPerAxis" << std::endl;
        return EXIT_FAILURE;
    }
    MeshFileReader reader(argv[1]);
    Mesh& mesh = (Mesh&)reader.GetMesh();
    //mesh.RemoveUselessVertices();
    //mesh.BuildAllConnectivities();
    //mesh.ExtractBoundary();
    const auto& V = mesh.V;
    const auto& C = mesh.C;
    double rx[3][3] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    double ry[3][3] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    double rz[3][3] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    for (int i = 0; i < num; i++) {
        float theta_x = 2.0 * PI * i / num;
        for (int j = 0; j < num; j++) {
            float theta_y = 2.0 * PI * j / num;
            for (int k = 0; k < num; k++) {
                float theta_z = 2.0 * PI * k / num;
                SetRotMatrix(rx, ry, rz, theta_x, theta_y, theta_z);
                std::vector<Vertex> Vx(V);
                for (size_t l = 0; l < Vx.size(); l++) {
                    Vx[l].x = rx[0][0] * V[l].x + rx[0][1] * V[l].y + rx[0][2] * V[l].z;
                    Vx[l].y = rx[1][0] * V[l].x + rx[1][1] * V[l].y + rx[1][2] * V[l].z;
                    Vx[l].z = rx[2][0] * V[l].x + rx[2][1] * V[l].y + rx[2][2] * V[l].z;

                    Vx[l].x = ry[0][0] * V[l].x + ry[0][1] * V[l].y + ry[0][2] * V[l].z;
                    Vx[l].y = ry[1][0] * V[l].x + ry[1][1] * V[l].y + ry[1][2] * V[l].z;
                    Vx[l].z = ry[2][0] * V[l].x + ry[2][1] * V[l].y + ry[2][2] * V[l].z;

                    Vx[l].x = rz[0][0] * V[l].x + rz[0][1] * V[l].y + rz[0][2] * V[l].z;
                    Vx[l].y = rz[1][0] * V[l].x + rz[1][1] * V[l].y + rz[1][2] * V[l].z;
                    Vx[l].z = rz[2][0] * V[l].x + rz[2][1] * V[l].y + rz[2][2] * V[l].z;
                }
                std::string strx = std::to_string(360 * i / num);
                std::string stry = std::to_string(360 * j / num);
                std::string strz = std::to_string(360 * k / num);
                const std::string str = std::string(argv[1]).substr(0, strlen(argv[1]) - 4) + std::string("_") + std::string(strx.c_str()) + std::string("_")
                        + std::string(stry.c_str()) + std::string("_") + std::string(strz.c_str()) + std::string(".vtk");
                MeshFileWriter meshWriter(Vx, C, str.c_str(), TETRAHEDRA);
                meshWriter.WriteFile();

            }
        }
    }

    return 0;
}
