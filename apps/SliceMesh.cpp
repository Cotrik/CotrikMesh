/*
 * SliceMesh.cpp
 *
 *  Created on: Dec 3, 2017
 *      Author: cotrik
 */

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "ArgumentManager.h"
#include <iostream>
#include <sstream>
#include <iomanip>

#include <vtkSmartPointer.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPlane.h>
#include <vtkClipDataSet.h>
#include <vtkCutter.h>

void getMeshSizeInfo(const Mesh& mesh, glm::vec3& mn, glm::vec3& mx, glm::vec3& s) {
    for (const auto& v : mesh.V) {
        if (v.x < mn.x) mn.x = v.x;
        if (v.y < mn.y) mn.y = v.y;
        if (v.z < mn.z) mn.z = v.z;

        if (v.x > mx.x) mx.x = v.x;
        if (v.y > mx.y) mx.y = v.y;
        if (v.z > mx.z) mx.z = v.z;
    }
    s = mx - mn;
    std::cout << "######## mesh size info ########\n";
    std::cout << "x : " << mn.x << " ~ " << mx.x << " size : " << s.x << "\n";
    std::cout << "y : " << mn.y << " ~ " << mx.y << " size : " << s.y << "\n";
    std::cout << "z : " << mn.z << " ~ " << mx.z << " size : " << s.z << "\n";
}

void getProj(const glm::vec3& mn, const glm::vec3& mx, const glm::vec3& s,
        const glm::vec3& proj_mn, const glm::vec3& proj_mx, const glm::vec3& proj_s, glm::vec3& normal) {
}
int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cout << "Usage: SliceMesh <file> <out.vtk> orig=\"0.0 0.0 0.0\" normal=\"1.0 0.0 0.0\" slice=1\n";
        return -1;
    }
    ArgumentManager am(argc, argv);
    glm::vec3 orig(0.0, 0.0, 0.0);
    glm::vec3 normal(1.0, 0.0, 0.0);
    int slice = 1;
    {
        string strorig = am.get("orig");
        if (!strorig.empty()) {
            std::stringstream ss(strorig);
            ss >> orig.x >> orig.y >> orig.z;
        }
        string strnormal = am.get("normal");
        if (!strnormal.empty()) {
            std::stringstream ss(strnormal);
            ss >> normal.x >> normal.y >> normal.z;
        }
        string strslice = am.get("slice");
        if (!strslice.empty()) slice = std::stoi(strslice);

        std::cout << "input = " << argv[1] << "\n";
        std::cout << "output = " << argv[2] << "\n";
        std::cout << "orig = (" << orig.x << ", " << orig.y << ", " << orig.z << ")\n";
        std::cout << "normal = (" << normal.x << ", " << normal.y << ", " << normal.z << ")\n";
        std::cout << "#slice = " << slice << "\n";
    }

    if (slice == 1) {
        vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
        reader->SetFileName(argv[1]);
        reader->Update();

        vtkSmartPointer<vtkUnstructuredGrid> model = reader->GetOutput();
        vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();

        plane->SetOrigin(orig.x, orig.y, orig.z);
        plane->SetNormal(normal.x, normal.y, normal.z);

        vtkSmartPointer<vtkCutter> clipDataSet = vtkSmartPointer<vtkCutter>::New();
        clipDataSet->SetCutFunction(plane);
        clipDataSet->SetInputConnection(reader->GetOutputPort());
        //clipDataSet->InsideOutOn();
        clipDataSet->GenerateTrianglesOn();
        clipDataSet->Update();
        //clipDataSet->GetOutput()->Print(std::cout);

        vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
        writer->SetFileName(argv[2]);
        writer->SetInputData(clipDataSet->GetOutput());
        writer->Update();
    } else {
        glm::vec3 mn(INT_MAX, INT_MAX, INT_MAX);
        glm::vec3 mx(INT_MIN, INT_MIN, INT_MIN);
        glm::vec3 s;
        {
            MeshFileReader reader(argv[1]);
            Mesh& mesh = (Mesh&) reader.GetMesh();
            getMeshSizeInfo(mesh, mn, mx, s);
        }
        glm::vec3 proj_mn(INT_MAX, INT_MAX, INT_MAX);
        glm::vec3 proj_mx(INT_MIN, INT_MIN, INT_MIN);
        glm::vec3 proj_s;
        //getProj_s(proj_mn, proj_mx, proj_s);

        vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
        reader->SetFileName(argv[1]);
        reader->Update();

        vtkSmartPointer<vtkUnstructuredGrid> model = reader->GetOutput();
        vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();

        double d = s.x / slice;
        for (int i = 0; i < slice; ++i) {
            plane->SetOrigin(mn.x + (i + 1) * d, 0, 0);
            plane->SetNormal(normal.x, normal.y, normal.z);

            vtkSmartPointer<vtkCutter> clipDataSet = vtkSmartPointer<vtkCutter>::New();
            clipDataSet->SetCutFunction(plane);
            clipDataSet->SetInputConnection(reader->GetOutputPort());
            //clipDataSet->InsideOutOn();
            clipDataSet->GenerateTrianglesOn();
            clipDataSet->Update();
            //clipDataSet->GetOutput()->Print(std::cout);

            vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
            std::string filename = "Slice" + std::to_string(i) + ".vtk";
            writer->SetFileName(filename.c_str());
            writer->SetInputData(clipDataSet->GetOutput());
            writer->Update();
        }
    }
    return 0;
}

