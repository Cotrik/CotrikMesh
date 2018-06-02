/*
 * ClipMesh.cpp
 *
 *  Created on: Dec 3, 2017
 *      Author: cotrik
 */

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "ArgumentManager.h"
#include <iostream>
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

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cout << "Usage: ClipMesh <file> <out.vtk> orig=\"0.0, 0.0, 0.0\" normal=\"1.0, 0.0, 0.0\"\n";
        return -1;
    }
    ArgumentManager argumentManager(argc, argv);

    vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName(argv[1]);
    reader->Update();

    vtkSmartPointer<vtkUnstructuredGrid> model = reader->GetOutput();
    vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
    plane->SetOrigin(0.419349824500387, 0.189770378710, 0.195257135783322); plane->SetNormal(1,0,0);

    // Clip data
    vtkSmartPointer<vtkClipDataSet> clipDataSet = vtkSmartPointer<vtkClipDataSet>::New();
    clipDataSet->SetClipFunction(plane);
    clipDataSet->SetInputConnection(reader->GetOutputPort());
    clipDataSet->InsideOutOn();
    clipDataSet->GenerateClippedOutputOn();
    clipDataSet->Update();
//    clipDataSet->GetOutput()->Print(std::cout);

    vtkSmartPointer<vtkUnstructuredGridWriter> writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    writer->SetFileName(argv[2]);
    writer->SetInputData(clipDataSet->GetOutput());
    writer->Update();

    return 0;
}

