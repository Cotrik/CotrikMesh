/*
 * smooth.cpp
 *
 *  Created on: Nov 8, 2016
 *      Author: cotrik
 */



#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkDataSet.h>
#include <vtkMeshQuality.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkTriangleFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkProperty.h>
#include <vtkScalarBarActor.h>
#include <vtkLookupTable.h>
#include <vtkColorTransferFunction.h>
#include <vtkThreshold.h>
#include <vtkExtractEdges.h>

#include <algorithm>
#include "MeshFileReader.h"
#include "MeshFileWriter.h"

void ReadIds(const char* filename, std::vector<size_t>& Ids)
{
    std::ifstream ifs(filename);
    unsigned long id;
    while (ifs >> id)
        Ids.push_back(id);

    std::sort(Ids.begin(), Ids.end());
    std::vector<size_t>::iterator iter = std::unique(Ids.begin(), Ids.end());
    Ids.resize(std::distance(Ids.begin(), iter));
    std::sort(Ids.begin(), Ids.end());
}

void GetVIds(const Mesh& mesh, const std::vector<size_t>& cellIds, std::vector<size_t>& vIds)
{
    for (size_t i = 0; i < cellIds.size(); i++)
    {
        const std::vector<size_t>& vids = mesh.C.at(cellIds.at(i)).Vids;
        std::copy(vids.begin(), vids.end(), std::back_inserter(vIds));
    }

    std::sort(vIds.begin(), vIds.end());
    std::vector<size_t>::iterator iter = std::unique(vIds.begin(), vIds.end());
    vIds.resize(std::distance(vIds.begin(), iter));
    std::sort(vIds.begin(), vIds.end());
}

void Smooth(Mesh& mesh, const std::vector<int>& fixedVertexData)
{
    for (size_t i = 0; i < mesh.V.size(); i++)
        if (fixedVertexData[i] == 0)
        {
            if (!mesh.V.at(i).isBoundary)
            {
                glm::vec3 sum(0.0, 0.0, 0.0);
                const size_t ns = mesh.V[i].N_Vids.size();
                for (size_t j = 0; j < ns; j++)
                {
                    sum.x += mesh.V[mesh.V[i].N_Vids[j]].x;
                    sum.y += mesh.V[mesh.V[i].N_Vids[j]].y;
                    sum.z += mesh.V[mesh.V[i].N_Vids[j]].z;
                }
                Vertex& v = mesh.V[i];

                v.x = sum.x / ns;
                v.y = sum.y / ns;
                v.z = sum.z / ns;
            }
            else
            {
                glm::vec3 sum(0.0, 0.0, 0.0);
                const size_t ns = mesh.V[i].N_Vids.size();
                int sn = 0;
                for (size_t j = 0; j < ns; j++)
                {
                    if (mesh.V[mesh.V[i].N_Vids[j]].isBoundary)
                    {
                        sum.x += mesh.V[mesh.V[i].N_Vids[j]].x;
                        sum.y += mesh.V[mesh.V[i].N_Vids[j]].y;
                        sum.z += mesh.V[mesh.V[i].N_Vids[j]].z;
                        sn++;
                    }
                }
                Vertex& v = mesh.V[i];

                v.x = sum.x / sn;
                v.y = sum.y / sn;
                v.z = sum.z / sn;
            }
        }
}

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " Filename(.vtk) fix_threshold iters cellIds" << std::endl;
        return EXIT_FAILURE;
    }

    // Get the filename from the command line
    std::string inputFilename = argv[1];

    // Get all data from the file
    vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName(inputFilename.c_str());
    reader->Update();

    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGridMesh = reader->GetOutput();
    vtkSmartPointer<vtkMeshQuality> qualityFilter = vtkSmartPointer<vtkMeshQuality>::New();
#if VTK_MAJOR_VERSION <= 5
    qualityFilter->SetInputConnection(mesh->GetProducerPort());
#else
    qualityFilter->SetInputConnection(reader->GetOutputPort());
#endif
    double minValue = 0;
    double maxValue = 0;
    double avgValue = 0;

    qualityFilter->SetHexQualityMeasureToScaledJacobian();
    qualityFilter->Update();

    vtkSmartPointer<vtkDoubleArray> qualityArray = vtkDoubleArray::SafeDownCast(qualityFilter->GetOutput()->GetCellData()->GetArray("Quality"));
    minValue = qualityArray->GetValue(0);
    maxValue = qualityArray->GetValue(0);
    avgValue = 0;
    float fix_threshold = 0.5;
    if (argc == 3)
        fix_threshold = std::stof(argv[2]);
    size_t numOfInvertedElements = 0;
    std::vector<int> fix(qualityArray->GetNumberOfTuples(), 0);
    for (vtkIdType i = 0; i < qualityArray->GetNumberOfTuples(); i++)
    {
        double val = qualityArray->GetValue(i);
        if (minValue > val) minValue = val;
        if (maxValue < val) maxValue = val;
        avgValue += val;

        if (val < 0) numOfInvertedElements++;
        if (val > fix_threshold && val <= 1)
            fix[i] = 1;
        //std::cout << "value " << i << " : " << val << std::endl;
    }
    avgValue /= qualityArray->GetNumberOfTuples();
    std::cout << "min value = " << minValue << std::endl;
    std::cout << "max value = " << maxValue << std::endl;
    std::cout << "avg value = " << avgValue << std::endl;
    std::cout << "#InvertedElements = " << numOfInvertedElements << std::endl;

    MeshFileReader myreader(argv[1]);
    Mesh& mesh = (Mesh&)myreader.GetMesh();
    mesh.BuildAllConnectivities();
    mesh.ExtractBoundary();
    std::vector<int> fixedVertexData(mesh.V.size(), 1);
    for (size_t i = 0; i < fix.size(); i++)
        if (fix[i] == 0)
        for (size_t j = 0; j < mesh.C[i].Vids.size(); j++)
            fixedVertexData[mesh.C[i].Vids[j]] = 0;

    int iter = 1;
    if (argc == 4)
        iter = std::stoi(argv[3]);
    if (argc == 5)
    {
        std::vector<size_t> cellIds;
        std::vector<size_t> vIds;
        ReadIds(argv[4], cellIds);
        GetVIds(mesh, cellIds, vIds);
        fixedVertexData.resize(fixedVertexData.size(), 1);
        for (size_t i = 0; i < vIds.size(); i++)
            fixedVertexData.at(vIds.at(i)) = 0;
    }
    while (iter--)
        Smooth(mesh, fixedVertexData);

    std::string fixFileName = inputFilename.substr(0, inputFilename.size() - 4) + "_smoothed.vtk";
    MeshFileWriter fixFileWriter(mesh.V, mesh.C, fixFileName.c_str(), HEXAHEDRA);
    fixFileWriter.WriteFile();
    //fixFileWriter.WritePointData(fixedVertexData, "fixed");

//    if (argc == 2)
        return 0;
}



