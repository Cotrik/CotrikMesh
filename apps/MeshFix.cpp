/*
 * meshfix.cpp
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

#include "MeshFileReader.h"
#include "MeshFileWriter.h"

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " Filename(.vtk)" << std::endl;
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
    //qualityFilter->SetTriangleQualityMeasureToArea();
    //qualityFilter->SetHexQualityMeasureToEdgeRatio()
    //qualityFilter->SetHexQualityMeasureToMedAspectFrobenius()
    //qualityFilter->SetHexQualityMeasureToMaxAspectFrobenius()
    //qualityFilter->SetHexQualityMeasureToMaxEdgeRatios()
    //qualityFilter->SetHexQualityMeasureToSkew()
    //qualityFilter->SetHexQualityMeasureToTaper()
    //qualityFilter->SetHexQualityMeasureToVolume()
    //qualityFilter->SetHexQualityMeasureToStretch()
    //qualityFilter->SetHexQualityMeasureToDiagonal()
    //qualityFilter->SetHexQualityMeasureToDimension()
    //qualityFilter->SetHexQualityMeasureToOddy()
    //qualityFilter->SetHexQualityMeasureToCondition()
    //qualityFilter->SetHexQualityMeasureToJacobian()
    qualityFilter->SetHexQualityMeasureToScaledJacobian();
    //qualityFilter->SetHexQualityMeasureToShear()
    //qualityFilter->SetHexQualityMeasureToShape()
    //qualityFilter->SetHexQualityMeasureToRelativeSizeSquared()
    //qualityFilter->SetHexQualityMeasureToShapeAndSize()
    //qualityFilter->SetHexQualityMeasureToShearAndSize()
    //qualityFilter->SetHexQualityMeasureToDistortion()
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
        if (val > fix_threshold)
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
    std::string fixFileName = inputFilename.substr(0, inputFilename.size() - 4) + "_fixed.vtk";
    MeshFileWriter innerCellsWriter(mesh.V, mesh.C, fixFileName.c_str(), HEXAHEDRA);
    innerCellsWriter.WriteFile();
    std::vector<int> fixedVertexData(mesh.V.size(), 0);
    for (size_t i = 0; i < fix.size(); i++)
        if (fix[i] == 1)
        for (size_t j = 0; j < mesh.C[i].Vids.size(); j++)
            fixedVertexData[mesh.C[i].Vids[j]] = 1;
    innerCellsWriter.WritePointData(fixedVertexData, "fixed");

//    if (argc == 2)
        return 0;


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Hightlight inverted elements
    vtkDataSet* qualityMesh = qualityFilter->GetOutput();
    vtkSmartPointer<vtkThreshold> selectCells = vtkSmartPointer<vtkThreshold>::New();
    selectCells->ThresholdByLower(0.00);
    selectCells->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, vtkDataSetAttributes::SCALARS);
#if VTK_MAJOR_VERSION <= 5
    selectCells->SetInput(qualityMesh);
#else
    selectCells->SetInputData(qualityMesh);
#endif
    selectCells->Update();
    vtkUnstructuredGrid* ug = selectCells->GetOutput();
    // Create a mapper and actor
    vtkSmartPointer<vtkDataSetMapper> inverted_mapper = vtkSmartPointer<vtkDataSetMapper>::New();
#if VTK_MAJOR_VERSION <= 5
    inverted_mapper->SetInput(ug);
#else
    inverted_mapper->SetInputData(ug);
#endif
    vtkSmartPointer<vtkActor> inverted_actor = vtkSmartPointer<vtkActor>::New();
    inverted_actor->SetMapper(inverted_mapper);
    inverted_actor->GetProperty()->SetColor(1.0, 0.0, 1.0);
    inverted_actor->GetProperty()->EdgeVisibilityOn();
    inverted_actor->GetProperty()->SetRepresentationToSurface();
    inverted_actor->GetProperty()->SetLineWidth(5);
    inverted_actor->GetProperty()->SetEdgeColor(1.0, 0, 1.0);
    //inverted_actor->GetProperty()->SetOpacity(0.4);

    // ------------------------------------
    // Display Edges
    vtkSmartPointer<vtkExtractEdges> extractEdges = vtkSmartPointer<vtkExtractEdges>::New();
    extractEdges->SetInputData(unstructuredGridMesh);
    extractEdges->Update();
    // Create a mapper and actor
    vtkSmartPointer<vtkPolyDataMapper> edges_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    edges_mapper->SetInputConnection(extractEdges->GetOutputPort());
    vtkSmartPointer<vtkActor> edges_actor = vtkSmartPointer<vtkActor>::New();
    edges_actor->SetMapper(edges_mapper);

    // ------------------------------------
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    unstructuredGrid->ShallowCopy(qualityFilter->GetOutput());

    // Visualize
    vtkSmartPointer<vtkDataSetMapper> mapper = vtkSmartPointer<vtkDataSetMapper>::New();
#if VTK_MAJOR_VERSION <= 5
    mapper->SetInputConnection(unstructuredGrid->GetProducerPort());
#else
    mapper->SetInputData(unstructuredGrid);
#endif
    mapper->SetScalarRange(unstructuredGrid->GetScalarRange());

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(0.9, 0.9, 0.9);
    actor->GetProperty()->EdgeVisibilityOn();
    actor->GetProperty()->SetRepresentationToSurface();
    actor->GetProperty()->SetEdgeColor(0, 0, 0);
    actor->GetProperty()->SetOpacity(0.4);

    vtkSmartPointer<vtkScalarBarActor> scalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
    scalarBar->SetLookupTable(mapper->GetLookupTable());
    scalarBar->SetTitle("Scaled Jacobian");
    scalarBar->SetNumberOfLabels(4);

    // Create a lookup table to share between the mapper and the scalarbar
//    vtkSmartPointer<vtkLookupTable> hueLut = vtkSmartPointer<vtkLookupTable>::New();
//    hueLut->SetTableRange(0, 1);
//    hueLut->SetHueRange(0, 1);
//    hueLut->SetSaturationRange(1, 1);
//    hueLut->SetValueRange(1, 1);
//    hueLut->Build();
//
//    mapper->SetLookupTable(hueLut);
//    scalarBar->SetLookupTable(hueLut);


    vtkSmartPointer<vtkColorTransferFunction> ctf = vtkSmartPointer<vtkColorTransferFunction>::New();
    if (minValue < 0 && maxValue > 0)
    {
        ctf->AddRGBPoint(minValue, 0, 0, 1);
        ctf->AddRGBPoint(0, 1, 1, 1);
        ctf->AddRGBPoint(maxValue, 1, 0, 0);
    }
    else if (minValue > 0 && maxValue > 0)
    {
        ctf->AddRGBPoint(minValue, 1, 1, 1);
        ctf->AddRGBPoint(maxValue, 1, 0, 0);
    }
    else if (minValue < 0 && maxValue < 0)
    {
        ctf->AddRGBPoint(minValue, 0, 0, 1);
        ctf->AddRGBPoint(maxValue, 1, 1, 1);
    }
    ctf->SetColorSpaceToDiverging();
    ctf->Build();
    mapper->SetLookupTable(ctf);
    scalarBar->SetLookupTable(ctf);


    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);
    vtkSmartPointer<vtkInteractorStyleTrackballCamera>  InteractorStyleTrackballCamera = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
    renderWindowInteractor->SetInteractorStyle(InteractorStyleTrackballCamera);

    renderer->AddActor(actor);
    renderer->AddActor2D(scalarBar);
    renderer->AddActor(inverted_actor);
    renderer->AddActor(edges_actor);

    renderer->SetBackground(1, 1, 1); // Background color green

    renderWindow->Render();
    renderWindowInteractor->Start();

    return EXIT_SUCCESS;
}



