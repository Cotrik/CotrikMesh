#include <iostream>
#include "MeshFileReader.h"
#include "MeshFileWriter.h"

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
#include <vtkTextProperty.h>
#include <vtkLight.h>
#include <vtkLightActor.h>

void GetMetrics(vtkMeshQuality* qualityFilter, double& minValue, double& maxValue, double& avgValue) {
    vtkSmartPointer<vtkDoubleArray> qualityArray = vtkDoubleArray::SafeDownCast(qualityFilter->GetOutput()->GetCellData()->GetArray("Quality"));
    minValue = qualityArray->GetValue(0);
    maxValue = qualityArray->GetValue(0);
    avgValue = 0;
    for (vtkIdType i = 0; i < qualityArray->GetNumberOfTuples(); i++) {
        double val = qualityArray->GetValue(i);
        if (minValue > val) minValue = val;
        if (maxValue < val) maxValue = val;
        avgValue += val;
        //std::cout << "value " << i << " : " << val << std::endl;
    }
    avgValue /= qualityArray->GetNumberOfTuples();
}

void OutputMetricFile(const char* vtkfilename, const char* filename = "metric.txt") {
    vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName(vtkfilename);
    reader->Update();
    std::cout << "#Vertices = " << reader->GetOutput()->GetNumberOfPoints() << std::endl;
    std::cout << "#Cells = " << reader->GetOutput()->GetNumberOfCells() << std::endl;
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGridMesh = reader->GetOutput();
    vtkSmartPointer<vtkMeshQuality> qualityFilter = vtkSmartPointer<vtkMeshQuality>::New();
#if VTK_MAJOR_VERSION <= 5
    qualityFilter->SetInputConnection(mesh->GetProducerPort());
#else
    ///qualityFilter->SetInputData(mesh);
    qualityFilter->SetInputConnection(reader->GetOutputPort());
#endif
    double minValue = 0;
    double maxValue = 0;
    double avgValue = 0;
    std::vector<double> metrics;
    //--------------------------
    qualityFilter->SetHexQualityMeasureToDiagonal();
    qualityFilter->Update();
    GetMetrics(qualityFilter, minValue, maxValue, avgValue);
    metrics.push_back(minValue);
    metrics.push_back(avgValue);
    //--------------------------
    qualityFilter->SetHexQualityMeasureToDimension();
    qualityFilter->Update();
    GetMetrics(qualityFilter, minValue, maxValue, avgValue);
    metrics.push_back(maxValue);
    metrics.push_back(avgValue);
    //--------------------------
    qualityFilter->SetHexQualityMeasureToDistortion();
    qualityFilter->Update();
    GetMetrics(qualityFilter, minValue, maxValue, avgValue);
    metrics.push_back(minValue);
    metrics.push_back(avgValue);
    //--------------------------
    qualityFilter->SetHexQualityMeasureToEdgeRatio();
    qualityFilter->Update();
    GetMetrics(qualityFilter, minValue, maxValue, avgValue);
    metrics.push_back(maxValue);
    metrics.push_back(avgValue);
    //--------------------------
    qualityFilter->SetHexQualityMeasureToJacobian();
    qualityFilter->Update();
    GetMetrics(qualityFilter, minValue, maxValue, avgValue);
    metrics.push_back(minValue);
    metrics.push_back(avgValue);
    //--------------------------
    qualityFilter->SetHexQualityMeasureToMaxEdgeRatios();
    qualityFilter->Update();
    GetMetrics(qualityFilter, minValue, maxValue, avgValue);
    metrics.push_back(maxValue);
    metrics.push_back(avgValue);
    //--------------------------
    qualityFilter->SetHexQualityMeasureToMaxAspectFrobenius();
    qualityFilter->Update();
    GetMetrics(qualityFilter, minValue, maxValue, avgValue);
    metrics.push_back(maxValue);
    metrics.push_back(avgValue);
    //--------------------------
    qualityFilter->SetHexQualityMeasureToMedAspectFrobenius();
    qualityFilter->Update();
    GetMetrics(qualityFilter, minValue, maxValue, avgValue);
    metrics.push_back(maxValue);
    metrics.push_back(avgValue);
    //--------------------------
    qualityFilter->SetHexQualityMeasureToOddy();
    qualityFilter->Update();
    GetMetrics(qualityFilter, minValue, maxValue, avgValue);
    metrics.push_back(maxValue);
    metrics.push_back(avgValue);
    //--------------------------
    qualityFilter->SetHexQualityMeasureToRelativeSizeSquared();
    qualityFilter->Update();
    GetMetrics(qualityFilter, minValue, maxValue, avgValue);
    metrics.push_back(minValue);
    metrics.push_back(avgValue);
    //--------------------------
    qualityFilter->SetHexQualityMeasureToScaledJacobian();
    qualityFilter->Update();
    GetMetrics(qualityFilter, minValue, maxValue, avgValue);
    metrics.push_back(minValue);
    metrics.push_back(avgValue);
    //--------------------------
    qualityFilter->SetHexQualityMeasureToShape();
    qualityFilter->Update();
    GetMetrics(qualityFilter, minValue, maxValue, avgValue);
    metrics.push_back(minValue);
    metrics.push_back(avgValue);
    //--------------------------
    qualityFilter->SetHexQualityMeasureToShapeAndSize();
    qualityFilter->Update();
    GetMetrics(qualityFilter, minValue, maxValue, avgValue);
    metrics.push_back(minValue);
    metrics.push_back(avgValue);
    //--------------------------
    qualityFilter->SetHexQualityMeasureToShear();
    qualityFilter->Update();
    GetMetrics(qualityFilter, minValue, maxValue, avgValue);
    metrics.push_back(minValue);
    metrics.push_back(avgValue);
    //--------------------------
    qualityFilter->SetHexQualityMeasureToShearAndSize();
    qualityFilter->Update();
    GetMetrics(qualityFilter, minValue, maxValue, avgValue);
    metrics.push_back(minValue);
    metrics.push_back(avgValue);
    //--------------------------
    qualityFilter->SetHexQualityMeasureToSkew();
    qualityFilter->Update();
    GetMetrics(qualityFilter, minValue, maxValue, avgValue);
    metrics.push_back(maxValue);
    metrics.push_back(avgValue);
    //--------------------------
    qualityFilter->SetHexQualityMeasureToStretch();
    qualityFilter->Update();
    GetMetrics(qualityFilter, minValue, maxValue, avgValue);
    metrics.push_back(minValue);
    metrics.push_back(avgValue);
    //--------------------------
    qualityFilter->SetHexQualityMeasureToTaper();
    qualityFilter->Update();
    GetMetrics(qualityFilter, minValue, maxValue, avgValue);
    metrics.push_back(maxValue);
    metrics.push_back(avgValue);
    //--------------------------
    qualityFilter->SetHexQualityMeasureToVolume();
    qualityFilter->Update();
    GetMetrics(qualityFilter, minValue, maxValue, avgValue);
    metrics.push_back(minValue);
    metrics.push_back(avgValue);

    std::ofstream ofs(filename);
    for (size_t i = 0; i < metrics.size(); i++)
        ofs << metrics.at(i) << std::endl;
}
int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " Filename(.vtk) outputFileName(.vtk) opacity metricFileName=<> " << std::endl;
        return EXIT_FAILURE;
    }
    OutputMetricFile(argv[1]);
    // Get the filename from the command line
    std::string inputFilename = argv[1];
    std::string outputFilename = argv[2];

    // Get all data from the file
    vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName(inputFilename.c_str());
    reader->Update();

//    vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
//    triangleFilter->SetInputConnection(reader->GetOutputPort());
//    triangleFilter->Update();
//
//    vtkPolyData* mesh = triangleFilter->GetOutput();
//    std::cout << "There are " << mesh->GetNumberOfCells() << " cells." << std::endl;
    std::cout << "#Vertices = " << reader->GetOutput()->GetNumberOfPoints() << std::endl;
    std::cout << "#Cells = " << reader->GetOutput()->GetNumberOfCells() << std::endl;
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGridMesh = reader->GetOutput();
    vtkSmartPointer<vtkMeshQuality> qualityFilter = vtkSmartPointer<vtkMeshQuality>::New();
#if VTK_MAJOR_VERSION <= 5
    qualityFilter->SetInputConnection(mesh->GetProducerPort());
#else
    ///qualityFilter->SetInputData(mesh);
    qualityFilter->SetInputConnection(reader->GetOutputPort());
#endif
    double minValue = 0;
    double maxValue = 0;
    double avgValue = 0;
    std::vector<double> metrics;
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
    // qualityFilter->SetHexQualityMeasureToScaledJacobian();
    qualityFilter->SetQuadQualityMeasureToScaledJacobian();
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
    size_t numOfInvertedElements = 0;
    // std::cout << "# cells: " << reader->GetOutput()->GetNumberOfCells() << "# jacobians: " << qualityArray->GetNumberOfTuples() << std::endl;
    std::vector<double> qualityValues;
    for (vtkIdType i = 0; i < qualityArray->GetNumberOfTuples(); i++) {
        double val = qualityArray->GetValue(i);
        qualityValues.push_back(val);
        if (minValue > val) minValue = val;
        if (maxValue < val) maxValue = val;
        avgValue += val;

        if (val < 0) numOfInvertedElements++;
        else if (val > 1) numOfInvertedElements++;
        // std::cout << "value " << i << " : " << val << std::endl;
    }
    avgValue /= qualityArray->GetNumberOfTuples();
    std::cout << "min value = " << minValue << std::endl;
    std::cout << "max value = " << maxValue << std::endl;
    std::cout << "avg value = " << avgValue << std::endl;
    std::cout << "#InvertedElements = " << numOfInvertedElements << std::endl;

    MeshFileReader meshReader(inputFilename.c_str());
    Mesh& mesh = (Mesh&)meshReader.GetMesh();
    for (size_t i = 0; i < mesh.C.size(); i++) {
        Cell& c = mesh.C.at(i);
        if (qualityValues.at(i) < 0 ) {
            c.qualityValue = -1;
        } else {
            c.qualityValue = 1;
        }
    }

    const size_t vnum = mesh.V.size();
    const size_t cnum = mesh.C.size();

    std::cout << outputFilename.c_str() << std::endl;

    std::ofstream ofs(outputFilename.c_str());
    ofs << "# vtk DataFile Version 3.0\n"
        << "jacobian" << "\n"
        << "ASCII\n\n"
        << "DATASET UNSTRUCTURED_GRID\n";
    ofs << "POINTS " << vnum << " double\n";
    // ofs << "POINTS " << vnum << " float\n";
    for (size_t i = 0; i < vnum; i++)
        ofs << mesh.V.at(i).x << " " << mesh.V.at(i).y << " " << mesh.V.at(i).z << "\n";
        // ofs << (float) V.at(i).x << " " << (float) V.at(i).y << " " << (float) V.at(i).z << "\n";
    ofs << "CELLS " << cnum << " " << cnum*5 << "\n";

    for (size_t i = 0; i < cnum; i++){
        ofs << mesh.C.at(i).Vids.size();
        for (size_t j = 0; j < mesh.C.at(i).Vids.size(); j++)
            ofs << " " << mesh.C.at(i).Vids.at(j);
        ofs << "\n";
    }
    ofs << "CELL_TYPES " << cnum << "\n";
    for (size_t i = 0; i < cnum; i++) {
        ofs << "9\n";
    }

    ofs << "CELL_DATA " << cnum << "\n";
    ofs << "SCALARS fixed double\n";
    ofs << "LOOKUP_TABLE default\n";

    for (auto& c: mesh.C) {
        ofs << c.qualityValue << "\n";
    }

    if (argc == 3) return 0;

    // Hightlight inverted elements
    vtkDataSet* qualityMesh = qualityFilter->GetOutput();
    vtkSmartPointer<vtkThreshold> selectCells = vtkSmartPointer<vtkThreshold>::New();
    if (minValue < 0) selectCells->ThresholdByLower(0.00);
    else selectCells->ThresholdByLower(minValue + 1e-2);
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
    inverted_actor->GetProperty()->SetColor(0.0, 1.0, 0.0);
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
    //if (numOfInvertedElements != 0)
    actor->GetProperty()->SetOpacity(std::stod(argv[2])/*0.4*/);

    vtkSmartPointer<vtkScalarBarActor> scalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
    scalarBar->SetLookupTable(mapper->GetLookupTable());
    scalarBar->SetTitle("Scaled Jac.");
    scalarBar->SetNumberOfLabels(21);
    scalarBar->GetLabelTextProperty()->SetColor(0, 0, 0);
    scalarBar->GetTitleTextProperty()->SetColor(0, 0, 0);

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
//    if (minValue < 0 && maxValue > 0)
//    {
//        ctf->AddRGBPoint(minValue, 0, 0, 1);
//        ctf->AddRGBPoint(0, 1, 1, 1);
//        ctf->AddRGBPoint(maxValue, 1, 0, 0);
//    }
//    else if (minValue > 0 && maxValue > 0)
//    {
//        ctf->AddRGBPoint(minValue, 1, 1, 1);
//        ctf->AddRGBPoint(maxValue, 1, 0, 0);
//    }
//    else if (minValue < 0 && maxValue < 0)
//    {
//        ctf->AddRGBPoint(minValue, 0, 0, 1);
//        ctf->AddRGBPoint(maxValue, 1, 1, 1);
//    }
    ctf->AddRGBPoint(-1, 1, 0, 0);
    ctf->AddRGBPoint(0, 1, 0, 0);
    ctf->AddRGBPoint(1e-6, 1, 1, 0);
    ctf->AddRGBPoint(0.2, 1, 1, 0);
    ctf->AddRGBPoint(0.2 + 1e-6, 0, 1, 0);
    ctf->AddRGBPoint(0.4, 0, 1, 0);
    ctf->AddRGBPoint(0.4 + 1e-6, 0, 1, 1);
    ctf->AddRGBPoint(0.6, 0, 1, 1);
    ctf->AddRGBPoint(0.6 + 1e-6, 1, 1, 1);
    ctf->AddRGBPoint(1, 1, 1, 1);

    ctf->SetColorSpaceToDiverging();
    ctf->Build();
    mapper->SetLookupTable(ctf);
    scalarBar->SetLookupTable(ctf);

    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);
    vtkSmartPointer<vtkInteractorStyleTrackballCamera> InteractorStyleTrackballCamera = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
    renderWindowInteractor->SetInteractorStyle(InteractorStyleTrackballCamera);

    renderer->AddActor(actor);
    renderer->AddActor2D(scalarBar);
    renderer->AddActor(inverted_actor);
    if (numOfInvertedElements != 0 && std::stod(argv[2]) != 1.0) renderer->AddActor(edges_actor);

    renderer->SetBackground(1, 1, 1);

    ////////////////////////////////////
//    double lightPosition[3] = {0, 0, 1};
//    double lightFocalPoint[3] = {0,0,0};
//    vtkSmartPointer<vtkLight> light = vtkSmartPointer<vtkLight>::New();
//    light->SetLightTypeToSceneLight();
//    light->SetPosition(lightPosition[0], lightPosition[1], lightPosition[2]);
//    light->SetPositional(true); // required for vtkLightActor below
//    light->SetConeAngle(10);
//    light->SetFocalPoint(lightFocalPoint[0], lightFocalPoint[1], lightFocalPoint[2]);
//    light->SetDiffuseColor(1,1,1);
//    light->SetAmbientColor(1,1,1);
//    light->SetSpecularColor(0,0,0);
//    light->SetIntensity(1);
//    vtkSmartPointer<vtkLightActor> lightActor = vtkSmartPointer<vtkLightActor>::New();
//    lightActor->SetLight(light);
//    renderer->AddViewProp(lightActor);
    ////////////////////////////////////

    renderWindow->Render();
    std::string name = std::string("Visualization ToolKit - OpenGL - ") + inputFilename;
    renderWindow->SetWindowName(name.c_str());
    renderWindow->SetSize(1920, 1080);
    renderWindowInteractor->Start();

    return EXIT_SUCCESS;
}
