/*
* FeatureExtractor.h
*
*  Created on: November 16, 2021
*      Author: https://github.com/naeem014
*/

#ifndef FEATURE_EXTRACTOR_H_
#define FEATURE_EXTRACTOR_H_

#include <thread>
#include <vector>

#include <vtkPolyData.h>
#include <glm/glm.hpp>


#include <vtkFeatureEdges.h>
#include <vtkPointData.h>
#include <vtkPointLocator.h>
#include <vtkImplicitPolyDataDistance.h>

#include "Mesh.h"
#include "MeshUtil.h"
// #include "ParallelFor.h"

#define NUM_THREADS std::thread::hardware_concurrency() == 0 ? 8 : std::thread::hardware_concurrency()

class FeatureExtractor {
    public:
        FeatureExtractor();
        FeatureExtractor(Mesh& mesh_, double angle_threshold_);
        ~FeatureExtractor();

        void SetMesh(Mesh& mesh_);
        void Extract();
    private:
        void SetFeatures(int i, vtkSmartPointer<vtkPoints> res, vtkSmartPointer<vtkPointLocator> pl, bool boundary);
        // void SetFeatures(int i, vtkSmartPointer<vtkPoints> points, bool boundary);

        Mesh& mesh = Mesh();
        MeshUtil mu;
        vtkSmartPointer<vtkPolyData> mesh_polyData = vtkSmartPointer<vtkPolyData>::New();
        vtkSmartPointer<vtkPolyData> GetPolyDataFromMesh();
        double angle_threshold = 15.0;
};

#endif