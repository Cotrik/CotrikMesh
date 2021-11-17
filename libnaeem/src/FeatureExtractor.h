/*
* FeatureExtractor.h
*
*  Created on: November 16, 2021
*      Author: https://github.com/naeem014
*/

#ifndef FEATURE_EXTRACTOR_H_
#define FEATURE_EXTRACTOR_H_

#include <vtkPolyData.h>
#include <glm/glm.hpp>

#include "Mesh.h"
#include "MeshUtil.h"

class FeatureExtractor {
    public:
        FeatureExtractor();
        FeatureExtractor(Mesh& mesh_, double angle_threshold_);
        ~FeatureExtractor();

        void SetMesh(Mesh& mesh_);
        void Extract();
    private:
        Mesh& mesh = Mesh();
        MeshUtil mu;
        vtkSmartPointer<vtkPolyData> mesh_polyData = vtkSmartPointer<vtkPolyData>::New();
        vtkSmartPointer<vtkPolyData> GetPolyDataFromMesh();
        double angle_threshold = 15.0;
};

#endif