/*
* QuadSurfaceMapper.h
*
*  Created on: September 24, 2021
*      Author: https://github.com/naeem014
*/

#ifndef QUAD_SURFACE_MAPPER_H_
#define QUAD_SURFACE_MAPPER_H_

#include <vtkPolyData.h>
#include <vtkCellLocator.h>
#include <glm/glm.hpp>

#include "Mesh.h"
#include "MeshUtil.h"

class SurfaceMapper {
    public:
        SurfaceMapper();
        SurfaceMapper(Mesh& source_, Mesh& target_);
        SurfaceMapper(Mesh& target_);
        ~SurfaceMapper();

        void SetSource(Mesh& mesh);
        void SetTarget(Mesh& mesh);
        glm::dvec3 GetClosestPoint(glm::dvec3 p);
        
    private:
        Mesh& source = Mesh();
        Mesh& target = Mesh();
        MeshUtil mu;
        vtkSmartPointer<vtkPolyData> source_polyData = vtkSmartPointer<vtkPolyData>::New();
        vtkSmartPointer<vtkPolyData> target_polyData = vtkSmartPointer<vtkPolyData>::New();
        vtkSmartPointer<vtkPolyData> target_cellPolyData = vtkSmartPointer<vtkPolyData>::New();
        vtkSmartPointer<vtkCellLocator> cell_locator = vtkSmartPointer<vtkCellLocator>::New();
        
        vtkSmartPointer<vtkPolyData> GetPolyDataFromMesh(Mesh& mesh);
        void SetCellLocator();
};

#endif