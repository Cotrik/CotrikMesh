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
#include <vtkImplicitPolyDataDistance.h>
#include <vtkDistancePolyDataFilter.h>
#include <vtkTriangleFilter.h>
// #include <vtkHausdorffDistancePointSetFilter.h>
#include <glm/glm.hpp>

#include "Mesh.h"
#include "MeshUtil.h"

class SurfaceMapper {
    public:
        SurfaceMapper();
        SurfaceMapper(Mesh& source_, Mesh& target_, MeshUtil& mu_);
        SurfaceMapper(Mesh& target_, MeshUtil& mu_);
        ~SurfaceMapper();

        void SetSource(Mesh& mesh);
        void SetTarget(Mesh& mesh);
        glm::dvec3 GetClosestPoint(glm::dvec3 p);
        void ExecutePolyDistanceFilter(Mesh& mesh);
        // double ExecuteHaursdorffDsitanceFilter(Mesh& m1, Mesh& m2);
        void RemapVertex(Mesh& mesh_, size_t vid, glm::dvec3 c);
        void SetLocator(Mesh& mesh_, std::vector<size_t> V);
    private:
        // Mesh& source = Mesh(std::vector<Vertex>{}, std::vector<Cell>{}, QUAD);
        // Mesh& target = Mesh();
        Mesh* source;
        Mesh* target;
        MeshUtil* mu;
        vtkSmartPointer<vtkPolyData> source_polyData = vtkSmartPointer<vtkPolyData>::New();
        vtkSmartPointer<vtkPolyData> target_polyData = vtkSmartPointer<vtkPolyData>::New();
        vtkSmartPointer<vtkPolyData> target_cellPolyData = vtkSmartPointer<vtkPolyData>::New();
        vtkSmartPointer<vtkCellLocator> cell_locator = vtkSmartPointer<vtkCellLocator>::New();
        vtkSmartPointer<vtkImplicitPolyDataDistance> point_finder = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
        vtkSmartPointer<vtkDistancePolyDataFilter> polyDistanceFilter = vtkSmartPointer<vtkDistancePolyDataFilter>::New();
        // vtkSmartPointer<vtkHausdorffDistancePointSetFilter> haursdorffDistanceFilter = vtkSmartPointer<vtkHausdorffDistancePointSetFilter>::New();

        vtkSmartPointer<vtkPolyData> poly = vtkSmartPointer<vtkPolyData>::New();
        vtkNew<vtkPoints> points;
        vtkNew<vtkCellArray> cells;
        std::vector<size_t> vertices;
        std::vector<size_t> faces;
        vtkNew<vtkPolygon> p;
        // vtkSmartPointer<vtkTriangleFilter> tf = vtkSmartPointer<vtkTriangleFilter>::New();
        vtkSmartPointer<vtkCellLocator> cl = vtkSmartPointer<vtkCellLocator>::New();
        vtkSmartPointer<vtkPolyData> tp = vtkSmartPointer<vtkPolyData>::New();

        vtkSmartPointer<vtkPolyData> GetPolyDataFromMesh(Mesh& mesh);
        void SetCellLocator();
};

#endif