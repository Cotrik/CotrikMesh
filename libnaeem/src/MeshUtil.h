/*
* MeshUtil.h
*
*  Created on: October 18, 2021
*      Author: https://github.com/naeem014
*/

#ifndef MESH_UTIL_H_
#define MESH_UTIL_H_

#include <vtkPolyData.h>
#include <glm/glm.hpp>

#include "Mesh.h"

class MeshUtil {
    public:
        // Constructors and Destructor
        MeshUtil();
        MeshUtil(const MeshUtil& r);
        MeshUtil(Mesh& mesh_);
        ~MeshUtil();

        MeshUtil& operator = (const MeshUtil&);

        // MeshUtil setters and getters
        void SetMembers(Mesh& mesh_);

        // Mesh Utils
        vtkSmartPointer<vtkPolyData> GetPolyData();
        vtkSmartPointer<vtkPolyData> GetPolyData(Mesh& mesh_, size_t vid);
        vtkSmartPointer<vtkPolyData> GetPolyData(Mesh& mesh_, std::vector<size_t> V);
        void SetMeshArea();
        double GetMeshArea();

        // Face Utils
        double GetFaceArea(int fid);

        // Vertex Utils
        double GetVertexEnergy(int vid);
        double GetInteriorAngleAtEdge(int vid, int eid);

        // Misc
        void AddContents(std::vector<size_t>& a, std::vector<size_t> b);
        void UpdateContents(std::vector<size_t>& a, std::vector<size_t> b);
        std::vector<size_t> GetDifference(std::vector<size_t> a, std::vector<size_t> b);
        std::vector<size_t> GetUnion(std::vector<size_t> a, std::vector<size_t> b);
        std::vector<size_t> GetIntersection(std::vector<size_t> a, std::vector<size_t> b);
        std::vector<size_t> GetIntersectionParallel(std::vector<size_t> a, std::vector<size_t> b);
        
    private:
        Mesh* mesh;

        void CheckValidity();
};

#endif