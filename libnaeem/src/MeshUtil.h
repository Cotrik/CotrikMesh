/*
* MeshUtil.h
*
*  Created on: October 18, 2021
*      Author: https://github.com/naeem014
*/

#ifndef MESH_UTIL_H_
#define MESH_UTIL_H_

#include <glm/glm.hpp>

#include "Mesh.h"

class MeshUtil {
    public:
        // Constructors and Destructor
        MeshUtil();
        MeshUtil(Mesh& mesh_);
        ~MeshUtil();

        MeshUtil& operator = (const MeshUtil&);

        // MeshUtil setters and getters
        void SetMesh(Mesh& mesh_);

        // Mesh Utils
        void SetMeshArea();
        double GetMeshArea();

        // Face Utils
        double GetFaceArea(int fid);

        // Vertex Utils
        double GetVertexEnergy(int vid);
        double GetInteriorAngleAtEdge(int vid, int eid);
    private:
        Mesh& mesh = Mesh();
        
        void CheckValidity();
};

#endif