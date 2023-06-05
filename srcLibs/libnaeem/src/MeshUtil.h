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


struct VertexScore {
    // glm::dmat4 Q;
    int same_singularity_min_dist = 3;
    int boundary_min_dist = 3;
    int diff_singularity_max_dist = std::numeric_limits<int>::max();
    int diff_count = 0;
    int same_count = 0;
    double boundary_score = 1;
    double same_singularity_score = 1;
    double diff_singularity_score = 0;
    double qem_score = 0;
};

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
        bool IsSharpFeature(size_t vid);

        size_t GetEdgeId(size_t vid, size_t vid2);
        size_t GetCCedgeAt(size_t vid, size_t eid, int counter);
        int GetVIdx(size_t vid, size_t fid);
        size_t GetFaceV(size_t fid, int idx, int offset);
        size_t GettCCFaceAt(size_t vid, size_t eid, int counter);
        
        void SetV_Scores();
        void SetQEMs();
        double CalculateQEM(size_t vid);
        double GetQEMcost(size_t vid);
        VertexScore* GetVertexScore(size_t vid, bool calculate_qem = false);
        std::vector<size_t> GetVertexDirections(size_t vid, std::vector<size_t> exceptions);



        double Q_A = 0.0;
        double delta = 0.0;
    private:
        Mesh* mesh;
        std::vector<VertexScore> v_scores;

        void CheckValidity();
};

#endif