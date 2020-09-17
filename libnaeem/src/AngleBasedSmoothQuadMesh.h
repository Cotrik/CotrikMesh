#ifndef ANGLE_BASED_SMOOTH_QUAD_MESH_H_
#define ANGLE_BASED_SMOOTH_QUAD_MESH_H_

#include "Mesh.h"
#include <glm/glm.hpp>

class SmoothAlgorithm {

private:
    Mesh& mesh;
    Mesh& boundary_mesh;
    std::vector<glm::dvec3> delta_coords;
    int iters;
    double lambda;
    double tau = 1e-4;
    float min_displacement_limit = 1e-6;
    double interiorE = 0.0;
    double boundaryE = 0.0;

public:
    SmoothAlgorithm(Mesh& mesh, Mesh& boundary_mesh, int it, double l_r);
    ~SmoothAlgorithm();

    std::vector<Vertex> original_vertices;
    std::vector<Vertex> input_vertices;
    std::vector<std::vector<int>> input_boundary; 
    std::vector<std::vector<int>> original_boundary; 

    void setOriginalVertices();
    double getMinEdgeLength();
    void resampleBoundaryVertices();
    void resampleBoundaryVertices1();
    void remapBoundaryVertices();
    void remapToOriginalBoundary();
    void extractBoundary();
    void smoothLaplacianSimple();
    void smoothLaplacianScaleBased();
    void smoothLaplacianCotangentBased();
    void SmoothAngleBased();
    void calculateMeshAngles();
    void angleBasedSmoothing();
    void smoothMesh();
    bool isMeshNonManifold();
    void setBoundaryVerticesMovable();
    void findNegativeElements();
    void fixCrossQuad(int fid);
    double getMeshEnergy(bool consider_boundary);
    double getVertexEnergy(int vid);
    bool isFaceNegative(int fid, int vid, glm::dvec3 false_coord);
};

#endif