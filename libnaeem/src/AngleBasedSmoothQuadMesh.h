#ifndef ANGLE_BASED_SMOOTH_QUAD_MESH_H_
#define ANGLE_BASED_SMOOTH_QUAD_MESH_H_

#include "Mesh.h"
#include <glm/glm.hpp>

class SmoothAlgorithm {

private:
    Mesh& mesh;
    std::vector<glm::dvec3> delta_coords;
    int iters;
    double lambda;
    float min_displacement_limit = 1e-6;

public:
    SmoothAlgorithm(Mesh& mesh, int it, double l_r);
    ~SmoothAlgorithm();

    std::vector<Vertex> original_vertices;


    void setOriginalVertices();
    double getMinEdgeLength();
    void smoothBoundary();
    void resampleBoundaryVertices();
    void smoothLaplacianSimple();
    void smoothLaplacianScaleBased();
    void smoothLaplacianCotangentBased();
    void SmoothAngleBased();
    void calculateMeshAngles();
    void angleBasedSmoothing();
    void smoothMesh();
    void remapBoundaryVertices();
};

#endif