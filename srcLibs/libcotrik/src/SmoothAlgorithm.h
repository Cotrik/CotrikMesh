/*
 * SmoothAlgorithm.h
 *
 *  Created on: Mar 29, 2017
 *      Author: cotrik
 */

#ifndef SMOOTH_ALGORITHM_H_
#define SMOOTH_ALGORITHM_H_

#include "Mesh.h"

enum Smooth_Algorithm {
    LAPLACIAN = 0,
    SCALED_JACOBIAN,
    ANGLE_BASED
};

class SmoothAlgorithm
{
public:
    SmoothAlgorithm(Mesh& mesh, const Smooth_Algorithm smoothAlgorithm= LAPLACIAN);
    virtual ~SmoothAlgorithm();

private:
    SmoothAlgorithm();
    SmoothAlgorithm(const SmoothAlgorithm&);
    SmoothAlgorithm& operator = (const SmoothAlgorithm&) const;

public:
    virtual int Run(const size_t iters = 1e2, const double eps = 1e-4);

private:
    double SmoothVolume(const Smooth_Algorithm smoothMethod = LAPLACIAN);
    glm::dvec3 LapLace(const Vertex& v);
    glm::dvec3 ScaledJacobian(const Vertex& v);

private:
    Mesh& mesh;
    Smooth_Algorithm m_smoothAlgorithm = LAPLACIAN;
};

#endif /* SMOOTH_ALGORITHM_H_ */
