/*
 * MeshOpt.h
 *
 *  Created on: Nov 25, 2016
 *      Author: cotrik
 */

#ifndef AUTO_MESH_OPT_H_
#define AUTO_MESH_OPT_H_

#include "Mesh.h"
#include "LocalMeshOpt.h"

class AutoMeshOpt : public LocalMeshOpt
{
public:
    AutoMeshOpt(const Mesh& mesh, const Mesh& origMesh);
    virtual ~AutoMeshOpt();

protected:
    AutoMeshOpt();
    AutoMeshOpt(const AutoMeshOpt& meshOpt);

public:
    virtual void Run();
    void SetHausdorffError(const double value = 1e-2);

protected:
    double GetHausdorffError() const;
    void InversionFreeDeformToTargetMesh();
protected:
    const Mesh& origMesh;
    double m_hausdorffError;
};

#endif /* AUTO_MESH_OPT_H_ */
