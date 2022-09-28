/*
 * MeshASJOpt.h
 *
 *  Created on: March 23, 2017
 *      Author: cotrik
 */

#ifndef MESH_ASJ_OPT_H_
#define MESH_ASJ_OPT_H_

#include "MeshOpt.h"
class MeshASJOpt: public MeshOpt
{
public:
    MeshASJOpt(const Mesh& mesh);
    virtual ~MeshASJOpt();

protected:
    MeshASJOpt();
    MeshASJOpt(const MeshASJOpt& meshOpt);

public:
    virtual size_t Run(const size_t iters = 1);
    bool Optimize(double& energy);

private:
    double currentMSJ = 0.0;
    double currentASJ = 0.0;
};

#endif /* MESH_ASJ_OPT_H_ */
