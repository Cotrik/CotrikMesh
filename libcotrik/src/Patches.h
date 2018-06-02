/*
 * Patches.h
 *
 *  Created on: May 17, 2017
 *      Author: cotrik
 */

#ifndef LIBCOTRIK_SRC_PATCHES_H_
#define LIBCOTRIK_SRC_PATCHES_H_
#include "Mesh.h"

class Patch
{
public:
    Patch(Mesh& mesh): mesh(mesh) {}
    Patch(const Patch& patch): mesh(patch.mesh), faceIds(patch.faceIds), edgeIds(patch.edgeIds), vertexIds(patch.vertexIds) {}
    Patch& operator=(const Patch& patch) {
        mesh = patch.mesh;
        faceIds = patch.faceIds;
        edgeIds = patch.edgeIds;
        vertexIds = patch.vertexIds;
    }
    ~Patch(){}

    void LabelFace(Face& initialFace, Face& face, size_t& label);
    void SetGlobalCosAngle(const double value = 0.0) {cosangle = value;}
    void WriteMeshFile(const char* filename) const;
private:
    Patch();
public:
    Mesh& mesh;
    std::vector<size_t> faceIds;
    std::vector<size_t> edgeIds;
    std::vector<size_t> vertexIds;
    double cosangle = 0;
};

class Patches
{
public:
    Patches(Mesh& mesh);
    virtual ~Patches();
public:
    void Extract();
    void Extract2();
    void SetGlobalCosAngle(const double value = 0.0) {cosangle = value;}
    void WriteMeshFile(const char* filename_prefix) const;
private:
    Patches();
    void LabelSurface();
    void LabelFace(Face& initialFace, Face& face, size_t& label);
    void ExtractPatches();
public:
    Mesh& mesh;
    std::vector<Patch> patches;
    double cosangle = 0;
};

#endif /* LIBCOTRIK_SRC_PATCHES_H_ */
