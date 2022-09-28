#ifndef GEODESIC_DISTANCE_H
#define GEODESIC_DISTANCE_H

#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "hmTriDistance.h"
#include "hmContext.h"
#include "hmUtility.h"
#include "hmVectorSizeT.h"

class GeodesicDistance{
public:
    GeodesicDistance();
    ~GeodesicDistance();
    int parseCommandLineArguments( int argc, char** argv, char** inputMesh, char** inputSources, char** outputMeshPrefix, char** outputDistancePrefix, char** inputReference, double* smoothness, double* boundaryConditions, char* verbose );
    void printHelp( int argc, char** argv );
    void printVersion( int argc, char** argv );
    void readMesh( hmTriDistance* distance, hmTriMesh* mesh, const char* filename );
    void readSources( int* nSourceSets, hmVectorSizeT** sourceSets, const char* filename );
    void readReferenceValues( hmTriDistance* distance, int* nReferenceColumns, double*** referenceValues, char*** referenceNames, const char* filename );
    void destroySourceSets( int nSourceSets, hmVectorSizeT** sourceSets );
    void destroyReferenceValues( int nReferenceColumns, double*** referenceValues, char*** referenceNames );
    void writeMesh( hmTriDistance* distance, const char* filename, int index );
    void writeDistances( hmTriDistance* distance, const char* prefix, int index );
    void setSources( hmTriDistance* distance, hmVectorSizeT* sourceSets, int index );
    void compareDistance( hmTriDistance* distance, int nReferenceColumns, double** referenceValues, char** referenceNames, const char* prefix );

    void Init();
    void Destroy();
    double GetGeodesicDistance(size_t vid1, size_t vid2);

    /* main data */
    hmContext context;
    hmTriMesh surface;
    hmTriDistance distance;

    size_t nSourceSets;// = 0;
    hmVectorSizeT* sourceSets;// = NULL;

    int nReferenceColumns;// = 0;
    double** referenceValues;// = NULL;
    char** referenceNames;// = NULL;

    /* parameters parsed from command-line arguments */
    char* inputMesh;// = NULL;
    char* inputSources;// = NULL;
    char* inputReference;// = NULL;
    char* outputMeshPrefix;// = NULL;
    char* outputDistancePrefix;// = NULL;
    double smoothness;// = -1.;
    double boundaryConditions;// = -1.;
    char verbose;// = 0;
};

//extern GeodesicDistance gd;

#endif //GEODESIC_DISTANCE_H
