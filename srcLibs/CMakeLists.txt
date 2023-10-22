include(CheckCXXCompilerFlag)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

find_package(VTK REQUIRED)
find_package(CGAL REQUIRED)
include(${CGAL_USE_FILE})
# include(${VTK_USE_FILE})   
include_directories(
    ${PROJECT_SOURCE_DIR}/srcLibs/libcotrik/src/ 
    ${PROJECT_SOURCE_DIR}/srcLibs/libnaeem/src/
    ${PROJECT_SOURCE_DIR}/eigen3/
    ${PROJECT_SOURCE_DIR}/verdict/vtkverdict/
    ${PROJECT_SOURCE_DIR}/libigl/include/
    ../
    ${VTK_INCLUDE_DIRS}
    ${CGAL_INCLUDE_DIRS}
)

set(ALL_LIBS
    ${VTK_LIBRARIES}
    ${CGAL_LIBRARIES}
)

add_library(srcLibs SHARED STATIC

    # Source files from libcotrik
    libcotrik/src/BaseComplex.cpp
    libcotrik/src/BaseComplex.h
    libcotrik/src/BaseComplexQuad.cpp
    libcotrik/src/BaseComplexQuad.h
    libcotrik/src/BaseComplexEditor.cpp
    libcotrik/src/BaseComplexEditor.h
    libcotrik/src/BaseComplexSheet.cpp
    libcotrik/src/BaseComplexSheet.h
    libcotrik/src/BaseComplexSheetQuad.cpp
    libcotrik/src/BaseComplexSheetQuad.h
    libcotrik/src/BaseComplexChord.cpp
    libcotrik/src/BaseComplexChord.h
    libcotrik/src/BaseComplexCleaner.cpp
    libcotrik/src/BaseComplexCleaner.h
    libcotrik/src/SingularityGraph.cpp
    libcotrik/src/SingularityGraph.h
    libcotrik/src/Component.cpp
    libcotrik/src/Component.h
    libcotrik/src/Dual.cpp
    libcotrik/src/Dual.h
    libcotrik/src/RefinedDual.cpp
    libcotrik/src/RefinedDual.h
    libcotrik/src/RefinedDualQuad.cpp
    libcotrik/src/RefinedDualQuad.h
    libcotrik/src/EdgeLines.cpp
    libcotrik/src/EdgeLines.h
    libcotrik/src/FeatureLine.cpp
    libcotrik/src/FeatureLine.h
    libcotrik/src/Frame.cpp
    libcotrik/src/FrameField.cpp
    libcotrik/src/FrameField.h
    libcotrik/src/Frame.h
    libcotrik/src/FrameOpt.cpp
    libcotrik/src/FrameOpt.h
    libcotrik/src/MeshOpt.cpp
    libcotrik/src/MeshOpt.h
    libcotrik/src/LocalMeshOptASJ.cpp
    libcotrik/src/LocalMeshOptASJ.h
    libcotrik/src/LocalMeshOpt.cpp
    libcotrik/src/LocalMeshOpt.h
    libcotrik/src/LocalMeshOptFixBoundary.cpp
    libcotrik/src/LocalMeshOptFixBoundary.h
    libcotrik/src/MeshASJOpt.cpp
    libcotrik/src/MeshASJOpt.h
    libcotrik/src/Mesh.cpp
    libcotrik/src/Mesh.h
    libcotrik/src/MeshFileReader.cpp
    libcotrik/src/MeshFileReader.h
    libcotrik/src/MeshFileWriter.cpp
    libcotrik/src/MeshFileWriter.h
    libcotrik/src/UnstructuredVTKWriter.cpp
    libcotrik/src/UnstructuredVTKWriter.h
    libcotrik/src/MeshOptFixBoundary.cpp
    libcotrik/src/MeshOptFixBoundary.h
    libcotrik/src/MeshQuality.cpp
    libcotrik/src/MeshQuality.h
    libcotrik/src/Patches.cpp
    libcotrik/src/Patches.h
    libcotrik/src/PolyLine.cpp
    libcotrik/src/PolyLine.h
    libcotrik/src/PhongProjector.cpp
    libcotrik/src/PhongProjector.h
    libcotrik/src/SmoothAlgorithm.cpp
    libcotrik/src/SmoothAlgorithm.h
    libcotrik/src/SurfaceMeshOpt.cpp
    libcotrik/src/SurfaceMeshOpt.h
    libcotrik/src/WeightedMeshOptFixBoundary.cpp
    libcotrik/src/WeightedMeshOptFixBoundary.h
    libcotrik/src/UnstructuredVTKWriter.h
    libcotrik/src/UnstructuredVTKWriter.cpp
    libcotrik/src/DualMesh.h
    libcotrik/src/DualMesh.cpp
    libcotrik/src/MST.h
    libcotrik/src/MST.cpp
    libcotrik/src/Simplifier.h
    libcotrik/src/Simplifier.cpp
    libcotrik/src/CollapseSimplifier.h
    libcotrik/src/CollapseSimplifier.cpp
    libcotrik/src/SplitSimplifier.h
    libcotrik/src/SplitSimplifier.cpp
    libcotrik/src/SheetSimplifier.h
    libcotrik/src/SheetSimplifier.cpp
    libcotrik/src/GlobalSheetSimplifier.h
    libcotrik/src/GlobalSheetSimplifier.cpp
    libcotrik/src/SingleSheetSimplifier.h
    libcotrik/src/SingleSheetSimplifier.cpp
    libcotrik/src/SheetSplitSimplifier.h
    libcotrik/src/SheetSplitSimplifier.cpp
    libcotrik/src/DoubletSimplifier.h
    libcotrik/src/DoubletSimplifier.cpp
    libcotrik/src/EdgeRotateSimplifier.h
    libcotrik/src/EdgeRotateSimplifier.cpp
    libcotrik/src/DiagnalCollapseSimplifier.h
    libcotrik/src/DiagnalCollapseSimplifier.cpp
    libcotrik/src/TriangleSimplifier.h
    libcotrik/src/TriangleSimplifier.cpp
    libcotrik/src/PatchSimplifier.h
    libcotrik/src/PatchSimplifier.cpp
    libcotrik/src/Util.h
    libcotrik/src/Util.cpp
    libcotrik/src/PatchSimplifierOperations.h
    libcotrik/src/LocalSimplifier.h
    libcotrik/src/LocalSimplifier.cpp

    # Source files from libnaeem
    libnaeem/src/AngleBasedSmoothQuadMesh.h
    libnaeem/src/AngleBasedSmoothQuadMesh.cpp
    libnaeem/src/ChordCollapse.h
    libnaeem/src/ChordCollapse.cpp
    libnaeem/src/ChordExtractor.h
    libnaeem/src/ChordExtractor.cpp
    libnaeem/src/DiagonalCollapse.h
    libnaeem/src/DiagonalCollapse.cpp
    libnaeem/src/DiagonalThreeFivePair.h
    libnaeem/src/DiagonalThreeFivePair.cpp
    libnaeem/src/DirectSeparatrixCollapse.h
    libnaeem/src/DirectSeparatrixCollapse.cpp
    libnaeem/src/EdgeCollapse.h
    libnaeem/src/EdgeCollapse.cpp
    libnaeem/src/EdgeRotation.h
    libnaeem/src/EdgeRotation.cpp
    libnaeem/src/FeatureExtractor.h
    libnaeem/src/FeatureExtractor.cpp
    libnaeem/src/KDTree.h
    libnaeem/src/KDTree.cpp
    libnaeem/src/Memento.h
    libnaeem/src/MeshUtil.h
    libnaeem/src/MeshUtil.cpp
    libnaeem/src/ParallelFor.h
    libnaeem/src/PQueue.h
    libnaeem/src/PQueue.cpp
    libnaeem/src/QuadSplit.h
    libnaeem/src/QuadSplit.cpp
    libnaeem/src/QuadSurfaceMapper.h
    libnaeem/src/QuadSurfaceMapper.cpp
    libnaeem/src/Renderer.h
    libnaeem/src/SemiGlobalSimplifier.h
    libnaeem/src/SemiGlobalSimplifier.cpp
    libnaeem/src/SeparatrixCollapse.h
    libnaeem/src/SeparatrixCollapse.cpp
    libnaeem/src/SimplificationOperation.h
    libnaeem/src/SimplificationOperation.cpp
    libnaeem/src/SingularityPair.h
    libnaeem/src/Smooth.h
    libnaeem/src/Smooth.cpp
    libnaeem/src/ThreadPool.h
    libnaeem/src/ThreeFivePair.h
    libnaeem/src/ThreeFivePair.cpp
    libnaeem/src/VertexRotation.h
    libnaeem/src/VertexRotation.cpp
    libnaeem/src/VertexSplit.h
    libnaeem/src/VertexSplit.cpp
)
target_link_libraries(srcLibs ${ALL_LIBS})
#install(FILES cotrik DESTINATION D:/Program\ Files/CotrikMesh/)
