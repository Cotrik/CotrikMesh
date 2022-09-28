# include(CheckCXXCompilerFlag)

# set(CMAKE_CXX_STANDARD 14)
# set(CMAKE_CXX_STANDARD_REQUIRED ON)

# # CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
# # CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
# # if(COMPILER_SUPPORTS_CXX11)
# #     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
# # elseif(COMPILER_SUPPORTS_CXX0X)
# #     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
# # else()
# #         message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
# # endif()

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

# #find_package(OpenGL REQUIRED)
# #find_package(Boost QUIET COMPONENTS program_options)
# #if(Boost_PROGRAM_OPTIONS_FOUND)
# #    if( CGAL_AUTO_LINK_ENABLED )
# #        message( STATUS "Boost.ProgramOptions library: found" )
# #    else()
# #        message( STATUS "Boost.ProgramOptions library: ${Boost_PROGRAM_OPTIONS_LIBRARY}" )
# #    endif()
# #    add_definitions( "-DCGAL_USE_BOOST_PROGRAM_OPTIONS" )
# #    list(APPEND CGAL_3RD_PARTY_LIBRARIES ${Boost_LIBRARIES})
# #endif()

# find_package(VTK REQUIRED)
# include(${VTK_USE_FILE})   
# include_directories(
#     ${PROJECT_SOURCE_DIR}/eigen3/
#     ${PROJECT_SOURCE_DIR}/libnaeem/src/
#     #/usr/local/include/
#     #/usr/local/include/vtk-6.3/
#     ${PROJECT_SOURCE_DIR}/verdict/vtkverdict/
#     ${PROJECT_SOURCE_DIR}/libigl/include/
#     ../
#     #${PROJECT_SOURCE_DIR}/libgeodesic_heart/include
#     ${VTK_INCLUDE_DIRS}
# )

# set(ALL_LIBS
# #       ${OPENGL_LIBRARY}
# #       GL
# #    boost_program_options
# #    pthread
# #    gmp
#     ${VTK_LIBRARIES}
# #    viennacl
# )

# add_definitions(
# #       -D_CRT_SECURE_NO_WARNINGS
# )
# add_library(cotrik SHARED STATIC 
#     #src/AutoMeshOpt.cpp
#     #src/AutoMeshOpt.h
#     src/BaseComplex.cpp
#     src/BaseComplex.h
#     src/BaseComplexQuad.cpp
#     src/BaseComplexQuad.h
#     src/BaseComplexEditor.cpp
#     src/BaseComplexEditor.h
#     src/BaseComplexSheet.cpp
#     src/BaseComplexSheet.h
#     src/BaseComplexSheetQuad.cpp
#     src/BaseComplexSheetQuad.h
#     src/BaseComplexChord.cpp
#     src/BaseComplexChord.h
#     src/BaseComplexCleaner.cpp
#     src/BaseComplexCleaner.h
#     src/SingularityGraph.cpp
#     src/SingularityGraph.h
#     src/Component.cpp
#     src/Component.h
#     src/Dual.cpp
#     src/Dual.h
#     src/RefinedDual.cpp
#     src/RefinedDual.h
#     src/RefinedDualQuad.cpp
#     src/RefinedDualQuad.h
#     src/EdgeLines.cpp
#     src/EdgeLines.h
#     src/FeatureLine.cpp
#     src/FeatureLine.h
#     src/Frame.cpp
#     src/FrameField.cpp
#     src/FrameField.h
#     src/Frame.h
#     src/FrameOpt.cpp
#     src/FrameOpt.h
#     #src/GeodesicDistance.cpp
#     #src/GeodesicDistance.h
#     #src/LayerOpt.cpp
#     #src/LayerOpt.h
#     src/MeshOpt.cpp
#     src/MeshOpt.h
#     src/LocalMeshOptASJ.cpp
#     src/LocalMeshOptASJ.h
#     src/LocalMeshOpt.cpp
#     src/LocalMeshOpt.h
#     src/LocalMeshOptFixBoundary.cpp
#     src/LocalMeshOptFixBoundary.h
#     src/MeshASJOpt.cpp
#     src/MeshASJOpt.h
#     src/Mesh.cpp
#     src/Mesh.h
#     src/MeshFileReader.cpp
#     src/MeshFileReader.h
#     src/MeshFileWriter.cpp
#     src/MeshFileWriter.h
#     src/UnstructuredVTKWriter.cpp
#     src/UnstructuredVTKWriter.h
#     src/MeshOptFixBoundary.cpp
#     src/MeshOptFixBoundary.h
#     src/MeshQuality.cpp
#     src/MeshQuality.h
#     src/Patches.cpp
#     src/Patches.h
#     src/PolyLine.cpp
#     src/PolyLine.h
#     src/PhongProjector.cpp
#     src/PhongProjector.h
#     src/SmoothAlgorithm.cpp
#     src/SmoothAlgorithm.h
#     src/SurfaceMeshOpt.cpp
#     src/SurfaceMeshOpt.h
#     src/WeightedMeshOptFixBoundary.cpp
#     src/WeightedMeshOptFixBoundary.h
#     src/UnstructuredVTKWriter.h
#     src/UnstructuredVTKWriter.cpp
# 	src/DualMesh.h
# 	src/DualMesh.cpp
# 	src/MST.h
# 	src/MST.cpp
# 	src/Simplifier.h
# 	src/Simplifier.cpp
# 	src/CollapseSimplifier.h
# 	src/CollapseSimplifier.cpp
# 	src/SplitSimplifier.h
# 	src/SplitSimplifier.cpp
# 	src/SheetSimplifier.h
# 	src/SheetSimplifier.cpp
# 	src/GlobalSheetSimplifier.h
# 	src/GlobalSheetSimplifier.cpp
# 	src/SingleSheetSimplifier.h
# 	src/SingleSheetSimplifier.cpp
# 	src/SheetSplitSimplifier.h
# 	src/SheetSplitSimplifier.cpp
# 	src/DoubletSimplifier.h
# 	src/DoubletSimplifier.cpp
# 	src/EdgeRotateSimplifier.h
# 	src/EdgeRotateSimplifier.cpp
# 	src/DiagnalCollapseSimplifier.h
# 	src/DiagnalCollapseSimplifier.cpp
# 	src/TriangleSimplifier.h
# 	src/TriangleSimplifier.cpp
# 	src/PatchSimplifier.h
# 	src/PatchSimplifier.cpp
# 	src/Util.h
# 	src/Util.cpp
#     src/PatchSimplifierOperations.h
#     src/LocalSimplifier.h
#     src/LocalSimplifier.cpp
# )
# target_link_libraries(cotrik ${ALL_LIBS})
# #install(FILES cotrik DESTINATION D:/Program\ Files/CotrikMesh/)
