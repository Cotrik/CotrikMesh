include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

add_definitions(
       -D_CRT_SECURE_NO_WARNINGS
	   -D_USE_MATH_DEFINES
)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})  

include_directories(
    ${PROJECT_SOURCE_DIR}/libcotrik/src/ 
    ${PROJECT_SOURCE_DIR}/libnaeem/src/ 
    ${PROJECT_SOURCE_DIR}/eigen3
    ${PROJECT_SOURCE_DIR}/libigl/include
    ${PROJECT_SOURCE_DIR}/libigl/external
    ../
    )
    
link_directories(/usr/local/lib)    

foreach(app 
    Align
    #AutoMeshOpt
    CheckMesh
    ClipMesh
    ConvertFileFormat
    SliceMesh
    MeasureMesh
    VerifyHexMesh
    CheckManifold
    CleanStructure
    ColorHexMesh
    ConvertQuadToTriangle
    #CotrikClustering
    DisplayFrameField
    ExtractBaseMesh
    ExtractBaseComplex
    ExtractBaseComplexQuad
    ExtractSheetDecompositions
    ExtractSheetDecompositionsQuad
    ExtractMainSheets
    RemovePadding
    CollapseComponent
    ExtractFaceSegments
    ExtractFrameField
    ExtractGaussianCurvature
    ExtractLayers
    ExtractPatches
    ExtractQuadMesh
    ExtractSharpEdges
    ExtractSurface
    ExtractSingularities
    ExtractSingularitiesGraphs
    Fix2DMeshVertices
    HexRefine
    HexSmooth
    #HexGen
    #LayerOpt
    LocalMeshOpt
    LocalMeshOptFixBoundary
    MeshAsjOpt
    #MeshChecker
    #MeshFix
    MeshOpt
    MeshOptFixBoundary
    OpenQuadMeshOpt
	QuadMeshOpt2D
    MeshProject
    MeshQuality
    MeshScale
    MeshShrink
	MeshMap
    #MIQ
    PerturbInnerVertices
    #QuadGen
    ReadMesh
    RemapMeshToTargetMesh
    #Slim
    #SLIM
    SlimHex
    SlimTet
    SlimTri
    #Smooth
    SmoothMap
    SmoothMesh
    SurfaceMeshOpt
    WeightedMeshOptFixBoundary
    tet2hex
    TriangleSmooth
    #Triangulate
    #IsotropicRemesh
    Test
    Rotate
    FixedSmooth
    PolycubeOpt
    PolycubeProcess
    QuadTriMeshRefine
    QuadMeshSimplify
    QuadMeshLocalSimplify
	QuadMeshLocalSimplify2
	QuadMeshLocalSimplify3
	QuadPatchSimplify
    QuadPatchLocalSimplify
	QuadPatchSimplify_refine
	ClosedQuadMeshSimplify
    Vol2Vtk
    ExtractLocalRegion
	HexMeshLocalSimplify
	ExtractBaseComplexSeparatedSurfaces
    KillHighValenceSingularities
	m2vtk
	uv2vtk
	vtk2m
	Quad2Tri
	Tri2Quad
	MapLoops
	ExtractMinSpanTree
	ExtractDualMesh
	ExtractLoops
	ExtractBaseHomology
	ExtractBaseHomologyQuad
	ExtractSurfaceFeature
	ExtractCoverSpace
	ExtractCoverSpace2
	ExtractCoverSpace3
	ExtractCutGraph
	ExtractHomotopyLoops
	ExtractHomotopyBases3
	ConvexCut
)
    add_executable(${app} ${app}.cpp)
    target_link_libraries(${app} cotrik naeem
    #CoMISo pthread
    )	#install(FILES ${app} DESTINATION D:/Program\ Files/CotrikMesh/)
endforeach()

#option(LIBIGL_WITH_COMISO           "Use CoMiso"         ON)
##Compile CoMISo
## NOTE: this cmakefile works only with the
## comiso available here: https://github.com/libigl/CoMISo
#if(LIBIGL_WITH_COMISO)
#  set(COMISO_DIR "${PROJECT_SOURCE_DIR}/../libigl/external/CoMISo")
#  set(COMISO_INCLUDE_DIRS
#    "${COMISO_DIR}/ext/gmm-4.2/include"
#    "${COMISO_DIR}/../")
#  list(APPEND LIBIGL_INCLUDE_DIRS ${COMISO_INCLUDE_DIRS})
#  #add_definitions(-DINCLUDE_TEMPLATES) (what need this?)
#  list(APPEND LIBIGL_DEFINITIONS "-DINCLUDE_TEMPLATES")
#  if(APPLE)
#    find_library(accelerate_library Accelerate)
#    list(APPEND LIBIGL_COMISO_EXTRA_LIBRARIES "CoMISo" ${accelerate_library})
#  elseif(UNIX)
#    find_package(BLAS REQUIRED)
#    list(APPEND LIBIGL_COMISO_EXTRA_LIBRARIES "CoMISo" ${BLAS_LIBRARIES})
#  endif(APPLE)
#  if(MSVC)
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /D_SCL_SECURE_NO_DEPRECATE")
#    #link_directories("${COMISO_ROOT}/CoMISo/ext/OpenBLAS-v0.2.14-Win64-int64/lib/")
#    list(APPEND LIBIGL_COMISO_EXTRA_LIBRARIES "CoMISo" "${COMISO_DIR}/ext/OpenBLAS-v0.2.14-Win64-int64/lib/libopenblas.dll.a.lib")
#  endif()
#  list(APPEND LIBIGL_EXTRA_LIBRARIES ${LIBIGL_COMISO_EXTRA_LIBRARIES})
#  add_subdirectory("${COMISO_DIR}" "CoMISo")
#  if(MSVC)
#    # Copy the dll
#    add_custom_target(Copy-CoMISo-DLL ALL       # Adds a post-build event to MyTest
#    COMMAND ${CMAKE_COMMAND} -E copy_if_different
#        "${COMISO_DIR}/ext/OpenBLAS-v0.2.14-Win64-int64/bin/libopenblas.dll"
#        "${CMAKE_CURRENT_BINARY_DIR}/../libopenblas.dll"
#    COMMAND ${CMAKE_COMMAND} -E copy_if_different
#        "${COMISO_DIR}/ext/OpenBLAS-v0.2.14-Win64-int64/bin/libgcc_s_seh-1.dll"
#        "${CMAKE_CURRENT_BINARY_DIR}/../libgcc_s_seh-1.dll"
#    COMMAND ${CMAKE_COMMAND} -E copy_if_different
#        "${COMISO_DIR}/ext/OpenBLAS-v0.2.14-Win64-int64/bin/libgfortran-3.dll"
#        "${CMAKE_CURRENT_BINARY_DIR}/../libgfortran-3.dll"
#    COMMAND ${CMAKE_COMMAND} -E copy_if_different
#        "${COMISO_DIR}/ext/OpenBLAS-v0.2.14-Win64-int64/bin/libquadmath-0.dll"
#        "${CMAKE_CURRENT_BINARY_DIR}/../libquadmath-0.dll")
#  endif()
#  if(LIBIGL_USE_STATIC_LIBRARY)
#    CompileIGL_Module_Copyleft("comiso")
#    target_include_directories(iglcomiso PRIVATE ${COMISO_INCLUDE_DIRS})
#    target_compile_definitions(iglcomiso PRIVATE -DINCLUDE_TEMPLATES)
#  endif()
#  
#  include_directories(
#    #${PROJECT_SOURCE_DIR}/libcotrik/src/ 
#    #/usr/include/eigen3
#    #../../libigl/include
#    #../
#    PRIVATE ${COMISO_INCLUDE_DIRS}
#    )
#    
#  foreach(app 
#    QuadGen
#  )
#    add_executable(${app} ${app}.cpp)
#    target_link_libraries(${app} cotrik ${LIBIGL_LIBRARIES} ${LIBIGL_COMISO_EXTRA_LIBRARIES})
#  endforeach()
#endif()


