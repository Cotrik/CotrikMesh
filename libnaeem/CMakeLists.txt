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

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})   
include_directories(
    ${PROJECT_SOURCE_DIR}/eigen3
    ${PROJECT_SOURCE_DIR}/verdict/vtkverdict
    ${PROJECT_SOURCE_DIR}/libigl/include
    ${PROJECT_SOURCE_DIR}/libcotrik/src
    ../
    ${VTK_INCLUDE_DIRS}
)

set(ALL_LIBS
    ${VTK_LIBRARIES}
)


add_library(naeem SHARED STATIC 
    src/AngleBasedSmoothQuadMesh.h
    src/AngleBasedSmoothQuadMesh.cpp
    src/ChordCollapse.h
    src/ChordCollapse.cpp
    src/ChordExtractor.h
    src/ChordExtractor.cpp
    src/DiagonalCollapse.h
    src/DiagonalCollapse.cpp
    src/DirectSeparatrixCollapse.h
    src/DirectSeparatrixCollapse.cpp
    src/EdgeCollapse.h
    src/EdgeCollapse.cpp
    src/EdgeRotation.h
    src/EdgeRotation.cpp
    src/FeatureExtractor.h
    src/FeatureExtractor.cpp
    src/LocalSimplifier.h
    src/LocalSimplifier.cpp
    src/MeshUtil.h
    src/MeshUtil.cpp
    src/ParallelFor.h
    src/PQueue.h
    src/PQueue.cpp
    src/QuadSurfaceMapper.h
    src/QuadSurfaceMapper.cpp
    src/SemiGlobalSimplifier.h
    src/SemiGlobalSimplifier.cpp
    src/SeparatrixCollapse.h
    src/SeparatrixCollapse.cpp
    src/SimplificationOperation.h
    src/SimplificationOperation.cpp
    src/Smooth.h
    src/Smooth.cpp
    src/VertexRotation.h
    src/VertexRotation.cpp
    src/VertexSplit.h
    src/VertexSplit.cpp
)
target_link_libraries(naeem cotrik ${ALL_LIBS})
