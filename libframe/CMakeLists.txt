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
    /usr/include/
    /usr/include/eigen3
)

set(ALL_LIBS
    pthread
    #${VTK_LIBRARIES}
)

add_definitions(
#       -D_CRT_SECURE_NO_WARNINGS
)
add_library(frame SHARED 
    src/component_class.cpp
    src/component_class.h
    src/constants.h
    #src/extract_sheet.cpp
    src/frame_of_mesh.cpp
    src/frame_of_mesh.h
    src/global_functions.cpp
    src/global_functions.h
    src/global_types.cpp
    src/global_types.h
    src/io.cpp
    src/io.h
    src/laplacian_smoothing.cpp
    src/laplacian_smoothing.h
    #src/main.cpp
    src/parameterization.cpp
    src/parameterization.h
    src/parametric_optimization.cpp
    src/parametric_optimization.h
    src/tetgen.h
)
target_link_libraries(frame ${ALL_LIBS})