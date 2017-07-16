if (ENABLE_CUDA)
   find_package(CUDA REQUIRED)
   set(CUDA_ARCH_FLAG "-arch=sm_30" CACHE STRING "Use one out of sm_13, sm_20, sm_30, ...")
   set(CUDA_NVCC_FLAGS "-std=c++11" "${CUDA_NVCC_FLAGS}" "${CUDA_ARCH_FLAG}" "-DVIENNACL_WITH_CUDA")
endif(ENABLE_CUDA)

# Export
########

configure_file(${PROJECT_SOURCE_DIR}/ViennaCL-1.7.1/cmake/FindOpenCL.cmake
   "${PROJECT_BINARY_DIR}/FindOpenCL.cmake" COPYONLY)

configure_file(${PROJECT_SOURCE_DIR}/ViennaCL-1.7.1/cmake/ViennaCLConfig.cmake.in
   "${PROJECT_BINARY_DIR}/ViennaCLConfig.cmake" @ONLY)

configure_file(${PROJECT_SOURCE_DIR}/ViennaCL-1.7.1/cmake/ViennaCLConfigVersion.cmake.in
   "${PROJECT_BINARY_DIR}/ViennaCLConfigVersion.cmake" @ONLY)

if (CMAKE_MINOR_VERSION GREATER 6)  # export(PACKAGE ...) introduced with CMake 2.8.0
  export(PACKAGE ViennaCL)
endif()

if (ENABLE_CUDA)
  if (ENABLE_UBLAS)
    include_directories(${PROJECT_SOURCE_DIR} ${Boost_INCLUDE_DIRS}
    ${PROJECT_SOURCE_DIR}/libcotrik/src/
    ${PROJECT_SOURCE_DIR}/ViennaCL-1.7.1
    ${PROJECT_SOURCE_DIR}/ViennaCL-1.7.1/viennacl
    ${PROJECT_SOURCE_DIR}/ViennaCL-1.7.1/libviennacl/include
    ${PROJECT_SOURCE_DIR}/ViennaCL-1.7.1/examples/tutorial
    /usr/include/eigen3)
    foreach(tut iterative)
        cuda_add_executable(${tut}-cotrik-cuda ${tut}.cu)
        target_link_libraries(${tut}-cotrik-cuda ${Boost_LIBRARIES} cotrik)
    endforeach()
  endif ()
endif (ENABLE_CUDA)