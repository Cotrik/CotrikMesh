# CMake entry point
cmake_minimum_required (VERSION 2.6)
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

project (CotrikMesh)

# Set release build:
IF(DEFINED CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}")
ELSE()
  SET (CMAKE_BUILD_TYPE Release)
ENDIF()

add_subdirectory(libcotrik)
add_subdirectory(libnaeem)
#add_subdirectory(libframe)
add_subdirectory(apps)
#add_subdirectory(libigl/external/CoMISo)
#add_subdirectory(cuda-apps)
#add_subdirectory(ViennaCL-1.7.1)

#set (CMAKE_PREFIX_PATH "D:\\Program Files\\CotrikMesh")