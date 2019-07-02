# Install script for directory: /home/cotrik/svn/CotrikMesh/libigl/external/cgal/src/CGAL_Project

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/cotrik/svn/CotrikMesh/libigl/external/cgal")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/doc/CGAL-4.11" TYPE FILE FILES
    "/home/cotrik/svn/CotrikMesh/libigl/external/cgal/src/CGAL_Project/AUTHORS"
    "/home/cotrik/svn/CotrikMesh/libigl/external/cgal/src/CGAL_Project/CHANGES"
    "/home/cotrik/svn/CotrikMesh/libigl/external/cgal/src/CGAL_Project/LICENSE"
    "/home/cotrik/svn/CotrikMesh/libigl/external/cgal/src/CGAL_Project/LICENSE.FREE_USE"
    "/home/cotrik/svn/CotrikMesh/libigl/external/cgal/src/CGAL_Project/LICENSE.GPL"
    "/home/cotrik/svn/CotrikMesh/libigl/external/cgal/src/CGAL_Project/LICENSE.LGPL"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE DIRECTORY FILES "/home/cotrik/svn/CotrikMesh/libigl/external/cgal/src/CGAL_Project/include/CGAL" REGEX "/\\.svn$" EXCLUDE)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE DIRECTORY FILES "/home/cotrik/svn/CotrikMesh/libigl/external/cgal/src/CGAL_Project-build/include/CGAL" REGEX "/\\.svn$" EXCLUDE)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE PROGRAM FILES
    "/home/cotrik/svn/CotrikMesh/libigl/external/cgal/src/CGAL_Project/scripts/cgal_create_CMakeLists"
    "/home/cotrik/svn/CotrikMesh/libigl/external/cgal/src/CGAL_Project/scripts/cgal_create_cmake_script"
    "/home/cotrik/svn/CotrikMesh/libigl/external/cgal/src/CGAL_Project/scripts/cgal_make_macosx_app"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/CGAL" TYPE DIRECTORY FILES "/home/cotrik/svn/CotrikMesh/libigl/external/cgal/src/CGAL_Project/cmake/modules/")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/CGAL" TYPE FILE FILES "/home/cotrik/svn/CotrikMesh/libigl/external/cgal/src/CGAL_Project/cmake/modules/UseCGAL.cmake")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/CGAL" TYPE FILE FILES
    "/home/cotrik/svn/CotrikMesh/libigl/external/cgal/src/CGAL_Project-build/CGALConfigVersion.cmake"
    "/home/cotrik/svn/CotrikMesh/libigl/external/cgal/src/CGAL_Project-build/config/CGALConfig.cmake"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/man/man1" TYPE FILE FILES "/home/cotrik/svn/CotrikMesh/libigl/external/cgal/src/CGAL_Project/auxiliary/cgal_create_cmake_script.1")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/cotrik/svn/CotrikMesh/libigl/external/cgal/src/CGAL_Project-build/src/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/cotrik/svn/CotrikMesh/libigl/external/cgal/src/CGAL_Project-build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
