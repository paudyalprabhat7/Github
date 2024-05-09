# Install script for directory: /home/ppaudya3/chrono_DEM_b/DEM-Engine/src/DEM

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/ppaudya3/chrono_DEM_b/DEM_install")
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
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/DEM" TYPE FILE FILES
    "/home/ppaudya3/chrono_DEM_b/DEM-Engine/src/DEM/kT.h"
    "/home/ppaudya3/chrono_DEM_b/DEM-Engine/src/DEM/dT.h"
    "/home/ppaudya3/chrono_DEM_b/DEM-Engine/src/DEM/API.h"
    "/home/ppaudya3/chrono_DEM_b/DEM-Engine/src/DEM/Defines.h"
    "/home/ppaudya3/chrono_DEM_b/DEM-Engine/src/DEM/Structs.h"
    "/home/ppaudya3/chrono_DEM_b/DEM-Engine/src/DEM/Models.h"
    "/home/ppaudya3/chrono_DEM_b/DEM-Engine/src/DEM/VariableTypes.h"
    "/home/ppaudya3/chrono_DEM_b/DEM-Engine/src/DEM/BdrsAndObjs.h"
    "/home/ppaudya3/chrono_DEM_b/DEM-Engine/src/DEM/HostSideHelpers.hpp"
    "/home/ppaudya3/chrono_DEM_b/DEM-Engine/src/DEM/utils/Samplers.hpp"
    "/home/ppaudya3/chrono_DEM_b/DEM-Engine/src/DEM/AuxClasses.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/DEM" TYPE DIRECTORY FILES "/home/ppaudya3/chrono_DEM_b/DEM-Engine/src/DEM/utils" FILES_MATCHING REGEX "/utils\\/[^/]*\\.h$" REGEX "/utils\\/[^/]*\\.hpp$")
endif()

