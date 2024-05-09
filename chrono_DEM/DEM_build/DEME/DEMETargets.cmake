# Generated by CMake

if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" LESS 2.6)
   message(FATAL_ERROR "CMake >= 2.6.0 required")
endif()
cmake_policy(PUSH)
cmake_policy(VERSION 2.6...3.21)
#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Protect against multiple inclusion, which would fail when already imported targets are added once more.
set(_targetsDefined)
set(_targetsNotDefined)
set(_expectedTargets)
foreach(_expectedTarget DEME::simulator_multi_gpu DEME::DEMERuntimeDataHelper DEME::DEMERuntimeDataHelper_install)
  list(APPEND _expectedTargets ${_expectedTarget})
  if(NOT TARGET ${_expectedTarget})
    list(APPEND _targetsNotDefined ${_expectedTarget})
  endif()
  if(TARGET ${_expectedTarget})
    list(APPEND _targetsDefined ${_expectedTarget})
  endif()
endforeach()
if("${_targetsDefined}" STREQUAL "${_expectedTargets}")
  unset(_targetsDefined)
  unset(_targetsNotDefined)
  unset(_expectedTargets)
  set(CMAKE_IMPORT_FILE_VERSION)
  cmake_policy(POP)
  return()
endif()
if(NOT "${_targetsDefined}" STREQUAL "")
  message(FATAL_ERROR "Some (but not all) targets in this export set were already defined.\nTargets Defined: ${_targetsDefined}\nTargets not yet defined: ${_targetsNotDefined}\n")
endif()
unset(_targetsDefined)
unset(_targetsNotDefined)
unset(_expectedTargets)


# Create imported target DEME::simulator_multi_gpu
add_library(DEME::simulator_multi_gpu STATIC IMPORTED)

set_target_properties(DEME::simulator_multi_gpu PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "/home/ppaudya3/chrono_DEM_b/DEM-Engine/src;/home/ppaudya3/chrono_DEM_b/DEM_build/src;/home/ppaudya3/chrono_DEM_b/DEM-Engine/thirdparty/nvidia_helper_math"
  INTERFACE_LINK_LIBRARIES "CUDA::cudart;CUDA::nvrtc;CUDA::cuda_driver;DEME::DEMERuntimeDataHelper"
)

# Create imported target DEME::DEMERuntimeDataHelper
add_library(DEME::DEMERuntimeDataHelper SHARED IMPORTED)

set_target_properties(DEME::DEMERuntimeDataHelper PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "/home/ppaudya3/chrono_DEM_b/DEM-Engine/src"
  INTERFACE_SOURCES "/home/ppaudya3/chrono_DEM_b/DEM-Engine/src/core/utils/RuntimeData.h"
)

# Create imported target DEME::DEMERuntimeDataHelper_install
add_library(DEME::DEMERuntimeDataHelper_install SHARED IMPORTED)

set_target_properties(DEME::DEMERuntimeDataHelper_install PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "/home/ppaudya3/chrono_DEM_b/DEM-Engine/src"
  INTERFACE_LINK_OPTIONS "LINKER:-soname,libDEMERuntimeDataHelper.so"
  INTERFACE_SOURCES "/home/ppaudya3/chrono_DEM_b/DEM-Engine/src/core/utils/RuntimeData.h"
)

# Import target "DEME::simulator_multi_gpu" for configuration "Release"
set_property(TARGET DEME::simulator_multi_gpu APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(DEME::simulator_multi_gpu PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CUDA;CXX"
  IMPORTED_LOCATION_RELEASE "/home/ppaudya3/chrono_DEM_b/DEM_build/libsimulator_multi_gpu.a"
  )

# Import target "DEME::DEMERuntimeDataHelper" for configuration "Release"
set_property(TARGET DEME::DEMERuntimeDataHelper APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(DEME::DEMERuntimeDataHelper PROPERTIES
  IMPORTED_LOCATION_RELEASE "/home/ppaudya3/chrono_DEM_b/DEM_build/src/core/libDEMERuntimeDataHelper.so"
  IMPORTED_SONAME_RELEASE "libDEMERuntimeDataHelper.so"
  )

# Import target "DEME::DEMERuntimeDataHelper_install" for configuration "Release"
set_property(TARGET DEME::DEMERuntimeDataHelper_install APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(DEME::DEMERuntimeDataHelper_install PROPERTIES
  IMPORTED_LOCATION_RELEASE "/home/ppaudya3/chrono_DEM_b/DEM_build/lib_install/libDEMERuntimeDataHelper.so"
  IMPORTED_SONAME_RELEASE "libDEMERuntimeDataHelper.so"
  )

# This file does not depend on other imported targets which have
# been exported from the same project but in a separate export set.

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
cmake_policy(POP)
