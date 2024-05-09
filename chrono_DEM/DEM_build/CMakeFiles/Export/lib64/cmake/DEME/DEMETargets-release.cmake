#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "DEME::simulator_multi_gpu" for configuration "Release"
set_property(TARGET DEME::simulator_multi_gpu APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(DEME::simulator_multi_gpu PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CUDA;CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib64/libsimulator_multi_gpu.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS DEME::simulator_multi_gpu )
list(APPEND _IMPORT_CHECK_FILES_FOR_DEME::simulator_multi_gpu "${_IMPORT_PREFIX}/lib64/libsimulator_multi_gpu.a" )

# Import target "DEME::DEMERuntimeDataHelper" for configuration "Release"
set_property(TARGET DEME::DEMERuntimeDataHelper APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(DEME::DEMERuntimeDataHelper PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib64/libDEMERuntimeDataHelper.so"
  IMPORTED_SONAME_RELEASE "libDEMERuntimeDataHelper.so"
  )

list(APPEND _IMPORT_CHECK_TARGETS DEME::DEMERuntimeDataHelper )
list(APPEND _IMPORT_CHECK_FILES_FOR_DEME::DEMERuntimeDataHelper "${_IMPORT_PREFIX}/lib64/libDEMERuntimeDataHelper.so" )

# Import target "DEME::DEMERuntimeDataHelper_install" for configuration "Release"
set_property(TARGET DEME::DEMERuntimeDataHelper_install APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(DEME::DEMERuntimeDataHelper_install PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib64/libDEMERuntimeDataHelper.so"
  IMPORTED_SONAME_RELEASE "libDEMERuntimeDataHelper.so"
  )

list(APPEND _IMPORT_CHECK_TARGETS DEME::DEMERuntimeDataHelper_install )
list(APPEND _IMPORT_CHECK_FILES_FOR_DEME::DEMERuntimeDataHelper_install "${_IMPORT_PREFIX}/lib64/libDEMERuntimeDataHelper.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
