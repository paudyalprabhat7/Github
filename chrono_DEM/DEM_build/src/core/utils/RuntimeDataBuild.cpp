//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//  SPDX-License-Identifier: BSD-3-Clause

#include <core/utils/RuntimeData.h>

#define RUNTIME_DATA_DIRECTORY "/home/ppaudya3/chrono_DEM_b/DEM_build"
#define RUNTIME_INCLUDE_DIRECTORY "/home/ppaudya3/chrono_DEM_b/DEM_build"

// It's called data_path, but it includes both data files and kernel files
std::filesystem::path RuntimeDataHelper::data_path = std::filesystem::path(RUNTIME_DATA_DIRECTORY);

// This include path is for the included headers in the kernels that need jitification
std::filesystem::path RuntimeDataHelper::include_path = std::filesystem::path(RUNTIME_INCLUDE_DIRECTORY);

