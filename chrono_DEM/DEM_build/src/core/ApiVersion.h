//	Copyright (c) 2021, SBEL GPU Development Team
//	Copyright (c) 2021, University of Wisconsin - Madison
//
//	SPDX-License-Identifier: BSD-3-Clause

#ifndef DEME_APIVER_H
#define DEME_APIVER_H

// Project Version Macros
#define DEME_VERSION_MAJOR 0
#define DEME_VERSION_MINOR 0
#define DEME_VERSION_PATCH 0

// The project version number expressed in the form 0xMMMMmmPP (for easy numerical comparisons)
#define DEME_API_VERSION ((DEME_VERSION_MAJOR << 16) | (DEME_VERSION_MINOR << 8) | (DEME_VERSION_PATCH))


// C++ Standard Macros
#define STD_AUTODETECT (__cplusplus)
#define STD_CXX98 199711L
#define STD_CXX11 201103L
#define STD_CXX14 201402L
#define STD_CXX17 201703L
#define STD_CXX20 202002L

// The C++ Standard Version targeted by the library
#define CXX_TARGET STD_AUTODETECT

// C++ Standard Comparisons
#define CXX_EQUAL(x)    (CXX_TARGET == x)
#define CXX_NEWER(x)    (CXX_TARGET >  x)
#define CXX_OLDER(x)    (CXX_TARGET <  x)

// C++ Standard Composite Comparisons
#define CXX_EQ_NEWER(x)	(CXX_EQUAL(x) || CXX_NEWER(x))
#define CXX_EQ_OLDER(x) (CXX_EQUAL(x) || CXX_OLDER(x))

#define DEME_CUDA_TOOLKIT_HEADERS "/packages/apps/spack/18/opt/spack/gcc-11.2.0/cuda-11.7.0-bhi/include"

#endif
