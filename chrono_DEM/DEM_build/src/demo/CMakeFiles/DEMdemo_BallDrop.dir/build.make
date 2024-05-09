# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.23

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /packages/apps/spack/18/opt/spack/gcc-11.2.0/cmake-3.23.1-okf/bin/cmake

# The command to remove a file.
RM = /packages/apps/spack/18/opt/spack/gcc-11.2.0/cmake-3.23.1-okf/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/ppaudya3/chrono_DEM_b/DEM-Engine

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ppaudya3/chrono_DEM_b/DEM_build

# Include any dependencies generated for this target.
include src/demo/CMakeFiles/DEMdemo_BallDrop.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/demo/CMakeFiles/DEMdemo_BallDrop.dir/compiler_depend.make

# Include the progress variables for this target.
include src/demo/CMakeFiles/DEMdemo_BallDrop.dir/progress.make

# Include the compile flags for this target's objects.
include src/demo/CMakeFiles/DEMdemo_BallDrop.dir/flags.make

src/demo/CMakeFiles/DEMdemo_BallDrop.dir/DEMdemo_BallDrop.cpp.o: src/demo/CMakeFiles/DEMdemo_BallDrop.dir/flags.make
src/demo/CMakeFiles/DEMdemo_BallDrop.dir/DEMdemo_BallDrop.cpp.o: /home/ppaudya3/chrono_DEM_b/DEM-Engine/src/demo/DEMdemo_BallDrop.cpp
src/demo/CMakeFiles/DEMdemo_BallDrop.dir/DEMdemo_BallDrop.cpp.o: src/demo/CMakeFiles/DEMdemo_BallDrop.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ppaudya3/chrono_DEM_b/DEM_build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/demo/CMakeFiles/DEMdemo_BallDrop.dir/DEMdemo_BallDrop.cpp.o"
	cd /home/ppaudya3/chrono_DEM_b/DEM_build/src/demo && /packages/apps/spack/18/opt/spack/gcc-11.2.0/gcc-11.2.0-244/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/demo/CMakeFiles/DEMdemo_BallDrop.dir/DEMdemo_BallDrop.cpp.o -MF CMakeFiles/DEMdemo_BallDrop.dir/DEMdemo_BallDrop.cpp.o.d -o CMakeFiles/DEMdemo_BallDrop.dir/DEMdemo_BallDrop.cpp.o -c /home/ppaudya3/chrono_DEM_b/DEM-Engine/src/demo/DEMdemo_BallDrop.cpp

src/demo/CMakeFiles/DEMdemo_BallDrop.dir/DEMdemo_BallDrop.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DEMdemo_BallDrop.dir/DEMdemo_BallDrop.cpp.i"
	cd /home/ppaudya3/chrono_DEM_b/DEM_build/src/demo && /packages/apps/spack/18/opt/spack/gcc-11.2.0/gcc-11.2.0-244/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ppaudya3/chrono_DEM_b/DEM-Engine/src/demo/DEMdemo_BallDrop.cpp > CMakeFiles/DEMdemo_BallDrop.dir/DEMdemo_BallDrop.cpp.i

src/demo/CMakeFiles/DEMdemo_BallDrop.dir/DEMdemo_BallDrop.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DEMdemo_BallDrop.dir/DEMdemo_BallDrop.cpp.s"
	cd /home/ppaudya3/chrono_DEM_b/DEM_build/src/demo && /packages/apps/spack/18/opt/spack/gcc-11.2.0/gcc-11.2.0-244/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ppaudya3/chrono_DEM_b/DEM-Engine/src/demo/DEMdemo_BallDrop.cpp -o CMakeFiles/DEMdemo_BallDrop.dir/DEMdemo_BallDrop.cpp.s

# Object files for target DEMdemo_BallDrop
DEMdemo_BallDrop_OBJECTS = \
"CMakeFiles/DEMdemo_BallDrop.dir/DEMdemo_BallDrop.cpp.o"

# External object files for target DEMdemo_BallDrop
DEMdemo_BallDrop_EXTERNAL_OBJECTS =

bin/DEMdemo_BallDrop: src/demo/CMakeFiles/DEMdemo_BallDrop.dir/DEMdemo_BallDrop.cpp.o
bin/DEMdemo_BallDrop: src/demo/CMakeFiles/DEMdemo_BallDrop.dir/build.make
bin/DEMdemo_BallDrop: libsimulator_multi_gpu.a
bin/DEMdemo_BallDrop: /packages/apps/spack/18/opt/spack/gcc-11.2.0/cuda-11.7.0-bhi/lib64/libcudart.so
bin/DEMdemo_BallDrop: /packages/apps/spack/18/opt/spack/gcc-11.2.0/cuda-11.7.0-bhi/lib64/libnvrtc.so
bin/DEMdemo_BallDrop: /packages/apps/spack/18/opt/spack/gcc-11.2.0/cuda-11.7.0-bhi/lib64/stubs/libcuda.so
bin/DEMdemo_BallDrop: src/core/libDEMERuntimeDataHelper.so
bin/DEMdemo_BallDrop: src/demo/CMakeFiles/DEMdemo_BallDrop.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ppaudya3/chrono_DEM_b/DEM_build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../bin/DEMdemo_BallDrop"
	cd /home/ppaudya3/chrono_DEM_b/DEM_build/src/demo && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/DEMdemo_BallDrop.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/demo/CMakeFiles/DEMdemo_BallDrop.dir/build: bin/DEMdemo_BallDrop
.PHONY : src/demo/CMakeFiles/DEMdemo_BallDrop.dir/build

src/demo/CMakeFiles/DEMdemo_BallDrop.dir/clean:
	cd /home/ppaudya3/chrono_DEM_b/DEM_build/src/demo && $(CMAKE_COMMAND) -P CMakeFiles/DEMdemo_BallDrop.dir/cmake_clean.cmake
.PHONY : src/demo/CMakeFiles/DEMdemo_BallDrop.dir/clean

src/demo/CMakeFiles/DEMdemo_BallDrop.dir/depend:
	cd /home/ppaudya3/chrono_DEM_b/DEM_build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ppaudya3/chrono_DEM_b/DEM-Engine /home/ppaudya3/chrono_DEM_b/DEM-Engine/src/demo /home/ppaudya3/chrono_DEM_b/DEM_build /home/ppaudya3/chrono_DEM_b/DEM_build/src/demo /home/ppaudya3/chrono_DEM_b/DEM_build/src/demo/CMakeFiles/DEMdemo_BallDrop.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/demo/CMakeFiles/DEMdemo_BallDrop.dir/depend

