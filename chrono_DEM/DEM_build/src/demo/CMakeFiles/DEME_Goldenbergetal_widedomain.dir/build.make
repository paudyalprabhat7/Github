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
include src/demo/CMakeFiles/DEME_Goldenbergetal_widedomain.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/demo/CMakeFiles/DEME_Goldenbergetal_widedomain.dir/compiler_depend.make

# Include the progress variables for this target.
include src/demo/CMakeFiles/DEME_Goldenbergetal_widedomain.dir/progress.make

# Include the compile flags for this target's objects.
include src/demo/CMakeFiles/DEME_Goldenbergetal_widedomain.dir/flags.make

src/demo/CMakeFiles/DEME_Goldenbergetal_widedomain.dir/DEME_Goldenbergetal_widedomain.cpp.o: src/demo/CMakeFiles/DEME_Goldenbergetal_widedomain.dir/flags.make
src/demo/CMakeFiles/DEME_Goldenbergetal_widedomain.dir/DEME_Goldenbergetal_widedomain.cpp.o: /home/ppaudya3/chrono_DEM_b/DEM-Engine/src/demo/DEME_Goldenbergetal_widedomain.cpp
src/demo/CMakeFiles/DEME_Goldenbergetal_widedomain.dir/DEME_Goldenbergetal_widedomain.cpp.o: src/demo/CMakeFiles/DEME_Goldenbergetal_widedomain.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ppaudya3/chrono_DEM_b/DEM_build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/demo/CMakeFiles/DEME_Goldenbergetal_widedomain.dir/DEME_Goldenbergetal_widedomain.cpp.o"
	cd /home/ppaudya3/chrono_DEM_b/DEM_build/src/demo && /packages/apps/spack/18/opt/spack/gcc-11.2.0/gcc-11.2.0-244/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/demo/CMakeFiles/DEME_Goldenbergetal_widedomain.dir/DEME_Goldenbergetal_widedomain.cpp.o -MF CMakeFiles/DEME_Goldenbergetal_widedomain.dir/DEME_Goldenbergetal_widedomain.cpp.o.d -o CMakeFiles/DEME_Goldenbergetal_widedomain.dir/DEME_Goldenbergetal_widedomain.cpp.o -c /home/ppaudya3/chrono_DEM_b/DEM-Engine/src/demo/DEME_Goldenbergetal_widedomain.cpp

src/demo/CMakeFiles/DEME_Goldenbergetal_widedomain.dir/DEME_Goldenbergetal_widedomain.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DEME_Goldenbergetal_widedomain.dir/DEME_Goldenbergetal_widedomain.cpp.i"
	cd /home/ppaudya3/chrono_DEM_b/DEM_build/src/demo && /packages/apps/spack/18/opt/spack/gcc-11.2.0/gcc-11.2.0-244/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ppaudya3/chrono_DEM_b/DEM-Engine/src/demo/DEME_Goldenbergetal_widedomain.cpp > CMakeFiles/DEME_Goldenbergetal_widedomain.dir/DEME_Goldenbergetal_widedomain.cpp.i

src/demo/CMakeFiles/DEME_Goldenbergetal_widedomain.dir/DEME_Goldenbergetal_widedomain.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DEME_Goldenbergetal_widedomain.dir/DEME_Goldenbergetal_widedomain.cpp.s"
	cd /home/ppaudya3/chrono_DEM_b/DEM_build/src/demo && /packages/apps/spack/18/opt/spack/gcc-11.2.0/gcc-11.2.0-244/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ppaudya3/chrono_DEM_b/DEM-Engine/src/demo/DEME_Goldenbergetal_widedomain.cpp -o CMakeFiles/DEME_Goldenbergetal_widedomain.dir/DEME_Goldenbergetal_widedomain.cpp.s

# Object files for target DEME_Goldenbergetal_widedomain
DEME_Goldenbergetal_widedomain_OBJECTS = \
"CMakeFiles/DEME_Goldenbergetal_widedomain.dir/DEME_Goldenbergetal_widedomain.cpp.o"

# External object files for target DEME_Goldenbergetal_widedomain
DEME_Goldenbergetal_widedomain_EXTERNAL_OBJECTS =

bin/DEME_Goldenbergetal_widedomain: src/demo/CMakeFiles/DEME_Goldenbergetal_widedomain.dir/DEME_Goldenbergetal_widedomain.cpp.o
bin/DEME_Goldenbergetal_widedomain: src/demo/CMakeFiles/DEME_Goldenbergetal_widedomain.dir/build.make
bin/DEME_Goldenbergetal_widedomain: libsimulator_multi_gpu.a
bin/DEME_Goldenbergetal_widedomain: /packages/apps/spack/18/opt/spack/gcc-11.2.0/cuda-11.7.0-bhi/lib64/libcudart.so
bin/DEME_Goldenbergetal_widedomain: /packages/apps/spack/18/opt/spack/gcc-11.2.0/cuda-11.7.0-bhi/lib64/libnvrtc.so
bin/DEME_Goldenbergetal_widedomain: /packages/apps/spack/18/opt/spack/gcc-11.2.0/cuda-11.7.0-bhi/lib64/stubs/libcuda.so
bin/DEME_Goldenbergetal_widedomain: src/core/libDEMERuntimeDataHelper.so
bin/DEME_Goldenbergetal_widedomain: src/demo/CMakeFiles/DEME_Goldenbergetal_widedomain.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ppaudya3/chrono_DEM_b/DEM_build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../bin/DEME_Goldenbergetal_widedomain"
	cd /home/ppaudya3/chrono_DEM_b/DEM_build/src/demo && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/DEME_Goldenbergetal_widedomain.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/demo/CMakeFiles/DEME_Goldenbergetal_widedomain.dir/build: bin/DEME_Goldenbergetal_widedomain
.PHONY : src/demo/CMakeFiles/DEME_Goldenbergetal_widedomain.dir/build

src/demo/CMakeFiles/DEME_Goldenbergetal_widedomain.dir/clean:
	cd /home/ppaudya3/chrono_DEM_b/DEM_build/src/demo && $(CMAKE_COMMAND) -P CMakeFiles/DEME_Goldenbergetal_widedomain.dir/cmake_clean.cmake
.PHONY : src/demo/CMakeFiles/DEME_Goldenbergetal_widedomain.dir/clean

src/demo/CMakeFiles/DEME_Goldenbergetal_widedomain.dir/depend:
	cd /home/ppaudya3/chrono_DEM_b/DEM_build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ppaudya3/chrono_DEM_b/DEM-Engine /home/ppaudya3/chrono_DEM_b/DEM-Engine/src/demo /home/ppaudya3/chrono_DEM_b/DEM_build /home/ppaudya3/chrono_DEM_b/DEM_build/src/demo /home/ppaudya3/chrono_DEM_b/DEM_build/src/demo/CMakeFiles/DEME_Goldenbergetal_widedomain.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/demo/CMakeFiles/DEME_Goldenbergetal_widedomain.dir/depend

