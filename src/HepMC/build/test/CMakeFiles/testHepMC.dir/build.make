# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/build

# Include any dependencies generated for this target.
include test/CMakeFiles/testHepMC.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/testHepMC.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/testHepMC.dir/flags.make

test/CMakeFiles/testHepMC.dir/testHepMC.cc.o: test/CMakeFiles/testHepMC.dir/flags.make
test/CMakeFiles/testHepMC.dir/testHepMC.cc.o: test/testHepMC.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/testHepMC.dir/testHepMC.cc.o"
	cd /scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/build/test && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/testHepMC.dir/testHepMC.cc.o -c /scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/build/test/testHepMC.cc

test/CMakeFiles/testHepMC.dir/testHepMC.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/testHepMC.dir/testHepMC.cc.i"
	cd /scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/build/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/build/test/testHepMC.cc > CMakeFiles/testHepMC.dir/testHepMC.cc.i

test/CMakeFiles/testHepMC.dir/testHepMC.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/testHepMC.dir/testHepMC.cc.s"
	cd /scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/build/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/build/test/testHepMC.cc -o CMakeFiles/testHepMC.dir/testHepMC.cc.s

test/CMakeFiles/testHepMC.dir/testHepMC.cc.o.requires:

.PHONY : test/CMakeFiles/testHepMC.dir/testHepMC.cc.o.requires

test/CMakeFiles/testHepMC.dir/testHepMC.cc.o.provides: test/CMakeFiles/testHepMC.dir/testHepMC.cc.o.requires
	$(MAKE) -f test/CMakeFiles/testHepMC.dir/build.make test/CMakeFiles/testHepMC.dir/testHepMC.cc.o.provides.build
.PHONY : test/CMakeFiles/testHepMC.dir/testHepMC.cc.o.provides

test/CMakeFiles/testHepMC.dir/testHepMC.cc.o.provides.build: test/CMakeFiles/testHepMC.dir/testHepMC.cc.o


test/CMakeFiles/testHepMC.dir/testHepMCMethods.cc.o: test/CMakeFiles/testHepMC.dir/flags.make
test/CMakeFiles/testHepMC.dir/testHepMCMethods.cc.o: /scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/test/testHepMCMethods.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object test/CMakeFiles/testHepMC.dir/testHepMCMethods.cc.o"
	cd /scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/build/test && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/testHepMC.dir/testHepMCMethods.cc.o -c /scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/test/testHepMCMethods.cc

test/CMakeFiles/testHepMC.dir/testHepMCMethods.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/testHepMC.dir/testHepMCMethods.cc.i"
	cd /scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/build/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/test/testHepMCMethods.cc > CMakeFiles/testHepMC.dir/testHepMCMethods.cc.i

test/CMakeFiles/testHepMC.dir/testHepMCMethods.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/testHepMC.dir/testHepMCMethods.cc.s"
	cd /scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/build/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/test/testHepMCMethods.cc -o CMakeFiles/testHepMC.dir/testHepMCMethods.cc.s

test/CMakeFiles/testHepMC.dir/testHepMCMethods.cc.o.requires:

.PHONY : test/CMakeFiles/testHepMC.dir/testHepMCMethods.cc.o.requires

test/CMakeFiles/testHepMC.dir/testHepMCMethods.cc.o.provides: test/CMakeFiles/testHepMC.dir/testHepMCMethods.cc.o.requires
	$(MAKE) -f test/CMakeFiles/testHepMC.dir/build.make test/CMakeFiles/testHepMC.dir/testHepMCMethods.cc.o.provides.build
.PHONY : test/CMakeFiles/testHepMC.dir/testHepMCMethods.cc.o.provides

test/CMakeFiles/testHepMC.dir/testHepMCMethods.cc.o.provides.build: test/CMakeFiles/testHepMC.dir/testHepMCMethods.cc.o


# Object files for target testHepMC
testHepMC_OBJECTS = \
"CMakeFiles/testHepMC.dir/testHepMC.cc.o" \
"CMakeFiles/testHepMC.dir/testHepMCMethods.cc.o"

# External object files for target testHepMC
testHepMC_EXTERNAL_OBJECTS =

test/testHepMC: test/CMakeFiles/testHepMC.dir/testHepMC.cc.o
test/testHepMC: test/CMakeFiles/testHepMC.dir/testHepMCMethods.cc.o
test/testHepMC: test/CMakeFiles/testHepMC.dir/build.make
test/testHepMC: lib/libHepMC.so.4.0.0
test/testHepMC: test/CMakeFiles/testHepMC.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable testHepMC"
	cd /scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/build/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/testHepMC.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/testHepMC.dir/build: test/testHepMC

.PHONY : test/CMakeFiles/testHepMC.dir/build

test/CMakeFiles/testHepMC.dir/requires: test/CMakeFiles/testHepMC.dir/testHepMC.cc.o.requires
test/CMakeFiles/testHepMC.dir/requires: test/CMakeFiles/testHepMC.dir/testHepMCMethods.cc.o.requires

.PHONY : test/CMakeFiles/testHepMC.dir/requires

test/CMakeFiles/testHepMC.dir/clean:
	cd /scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/build/test && $(CMAKE_COMMAND) -P CMakeFiles/testHepMC.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/testHepMC.dir/clean

test/CMakeFiles/testHepMC.dir/depend:
	cd /scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09 /scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/test /scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/build /scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/build/test /scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/build/test/CMakeFiles/testHepMC.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/testHepMC.dir/depend

