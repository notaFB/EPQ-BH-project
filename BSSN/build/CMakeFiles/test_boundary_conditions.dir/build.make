# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_SOURCE_DIR = /home/andres/Desktop/not/BSSN01/bssn_simulation

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/andres/Desktop/not/BSSN01/bssn_simulation/build

# Include any dependencies generated for this target.
include CMakeFiles/test_boundary_conditions.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/test_boundary_conditions.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/test_boundary_conditions.dir/flags.make

CMakeFiles/test_boundary_conditions.dir/test/test_boundary_conditions.cpp.o: CMakeFiles/test_boundary_conditions.dir/flags.make
CMakeFiles/test_boundary_conditions.dir/test/test_boundary_conditions.cpp.o: ../test/test_boundary_conditions.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/andres/Desktop/not/BSSN01/bssn_simulation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/test_boundary_conditions.dir/test/test_boundary_conditions.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/test_boundary_conditions.dir/test/test_boundary_conditions.cpp.o -c /home/andres/Desktop/not/BSSN01/bssn_simulation/test/test_boundary_conditions.cpp

CMakeFiles/test_boundary_conditions.dir/test/test_boundary_conditions.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_boundary_conditions.dir/test/test_boundary_conditions.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/andres/Desktop/not/BSSN01/bssn_simulation/test/test_boundary_conditions.cpp > CMakeFiles/test_boundary_conditions.dir/test/test_boundary_conditions.cpp.i

CMakeFiles/test_boundary_conditions.dir/test/test_boundary_conditions.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_boundary_conditions.dir/test/test_boundary_conditions.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/andres/Desktop/not/BSSN01/bssn_simulation/test/test_boundary_conditions.cpp -o CMakeFiles/test_boundary_conditions.dir/test/test_boundary_conditions.cpp.s

CMakeFiles/test_boundary_conditions.dir/src/boundary_conditions.cpp.o: CMakeFiles/test_boundary_conditions.dir/flags.make
CMakeFiles/test_boundary_conditions.dir/src/boundary_conditions.cpp.o: ../src/boundary_conditions.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/andres/Desktop/not/BSSN01/bssn_simulation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/test_boundary_conditions.dir/src/boundary_conditions.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/test_boundary_conditions.dir/src/boundary_conditions.cpp.o -c /home/andres/Desktop/not/BSSN01/bssn_simulation/src/boundary_conditions.cpp

CMakeFiles/test_boundary_conditions.dir/src/boundary_conditions.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_boundary_conditions.dir/src/boundary_conditions.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/andres/Desktop/not/BSSN01/bssn_simulation/src/boundary_conditions.cpp > CMakeFiles/test_boundary_conditions.dir/src/boundary_conditions.cpp.i

CMakeFiles/test_boundary_conditions.dir/src/boundary_conditions.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_boundary_conditions.dir/src/boundary_conditions.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/andres/Desktop/not/BSSN01/bssn_simulation/src/boundary_conditions.cpp -o CMakeFiles/test_boundary_conditions.dir/src/boundary_conditions.cpp.s

CMakeFiles/test_boundary_conditions.dir/src/config.cpp.o: CMakeFiles/test_boundary_conditions.dir/flags.make
CMakeFiles/test_boundary_conditions.dir/src/config.cpp.o: ../src/config.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/andres/Desktop/not/BSSN01/bssn_simulation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/test_boundary_conditions.dir/src/config.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/test_boundary_conditions.dir/src/config.cpp.o -c /home/andres/Desktop/not/BSSN01/bssn_simulation/src/config.cpp

CMakeFiles/test_boundary_conditions.dir/src/config.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_boundary_conditions.dir/src/config.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/andres/Desktop/not/BSSN01/bssn_simulation/src/config.cpp > CMakeFiles/test_boundary_conditions.dir/src/config.cpp.i

CMakeFiles/test_boundary_conditions.dir/src/config.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_boundary_conditions.dir/src/config.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/andres/Desktop/not/BSSN01/bssn_simulation/src/config.cpp -o CMakeFiles/test_boundary_conditions.dir/src/config.cpp.s

CMakeFiles/test_boundary_conditions.dir/src/evolution.cpp.o: CMakeFiles/test_boundary_conditions.dir/flags.make
CMakeFiles/test_boundary_conditions.dir/src/evolution.cpp.o: ../src/evolution.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/andres/Desktop/not/BSSN01/bssn_simulation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/test_boundary_conditions.dir/src/evolution.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/test_boundary_conditions.dir/src/evolution.cpp.o -c /home/andres/Desktop/not/BSSN01/bssn_simulation/src/evolution.cpp

CMakeFiles/test_boundary_conditions.dir/src/evolution.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_boundary_conditions.dir/src/evolution.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/andres/Desktop/not/BSSN01/bssn_simulation/src/evolution.cpp > CMakeFiles/test_boundary_conditions.dir/src/evolution.cpp.i

CMakeFiles/test_boundary_conditions.dir/src/evolution.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_boundary_conditions.dir/src/evolution.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/andres/Desktop/not/BSSN01/bssn_simulation/src/evolution.cpp -o CMakeFiles/test_boundary_conditions.dir/src/evolution.cpp.s

CMakeFiles/test_boundary_conditions.dir/src/fieldData.cpp.o: CMakeFiles/test_boundary_conditions.dir/flags.make
CMakeFiles/test_boundary_conditions.dir/src/fieldData.cpp.o: ../src/fieldData.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/andres/Desktop/not/BSSN01/bssn_simulation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/test_boundary_conditions.dir/src/fieldData.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/test_boundary_conditions.dir/src/fieldData.cpp.o -c /home/andres/Desktop/not/BSSN01/bssn_simulation/src/fieldData.cpp

CMakeFiles/test_boundary_conditions.dir/src/fieldData.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_boundary_conditions.dir/src/fieldData.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/andres/Desktop/not/BSSN01/bssn_simulation/src/fieldData.cpp > CMakeFiles/test_boundary_conditions.dir/src/fieldData.cpp.i

CMakeFiles/test_boundary_conditions.dir/src/fieldData.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_boundary_conditions.dir/src/fieldData.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/andres/Desktop/not/BSSN01/bssn_simulation/src/fieldData.cpp -o CMakeFiles/test_boundary_conditions.dir/src/fieldData.cpp.s

CMakeFiles/test_boundary_conditions.dir/src/grid.cpp.o: CMakeFiles/test_boundary_conditions.dir/flags.make
CMakeFiles/test_boundary_conditions.dir/src/grid.cpp.o: ../src/grid.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/andres/Desktop/not/BSSN01/bssn_simulation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/test_boundary_conditions.dir/src/grid.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/test_boundary_conditions.dir/src/grid.cpp.o -c /home/andres/Desktop/not/BSSN01/bssn_simulation/src/grid.cpp

CMakeFiles/test_boundary_conditions.dir/src/grid.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_boundary_conditions.dir/src/grid.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/andres/Desktop/not/BSSN01/bssn_simulation/src/grid.cpp > CMakeFiles/test_boundary_conditions.dir/src/grid.cpp.i

CMakeFiles/test_boundary_conditions.dir/src/grid.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_boundary_conditions.dir/src/grid.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/andres/Desktop/not/BSSN01/bssn_simulation/src/grid.cpp -o CMakeFiles/test_boundary_conditions.dir/src/grid.cpp.s

CMakeFiles/test_boundary_conditions.dir/src/initial_conditions.cpp.o: CMakeFiles/test_boundary_conditions.dir/flags.make
CMakeFiles/test_boundary_conditions.dir/src/initial_conditions.cpp.o: ../src/initial_conditions.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/andres/Desktop/not/BSSN01/bssn_simulation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/test_boundary_conditions.dir/src/initial_conditions.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/test_boundary_conditions.dir/src/initial_conditions.cpp.o -c /home/andres/Desktop/not/BSSN01/bssn_simulation/src/initial_conditions.cpp

CMakeFiles/test_boundary_conditions.dir/src/initial_conditions.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_boundary_conditions.dir/src/initial_conditions.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/andres/Desktop/not/BSSN01/bssn_simulation/src/initial_conditions.cpp > CMakeFiles/test_boundary_conditions.dir/src/initial_conditions.cpp.i

CMakeFiles/test_boundary_conditions.dir/src/initial_conditions.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_boundary_conditions.dir/src/initial_conditions.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/andres/Desktop/not/BSSN01/bssn_simulation/src/initial_conditions.cpp -o CMakeFiles/test_boundary_conditions.dir/src/initial_conditions.cpp.s

CMakeFiles/test_boundary_conditions.dir/src/main.cpp.o: CMakeFiles/test_boundary_conditions.dir/flags.make
CMakeFiles/test_boundary_conditions.dir/src/main.cpp.o: ../src/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/andres/Desktop/not/BSSN01/bssn_simulation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/test_boundary_conditions.dir/src/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/test_boundary_conditions.dir/src/main.cpp.o -c /home/andres/Desktop/not/BSSN01/bssn_simulation/src/main.cpp

CMakeFiles/test_boundary_conditions.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_boundary_conditions.dir/src/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/andres/Desktop/not/BSSN01/bssn_simulation/src/main.cpp > CMakeFiles/test_boundary_conditions.dir/src/main.cpp.i

CMakeFiles/test_boundary_conditions.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_boundary_conditions.dir/src/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/andres/Desktop/not/BSSN01/bssn_simulation/src/main.cpp -o CMakeFiles/test_boundary_conditions.dir/src/main.cpp.s

CMakeFiles/test_boundary_conditions.dir/src/utilities.cpp.o: CMakeFiles/test_boundary_conditions.dir/flags.make
CMakeFiles/test_boundary_conditions.dir/src/utilities.cpp.o: ../src/utilities.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/andres/Desktop/not/BSSN01/bssn_simulation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object CMakeFiles/test_boundary_conditions.dir/src/utilities.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/test_boundary_conditions.dir/src/utilities.cpp.o -c /home/andres/Desktop/not/BSSN01/bssn_simulation/src/utilities.cpp

CMakeFiles/test_boundary_conditions.dir/src/utilities.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_boundary_conditions.dir/src/utilities.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/andres/Desktop/not/BSSN01/bssn_simulation/src/utilities.cpp > CMakeFiles/test_boundary_conditions.dir/src/utilities.cpp.i

CMakeFiles/test_boundary_conditions.dir/src/utilities.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_boundary_conditions.dir/src/utilities.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/andres/Desktop/not/BSSN01/bssn_simulation/src/utilities.cpp -o CMakeFiles/test_boundary_conditions.dir/src/utilities.cpp.s

# Object files for target test_boundary_conditions
test_boundary_conditions_OBJECTS = \
"CMakeFiles/test_boundary_conditions.dir/test/test_boundary_conditions.cpp.o" \
"CMakeFiles/test_boundary_conditions.dir/src/boundary_conditions.cpp.o" \
"CMakeFiles/test_boundary_conditions.dir/src/config.cpp.o" \
"CMakeFiles/test_boundary_conditions.dir/src/evolution.cpp.o" \
"CMakeFiles/test_boundary_conditions.dir/src/fieldData.cpp.o" \
"CMakeFiles/test_boundary_conditions.dir/src/grid.cpp.o" \
"CMakeFiles/test_boundary_conditions.dir/src/initial_conditions.cpp.o" \
"CMakeFiles/test_boundary_conditions.dir/src/main.cpp.o" \
"CMakeFiles/test_boundary_conditions.dir/src/utilities.cpp.o"

# External object files for target test_boundary_conditions
test_boundary_conditions_EXTERNAL_OBJECTS =

bin/test_boundary_conditions: CMakeFiles/test_boundary_conditions.dir/test/test_boundary_conditions.cpp.o
bin/test_boundary_conditions: CMakeFiles/test_boundary_conditions.dir/src/boundary_conditions.cpp.o
bin/test_boundary_conditions: CMakeFiles/test_boundary_conditions.dir/src/config.cpp.o
bin/test_boundary_conditions: CMakeFiles/test_boundary_conditions.dir/src/evolution.cpp.o
bin/test_boundary_conditions: CMakeFiles/test_boundary_conditions.dir/src/fieldData.cpp.o
bin/test_boundary_conditions: CMakeFiles/test_boundary_conditions.dir/src/grid.cpp.o
bin/test_boundary_conditions: CMakeFiles/test_boundary_conditions.dir/src/initial_conditions.cpp.o
bin/test_boundary_conditions: CMakeFiles/test_boundary_conditions.dir/src/main.cpp.o
bin/test_boundary_conditions: CMakeFiles/test_boundary_conditions.dir/src/utilities.cpp.o
bin/test_boundary_conditions: CMakeFiles/test_boundary_conditions.dir/build.make
bin/test_boundary_conditions: CMakeFiles/test_boundary_conditions.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/andres/Desktop/not/BSSN01/bssn_simulation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Linking CXX executable bin/test_boundary_conditions"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_boundary_conditions.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/test_boundary_conditions.dir/build: bin/test_boundary_conditions

.PHONY : CMakeFiles/test_boundary_conditions.dir/build

CMakeFiles/test_boundary_conditions.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/test_boundary_conditions.dir/cmake_clean.cmake
.PHONY : CMakeFiles/test_boundary_conditions.dir/clean

CMakeFiles/test_boundary_conditions.dir/depend:
	cd /home/andres/Desktop/not/BSSN01/bssn_simulation/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/andres/Desktop/not/BSSN01/bssn_simulation /home/andres/Desktop/not/BSSN01/bssn_simulation /home/andres/Desktop/not/BSSN01/bssn_simulation/build /home/andres/Desktop/not/BSSN01/bssn_simulation/build /home/andres/Desktop/not/BSSN01/bssn_simulation/build/CMakeFiles/test_boundary_conditions.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/test_boundary_conditions.dir/depend

