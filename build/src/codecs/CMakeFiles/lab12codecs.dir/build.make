# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

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
CMAKE_SOURCE_DIR = /home/sgfairbro/moos-ivp/moos-ivp-extend-acomms

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/sgfairbro/moos-ivp/moos-ivp-extend-acomms/build

# Include any dependencies generated for this target.
include src/codecs/CMakeFiles/lab12codecs.dir/depend.make

# Include the progress variables for this target.
include src/codecs/CMakeFiles/lab12codecs.dir/progress.make

# Include the compile flags for this target's objects.
include src/codecs/CMakeFiles/lab12codecs.dir/flags.make

src/codecs/CMakeFiles/lab12codecs.dir/varint.cpp.o: src/codecs/CMakeFiles/lab12codecs.dir/flags.make
src/codecs/CMakeFiles/lab12codecs.dir/varint.cpp.o: ../src/codecs/varint.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/sgfairbro/moos-ivp/moos-ivp-extend-acomms/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/codecs/CMakeFiles/lab12codecs.dir/varint.cpp.o"
	cd /home/sgfairbro/moos-ivp/moos-ivp-extend-acomms/build/src/codecs && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/lab12codecs.dir/varint.cpp.o -c /home/sgfairbro/moos-ivp/moos-ivp-extend-acomms/src/codecs/varint.cpp

src/codecs/CMakeFiles/lab12codecs.dir/varint.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lab12codecs.dir/varint.cpp.i"
	cd /home/sgfairbro/moos-ivp/moos-ivp-extend-acomms/build/src/codecs && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/sgfairbro/moos-ivp/moos-ivp-extend-acomms/src/codecs/varint.cpp > CMakeFiles/lab12codecs.dir/varint.cpp.i

src/codecs/CMakeFiles/lab12codecs.dir/varint.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lab12codecs.dir/varint.cpp.s"
	cd /home/sgfairbro/moos-ivp/moos-ivp-extend-acomms/build/src/codecs && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/sgfairbro/moos-ivp/moos-ivp-extend-acomms/src/codecs/varint.cpp -o CMakeFiles/lab12codecs.dir/varint.cpp.s

src/codecs/CMakeFiles/lab12codecs.dir/varint.cpp.o.requires:
.PHONY : src/codecs/CMakeFiles/lab12codecs.dir/varint.cpp.o.requires

src/codecs/CMakeFiles/lab12codecs.dir/varint.cpp.o.provides: src/codecs/CMakeFiles/lab12codecs.dir/varint.cpp.o.requires
	$(MAKE) -f src/codecs/CMakeFiles/lab12codecs.dir/build.make src/codecs/CMakeFiles/lab12codecs.dir/varint.cpp.o.provides.build
.PHONY : src/codecs/CMakeFiles/lab12codecs.dir/varint.cpp.o.provides

src/codecs/CMakeFiles/lab12codecs.dir/varint.cpp.o.provides.build: src/codecs/CMakeFiles/lab12codecs.dir/varint.cpp.o

src/codecs/CMakeFiles/lab12codecs.dir/arithmetic.cpp.o: src/codecs/CMakeFiles/lab12codecs.dir/flags.make
src/codecs/CMakeFiles/lab12codecs.dir/arithmetic.cpp.o: ../src/codecs/arithmetic.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/sgfairbro/moos-ivp/moos-ivp-extend-acomms/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/codecs/CMakeFiles/lab12codecs.dir/arithmetic.cpp.o"
	cd /home/sgfairbro/moos-ivp/moos-ivp-extend-acomms/build/src/codecs && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/lab12codecs.dir/arithmetic.cpp.o -c /home/sgfairbro/moos-ivp/moos-ivp-extend-acomms/src/codecs/arithmetic.cpp

src/codecs/CMakeFiles/lab12codecs.dir/arithmetic.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lab12codecs.dir/arithmetic.cpp.i"
	cd /home/sgfairbro/moos-ivp/moos-ivp-extend-acomms/build/src/codecs && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/sgfairbro/moos-ivp/moos-ivp-extend-acomms/src/codecs/arithmetic.cpp > CMakeFiles/lab12codecs.dir/arithmetic.cpp.i

src/codecs/CMakeFiles/lab12codecs.dir/arithmetic.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lab12codecs.dir/arithmetic.cpp.s"
	cd /home/sgfairbro/moos-ivp/moos-ivp-extend-acomms/build/src/codecs && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/sgfairbro/moos-ivp/moos-ivp-extend-acomms/src/codecs/arithmetic.cpp -o CMakeFiles/lab12codecs.dir/arithmetic.cpp.s

src/codecs/CMakeFiles/lab12codecs.dir/arithmetic.cpp.o.requires:
.PHONY : src/codecs/CMakeFiles/lab12codecs.dir/arithmetic.cpp.o.requires

src/codecs/CMakeFiles/lab12codecs.dir/arithmetic.cpp.o.provides: src/codecs/CMakeFiles/lab12codecs.dir/arithmetic.cpp.o.requires
	$(MAKE) -f src/codecs/CMakeFiles/lab12codecs.dir/build.make src/codecs/CMakeFiles/lab12codecs.dir/arithmetic.cpp.o.provides.build
.PHONY : src/codecs/CMakeFiles/lab12codecs.dir/arithmetic.cpp.o.provides

src/codecs/CMakeFiles/lab12codecs.dir/arithmetic.cpp.o.provides.build: src/codecs/CMakeFiles/lab12codecs.dir/arithmetic.cpp.o

src/codecs/CMakeFiles/lab12codecs.dir/codec_load.cpp.o: src/codecs/CMakeFiles/lab12codecs.dir/flags.make
src/codecs/CMakeFiles/lab12codecs.dir/codec_load.cpp.o: ../src/codecs/codec_load.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/sgfairbro/moos-ivp/moos-ivp-extend-acomms/build/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/codecs/CMakeFiles/lab12codecs.dir/codec_load.cpp.o"
	cd /home/sgfairbro/moos-ivp/moos-ivp-extend-acomms/build/src/codecs && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/lab12codecs.dir/codec_load.cpp.o -c /home/sgfairbro/moos-ivp/moos-ivp-extend-acomms/src/codecs/codec_load.cpp

src/codecs/CMakeFiles/lab12codecs.dir/codec_load.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lab12codecs.dir/codec_load.cpp.i"
	cd /home/sgfairbro/moos-ivp/moos-ivp-extend-acomms/build/src/codecs && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/sgfairbro/moos-ivp/moos-ivp-extend-acomms/src/codecs/codec_load.cpp > CMakeFiles/lab12codecs.dir/codec_load.cpp.i

src/codecs/CMakeFiles/lab12codecs.dir/codec_load.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lab12codecs.dir/codec_load.cpp.s"
	cd /home/sgfairbro/moos-ivp/moos-ivp-extend-acomms/build/src/codecs && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/sgfairbro/moos-ivp/moos-ivp-extend-acomms/src/codecs/codec_load.cpp -o CMakeFiles/lab12codecs.dir/codec_load.cpp.s

src/codecs/CMakeFiles/lab12codecs.dir/codec_load.cpp.o.requires:
.PHONY : src/codecs/CMakeFiles/lab12codecs.dir/codec_load.cpp.o.requires

src/codecs/CMakeFiles/lab12codecs.dir/codec_load.cpp.o.provides: src/codecs/CMakeFiles/lab12codecs.dir/codec_load.cpp.o.requires
	$(MAKE) -f src/codecs/CMakeFiles/lab12codecs.dir/build.make src/codecs/CMakeFiles/lab12codecs.dir/codec_load.cpp.o.provides.build
.PHONY : src/codecs/CMakeFiles/lab12codecs.dir/codec_load.cpp.o.provides

src/codecs/CMakeFiles/lab12codecs.dir/codec_load.cpp.o.provides.build: src/codecs/CMakeFiles/lab12codecs.dir/codec_load.cpp.o

# Object files for target lab12codecs
lab12codecs_OBJECTS = \
"CMakeFiles/lab12codecs.dir/varint.cpp.o" \
"CMakeFiles/lab12codecs.dir/arithmetic.cpp.o" \
"CMakeFiles/lab12codecs.dir/codec_load.cpp.o"

# External object files for target lab12codecs
lab12codecs_EXTERNAL_OBJECTS =

../lib/liblab12codecs.so: src/codecs/CMakeFiles/lab12codecs.dir/varint.cpp.o
../lib/liblab12codecs.so: src/codecs/CMakeFiles/lab12codecs.dir/arithmetic.cpp.o
../lib/liblab12codecs.so: src/codecs/CMakeFiles/lab12codecs.dir/codec_load.cpp.o
../lib/liblab12codecs.so: src/codecs/CMakeFiles/lab12codecs.dir/build.make
../lib/liblab12codecs.so: /usr/lib/x86_64-linux-gnu/libboost_signals.so
../lib/liblab12codecs.so: /usr/lib/x86_64-linux-gnu/libboost_system.so
../lib/liblab12codecs.so: /home/sgfairbro/goby/lib/libgoby_acomms.so
../lib/liblab12codecs.so: /home/sgfairbro/goby/lib/libgoby_moos.so
../lib/liblab12codecs.so: /home/sgfairbro/goby/lib/libgoby_common.so
../lib/liblab12codecs.so: /home/sgfairbro/goby/lib/libgoby_util.so
../lib/liblab12codecs.so: /usr/lib/x86_64-linux-gnu/libprotobuf.so
../lib/liblab12codecs.so: src/codecs/CMakeFiles/lab12codecs.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX shared library ../../../lib/liblab12codecs.so"
	cd /home/sgfairbro/moos-ivp/moos-ivp-extend-acomms/build/src/codecs && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/lab12codecs.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/codecs/CMakeFiles/lab12codecs.dir/build: ../lib/liblab12codecs.so
.PHONY : src/codecs/CMakeFiles/lab12codecs.dir/build

src/codecs/CMakeFiles/lab12codecs.dir/requires: src/codecs/CMakeFiles/lab12codecs.dir/varint.cpp.o.requires
src/codecs/CMakeFiles/lab12codecs.dir/requires: src/codecs/CMakeFiles/lab12codecs.dir/arithmetic.cpp.o.requires
src/codecs/CMakeFiles/lab12codecs.dir/requires: src/codecs/CMakeFiles/lab12codecs.dir/codec_load.cpp.o.requires
.PHONY : src/codecs/CMakeFiles/lab12codecs.dir/requires

src/codecs/CMakeFiles/lab12codecs.dir/clean:
	cd /home/sgfairbro/moos-ivp/moos-ivp-extend-acomms/build/src/codecs && $(CMAKE_COMMAND) -P CMakeFiles/lab12codecs.dir/cmake_clean.cmake
.PHONY : src/codecs/CMakeFiles/lab12codecs.dir/clean

src/codecs/CMakeFiles/lab12codecs.dir/depend:
	cd /home/sgfairbro/moos-ivp/moos-ivp-extend-acomms/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/sgfairbro/moos-ivp/moos-ivp-extend-acomms /home/sgfairbro/moos-ivp/moos-ivp-extend-acomms/src/codecs /home/sgfairbro/moos-ivp/moos-ivp-extend-acomms/build /home/sgfairbro/moos-ivp/moos-ivp-extend-acomms/build/src/codecs /home/sgfairbro/moos-ivp/moos-ivp-extend-acomms/build/src/codecs/CMakeFiles/lab12codecs.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/codecs/CMakeFiles/lab12codecs.dir/depend

