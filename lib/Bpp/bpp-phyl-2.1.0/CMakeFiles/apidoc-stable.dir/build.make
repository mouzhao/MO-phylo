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

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/mophylo/MO-Phylogenetics/lib/Bpp/bpp-phyl-2.1.0

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/mophylo/MO-Phylogenetics/lib/Bpp/bpp-phyl-2.1.0

# Utility rule file for apidoc-stable.

# Include the progress variables for this target.
include CMakeFiles/apidoc-stable.dir/progress.make

CMakeFiles/apidoc-stable:
	cp Doxyfile /home/mophylo/MO-Phylogenetics/lib/Bpp/bpp-phyl-2.1.0/Doxyfile-stable
	echo OUTPUT_DIRECTORY=/home/mophylo/MO-Phylogenetics/lib/Bpp/bpp-phyl-2.1.0 >> /home/mophylo/MO-Phylogenetics/lib/Bpp/bpp-phyl-2.1.0/Doxyfile-stable
	echo HTML_HEADER=header.html >> /home/mophylo/MO-Phylogenetics/lib/Bpp/bpp-phyl-2.1.0/Doxyfile-stable
	/usr/bin/doxygen /home/mophylo/MO-Phylogenetics/lib/Bpp/bpp-phyl-2.1.0/Doxyfile-stable

apidoc-stable: CMakeFiles/apidoc-stable
apidoc-stable: CMakeFiles/apidoc-stable.dir/build.make
.PHONY : apidoc-stable

# Rule to build all files generated by this target.
CMakeFiles/apidoc-stable.dir/build: apidoc-stable
.PHONY : CMakeFiles/apidoc-stable.dir/build

CMakeFiles/apidoc-stable.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/apidoc-stable.dir/cmake_clean.cmake
.PHONY : CMakeFiles/apidoc-stable.dir/clean

CMakeFiles/apidoc-stable.dir/depend:
	cd /home/mophylo/MO-Phylogenetics/lib/Bpp/bpp-phyl-2.1.0 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/mophylo/MO-Phylogenetics/lib/Bpp/bpp-phyl-2.1.0 /home/mophylo/MO-Phylogenetics/lib/Bpp/bpp-phyl-2.1.0 /home/mophylo/MO-Phylogenetics/lib/Bpp/bpp-phyl-2.1.0 /home/mophylo/MO-Phylogenetics/lib/Bpp/bpp-phyl-2.1.0 /home/mophylo/MO-Phylogenetics/lib/Bpp/bpp-phyl-2.1.0/CMakeFiles/apidoc-stable.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/apidoc-stable.dir/depend

