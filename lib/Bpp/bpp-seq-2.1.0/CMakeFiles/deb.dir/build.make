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
CMAKE_SOURCE_DIR = /home/mophylo/MO-Phylogenetics/lib/Bpp/bpp-seq-2.1.0

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/mophylo/MO-Phylogenetics/lib/Bpp/bpp-seq-2.1.0

# Utility rule file for deb.

# Include the progress variables for this target.
include CMakeFiles/deb.dir/progress.make

CMakeFiles/deb:
	dpkg-buildpackage -uc -us -ibpp-seq-2.1.0.tar.gz

deb: CMakeFiles/deb
deb: CMakeFiles/deb.dir/build.make
.PHONY : deb

# Rule to build all files generated by this target.
CMakeFiles/deb.dir/build: deb
.PHONY : CMakeFiles/deb.dir/build

CMakeFiles/deb.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/deb.dir/cmake_clean.cmake
.PHONY : CMakeFiles/deb.dir/clean

CMakeFiles/deb.dir/depend:
	cd /home/mophylo/MO-Phylogenetics/lib/Bpp/bpp-seq-2.1.0 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/mophylo/MO-Phylogenetics/lib/Bpp/bpp-seq-2.1.0 /home/mophylo/MO-Phylogenetics/lib/Bpp/bpp-seq-2.1.0 /home/mophylo/MO-Phylogenetics/lib/Bpp/bpp-seq-2.1.0 /home/mophylo/MO-Phylogenetics/lib/Bpp/bpp-seq-2.1.0 /home/mophylo/MO-Phylogenetics/lib/Bpp/bpp-seq-2.1.0/CMakeFiles/deb.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/deb.dir/depend

