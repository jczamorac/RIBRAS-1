# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.13

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
CMAKE_SOURCE_DIR = /home/leo/Desktop/RIBRAS

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/leo/Desktop/RIBRAS/build

# Include any dependencies generated for this target.
include CMakeFiles/arquivobinario.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/arquivobinario.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/arquivobinario.dir/flags.make

CMakeFiles/arquivobinario.dir/Example.cc.o: CMakeFiles/arquivobinario.dir/flags.make
CMakeFiles/arquivobinario.dir/Example.cc.o: ../Example.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/leo/Desktop/RIBRAS/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/arquivobinario.dir/Example.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/arquivobinario.dir/Example.cc.o -c /home/leo/Desktop/RIBRAS/Example.cc

CMakeFiles/arquivobinario.dir/Example.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/arquivobinario.dir/Example.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/leo/Desktop/RIBRAS/Example.cc > CMakeFiles/arquivobinario.dir/Example.cc.i

CMakeFiles/arquivobinario.dir/Example.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/arquivobinario.dir/Example.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/leo/Desktop/RIBRAS/Example.cc -o CMakeFiles/arquivobinario.dir/Example.cc.s

CMakeFiles/arquivobinario.dir/src/BeamTestPrimaryGeneratorAction.cc.o: CMakeFiles/arquivobinario.dir/flags.make
CMakeFiles/arquivobinario.dir/src/BeamTestPrimaryGeneratorAction.cc.o: ../src/BeamTestPrimaryGeneratorAction.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/leo/Desktop/RIBRAS/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/arquivobinario.dir/src/BeamTestPrimaryGeneratorAction.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/arquivobinario.dir/src/BeamTestPrimaryGeneratorAction.cc.o -c /home/leo/Desktop/RIBRAS/src/BeamTestPrimaryGeneratorAction.cc

CMakeFiles/arquivobinario.dir/src/BeamTestPrimaryGeneratorAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/arquivobinario.dir/src/BeamTestPrimaryGeneratorAction.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/leo/Desktop/RIBRAS/src/BeamTestPrimaryGeneratorAction.cc > CMakeFiles/arquivobinario.dir/src/BeamTestPrimaryGeneratorAction.cc.i

CMakeFiles/arquivobinario.dir/src/BeamTestPrimaryGeneratorAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/arquivobinario.dir/src/BeamTestPrimaryGeneratorAction.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/leo/Desktop/RIBRAS/src/BeamTestPrimaryGeneratorAction.cc -o CMakeFiles/arquivobinario.dir/src/BeamTestPrimaryGeneratorAction.cc.s

CMakeFiles/arquivobinario.dir/src/BeamTestPrimaryGeneratorMessenger.cc.o: CMakeFiles/arquivobinario.dir/flags.make
CMakeFiles/arquivobinario.dir/src/BeamTestPrimaryGeneratorMessenger.cc.o: ../src/BeamTestPrimaryGeneratorMessenger.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/leo/Desktop/RIBRAS/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/arquivobinario.dir/src/BeamTestPrimaryGeneratorMessenger.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/arquivobinario.dir/src/BeamTestPrimaryGeneratorMessenger.cc.o -c /home/leo/Desktop/RIBRAS/src/BeamTestPrimaryGeneratorMessenger.cc

CMakeFiles/arquivobinario.dir/src/BeamTestPrimaryGeneratorMessenger.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/arquivobinario.dir/src/BeamTestPrimaryGeneratorMessenger.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/leo/Desktop/RIBRAS/src/BeamTestPrimaryGeneratorMessenger.cc > CMakeFiles/arquivobinario.dir/src/BeamTestPrimaryGeneratorMessenger.cc.i

CMakeFiles/arquivobinario.dir/src/BeamTestPrimaryGeneratorMessenger.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/arquivobinario.dir/src/BeamTestPrimaryGeneratorMessenger.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/leo/Desktop/RIBRAS/src/BeamTestPrimaryGeneratorMessenger.cc -o CMakeFiles/arquivobinario.dir/src/BeamTestPrimaryGeneratorMessenger.cc.s

CMakeFiles/arquivobinario.dir/src/DetectorConstruction.cc.o: CMakeFiles/arquivobinario.dir/flags.make
CMakeFiles/arquivobinario.dir/src/DetectorConstruction.cc.o: ../src/DetectorConstruction.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/leo/Desktop/RIBRAS/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/arquivobinario.dir/src/DetectorConstruction.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/arquivobinario.dir/src/DetectorConstruction.cc.o -c /home/leo/Desktop/RIBRAS/src/DetectorConstruction.cc

CMakeFiles/arquivobinario.dir/src/DetectorConstruction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/arquivobinario.dir/src/DetectorConstruction.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/leo/Desktop/RIBRAS/src/DetectorConstruction.cc > CMakeFiles/arquivobinario.dir/src/DetectorConstruction.cc.i

CMakeFiles/arquivobinario.dir/src/DetectorConstruction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/arquivobinario.dir/src/DetectorConstruction.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/leo/Desktop/RIBRAS/src/DetectorConstruction.cc -o CMakeFiles/arquivobinario.dir/src/DetectorConstruction.cc.s

CMakeFiles/arquivobinario.dir/src/DetectorMessenger.cc.o: CMakeFiles/arquivobinario.dir/flags.make
CMakeFiles/arquivobinario.dir/src/DetectorMessenger.cc.o: ../src/DetectorMessenger.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/leo/Desktop/RIBRAS/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/arquivobinario.dir/src/DetectorMessenger.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/arquivobinario.dir/src/DetectorMessenger.cc.o -c /home/leo/Desktop/RIBRAS/src/DetectorMessenger.cc

CMakeFiles/arquivobinario.dir/src/DetectorMessenger.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/arquivobinario.dir/src/DetectorMessenger.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/leo/Desktop/RIBRAS/src/DetectorMessenger.cc > CMakeFiles/arquivobinario.dir/src/DetectorMessenger.cc.i

CMakeFiles/arquivobinario.dir/src/DetectorMessenger.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/arquivobinario.dir/src/DetectorMessenger.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/leo/Desktop/RIBRAS/src/DetectorMessenger.cc -o CMakeFiles/arquivobinario.dir/src/DetectorMessenger.cc.s

CMakeFiles/arquivobinario.dir/src/EventAction.cc.o: CMakeFiles/arquivobinario.dir/flags.make
CMakeFiles/arquivobinario.dir/src/EventAction.cc.o: ../src/EventAction.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/leo/Desktop/RIBRAS/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/arquivobinario.dir/src/EventAction.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/arquivobinario.dir/src/EventAction.cc.o -c /home/leo/Desktop/RIBRAS/src/EventAction.cc

CMakeFiles/arquivobinario.dir/src/EventAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/arquivobinario.dir/src/EventAction.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/leo/Desktop/RIBRAS/src/EventAction.cc > CMakeFiles/arquivobinario.dir/src/EventAction.cc.i

CMakeFiles/arquivobinario.dir/src/EventAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/arquivobinario.dir/src/EventAction.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/leo/Desktop/RIBRAS/src/EventAction.cc -o CMakeFiles/arquivobinario.dir/src/EventAction.cc.s

CMakeFiles/arquivobinario.dir/src/MagneticField.cc.o: CMakeFiles/arquivobinario.dir/flags.make
CMakeFiles/arquivobinario.dir/src/MagneticField.cc.o: ../src/MagneticField.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/leo/Desktop/RIBRAS/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/arquivobinario.dir/src/MagneticField.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/arquivobinario.dir/src/MagneticField.cc.o -c /home/leo/Desktop/RIBRAS/src/MagneticField.cc

CMakeFiles/arquivobinario.dir/src/MagneticField.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/arquivobinario.dir/src/MagneticField.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/leo/Desktop/RIBRAS/src/MagneticField.cc > CMakeFiles/arquivobinario.dir/src/MagneticField.cc.i

CMakeFiles/arquivobinario.dir/src/MagneticField.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/arquivobinario.dir/src/MagneticField.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/leo/Desktop/RIBRAS/src/MagneticField.cc -o CMakeFiles/arquivobinario.dir/src/MagneticField.cc.s

CMakeFiles/arquivobinario.dir/src/MagneticFieldMessenger.cc.o: CMakeFiles/arquivobinario.dir/flags.make
CMakeFiles/arquivobinario.dir/src/MagneticFieldMessenger.cc.o: ../src/MagneticFieldMessenger.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/leo/Desktop/RIBRAS/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/arquivobinario.dir/src/MagneticFieldMessenger.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/arquivobinario.dir/src/MagneticFieldMessenger.cc.o -c /home/leo/Desktop/RIBRAS/src/MagneticFieldMessenger.cc

CMakeFiles/arquivobinario.dir/src/MagneticFieldMessenger.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/arquivobinario.dir/src/MagneticFieldMessenger.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/leo/Desktop/RIBRAS/src/MagneticFieldMessenger.cc > CMakeFiles/arquivobinario.dir/src/MagneticFieldMessenger.cc.i

CMakeFiles/arquivobinario.dir/src/MagneticFieldMessenger.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/arquivobinario.dir/src/MagneticFieldMessenger.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/leo/Desktop/RIBRAS/src/MagneticFieldMessenger.cc -o CMakeFiles/arquivobinario.dir/src/MagneticFieldMessenger.cc.s

CMakeFiles/arquivobinario.dir/src/PhysicsList.cc.o: CMakeFiles/arquivobinario.dir/flags.make
CMakeFiles/arquivobinario.dir/src/PhysicsList.cc.o: ../src/PhysicsList.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/leo/Desktop/RIBRAS/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object CMakeFiles/arquivobinario.dir/src/PhysicsList.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/arquivobinario.dir/src/PhysicsList.cc.o -c /home/leo/Desktop/RIBRAS/src/PhysicsList.cc

CMakeFiles/arquivobinario.dir/src/PhysicsList.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/arquivobinario.dir/src/PhysicsList.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/leo/Desktop/RIBRAS/src/PhysicsList.cc > CMakeFiles/arquivobinario.dir/src/PhysicsList.cc.i

CMakeFiles/arquivobinario.dir/src/PhysicsList.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/arquivobinario.dir/src/PhysicsList.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/leo/Desktop/RIBRAS/src/PhysicsList.cc -o CMakeFiles/arquivobinario.dir/src/PhysicsList.cc.s

CMakeFiles/arquivobinario.dir/src/PhysicsListMessenger.cc.o: CMakeFiles/arquivobinario.dir/flags.make
CMakeFiles/arquivobinario.dir/src/PhysicsListMessenger.cc.o: ../src/PhysicsListMessenger.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/leo/Desktop/RIBRAS/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object CMakeFiles/arquivobinario.dir/src/PhysicsListMessenger.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/arquivobinario.dir/src/PhysicsListMessenger.cc.o -c /home/leo/Desktop/RIBRAS/src/PhysicsListMessenger.cc

CMakeFiles/arquivobinario.dir/src/PhysicsListMessenger.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/arquivobinario.dir/src/PhysicsListMessenger.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/leo/Desktop/RIBRAS/src/PhysicsListMessenger.cc > CMakeFiles/arquivobinario.dir/src/PhysicsListMessenger.cc.i

CMakeFiles/arquivobinario.dir/src/PhysicsListMessenger.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/arquivobinario.dir/src/PhysicsListMessenger.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/leo/Desktop/RIBRAS/src/PhysicsListMessenger.cc -o CMakeFiles/arquivobinario.dir/src/PhysicsListMessenger.cc.s

CMakeFiles/arquivobinario.dir/src/Reaction.cc.o: CMakeFiles/arquivobinario.dir/flags.make
CMakeFiles/arquivobinario.dir/src/Reaction.cc.o: ../src/Reaction.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/leo/Desktop/RIBRAS/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object CMakeFiles/arquivobinario.dir/src/Reaction.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/arquivobinario.dir/src/Reaction.cc.o -c /home/leo/Desktop/RIBRAS/src/Reaction.cc

CMakeFiles/arquivobinario.dir/src/Reaction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/arquivobinario.dir/src/Reaction.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/leo/Desktop/RIBRAS/src/Reaction.cc > CMakeFiles/arquivobinario.dir/src/Reaction.cc.i

CMakeFiles/arquivobinario.dir/src/Reaction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/arquivobinario.dir/src/Reaction.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/leo/Desktop/RIBRAS/src/Reaction.cc -o CMakeFiles/arquivobinario.dir/src/Reaction.cc.s

CMakeFiles/arquivobinario.dir/src/RootSaver.cc.o: CMakeFiles/arquivobinario.dir/flags.make
CMakeFiles/arquivobinario.dir/src/RootSaver.cc.o: ../src/RootSaver.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/leo/Desktop/RIBRAS/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building CXX object CMakeFiles/arquivobinario.dir/src/RootSaver.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/arquivobinario.dir/src/RootSaver.cc.o -c /home/leo/Desktop/RIBRAS/src/RootSaver.cc

CMakeFiles/arquivobinario.dir/src/RootSaver.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/arquivobinario.dir/src/RootSaver.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/leo/Desktop/RIBRAS/src/RootSaver.cc > CMakeFiles/arquivobinario.dir/src/RootSaver.cc.i

CMakeFiles/arquivobinario.dir/src/RootSaver.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/arquivobinario.dir/src/RootSaver.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/leo/Desktop/RIBRAS/src/RootSaver.cc -o CMakeFiles/arquivobinario.dir/src/RootSaver.cc.s

CMakeFiles/arquivobinario.dir/src/RunAction.cc.o: CMakeFiles/arquivobinario.dir/flags.make
CMakeFiles/arquivobinario.dir/src/RunAction.cc.o: ../src/RunAction.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/leo/Desktop/RIBRAS/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Building CXX object CMakeFiles/arquivobinario.dir/src/RunAction.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/arquivobinario.dir/src/RunAction.cc.o -c /home/leo/Desktop/RIBRAS/src/RunAction.cc

CMakeFiles/arquivobinario.dir/src/RunAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/arquivobinario.dir/src/RunAction.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/leo/Desktop/RIBRAS/src/RunAction.cc > CMakeFiles/arquivobinario.dir/src/RunAction.cc.i

CMakeFiles/arquivobinario.dir/src/RunAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/arquivobinario.dir/src/RunAction.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/leo/Desktop/RIBRAS/src/RunAction.cc -o CMakeFiles/arquivobinario.dir/src/RunAction.cc.s

CMakeFiles/arquivobinario.dir/src/SensitiveDetector.cc.o: CMakeFiles/arquivobinario.dir/flags.make
CMakeFiles/arquivobinario.dir/src/SensitiveDetector.cc.o: ../src/SensitiveDetector.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/leo/Desktop/RIBRAS/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_14) "Building CXX object CMakeFiles/arquivobinario.dir/src/SensitiveDetector.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/arquivobinario.dir/src/SensitiveDetector.cc.o -c /home/leo/Desktop/RIBRAS/src/SensitiveDetector.cc

CMakeFiles/arquivobinario.dir/src/SensitiveDetector.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/arquivobinario.dir/src/SensitiveDetector.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/leo/Desktop/RIBRAS/src/SensitiveDetector.cc > CMakeFiles/arquivobinario.dir/src/SensitiveDetector.cc.i

CMakeFiles/arquivobinario.dir/src/SensitiveDetector.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/arquivobinario.dir/src/SensitiveDetector.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/leo/Desktop/RIBRAS/src/SensitiveDetector.cc -o CMakeFiles/arquivobinario.dir/src/SensitiveDetector.cc.s

CMakeFiles/arquivobinario.dir/src/SiHit.cc.o: CMakeFiles/arquivobinario.dir/flags.make
CMakeFiles/arquivobinario.dir/src/SiHit.cc.o: ../src/SiHit.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/leo/Desktop/RIBRAS/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_15) "Building CXX object CMakeFiles/arquivobinario.dir/src/SiHit.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/arquivobinario.dir/src/SiHit.cc.o -c /home/leo/Desktop/RIBRAS/src/SiHit.cc

CMakeFiles/arquivobinario.dir/src/SiHit.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/arquivobinario.dir/src/SiHit.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/leo/Desktop/RIBRAS/src/SiHit.cc > CMakeFiles/arquivobinario.dir/src/SiHit.cc.i

CMakeFiles/arquivobinario.dir/src/SiHit.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/arquivobinario.dir/src/SiHit.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/leo/Desktop/RIBRAS/src/SiHit.cc -o CMakeFiles/arquivobinario.dir/src/SiHit.cc.s

CMakeFiles/arquivobinario.dir/src/StepMax.cc.o: CMakeFiles/arquivobinario.dir/flags.make
CMakeFiles/arquivobinario.dir/src/StepMax.cc.o: ../src/StepMax.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/leo/Desktop/RIBRAS/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_16) "Building CXX object CMakeFiles/arquivobinario.dir/src/StepMax.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/arquivobinario.dir/src/StepMax.cc.o -c /home/leo/Desktop/RIBRAS/src/StepMax.cc

CMakeFiles/arquivobinario.dir/src/StepMax.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/arquivobinario.dir/src/StepMax.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/leo/Desktop/RIBRAS/src/StepMax.cc > CMakeFiles/arquivobinario.dir/src/StepMax.cc.i

CMakeFiles/arquivobinario.dir/src/StepMax.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/arquivobinario.dir/src/StepMax.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/leo/Desktop/RIBRAS/src/StepMax.cc -o CMakeFiles/arquivobinario.dir/src/StepMax.cc.s

CMakeFiles/arquivobinario.dir/src/StepMaxMessenger.cc.o: CMakeFiles/arquivobinario.dir/flags.make
CMakeFiles/arquivobinario.dir/src/StepMaxMessenger.cc.o: ../src/StepMaxMessenger.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/leo/Desktop/RIBRAS/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_17) "Building CXX object CMakeFiles/arquivobinario.dir/src/StepMaxMessenger.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/arquivobinario.dir/src/StepMaxMessenger.cc.o -c /home/leo/Desktop/RIBRAS/src/StepMaxMessenger.cc

CMakeFiles/arquivobinario.dir/src/StepMaxMessenger.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/arquivobinario.dir/src/StepMaxMessenger.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/leo/Desktop/RIBRAS/src/StepMaxMessenger.cc > CMakeFiles/arquivobinario.dir/src/StepMaxMessenger.cc.i

CMakeFiles/arquivobinario.dir/src/StepMaxMessenger.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/arquivobinario.dir/src/StepMaxMessenger.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/leo/Desktop/RIBRAS/src/StepMaxMessenger.cc -o CMakeFiles/arquivobinario.dir/src/StepMaxMessenger.cc.s

# Object files for target arquivobinario
arquivobinario_OBJECTS = \
"CMakeFiles/arquivobinario.dir/Example.cc.o" \
"CMakeFiles/arquivobinario.dir/src/BeamTestPrimaryGeneratorAction.cc.o" \
"CMakeFiles/arquivobinario.dir/src/BeamTestPrimaryGeneratorMessenger.cc.o" \
"CMakeFiles/arquivobinario.dir/src/DetectorConstruction.cc.o" \
"CMakeFiles/arquivobinario.dir/src/DetectorMessenger.cc.o" \
"CMakeFiles/arquivobinario.dir/src/EventAction.cc.o" \
"CMakeFiles/arquivobinario.dir/src/MagneticField.cc.o" \
"CMakeFiles/arquivobinario.dir/src/MagneticFieldMessenger.cc.o" \
"CMakeFiles/arquivobinario.dir/src/PhysicsList.cc.o" \
"CMakeFiles/arquivobinario.dir/src/PhysicsListMessenger.cc.o" \
"CMakeFiles/arquivobinario.dir/src/Reaction.cc.o" \
"CMakeFiles/arquivobinario.dir/src/RootSaver.cc.o" \
"CMakeFiles/arquivobinario.dir/src/RunAction.cc.o" \
"CMakeFiles/arquivobinario.dir/src/SensitiveDetector.cc.o" \
"CMakeFiles/arquivobinario.dir/src/SiHit.cc.o" \
"CMakeFiles/arquivobinario.dir/src/StepMax.cc.o" \
"CMakeFiles/arquivobinario.dir/src/StepMaxMessenger.cc.o"

# External object files for target arquivobinario
arquivobinario_EXTERNAL_OBJECTS =

arquivobinario: CMakeFiles/arquivobinario.dir/Example.cc.o
arquivobinario: CMakeFiles/arquivobinario.dir/src/BeamTestPrimaryGeneratorAction.cc.o
arquivobinario: CMakeFiles/arquivobinario.dir/src/BeamTestPrimaryGeneratorMessenger.cc.o
arquivobinario: CMakeFiles/arquivobinario.dir/src/DetectorConstruction.cc.o
arquivobinario: CMakeFiles/arquivobinario.dir/src/DetectorMessenger.cc.o
arquivobinario: CMakeFiles/arquivobinario.dir/src/EventAction.cc.o
arquivobinario: CMakeFiles/arquivobinario.dir/src/MagneticField.cc.o
arquivobinario: CMakeFiles/arquivobinario.dir/src/MagneticFieldMessenger.cc.o
arquivobinario: CMakeFiles/arquivobinario.dir/src/PhysicsList.cc.o
arquivobinario: CMakeFiles/arquivobinario.dir/src/PhysicsListMessenger.cc.o
arquivobinario: CMakeFiles/arquivobinario.dir/src/Reaction.cc.o
arquivobinario: CMakeFiles/arquivobinario.dir/src/RootSaver.cc.o
arquivobinario: CMakeFiles/arquivobinario.dir/src/RunAction.cc.o
arquivobinario: CMakeFiles/arquivobinario.dir/src/SensitiveDetector.cc.o
arquivobinario: CMakeFiles/arquivobinario.dir/src/SiHit.cc.o
arquivobinario: CMakeFiles/arquivobinario.dir/src/StepMax.cc.o
arquivobinario: CMakeFiles/arquivobinario.dir/src/StepMaxMessenger.cc.o
arquivobinario: CMakeFiles/arquivobinario.dir/build.make
arquivobinario: /home/leo/G4/geant4-install/lib/libG4Tree.so
arquivobinario: /home/leo/G4/geant4-install/lib/libG4GMocren.so
arquivobinario: /home/leo/G4/geant4-install/lib/libG4visHepRep.so
arquivobinario: /home/leo/G4/geant4-install/lib/libG4RayTracer.so
arquivobinario: /home/leo/G4/geant4-install/lib/libG4VRML.so
arquivobinario: /home/leo/G4/geant4-install/lib/libG4OpenGL.so
arquivobinario: /home/leo/G4/geant4-install/lib/libG4gl2ps.so
arquivobinario: /home/leo/G4/geant4-install/lib/libG4interfaces.so
arquivobinario: /home/leo/G4/geant4-install/lib/libG4persistency.so
arquivobinario: /home/leo/G4/geant4-install/lib/libG4error_propagation.so
arquivobinario: /home/leo/G4/geant4-install/lib/libG4readout.so
arquivobinario: /home/leo/G4/geant4-install/lib/libG4physicslists.so
arquivobinario: /home/leo/G4/geant4-install/lib/libG4parmodels.so
arquivobinario: /home/leo/G4/geant4-install/lib/libG4FR.so
arquivobinario: /home/leo/G4/geant4-install/lib/libG4vis_management.so
arquivobinario: /home/leo/G4/geant4-install/lib/libG4modeling.so
arquivobinario: /usr/lib/x86_64-linux-gnu/libGLU.so
arquivobinario: /usr/lib/x86_64-linux-gnu/libGL.so
arquivobinario: /usr/lib/x86_64-linux-gnu/libXmu.so
arquivobinario: /usr/lib/x86_64-linux-gnu/libXext.so
arquivobinario: /usr/lib/x86_64-linux-gnu/libXt.so
arquivobinario: /usr/lib/x86_64-linux-gnu/libSM.so
arquivobinario: /usr/lib/x86_64-linux-gnu/libICE.so
arquivobinario: /usr/lib/x86_64-linux-gnu/libX11.so
arquivobinario: /usr/lib/x86_64-linux-gnu/libQtOpenGL.so
arquivobinario: /usr/lib/x86_64-linux-gnu/libQtGui.so
arquivobinario: /usr/lib/x86_64-linux-gnu/libQtCore.so
arquivobinario: /home/leo/G4/geant4-install/lib/libG4run.so
arquivobinario: /home/leo/G4/geant4-install/lib/libG4event.so
arquivobinario: /home/leo/G4/geant4-install/lib/libG4tracking.so
arquivobinario: /home/leo/G4/geant4-install/lib/libG4processes.so
arquivobinario: /home/leo/G4/geant4-install/lib/libG4analysis.so
arquivobinario: /usr/lib/x86_64-linux-gnu/libz.so
arquivobinario: /usr/lib/x86_64-linux-gnu/libexpat.so
arquivobinario: /home/leo/G4/geant4-install/lib/libG4digits_hits.so
arquivobinario: /home/leo/G4/geant4-install/lib/libG4track.so
arquivobinario: /home/leo/G4/geant4-install/lib/libG4particles.so
arquivobinario: /home/leo/G4/geant4-install/lib/libG4geometry.so
arquivobinario: /home/leo/G4/geant4-install/lib/libG4materials.so
arquivobinario: /home/leo/G4/geant4-install/lib/libG4graphics_reps.so
arquivobinario: /home/leo/G4/geant4-install/lib/libG4intercoms.so
arquivobinario: /home/leo/G4/geant4-install/lib/libG4global.so
arquivobinario: /home/leo/CLHEP/CLHEP-install/lib/libCLHEP-Evaluator-2.4.1.3.so
arquivobinario: /home/leo/CLHEP/CLHEP-install/lib/libCLHEP-Geometry-2.4.1.3.so
arquivobinario: /home/leo/CLHEP/CLHEP-install/lib/libCLHEP-Vector-2.4.1.3.so
arquivobinario: /home/leo/CLHEP/CLHEP-install/lib/libCLHEP-Random-2.4.1.3.so
arquivobinario: CMakeFiles/arquivobinario.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/leo/Desktop/RIBRAS/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_18) "Linking CXX executable arquivobinario"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/arquivobinario.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/arquivobinario.dir/build: arquivobinario

.PHONY : CMakeFiles/arquivobinario.dir/build

CMakeFiles/arquivobinario.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/arquivobinario.dir/cmake_clean.cmake
.PHONY : CMakeFiles/arquivobinario.dir/clean

CMakeFiles/arquivobinario.dir/depend:
	cd /home/leo/Desktop/RIBRAS/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/leo/Desktop/RIBRAS /home/leo/Desktop/RIBRAS /home/leo/Desktop/RIBRAS/build /home/leo/Desktop/RIBRAS/build /home/leo/Desktop/RIBRAS/build/CMakeFiles/arquivobinario.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/arquivobinario.dir/depend

