#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(Example)

#set(CMAKE_CXX_STANDARD 11)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -Wall -Wextra")


# GEANT4 core
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)


#Find ROOT
find_package(ROOT REQUIRED)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}  -pthread -Wno-deprecated-declarations")
include_directories(${ROOT_INCLUDE_DIR}  ${CMAKE_SOURCE_DIR})
link_directories(${ROOT_LIBRARY_DIR})
message( ${ROOT_LIBRARIES} )


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build")


## User code
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)


#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#

#set(ROOT_LIBRARY_PROPERTIES ${ROOT_LIBRARY_PROPERTIES}) ## test DFlechas
add_executable(arquivobinario Example.cc ${sources} ${headers})
target_link_libraries(arquivobinario ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
#target_link_libraries(USP_RIBRAS ${cadmesh_LIBRARIES})
#target_link_libraries(USP_RIBRAS ${Geant4_LIBRARIES} ${ROOT_LIBRARIES}) ## test DFlechas


#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS arquivobinario DESTINATION bin)

