# ===============================================================
# = CMake file for Molosov problem - Jess Robertson, 2010-05-20 =
# ===============================================================
project(mosolov) 
include("modules/LibFindMacros.cmake")

# Configure cmake build
cmake_minimum_required(VERSION 2.8.1 FATAL_ERROR)        
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "modules") 
set(CMAKE_VERBOSE_MAKEFILE false)
set(CMAKE_C_COMPILER gcc) 
set(CMAKE_CXX_COMPILER g++) 

# Set source file directory
set(SOURCE_DIRECTORY .)

# Find and add libraries and headers to build script 
# Set the include dir variables and the libraries and let libfind_process do
# the rest. NOTE: Singular variables for this library, plural for libraries
# this this lib depends on. 
find_package(MULTIGRID REQUIRED) 
set(${PROJECT_NAME}_PROCESS_INCLUDES MULTIGRID_INCLUDE_DIRS)
set(${PROJECT_NAME}_PROCESS_LIBS MULTIGRID_LIBRARIES)
libfind_process(${PROJECT_NAME})

# Set up project    
include_directories(${SOURCE_DIR})
file(GLOB SRC ${SOURCE_DIRECTORY}/*.cpp)    # Glob for source files 
add_executable(${PROJECT_NAME} ${SRC})   

include_directories(${${PROJECT_NAME}_INCLUDE_DIRS})
link_directories(${${PROJECT_NAME}_LIBRARY_DIRS})
target_link_libraries(${PROJECT_NAME} ${${PROJECT_NAME}_LIBRARIES})
set_target_properties(${PROJECT_NAME} 
    PROPERTIES COMPILER_FLAGS "-fast -m64 -arch x86_64 -msse -Wall -pedantic"
               LINKER_FLAGS "-fast -m64 -arch x86_64 -msse")