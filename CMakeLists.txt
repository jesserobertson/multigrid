# =================================================================
# = CMake file for multigrid library - Jess Robertson, 2010-05-20 =
# =================================================================
project(multigrid)      

# Configure cmake build
cmake_minimum_required(VERSION 2.8.1 FATAL_ERROR)        
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "modules") 
set(CMAKE_VERBOSE_MAKEFILE false)
set(CMAKE_C_COMPILER gcc) 
set(CMAKE_CXX_COMPILER g++) 

# Set source file directory  
set(source_directory .)   

# Find and add libraries and headers to build script                     
find_path(BOOST_INCLUDE_DIR NAMES boost/foreach.hpp)
find_package(BLITZ REQUIRED)
find_package(NETCDF_CPP REQUIRED)

# Decide what to build
set(build_library true) 
IF(${build_library}) 
    # Get headers & sources but not main.cpp  
    file(GLOB headers ${source_directory}/*.hpp)
    file(GLOB sources 
        ${source_directory}/boundary_conditions.cpp
        ${source_directory}/fdarray.cpp
        ${source_directory}/fdbase.cpp
        ${source_directory}/fdvecarray.cpp
        ${source_directory}/multigrid_base.cpp
        ${source_directory}/multigrid_linear.cpp
        ${source_directory}/multigrid_nonlinear.cpp
        ${source_directory}/stack.cpp
        ${source_directory}/settings.cpp
        ${source_directory}/boundary_conditions.cpp)    
    
    # Set up library    
    add_library(${PROJECT_NAME} ${sources})
    include_directories(${INCLUDES} ${source_directory} 
        ${BOOST_INCLUDE_DIR} ${BLITZ_INCLUDE_DIRS} ${NETCDF_INCLUDE_DIRS})
    set_target_properties(${PROJECT_NAME} 
        PROPERTIES COMPILER_FLAGS "-fast -m64 -arch x86_64 -msse -Wall -pedantic") 

    # Install commands
    install(TARGETS multigrid
        LIBRARY DESTINATION /usr/local/lib
        ARCHIVE DESTINATION /usr/local/lib
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                    GROUP_READ GROUP_EXECUTE)    
    install(FILES ${headers}
        DESTINATION /usr/local/include/multigrid
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                    GROUP_READ GROUP_EXECUTE)
ELSE(${build_library}) 
    # Get sources
    file(GLOB sources ${source_directory}/*.cpp)
    
    # Set up executable 
    add_executable(${PROJECT_NAME} ${sources})
    include_directories(${INCLUDES} ${source_directory} 
        ${BOOST_INCLUDE_DIR} ${BLITZ_INCLUDE_DIRS} ${NETCDF_INCLUDE_DIRS}) 
    target_link_libraries(${PROJECT_NAME} ${BLITZ_LIBRARIES} ${NETCDF_CPP_LIBRARIES})
    set_target_properties(${PROJECT_NAME} 
        PROPERTIES COMPILER_FLAGS "-g -m64 -arch x86_64 -msse -Wall -pedantic"
                   LINKER_FLAGS "-m64 -arch x86_64 -msse")
ENDIF(${build_library})          