# Find multigrid include directories and libraries
#
# MULTIGRID_INCLUDE_DIRECTORIES - where to find netcdf.h
# MULTIGRID_LIBRARIES - list of libraries to link against when using NetCDF
# MULTIGRID_FOUND - Do not attempt to use NetCDF if "no", "0", or undefined.

include(LibFindMacros)

# Dependencies 
set(MULTIGRID_PREFIX "/usr/local/" 
    CACHE PATH "Path to search for Multigrid header and library files" )

# ==============
# = Find Blitz =
# ==============       
# Dependencies 
set( BLITZ_PREFIX "/usr/local/" 
    CACHE PATH "Path to search for Blitz++ header and library files" ) 

# Include dir
find_path(BLITZ_INCLUDE_DIR NAMES blitz/blitz.h PATHS ${BLITZ_PREFIX})   

# Finally the library itself
find_library(BLITZ_LIBRARY NAMES blitz PATHS ${BLITZ_PREFIX})

# Set the include dir variables and the libraries and let libfind_process do
# the rest. NOTE: Singular variables for this library, plural for libraries
# this this lib depends on.
set(BLITZ_PROCESS_INCLUDES BLITZ_INCLUDE_DIR)
set(BLITZ_PROCESS_LIBS BLITZ_LIBRARY)
libfind_process(BLITZ)

# ===============
# = Find NetCDF =
# ===============
# Dependencies 
set(NETCDF_CPP_PREFIX "/usr/local" 
    CACHE PATH "Path to search for NetCDF header and library files" )
set( NETCDF_PREFIX "/usr/local" 
    CACHE PATH "Path to search for NetCDF header and library files" ) 
find_package(CURL REQUIRED)
find_package(HDF5 REQUIRED)
find_package(ZLIB REQUIRED)

# Include dir
find_path(NETCDF_INCLUDE_DIR NAMES netcdf.h PATHS ${NETCDF_PREFIX})

# Finally the library itself
find_library(NETCDF_LIBRARY NAMES netcdf PATHS ${NETCDF_PREFIX})

# Set the include dir variables and the libraries and let libfind_process do
# the rest. NOTE: Singular variables for this library, plural for libraries
# this this lib depends on.
set(NETCDF_PROCESS_INCLUDES 
    NETCDF_INCLUDE_DIR CURL_INCLUDE_DIRS HDF5_INCLUDE_DIRS ZLIB_INCLUDE_DIRS)
set(NETCDF_PROCESS_LIBS 
    NETCDF_LIBRARY CURL_LIBRARIES HDF5_LIBRARIES ZLIB_LIBRARIES)
libfind_process(NETCDF)   

# Include dir
find_path(NETCDF_CPP_INCLUDE_DIR NAMES netcdf.h PATHS ${NETCDF_CPP_PREFIX})

# Finally the library itself
find_library(NETCDF_CPP_LIBRARY NAMES netcdf_c++ PATHS ${NETCDF_CPP_PREFIX})

# Set the include dir variables and the libraries and let libfind_process do
# the rest. NOTE: Singular variables for this library, plural for libraries
# this this lib depends on.
set(NETCDF_CPP_PROCESS_INCLUDES NETCDF_CPP_INCLUDE_DIR NETCDF_INCLUDE_DIRS)
set(NETCDF_CPP_PROCESS_LIBS NETCDF_CPP_LIBRARY NETCDF_LIBRARIES)
libfind_process(NETCDF_CPP)  

# ==================
# = Find Multigrid =
# ==================
# Include dir
find_path(MULTIGRID_INCLUDE_DIR NAMES multigrid/multigrid.hpp 
    PATHS ${MULTIGRID_PREFIX}/include)

# Finally the library itself
find_library(MULTIGRID_LIBRARY NAMES multigrid PATHS ${MULTIGRID_PREFIX}/lib)

# Set the include dir variables and the libraries and let libfind_process do
# the rest. NOTE: Singular variables for this library, plural for libraries
# this this lib depends on.
set(MULTIGRID_PROCESS_INCLUDES 
    MULTIGRID_INCLUDE_DIR 
    BLITZ_INCLUDE_DIRS
    NETCDF_CPP_INCLUDE_DIRS)
set(MULTIGRID_PROCESS_LIBS 
    MULTIGRID_LIBRARY 
    BLITZ_LIBRARIES 
    NETCDF_CPP_LIBRARIES)
libfind_process(MULTIGRID)