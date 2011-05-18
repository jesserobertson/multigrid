# Find netcdf include directories and libraries
#
# NETCDF_INCLUDE_DIRECTORIES - where to find netcdf.h
# NETCDF_LIBRARIES - list of libraries to link against when using NetCDF
# NETCDF_FOUND - Do not attempt to use NetCDF if "no", "0", or undefined.

include(LibFindMacros)

# Dependencies 
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