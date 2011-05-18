# Find netcdf_cpp include directories and libraries
#
# NETCDF_INCLUDE_DIRECTORIES - where to find netcdf.h
# NETCDF_LIBRARIES - list of libraries to link against when using NetCDF
# NETCDF_FOUND - Do not attempt to use NetCDF if "no", "0", or undefined.

include(LibFindMacros)

# Dependencies 
set(NETCDF_CPP_PREFIX "/usr/local" 
    CACHE PATH "Path to search for NetCDF header and library files" )
libfind_package(NETCDF_CPP NETCDF)

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