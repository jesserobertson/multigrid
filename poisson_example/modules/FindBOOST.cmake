# Find boost include directories and libraries
#
# BOOST_INCLUDE_DIRECTORIES - where to find netcdf.h
# BOOST_LIBRARIES - list of libraries to link against when using NetCDF
# BOOST_FOUND - Do not attempt to use NetCDF if "no", "0", or undefined.

include(LibFindMacros)

# Dependencies 
set(BOOST_PREFIX "/usr/local" 
    CACHE PATH "Path to search for boost header and library files" )

# Find include directories
find_path(BOOST_PROGRAM_OPTIONS_INCLUDE_DIR 
    NAMES boost/program_options.hpp PATHS ${BOOST_PREFIX}) 

# Finally the libraries themselves
find_library(BOOST_PROGRAM_OPTIONS_LIBRARY 
    NAMES boost_program_options-mt PATHS ${BOOST_PREFIX})

# Set the include dir variables and the libraries and let libfind_process do
# the rest. NOTE: Singular variables for this library, plural for libraries
# this this lib depends on.
set(BOOST_PROCESS_INCLUDES BOOST_PROGRAM_OPTIONS_INCLUDE_DIR)
set(BOOST_PROCESS_LIBS BOOST_PROGRAM_OPTIONS_LIBRARY)
libfind_process(BOOST)