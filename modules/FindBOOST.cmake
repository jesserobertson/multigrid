# Find boost include directories and libraries
#
# BOOST_INCLUDE_DIRECTORIES - where to find boost headers
# BOOST_FOUND - Do not attempt to use NetCDF if "no", "0", or undefined. 
#
# N.B. Only boost headers are required for this library.

include(LibFindMacros)

# Dependencies 
set(BOOST_PREFIX "/usr/local" 
    CACHE PATH "Path to search for boost header and library files" )

# Find include directories
find_path(BOOST_INCLUDE_DIR
    NAMES boost/foreach.hpp PATHS ${BOOST_PREFIX})

# Set the include dir variables and the libraries and let libfind_process do
# the rest. NOTE: Singular variables for this library, plural for libraries
# this this lib depends on.
set(BOOST_PROCESS_INCLUDES BOOST_INCLUDE_DIR)
set(BOOST_PROCESS_LIBS )
libfind_process(BOOST)