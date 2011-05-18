# Find Blitz++ include directories
#
# BLITZ_INCLUDE_DIRECTORIES - where to find blitz.hpp
# BLITZ_FOUND - Do not attempt to use blitz++ if "no", "0", or undefined.

include(LibFindMacros)

# Dependencies 
set( BLITZ_PREFIX "/usr/local" 
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