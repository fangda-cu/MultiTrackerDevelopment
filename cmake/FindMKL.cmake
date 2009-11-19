# - Find MKL
# Find the MKL headers and libraries
#
# MKL_INCLUDE_DIR - include path for headers
# MKL_LIBRARIES   - libraries to include when linking
# MKL_FOUND       - True if MKL is found

if (MKL_INCLUDE_DIR AND MKL_LIBRARIES)
  # already in cache, be silent
  set (MKL_FIND_QUIETLY TRUE)
endif (MKL_INCLUDE_DIR AND MKL_LIBRARIES)


##############################################################
##
## Search for the header location. There are hundreds of
## headers so just search for one we know must be there and 
## assume it is a valid install.
##
##############################################################
find_path (MKL_INCLUDE_PATH
  /include/mkl.h
  HINTS ENV MKLROOT
  )

##############################################################
##
## Search for all the libraries to make sure they can be found
##
##############################################################
set (MKL_LIBRARIES_FOUND)
set (MKL_LIBRARIES_MISSING)
set (MKL_LIBS libmkl_intel_lp64.so libmkl_intel_thread.so libmkl_core.so libiomp5.so libmkl_lapack.so)
foreach (MKLLIB ${MKL_LIBS})
    set (MKL_SEARCH_LIB "MKL_SEARCH_LIB-NOTFOUND" CACHE FILEPATH "Cleared" FORCE)
    find_library (MKL_SEARCH_LIB ${MKLLIB} PATHS $ENV{MKLROOT}/lib/em64t/)
    if (MKL_SEARCH_LIB)
        list (APPEND MKL_LIBRARIES_FOUND ${MKL_SEARCH_LIB})        
    else (MKL_SEARCH_LIB)
        list (APPEND MKL_LIBRARIES_MISSING ${MKL_SEARCH_LIB})
        message (SEND_ERROR "Unable to find MKL library ${MKLLIB}")
    endif (MKL_SEARCH_LIB)
endforeach (MKLLIB)

set (MKL_LIBRARY ${MKL_LIBRARIES_FOUND} CACHE STRING "MKL libraries" FORCE)

# handle the QUIETLY and REQUIRED arguments and set
# MKL_FOUND to TRUE if all listed variables are TRUE
include (FindPackageHandleStandardArgs)
#find_package_handle_standard_args (MKL "MKL not found. Set MKL_LOCATION " MKL_INCLUDE_PATH MKL_LIBRARIES)
find_package_handle_standard_args (MKL "MKL not found. Set MKL_LOCATION " MKL_INCLUDE_PATH MKL_LIBRARY)

if (MKL_FOUND)
  set (MKL_INCLUDE_DIR ${MKL_INCLUDE_PATH}/include)
  set (MKL_LIBRARIES ${MKL_LIBRARY})
else (MKL_FOUND)
  set (MKL_INCLUDE_DIR)
  set (MKL_LIBRARIES)
endif (MKL_FOUND)