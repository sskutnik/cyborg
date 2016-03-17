#  SCALE_FOUND - system has the SCALE shared libraries
#  SCALE_INCLUDE_DIR - the SCALE include directory
#  SCALE_LIBRARIES - The libraries needed to use the ORIGEN/SCALE API

# Check if we have an environment variable to SCALE
IF(DEFINED ENV{SCALE})
    IF(NOT DEFINED SCALE_ROOT_DIR)
        SET(SCALE_ROOT_DIR "$ENV{SCALE}")
    ELSE(NOT DEFINED SCALE_ROOT_DIR)
        # Yell if both exist
        MESSAGE(STATUS "\tTwo SCALE_ROOT_DIRs have been found:")
        MESSAGE(STATUS "\t\tThe defined cmake variable SCALE_ROOT_DIR: ${SCALE_ROOT_DIR}")
        MESSAGE(STATUS "\t\tThe environment variable SCALE_ROOT_DIR: $ENV{SCALE}")
    ENDIF(NOT DEFINED SCALE_ROOT_DIR)
ENDIF(DEFINED ENV{SCALE})

# Let the user know if we're using a hint
MESSAGE(STATUS "Using ${SCALE_ROOT_DIR} as SCALE_ROOT_DIR.")

# Set the include dir, this will be the future basis for other
# defined dirs
FIND_PATH(ORIGEN_INCLUDE_DIR origen_config.h
    HINTS ${SCALE_ROOT_DIR}
    "${SCALE_ROOT_DIR}/${CMAKE_SYSTEM_NAME}_${CMAKE_SYSTEM_PROCESSOR}"
    PATHS /opt/scale6.2b5
    PATH_SUFFIXES include include/Origen/Core )

# Add the root dir to the hints
SET(SCALE_INCLUDE_DIR "${ORIGEN_INCLUDE_DIR}/../..")

# Look for the library
FIND_LIBRARY(ORIGEN_CORE_LIBRARY NAMES OrigenCore
    HINTS ${SCALE_ROOT_DIR}
    "${SCALE_ROOT_DIR}/${CMAKE_SYSTEM_NAME}_${CMAKE_SYSTEM_PROCESSOR}"
    PATHS /opt/scale6.2b5
    PATH_SUFFIXES lib)

# Copy the results to the output variables.
IF(SCALE_INCLUDE_DIR AND ORIGEN_CORE_LIBRARY)
    SET(SCALE_FOUND 1)
    SET(SCALE_LIBRARIES ${ORIGEN_CORE_LIBRARY})
    SET(SCALE_INCLUDE_DIRS ${SCALE_INCLUDE_DIR})
ELSE()
    SET(SCALE_FOUND 0)
    SET(SCALE_LIBRARIES)
    SET(SCALE_INCLUDE_DIRS)
ENDIF()

# Report the results.
IF(SCALE_FOUND)
    SET(SCALE_DIR_MESSAGE "Found SCALE Headers : "
        ${SCALE_INCLUDE_DIRS} )
    SET(SCALE_LIB_MESSAGE "Found SCALE Libraries : "
        ${SCALE_LIBRARIES} )
    MESSAGE(STATUS ${SCALE_DIR_MESSAGE})
    MESSAGE(STATUS ${SCALE_LIB_MESSAGE})
ELSE(SCALE_FOUND)
    SET(SCALE_DIR_MESSAGE
        "SCALE was not found. Make sure SCALE_LIBRARY and SCALE_INCLUDE_DIR are set.")
    IF(NOT SCALE_FIND_QUIETLY)
        MESSAGE(STATUS "${SCALE_DIR_MESSAGE}")
        MESSAGE(STATUS "SCALE_INCLUDE_DIR was set to : ${SCALE_INCLUDE_DIR}")
        MESSAGE(STATUS "SCALE_LIBRARY was set to : ${SCALE_LIBRARIES}")
        IF(SCALE_FIND_REQUIRED)
            MESSAGE(FATAL_ERROR "${SCALE_DIR_MESSAGE}")
        ENDIF(SCALE_FIND_REQUIRED)
    ENDIF(NOT SCALE_FIND_QUIETLY)

ENDIF(SCALE_FOUND)

MARK_AS_ADVANCED(
    SCALE_INCLUDE_DIR
    SCALE_LIBRARY
    )
