#################################################################
######### build libOrigenInterface ##############################
#################################################################

# Configure the default ORIGEN library directory
IF(DEFINED USER_ORGLIB_DIR) 
   SET(ORGLIB_DIR "${USER_ORGLIB_DIR}")
ELSE() 
   SET(ORGLIB_DIR "$ENV{DATA}/arplibs/")
ENDIF()
MESSAGE("-- Setting ORGLIB_DIR to: ${ORGLIB_DIR}")
CONFIGURE_FILE(orglib_default_location.h.in "${CMAKE_CURRENT_SOURCE_DIR}/orglib_default_location.h" @ONLY)

INCLUDE_DIRECTORIES(
    "${CMAKE_CURRENT_SOURCE_DIR}"
    )

FILE(GLOB ORIGEN_INTERFACE_SRC "${CMAKE_CURRENT_SOURCE_DIR}/*origen*.cpp")

ADD_LIBRARY(OrigenInterface SHARED ${ORIGEN_INTERFACE_SRC})

SET_TARGET_PROPERTIES(OrigenInterface
    PROPERTIES
    INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib"
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib"
    LINKER_LANGUAGE CXX
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED ON
    )
TARGET_LINK_LIBRARIES(OrigenInterface ${LIBS} ${SCALE_LIBRARIES})

INSTALL(
    TARGETS OrigenInterface
    LIBRARY DESTINATION lib
    COMPONENT OrigenInterface
    )

#SET(ORIGEN_LIB OrigenInterface PARENT_SCOPE)
#SET(LIBS_NO_SCALE ${LIBS} PARENT_SCOPE)
SET(LIBS ${LIBS} OrigenInterface)


#################################################################

### DO NOT DELETE THIS COMMENT: INSERT_ARCHETYPES_HERE ###

SET(CYCLUS_CUSTOM_HEADERS "cyclus_origen_interface.h" "orglib_default_location.h")

# Want to install OrigenInterface *only* for Cyborg archetypes and NOT unit test driver

#SET(LIBS ${LIBS} ${OrigenInterface} PARENT_SCOPE)
#MESSAGE("Building cyborg: LIBS = ${LIBS}")
USE_CYCLUS("cyborg" "reactor")

INSTALL_CYCLUS_MODULE("cyborg" "" "NONE")
SET(CYBORG_LIB ${cyborg_LIB} PARENT_SCOPE)

SET(TestSource ${cyborg_TEST_CC} PARENT_SCOPE)

# install header files
FILE(GLOB h_files "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
INSTALL(FILES ${h_files} DESTINATION include/cyborg COMPONENT cyborg)
