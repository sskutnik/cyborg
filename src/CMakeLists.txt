#################################################################
######### build libOrigenInterface ##############################
#################################################################

INCLUDE_DIRECTORIES(
    "${CMAKE_CURRENT_SOURCE_DIR}"
    )


FILE(GLOB origen_interface_src "${CMAKE_CURRENT_SOURCE_DIR}/[^_]origen*.cpp")

ADD_LIBRARY(OrigenInterface ${ORIGEN_INTERFACE_SRC})

SET_TARGET_PROPERTIES(OrigenInterface
    PROPERTIES
    INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib"
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib"
    )

TARGET_LINK_LIBRARIES(OrigenInterface ${LIBS})
set_target_properties(OrigenInterface PROPERTIES LINKER_LANGUAGE CXX)

INSTALL(
    TARGETS OrigenInterface
    LIBRARY DESTINATION lib
    COMPONENT OrigenInterface
    )


SET(LIBS ${LIBS} OrigenInterface PARENT_SCOPE)


#################################################################

### DO NOT DELETE THIS COMMENT: INSERT_ARCHETYPES_HERE ###
SET(CYCLUS_CUSTOM_HEADERS "cyclus_origen_interface.h")
USE_CYCLUS("cyborg" "reactor")

INSTALL_CYCLUS_MODULE("cyborg" "" "NONE")

# install header files
FILE(GLOB h_files "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
INSTALL(FILES ${h_files} DESTINATION include/cyborg COMPONENT cyborg)

