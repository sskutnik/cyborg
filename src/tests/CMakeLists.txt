# TODO: GTest-ify this test & add a gtest dependency

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} "${PROJECT_SOURCE_DIR}/src")

add_executable(tstCyclusOrigenInterface tst_cyclus_origen_interface.cpp )
SET_TARGET_PROPERTIES( tstCyclusOrigenInterface
    PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib"
    LINKER_LANGUAGE CXX
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED ON
    )
TARGET_LINK_LIBRARIES( tstCyclusOrigenInterface OrigenInterface ${LIBS} )
MESSAGE("TARGETING ${LIBS}")
