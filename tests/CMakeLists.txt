###############################################################################
################################### begin cyborg testing ######################
###############################################################################

  FIND_PACKAGE(Threads)
  IF(CMAKE_USE_PTHREADS_INIT)  # The pthreads library is available.
      SET(cxx_base_flags "${cxx_base_flags} -DGTEST_HAS_PTHREAD=1")
  ENDIF()

###############################################################################
# CyBORG unit tests
###############################################################################

  include_directories( ${CMAKE_CURRENT_SOURCE_DIR} "${PROJECT_SOURCE_DIR}/src" ${CYBORG_INCLUDE_DIRS})

  # Build cyborg_unit_tests
  ADD_EXECUTABLE(cyborg_unit_tests
      cyborg_unit_test_driver.cc
      ${TestSource}
  )
  MESSAGE("Cyborg unit tests: ${TestSource}")    
  TARGET_LINK_LIBRARIES(cyborg_unit_tests
      dl ${LIBS} ${CYBORG_LIB} ${CYCLUS_TEST_LIBRARIES} )

  INSTALL(TARGETS cyborg_unit_tests
      RUNTIME DESTINATION bin
      COMPONENT testing
  )

###############################################################################
# OrigenInterface unit tests
###############################################################################

  # Add source for OrigenInterface tests
  FILE(GLOB OrigenInterface_TEST_CC "${CMAKE_CURRENT_SOURCE_DIR}/*origen*.cpp")
   
  # Build cyborg_unit_tests
  ADD_EXECUTABLE(OrigenInterface_unit_tests
      cyborg_unit_test_driver.cc
      ${OrigenInterface_TEST_CC} )

  TARGET_INCLUDE_DIRECTORIES(OrigenInterface_unit_tests PUBLIC ${SCALE_INCLUDE_DIRS} ) 

  TARGET_LINK_LIBRARIES(OrigenInterface_unit_tests
      dl OrigenInterface ${LIBS} ${CYCLUS_TEST_LIBRARIES} )

  INSTALL(TARGETS OrigenInterface_unit_tests
      RUNTIME DESTINATION bin
      COMPONENT testing
  )

  # Add OrigenInterface test files back to TestSource
  SET(TestSource ${OrigenInterface_TEST_CC} PARENT_SCOPE)

###############################################################################
#################################### end cyborg testing #######################
###############################################################################
