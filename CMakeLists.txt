CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
# taken from http://geant4.cern.ch/support/source/geant4/CMakeLists.txt
IF(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(STATUS "Stub requires an out-of-source build.")
    message(STATUS "Please remove these files from ${CMAKE_BINARY_DIR} first:")
    message(STATUS "CMakeCache.txt")
    message(STATUS "CMakeFiles")
    message(STATUS "Once these files are removed, create a separate directory")
    message(STATUS "and run CMake from there")
    message(FATAL_ERROR "in-source build detected")
ENDIF()

PROJECT(CYBORG)

# quiets fortify_source warnings when not compiling with optimizations
# in linux distros where compilers were compiled with fortify_source enabled by
# default (e.g. Arch linux).

# ORIGEN requires C++11
SET(CMAKE_CXX_STANDARD 11)
SET(CXX_STANDARD_REQUIRED ON )



MESSAGE("-- CMake Build type: ${CMAKE_BUILD_TYPE}")
STRING(TOLOWER "${CMAKE_BUILD_TYPE}" BUILD_TYPE)
IF(NOT BUILD_TYPE STREQUAL "release")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -U_FORTIFY_SOURCE -fno-stack-protector -g -O0")
    # embed entire source code into binaries for debug builds
    IF(BUILD_TYPE STREQUAL "debug")
       SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g3")
    ENDIF()
ENDIF()

# suppress warnings regarding deprecated #_BSD_SOURCE and #_SVID_SOURCE 
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_DEFAULT_SOURCE")

# no overflow warnings because of silly coin-ness
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-overflow")

# Direct any out-of-source builds to this directory
SET(SOURCE_DIR ${CMAKE_SOURCE_DIR})

IF(NOT CYCLUS_DOC_ONLY)
    # try excluding external symbols from libs?
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,--no-export-dynamic")

    # Direct any binary installation paths to this directory
    SET(BINARY_DIR ${CMAKE_BINARY_DIR})

    # This makes all the libraries build as SHARED
    SET(BUILD_SHARED_LIBS true)

    # Setup build locations.
    IF(NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
        SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BINARY_DIR}/bin)
    endif()
    IF(NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
        SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${BINARY_DIR}/lib)
    endif()
    IF(NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
        SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${BINARY_DIR}/lib)
    ENDIF()

    SET(CYBORG_EXECUTABLE_DIR ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

    # use, i.e. don't skip the full RPATH for the build tree
    SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

    # when building, don't use the install RPATH already
    # (but later on when installing)
    SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
    SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

    # add the automatically determined parts of the RPATH
    # which point to directories outside the build tree to the install RPATH
    SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

    # the RPATH to be used when installing, but only if it's not a system directory
    LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
    IF("${isSystemDir}" STREQUAL "-1")
       SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
    ENDIF("${isSystemDir}" STREQUAL "-1")

    # Tell CMake to use the library type it can find (debug/release)
    INCLUDE(SelectLibraryConfigurations)

    # Tell CMake where the modules are
    # Use the "new" CMake behavior where user-provided modules are preferred
    cmake_policy(SET CMP0017 NEW)
    SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_ROOT}/share/Modules")
    SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CYBORG_SOURCE_DIR}/cmake)
    SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CYCLUS_CORE_SHARE_DIR}/cmake)

    # Find cyclus
    FIND_PACKAGE(Cyclus REQUIRED)
    SET(CYBORG_INCLUDE_DIRS ${CYBORG_INCLUDE_DIRS} ${CYCLUS_CORE_INCLUDE_DIR})
    SET(LIBS ${LIBS} ${CYCLUS_CORE_LIBRARIES})

    # Include macros
    INCLUDE(UseCyclus)

    # Find LibXML++ and dependencies
    FIND_PACKAGE(LibXML++ REQUIRED)
    SET(CYBORG_INCLUDE_DIRS ${CYBORG_INCLUDE_DIRS} ${LibXML++_INCLUDE_DIR} ${Glibmm_INCLUDE_DIRS} ${LibXML++Config_INCLUDE_DIR})
    SET(LIBS ${LIBS} ${LibXML++_LIBRARIES})

    MESSAGE("--LIBS: ${LIBS}")
    MESSAGE("--LD_LIBRARY_PATH: $ENV{LD_LIBRARY_PATH}")

    # Include the boost header files, system, and filesystem libraries
    SET(Boost_USE_STATIC_LIBS       OFF)
    SET(Boost_USE_STATIC_RUNTIME    OFF)
    FIND_PACKAGE(Boost COMPONENTS filesystem system REQUIRED)
    SET(CYBORG_INCLUDE_DIRS ${CYBORG_INCLUDE_DIRS} ${Boost_INCLUDE_DIR})
    SET(LIBS ${LIBS} ${Boost_FILESYSTEM_LIBRARY})
    SET(LIBS ${LIBS} ${Boost_SYSTEM_LIBRARY})

    # find lapack and link to it
    # note there is no include directory variable:
    # http://www.cmake.org/cmake/help/v3.0/module/FindLAPACK.html
    FIND_PACKAGE(LAPACK REQUIRED)
    set(LIBS ${LIBS} ${LAPACK_LIBRARIES})
    MESSAGE("\tFound LAPACK Libraries: ${LAPACK_LIBRARIES}")

    # Find HDF5
    FIND_PACKAGE(HDF5 REQUIRED)
    ADD_DEFINITIONS(${HDF5_DEFINITIONS})
    SET(CYBORG_INCLUDE_DIRS ${CYBORG_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})
    set(LIBS ${LIBS} ${HDF5_LIBRARIES})
    MESSAGE("--    HDF5 Root: ${HDF5_ROOT}")
    MESSAGE("--    HDF5 Include directory: ${HDF5_INCLUDE_DIRS}")
    MESSAGE("--    HDF5 Library directories: ${HDF5_LIBRARY_DIRS}")
    MESSAGE("--    HDF5 Libraries: ${HDF5_LIBRARIES}")

    # find coin and link to it
#    FIND_PACKAGE(COIN REQUIRED)
#    SET(CYBORG_INCLUDE_DIRS ${CYBORG_INCLUDE_DIRS} ${COIN_INCLUDE_DIRS})
#    set(LIBS ${LIBS} ${COIN_LIBRARIES})
#    MESSAGE("--    COIN Root: ${COIN_ROOT}")
#    MESSAGE("--    COIN Include directories: ${COIN_INCLUDE_DIRS}")
#    MESSAGE("--    COIN Libraries: ${COIN_LIBRARIES}")

    FIND_PACKAGE( Sqlite3 REQUIRED )
    SET(CYBORG_INCLUDE_DIRS ${CYBORG_INCLUDE_DIRS} ${SQLITE3_INCLUDE_DIR})
    SET(LIBS ${LIBS} ${SQLITE3_LIBRARIES})
    MESSAGE("--    SQLITE3 Include directories: ${SQLITE3_INCLUDE_DIR}")
    MESSAGE("--    SQLITE3 Libraries: ${SQLITE3_LIBRARIES}")

    # Find Qt libraries (required for ORIGEN)
    # Note that ${QT_USE_FILE} sets up required libraries for Qt
    # See https://cmake.org/cmake/help/cmake2.6docs.html#module%3aFindQt4
    FIND_PACKAGE( Qt COMPONENTS QtCore REQUIRED )
    include(${QT_USE_FILE})
    SET(CYBORG_INCLUDE_DIRS ${CYBORG_INCLUDE_DIRS} ${QT_INCLUDE_DIR} )

    # find SCALE libraries (for ORIGEN)
    # Link ONLY as necessary (i.e., OrigenInterface) to avoid gtest conflict via Nemesis
    FIND_PACKAGE(SCALE REQUIRED)
    SET(CYBORG_INCLUDE_DIRS ${CYBORG_INCLUDE_DIRS} ${SCALE_INCLUDE_DIRS})

    # include all the directories we just found
    INCLUDE_DIRECTORIES(${CYBORG_INCLUDE_DIRS})

    # ------------------------- Add the Agents -----------------------------------
    ADD_SUBDIRECTORY(src)

    # ------------------------- Google Test -----------------------------------

    # Be sure to clear these each time
    SET(GENERATED_TEST_LIST "" CACHE INTERNAL "")
    SET(GENERATED_FILTER_LIST "" CACHE INTERNAL "")
    SET(GENERATED_TEST_SOURCE_LIST "" CACHE INTERNAL "")

    # add tests for the ORIGEN interface
    OPTION(USE_TESTING "Build testing" ON)
    IF(USE_TESTING)
        enable_testing()
        INCLUDE(CTest)
    ENDIF()

    SET(TestSource
        ${TestSource}
        ${CYCLUS_TEST_CORE} )
    LIST(REMOVE_DUPLICATES TestSource)

    SET(ROOT_DIR ${CYCLUS_ROOT_DIR})

    ADD_SUBDIRECTORY(tests)
    #ADD_TEST( NAME tstOrigenInterface COMMAND tstCyclusOrigenInterface )
    ADD_TEST( NAME testCyborg COMMAND cyborg_unit_tests )

    SET(PROJ_DIR ${PROJECT_BINARY_DIR})
    SET(CORE_SHARE_DIR ${CYCLUS_CORE_SHARE_DIR})

ENDIF(NOT CYCLUS_DOC_ONLY)

# This is the directory that holds the doxygen doxyfile template (doxy.conf.in)
SET(DOC_INPUT_DIR ${CYBORG_SOURCE_DIR}/doc)

# If doxygen exists, use the doc/CMakeLists.txt file for further instructions.
FIND_PACKAGE(Doxygen)
IF(DOXYGEN_FOUND)
    #ADD_SUBDIRECTORY(doc)
    #SET(DOC_OUTPUT_DIR ${CMAKE_BINARY_DIR}/doc)
ELSE(DOXYGEN_FOUND)
    MESSAGE(STATUS "WARNING: Doxygen not found - doc won't be created")
ENDIF(DOXYGEN_FOUND)

# uninstall target
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
