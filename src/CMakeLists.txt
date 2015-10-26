### DO NOT DELETE THIS COMMENT: INSERT_ARCHETYPES_HERE ###
USE_CYCLUS("cyborg" "reactor")

INSTALL_CYCLUS_MODULE("cyborg" "/home/jenny/cyclus/cyborg")

# install header files
FILE(GLOB h_files "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
INSTALL(FILES ${h_files} DESTINATION include/cyborg COMPONENT cyborg)
