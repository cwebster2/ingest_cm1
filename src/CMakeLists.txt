set (SOURCES cm1_base.f90 cm1_grads.f90 cm1_grads_mpi.f90 cm1_grads_single.f90 cm1_hdf5.F90 ingest_cm1.F90)

add_library(ingest_cm1 ${SOURCES})

set_target_properties(ingest_cm1
  PROPERTIES
  VERSION 1.0.0
  SOVERSION 1)

target_link_libraries(ingest_cm1 ${LINK_LIBS})

#target_link_libraries(ingest_cm1 hdf-libs)

INSTALL(
  TARGETS ingest_cm1
  EXPORT ingest_cm1
  DESTINATION lib
  PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
)

install(
  EXPORT ingest_cm1
  NAMESPACE ingest_cm1
  DESTINATION lib
)

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/
  DESTINATION include
  FILES_MATCHING PATTERN "*.mod"
  PATTERN "CMakeFiles" EXCLUDE
  PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
)              
