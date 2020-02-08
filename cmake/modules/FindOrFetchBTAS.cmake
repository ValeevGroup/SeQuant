find_package(BTAS 1.0.0 QUIET)

if (NOT TARGET BTAS::btas)

  include(DownloadProject)
  download_project(PROJ                BTAS
    GIT_REPOSITORY      https://github.com/BTAS/btas.git
    GIT_TAG             ${SEQUANT_TRACKED_BTAS_TAG}
    UPDATE_DISCONNECTED 1
    )

  add_subdirectory(${BTAS_SOURCE_DIR} ${BTAS_BINARY_DIR})

endif(NOT TARGET BTAS::btas)