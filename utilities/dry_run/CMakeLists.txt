include(FindOrFetchCLI11)
include(FindOrFetchJSON)

add_executable(dry_run dry_run.cpp)
set_target_properties(dry_run PROPERTIES
    CXX_SCAN_FOR_MODULES OFF
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
target_link_libraries(dry_run PRIVATE
    SeQuant::mbpt nlohmann_json::nlohmann_json CLI11::CLI11)
target_set_compiler_flags(dry_run)
target_set_optimization_flags(dry_run)

# Each example is a driver written by MPQC (`plan_sizes:save`) for a
# Kramers-restricted PNS-CCD particle-particle-ladder residual, with the
# meter report MPQC's own dry run produced recorded in it (and the tree
# fingerprints, which are hash values, taken out). The utility exits with
# status 2 unless its replay reproduces that report.
foreach(NAME IN ITEMS "pns_ccd_ppl" "pns_ccd_ppl_large")
    add_test(
        NAME "sequant/dry_run/${NAME}"
        COMMAND dry_run --driver "${CMAKE_CURRENT_LIST_DIR}/examples/${NAME}.json")
endforeach()
