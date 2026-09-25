# Compiles and runs C code produced by the TAPP export generator
enable_language(C)

add_executable(export_tapp_generate ${BUILD_BY_DEFAULT} "generate.cpp")
set_target_properties(export_tapp_generate PROPERTIES CXX_SCAN_FOR_MODULES OFF)
target_link_libraries(export_tapp_generate PRIVATE SeQuant::export)
target_set_warning_flags(export_tapp_generate)

set(generated_sources "")
foreach(kind IN ITEMS real complex)
    set(generated "${CMAKE_CURRENT_BINARY_DIR}/tapp_${kind}.c")
    add_custom_command(
        OUTPUT "${generated}"
        COMMAND export_tapp_generate "${generated}" ${kind}
        DEPENDS export_tapp_generate
        COMMENT "Generating TAPP C code (${kind})"
    )
    list(APPEND generated_sources "${generated}")
endforeach()
# The generated files are #included by the driver
set_source_files_properties(${generated_sources} PROPERTIES HEADER_FILE_ONLY ON)

add_executable(export_tapp ${BUILD_BY_DEFAULT} "driver.c" ${generated_sources})
set_target_properties(export_tapp PROPERTIES
    C_STANDARD 11
    C_STANDARD_REQUIRED ON
    C_EXTENSIONS OFF
)
target_include_directories(export_tapp PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")
target_link_libraries(export_tapp PRIVATE tapp::reference)
target_set_warning_flags(export_tapp)

set(test_name "sequant/integration/export_tapp")
add_test(NAME "${test_name}" COMMAND export_tapp)
build_test_as_needed(export_tapp "${test_name}" test_name)
