set(EXECUTABLE_SOURCES
    "tensor_network_graphs.cpp"
    "sequant-schedule.cpp"
)

foreach(current IN LISTS EXECUTABLE_SOURCES)
    string(REGEX REPLACE "\\.cpp$" "" name "${current}")
    add_executable(${name} ${current})
    set_target_properties(${name} PROPERTIES CXX_SCAN_FOR_MODULES OFF)
    target_link_libraries(${name} SeQuant)
    
    target_set_compiler_flags("${name}")
    target_set_optimization_flags("${name}")
endforeach()

target_compile_definitions(sequant-schedule PRIVATE SEQUANT_SCHEDULE_VIEW_TEMPLATE="${CMAKE_CURRENT_SOURCE_DIR}/schedule-view/schedule-view.html")

add_subdirectory(external-interface)
add_subdirectory(cost_analysis)
