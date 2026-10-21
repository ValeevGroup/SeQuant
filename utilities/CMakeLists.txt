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

# the viewer page is compiled into sequant-schedule, so `view` needs no file
# beside the binary; --template substitutes another copy
set(schedule_view_html "${CMAKE_CURRENT_SOURCE_DIR}/schedule-view/schedule-view.html")
set(schedule_view_cpp "${CMAKE_CURRENT_BINARY_DIR}/schedule_view_template.cpp")
add_custom_command(
    OUTPUT "${schedule_view_cpp}"
    COMMAND ${CMAKE_COMMAND} -DINPUT=${schedule_view_html} -DOUTPUT=${schedule_view_cpp} -P "${CMAKE_CURRENT_SOURCE_DIR}/schedule-view/embed.cmake"
    DEPENDS "${schedule_view_html}" "${CMAKE_CURRENT_SOURCE_DIR}/schedule-view/embed.cmake"
    COMMENT "Embedding schedule-view.html into sequant-schedule"
    VERBATIM
)
target_sources(sequant-schedule PRIVATE "${schedule_view_cpp}")

add_subdirectory(external-interface)
add_subdirectory(cost_analysis)

add_test(
    NAME "sequant/schedule-view/validate-example"
    COMMAND sequant-schedule validate "${CMAKE_CURRENT_SOURCE_DIR}/schedule-view/example.json"
)
add_test(
    NAME "sequant/schedule-view/view-example"
    COMMAND sequant-schedule view "${CMAKE_CURRENT_SOURCE_DIR}/schedule-view/example.json" -o "${CMAKE_CURRENT_BINARY_DIR}/example.html"
)
