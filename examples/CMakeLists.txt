macro(prepare_executable target source)
    add_executable(${target} ${source})
    target_link_libraries(${target} m cblas gsl pthread model_static)
    target_include_directories(${target} PUBLIC ${MODEL_STATIC_INCLUDE_DIR})
endmacro(prepare_executable)

prepare_executable(example_simple   example_simple.c)
prepare_executable(example_data     example_data.c)
prepare_executable(example_graphics example_graphics.c)

add_custom_command(
        TARGET example_simple POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_SOURCE_DIR}/default.toml
                ${CMAKE_CURRENT_BINARY_DIR}/default.toml)
