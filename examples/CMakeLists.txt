prepare_executable(example_simple   example_simple.c)
prepare_executable(example_data     example_data.c)
prepare_executable(example_graphics example_graphics.c)
prepare_executable(example_spheres  example_spheres.c.c)

add_custom_command(
        TARGET example_simple POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_SOURCE_DIR}/default.toml
                ${CMAKE_CURRENT_BINARY_DIR}/default.toml)
