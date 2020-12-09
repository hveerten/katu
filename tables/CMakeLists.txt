macro(prepare_executable target source)
    add_executable(${target} ${source})
    target_link_libraries(${target} m gsl gslcblas pthread model_static)
    target_include_directories(${target} PUBLIC ${MODEL_STATIC_INCLUDE_DIR})
endmacro(prepare_executable)

prepare_executable(generate_table_BH_ka generate_table_BH_ka.c)
