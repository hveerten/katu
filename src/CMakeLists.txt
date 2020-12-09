add_library(
    model_static
    STATIC
    external_libs/toml.c
    utils/tpool.c
	acceleration.c
	bethe_heitler.c
	bethe_heitler_kelner_aharonian.c
	config.c
	distribution.c
	decay_and_escape.c
	escape.c
	inverse_compton.c
	muon_decay.c
	pair_production.c
	pair_annihilation.c
	pion_decay.c
	pion_production.c
	population.c
	state.c
	state_step.c
	state_step_common.c
	synchrotron.c
    )

target_include_directories(model_static PUBLIC "./external_libs")
target_include_directories(model_static PUBLIC "./utils")
target_include_directories(model_static PUBLIC ".")

target_link_libraries(model_static m)
target_link_libraries(model_static cblas)
target_link_libraries(model_static gsl)
target_link_libraries(model_static pthread)

set_property(TARGET model_static PROPERTY INTERPROCEDURAL_OPTIMIZATION True)

set(MODEL_STATIC_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
