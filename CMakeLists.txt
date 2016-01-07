cmake_minimum_required (VERSION 2.6)

project (cosmos)

set (resources ${PROJECT_SOURCE_DIR}/resources)

# resource compiler

set (generated ${PROJECT_BINARY_DIR}/generated)
make_directory (${generated})
message (STATUS "generated : ${generated}")
add_executable (gravity_compiler
	src/util/gravity_compiler.c
)
target_link_libraries (gravity_compiler m)
include_directories (src test src/util ${generated})

# gem4

add_custom_command (
	OUTPUT ${generated}/gravity_model_gem4.c
	COMMAND gravity_compiler ${generated}/gravity_model_gem4.c ${resources}/gem4/gem4.gfc
	DEPENDS gravity_compiler
)
add_library (gravity_model_gem4 ${generated}/gravity_model_gem4.c)
add_dependencies (gravity_model_gem4 gravity_compiler)

# main library

add_library (gravity src/gravity.c)

# tests

enable_testing ()

add_executable (potential44 test/potential44.c)
target_link_libraries (potential44 m gravity gravity_model_gem4)
add_test (NAME gravity.potential44 COMMAND potential44)

add_executable (gdf test/gdf.c)
target_link_libraries (gdf m gravity gravity_model_gem4)
add_test (NAME gravity.gdf.potential.gem4 COMMAND gdf ${resources}/gem4/gem4-9280.gdf)

