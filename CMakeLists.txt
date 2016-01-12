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

add_custom_command (
	OUTPUT ${generated}/gravity_model_egm96.c
	COMMAND gravity_compiler ${generated}/gravity_model_egm96.c ${resources}/egm96/egm96.gfc
	DEPENDS gravity_compiler
)
add_library (gravity_model_egm96 ${generated}/gravity_model_egm96.c)

# main library

add_library (gravity src/gravity.c)

# tests

enable_testing ()

make_directory (${generated}/test)
set (test_resources ${PROJECT_SOURCE_DIR}/test/resources)
add_custom_command (
	OUTPUT ${generated}/test/test_model.c
	COMMAND gravity_compiler ${generated}/test/test_model.c ${test_resources}/test_model.gfc
	DEPENDS gravity_compiler ${test_resources}/test_model.gfc
)
add_library (test_model ${generated}/test/test_model.c)

add_executable (potential44 test/potential44.c)
target_link_libraries (potential44 m gravity test_model)
add_test (NAME gravity.potential44 COMMAND potential44)

add_executable(gdf.gem4 test/gdf.c)
target_link_libraries (gdf.gem4 m gravity gravity_model_gem4)
add_test (NAME gravity.gdf.potential.gem4 COMMAND gdf.gem4 ${resources}/gem4/gem4-9280.gdf)

add_executable(gdf.egm96 test/gdf.c)
target_link_libraries (gdf.egm96 m gravity gravity_model_egm96)
add_test (NAME gravity.gdf.potential.egm96 COMMAND gdf.egm96 ${resources}/egm96/egm96-9287.gdf)
