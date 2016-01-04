cmake_minimum_required (VERSION 2.6)

project (cosmos)

set (resources ${PROJECT_SOURCE_DIR}/resources)

# compile resources

set (generated ${PROJECT_BINARY_DIR}/generated)
make_directory (${generated})
message (STATUS "generated : ${generated}")

add_executable (gravity_compiler
	src/util/gravity_compiler.c
)
target_link_libraries (gravity_compiler m)

add_custom_command (
	OUTPUT ${generated}/gravity_data.c
	COMMAND gravity_compiler ${generated}/gravity_data.c ${resources}/gem4/gem4.gfc
	DEPENDS gravity_compiler
)

# comile main library

include_directories (src test src/util ${generated})

add_library (cosmos
	src/gravity.c
	${generated}/gravity_data.c
)

add_dependencies(cosmos gravity_compiler)

# test

enable_testing ()

add_executable (potential44 test/potential44.c)
target_link_libraries (potential44 m cosmos)
add_test (NAME gravity.potential44 COMMAND potential44)
