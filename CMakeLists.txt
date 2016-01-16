cmake_minimum_required (VERSION 2.6)

project (cosmos)

set (resources ${PROJECT_SOURCE_DIR}/resources)
set (test      ${PROJECT_SOURCE_DIR}/test)
set (generated ${PROJECT_BINARY_DIR}/generated)
make_directory (${generated})

include_directories (src test src/util ${generated})

# comile resources compiler

add_executable (gravity_compiler src/util/gravity_compiler.c)
target_link_libraries (gravity_compiler m)

# compile gem4

add_custom_command (
	OUTPUT ${generated}/gravity_model_gem4.c
	COMMAND gravity_compiler 
		${generated}/gravity_model_gem4.c 
		${resources}/gem4/gem4.gfc
	DEPENDS gravity_compiler
)
add_library (gravity_model_gem4 ${generated}/gravity_model_gem4.c)

# compile egm96

add_custom_command (
	OUTPUT ${generated}/gravity_model_egm96.c
	COMMAND gravity_compiler 
		${generated}/gravity_model_egm96.c 
		${resources}/egm96/egm96.gfc
	DEPENDS gravity_compiler
)
add_library (gravity_model_egm96 ${generated}/gravity_model_egm96.c)

# main library

add_library (gravity src/gravity.c)

# tests

enable_testing ()
set(gravity.potential44 ON)
set(gravity.gdf.gem4    ON)
#set(gravity.gdf.egm96   ON)

make_directory (${generated}/test)
make_directory (${generated}/test/resources)

if(gravity.potential44)
	file(READ ${test}/resources/constants.txt gravity.constants)
	foreach(c ${gravity.constants})
		string(REPLACE "," ";" cc ${c})
		list(GET cc 0 variant)
		list(GET cc 1 earth_gravity_constant)
		list(GET cc 2 radius)
		foreach(m RANGE 0 4)
			foreach (n RANGE ${m} 4)
				set (test_model test_model_${variant}_${n}_${m})
				configure_file (${test}/resources/test_model.gfc
					${generated}/test/resources/${test_model}.gfc)
				add_custom_command (
					OUTPUT ${generated}/test/${test_model}.c
					COMMAND gravity_compiler 
						${generated}/test/${test_model}.c 
						${generated}/test/resources/${test_model}.gfc
					DEPENDS gravity_compiler 
						${generated}/test/resources/${test_model}.gfc
						${test}/resources/test_model.gfc
						${test}/resources/constants.txt)
				add_library (${test_model} ${generated}/test/${test_model}.c)

				set (potential44 potential44_${variant}_${n}_${m})
				add_executable (${potential44} test/potential44.c)
				target_link_libraries (${potential44} m gravity ${test_model})
				add_test (NAME gravity.${potential44} 
					COMMAND ${potential44} 1e-10 33 33)
			endforeach(n)
		endforeach(m)
	endforeach(c)
endif(gravity.potential44)

if(gravity.gdf.gem4)
add_executable(gdf.gem4 test/gdf.c)
	target_link_libraries (gdf.gem4 m gravity gravity_model_gem4)
	add_test (NAME gdf.gem4 
		COMMAND gdf.gem4 ${resources}/gem4/gem4-9280.gdf)
endif(gravity.gdf.gem4)

if(gravity.gdf.egm96)
	add_executable(gdf.egm96 test/gdf.c)
	target_link_libraries (gdf.egm96 m gravity gravity_model_egm96)
	add_test (NAME gdf.egm96 
		COMMAND gdf.egm96 ${resources}/egm96/egm96-9287.gdf)
endif(gravity.gdf.egm96)

