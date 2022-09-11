project(HopfieldNetwork)
cmake_minimum_required(VERSION 3.0.0)
find_package(OpenMP REQUIRED)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ../build/)
add_executable(example ./example.C)

if(NOT DEFINED N_PARALLEL_THREADS)
	execute_process (COMMAND bash -c "nproc" OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE NPROC)
	set(N_PARALLEL_THREADS ${NPROC})
	message("-- WARNING: No value explicitly passed for N_PARALLEL_THREADS, compiling with N_PARALLEL_THREADS = nproc = " ${N_PARALLEL_THREADS})
else()
	message("-- Compiling with N_PARALLEL_THREADS = " ${N_PARALLEL_THREADS})
endif()
target_compile_options(example PRIVATE ${OpenMP_CXX_FLAGS} -DN_PARALLEL_THREADS=${N_PARALLEL_THREADS} -Ofast -Wall -Wextra )
target_link_options(example PRIVATE ${OpenMP_CXX_FLAGS})

include(CTest)

add_executable(pushback_wrong_size ../test/pushback.test.C)
add_executable(uninitialised_patterns ../test/uninitpatterns.test.C)
add_executable(uninitialised_weights ../test/uninitweights.test.C)

add_test(NAME pushback_wrong_size COMMAND pushback_wrong_size)
add_test(NAME uninitialised_patterns COMMAND uninitialised_patterns)
add_test(NAME uninitialised_weights COMMAND uninitialised_weights)

set_tests_properties(pushback_wrong_size uninitialised_weights uninitialised_patterns PROPERTIES WILL_FAIL True)