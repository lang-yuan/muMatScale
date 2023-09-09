find_package(PythonInterp 3.2 REQUIRED)

include_directories( ${CMAKE_SOURCE_DIR}/src)

add_executable(testPacking ${CMAKE_SOURCE_DIR}/tests/Packing/main.c
               ${CMAKE_SOURCE_DIR}/src/globals.c
               ${CMAKE_SOURCE_DIR}/tests/Packing/functions.c
               ${CMAKE_SOURCE_DIR}/tests/Packing/read_ctrl.c
               ${CMAKE_SOURCE_DIR}/src/packing.c)
target_link_libraries(testPacking PRIVATE m)
if(${BUILD_FOR_GPU})
  target_compile_definitions(testPacking PRIVATE GPU_PACK)
  target_link_libraries(testPacking PRIVATE OpenMP::OpenMP_C)
  #set_target_properties(testPacking PROPERTIES COMPILE_FLAGS ${EXTRA_OpenMP_C_FLAGS})
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${EXTRA_OpenMP_C_FLAGS}")
  target_link_libraries(testPacking PRIVATE ${EXTRA_OpenMP_C_FLAGS})
endif()

add_test(NAME testPacking
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 6 ${MPIEXEC_PREFLAGS}
	 ${CMAKE_CURRENT_BINARY_DIR}/testPacking
         ${CMAKE_CURRENT_SOURCE_DIR}/Packing/pack.in)

#regression tests
add_test(NAME OneGrain
         COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testOneGrain.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 6 ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/../muMatScale
         ${CMAKE_CURRENT_SOURCE_DIR}/one_grain.in
         )
add_test(NAME MultiGrains
         COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testMultiGrains.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 6 ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/../muMatScale
         ${CMAKE_CURRENT_SOURCE_DIR}/ca6mpi.in
         )
