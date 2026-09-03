find_package(PythonInterp 3.2 REQUIRED)

include_directories( ${CMAKE_SOURCE_DIR}/src)

add_executable(testHaloFill ${CMAKE_SOURCE_DIR}/tests/HaloFill/main.c
               ${CMAKE_SOURCE_DIR}/src/globals.c
               ${CMAKE_SOURCE_DIR}/src/calculate.c
               ${CMAKE_SOURCE_DIR}/src/growth.c
               ${CMAKE_SOURCE_DIR}/src/grain.c
               ${CMAKE_SOURCE_DIR}/src/face_util.c
               ${CMAKE_SOURCE_DIR}/src/functions.c
               ${CMAKE_SOURCE_DIR}/src/temperature.c
               ${CMAKE_SOURCE_DIR}/src/ll.c
               ${CMAKE_SOURCE_DIR}/src/distribute.c
               ${CMAKE_SOURCE_DIR}/tests/HaloFill/read_ctrl.c
               ${CMAKE_SOURCE_DIR}/src/profiler.c
               ${CMAKE_SOURCE_DIR}/src/stopwatch.c
               ${CMAKE_SOURCE_DIR}/src/packing.c)
target_link_libraries(testHaloFill PRIVATE m)
target_link_libraries(testHaloFill PRIVATE ${HDF5_LIBRARIES})
if(${BUILD_FOR_GPU})
  message(STATUS "Build testHaloFill with OpenMP offload")
  target_compile_definitions(testHaloFill PRIVATE GPU_OMP)
  target_compile_definitions(testHaloFill PRIVATE GPU_PACK)
  target_link_libraries(testHaloFill PRIVATE OpenMP::OpenMP_C)
  #set_target_properties(testHaloFill PROPERTIES COMPILE_FLAGS ${EXTRA_OpenMP_C_FLAGS})
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${EXTRA_OpenMP_C_FLAGS}")
  target_link_libraries(testHaloFill PRIVATE ${EXTRA_OpenMP_C_FLAGS})
  if(NOT ${BUILD_WITH_GPU_AWARE_MPI})
    target_compile_definitions(muMatScale PRIVATE CPU_MPI)
  endif()
endif()

target_include_directories(testHaloFill PRIVATE ${HDF5_INCLUDE_DIRS})
target_compile_definitions(testHaloFill PRIVATE NUC_SEP)

add_test(NAME testHaloFill
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 8 ${MPIEXEC_PREFLAGS}
	 ${CMAKE_CURRENT_BINARY_DIR}/testHaloFill
         ${CMAKE_CURRENT_SOURCE_DIR}/HaloFill/pack.in)

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

add_test(NAME AMtrack
         COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testAMtrack.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 6 ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/../muMatScale
         ${CMAKE_CURRENT_SOURCE_DIR}/am_singletrack.in
         )
