cmake_minimum_required(VERSION 3.16)
project(muMatScale C)

set(HDF5_STATIC FALSE CACHE BOOL "Whether to link with static hdf5")
if(${HDF5_STATIC})
  set(HDF5_USE_STATIC_LIBRARIES TRUE)
endif()

find_package(MPI REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS C)
find_package(ZLIB REQUIRED)
find_package(OpenMP REQUIRED)

if (${HDF5_FOUND})
  message(STATUS "HDF5_LIBRARIES after search: ${HDF5_LIBRARIES}")
  message(STATUS "HDF5_INCLUDE_DIRS: ${HDF5_INCLUDE_DIRS}")
else(${HDF5_FOUND})
  message(FATAL_ERROR "Required HDF5 package not found.")
endif (${HDF5_FOUND})

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
  message(STATUS "No user specified build type, using default...")
endif()
message(STATUS "Build type is ${CMAKE_BUILD_TYPE}")

set(src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/packing.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/calculate.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ca_main.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/checkpoint.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/curvature.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/distribute.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/face_util.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/file_io.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/globals.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/growth.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/functions.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/grain.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/io_vtk.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/io_xdmf.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ll.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/loop.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/profiler.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/stopwatch.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/read_ctrl.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/setup.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/temperature.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/vtkCompress.c
  )

add_executable(muMatScale ${src})

set_property(TARGET muMatScale PROPERTY C_STANDARD 99)

target_link_libraries(muMatScale PRIVATE hdf5::hdf5)
target_link_libraries(muMatScale PRIVATE MPI::MPI_C)
target_link_libraries(muMatScale PRIVATE ZLIB::ZLIB)
target_link_libraries(muMatScale PRIVATE OpenMP::OpenMP_C)
target_link_libraries(muMatScale PRIVATE m)

target_include_directories(muMatScale PRIVATE ${HDF5_INCLUDE_DIRS})

target_compile_definitions(muMatScale PRIVATE INDEX_SEP)
target_compile_definitions(muMatScale PRIVATE NUC_SEP)
target_compile_definitions(muMatScale PRIVATE EBUG=0)

set(BUILD_FOR_GPU FALSE CACHE BOOL "Whether to build for GPU")
if(${BUILD_FOR_GPU})
    if(OpenMP_C_VERSION LESS "4.5")
        message(SEND_ERROR "Requires OpenMP4.5 or newer for GPU offload")
    endif()
    message(STATUS "Build with OpenMP offload")
    set(EXTRA_OpenMP_C_FLAGS "" CACHE STRING "Extra flags for OpenMP")
    if(CMAKE_C_COMPILER_ID STREQUAL "XLClang")
        set(EXTRA_OpenMP_C_FLAGS "-qsmp=omp -qoffload")
    endif()
    message(STATUS "OpenMP_C_FLAGS: ${OpenMP_C_FLAGS}")
    message(STATUS "EXTRA_OpenMP_C_FLAGS: ${EXTRA_OpenMP_C_FLAGS}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${EXTRA_OpenMP_C_FLAGS}")
    target_compile_definitions(muMatScale PRIVATE GPU_OMP)
    target_compile_definitions(muMatScale PRIVATE GPU_PACK)
    target_compile_definitions(muMatScale PRIVATE GPU_OMP_NUC)
    target_link_libraries(muMatScale PRIVATE ${EXTRA_OpenMP_C_FLAGS})
endif()

enable_testing()

add_subdirectory(tests)

install(TARGETS muMatScale DESTINATION bin)
