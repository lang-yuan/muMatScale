#!/bin/bash -l
module load PrgEnv-amd/8.7.0
module load rocm
module load craype-accel-amd-gfx90a
module load cray-hdf5
module load cmake
module load cray-pmi

rm -rf build
mkdir build
cd build

cmake -DCMAKE_C_COMPILER=cc -DBUILD_FOR_GPU=yes \
      -DBUILD_WITH_GPU_AWARE_MPI=yes \
      -DMPIEXEC_EXECUTABLE="/usr/bin/srun" \
      -DMPIEXEC_NUMPROC_FLAG="-n" \
      -DMPIEXEC_PREFLAGS="-c7;--gpus-per-task=1;--gpu-bind=closest" \
      -DCMAKE_C_FLAGS="-fopenmp -fopenmp-assume-no-thread-state -fopenmp-targets=amdgcn-amd-amdhsa -Xopenmp-target=amdgcn-amd-amdhsa -march=gfx90a" \
      -DCMAKE_EXE_LINKER_FLAGS="-lm -fopenmp" \
      -DCMAKE_EXE_LINKER_FLAGS="${PE_MPICH_GTL_DIR_amd_gfx90a} -Wl,--no-as-needed -lmpi_gtl_hsa" \
      -DCMAKE_BUILD_TYPE=Release \
      ..

make -j5

