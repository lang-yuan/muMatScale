#!/bin/csh
module load PrgEnv-amd/8.4.0
module load craype-accel-amd-gfx90a
module load hdf5
module load cmake

rm -rf build
mkdir build
cd build

cmake -DCMAKE_C_COMPILER=cc -DBUILD_FOR_GPU=yes \
      -DMPIEXEC_EXECUTABLE="/usr/bin/srun" \
      -DMPIEXEC_NUMPROC_FLAG="-n" \
      -DMPIEXEC_PREFLAGS="-c7;--gpus-per-task=1;--gpu-bind=closest" \
      -DCMAKE_C_FLAGS="-fopenmp -fopenmp-assume-no-thread-state -fopenmp-targets=amdgcn-amd-amdhsa -Xopenmp-target=amdgcn-amd-amdhsa -march=gfx90a" \
      -DCMAKE_EXE_LINKER_FLAGS="-lm -fopenmp" \
      -DCMAKE_BUILD_TYPE=Release \
      ..

make

