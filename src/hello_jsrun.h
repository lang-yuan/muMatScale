/***************************************************************/
/* Copyright (c) 2023, Lang Yuan, Univeristy of South Carolina */
/* All rights reserved.                                        */
/* This file is part of muMatScale.                            */
/* See the top-level LICENSE file for details.                 */
/***************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <sched.h>
#include <nvml.h>
#include <omp.h>
#include <cuda.h>
#include <cuda_runtime.h>

int hello_jsrun(
    int size,
    int rank);
