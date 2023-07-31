/***************************************************************/
/* Copyright (c) 2023, Lang Yuan, Univeristy of South Carolina */
/* All rights reserved.                                        */
/* This file is part of muMatScale.                            */
/* See the top-level LICENSE file for details.                 */
/***************************************************************/
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include "globals.h"
#include "functions.h"
#include "grain.h"
#include "debug.h"
#include "temperature.h"
#include "calculate.h"

#define GAMMA 0                 // 0 CENTRAL DIFFERENCE; 1 UPWIND
#define OMEGA 1.7               //Relaxation parameter for SOR iteration
#define ITERMAX 300             //itermax: maximal number of pressure iterations in one time step
#define EPS 0.00001             //eps : stopping tolerance for pressure iteration
#define TOTALITER 200000
