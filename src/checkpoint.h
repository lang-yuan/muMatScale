/***************************************************************/
/* Copyright (c) 2023, Lang Yuan, Univeristy of South Carolina */
/* All rights reserved.                                        */
/* This file is part of muMatScale.                            */
/* See the top-level LICENSE file for details.                 */
/***************************************************************/

#ifndef __CHECKPOINT_H_
#define __CHECKPOINT_H_
#include <inttypes.h>
#include "hdf5.h"
//#include "/usr/lib/x86_64-linux-gnu/hdf5/serial/include/hdf5.h"

void writeTaskCheckpoint(
    uint64_t timestep);
void writeMainCheckpoint(
    );
void mainRestart(
    uint64_t restart);
void taskRestart(
    uint64_t restart);
void clean_checkpoint(
    );

typedef struct taskSub
{
    int rank;
    hvl_t assigned_sb;

} taskData;


#endif
