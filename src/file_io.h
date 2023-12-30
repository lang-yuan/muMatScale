/***************************************************************/
/* Copyright (c) 2023, Lang Yuan, Univeristy of South Carolina */
/* All rights reserved.                                        */
/* This file is part of muMatScale.                            */
/* See the top-level LICENSE file for details.                 */
/***************************************************************/

#ifndef __FILE_IO_H__
#define __FILE_IO_H__

#include <inttypes.h>

#include "globals.h"


/**
 * Returns the name of the type of IO that will be done
 * \return a null-terminated string with the type of File-IO to be done
 */
const char *getIOTypeName(
    );

/**
 * Prepares the IO subsystem
 * Called by both Main and Tasks
 */
void prepareIO(
    );

/**
 * Writes out any data that the Main needs to write
 */
void writeMain(
    int timestep);

/**
 * Prepares the IO subsystem for a new Subblock
 */
void prepareIOSubblock(
    SB_struct * sb);

/**
 * Writes out any data that the Tasks need to write for their subblocks
 */
void writeSubblocks(
    int timestep);

/**
 * Finialize any IO for the simulation.
 */
void closeIO(
    );

#endif
