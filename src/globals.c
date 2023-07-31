/***************************************************************/
/* Copyright (c) 2023, Lang Yuan, Univeristy of South Carolina */
/* All rights reserved.                                        */
/* This file is part of muMatScale.                            */
/* See the top-level LICENSE file for details.                 */
/***************************************************************/
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "globals.h"

/** \file
 * Contains definitions of global data structures.
 */

int iproc;                      ///< My process rank
int nproc;                      ///< Number of processors

MPI_Comm mpi_comm_new;

MSB_struct *gmsp;               ///< list of all subblocks

#ifdef GPU_OMP
#pragma omp declare target
#endif

BB_struct *bp;                  ///< Main configuration structure

#ifdef GPU_OMP
#pragma omp end declare target
#endif


/// Lookup table to "reverse" a face (top->bottom, etc)
const int gface_reverse[NUM_NEIGHBORS] = {
    FACE_BOTTOM, FACE_TOP,
    FACE_RIGHT, FACE_LEFT,
    FACE_BACK, FACE_FRONT
};


int gMPIThreadLevel = 0;
MPI_Datatype MPI_Type_Neighbor;

/**
 * Copies data from a MSB_struct to an SB_struct as part
 * of the setup phase of a new subblock.
 * \param[in] m  Incoming \c MSB_struct
 * \param[out] s Outgoing \c SB_struct
 */
void
CopyMSBtoSB(
    MSB_struct * m,
    SB_struct * s)
{
    s->procid = m->procid;
    s->subblockid = m->subblockid;
    s->coords = m->coords;
    s->nnuc = m->nnuc;
    memcpy(s->neighbors, m->neighbors,
           NUM_NEIGHBORS * NBOR_INFO * sizeof(int));
}
