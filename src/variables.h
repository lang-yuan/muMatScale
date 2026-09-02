/***************************************************************/
/* Copyright (c) 2023, Lang Yuan, Univeristy of South Carolina */
/* All rights reserved.                                        */
/* This file is part of muMatScale.                            */
/* See the top-level LICENSE file for details.                 */
/***************************************************************/

#ifndef __VARIABLES_H__
#define __VARIABLES_H__

#include <mpi.h>

int registerCommInfo(size_t datasize);
void FinishExchangeForVar(int variable_key, void* data);
void ExchangeFacesForVar(int variable_key, void* d);

typedef struct variable_registration
{
    MPI_Request *reqs;
    int nreq;
    size_t datasize;
    // pointers to halo exchange buffer allocations
    char *rbuf_base;
    char *sbuf_base;
    // max. number of elements in halo exchange buffer
    int buffer_slot_cells;
    // max. size in bytes for 6 halo exchange buffers
    size_t buffer_slot_bytes;
    void *rbuf[6];
    void *sbuf[6];
} variable_registration;

#endif
