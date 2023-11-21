/***************************************************************/
/* Copyright (c) 2023, Lang Yuan, Univeristy of South Carolina */
/* All rights reserved.                                        */
/* This file is part of muMatScale.                            */
/* See the top-level LICENSE file for details.                 */
/***************************************************************/

#include <sys/types.h>
#include <sys/time.h>
#include <mpi.h>
#ifndef __USE_POSIX             // For 'fileno' in stdio.h
#define __USE_POSIX
#endif
#include <values.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <inttypes.h>
#include <time.h>
#include <math.h>

#include "distribute.h"
#include "functions.h"
#include "globals.h"
#include "xmalloc.h"
#include "debug.h"
#include "profiler.h"
#include "ll.h"

extern SB_struct *lsp;

// list of all subblocks
ll_t *gsubblock_list;

// for MPI sned/recv requests
MPI_Request *request;
int nrequests;

// Initialize the subblock array
void
init_gmsp(
    )
{
    size_t gnsb = bp->gnsbx * bp->gnsby * bp->gnsbz;

    // Allocate array of empty subblocks
    printf("Initializing %lu Subblocks...\n", gnsb);
    xmalloc(gmsp, MSB_struct, gnsb);

    for (int id = 0; id < (bp->gnsbx * bp->gnsby * bp->gnsbz); id++)
    {
        gmsp[id].subblockid = id;
        gmsp[id].nnuc = 0;
    }

    for (int id = 0; id < (bp->gnsbx * bp->gnsby * bp->gnsbz); id++)
        ll_append(gsubblock_list, &gmsp[id]);
}


void
sendSubblockInfo(
    )
{
    // Send the subblocks to each processor
    int count = ll_count(gsubblock_list);
    request = malloc(2 * count * sizeof(MPI_Request));

    // count number of requests sent
    nrequests = 0;

    lli_t *llp = gsubblock_list->head;
    while (llp != NULL)
    {
        MSB_struct *msb = (MSB_struct *) llp->data;

        timing(COMPUTATION, timer_elapsed());

        MPI_Isend(msb, sizeof(MSB_struct), MPI_BYTE, msb->subblockid, SYNCTAG,
                  mpi_comm_new, request + nrequests);
        nrequests++;
        if (msb->nnuc)
        {
            MPI_Isend(msb->nuc_pt, sizeof(nucleation_t) * msb->nnuc, MPI_BYTE,
                      msb->subblockid, SYNCTAG, mpi_comm_new,
                      request + nrequests);
            nrequests++;
        }
        timing(COMMUNICATION, timer_elapsed());

        llp = llp->next;
    }
}

void
completeSendSubblockInfo(
    )
{
    MPI_Waitall(nrequests, request, MPI_STATUSES_IGNORE);
    free(request);
}

/**
 * Receive subblock info from task 0
 */
void
recvSubblockInfo(
    )
{
    timing(COMPUTATION, timer_elapsed());

    MSB_struct msb;
    MPI_Status status;
    MPI_Recv(&msb, sizeof(MSB_struct), MPI_BYTE, 0, SYNCTAG,
             mpi_comm_new, &status);

    timing(COMMUNICATION, timer_elapsed());

    // allocate structure for subblock data
    xmalloc(lsp, SB_struct, 1);

    CopyMSBtoSB(&msb, lsp);

    if (lsp->nnuc)
    {
        xmalloc(lsp->nuc_pt, nucleation_t, lsp->nnuc);

        timing(COMPUTATION, timer_elapsed());
        MPI_Recv(lsp->nuc_pt, sizeof(nucleation_t) * lsp->nnuc, MPI_BYTE,
                 0, SYNCTAG, mpi_comm_new, &status);
        timing(COMMUNICATION, timer_elapsed());
    }
}

/**
 * Frees memory associated with subblocks
 */
void
clearSubblocks(
    )
{
    SB_struct *s = lsp;
    if (s->mold)
        xfree(s->mold);
    if (s->fs)
        xfree(s->fs);
    if (s->nuc_threshold)
        xfree(s->nuc_threshold);
    if (s->ce)
        xfree(s->ce);
    if (s->oce)
        xfree(s->oce);
    if (s->cl)
        xfree(s->cl);
    if (s->diff_id)
        xfree(s->diff_id);
    if (s->fs_id)
        xfree(s->fs_id);
    if (s->nuc_id)
        xfree(s->nuc_id);
    if (s->nuc_id2)
        xfree(s->nuc_id2);
    if (s->dc)
        xfree(s->dc);
    if (s->d)
        xfree(s->d);
    if (s->gr)
        xfree(s->gr);
    if (s->ogr)
        xfree(s->ogr);
    if (s->temperature)
        xfree(s->temperature);
    //ADD_CURV_LY
    if (bp->tip_curv == 1 || s->curv)
        xfree(s->curv);

    //Fluid flow
    if (s->cell_u)
        xfree(s->cell_u);
    if (s->cell_v)
        xfree(s->cell_v);
    if (s->cell_w)
        xfree(s->cell_w);

    // AM
    if (s->lsindex)
        xfree(s->lsindex);
}
