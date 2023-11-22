/***************************************************************/
/* Copyright (c) 2023, Lang Yuan, Univeristy of South Carolina */
/* All rights reserved.                                        */
/* This file is part of muMatScale.                            */
/* See the top-level LICENSE file for details.                 */
/***************************************************************/

#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <stddef.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "globals.h"
#include "functions.h"
#include "setup.h"
#include "xmalloc.h"
#include "debug.h"
#include "growth.h"
#include "temperature.h"
#include "grain.h"
#include "checkpoint.h"
#include "profiler.h"
#include "face_util.h"
#include "calculate.h"
#include "packing.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef USE_MPE
#include <mpe.h>
extern int gLogStateBeginID;
extern int gLogStateEndID;
#endif

#include "curvature.h"

SB_struct *lsp;


typedef struct variable_registration
{
    MPI_Request *reqs;
    int nreq;
    size_t datasize;
    void *rbuf[6];
    void *sbuf[6];
} variable_registration;

static int var_count = 0;
static variable_registration var_regs[1 << TAG_DATA_KEY_SHIFT];

static int
registerCommInfo(
    size_t datasize)
{
    assert(var_count < (1 << TAG_DATA_KEY_SHIFT));
    int varnum = var_count++;

    variable_registration *v = &var_regs[varnum];
    v->reqs = NULL;
    v->nreq = 0;
    v->datasize = datasize;
    for (int i = 0; i < 6; i++)
    {
        v->rbuf[i] = NULL;
        v->sbuf[i] = NULL;
    }
    return varnum;
}


static void
FinishExchangeForVar(
    int variable_key, void* data)
{
    variable_registration *v = &var_regs[variable_key];

    dwrite(DEBUG_MPI, "%d: Just before Waitall\n", iproc);
    timing(COMPUTATION, timer_elapsed());

    MPI_Waitall(v->nreq, v->reqs, MPI_STATUSES_IGNORE);
    timing(COMMUNICATION, timer_elapsed());

    dwrite(DEBUG_MPI, "%d: Waitall returned\n", iproc);
    profile(FACE_EXCHNG_REMOTE_WAIT);

    // unpack received data
    SB_struct *s = lsp;
    for (int face = 0; face < NUM_NEIGHBORS; face++)
    {
        int rank = s->neighbors[face][0];
        // A rank of less than 0 means that it isn't assigned
        if (rank >= 0 && rank != iproc)
        {
            unpack_plane(data, v->datasize, face, v->rbuf);
        }
    }

}

static void
ExchangeFacesForVar(
    int variable_key, void* d)
{
    /* Allocate list of MPI Requests.  One for each potential Send and one for
     * each potential Recv that we might do. */
    variable_registration *v = &var_regs[variable_key];
    size_t req_len = 2 * NUM_NEIGHBORS;

    double* dbuf;
    int* ibuf;

    if (v->reqs == NULL)
    {
        xrealloc(v->reqs, MPI_Request, req_len);
    }
    memset(v->reqs, '\0', sizeof(MPI_Request) * req_len);
    if (v->rbuf[0] == NULL)
    {
        int dimx = bp->gsdimx;
        int dimy = bp->gsdimy;
        int dimz = bp->gsdimz;

        int nxy = (dimx + 2) * (dimy + 2);
        int nxz = (dimx + 2) * (dimz + 2);
        int nyz = (dimy + 2) * (dimz + 2);
        int n2 = nxy;
        if (nxz > n2)
            n2 = nxz;
        if (nyz > n2)
            n2 = nyz;
        for (int face = 0; face < NUM_NEIGHBORS; face++)
        {
            v->rbuf[face] = malloc(v->datasize * n2);
            memset(v->rbuf[face], 1, v->datasize * n2);
            v->sbuf[face] = malloc(v->datasize * n2);
            memset(v->sbuf[face], 1, v->datasize * n2);
#ifdef GPU_PACK
switch( v->datasize)
{
    case 8:
        dbuf = (double*)v->rbuf[face];
#pragma omp target enter data map(to:dbuf[:n2])
        dbuf = (double*)v->sbuf[face];
#pragma omp target enter data map(to:dbuf[:n2])
        break;
   case 4:
       ibuf = (int*)v->rbuf[face];
#pragma omp target enter data map(to:ibuf[:n2])
       ibuf = (int*)v->sbuf[face];
#pragma omp target enter data map(to:ibuf[:n2])
       break;
   case 24:
        dbuf = (double*)v->rbuf[face];
#pragma omp target enter data map(to:dbuf[:3*n2])
        dbuf = (double*)v->sbuf[face];
#pragma omp target enter data map(to:dbuf[:3*n2])
   default:
       break;
}
#endif
        }
    }

    timing(COMPUTATION, timer_elapsed());

    // Post Receives & Sends
    {
        SB_struct *s = lsp;
        dwrite(DEBUG_TASK_CTRL, "Updating variable for MPI task %d\n",
               iproc);

        timing(COMPUTATION, timer_elapsed());

        v->nreq = SendRecvHalosNB(d, variable_key, v->datasize,
                                  s->neighbors, v->sbuf, v->rbuf,
                                  &v->reqs[0]);
    }
    profile(FACE_EXCHNG_REMOTE_SEND);
    timing(COMPUTATION, timer_elapsed());

    profile(FACE_EXCHNG_LOCAL);

}

/**
 * Main "do work" function
 */
double
doiteration(
    )
{
    static int grain_var = -1;
    static int d_var = -1;
    static int cl_var = -1;
    static int fs_var = 1;
    static int dc_var = -1;
    static int first_time = 1;

    static int cell_u_var = -1;
    static int cell_v_var = -1;
    static int cell_w_var = -1;


    if (first_time)
    {
        grain_var = registerCommInfo(sizeof(int));
        d_var = registerCommInfo(sizeof(double));
        cl_var = registerCommInfo(sizeof(double));
        fs_var = registerCommInfo(sizeof(double));
        dc_var = registerCommInfo(3 * sizeof(double));

        if (bp->fluidflow)
        {
            cell_u_var = registerCommInfo(sizeof(double));
            cell_v_var =
                registerCommInfo(sizeof(double));
            cell_w_var =
                registerCommInfo(sizeof(double));
        }

    }

    {
        // start communication for cl
        {
            ExchangeFacesForVar(cl_var, lsp->cl);
        }

        // start communication for fs
        {
            ExchangeFacesForVar(fs_var, lsp->fs);
        }

        // start communication for dc
        {
            ExchangeFacesForVar(dc_var, lsp->dc);
        }

        // Produces:  temperature
        {
            tempUpdate(false);
        }

        // Uses No Halo:  mold, gr, nuc_threshold, temperature
        // Uses w/ Halo:
        // Produces:  gr
        {
            cell_nucleation(lsp, NULL);
        }
        {
            activateNewGrains();
        }

        {
            ExchangeFacesForVar(grain_var, lsp->gr);
        }

        {
            FinishExchangeForVar(cl_var, lsp->cl);
        }
        {
            FinishExchangeForVar(fs_var, lsp->fs);
        }

        {
            FinishExchangeForVar(grain_var, lsp->gr);
        }

        // Uses No Halo: ce
        // Uses w/ Halo: cl, fs, mold
        // Produces: ce, oce
        {
            sb_diffuse_alloy_decentered(lsp, NULL);
            profile(CALC_DIFFUSE_ALLOY);
        }


        // Uses No Halo: gr, cl, ce, oce, fs, temperature
        // Uses w/ Halo:
        // Produces: fs, cl, gr
        {
            fs_change_diffuse(lsp, NULL);
            profile(CALC_FS_CHANGE);
        }

        // Uses No Halo: gr, fs, dc
        // Uses w/ Halo:
        // Produces: d, fs
        {
            grow_octahedra(lsp, NULL);
            profile(CALC_GROW_OCTAHEDRA);
        }

        // start communication for d
        {
            ExchangeFacesForVar(d_var, lsp->d);
        }

        // finish communications for dc
        {
            FinishExchangeForVar(dc_var, lsp->dc);
        }
        // finish communications for d
        {
            FinishExchangeForVar(d_var, lsp->d);
        }

#ifdef INDEX_SEP
        // Uses w/ Halo: gr
        {
            grow_cell_reduction(lsp, NULL);
            profile(CALC_CELL_INDEX);
        }
#endif

        // Uses No Halo: mold, fs, cl, ce, diff_id, ogr
        // Uses w/ Halo: gr, dc, d
        // Produces: gr, dc, fs, cl
        {
            capture_octahedra_diffuse(lsp);
            profile(CALC_CAPTURE_OCTAHEDRA);
            timing(COMPUTATION, timer_elapsed());
        }
    }

    timing(COMPUTATION, timer_elapsed());

    //                all tasks - synchronization point for each time step
    double gsolid_volume = solid_volume(lsp);
    double svol;
    MPI_Reduce(&gsolid_volume, &svol, 1, MPI_DOUBLE,
               MPI_SUM, 0, mpi_comm_new);
    profile(REDUCE_FS);
    timing(COMMUNICATION, timer_elapsed());

    first_time = 0;

    return svol;
}
