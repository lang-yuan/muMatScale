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

#include "file_io.h"
#include "globals.h"
#include "functions.h"
#include "setup.h"
#include "face_util.h"
#include "xmalloc.h"
#include "debug.h"
#include "temperature.h"
#include "grain.h"
#include "profiler.h"

#ifdef USE_MPE
#include <mpe.h>
extern int gLogStateBeginID;
extern int gLogStateEndID;
#endif

extern SB_struct *lsp;


/**
 * Sets up the initial conditions and initializes a new subblock
 * \param lsp Pointer to subblock to be set up
 */
static void
SetupInitialConditions(
    SB_struct * lsp)
{
    int dimx = bp->gsdimx;
    int dimy = bp->gsdimy;
    int dimz = bp->gsdimz;
    int full_volume = (dimx + 2) * (dimy + 2) * (dimz + 2);
    dwrite(DEBUG_TASK_CTRL, "%d: Initializing subblock %d (%d x %d x %d)\n",
           iproc, lsp->subblockid, dimx, dimy, dimz);

    // Create output directory for this subblock
    prepareIOSubblock(lsp);

    // Set initial temp
    tempPrepareSB(lsp);

    // Initialize the grains to 0 and fraction solid to 0
    for (int i = 0; i < full_volume; i++)
    {
        lsp->nuc_threshold[i] = INFINITY;

        if (bp->calc_type == DIFFUSION)
        {
            lsp->ce[i] = lsp->cl[i] = bp->Cinit;
            lsp->nuc_id[i] = i;
            //ADD_CURV_LY
            if (bp->tip_curv == 1)
                lsp->curv[i] = 0.0;

            if (bp->fluidflow)
            {
                lsp->cell_u[i] = bp->initvelo[0];
                lsp->cell_v[i] = bp->initvelo[1];
                lsp->cell_w[i] = bp->initvelo[2];

            }
        }
    }

    lsp->totaldim = (dimx + 2) * (dimy + 2) * (dimz + 2);
    lsp->gindex = 0;
    lsp->nindex = lsp->totaldim;

    double low_x, low_y, low_z;
    double high_x, high_y, high_z;
    getSBRealBounds(lsp->coords.x, lsp->coords.y, lsp->coords.z, &low_x,
                    &low_y, &low_z, &high_x, &high_y, &high_z);
    dwrite(DEBUG_TASK_CTRL,
           "SB %d checking for nuc points in box [%g %g] [%g %g] [%g %g]\n",
           lsp->subblockid, low_x, high_x, low_y, high_y, low_z, high_z);

    init_sb_nucleation(lsp);
}

// Calculate the number of cells that are mold
static void
countMold(
    SB_struct * lsp)
{
    int dimx = bp->gsdimx;
    int dimy = bp->gsdimy;
    int dimz = bp->gsdimz;

    for (int z = 1; z <= dimz; z++)
        for (int y = 1; y <= dimy; y++)
            for (int x = 1; x <= dimx; x++)
                if (lsp->mold[SBIDX(x, y, z)])
                    lsp->mold_size++;
}

void
allocateFields(
    )
{
    allocate_byte(&lsp->mold);
    allocate_float(&(lsp->fs));
    allocate_float(&(lsp->temperature));

    xmalloc(lsp->nuc_threshold, float,
              (bp->gsdimx + 2) * (bp->gsdimy + 2) * (bp->gsdimz + 2));

    // Allocate arrays for diffusion - should only do if necessary
    if (bp->calc_type == DIFFUSION)
    {
        allocate_float(&(lsp->ce));
        allocate_float(&(lsp->oce));
        allocate_float(&(lsp->cl));
        //ADD_CURV_LY
        if (bp->tip_curv == 1)
            allocate_float(&(lsp->curv));

        allocate_int(&(lsp->diff_id));
        allocate_int(&(lsp->fs_id));
        allocate_int(&(lsp->nuc_id));
        allocate_int(&(lsp->nuc_id2));

        //Fluidflow
        if (bp->fluidflow)
        {
            allocate_float(&(lsp->cell_u));
            allocate_float(&(lsp->cell_v));
            allocate_float(&(lsp->cell_w));
        }
    }

    // Allocate arrays for decentered octahedron information
    allocate_decentered(&(lsp->dc));
    allocate_float(&(lsp->d));

    allocate_int(&(lsp->gr));
    allocate_int(&(lsp->ogr));

    allocate_int(&(lsp->lsindex));

    // allocate arrays for IO
    allocate_sp(&(lsp->fs_io));
    allocate_sp(&(lsp->temperature_io));
    allocate_int(&(lsp->gr_io));
    if (bp->calc_type == DIFFUSION)
    {
        allocate_sp(&(lsp->ce_io));
    }

    int totaldim = (bp->gsdimx + 2) * (bp->gsdimy + 2) * (bp->gsdimz + 2);

    memset(lsp->cl, 0, totaldim*sizeof(double));
    memset(lsp->d, 0, totaldim*sizeof(double));
    memset(lsp->ogr, 0, totaldim*sizeof(int));

#ifdef GPU_OMP

#pragma omp target enter data map(to:bp[0:1])

    double* d = lsp->d;
#pragma omp target enter data map(to:d[0:totaldim]) nowait
    double* fs = lsp->fs;
#pragma omp target enter data map(to:fs[0:totaldim]) nowait
    int* gr = lsp->gr;
#pragma omp target enter data map(to:gr[0:totaldim]) nowait
    int *ogr = lsp->ogr;
#pragma omp target enter data map(to:ogr[0:totaldim]) nowait
    double *temperature = lsp->temperature;
#pragma omp target enter data map(to:temperature[0:totaldim])
    int* lsindex = lsp->lsindex;
#pragma omp target enter data map(to:lsindex[0:totaldim])

#endif

    SetupInitialConditions(lsp);
    dprintf(stderr, "complete initial Rank %d \n", iproc);

#ifdef GPU_OMP
    double* cl = lsp->cl;
#pragma omp target enter data map(to:cl[0:totaldim]) nowait
    decentered_t* dc = lsp->dc;
#pragma omp target enter data map(to:dc[0:totaldim]) nowait

#ifdef GPU_OMP_NUC
    float *nuc_threshold = lsp->nuc_threshold;
#pragma omp target enter data map(to:nuc_threshold[0:totaldim])
#endif

    int *fs_id = lsp->fs_id;
#pragma omp target enter data map(to:fs_id[0:totaldim])
    int *nuc_id = lsp->nuc_id;
#pragma omp target enter data map(to:nuc_id[0:totaldim])
    int *nuc_id2 = lsp->nuc_id2;
#pragma omp target enter data map(to:nuc_id2[0:totaldim])
    int *diff_id = lsp->diff_id;
#pragma omp target enter data map(to:diff_id[0:totaldim])
    int8_t *mold = lsp->mold;
#pragma omp target enter data map(to:mold[0:totaldim]) nowait
    double *oce = lsp->oce;
#pragma omp target enter data map(to:oce[0:totaldim]) nowait
    double *ce = lsp->ce;
#pragma omp target enter data map(to:ce[0:totaldim]) nowait

// IO arrays
    float* fs_io = lsp->fs_io;
#pragma omp target enter data map(to:fs_io[0:totaldim]) nowait
    float *ce_io = lsp->ce_io;
#pragma omp target enter data map(to:ce_io[0:totaldim]) nowait
    float *temperature_io = lsp->temperature_io;
#pragma omp target enter data map(to:temperature_io[0:totaldim])
    int* gr_io = lsp->gr_io;
#pragma omp target enter data map(to:gr_io[0:totaldim]) nowait

#endif

    countMold(lsp);

    dwrite(DEBUG_TASK_CTRL,
           "%d: Allocated memory for 1 subblocks of size %d x %d x %d\n",
           iproc, bp->gsdimx, bp->gsdimy, bp->gsdimz);
}
