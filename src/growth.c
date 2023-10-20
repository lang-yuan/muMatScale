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
#include "profiler.h"
#include "fluidflow.h"
#include "omp.h"
#include "xmalloc.h"

double* fs;
double* cl;
double* d;
int* gr;
decentered_t *dc;

void
fs_dataexchange_to(
    SB_struct * lsp,
    void * __attribute__ ((__unused__)) __unused)
{
#ifdef GPU_OMP
#pragma omp target update to(fs[0:lsp->totaldim])  //nowait
#endif
}

void
fs_dataexchange_from(
    SB_struct * lsp,
    void * __attribute__ ((__unused__)) __unused)
{
#ifdef GPU_OMP
#pragma omp target update from(fs[0:lsp->totaldim])        //nowait
#endif
}

void
cl_dataexchange_to(
    SB_struct * lsp,
    void * __attribute__ ((__unused__)) __unused)
{
#ifdef GPU_OMP
#pragma omp target update to(cl[0:lsp->totaldim])  //nowait
#endif
}

void
cl_dataexchange_from(
    SB_struct * lsp,
    void * __attribute__ ((__unused__)) __unused)
{
#ifdef GPU_OMP
#pragma omp target update from(cl[0:lsp->totaldim])        //nowait
#endif
}

void
gr_dataexchange_to(
    SB_struct * lsp,
    void * __attribute__ ((__unused__)) __unused)
{
#ifdef GPU_OMP
#pragma omp target update to(gr[0:lsp->totaldim])  //nowait
#endif
}

void
gr_dataexchange_from(
    SB_struct * lsp,
    void * __attribute__ ((__unused__)) __unused)
{
#ifdef GPU_OMP
#pragma omp target update from(gr[0:lsp->totaldim])        //nowait
#endif
}

void
d_dataexchange_to(
    SB_struct * lsp,
    void * __attribute__ ((__unused__)) __unused)
{
#ifdef GPU_OMP
#pragma omp target update to(d[0:lsp->totaldim])   //nowait
#endif
}

void
d_dataexchange_from(
    SB_struct * lsp,
    void * __attribute__ ((__unused__)) __unused)
{
#ifdef GPU_OMP
#pragma omp target update from(d[0:lsp->totaldim]) //nowait
#endif
}

void
dc_dataexchange_to(
    SB_struct * lsp,
    void * __attribute__ ((__unused__)) __unused)
{
#ifdef GPU_OMP
#pragma omp target update to(dc[0:lsp->totaldim])  //nowait
#endif
}

void
dc_dataexchange_from(
    SB_struct * lsp,
    void * __attribute__ ((__unused__)) __unused)
{
#ifdef GPU_OMP
#pragma omp target update from(dc[0:lsp->totaldim])        //nowait
#endif
}

/**
 * Calculate the diffusion for each cell in the subblock
 * \param[in] vlsp Pointer to the subblock
 */
void
sb_diffuse_alloy_decentered(
    SB_struct * lsp,
    void * __attribute__ ((__unused__)) __unused)
{
    int dimx = bp->gsdimx;
    int dimy = bp->gsdimy;
    int dimz = bp->gsdimz;

    double dtx = bp->ts_delt / (bp->cellSize * bp->cellSize);
    double rs = (bp->Dsol * dtx);
    double rl = (bp->Dliq * dtx);

    /* check courant stability */
    if (rs > COURANT_LIMIT || rl > COURANT_LIMIT)
    {
        fprintf(stderr, "Courant Criterion: Solid: , %1.2e, Liquid: %1.2e\n",
                rs, rl);
        error
            ("Instability by Courant criterion, lower calculation time step: %g \n!",
             bp->ts_delt);
    }

#ifndef GPU_OMP
    memcpy(lsp->ogr, gr,
           (dimx + 2) * (dimy + 2) * (dimz + 2) * sizeof(int));
#endif

    int totaldim = (dimx + 2) * (dimy + 2) * (dimz + 2);

// copy ce into oce before updating ce
    double *oce = lsp->oce;
    double *ce = lsp->ce;
#ifdef GPU_OMP
#pragma omp target teams distribute parallel for
    for (int idx = 0; idx < totaldim; idx++)
    {
        oce[idx] = ce[idx];
    }
#else
    memcpy(oce, ce, totaldim * sizeof(double));
#endif

    int8_t *mold = lsp->mold;
#ifdef GPU_OMP
#pragma omp target teams distribute
#endif
    for (int k = 1; k <= dimz; k++)
    {
#ifdef GPU_OMP
#pragma omp parallel for collapse(2) schedule(static,1)
#endif
        for (int j = 1; j <= dimy; j++)
        {
            for (int i = 1; i <= dimx; i++)
            {
                int idx = k * (dimy + 2) * (dimx + 2) + j * (dimx + 2) + i;
                double nbsum = 0;
                double conc = cl[idx];
                double lspfs = fs[idx];

                int idxp = idx + 1;
                int idxm = idx - 1;
                int idyp = idx + dimx + 2;
                int idym = idx - (dimx + 2);
                int idzp = idx + (dimy + 2) * (dimx + 2);
                int idzm = idx - (dimy + 2) * (dimx + 2);
                // Calculate average of neighbors fraction solid
                double fs_av, r;


                fs_av = 0.5 * (lspfs + fs[idxp]);
                r = (rs * fs_av + rl * (1 - fs_av)) * (1 - mold[idxp]);
                double concxp = cl[idxp];
                nbsum += r * (concxp - conc);

                fs_av = 0.5 * (lspfs + fs[idxm]);
                r = (rs * fs_av + rl * (1 - fs_av)) * (1 - mold[idxm]);
                double concxm = cl[idxm];
                nbsum += r * (concxm - conc);

                fs_av = 0.5 * (lspfs + fs[idyp]);
                r = (rs * fs_av + rl * (1 - fs_av)) * (1 - mold[idyp]);
                double concyp = cl[idyp];
                nbsum += r * (concyp - conc);

                fs_av = 0.5 * (lspfs + fs[idym]);
                r = (rs * fs_av + rl * (1 - fs_av)) * (1 - mold[idym]);
                double concym = cl[idym];
                nbsum += r * (concym - conc);

                fs_av = 0.5 * (lspfs + fs[idzp]);
                r = (rs * fs_av + rl * (1 - fs_av)) * (1 - mold[idzp]);
                double conczp = cl[idzp];
                nbsum += r * (conczp - conc);

                fs_av = 0.5 * (lspfs + fs[idzm]);
                r = (rs * fs_av + rl * (1 - fs_av)) * (1 - mold[idzm]);
                double conczm = cl[idzm];
                nbsum += r * (conczm - conc);

/*
                if (bp->fluidflow)
                {

                    double cell_u = lsp->cell_u[idx];
                    double dcll = (conc - concxm) / bp->cellSize;
                    double dclr = (concxp - conc) / bp->cellSize;
                    double sgu = (cell_u >= 0 ? 1.0 : -1.0);
                    double dcldx =
                        (dcll + dclr + GAMMA * sgu * (dcll - dclr)) / 2.0;

                    double cell_v = lsp->cell_v[idx];
                    double dclb = (conc - concym) / bp->cellSize;
                    double dclf = (concyp - conc) / bp->cellSize;
                    sgu = (cell_v >= 0 ? 1.0 : -1.0);
                    double dcldy =
                        (dclb + dclf + GAMMA * sgu * (dclb - dclf)) / 2.0;

                    double cell_w = lsp->cell_w[idx];
                    double dcld = (conc - conczm) / bp->cellSize;
                    double dclu = (conczp - conc) / bp->cellSize;
                    sgu = (cell_w >= 0 ? 1.0 : -1.0);
                    double dcldz =
                        (dcld + dclu + GAMMA * sgu * (dcld - dclu)) / 2.0;

                    double vtoce =
                        (cell_u * dcldx + cell_v * dcldy +
                         cell_w * dcldz) * (1 - fs[idx]);
                    vtoce = vtoce * bp->ts_delt;

                    nbsum = nbsum - vtoce;
                }
*/
/*
     for (int n = 0; n < nn_num; n++)
	 {
	 	int ndx = (k + nn[n][2]) * (dimy+2)*(dimx+2) +
	 			(j + nn[n][1]) * (dimx+2) +
	 			(i + nn[n][0]);
         double fs_av = 0.5 * (lspfs + fs[ndx]);
		 double r = rs * fs_av + rl * (1 - fs_av);
		 double nbconc = cl[ndx];
		 nbsum += r * (nbconc - conc);
	  }
*/
                //Update the final ce value

                ce[idx] += nbsum;

            }
        }
    }
// profile(TIMER_4);
}

void
pre_cell_reduction(
    SB_struct * lsp,
    void * __attribute__ ((__unused__)) __unused)
{
    int dimx = bp->gsdimx;
    int dimy = bp->gsdimy;
    int dimz = bp->gsdimz;
    int gindex = 0;
    int nindex = 0;
    int fsindex = 0;

    double k_inv = 1.0 / bp->part_coef;

    int8_t *mold = lsp->mold;
    double *ce = lsp->ce;
    int *diff_id = lsp->diff_id;
#if defined(GPU_OMP)
#pragma omp target map(tofrom:gindex, nindex, fsindex)
#pragma omp teams distribute
#endif
    for (int k = 1; k <= dimz; k++)
    {
#if defined(GPU_OMP)
#pragma omp parallel for collapse(2) schedule(static,1)
#endif
        for (int j = 1; j <= dimy; j++)
        {
            for (int i = 1; i <= dimx; i++)
            {
                int idx = k * (dimy + 2) * (dimx + 2) + j * (dimx + 2) + i;
                int idxp = idx + 1;
                int idxm = idx - 1;
                int idyp = idx + dimx + 2;
                int idym = idx - (dimx + 2);
                int idzp = idx + (dimy + 2) * (dimx + 2);
                int idzm = idx - (dimy + 2) * (dimx + 2);

                if (mold[idx])
                    continue;

#ifdef SEP_LISTS
                if (gr[idx] <= 0)  // liquid cell
                {
                    cl[idx] = ce[idx];
                    int nbindex =
                        gr[idxp] + gr[idxm] + gr[idyp] +
                        gr[idym] + gr[idzp] + gr[idzm];
                    if (nbindex > 0)
                    {
                        int nn = 0;
#pragma omp atomic capture
                        {
                            gindex = gindex + 1;
                            nn = gindex;
                        }
                        diff_id[nn] = idx;

                    }
                    continue;
                }

                if (fs[idx] == 1)
                {
                    cl[idx] = ce[idx] * k_inv;
                    continue;
                }

                if (fs[idx] < 1.0)
                {

                    int ff = 0;
#pragma omp atomic capture
                    {
                        fsindex = fsindex + 1;
                        ff = fsindex;
                    }
                    fs_id[ff] = idx;

                }
#else // NOT SEP_LISTS
                if (fs[idx] >= 1.0)
                {
                    cl[idx] = ce[idx] * k_inv;
                    continue;
                }

                if (gr[idx] <= 0)  // liquid cell
                    cl[idx] = ce[idx];

                int nbindex =
                    gr[idx] + gr[idxp] + gr[idxm] +
                    gr[idyp] + gr[idym] + gr[idzp] +
                    gr[idzm];
                if (nbindex > 0)
                {
                    int nn = 0;
#pragma omp atomic capture
                    {
                        gindex = gindex + 1;
                        nn = gindex;
                    }
                    diff_id[nn] = idx;
                }

#endif
            }
        }
    }

#ifdef SEP_LISTS
    lsp->gindex = gindex;
    lsp->fsindex = fsindex;
    lsp->growindex = fsindex;
#else
    lsp->gindex = gindex;
    lsp->fsindex = gindex;
    lsp->growindex = gindex;

#endif

//#ifdef NUC_PRELIST
//    lsp->nindex = nindex;
//#endif

}

/**
 * Update the fraction solid across the subblock
 * \param[in] vlsp Pointer to the subblock
 */
void
fs_change_diffuse(
    SB_struct * lsp,
    void * __attribute__ ((__unused__)) __unused)
{
    const double liq_slope = bp->m_solute0;

    const double k_inv = 1.0 / bp->part_coef;
    const double km = 1.0 - bp->part_coef;
    const double km_inv = 1.0 / km;
    const double m_inv = 1.0 / liq_slope;

    const double D_x = bp->Dliq / bp->cellSize / bp->cellSize;
    const double D_xt = D_x * bp->ts_delt;

    const int dimx = bp->gsdimx;
    const int dimy = bp->gsdimy;
    const int dimz = bp->gsdimz;

    int totaldim = (dimx + 2) * (dimy + 2) * (dimz + 2);
    const double cinit = bp->Cinit;
    const double Tliq = bp->liquidusTemp;

    double *oce = lsp->oce;
    double *ce = lsp->ce;
    double *temperature = lsp->temperature;
#ifndef FSSEP

#if defined(GPU_OMP)
#pragma omp target teams distribute
#endif
    for (int k = 1; k <= dimz; k++)
    {
#if defined(GPU_OMP)
#pragma omp parallel for collapse(2) schedule(static,1)
#endif
        for (int j = 1; j <= dimy; j++)
        {
            for (int i = 1; i <= dimx; i++)
            {
                int idx = k * (dimy + 2) * (dimx + 2) + j * (dimx + 2) + i;
                double Tcell = temperature[idx];
                double c_int = cinit + m_inv * (Tcell - Tliq);
                double Tunder = Tliq - Tcell;

                Tunder += liq_slope * (cl[idx] - cinit);


                if (gr[idx] <= 0 || Tcell < 1.0 || Tcell > Tliq)
                {
                    cl[idx] = ce[idx];
                    lsp->ogr[idx] = 0;
                    gr[idx] = 0;
                    fs[idx] = 0;
                    continue;
                }

                if (fs[idx] == 1.0 && (Tunder >= 0.0))
                {
                    cl[idx] = ce[idx] * k_inv;
                    continue;
                }

                if (bp->tip_curv)
                {
                    c_int =
                        cinit + m_inv * (Tcell - Tliq +
                                         bp->gibbs_coef * lsp->curv[idx]);
                    Tunder = Tliq - Tcell - bp->gibbs_coef * lsp->curv[idx];
                }

                // Note: The following line depends on the memcpy
                // from the previous function
                // Shouldnt need to be exchanged since not looking
                // at neighbors
                double dce = ce[idx] - oce[idx];
                double dcl =
                    (c_int - cl[idx]) * D_xt * 2 + dce / (1 -
                                                               km *
                                                               fs[idx]);
                double dfs =
                    (dcl * (1. - km * fs[idx]) -
                     dce) * km_inv / cl[idx];
                cl[idx] += dcl;
                fs[idx] += dfs;
                fs[idx] = MIN(1.0, fs[idx]);

                if (fs[idx] < 0.)
                {
                    fs[idx] = 0.;
                    gr[idx] = 0;
                    lsp->ogr[idx] = 0;
                    dfs = 0;
                    cl[idx] = ce[idx];
                }

            }
        }
    }

#else

    int fsindex = fsindex;
    double *oce = lsp->oce;
#if defined(GPU_OMP)
#pragma omp target teams distribute parallel for schedule(static,1)
#endif
    for (int i = 1; i <= fsindex; i++)
    {
#ifdef SEP_LISTS
        int idx = fs_id[i];
#else
        int idx = diff_id[i];

        if (gr[idx] <= 0)
        {
            cl[idx] = ce[idx];
            continue;
        }

        if (fs[idx] == 1.0)
        {
            cl[idx] = ce[idx] * k_inv;
            continue;
        }

#endif
        double Tcell = temperature[idx];
        double fs = fs[idx];

        double c_int = cinit + m_inv * (Tcell - Tliq);

        double dce = ce[idx] - oce[idx];
        double dcl = (c_int - cl[idx]) * D_xt * 2 + dce / (1 - km * fs);
        double dfs = (dcl * (1. - km * fs) - dce) * km_inv / cl[idx];
        cl[idx] += dcl;
        fs[idx] += dfs;
        fs[idx] = MIN(1.0, fs[idx]);

        if (fs[idx] < 0.)
        {
            fs[idx] = 0.;
            gr[idx] = 0;
            cl[idx] = ce[idx];
        }

    }

#endif //FSSEP

//    profile(CALC_FS_CHANGE);
}

/**
 * Based on the growth rate, increase the size of the octahedral grains
 * by increasing the diameter
 * \param[in] vlsp A pointer to the subblock
 */


void
grow_octahedra(
    SB_struct * lsp,
    void * __attribute__ ((__unused__)) __unused)
{
    const int dimx = bp->gsdimx;
    const int dimy = bp->gsdimy;
    const int dimz = bp->gsdimz;

    const double fsgrow = bp->rev_fsgrow;

    int totaldim = (dimx + 2) * (dimy + 2) * (dimz + 2);

#ifndef GROWSEP

#if defined(GPU_OMP)
#pragma omp target teams distribute
#endif
    for (int k = 1; k <= dimz; k++)
    {
#if defined(GPU_OMP)
#pragma omp parallel for collapse(2) schedule(static,1)
#endif

        for (int j = 1; j <= dimy; j++)
        {
            for (int i = 1; i <= dimx; i++)
            {
                int idx = k * (dimy + 2) * (dimx + 2) + j * (dimx + 2) + i;
                int grid = gr[idx];

                // if cell part of a grain, grow it
                if (grid > 0 && grid < bp->maxTotalGrains)
                {
                    grain_t grain = grain_cache[grid];

                    double g00 = grain.rotmat[0][0];
                    double g10 = grain.rotmat[1][0];
                    double g20 = grain.rotmat[2][0];

                    double g01 = grain.rotmat[0][1];
                    double g11 = grain.rotmat[1][1];
                    double g21 = grain.rotmat[2][1];

                    double g02 = grain.rotmat[0][2];
                    double g12 = grain.rotmat[1][2];
                    double g22 = grain.rotmat[2][2];

                    double rx0 = 0.0 - dc[idx].x;
                    double ry0 = 0.0 - dc[idx].y;
                    double rz0 = 0.0 - dc[idx].z;

                    double gx0 = g00 * rx0 + g10 * ry0 + g20 * rz0;
                    double gy0 = g01 * rx0 + g11 * ry0 + g21 * rz0;
                    double gz0 = g02 * rx0 + g12 * ry0 + g22 * rz0;

                    double ds1 = ABS(gx0) + ABS(gy0) + ABS(gz0);
                    ds1 = MIN(ds1, 1.5);

                    d[idx] = ds1 + fs[idx] * fsgrow;
                }
            }
        }
    }

#else // GROWSEP

    int growindex = growindex;

#if defined(GPU_OMP)
#pragma omp target teams distribute parallel for schedule(static,1)
#endif
    for (int i = 1; i <= growindex; i++)
    {
#ifdef SEP_LISTS
        int idx = fs_id[i];
#else
        int idx = diff_id[i];
        if (gr[idx] <= 0)
            continue;
#endif
        //int idx = grow_id[i];
        int gid = gr[idx];
        grain_t grain = grain_cache[gid];

        double g00 = grain.rotmat[0][0];
        double g10 = grain.rotmat[1][0];
        double g20 = grain.rotmat[2][0];

        double g01 = grain.rotmat[0][1];
        double g11 = grain.rotmat[1][1];
        double g21 = grain.rotmat[2][1];

        double g02 = grain.rotmat[0][2];
        double g12 = grain.rotmat[1][2];
        double g22 = grain.rotmat[2][2];

        double rx0 = 0.0 - dc[idx].x;
        double ry0 = 0.0 - dc[idx].y;
        double rz0 = 0.0 - dc[idx].z;

        double gx0 = g00 * rx0 + g10 * ry0 + g20 * rz0;
        double gy0 = g01 * rx0 + g11 * ry0 + g21 * rz0;
        double gz0 = g02 * rx0 + g12 * ry0 + g22 * rz0;

        double ds1 = ABS(gx0) + ABS(gy0) + ABS(gz0);
        ds1 = MIN(ds1, 1.5);

        d[idx] = ds1 + fs[idx] * fsgrow;
    }
#endif //GROWSEP

    int *ogr = lsp->ogr;
#if defined(GPU_OMP)
#pragma omp target teams distribute parallel for schedule(static, 1)
#endif
    for (int i = 0; i < totaldim; i++)
    {
        ogr[i] = gr[i];
    }

}

void
grow_cell_reduction(
    SB_struct * lsp,
    void * __attribute__ ((__unused__)) __unused)
{
    int dimx = bp->gsdimx;
    int dimy = bp->gsdimy;
    int dimz = bp->gsdimz;
    int gindex = 0;
    int nindex = 0;

    int8_t *mold = lsp->mold;
    int *diff_id = lsp->diff_id;
    int *nuc_id = lsp->nuc_id;
#if defined(GPU_OMP)
#pragma omp target map(tofrom:gindex, nindex)
#pragma omp teams distribute
#endif
    for (int k = 1; k <= dimz; k++)
    {
#if defined(GPU_OMP)
#pragma omp parallel for collapse(2) schedule(static,1)
#endif
        for (int j = 1; j <= dimy; j++)
        {
            for (int i = 1; i <= dimx; i++)
            {
                int idx = k * (dimy + 2) * (dimx + 2) + j * (dimx + 2) + i;
                int idxp = idx + 1;
                int idxm = idx - 1;
                int idyp = idx + dimx + 2;
                int idym = idx - (dimx + 2);
                int idzp = idx + (dimy + 2) * (dimx + 2);
                int idzm = idx - (dimy + 2) * (dimx + 2);

                if (mold[idx])
                    continue;

                if (gr[idx] <= 0)  // liquid cell
                {
                    int nbindex =
                        gr[idxp] + gr[idxm] + gr[idyp] +
                        gr[idym] + gr[idzp] + gr[idzm];
                    if (nbindex > 0)
                    {
                        int nn = 0;
#pragma omp atomic capture
                        {
                            gindex = gindex + 1;
                            nn = gindex;
                        }
                        diff_id[nn] = idx;

                    }
#ifdef NUC_PRELIST
                    else
                    {
                        int mm = 0;
#pragma omp atomic capture
                        {
                            nindex = nindex + 1;
                            mm = nindex;
                        }
                        nuc_id[mm] = idx;

                    }
#endif // NUC_PRELIST
                }

            }
        }
    }

    lsp->gindex = gindex;
#ifdef NUC_PRELIST
    lsp->nindex = nindex;
#endif
}


/**
 * Determine if a growing octahedral grain captures a cell
 * \param[in] vlsp A pointer to the subblock
 * \param[out] vnum_solid The number of solid cells in the subblock
 */
void
capture_octahedra_diffuse(
    SB_struct * lsp,
    void *vsolid_volume)
{
    double *solid_volume = (double *) vsolid_volume;
    double sum_fs = 0.0;
    uint64_t solid_count = 0;


    int dimx = bp->gsdimx;
    int dimy = bp->gsdimy;
    int dimz = bp->gsdimz;

    static decentered_t *dc_old = NULL;
    if (dc_old == NULL)
    {
        allocate_decentered(&dc_old);
    }

#ifndef GPU_OMP
    int sizedc =
        sizeof(decentered_t) * ((bp->gsdimx + 2) * (bp->gsdimy + 2) *
                                (bp->gsdimz + 2));
    memcpy(dc_old, dc, sizedc);
#endif

    dwrite(DEBUG_TASK_CTRL,
           "%d: Calculating for subblock %d (%d x %d x %d)\n", iproc,
           lsp->subblockid, dimx, dimy, dimz);

    // Neighbor mapping x, y, z
    // NOTE: This creates a subtle bias based on neighbor ordering. The ordering is this way to match umatic as of 4/13/2011
    const short nn[6][3] = { {1, 0, 0},
    {0, 1, 0},
    {-1, 0, 0},
    {0, -1, 0},
    {0, 0, 1},
    {0, 0, -1}
    };
    const int nn_num = 6;

//General random number locally
#ifdef RAND_LOCAL
    double *randnum;
    allocate_float(&randnum);
    int totaldim = (dimx + 2) * (dimy + 2) * (dimz + 2);
    for (int i = 0; i <= totaldim; i++)
    {
        randnum[i] = (float) rand() / (float) RAND_MAX;
    }
#endif

    double *ce = lsp->ce;
    int *ogr = lsp->ogr;

#ifdef INDEX_SEP

    int gindex = lsp->gindex;

    static int dc_tmp_size = 0;
    static decentered_t *dc_tmp;
    if (gindex > dc_tmp_size)
    {
        if (dc_tmp_size > 0)
        {
#if defined(GPU_OMP)
#pragma omp target exit data map(release: dc_tmp[0:dc_tmp_size])
#endif
            xfree(dc_tmp);
        }
        dc_tmp_size = gindex + 100;
        xmalloc(dc_tmp, decentered_t, dc_tmp_size);
#if defined(GPU_OMP)
#pragma omp target enter data map(alloc:dc_tmp[0:dc_tmp_size])
#endif
    }

    int *diff_id = lsp->diff_id;
#if defined(GPU_OMP)
#pragma omp target teams distribute parallel for schedule(static,1)
#endif
    for (int i = 1; i <= gindex; i++)
    {
        int idx = diff_id[i];

#ifndef SEP_LISTS
        if (ogr[idx] > 0)
            continue;
#endif
        const int ndxs[6] =
            { idx + 1, idx + dimx + 2, idx - 1, idx - (dimx + 2),
            idx + (dimy + 2) * (dimx + 2), idx - (dimy + 2) * (dimx + 2)
        };

        double g00;
        double g10;
        double g20;
        double g01;
        double g11;
        double g21;
        double g02;
        double g12;
        double g22;
        double dx;
        double dy;
        double dz;

        int ndx;
        short found = -1;

        // loop over 6 neighbors
        for (short n = 0; n < 6; n++)
        {
            ndx = ndxs[n];
            int gid = ogr[ndx];
            if (gid > 0 && gid < bp->maxTotalGrains)
            {
                grain_t grain = grain_cache[gid];

                // Determine if my neighbor captures me
                g00 = grain.rotmat[0][0];
                g10 = grain.rotmat[1][0];
                g20 = grain.rotmat[2][0];

                g01 = grain.rotmat[0][1];
                g11 = grain.rotmat[1][1];
                g21 = grain.rotmat[2][1];

                g02 = grain.rotmat[0][2];
                g12 = grain.rotmat[1][2];
                g22 = grain.rotmat[2][2];

                double rx1 = -1.0 * nn[n][0] - dc[ndx].x;
                double ry1 = -1.0 * nn[n][1] - dc[ndx].y;
                double rz1 = -1.0 * nn[n][2] - dc[ndx].z;

                dx = g00 * rx1 + g10 * ry1 + g20 * rz1;
                dy = g01 * rx1 + g11 * ry1 + g21 * rz1;
                dz = g02 * rx1 + g12 * ry1 + g22 * rz1;

                double dfs = d[ndx] - (ABS(dx) + ABS(dy) + ABS(dz));

#ifdef RAND_LOCAL
                dfs += 0.05 * randnum[ndx];
#endif

                if (dfs > 0)
                {
                    found = n;
                    break;
                }
            }
        }

        if (found > -1)
        {
            double ncx = 0.;
            double ncz = 0.;
            double ncy = 0.;

            double ncd = MIN(d[ndx], 1.0);
            double dd = d[ndx] - ncd;

            /* determine which corner the new capture cell is close to */
            if (dx + dy >= 0 && dx - dy > 0 && dx + dz > 0 && dx - dz >= 0)
            {
                ncx = dd;
            }
            else if (dx + dy <= 0 && dx - dy < 0 && dx + dz < 0
                     && dx - dz <= 0)
            {
                ncx = -1. * dd;
            }
            else if (dy + dx > 0 && dy - dx >= 0 && dy + dz >= 0
                     && dy - dz > 0)
            {
                ncy = dd;
            }
            else if (dy + dx < 0 && dy - dx <= 0 && dy + dz <= 0
                     && dy - dz < 0)
            {
                ncy = -1. * dd;
            }
            else if (dz + dx >= 0 && dz - dx > 0 && dz + dy > 0
                     && dz - dy >= 0)
            {
                ncz = dd;
            }
            else
            {
                ncz = -1. * dd;
            }

            // set dc_tmp based on value of dc at neighboring cell
            dc_tmp[i].x =
                dc[ndx].x + nn[found][0] + g00 * ncx + g01 * ncy +
                g02 * ncz;
            dc_tmp[i].y =
                dc[ndx].y + nn[found][1] + g10 * ncx + g11 * ncy +
                g12 * ncz;
            dc_tmp[i].z =
                dc[ndx].z + nn[found][2] + g20 * ncx + g21 * ncy +
                g22 * ncz;

            gr[idx] = ogr[ndx];
            fs[idx] = 0.;  // Just got captured, set my fraction solid to 0
            cl[idx] = ce[idx];

        }
    }

    // update dc
#if defined(GPU_OMP)
#pragma omp target teams distribute parallel for schedule(static,1)
#endif
    for (int i = 1; i <= gindex; i++)
    {
        int idx = diff_id[i];

        dc[idx].x = dc_tmp[i].x;
        dc[idx].y = dc_tmp[i].y;
        dc[idx].z = dc_tmp[i].z;
    }

#if defined(GPU_OMP)
#pragma omp target teams distribute parallel for collapse(3) reduction(+:sum_fs) schedule(static,1)
#endif
    for (int k = 1; k <= dimz; k++)
    {
        for (int j = 1; j <= dimy; j++)
        {
            for (int i = 1; i <= dimx; i++)
            {
                int idx = k * (dimy + 2) * (dimx + 2) + j * (dimx + 2) + i;
                sum_fs += MIN(1.0, fs[idx]);
            }
        }
    }

#else // orignal code without INDEX_SEP


    for (int k = 1; k <= dimz; k++)
    {
        for (int j = 1; j <= dimy; j++)
        {
            for (int i = 1; i <= dimx; i++)
            {
                int idx = k * (dimy + 2) * (dimx + 2) + j * (dimx + 2) + i;
                if (mold[idx])
                {
                    continue;
                }

                if (ogr[idx] <= 0)      // liquid cell
                {
                    int nbindex = 0;
                    for (int n = 0; n < nn_num; n++)
                    {
                        int ndx = (k + nn[n][2]) * (dimy + 2) * (dimx + 2) +
                            (j + nn[n][1]) * (dimx + 2) + (i + nn[n][0]);
                        nbindex += ogr[ndx];
                    }

                    if (nbindex ? 1 : 0)
                    {

                        int capt = 0;   // Cell not captured -- should be a boolean

                        //Check my neighbors to see if any of them are solid
                        for (int n = 0; n < nn_num && capt == 0; n++)
                        {
                            int ndx =
                                (k + nn[n][2]) * (dimy + 2) * (dimx + 2) +
                                (j + nn[n][1]) * (dimx + 2) + (i + nn[n][0]);
                            int gid = ogr[ndx];

                            if (gid > 0)
                            {
                                grain_t grain = grain_cache[gid];

                                // Determine if my neighbor captures me
                                double g00 = grain.rotmat[0][0];
                                double g10 = grain.rotmat[1][0];
                                double g20 = grain.rotmat[2][0];

                                double g01 = grain.rotmat[0][1];
                                double g11 = grain.rotmat[1][1];
                                double g21 = grain.rotmat[2][1];

                                double g02 = grain.rotmat[0][2];
                                double g12 = grain.rotmat[1][2];
                                double g22 = grain.rotmat[2][2];

                                double rx1 = -1.0 * nn[n][0] - dc_old[ndx].x;
                                double ry1 = -1.0 * nn[n][1] - dc_old[ndx].y;
                                double rz1 = -1.0 * nn[n][2] - dc_old[ndx].z;

                                double dx = g00 * rx1 + g10 * ry1 + g20 * rz1;
                                double dy = g01 * rx1 + g11 * ry1 + g21 * rz1;
                                double dz = g02 * rx1 + g12 * ry1 + g22 * rz1;

                                double ds = ABS(dx) + ABS(dy) + ABS(dz);

                                double dfs = d[ndx] - ds;
                                //  dfs += 0.05 * getRandScale();

                                if (dfs > 0)
                                {
                                    double ncx, ncz, ncy;
                                    capt = 1;

                                    double ncd = MIN(d[ndx], 1.0);

                                    double dd = d[ndx] - ncd;

                                    /* determine which corner the new capture cell is close to */
                                    if (dx + dy >= 0 && dx - dy > 0
                                        && dx + dz > 0 && dx - dz >= 0)
                                    {
                                        ncx = dd;
                                        ncy = 0.;
                                        ncz = 0.;

                                    }
                                    else if (dx + dy <= 0 && dx - dy < 0
                                             && dx + dz < 0 && dx - dz <= 0)
                                    {
                                        ncx = -1. * dd;
                                        ncy = 0.;
                                        ncz = 0.;

                                    }
                                    else if (dy + dx > 0 && dy - dx >= 0
                                             && dy + dz >= 0 && dy - dz > 0)
                                    {
                                        ncx = 0.;
                                        ncy = dd;
                                        ncz = 0.;

                                    }
                                    else if (dy + dx < 0 && dy - dx <= 0
                                             && dy + dz <= 0 && dy - dz < 0)
                                    {
                                        ncx = 0.;
                                        ncy = -1. * dd;
                                        ncz = 0.;

                                    }
                                    else if (dz + dx >= 0 && dz - dx > 0
                                             && dz + dy > 0 && dz - dy >= 0)
                                    {
                                        ncx = 0.;
                                        ncy = 0.;
                                        ncz = dd;

                                    }
                                    else if (dz + dx <= 0 && dz - dx < 0
                                             && dz + dy < 0 && dz - dy <= 0)
                                    {
                                        ncx = 0.;
                                        ncy = 0.;
                                        ncz = -1. * dd;
                                    }
                                    else
                                    {
                                        ncx = 0.;
                                        ncy = 0.;
                                        ncz = 0.;
                                    }

                                    gr[idx] = ogr[ndx];

                                    dc[idx].x =
                                        dc_old[ndx].x + nn[n][0] + g00 * ncx +
                                        g01 * ncy + g02 * ncz;
                                    dc[idx].y =
                                        dc_old[ndx].y + nn[n][1] + g10 * ncx +
                                        g11 * ncy + g12 * ncz;
                                    dc[idx].z =
                                        dc_old[ndx].z + nn[n][2] + g20 * ncx +
                                        g21 * ncy + g22 * ncz;


                                    dfs = 0;
                                    fs[idx] = dfs; // Just got captured, set my fraction solid to 0
                                    cl[idx] = ce[idx];

                                }
                            }
                        }
                    }
                }

                sum_fs += MIN(1.0, fs[idx]);
            }
        }
    }

#endif //INDEX_SEP

    (*solid_volume) += sum_fs * pow(bp->cellSize, 3.0);
    dwrite(DEBUG_TASK_CTRL, "(sb: %d) Returning %lu solids\n",
           lsp->subblockid, solid_count);
}
