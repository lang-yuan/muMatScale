/***************************************************************/
/* Copyright (c) 2023, Lang Yuan, Univeristy of South Carolina */
/* All rights reserved.                                        */
/* This file is part of muMatScale.                            */
/* See the top-level LICENSE file for details.                 */
/***************************************************************/

#include "ll.h"
#include "globals.h"
#include "functions.h"
#include "temperature.h"

#include <math.h>
#include <assert.h>

extern SB_struct *lsp;

/****** INTERNAL TEMPERATURE FUNCTIONS ******/
#ifdef GPU_OMP
#pragma omp declare target
#endif
static double
direct_temp_calc(
    uint32_t sbx,
    uint32_t sby,
    uint32_t sbz,
    size_t x,
    size_t y,
    size_t z)
{
    internal_temp_ctrl_t *t = &bp->temp_ctrl;
    double sim_time = bp->timestep * bp->ts_delt;
    double T;

    if (bp->thermal_ana == 0 && t->gradient == 0)
    {
        T = t->initialTemperature - sim_time * t->cool_rate;
    }
    else
    {

        double rx, ry, rz;
        scoord2realcoord(sbx, sby, sbz, x, y, z, &rx, &ry, &rz);

        if (!bp->am)
        {
            double axis = ((bp->gdimz == 1) ? ry : rz) - bp->cellSize * 0.5;    // 2D grow y direction, 3D grow z direction

            double tgrad = t->gradient + t->grad_coef * sim_time;
            double v = t->velocity + t->velo_coef * sim_time;
            double vt = 0.5 * (v + t->velocity) * sim_time;
            //   double iso_coef2 = t->iso_coef2 + t->coef_iso2 * sim_time;
            //   double slope = t->grad_slope + t->slope_coef * sim_time;
            //   double T =
            //       t->initialTemperature + tgrad * (axis -
            //                                        vt) * cos(slope * M_PI / 180.0) +
            //       rad * sin(slope * M_PI / 180.0) * tgrad;
            //   T -= tgrad * (t->iso_coef1 * rad + iso_coef2 * rad * rad);
            T = t->initialTemperature + tgrad * (axis - vt);
        }
        else if (bp->thermal_ana == 1)
        {
            double lv = bp->lv;
            double ls = bp->lsd / 2 / 1.414;

            double xlen = bp->gnsbx * bp->gsdimx * bp->cellSize;        // lenght of x direction
            xlen = xlen + bp->lsd;
            double tline = xlen / lv;   //time for each line in x direction
            double tpause = 1.5e-4 / lv;
            tline = tline + tpause;
            //tline = 3.5e-3;
            double trackno = 0;

            int layno = bp->amlayer;    //number of layers
            double tlayer = 50e-6;      //layer thickness

            if (bp->amtracks > 1)
                trackno = floor(sim_time / tline);

            double xnow = bp->x0 + (sim_time - trackno * tline) * lv;
            double ynow = bp->y0 + trackno * bp->lhs;
            double znow = bp->z0;
            double lpow = bp->lpow;

            double lx = rx - xnow;
            double ly = ry - ynow;
            double lz = rz - znow;

            if (rz <= znow + bp->cellSize)
            {

                int nt = 400;
                double lt;
                double Tinteg = 0.0;
                double Tlocal, Tfront;
                double pi15 = 5.568328;
                double lalpha = bp->lthcon / bp->rho / bp->lcp;
                double tsteady = 20 * lalpha / lv / lv;
                double ldt = tsteady / nt;

                if (sim_time < 1e-4)
                    lpow = 0;   //define sometime for substrate to form grains

                for (int i = 0; i < nt; i++)
                {
                    lt = (0.5 + i) * ldt;
                    Tlocal = sqrt(lalpha * lt) * (ls * ls + 4 * lalpha * lt);
                    Tinteg +=
                        ldt *
                        exp(-((lx + lv * lt) * (lx + lv * lt) + ly * ly) /
                            (ls * ls + 4 * lalpha * lt) -
                            (lz * lz / 4 / lalpha / lt)) / Tlocal;
                }

                Tfront = bp->lab * lpow * lalpha / pi15 / bp->lthcon;
                T = bp->liniT + Tfront * Tinteg;
            }
            else
            {
                T = 0.0;
            }
        }                       // bp->thermal_ana ==1
    }
    return T;
}

#ifdef GPU_OMP
#pragma omp end declare target
#endif

static void
blank_moldface(
    int8_t * mold,
    size_t x0,
    size_t x1,
    size_t y0,
    size_t y1,
    size_t z0,
    size_t z1)
{
    for (size_t z = z0; z < z1; z++)
    {
        for (size_t y = y0; y < y1; y++)
        {
            for (size_t x = x0; x < x1; x++)
            {
                mold[SBIDX(x, y, z)] = 1;
            }
        }
    }
}


static void
tempPrepareSB_int(
    SB_struct * sb)
{
    size_t dimx = bp->gsdimx + 2;
    size_t dimy = bp->gsdimy + 2;
    size_t dimz = bp->gsdimz + 2;

    int totaldim = dimx * dimy * dimz;

    double* temperature = sb->temperature;
    int* lsindex = sb->lsindex;
    uint32_t sbx = sb->coords.x;
    uint32_t sby = sb->coords.y;
    uint32_t sbz = sb->coords.z;

#ifdef GPU_OMP
#pragma omp target teams distribute parallel for collapse(3) schedule(static,1)
#endif
    for (size_t z = 0; z < dimz; z++)
    {
        for (size_t y = 0; y < dimy; y++)
        {
            for (size_t x = 0; x < dimx; x++)
            {
                temperature[SBIDX(x, y, z)] =
                    direct_temp_calc(sbx, sby, sbz, x, y, z);

                if (bp->am)
                {
                    if (temperature[SBIDX(x, y, z)] >= bp->liquidusTemp)
                        lsindex[SBIDX(x, y, z)] = 1;
                    else
                        lsindex[SBIDX(x, y, z)] = 0;
                }
            }
        }
    }

#ifdef GPU_OMP
#pragma omp target update from(temperature[0:totaldim])
#pragma omp target update from(lsindex[0:totaldim])
#endif

    /* Set up the mold in the borderlands */
    if (sb->neighbors[FACE_TOP][0] == -1)
        blank_moldface(sb->mold, 0, bp->gsdimx + 2, 0, bp->gsdimy + 2,
                       bp->gsdimz + 1, bp->gsdimz + 2);
    if (sb->neighbors[FACE_BOTTOM][0] == -1)
        blank_moldface(sb->mold, 0, bp->gsdimx + 2, 0, bp->gsdimy + 2, 0, 1);
    if (sb->neighbors[FACE_LEFT][0] == -1)
        blank_moldface(sb->mold, 0, 1, 0, bp->gsdimy + 2, 0, bp->gsdimz + 2);
    if (sb->neighbors[FACE_RIGHT][0] == -1)
        blank_moldface(sb->mold, bp->gsdimx + 1, bp->gsdimx + 2, 0,
                       bp->gsdimy + 2, 0, bp->gsdimz + 2);
    if (sb->neighbors[FACE_FRONT][0] == -1)
        blank_moldface(sb->mold, 0, bp->gsdimx + 2, 0, 1, 0, bp->gsdimz + 2);
    if (sb->neighbors[FACE_BACK][0] == -1)
        blank_moldface(sb->mold, 0, bp->gsdimx + 2, bp->gsdimy + 1,
                       bp->gsdimy + 2, 0, bp->gsdimz + 2);
}


static void
tempUpdateSB_int(
    SB_struct * sb,
    void *vUpdateAll)
{
    assert(sb != NULL);

    int dimx = bp->gsdimx + 2;
    int dimy = bp->gsdimy + 2;
    int dimz = bp->gsdimz + 2;

    double *temperature = sb->temperature;
    int* lsindex = sb->lsindex;

    uint32_t sbx = sb->coords.x;
    uint32_t sby = sb->coords.y;
    uint32_t sbz = sb->coords.z;

#ifdef GPU_OMP

    int totaldim = dimx * dimy * dimz;
#pragma omp target update to(bp[0:1])
#pragma omp target teams distribute parallel for collapse(3) schedule(static,1)
#endif
    for (int z = 0; z < dimz; z++)
    {
        for (int y = 0; y < dimy; y++)
        {
            for (int x = 0; x < dimx; x++)
            {
                /* Note:  This function cannot auto-vectorize in GCC 4.7.0
                 * because __builtin_sincos() uses a double-complex exponent
                 * for it's calculations, and double-complex cannot, at this
                 * time, be vectorized.
                 */
                temperature[SBIDX(x, y, z)] =
                    direct_temp_calc(sbx, sby, sbz, x, y, z);

                if ((bp->thermal_ana)
                    && (temperature[SBIDX(x, y, z)] > bp->liquidusTemp))
                    lsindex[SBIDX(x, y, z)] = 1;
            }
        }
    }

#ifdef GPU_OMP
#ifndef GPU_OMP_NUC
#pragma omp target update from(temperature[0:totaldim])
#endif
#pragma omp target update from(lsindex[0:totaldim])
#endif

}

static void
tempUpdate_int(
    bool updateAll)
{
    tempUpdateSB_int(lsp, &updateAll);
}

static uint64_t
timeToReachTemp_int(
    volume_t * v,
    double target)
{
#if 0
    double t_time =
        (bp->temp_ctrl.initialTemperature - target) / fabs(bp->coolrate);
    if (t_time < 0)
        return 0;
    return (uint64_t) (t_time / bp->ts_delt);
#else
    return 0;
#endif
}


static double
totalNonMoldVolume_int(
    )
{
    return (bp->gdimx * bp->gdimy * bp->gdimz) * pow(bp->cellSize, 3);
}


/****** PUBLIC FUNCTIONS ******/

void
tempPrepareSB(
    SB_struct * sb)
{
    switch (bp->temp_type)
    {
        case INTERNAL:
            tempPrepareSB_int(sb);
            break;
    }
    return;
}

void
tempUpdate(
    bool updateAll)
{
    switch (bp->temp_type)
    {
        case INTERNAL:
            tempUpdate_int(updateAll);
            break;
    }
    return;
}

uint64_t
timeToReachTemp(
    volume_t * volume,
    double target)
{
    switch (bp->temp_type)
    {
        case INTERNAL:
            return timeToReachTemp_int(volume, target);
    }
    return TS_NEVER;
}


uint64_t
timeToReachTempSBID(
    size_t sbid,
    double target)
{
    volume_t vol;
    getSBIDRealBoundsVolume(sbid, &vol);
    return timeToReachTemp(&vol, target);
}


double
totalNonMoldVolume(
    )
{
    switch (bp->temp_type)
    {
        case INTERNAL:
            return totalNonMoldVolume_int();
    }
    return 1.0;
}
