/***************************************************************/
/* Copyright (c) 2023, Lang Yuan, Univeristy of South Carolina */
/* All rights reserved.                                        */
/* This file is part of muMatScale.                            */
/* See the top-level LICENSE file for details.                 */
/***************************************************************/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "globals.h"
#include "debug.h"
#include "functions.h"
#include "temperature.h"
#include "face_util.h"
#include "xmalloc.h"


/**
 * Calculates the subblock IDs of the neighboring subblocks
 * \param[in] id Subblock ID of the target block
 * \param[out] neighbors Output array of length NUM_NEIGHBORS with associated id numbers
 */
void
determine_3dneighbors(
    int id,
    int neighbors[NUM_NEIGHBORS])
{
    assert(id >= 0);
    assert(id < nproc);

    int nsbplane = bp->gnsbx * bp->gnsby;
    int nsbrow = bp->gnsbx;
    int my_x = id % nsbrow;
    int my_y = (id / nsbrow) % bp->gnsby;
    int my_z = id / nsbplane;

    neighbors[FACE_TOP] = sbID_from_sbcoords(my_x, my_y, my_z + 1);
    neighbors[FACE_BOTTOM] = sbID_from_sbcoords(my_x, my_y, my_z - 1);
    neighbors[FACE_BACK] = sbID_from_sbcoords(my_x, my_y + 1, my_z);
    neighbors[FACE_FRONT] = sbID_from_sbcoords(my_x, my_y - 1, my_z);
    neighbors[FACE_RIGHT] = sbID_from_sbcoords(my_x + 1, my_y, my_z);
    neighbors[FACE_LEFT] = sbID_from_sbcoords(my_x - 1, my_y, my_z);
}


/**
 * Translates inside-of-subblock indicies to global cell indicies
 * \param[in] s Reference Subblock
 * \param[in] i index
 * \param[in] j index
 * \param[in] k index
 * \param[out] x global index
 * \param[out] y global index
 * \param[out] z global index
 */

#ifdef GPU_OMP
#pragma omp declare target
#endif

void
scoord2globalcellcoord(
    uint32_t sbx,
    uint32_t sby,
    uint32_t sbz,
    int i,
    int j,
    int k,
    int *x,
    int *y,
    int *z)
{
    *x = sbx * bp->gsdimx + i;
    *y = sby * bp->gsdimy + j;
    *z = sbz * bp->gsdimz + k;
}


/**
 * Translates global coordinates to global cell indicies
 * \param[in] x coordinate
 * \param[in] y coordinate
 * \param[in] z coordinate
 * \param[out] i global index
 * \param[out] j global index
 * \param[out] k global index
 */
void
realcoord2scoord(
    double x,
    double y,
    double z,
    size_t * i,
    size_t * j,
    size_t * k)
{
    x = x - bp->origin_offset[0] + bp->cellSize * 0.1;
    y = y - bp->origin_offset[1] + bp->cellSize * 0.1;
    z = z - bp->origin_offset[2] + bp->cellSize * 0.1;

    *i = (size_t) (x / bp->cellSize);
    *j = (size_t) (y / bp->cellSize);
    *k = (size_t) (z / bp->cellSize);
}


/**
 * Translates inside-of-subblock indicies to global cell coordinates
 * \param[in] sbx,sby,sbz Reference Subblock
 * \param[in] i index
 * \param[in] j index
 * \param[in] k index
 * \param[out] x coordinate
 * \param[out] y coordinate
 * \param[out] z coordinate
 */

void
scoord2realcoord(
    uint32_t sbx,
    uint32_t sby,
    uint32_t sbz,
    int i,
    int j,
    int k,
    double *x,
    double *y,
    double *z)
{
    int cx = sbx * bp->gsdimx;
    int cy = sby * bp->gsdimy;
    int cz = sbz * bp->gsdimz;

    *x = CELLCOORD2REALCOORD0(cx + i);
    *y = CELLCOORD2REALCOORD1(cy + j);
    *z = CELLCOORD2REALCOORD2(cz + k);
}

/**
 * Calculates the coordinate boundaries of a subblock
 * \param[in] s Reference Subblock
 * \param[out] lx Low x coordinate
 * \param[out] ly Low y coordinate
 * \param[out] lz Low z coordinate
 * \param[out] ux High x coordinate
 * \param[out] uy High y coordinate
 * \param[out] uz High z coordinate
 */
void
getSBRealBounds(
    uint32_t sbx,
    uint32_t sby,
    uint32_t sbz,
    double *lx,
    double *ly,
    double *lz,
    double *ux,
    double *uy,
    double *uz)
{
    scoord2realcoord(sbx, sby, sbz, 0, 0, 0, lx, ly, lz);
    *lx -= bp->cellSize / 2.0;
    *ly -= bp->cellSize / 2.0;
    *lz -= bp->cellSize / 2.0;
    scoord2realcoord(sbx, sby, sbz, bp->gsdimx, bp->gsdimy, bp->gsdimz, ux,
                     uy, uz);
    *ux += bp->cellSize / 2.0;
    *uy += bp->cellSize / 2.0;
    *uz += bp->cellSize / 2.0;
}

#ifdef GPU_OMP
#pragma omp end declare target
#endif
/**
 * Calculates the coordinate boundaries of a subblock, given an ID
 * \param[in] sbid Reference subblock ID
 * \param[out] lx Low x coordinate
 * \param[out] ly Low y coordinate
 * \param[out] lz Low z coordinate
 * \param[out] ux High x coordinate
 * \param[out] uy High y coordinate
 * \param[out] uz High z coordinate
 */
void
getSBIDRealBounds(
    size_t sbid,
    double *lx,
    double *ly,
    double *lz,
    double *ux,
    double *uy,
    double *uz)
{
    size_t x, y, z;
    sbcoords_from_sbID(sbid, &x, &y, &z);
    x *= bp->gsdimx;
    y *= bp->gsdimy;
    z *= bp->gsdimz;
    *lx = CELLCOORD2REALCOORD0(x);
    *ly = CELLCOORD2REALCOORD1(y);
    *lz = CELLCOORD2REALCOORD2(z);
    *lx -= 0.5 * bp->cellSize;
    *ly -= 0.5 * bp->cellSize;
    *lz -= 0.5 * bp->cellSize;
    *ux = CELLCOORD2REALCOORD0(x + bp->gsdimx);
    *uy = CELLCOORD2REALCOORD1(y + bp->gsdimy);
    *uz = CELLCOORD2REALCOORD2(z + bp->gsdimz);
    *ux += 0.5 * bp->cellSize;
    *uy += 0.5 * bp->cellSize;
    *uz += 0.5 * bp->cellSize;
}


/**
 * Calculates the coordinate boundaries of a subblock, given an ID
 * \param[in] sbid Reference subblock ID
 * \param[out] v   Volume description of boundaries
 */
void
getSBIDRealBoundsVolume(
    size_t sbid,
    volume_t * v)
{
    getSBIDRealBounds(sbid, &v->minX, &v->minY, &v->minZ, &v->maxX, &v->maxY,
                      &v->maxZ);
}


/**
 * Calculates the coordinate boundaries of a subblock
 * \param[in]  s Reference subblock
 * \param[out] v Volume description of boundaries
 */
void
getSBRealBoundsVolume(
    SB_struct * s,
    volume_t * v)
{
    getSBRealBounds(s->coords.x, s->coords.y, s->coords.z, &v->minX, &v->minY,
                    &v->minZ, &v->maxX, &v->maxY, &v->maxZ);
}


/**
 * Calculates the coordinate boundaries of a subblock
 * \param[in]  s Reference subblock
 * \param[out] v Volume description of boundaries
 */
void
getMSBRealBoundsVolume(
    MSB_struct * s,
    volume_t * v)
{
    SB_struct *sb = (SB_struct *) s;
    getSBRealBounds(sb->coords.x, sb->coords.y, sb->coords.z, &v->minX,
                    &v->minY, &v->minZ, &v->maxX, &v->maxY, &v->maxZ);
}


/**
 * Calculate subblock ID from subblock indicies
 * \param[in] x subblock index
 * \param[in] y subblock index
 * \param[in] z subblock index
 * \return subblock ID
 */
int
sbID_from_sbcoords(
    int x,
    int y,
    int z)
{
    if (bp->padBndyX)
    {
        if (x >= bp->gnsbx)
            return -1;
        if (x < 0)
            return -1;
    }
    else
    {
        x = (x + bp->gnsbx) % bp->gnsbx;        // Implicity deal with negatives
    }
    if (bp->padBndyY)
    {
        if (y >= bp->gnsby)
            return -1;
        if (y < 0)
            return -1;
    }
    else
    {
        y = (y + bp->gnsby) % bp->gnsby;        // Implicity deal with negatives
    }
    if (bp->padBndyZ)
    {
        if (z >= bp->gnsbz)
            return -1;
        if (z < 0)
            return -1;
    }
    else
    {
        z = (z + bp->gnsbz) % bp->gnsbz;        // Implicity deal with negatives
    }
    return x + y * bp->gnsbx + z * bp->gnsbx * bp->gnsby;
}

/**
 * Calculates subblock indicies from subblock ID
 * \param[in] sbid Reference subblock ID
 * \param[out] x subblock index
 * \param[out] y subblock index
 * \param[out] z subblock index
 */
void
sbcoords_from_sbID(
    size_t sbid,
    size_t * x,
    size_t * y,
    size_t * z)
{
    *z = sbid / (bp->gnsbx * bp->gnsby);
    *y = (sbid / bp->gnsbx) % bp->gnsby;
    *x = sbid % bp->gnsbx;
}


/**
 * Calculates subblock ID from real coordinates
 * \param[in] x coordinate
 * \param[in] y coordinate
 * \param[in] z coordinate
 * \return subblock ID
 */
int
sbID_from_realcoords(
    double x,
    double y,
    double z)
{
    size_t ix, iy, iz;
    realcoord2scoord(x, y, z, &ix, &iy, &iz);
    ix = ix / bp->gsdimx;
    iy = iy / bp->gsdimy;
    iz = iz / bp->gsdimz;
    return sbID_from_sbcoords(ix, iy, iz);
}

void
create_MPIdtypes(
    )
{
    createCommTypes(MPI_INT, &bp->intCommTypes);
    createCommTypes(MPI_DOUBLE, &bp->floatCommTypes);
    MPI_Type_vector(1, 3, 1, MPI_DOUBLE, &bp->MPI_DECENTERED);
    MPI_Type_commit(&bp->MPI_DECENTERED);
    createCommTypes(bp->MPI_DECENTERED, &bp->decenteredCommTypes);
}

void
destroy_MPIdtypes(
    )
{
    MPI_Type_free(&bp->MPI_DECENTERED);
}

/**
 * Sets up MPI datatypes used
 */
void
setup_mpi_datatypes(
    )
{
    // MPI_Type_Neighbors = nbr_info_t
    const int nbrNum = 2;
    MPI_Datatype nbrTypes[] = { MPI_INT, MPI_INT };
    int nbrBlockLen[] = { 1, NUM_NEIGHBORS * NBOR_INFO };
    MPI_Aint nbrDisp[] = { 0, offsetof(nbr_info_t, nbors) };
    MPI_Type_create_struct(nbrNum, nbrBlockLen, nbrDisp, nbrTypes,
                           &MPI_Type_Neighbor);
    MPI_Type_commit(&MPI_Type_Neighbor);
}


/**
 * \return the next pseudorandom "normally" distributed double value
 * with a mean of 0.0 and a standard deviation of 1.0
 */
float
getRandomGaussian(
    void)
{
    static int have_one_stored = false;
    static float next_val = 0.0;

    if (have_one_stored)
    {
        have_one_stored = false;
        return next_val;
    }
    else
    {
        float val;
        float n0, n1, s;
        do
        {
            // Random double between -1 and 1.0
            n0 = ((float) rand() / (float) (RAND_MAX / 2.0)) - 1.0;
            n1 = ((float) rand() / (float) (RAND_MAX / 2.0)) - 1.0;
            s = n0 * n0 + n1 * n1;
        }
        while (s >= 1.0 || s == 0);
        float multiplier = sqrtf(-2.0 * log(s) / s);
        val = n0 * multiplier;
        next_val = n1 * multiplier;
        have_one_stored = true;
        return val;
    }
}


/**
 * Gets a random number with "normal" distribution in a range
 * \param[in] min Minimum value for the range
 * \param[in] max Maximum value for the range
 * \return Random number
 */
float
getRandomGaussianInRange(
    float min,
    float max)
{
    if (min > max)
        error("min %f > max %f\n", min, max);

    float scale = max - min;

    float val;
    do
    {
        // Assumption, 95% of the time, this will fire once, as we're
        // scaling withing 2 sigma of the mean (zero).
        val = (getRandomGaussian() + 1.0) / 4.0;        // Move from [-2 2] to [0 1]
    }
    while (val < 0 || val > 1);
    return min + scale * val;
}


/**
 * \return a uniformly distributed random number in the range 0 to 1
 */
float
getRandScale(
    )
{
    return (float) rand() / (float) RAND_MAX;
}

/**
 * Allocates memory, and tags for automatic freeing, an array
 * of type \c double with a size equal to a subblock, including HALO
 */
void
allocate_float(
    double **f)
{
    xmalloc(*f, double,
              (bp->gsdimx + 2) * (bp->gsdimy + 2) * (bp->gsdimz + 2));
}

/**
 * Allocates memory, and tags for automatic freeing, an array
 * of type \c int with a size equal to a subblock, including HALO
 */
void
allocate_int(
    int **f)
{
    xmalloc(*f, int,
              (bp->gsdimx + 2) * (bp->gsdimy + 2) * (bp->gsdimz + 2));
}

/**
 * Allocates memory, and tags for automatic freeing, an array
 * of type \c int8_t with a size equal to a subblock, including HALO
 */
void
allocate_byte(
    int8_t ** f)
{
    xmalloc(*f, int8_t,
            (bp->gsdimx + 2) * (bp->gsdimy + 2) * (bp->gsdimz + 2));
}

/**
 * Allocates memory, and tags for automatic freeing, an array
 * of type \c decentered_t with a size equal to a subblock, including HALO
 */
void
allocate_decentered(
    decentered_t ** f)
{
    xmalloc(*f, decentered_t,
            (bp->gsdimx + 2) * (bp->gsdimy + 2) * (bp->gsdimz + 2));
}


/**
 * Called during initialization by all ranks, uses the temperature
 * model to find the fast-forward time.
 */
void
findStartTime(
    void)
{
    uint64_t min_active_ts = TS_NEVER;
    size_t gnsb = bp->gnsbx * bp->gnsby * bp->gnsbz;
    uint64_t *timestamps;
    int *recvcnts;
    int *displacements;
    xmalloc(timestamps, uint64_t, gnsb);
    xmalloc(recvcnts, int,
            nproc);
    xmalloc(displacements, int,
            nproc);

    recvcnts[0] = gnsb / nproc;
    displacements[0] = 0;
    for (int i = 1; i < nproc; i++)
    {
        recvcnts[i] = gnsb / nproc;
        if (i <= (gnsb % nproc))
            recvcnts[i]++;
        displacements[i] = displacements[i - 1] + recvcnts[i - 1];
    }

#pragma omp parallel for
    for (int count = 0; count < recvcnts[iproc]; count++)
    {
        int64_t id = count + displacements[iproc];
        timestamps[id] = timeToReachTempSBID(id, bp->liquidusTemp);
    }

    MPI_Allgatherv(MPI_IN_PLACE, recvcnts[iproc], MPI_LONG_LONG_INT,
                   timestamps, recvcnts, displacements,
                   MPI_LONG_LONG_INT, mpi_comm_new);

    if (iproc == 0)
    {
        size_t mayactivate = 0;
        for (size_t id = 0; id < gnsb; id++)
        {
            gmsp[id].activate_ts = timestamps[id];
            min_active_ts = MIN(min_active_ts, gmsp[id].activate_ts);
            if (gmsp[id].activate_ts != TS_NEVER)
            {
                mayactivate++;
            }
        }
        printf("# of subblocks we expect to activate:  %zu\n", mayactivate);
    }
    for (size_t id = 0; id < gnsb; id++)
    {
        min_active_ts = MIN(min_active_ts, timestamps[id]);
    }
    xfree(timestamps);
    xfree(recvcnts);
    xfree(displacements);

    if (min_active_ts == TS_NEVER)
    {
        error("This simulation will never reach Liquidus Temperature!\n");
    }

    bp->timestep = min_active_ts;

}
