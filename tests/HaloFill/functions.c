
/***************************************************************/
/* Yuan, Univeristy of South Carolina		               */
/* 					                       */
/* All rights reserved.                                        */
/***************************************************************/

#include "globals.h"
#include "face_util.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#define CELLCOORD2REALCOORD0(i) (bp->origin_offset[0]+ (i+0.5) * bp->cellSize)
#define CELLCOORD2REALCOORD1(i) (bp->origin_offset[1]+ (i+0.5) * bp->cellSize)
#define CELLCOORD2REALCOORD2(i) (bp->origin_offset[2]+ (i+0.5) * bp->cellSize)

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
 * Translates inside-of-subblock indicies to global cell coordinates
 * \param[in] sbx,sby,sbz Reference Subblock
 * \param[in] i index
 * \param[in] j index
 * \param[in] k index
 * \param[out] x coordinate
 * \param[out] y coordinate
 * \param[out] z coordinate
 */

#ifdef GPU_OMP
#pragma omp declare target
#endif

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

#ifdef GPU_OMP
#pragma omp end declare target
#endif

