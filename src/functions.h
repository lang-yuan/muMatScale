/***************************************************************/
/* Copyright (c) 2023, Lang Yuan, Univeristy of South Carolina */
/* All rights reserved.                                        */
/* This file is part of muMatScale.                            */
/* See the top-level LICENSE file for details.                 */
/***************************************************************/

#ifndef __FUNCTIONS_H__
#define __FUNCTIONS_H__

#include "globals.h"

/**
 * Calculates the subblock IDs of the neighboring subblocks
 * \param[in] id Subblock ID of the target block
 * \param[out] neighbors Output array of length NUM_NEIGHBORS with associated id numbers
 */
void determine_3dneighbors(
    int id,
    int neighbors[NUM_NEIGHBORS]);

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
void scoord2globalcellcoord(
    uint32_t sx,
    uint32_t sy,
    uint32_t sz,
    int i,
    int j,
    int k,
    int *x,
    int *y,
    int *z);

/**
 * Translates global coordinates to global cell indicies
 * \param[in] x coordinate
 * \param[in] y coordinate
 * \param[in] z coordinate
 * \param[out] i global index
 * \param[out] j global index
 * \param[out] k global index
 */
void realcoord2scoord(
    double x,
    double y,
    double z,
    size_t * i,
    size_t * j,
    size_t * k);

/**
 * Translates inside-of-subblock indicies to global cell coordinates
 * \param[in] s Reference Subblock
 * \param[in] i index
 * \param[in] j index
 * \param[in] k index
 * \param[out] x coordinate
 * \param[out] y coordinate
 * \param[out] z coordinate
 */
void scoord2realcoord(
    uint32_t sx,
    uint32_t sy,
    uint32_t sz,
    int i,
    int j,
    int k,
    double *x,
    double *y,
    double *z);


/**
 * Calculates the coordinate boundaries of a subblock
 * \param[in] sbx,sby,sbz Reference Subblock
 * \param[out] lx Low x coordinate
 * \param[out] ly Low y coordinate
 * \param[out] lz Low z coordinate
 * \param[out] ux High x coordinate
 * \param[out] uy High y coordinate
 * \param[out] uz High z coordinate
 */
void getSBRealBounds(
    uint32_t sbx,
    uint32_t sby,
    uint32_t sbz,
    double *lx,
    double *ly,
    double *lz,
    double *ux,
    double *uy,
    double *uz);

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
void getSBIDRealBounds(
    size_t sbid,
    double *lx,
    double *ly,
    double *lz,
    double *ux,
    double *uy,
    double *uz);

/**
 * Calculates the coordinate boundaries of a subblock, given an ID
 * \param[in] sbid Reference subblock ID
 * \param[out] v   Volume description of boundaries
 */
void getSBIDRealBoundsVolume(
    size_t sbid,
    volume_t * volume);

/**
 * Calculates the coordinate boundaries of a subblock
 * \param[in]  s Reference subblock
 * \param[out] v Volume description of boundaries
 */
void getSBRealBoundsVolume(
    SB_struct * s,
    volume_t * volume);

/**
 * Calculates the coordinate boundaries of a subblock
 * \param[in]  s Reference subblock
 * \param[out] v Volume description of boundaries
 */
void getMSBRealBoundsVolume(
    MSB_struct * s,
    volume_t * volume);


/**
 * Calculate subblock ID from subblock indicies
 * \param[in] x subblock index
 * \param[in] y subblock index
 * \param[in] z subblock index
 * \return subblock ID
 */
int sbID_from_sbcoords(
    int x,
    int y,
    int z);

/**
 * Calculates subblock ID from real coordinates
 * \param[in] x coordinate
 * \param[in] y coordinate
 * \param[in] z coordinate
 * \return subblock ID
 */
int sbID_from_realcoords(
    double x,
    double y,
    double z);

/**
 * Calculates subblock indicies from subblock ID
 * \param[in] sbid Reference subblock ID
 * \param[out] x subblock index
 * \param[out] y subblock index
 * \param[out] z subblock index
 */
void sbcoords_from_sbID(
    size_t sbid,
    size_t * x,
    size_t * y,
    size_t * z);

void create_MPIdtypes(
    );
void destroy_MPIdtypes(
    );

/**
 * Sets up MPI datatypes used
 */
void setup_mpi_datatypes(
    void);



/**
 * \return the next pseudorandom "normally" distributed double value
 * with a mean of 0.0 and a standard deviation of 1.0
 */
float getRandomGaussian(
    void);

/**
 * Gets a random number with "normal" distribution in a range
 * \param[in] min Minimum value for the range
 * \param[in] max Maximum value for the range
 * \return Random number
 */
float getRandomGaussianInRange(
    float min,
    float max);

/**
 * \return a uniformly distributed random number in the range 0 to 1
 */
float getRandScale(
    void);







/**
 * Allocates memory, an array
 * of type \c double with a size equal to a subblock, including HALO
 */
void allocate_float(
    double **f);

/**
 * Allocates memory, an array
 * of type \c int with a size equal to a subblock, including HALO
 */
void allocate_int(
    int **f);

/**
 * Allocates memory, an array
 * of type \c int8_t with a size equal to a subblock, including HALO
 */
void allocate_decentered(
    decentered_t ** f);

/**
 * Allocates memory, an array
 * of type \c decentered_t with a size equal to a subblock, including HALO
 */
void allocate_byte(
    int8_t ** f);



/**
 * Called during initialization by all ranks, uses the temperature
 * model to find the fast-forward time.
 */
void findStartTime(
    void);


/**
 * Check to see if using an external temperature model such as PROCAST
 *
 * \return True if using an External Temperature model
 */


#if 0                           /* useful for debugging */
#include <execinfo.h>
#define dump_stack_trace() \
	do { \
		int _depth = 16; \
		void *_btarray[16]; \
		fprintf(stderr, "BACKTRACE FROM %s::%s:%d\n", __FUNCTION__, __FILE__, __LINE__); \
		_depth = backtrace(_btarray, _depth); \
		backtrace_symbols_fd(_btarray, _depth, 2); \
	} while (0)
#endif

// Translates global cell indices to global cell coordinates
#define CELLCOORD2REALCOORD0(i) (bp->origin_offset[0]+ (i+0.5) * bp->cellSize)
#define CELLCOORD2REALCOORD1(i) (bp->origin_offset[1]+ (i+0.5) * bp->cellSize)
#define CELLCOORD2REALCOORD2(i) (bp->origin_offset[2]+ (i+0.5) * bp->cellSize)

#endif
