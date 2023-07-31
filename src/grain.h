/***************************************************************/
/* Copyright (c) 2023, Lang Yuan, Univeristy of South Carolina */
/* All rights reserved.                                        */
/* This file is part of muMatScale.                            */
/* See the top-level LICENSE file for details.                 */
/***************************************************************/
#ifndef __GRAIN_H__
#define __GRAIN_H__

#include "globals.h"

/**
 * Initialize the grain module.
 * Called by all tasks during initialization
 */
void grainSetup(
    void);

/**
 * Deinitialize the grain module.
 * Called by all tasks during shutdown
 */
void grainShutdown(
    void);

/**
 * Organizes and activates any new grains, communicating between
 * all tasks
 */
void activateNewGrains(
    void);



/* Task 0 Only Functions */

/**
 * Takes a coordinate and orientation angle to specify a nucleation point.
 * Adds this information to a list to later be created as a nucleation point.
 * Called by the main during configuration reading.
 *
 * \see init_nucleation_sites
 *
 * \param[in] x  global coordinate
 * \param[in] y  global coordinate
 * \param[in] z  global coordinate
 * \param[in] ax Rotation Angle
 * \param[in] ay Rotation Angle
 * \param[in] az Rotation Angle
 */
void add_fixed_nuc(
    double x,
    double y,
    double z,
    double ax,
    double ay,
    double az,
    double tsh);

/**
 * Specifies a volume in which to seed nucleation points.
 * Stores this information for later processing.
 * \see init_nucleation_sites
 * \param[in] density  The density (as a percentage [0 to 1]) of nucleation sites in the volume
 * \param[in] minx     Minimum X coordinate for the volume
 * \param[in] miny     Minimum Y coordinate for the volume
 * \param[in] minz     Minimum Z coordinate for the volume
 * \param[in] maxx     Maximum X coordinate for the volume
 * \param[in] maxy     Maximum Y coordinate for the volume
 * \param[in] maxz     Maximum Z coordinate for the volume
 */
void add_nuc_volume(
    double density,
    double minx,
    double miny,
    double minz,
    double maxx,
    double maxy,
    double maxz,
    double tsh);

/**
 * Initializes all the nucleation sites.
 *
 * There are 3 ways to specify nucleation sites:
 *  \li FixedNuc: Fixed Nucleation Sites
 *  \li NucVolume: Random fill inside specific volume
 *  \li "Purely" Random 'GNGaussian' in the body or on the surface
 *
 * Each potential nucleation site will have an undercooling threshold, that
 * when reached, will cause a grain to form at that site (unless one has
 * already grown into the cell).
 *
 * This function takes all the specified sites and generates their info.
 *
 * \see add_fixed_nuc
 * \see add_nuc_volume
 */
void init_nucleation_sites(
    void);

/**
 * Creates a CSV file with information about every grain.
 * Called by all tasks
 */
void output_grains(
    SB_struct * sb);



/**
 * Called when a new subblock is being activated, this function fills in
 * the potential nucleation sites that reside within this subblock's
 * boundaries.
 *
 * Will assign both fixed nuclei and gaussian-distributed surface and
 * volume nucleation points.
 *
 * \param[in] sb  The subblock
 */
void init_sb_nucleation(
    SB_struct * sb);

/**
 * Given a grain number, retrieve the grain information.
 *
 * \param[in] grnum  Grain Identifying number
 * \return pointer to a grain structure
 */
grain_t *getGrain2(
    int grnum,
    BB_struct * bp);

/**
 * Determines if a cell is past its nucleation undercooling limit.
 * If so, will call for a new grain to be created.
 *
 * Designed to be called from a (potentially parallel \c ll_walk()
 *
 * \param[in] vsb  Void pointer to a Subblock
 */
void cell_nucleation(
    SB_struct * sb,
    void *__unused);

#endif /* __GRAIN_H__ */
