/***************************************************************/
/* Copyright (c) 2023, Lang Yuan, Univeristy of South Carolina */
/* All rights reserved.                                        */
/* This file is part of muMatScale.                            */
/* See the top-level LICENSE file for details.                 */
/***************************************************************/
#include "globals.h"
#include "file_io.h"
#include "temperature.h"

#include "io_vtk.h"
#include "io_xdmf.h"


struct io_switch
{
    const char *name;
    void (
    *prepareIO) (
    void);
    void (
    *writeMain) (
    void);
    void (
    *prepareIOSubblock) (
    SB_struct * sb);
    void (
    *writeSubblocks) (
    void);
    void (
    *closeIO) (
    void);
};


static struct io_switch g_format_types[] = {
    {"NONE", NULL, NULL, NULL, NULL, NULL},     // NONE type
    {"VTK", NULL, vtk_writeMain, vtk_prepareIOSubblock, vtk_writeSubblocks, NULL},      // VTK
    {"XDMF", xdmf_prepareIO, xdmf_writeMain, NULL, xdmf_writeSubblocks,
     xdmf_closeIO},
};


/**
 * Returns the name of the type of IO that will be done
 * \return a null-terminated string with the type of File-IO to be done
 */
const char *
getIOTypeName(
    )
{
    return g_format_types[bp->vis_format].name;
}


/**
 * Prepares the IO subsystem
 * Called by both Main and Tasks
 */
void
prepareIO(
    )
{
    if (g_format_types[bp->vis_format].prepareIO != NULL)
        g_format_types[bp->vis_format].prepareIO();
}


/**
 * Writes out any data that the Main needs to write
 */
void
writeMain(
    )
{
    if (g_format_types[bp->vis_format].writeMain != NULL)
        g_format_types[bp->vis_format].writeMain();
}


/**
 * Prepares the IO subsystem for a new Subblock
 */
void
prepareIOSubblock(
    SB_struct * sb)
{
    if (g_format_types[bp->vis_format].prepareIOSubblock != NULL)
        g_format_types[bp->vis_format].prepareIOSubblock(sb);
}


/**
 * Writes out any data that the Tasks need to write for their subblocks
 */
void
writeSubblocks(
    )
{
    tempUpdate(true);           // make sure temperature is up to date for ALL subblocks

    if (g_format_types[bp->vis_format].writeSubblocks != NULL)
        g_format_types[bp->vis_format].writeSubblocks();
}

/**
 * Finialize any IO for the simulation.
 */
void
closeIO(
    )
{
    if (g_format_types[bp->vis_format].closeIO != NULL)
        g_format_types[bp->vis_format].closeIO();
}
