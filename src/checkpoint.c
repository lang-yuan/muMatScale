/***************************************************************/
/* Copyright (c) 2023, Lang Yuan, Univeristy of South Carolina */
/* All rights reserved.                                        */
/* This file is part of muMatScale.                            */
/* See the top-level LICENSE file for details.                 */
/***************************************************************/
#define H5_USE_110_API

#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h>

#include "globals.h"
#include "checkpoint.h"
#include "ll.h"
#include "xmalloc.h"
#include "functions.h"
#include "file_io.h"
#include "temperature.h"
#include "grain.h"
#include "distribute.h"

#ifndef PATH_MAX
#define PATH_MAX 512
#endif

extern ll_t *gsubblock_list;
extern SB_struct *lsp;
tasks_t gtasks;

static const char *
add_subblock(
    SB_struct * sb)
{
    /*
     * INPUT: sb: a subblock
     * OUTPUT: name of the new group based on subblock coordinate, used in writeTaskCheckpoint function.
     */
    static char path[PATH_MAX] = { 0 };

    if (sb != NULL)
        snprintf(path, PATH_MAX - 1, "SubBlock_%dx%dx%d",
                 sb->coords.x, sb->coords.y, sb->coords.z);
    return path;
}



static const char *
read_subblock_cp_path(
    const char *name,
    const char *var)
{
    /*
     * INPUT: name: the name of the group
     * var: the name of the dataset
     * OUTPUT: the path of the dataset in the group, used in  restoreSubblock function
     */
    static char path[PATH_MAX] = { 0 };
    snprintf(path, PATH_MAX - 1, "%s/%s", name, var);
    return path;

}

static const char *
get_subblock_cp_path(
    SB_struct * sb,
    int rank,
    const char *var)
{
    /*
     * INPUT: sb: a subblock
     * rank: the rank of the process
     * var: name of the dataset
     * OUTPUT: the path of the dataset corresponding to the subblock coordinate,
     * used in writeMainCheckpoint, writeTaskCheckpoint, mainRestart, taskRestart
     */
    static char path[PATH_MAX] = { 0 };
    if (!rank)
    {
        snprintf(path, PATH_MAX - 1, "/%s", var);
    }
    else
    {
        if (sb != NULL)
            snprintf(path, PATH_MAX - 1, "/SubBlock_%dx%dx%d/%s",
                     sb->coords.x, sb->coords.y, sb->coords.z, var);
        else
        {
            snprintf(path, PATH_MAX - 1, "/%s", var);

        }
    }
    return path;
}

static const char *
get_rank_dir(
    int rank)
{
    /*
     * INPUT: rank: the rank of the process
     * OUTPUT: the DIRECTORY path, used in get_rank_cp_path
     */
    static char path[PATH_MAX] = { 0 };
    snprintf(path, PATH_MAX - 1, "%s_r%d", bp->basefilename, rank);
    return path;
}

static const char *
get_rank_cp_path(
    int rank,
    uint64_t ts)
{
    /*
     * INPUT: rank: the rank of the process
     * ts: checkpoint timestep
     * OUTPUT: the location of the checkpoint file, used in writeMainCheckpoint,
     * writeTaskCheckpoint, mainRestart, taskRestart, clean_checkpoint
     * NOTE: this function is to get the location of checkpoint file for MULTIPLE restart
     */
    static char path[PATH_MAX] = { 0 };
    if (!rank)
    {
        snprintf(path, PATH_MAX - 1, "cp_main_%05lu.h5", ts);
    }
    else
    {
        snprintf(path, PATH_MAX - 1, "%s/cp_task_%05lu.h5",
                 get_rank_dir(rank), ts);

    }

    return path;
}

static const char *
get_rank_cp_file(
    int rank)
{
    /*
     * INPUT: rank: the rank of the process
     * OUTPUT: the location of the checkpoint file, used in clean_checkpoint
     * NOTE: this function is to get the location of checkpoint file for ONE restart location
     */

    static char path[PATH_MAX] = { 0 };
    if (!rank)
    {
        snprintf(path, PATH_MAX - 1, "cp_main.h5");
    }
    else
    {
        snprintf(path, PATH_MAX - 1, "cp_task_%03d.h5", rank);

    }

    return path;
}


void
clean_checkpoint(
    )
{
    /*
     * INPUT: none
     * OUTPUT: checkpoint files for ONE restart location
     * NOTE: first, this function will save only the checkpoint of the last
     * timestep, remove the checkpoint of other timesteps.
     * (DEBUG_CHECKPOINT==1) will turn on debugging, umatgc will save
     * multiple checkpoints, clean_checkpoint will not be used.
     */
    rename(get_rank_cp_path(iproc, bp->timestep), get_rank_cp_file(iproc));
}





static void
prepareTemp(
    SB_struct * sb)
{
    /*
     * INPUT: sb: a subblock
     * OUTPUT: setup initial conditions for the subblock
     * USE: prepareIOSubblock, tempPrepareSB, getSBrealBounds, init_sb_nucleation
     */
    prepareIOSubblock(sb);
    tempPrepareSB(sb);
    double low_x, low_y, low_z;
    double high_x, high_y, high_z;
    getSBRealBounds(sb->coords.x, sb->coords.y, sb->coords.z, &low_x, &low_y,
                    &low_z, &high_x, &high_y, &high_z);
    init_sb_nucleation(sb);

}


static void
restoreSubblock(
    hid_t file_id,
    const char *name)
{
    /*
     * INPUT: file_id: file identifier of the checkpoint file
     * name: the group ID
     * OUTPUT:
     * The new subblock will have the information which is stored
     * in the group ID
     * This function is used in restoreTaskSubList
     * USE:   get_subblock_cp_path, prepareTemp, ll_append
     */


    /*To reduce the IO operations, I collect all of small variables into one large
     * variable. Then we only need to use one IO operation to save all these
     * variables.  */
    typedef struct cp_taskSubblock
    {
        COMMON_SUBBLOCK_DEFS;
        int last_active_ts;     // <- TODO: Remove?
        size_t mold_size;
    } taskSubblock;



    //Create the subblock size
    hsize_t start[3] = { 0, 0, 0 };
    hsize_t sbdims[3];
    sbdims[0] = bp->gsdimz + 2;
    sbdims[1] = bp->gsdimy + 2;
    sbdims[2] = bp->gsdimx + 2;

    //Create the memory space for subblock
    hid_t memspace = H5Screate_simple(3, sbdims, NULL);

    //Select the hyperslab
    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, start, NULL, sbdims, NULL);

    //Create the datatype for decentered_t
    hid_t dendriteType = H5Tcreate(H5T_COMPOUND, sizeof(decentered_t));
    H5Tinsert(dendriteType, "x", HOFFSET(decentered_t, x), H5T_NATIVE_DOUBLE);
    H5Tinsert(dendriteType, "y", HOFFSET(decentered_t, y), H5T_NATIVE_DOUBLE);
    H5Tinsert(dendriteType, "z", HOFFSET(decentered_t, z), H5T_NATIVE_DOUBLE);

    hid_t memtype, memCoord, memActiveState, memNeiBor;

    hsize_t neiBorDims[2] = { NUM_NEIGHBORS, NBOR_INFO };

    //Create the datatype for loc_t
    memCoord = H5Tcreate(H5T_COMPOUND, sizeof(loc_t));
    H5Tinsert(memCoord, "x", HOFFSET(loc_t, x), H5T_NATIVE_UINT32);
    H5Tinsert(memCoord, "y", HOFFSET(loc_t, y), H5T_NATIVE_UINT32);
    H5Tinsert(memCoord, "z", HOFFSET(loc_t, z), H5T_NATIVE_UINT32);

    //Create the datatype for activatestate_t
    activestate_t val;

    memActiveState = H5Tenum_create(H5T_NATIVE_INT);
    val = (activestate_t) INACTIVE;
    H5Tenum_insert(memActiveState, "inactive", &val);
    val = (activestate_t) ACTIVE;
    H5Tenum_insert(memActiveState, "active", &val);

    //Create the datatype for neighbors_t
    memNeiBor = H5Tarray_create(H5T_NATIVE_INT, 2, neiBorDims);

    //Create the datatype for taskSubblock
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(taskSubblock));
    H5Tinsert(memtype, "last_active",
              HOFFSET(taskSubblock, last_active_ts), H5T_NATIVE_INT);
    H5Tinsert(memtype, "coord", HOFFSET(taskSubblock, coords), memCoord);

    H5Tinsert(memtype, "subblockid",
              HOFFSET(taskSubblock, subblockid), H5T_NATIVE_INT);
    H5Tinsert(memtype, "procid",
              HOFFSET(taskSubblock, procid), H5T_NATIVE_INT);

    H5Tinsert(memtype, "neighbor", HOFFSET(taskSubblock, neighbors),
              memNeiBor);
    H5Tinsert(memtype, "nnuc", HOFFSET(taskSubblock, nnuc), H5T_NATIVE_INT);
    H5Tinsert(memtype, "mold_size",
              HOFFSET(taskSubblock, mold_size), H5T_NATIVE_UINT64);

    hid_t nucleationType = H5Tcreate(H5T_COMPOUND, sizeof(nucleation_t));
    H5Tinsert(nucleationType, "x", HOFFSET(nucleation_t, x),
              H5T_NATIVE_DOUBLE);
    H5Tinsert(nucleationType, "y", HOFFSET(nucleation_t, y),
              H5T_NATIVE_DOUBLE);
    H5Tinsert(nucleationType, "z", HOFFSET(nucleation_t, z),
              H5T_NATIVE_DOUBLE);
    H5Tinsert(nucleationType, "cx", HOFFSET(nucleation_t, cx),
              H5T_NATIVE_UINT64);
    H5Tinsert(nucleationType, "cy", HOFFSET(nucleation_t, cy),
              H5T_NATIVE_UINT64);
    H5Tinsert(nucleationType, "cz", HOFFSET(nucleation_t, cz),
              H5T_NATIVE_UINT64);
    H5Tinsert(nucleationType, "gx", HOFFSET(nucleation_t, gx),
              H5T_NATIVE_DOUBLE);
    H5Tinsert(nucleationType, "gy", HOFFSET(nucleation_t, gy),
              H5T_NATIVE_DOUBLE);
    H5Tinsert(nucleationType, "gz", HOFFSET(nucleation_t, gz),
              H5T_NATIVE_DOUBLE);
    H5Tinsert(nucleationType, "threshold", HOFFSET(nucleation_t, threshold),
              H5T_NATIVE_DOUBLE);

    /*
     * ALLOCATING the memory space of subblock
     */
    assert(lsp == NULL);
    SB_struct *sb = lsp;

    xmalloc(sb, SB_struct, 1);
    allocate_byte(&sb->mold);
    allocate_float(&(sb->fs));
    xmalloc(sb->nuc_threshold, float,
              (bp->gsdimx + 2) * (bp->gsdimy + 2) * (bp->gsdimz + 2));
    allocate_float(&(sb->temperature));

    if (bp->calc_type == DIFFUSION)
    {
        allocate_float(&(sb->ce));
        allocate_float(&(sb->oce));
        allocate_float(&(sb->cl));
        allocate_int(&(sb->diff_id));
        allocate_int(&(sb->nuc_id));
        allocate_int(&(sb->nuc_id2));
        // ADD_CURV_LY
        allocate_float(&(sb->curv));
    }

    // Allocate arrays for decentered octahedron information
    allocate_decentered(&(sb->dc));
    allocate_float(&(sb->d));
    allocate_int(&(sb->gr));
    allocate_int(&(sb->ogr));


    taskSubblock staticSub;
    //Read the staticSub in the dataset
    hid_t dataset = H5Dopen(file_id, read_subblock_cp_path(name, "CoordTemp"),
                            H5P_DEFAULT);
    H5Dread(dataset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &staticSub);

    //Copy the staticSub in the dataset to the new subblock "sb"
    sb->last_active_ts = staticSub.last_active_ts;
    sb->procid = staticSub.procid;
    sb->subblockid = staticSub.subblockid;
    sb->nnuc = staticSub.nnuc;
    sb->mold_size = staticSub.mold_size;

    xmalloc(sb->nuc_pt, nucleation_t, sb->nnuc);

    sb->coords.x = staticSub.coords.x;
    sb->coords.y = staticSub.coords.y;
    sb->coords.z = staticSub.coords.z;

    for (int j = 0; j < NUM_NEIGHBORS; j++)
    {
        for (int k = 0; k < NBOR_INFO; k++)
        {
            sb->neighbors[j][k] = staticSub.neighbors[j][k];
        }
    }

    //Set up the initial conditions for the new subblock
    prepareTemp(sb);


#define READ_FIELD(field, type) \
	do { \
		if ( sb->field ) { \
			hid_t dset = H5Dopen(file_id, read_subblock_cp_path(name, #field), H5P_DEFAULT); \
			H5Dread(dset, type, memspace, H5S_ALL, H5P_DEFAULT, sb->field); \
			H5Dclose(dset); \
		} \
	} while(0)

    READ_FIELD(temperature, H5T_NATIVE_DOUBLE);
    READ_FIELD(gr, H5T_NATIVE_INT);
    READ_FIELD(fs, H5T_NATIVE_DOUBLE);
    READ_FIELD(ce, H5T_NATIVE_DOUBLE);
    READ_FIELD(oce, H5T_NATIVE_DOUBLE);
    READ_FIELD(cl, H5T_NATIVE_DOUBLE);
    READ_FIELD(diff_id, H5T_NATIVE_CHAR);
    READ_FIELD(mold, H5T_NATIVE_CHAR);
    READ_FIELD(d, H5T_NATIVE_DOUBLE);
    READ_FIELD(nuc_threshold, H5T_NATIVE_FLOAT);
    READ_FIELD(dc, dendriteType);
    // ADD_CURV_LY
    READ_FIELD(curv, H5T_NATIVE_DOUBLE);

#undef READ_FIELD
    if (sb->nnuc > 0)
    {
        hsize_t len = sb->nnuc;
        xmalloc(sb->nuc_pt, nucleation_t, sb->nnuc);
        hid_t lenhdl = H5Screate_simple(1, &len, NULL);
        if (lenhdl < 0)
            error("Harumph!\n");

        hid_t nuc_pt_hdl =
            H5Dopen(file_id, read_subblock_cp_path(name, "nuc_pt"),
                    H5P_DEFAULT);
        H5Dread(nuc_pt_hdl, nucleationType, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                sb->nuc_pt);
        H5Dclose(nuc_pt_hdl);
        H5Sclose(lenhdl);
    }

    //Release all resources
    H5Dclose(dataset);
    H5Sclose(memspace);
    H5Tclose(memActiveState);
    H5Tclose(memCoord);
    H5Tclose(memNeiBor);
    H5Tclose(memtype);
    H5Tclose(nucleationType);
}

/************************************************************

  Operator function.

 ************************************************************/
static herr_t
restoreTaskSubList(
    hid_t loc_id,
    const char *name,
    const H5L_info_t * info,
    void *operator_data)
{
    /*
     * INPUT: loc_id: file identifier of the checkpoint file
     * name: the group ID
     * OUTPUT: info: type of the group
     * operator_data: NULL in this case
     * NOTE:   The data in the SUBBLOCK will be allocated
     * This function is used in taskRestart
     * USE: restoreSubblock
     */
    H5O_info_t infobuf;

    /*
     * Get type of the object and display its name and type.
     */


    H5Oget_info_by_name(loc_id, name, &infobuf, H5P_DEFAULT);
    switch (infobuf.type)
    {
        case H5O_TYPE_GROUP:
            //Get the data in the group and add it to the subblock
            restoreSubblock(loc_id, name);
            break;
        case H5O_TYPE_DATASET:
            printf("  Dataset: %s\n", name);
            break;
        case H5O_TYPE_NAMED_DATATYPE:
            printf("  Datatype: %s\n", name);
            break;
        default:
            printf("  Unknown: %s\n", name);
    }



    return 0;
}




static void
mainActivateSubblock(
    ll_t * gsubblock_list)
{
    /*
     * INPUT: gsubblock_list: empty gsubblock_list
     * OUTPUT: every elements in gmsp will be added in qsubblock_list of MAIN
     * number of active and deactive subblocks will be updated
     * NOTE: This function is used for main restart
     */
    int dimSize = bp->gnsbx * bp->gnsby * bp->gnsbz;
    for (int i = 0; i < dimSize; i++)
    {
        ll_append(gsubblock_list, &gmsp[i]);
    }

}

static void
mainAddSubblock(
    uint64_t restart,
    hid_t file)
{
    /*
     * INPUT: restart: restart timestep
     * file: file identifier of main checkpoint
     * OUTPUT: gmsp will be restored
     * NOTE: This function is used for main restart
     */

    hid_t memtype, memCoord, memActiveState, memNeiBor, space, dset;
    /* Handles */
    int dimSize = bp->gnsbx * bp->gnsby * bp->gnsbz;
    hsize_t dims[1] = { dimSize }, neiBorDims[2] =
    {
    NUM_NEIGHBORS, NBOR_INFO};

    //create datatype for loc_t
    memCoord = H5Tcreate(H5T_COMPOUND, sizeof(loc_t));
    H5Tinsert(memCoord, "x", HOFFSET(loc_t, x), H5T_NATIVE_UINT32);
    H5Tinsert(memCoord, "y", HOFFSET(loc_t, y), H5T_NATIVE_UINT32);
    H5Tinsert(memCoord, "z", HOFFSET(loc_t, z), H5T_NATIVE_UINT32);

    //create datatype for activestate_t
    activestate_t val;

    memActiveState = H5Tenum_create(H5T_NATIVE_INT);
    val = (activestate_t) INACTIVE;
    H5Tenum_insert(memActiveState, "inactive", &val);
    val = (activestate_t) ACTIVE;
    H5Tenum_insert(memActiveState, "active", &val);

    //create datatype for neighbors_t
    memNeiBor = H5Tarray_create(H5T_NATIVE_INT, 2, neiBorDims);

    hid_t dataspace_pl = H5Pcreate(H5P_DATASET_CREATE);

    H5Pset_chunk(dataspace_pl, 1, dims);
    H5Pset_fill_time(dataspace_pl, H5D_FILL_TIME_NEVER);
    H5Pset_shuffle(dataspace_pl);
    H5Pset_deflate(dataspace_pl, 9);

    // Create the compound datatype for gmsp.
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(MSB_struct));
    H5Tinsert(memtype, "subblockid",
              HOFFSET(MSB_struct, subblockid), H5T_NATIVE_INT);
    H5Tinsert(memtype, "procid", HOFFSET(MSB_struct, procid), H5T_NATIVE_INT);
    H5Tinsert(memtype, "activate_ts",
              HOFFSET(MSB_struct, activate_ts), H5T_NATIVE_UINT64);

    H5Tinsert(memtype, "coord", HOFFSET(MSB_struct, coords), memCoord);

    H5Tinsert(memtype, "neighbor", HOFFSET(MSB_struct, neighbors), memNeiBor);
    H5Tinsert(memtype, "nnuc", HOFFSET(MSB_struct, nnuc), H5T_NATIVE_INT);

    hid_t nucleationType = H5Tcreate(H5T_COMPOUND, sizeof(nucleation_t));
    H5Tinsert(nucleationType, "x", HOFFSET(nucleation_t, x),
              H5T_NATIVE_DOUBLE);
    H5Tinsert(nucleationType, "y", HOFFSET(nucleation_t, y),
              H5T_NATIVE_DOUBLE);
    H5Tinsert(nucleationType, "z", HOFFSET(nucleation_t, z),
              H5T_NATIVE_DOUBLE);
    H5Tinsert(nucleationType, "cx", HOFFSET(nucleation_t, cx),
              H5T_NATIVE_UINT64);
    H5Tinsert(nucleationType, "cy", HOFFSET(nucleation_t, cy),
              H5T_NATIVE_UINT64);
    H5Tinsert(nucleationType, "cz", HOFFSET(nucleation_t, cz),
              H5T_NATIVE_UINT64);
    H5Tinsert(nucleationType, "gx", HOFFSET(nucleation_t, gx),
              H5T_NATIVE_DOUBLE);
    H5Tinsert(nucleationType, "gy", HOFFSET(nucleation_t, gy),
              H5T_NATIVE_DOUBLE);
    H5Tinsert(nucleationType, "gz", HOFFSET(nucleation_t, gz),
              H5T_NATIVE_DOUBLE);
    H5Tinsert(nucleationType, "threshold", HOFFSET(nucleation_t, threshold),
              H5T_NATIVE_DOUBLE);

    //create dataspace for gmsp
    space = H5Screate_simple(1, dims, NULL);

    //read gmsp from main checkpoint file
    dset =
        H5Dopen(file, get_subblock_cp_path(NULL, iproc, "mainSubblock"),
                H5P_DEFAULT);

    H5Dread(dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, gmsp);

    //save nuc_pt
    for (int i = 0; i < dimSize; i++)
    {
        xmalloc(gmsp[i].nuc_pt, nucleation_t, gmsp[i].nnuc);
        gmsp[i].nuc_pt_len = gmsp[i].nnuc;

        char path[PATH_MAX] = { 0 };
        snprintf(path, PATH_MAX - 1, "/nuc_pts/sb_%d", i);

        if (gmsp[i].nnuc > 0)
        {
            hsize_t len = gmsp[i].nnuc;
            hid_t lenhdl = H5Screate_simple(1, &len, NULL);

            hid_t nuc_pt_hdl = H5Dopen(file, path, H5P_DEFAULT);

            H5Dread(nuc_pt_hdl, nucleationType, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                    gmsp[i].nuc_pt);
            H5Dclose(nuc_pt_hdl);
            H5Sclose(lenhdl);
        }
    }

    //release resources
    H5Dclose(dset);
    H5Tclose(memCoord);

    H5Tclose(memActiveState);
    H5Tclose(memNeiBor);

    H5Pclose(dataspace_pl);
    H5Tclose(memtype);
    H5Tclose(nucleationType);
    H5Sclose(space);
}


static void
storeGlobalGrainInfo(
    hid_t file)
{
    hsize_t romatDims[2] = { 3, 3 };
    hsize_t grainDims[1] = { bp->maxTotalGrains + 1 };

    hid_t memRotmat = H5Tarray_create(H5T_NATIVE_DOUBLE, 2, romatDims);

    hid_t memGrain = H5Tcreate(H5T_COMPOUND, sizeof(grain_t));
    H5Tinsert(memGrain, "nuc_x", HOFFSET(grain_t, nuc_x), H5T_NATIVE_DOUBLE);
    H5Tinsert(memGrain, "nuc_y", HOFFSET(grain_t, nuc_y), H5T_NATIVE_DOUBLE);
    H5Tinsert(memGrain, "nuc_z", HOFFSET(grain_t, nuc_z), H5T_NATIVE_DOUBLE);
    H5Tinsert(memGrain, "nuc_timestep",
              HOFFSET(grain_t, nuc_timestep), H5T_NATIVE_UINT64);
    H5Tinsert(memGrain, "rotmat", HOFFSET(grain_t, rotmat), memRotmat);

    hid_t grainSpace = H5Screate_simple(1, grainDims, NULL);

    hid_t dataspace_grain = H5Pcreate(H5P_DATASET_CREATE);

    H5Pset_chunk(dataspace_grain, 1, grainDims);
    H5Pset_fill_time(dataspace_grain, H5D_FILL_TIME_NEVER);
    H5Pset_shuffle(dataspace_grain);
    H5Pset_deflate(dataspace_grain, 9);

    hid_t grainSet = H5Dcreate(file,
                               get_subblock_cp_path(NULL, iproc,
                                                    "grainCache"), memGrain,
                               grainSpace, H5P_DEFAULT, dataspace_grain,
                               H5P_DEFAULT);

    H5Dwrite(grainSet, memGrain, H5S_ALL, H5S_ALL, H5P_DEFAULT,
             grain_cache);
    H5Dclose(grainSet);
    H5Pclose(dataspace_grain);
    H5Sclose(grainSpace);
    H5Tclose(memGrain);
    H5Tclose(memRotmat);


    //save bp->num_grains
    hsize_t bpDims[1] = { 1 };
    hid_t bpSpace = H5Screate_simple(1, bpDims, NULL);
    hid_t numGrain =
        H5Dcreate(file, get_subblock_cp_path(NULL, iproc, "numGrain"),
                  H5T_NATIVE_INT, bpSpace, H5P_DEFAULT, H5P_DEFAULT,
                  H5P_DEFAULT);
    H5Dwrite(numGrain, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
             &bp->num_grains);

    H5Dclose(numGrain);
    H5Sclose(bpSpace);

}


static void
addGrainCache(
    hid_t file)
{
    /*
     * INPUT: file: file identifier of main checkpoint
     * OUTPUT: grain_cached will be restored
     * bp->num_grains will also be restored
     * NOTE: This function is used for main restart
     */


    //int i,j;
    unsigned int i, j;          //[C++]
    hsize_t romatDims[2] = { 3, 3 };
    hsize_t grainDims[1] = { bp->maxTotalGrains + 1 };

    //create datatype for rotmat
    hid_t memRotmat = H5Tarray_create(H5T_NATIVE_DOUBLE, 2, romatDims);

    //create datype for grain_t
    hid_t memGrain = H5Tcreate(H5T_COMPOUND, sizeof(grain_t));
    H5Tinsert(memGrain, "nuc_x", HOFFSET(grain_t, nuc_x), H5T_NATIVE_DOUBLE);
    H5Tinsert(memGrain, "nuc_y", HOFFSET(grain_t, nuc_y), H5T_NATIVE_DOUBLE);
    H5Tinsert(memGrain, "nuc_z", HOFFSET(grain_t, nuc_z), H5T_NATIVE_DOUBLE);
    H5Tinsert(memGrain, "nuc_timestep",
              HOFFSET(grain_t, nuc_timestep), H5T_NATIVE_UINT64);
    H5Tinsert(memGrain, "rotmat", HOFFSET(grain_t, rotmat), memRotmat);

    hid_t grainSpace = H5Screate_simple(1, grainDims, NULL);

    grain_t *cpGrain;

    xmalloc(cpGrain, grain_t, bp->maxTotalGrains + 1);
    //read grain cache from task checkpoint files
    hid_t cp_grain =
        H5Dopen(file, get_subblock_cp_path(NULL, iproc, "grainCache"),
                H5P_DEFAULT);

    H5Dread(cp_grain, memGrain, H5S_ALL, H5S_ALL, H5P_DEFAULT, cpGrain);
    //restore grain_cache
    for (i = 0; i < bp->maxTotalGrains + 1; i++)
    {
        grain_cache[i].nuc_x = cpGrain[i].nuc_x;
        grain_cache[i].nuc_y = cpGrain[i].nuc_y;
        grain_cache[i].nuc_z = cpGrain[i].nuc_z;
        grain_cache[i].nuc_timestep = cpGrain[i].nuc_timestep;

        for (j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                grain_cache[i].rotmat[j][k] = cpGrain[i].rotmat[j][k];
            }
        }
    }

    xfree(cpGrain);

    //restore bp->num_grains
    hsize_t bpDims[1] = { 1 };
    hid_t bpSpace = H5Screate_simple(1, bpDims, NULL);
    hid_t numGrain =
        H5Dopen(file, get_subblock_cp_path(NULL, iproc, "numGrain"),
                H5P_DEFAULT);
    H5Dread(numGrain, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            &bp->num_grains);

    //release resource
    H5Sclose(bpSpace);
    H5Dclose(numGrain);

    H5Sclose(grainSpace);
    H5Tclose(memGrain);
    H5Tclose(memRotmat);

}

static void
mainAddSubblockMapping(
    uint64_t restart,
    hid_t file)
{
    /*
     * INPUT: restart: restart timestep
     * file: file identifier of main checkpoint
     * OUTPUT: gtasks will be restored
     * NOTE: This function is used for main restart
     */


    taskData *wdata;            // Pointer to vlen structures //
    hsize_t wDims[1] = { nproc };
    int *ptr;

    xmalloc(wdata, taskData, nproc);

    hid_t taskSpace = H5Screate_simple(1, wDims, NULL);
    hid_t taskAssignedSB = H5Tvlen_create(H5T_NATIVE_INT);

    //create datatype for gtasks
    hid_t taskType = H5Tcreate(H5T_COMPOUND, sizeof(taskData));
    H5Tinsert(taskType, "rank", HOFFSET(taskData, rank), H5T_NATIVE_INT);
    H5Tinsert(taskType, "assigned_sb", HOFFSET(taskData, assigned_sb),
              taskAssignedSB);

    //read gtasks data from main checkpoint file
    hid_t wSet =
        H5Dopen(file, get_subblock_cp_path(NULL, iproc, "taskSubMap"),
                H5P_DEFAULT);

    H5Dread(wSet, taskType, H5S_ALL, H5S_ALL, H5P_DEFAULT, wdata);

    for (int i = 0; i < nproc; i++)
    {
        gtasks.task[i].rank = wdata[i].rank;
        ptr = (int *) wdata[i].assigned_sb.p;
        //for (int j=0;j<wdata[i].assigned_sb.len;j++)
        {
            gtasks.task[i].assigned_sb = ptr[0];
        }
    }


    // Now update the global subblock map.
    lli_t *llp = gsubblock_list->head;
    while (llp != NULL)
    {
        MSB_struct *msb = (MSB_struct *) llp->data;
        // Figure out who my neighboring subblocks are
        int blockid = msb->subblockid;
        int neighbors[NUM_NEIGHBORS];

        determine_3dneighbors(blockid, neighbors);

        // Go through my neighbors, set the processor and subblock ids in my neighbor map
        for (int i = 0; i < NUM_NEIGHBORS; i++)
        {
            if (neighbors[i] >= 0)
            {
                MSB_struct *neighbor = &gmsp[neighbors[i]];
                msb->neighbors[i][0] = neighbor->procid;
                dwrite(DEBUG_MAIN_CTRL,
                       "sb %3d Setting neighbor[%d][0,1] = %3d, %3d\n",
                       msb->subblockid, i,
                       msb->neighbors[i][0], msb->neighbors[i][1]);
            }
        }
        llp = llp->next;
    }

    // Rebuild global (scatter) neighbor maps
    int nbufsz = ll_count(gsubblock_list);
    assert(nbufsz > 0);
    xrealloc(gtasks.nbrbuf, nbr_info_t, nbufsz);

    gtasks.nbrdisp[0] = 0;
    gtasks.nbrcnts[0] = 0;

    for (int p = 1, displacement = 0; p < nproc; p++)
    {
        task_t *w = &gtasks.task[p];

        gtasks.nbrdisp[p] = displacement;
        gtasks.nbrcnts[p] = 1;

        // Fill in neighbor info
        nbr_info_t *ni = &gtasks.nbrbuf[displacement];
        ni[0].id = w->assigned_sb;
        memcpy(&(ni[0].nbors), gmsp[w->assigned_sb].neighbors,
               sizeof(neighbors_t));

        // Advance to the next slot
        displacement += gtasks.nbrcnts[p];
    }

    //release resource
    H5Sclose(taskSpace);
    H5Tclose(taskAssignedSB);
    H5Dclose(wSet);
    H5Tclose(taskType);

    for (int i = 0; i < nproc; i++)
    {
        free(wdata[i].assigned_sb.p);
    }
    xfree(wdata);

}


void
mainRestart(
    uint64_t restart)
{
    /*
     * INPUT: restart: restart file ID at specified timestep if DEBUG_CHECKPOINT==1
     * restart file ID = cp_main.h5 if DEBUG_CHECKPOINT==0
     * OUTPUT: gmsp, gsubblock_list, gtasks will be restored
     * USE: mainAddSubblock,mainActivateSubblock,addGrainCache, mainAddSubblockMapping
     */
    hid_t file_id;

    file_id = H5Fopen(get_rank_cp_file(iproc), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id <= 0)
        printf("\n\n Main task COULD NOT open restart file: %s\n",
               get_rank_cp_file(iproc));
    assert(file_id > 0);


    //restore current timestep
    hid_t curTimeStep;
    curTimeStep =
        H5Dopen(file_id, get_subblock_cp_path(NULL, iproc, "curTimeStep"),
                H5P_DEFAULT);

    H5Dread(curTimeStep, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            &bp->timestep);
    H5Dclose(curTimeStep);

    //restore gmsp
    mainAddSubblock(restart, file_id);
    //add gmsp to gsubblock_list
    mainActivateSubblock(gsubblock_list);
    //restore grain cache
    addGrainCache(file_id);
    //restore gtasks
    mainAddSubblockMapping(restart, file_id);

    H5Fclose(file_id);

}

void
taskRestart(
    uint64_t restart)
{
    /*
     * INPUT: restart: restart file ID at specified timestep if DEBUG_CHECKPOINT==1
     * restart file ID = cp_task_###.h5, ### is task ID if DEBUG_CHECKPOINT==0
     * OUTPUT: gProfile.ts_idx, and subblock data will be restored
     * NOTE: This function is used for task restart
     * USE: H5Literate,restoreTaskSubList
     */

    hid_t file_id;

    file_id = H5Fopen(get_rank_cp_file(iproc), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id <= 0)
        printf("\n\n Task %d COULD NOT open restart file: %s\n", iproc,
               get_rank_cp_file(iproc));
    assert(file_id > 0);

    //restore timestep
    hid_t curTimeStep;
    curTimeStep =
        H5Dopen(file_id, get_subblock_cp_path(NULL, iproc, "curTimeStep"),
                H5P_DEFAULT);
    H5Dread(curTimeStep, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            &bp->timestep);
    H5Dclose(curTimeStep);

    hid_t group = H5Gopen(file_id, "SUBBLOCK", H5P_DEFAULT);

    //restore subblock
    H5Literate(group, H5_INDEX_CRT_ORDER, H5_ITER_INC, NULL,
               restoreTaskSubList, NULL);

    addGrainCache(file_id);

    printf("task %d done restarting at timestep %lu\n", iproc, bp->timestep);

    H5Fclose(file_id);
    H5Gclose(group);
}




void
writeTaskCheckpoint(
    uint64_t restart)
{
    /*
     * INPUT: restart: the checkpoint timestep (current timestep)
     * OUTPUT: current timestep, gProfile.ts_idx, and all data in subblock
     * will be saved in the task checkpoint file
     * USE: get_rank_cp_path, add_subblock
     */

    /* To reduce the IO operations, I collect all of static variables into one large
     * variable. Then we only need to use one IO operation to save all these
     * variables.
     */
    typedef struct cp_taskSubblock
    {
        COMMON_SUBBLOCK_DEFS;
        int last_active_ts;
        size_t mold_size;
    } taskSubblock;

    taskSubblock cpSubblock;

    hid_t file_id;
    //create file identifier for the checkpoint file
    file_id =
        H5Fcreate(get_rank_cp_path(iproc, restart), H5F_ACC_TRUNC,
                  H5P_DEFAULT, H5P_DEFAULT);

    //set the dimension size for the subblocks
    hsize_t start[3] = { 0, 0, 0 };
    hsize_t sbdims[3];
    sbdims[0] = bp->gsdimz + 2;
    sbdims[1] = bp->gsdimy + 2;
    sbdims[2] = bp->gsdimx + 2;

    //create the chunk size
    hsize_t chunkdims[3];
    chunkdims[0] = bp->gsdimz + 2;
    chunkdims[1] = bp->gsdimy + 2;
    chunkdims[2] = bp->gsdimx + 2;

    //create the dataspace
    hid_t dataspace = H5Screate_simple(3, sbdims, NULL);
    hid_t link_pl = H5Pcreate(H5P_LINK_CREATE);
    H5Pset_create_intermediate_group(link_pl, 1 /* Positive means YES */ );

    hid_t dataspace_pl = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(dataspace_pl, 3, chunkdims);
    H5Pset_fill_time(dataspace_pl, H5D_FILL_TIME_NEVER);
    H5Pset_shuffle(dataspace_pl);
    H5Pset_deflate(dataspace_pl, 9);

    hid_t memspace = H5Screate_simple(3, sbdims, NULL);
    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, start, NULL, sbdims, NULL);


    hid_t dendriteType = H5Tcreate(H5T_COMPOUND, sizeof(decentered_t));
    H5Tinsert(dendriteType, "x", HOFFSET(decentered_t, x), H5T_NATIVE_DOUBLE);
    H5Tinsert(dendriteType, "y", HOFFSET(decentered_t, y), H5T_NATIVE_DOUBLE);
    H5Tinsert(dendriteType, "z", HOFFSET(decentered_t, z), H5T_NATIVE_DOUBLE);


    hid_t memtype, memCoord, memActiveState, memNeiBor, space;
    hsize_t dims[1] = { 1 };
    hsize_t neiBorDims[2] = { NUM_NEIGHBORS, NBOR_INFO };
    hsize_t uint64Size[1] = { 1 };

    hid_t u64Space = H5Screate_simple(1, uint64Size, NULL);


    memCoord = H5Tcreate(H5T_COMPOUND, sizeof(loc_t));
    H5Tinsert(memCoord, "x", HOFFSET(loc_t, x), H5T_NATIVE_UINT32);
    H5Tinsert(memCoord, "y", HOFFSET(loc_t, y), H5T_NATIVE_UINT32);
    H5Tinsert(memCoord, "z", HOFFSET(loc_t, z), H5T_NATIVE_UINT32);

    activestate_t val;

    memActiveState = H5Tenum_create(H5T_NATIVE_INT);
    val = (activestate_t) INACTIVE;
    H5Tenum_insert(memActiveState, "inactive", &val);
    val = (activestate_t) ACTIVE;
    H5Tenum_insert(memActiveState, "active", &val);

    memNeiBor = H5Tarray_create(H5T_NATIVE_INT, 2, neiBorDims);

    //create datatype for taskSubblock
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(taskSubblock));
    H5Tinsert(memtype, "last_active",
              HOFFSET(taskSubblock, last_active_ts), H5T_NATIVE_INT);
    H5Tinsert(memtype, "coord", HOFFSET(taskSubblock, coords), memCoord);

    H5Tinsert(memtype, "subblockid",
              HOFFSET(taskSubblock, subblockid), H5T_NATIVE_INT);
    H5Tinsert(memtype, "procid",
              HOFFSET(taskSubblock, procid), H5T_NATIVE_INT);

    H5Tinsert(memtype, "neighbor", HOFFSET(taskSubblock, neighbors),
              memNeiBor);
    H5Tinsert(memtype, "nnuc", HOFFSET(taskSubblock, nnuc), H5T_NATIVE_INT);
    H5Tinsert(memtype, "mold_size",
              HOFFSET(taskSubblock, mold_size), H5T_NATIVE_UINT64);

    hid_t nucleationType = H5Tcreate(H5T_COMPOUND, sizeof(nucleation_t));
    H5Tinsert(nucleationType, "x", HOFFSET(nucleation_t, x),
              H5T_NATIVE_DOUBLE);
    H5Tinsert(nucleationType, "y", HOFFSET(nucleation_t, y),
              H5T_NATIVE_DOUBLE);
    H5Tinsert(nucleationType, "z", HOFFSET(nucleation_t, z),
              H5T_NATIVE_DOUBLE);
    H5Tinsert(nucleationType, "cx", HOFFSET(nucleation_t, cx),
              H5T_NATIVE_UINT64);
    H5Tinsert(nucleationType, "cy", HOFFSET(nucleation_t, cy),
              H5T_NATIVE_UINT64);
    H5Tinsert(nucleationType, "cz", HOFFSET(nucleation_t, cz),
              H5T_NATIVE_UINT64);
    H5Tinsert(nucleationType, "gx", HOFFSET(nucleation_t, gx),
              H5T_NATIVE_DOUBLE);
    H5Tinsert(nucleationType, "gy", HOFFSET(nucleation_t, gy),
              H5T_NATIVE_DOUBLE);
    H5Tinsert(nucleationType, "gz", HOFFSET(nucleation_t, gz),
              H5T_NATIVE_DOUBLE);
    H5Tinsert(nucleationType, "threshold", HOFFSET(nucleation_t, threshold),
              H5T_NATIVE_DOUBLE);


    /*
     * Create dataspace.  Setting maximum size to NULL sets the maximum
     * size to be the current size.
     */

    space = H5Screate_simple(1, dims, NULL);

    //write current timestep into checkpoint file
    hid_t curTimeStep;

    curTimeStep =
        H5Dcreate(file_id, get_subblock_cp_path(NULL, iproc, "curTimeStep"),
                  H5T_NATIVE_UINT32, u64Space, H5P_DEFAULT, H5P_DEFAULT,
                  H5P_DEFAULT);
    H5Dwrite(curTimeStep, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT,
             &bp->timestep);

    H5Dclose(curTimeStep);


    //checkpoint ts_idx in gProfile.ts_idx
    hid_t g_ts_idx;
    switch (bp->temp_type)
    {
        case INTERNAL:         //dummy data
            g_ts_idx =
                H5Dcreate(file_id,
                          get_subblock_cp_path(NULL, iproc, "g_ts_idx"),
                          H5T_NATIVE_UINT32, u64Space, H5P_DEFAULT,
                          H5P_DEFAULT, H5P_DEFAULT);

            H5Dwrite(g_ts_idx, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL,
                     H5P_DEFAULT, &bp->temp_type);
            H5Dclose(g_ts_idx);


            break;
    }


    //write all data of subblock into the checkpoint file by
    //creation order
    hid_t gcpl = H5Pcreate(H5P_GROUP_CREATE);
    H5Pset_link_creation_order(gcpl, H5P_CRT_ORDER_TRACKED |
                               H5P_CRT_ORDER_INDEXED);

    hid_t subblock =
        H5Gcreate(file_id, "SUBBLOCK", H5P_DEFAULT, gcpl, H5P_DEFAULT);

    {
        SB_struct *sb = lsp;

        hid_t subgroup =
            H5Gcreate(subblock, add_subblock(sb), H5P_DEFAULT, gcpl,
                      H5P_DEFAULT);

        //Collect all values of small variable into one compound variable
        cpSubblock.last_active_ts = sb->last_active_ts;
        cpSubblock.subblockid = sb->subblockid;
        cpSubblock.procid = sb->procid;
        cpSubblock.coords.x = sb->coords.x;
        cpSubblock.coords.y = sb->coords.y;
        cpSubblock.coords.z = sb->coords.z;
        cpSubblock.nnuc = sb->nnuc;
        cpSubblock.mold_size = sb->mold_size;

        for (int j = 0; j < NUM_NEIGHBORS; j++)
        {
            for (int k = 0; k < NBOR_INFO; k++)
            {
                cpSubblock.neighbors[j][k] = sb->neighbors[j][k];
            }
        }


        /*
         * Create the dataset and write the compound data to it. Use only one IO
         * operation.
         */
        hid_t dset = H5Dcreate(subgroup, "CoordTemp",
                               memtype, space, link_pl, H5P_DEFAULT,
                               H5P_DEFAULT);
        H5Dwrite(dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &cpSubblock);
        H5Dclose(dset);


#define WRITE_FIELD(name, type) \
		do { \
			/* Only write ones that exist */ \
			if ( sb->name ) \
			{ \
				hid_t dset = H5Dcreate(subgroup, #name, type, dataspace, link_pl, dataspace_pl, H5P_DEFAULT); \
				H5Dwrite(dset, type, memspace, H5S_ALL, H5P_DEFAULT, sb->name); \
				H5Dclose(dset); \
			} \
		} while(0)

        WRITE_FIELD(temperature, H5T_NATIVE_DOUBLE);
        WRITE_FIELD(gr, H5T_NATIVE_INT);
        WRITE_FIELD(fs, H5T_NATIVE_DOUBLE);
        WRITE_FIELD(ce, H5T_NATIVE_DOUBLE);
        WRITE_FIELD(oce, H5T_NATIVE_DOUBLE);
        WRITE_FIELD(cl, H5T_NATIVE_DOUBLE);
        WRITE_FIELD(diff_id, H5T_NATIVE_CHAR);
        WRITE_FIELD(mold, H5T_NATIVE_CHAR);
        WRITE_FIELD(d, H5T_NATIVE_DOUBLE);
        WRITE_FIELD(nuc_threshold, H5T_NATIVE_FLOAT);
        WRITE_FIELD(dc, dendriteType);
//ADD_CURV_LY
        WRITE_FIELD(curv, H5T_NATIVE_DOUBLE);

#undef WRITE_FIELD

        if (sb->nnuc > 0)
        {
            hsize_t len = sb->nnuc;
            hid_t lenhdl = H5Screate_simple(1, &len, NULL);
            if (lenhdl < 0)
                error("H5Screate_simple failed!\n");

            hid_t nuc_pt_hdl = H5Dcreate(subgroup, "nuc_pt", nucleationType,
                                         lenhdl, link_pl, H5P_DEFAULT,
                                         H5P_DEFAULT);
            if (nuc_pt_hdl < 0)
                error("H5Dcreate failed!\n");
            H5Dwrite(nuc_pt_hdl, nucleationType, H5S_ALL, H5S_ALL,
                     H5P_DEFAULT, sb->nuc_pt);
            H5Dclose(nuc_pt_hdl);
            H5Sclose(lenhdl);
        }

        /*sb->temp_info->deltas, sb->temp_info->lastTS */

        switch (bp->temp_type)
        {
            case INTERNAL:
                break;
        }

        H5Gclose(subgroup);

    }
    storeGlobalGrainInfo(file_id);
    H5Gclose(subblock);

    H5Sclose(u64Space);

    H5Fclose(file_id);

    H5Pclose(dataspace_pl);
    H5Pclose(link_pl);
    H5Sclose(dataspace);
    H5Sclose(memspace);
    H5Sclose(space);
    H5Tclose(memActiveState);
    H5Tclose(memCoord);
    H5Tclose(memNeiBor);

    H5Tclose(memtype);
    H5Tclose(nucleationType);
}




void
writeMainCheckpoint(
    )
{
    /*
     * INPUT: none
     * OUTPUT: bp->timestep,grain_cache, bp->num_grains, gtasks, gmsp will be saved in the main checkpoint
     * USE: get_rank_cp_path, get_subblock_cp_path
     */

    hid_t file, memtype, memCoord, memActiveState, memNeiBor, space, dset;
    /* Handles */
    int dimSize = bp->gnsbx * bp->gnsby * bp->gnsbz;
    hsize_t mDims[1] = { dimSize };
    hsize_t neiBorDims[2] = { NUM_NEIGHBORS, NBOR_INFO };

    //create datatype for loc_t
    memCoord = H5Tcreate(H5T_COMPOUND, sizeof(loc_t));
    H5Tinsert(memCoord, "x", HOFFSET(loc_t, x), H5T_NATIVE_UINT32);
    H5Tinsert(memCoord, "y", HOFFSET(loc_t, y), H5T_NATIVE_UINT32);
    H5Tinsert(memCoord, "z", HOFFSET(loc_t, z), H5T_NATIVE_UINT32);

    //create datatype for activestate_t
    activestate_t val;

    memActiveState = H5Tenum_create(H5T_NATIVE_INT);
    val = (activestate_t) INACTIVE;
    H5Tenum_insert(memActiveState, "inactive", &val);
    val = (activestate_t) ACTIVE;
    H5Tenum_insert(memActiveState, "active", &val);

    memNeiBor = H5Tarray_create(H5T_NATIVE_INT, 2, neiBorDims);

    hid_t dataspace_pl = H5Pcreate(H5P_DATASET_CREATE);

    H5Pset_chunk(dataspace_pl, 1, mDims);
    H5Pset_fill_time(dataspace_pl, H5D_FILL_TIME_NEVER);
    H5Pset_shuffle(dataspace_pl);
    H5Pset_deflate(dataspace_pl, 9);

    /*
     * Create a new file using the default properties.
     */
    file =
        H5Fcreate(get_rank_cp_path(iproc, bp->timestep), H5F_ACC_TRUNC,
                  H5P_DEFAULT, H5P_DEFAULT);

    hsize_t bpDims[1] = { 1 };
    hid_t bpSpace = H5Screate_simple(1, bpDims, NULL);

    //save current timestep
    hid_t curTimeStep =
        H5Dcreate(file, get_subblock_cp_path(NULL, iproc, "curTimeStep"),
                  H5T_NATIVE_INT, bpSpace, H5P_DEFAULT, H5P_DEFAULT,
                  H5P_DEFAULT);
    H5Dwrite(curTimeStep, H5T_NATIVE_UINT64, H5S_ALL, H5S_ALL, H5P_DEFAULT,
             &bp->timestep);
    H5Dclose(curTimeStep);
    H5Sclose(bpSpace);


    //create dataset for grain cache
    storeGlobalGrainInfo(file);


    //save gmsp
    /*
     * Create the compound datatype for memory.
     */
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(MSB_struct));
    H5Tinsert(memtype, "subblockid",
              HOFFSET(MSB_struct, subblockid), H5T_NATIVE_INT);
    H5Tinsert(memtype, "procid", HOFFSET(MSB_struct, procid), H5T_NATIVE_INT);
    H5Tinsert(memtype, "activate_ts",
              HOFFSET(MSB_struct, activate_ts), H5T_NATIVE_UINT64);

    H5Tinsert(memtype, "coord", HOFFSET(MSB_struct, coords), memCoord);

    H5Tinsert(memtype, "neighbor", HOFFSET(MSB_struct, neighbors), memNeiBor);
    H5Tinsert(memtype, "nnuc", HOFFSET(MSB_struct, nnuc), H5T_NATIVE_INT);

    hid_t nucleationType = H5Tcreate(H5T_COMPOUND, sizeof(nucleation_t));
    H5Tinsert(nucleationType, "x", HOFFSET(nucleation_t, x),
              H5T_NATIVE_DOUBLE);
    H5Tinsert(nucleationType, "y", HOFFSET(nucleation_t, y),
              H5T_NATIVE_DOUBLE);
    H5Tinsert(nucleationType, "z", HOFFSET(nucleation_t, z),
              H5T_NATIVE_DOUBLE);
    H5Tinsert(nucleationType, "cx", HOFFSET(nucleation_t, cx),
              H5T_NATIVE_UINT64);
    H5Tinsert(nucleationType, "cy", HOFFSET(nucleation_t, cy),
              H5T_NATIVE_UINT64);
    H5Tinsert(nucleationType, "cz", HOFFSET(nucleation_t, cz),
              H5T_NATIVE_UINT64);
    H5Tinsert(nucleationType, "gx", HOFFSET(nucleation_t, gx),
              H5T_NATIVE_DOUBLE);
    H5Tinsert(nucleationType, "gy", HOFFSET(nucleation_t, gy),
              H5T_NATIVE_DOUBLE);
    H5Tinsert(nucleationType, "gz", HOFFSET(nucleation_t, gz),
              H5T_NATIVE_DOUBLE);
    H5Tinsert(nucleationType, "threshold", HOFFSET(nucleation_t, threshold),
              H5T_NATIVE_DOUBLE);

    /*
     * Create dataspace.  Setting maximum size to NULL sets the maximum
     * size to be the current size.
     */
    space = H5Screate_simple(1, mDims, NULL);

    /*
     * Create the dataset and write the compound data to it.
     */
    dset =
        H5Dcreate(file, get_subblock_cp_path(NULL, iproc, "mainSubblock"),
                  memtype, space, H5P_DEFAULT, dataspace_pl, H5P_DEFAULT);
    H5Dwrite(dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, gmsp);

    hid_t link_pl = H5Pcreate(H5P_LINK_CREATE);
    H5Pset_create_intermediate_group(link_pl, 1 /* Positive means YES */ );

    //save nuc_pt
    for (int i = 0; i < dimSize; i++)
    {
        if (gmsp[i].nnuc > 0)
        {
            char path[PATH_MAX] = { 0 };
            snprintf(path, PATH_MAX - 1, "/nuc_pts/sb_%d", i);
            hsize_t len = gmsp[i].nnuc;
            hid_t lenhdl = H5Screate_simple(1, &len, NULL);
            hid_t nuc_pt_hdl = H5Dcreate(file, path, nucleationType,
                                         lenhdl, link_pl, dataspace_pl,
                                         H5P_DEFAULT);
            H5Dwrite(nuc_pt_hdl, nucleationType, H5S_ALL, H5S_ALL,
                     H5P_DEFAULT, gmsp[i].nuc_pt);
            H5Dclose(nuc_pt_hdl);
            H5Sclose(lenhdl);
        }
    }
    H5Pclose(link_pl);


    //save gtasks
    taskData wdata[nproc];      // Pointer to vlen structures //
    hsize_t wDims[1] = { nproc };

    int *ptr;

    for (int i = 0; i < nproc; i++)
    {
        wdata[i].rank = gtasks.task[i].rank;

        wdata[i].assigned_sb.len = 1;
        ptr = (int *) malloc(wdata[i].assigned_sb.len * sizeof(int));
        //for (int j=0;j<wdata[i].assigned_sb.len;j++)
        {
            ptr[0] = gtasks.task[i].assigned_sb;
        }
        wdata[i].assigned_sb.p = (void *) ptr;

    }

    hid_t taskSpace = H5Screate_simple(1, wDims, NULL);
    hid_t taskAssignedSB = H5Tvlen_create(H5T_NATIVE_INT);


    hid_t taskType = H5Tcreate(H5T_COMPOUND, sizeof(taskData));
    H5Tinsert(taskType, "rank", HOFFSET(taskData, rank), H5T_NATIVE_INT);
    H5Tinsert(taskType, "assigned_sb", HOFFSET(taskData, assigned_sb),
              taskAssignedSB);

    hid_t wSet =
        H5Dcreate(file, get_subblock_cp_path(NULL, iproc, "taskSubMap"),
                  taskType, taskSpace, H5P_DEFAULT, H5P_DEFAULT,
                  H5P_DEFAULT);


    H5Dwrite(wSet, taskType, H5S_ALL, H5S_ALL, H5P_DEFAULT, wdata);


    H5Sclose(taskSpace);
    H5Tclose(taskAssignedSB);
    H5Dclose(wSet);
    H5Tclose(taskType);

    for (int i = 0; i < nproc; i++)
    {
        free(wdata[i].assigned_sb.p);
    }

    /*
     * Close and release resources.
     */
    H5Dclose(dset);
    H5Sclose(space);
    H5Tclose(memtype);
    H5Tclose(memActiveState);
    H5Tclose(memNeiBor);
    H5Tclose(memCoord);
    H5Tclose(nucleationType);
    H5Fclose(file);
    H5Pclose(dataspace_pl);

}
