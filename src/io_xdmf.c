/***************************************************************/
/* Copyright (c) 2023, Lang Yuan, Univeristy of South Carolina */
/* All rights reserved.                                        */
/* This file is part of muMatScale.                            */
/* See the top-level LICENSE file for details.                 */
/***************************************************************/
#include <stdio.h>
#include <inttypes.h>
#include <errno.h>
#include <endian.h>
#include <limits.h>
#include <time.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <mpi.h>
#include <hdf5.h>
#include "globals.h"
#include "functions.h"
#include "debug.h"
#include "ll.h"
#include "grain.h"
#include "xmalloc.h"
#include "profiler.h"

#ifndef PATH_MAX
#define PATH_MAX 512
#endif

extern SB_struct *lsp;

static const char *
get_msubblock_path(
    int itask,
    const char *var)
{
    static char path[PATH_MAX] = { 0 };
    int coords[3];
    MPI_Cart_coords(mpi_comm_new, itask, nproc, coords);
    snprintf(path, PATH_MAX - 1, "/SubBlock_%dx%dx%d/%s",
             coords[2], coords[1], coords[0], var);
    return path;
}

static const char *
get_subblock_path(
    SB_struct * sb,
    const char *var)
{
    static char path[PATH_MAX] = { 0 };
    int coords[3];
    MPI_Cart_coords(mpi_comm_new, iproc, nproc, coords);
    snprintf(path, PATH_MAX - 1, "/SubBlock_%dx%dx%d/%s",
             coords[2], coords[1], coords[0], var);
    return path;
}

static const char *
get_rank_dir(
    int rank)
{
    static char path[PATH_MAX] = { 0 };
    snprintf(path, PATH_MAX - 1, "%s_r%d", bp->basefilename, rank);
    return path;
}

static const char *
get_rank_file(
    int rank,
    uint64_t ts)
{
    static char path[PATH_MAX] = { 0 };
    snprintf(path, PATH_MAX - 1, "%s/%05lu.h5", get_rank_dir(rank), ts);
    return path;
}

/**
 * This method is similar to get_rank_file, but gives a local directory name.
 * If the bp->basefilename includes a directory separator, we don't want to have
 * that directory referenced again in our XDMF file.
 */
static const char *
get_local_rank_file(
    int rank,
    uint64_t ts)
{
    static char path[PATH_MAX] = { 0 };
    char *fptr = strrchr(bp->basefilename, '/');
    if (!fptr)
        fptr = bp->basefilename;
    else
        fptr++;                 // Advance past the '/'
    snprintf(path, PATH_MAX - 1, "%s_r%d/%05lu.h5", fptr, rank, ts);
    return path;
}

static const char *
get_global_file(
    void)
{
    static char path[PATH_MAX] = { 0 };
    snprintf(path, PATH_MAX - 1, "%s_gbl.h5", bp->basefilename);
    return path;
}


/**
 * This method is similar to get_global_file, but gives a local directory name.
 * If the bp->basefilename includes a directory separator, we don't want to have
 * that directory referenced again in our XDMF file.
 */
static const char *
get_local_global_file(
    void)
{
    static char path[PATH_MAX] = { 0 };
    char *fptr = strrchr(bp->basefilename, '/');
    if (!fptr)
        fptr = bp->basefilename;
    else
        fptr++;                 // Advance past the '/'
    snprintf(path, PATH_MAX - 1, "%s_gbl.h5", fptr);
    return path;
}


#ifdef WRITE_TIMESTEPED_XMF
static const char *
get_file_path(
    uint64_t timestamp,
    const char *ext)
{
    static char path[PATH_MAX] = { 0 };
    snprintf(path, PATH_MAX - 1, "%s_%05lu.%s",
             bp->basefilename, timestamp, ext);
    return path;
}
#endif


static void
writeXDMFGrids(
    FILE * fp,
    uint64_t timestamp,
    const char *indent)
{
    fprintf(fp, "%s<Grid Name=\"Full Volume at time %g\" "
            "GridType=\"Collection\" CollectionType=\"Spatial\">\n",
            indent, timestamp * bp->ts_delt);
    fprintf(fp, "%s\t<Time Value=\"%g\" />\n", indent,
            timestamp * bp->ts_delt);
    // Write each block
    uint32_t num_sb = bp->gnsbx * bp->gnsby * bp->gnsbz;
    for (uint32_t i = 0; i < num_sb; i++)
    {
        int coords[3];
        MPI_Cart_coords(mpi_comm_new, i, num_sb, coords);

        fprintf(fp, "%s\t<Grid Name=\"SubBlock %dx%dx%d (%d)\" "
                "GridType=\"Uniform\">\n",
                indent,
                coords[2], coords[1], coords[0], i);
        fprintf(fp, "%s\t\t<Topology TopologyType=\"3DCoRectMesh\" "
                "Dimensions=\"%u %u %u\" />\n",
                // We do a +1 because we're going to be Cell data
                indent, bp->gsdimz + 1, bp->gsdimy + 1, bp->gsdimx + 1);
        fprintf(fp, "%s\t\t<Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n",
                indent);
        fprintf(fp, "%s\t\t\t<DataItem Format=\"XML\" Dimensions=\"3\">\n",
                indent);
        // Origin
        double rx, ry, rz, ux, uy, uz;
        getSBRealBounds(coords[2], coords[1], coords[0], &rx, &ry,
                        &rz, &ux, &uy, &uz);
        fprintf(fp, "%s\t\t\t\t%e %e %e\n", indent, rz, ry, rx);
        fprintf(fp, "%s\t\t\t</DataItem>\n", indent);
        fprintf(fp, "%s\t\t\t<DataItem Format=\"XML\" Dimensions=\"3\">\n",
                indent);
        // Size
        fprintf(fp, "%s\t\t\t\t%e %e %e\n",
                indent, bp->cellSize, bp->cellSize, bp->cellSize);
        fprintf(fp, "%s\t\t\t</DataItem>\n", indent);
        fprintf(fp, "%s\t\t</Geometry>\n", indent);

        // Now, for the actual *DATA*

                /***** TEMPERATURE *****/
        fprintf(fp,
                "%s\t\t<Attribute Name=\"Temperature\" Center=\"Cell\">\n",
                indent);
        fprintf(fp,
                "%s\t\t\t<DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"%lu\" Dimensions=\"%u %u %u\">\n",
                indent, sizeof(double), bp->gsdimz, bp->gsdimy, bp->gsdimx);
        fprintf(fp, "%s\t\t\t\tSERIAL:%s:%s\n", indent,
                get_local_rank_file(i, timestamp),
                get_msubblock_path(i, "Temperature"));
        fprintf(fp, "%s\t\t\t</DataItem>\n", indent);
        fprintf(fp, "%s\t\t</Attribute>\n", indent);



                /***** Solute Composition *****/
        if (bp->calc_type == DIFFUSION)
        {
            fprintf(fp,
                    "%s\t\t<Attribute Name=\"Solute Composition\" Center=\"Cell\">\n",
                    indent);
            fprintf(fp,
                    "%s\t\t\t<DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"%lu\" Dimensions=\"%u %u %u\">\n",
                    indent, sizeof(double), bp->gsdimz, bp->gsdimy,
                    bp->gsdimx);
            fprintf(fp, "%s\t\t\t\tSERIAL:%s:%s\n", indent,
                    get_local_rank_file(i, timestamp),
                    get_msubblock_path(i, "CE"));
            fprintf(fp, "%s\t\t\t</DataItem>\n", indent);
            fprintf(fp, "%s\t\t</Attribute>\n", indent);
        }

                /***** GRAIN *****/
        fprintf(fp, "%s\t\t<Attribute Name=\"Grain\" Center=\"Cell\">\n",
                indent);
        if (bp->randomize_grains)
        {
            fprintf(fp,
                    "%s\t\t\t<DataItem ItemType=\"Function\" NumberType=\"UInt\" Dimensions=\"%u %u %u\" Function=\"$0[$1]\">\n",
                    indent, bp->gsdimz, bp->gsdimy, bp->gsdimx);
            fprintf(fp,
                    "%s\t\t\t\t<DataItem Reference=\"XML\" Dimensions=\"%u %u %u\">\n",
                    indent, bp->gsdimz, bp->gsdimy, bp->gsdimx);
            fprintf(fp,
                    "%s\t\t\t\t\t/Xdmf/Domain/DataItem[@Name=\"Grain Renumber\"]\n",
                    indent);
            fprintf(fp, "%s\t\t\t\t</DataItem>\n", indent);
        }
        fprintf(fp,
                "%s\t\t\t\t<DataItem Format=\"HDF\" NumberType=\"UInt\" Dimensions=\"%u %u %u\">\n",
                indent, bp->gsdimz, bp->gsdimy, bp->gsdimx);
        fprintf(fp, "%s\t\t\t\t\tSERIAL:%s:%s\n", indent,
                get_local_rank_file(i, timestamp),
                get_msubblock_path(i, "Grain"));
        fprintf(fp, "%s\t\t\t\t</DataItem>\n", indent);
        if (bp->randomize_grains)
            fprintf(fp, "%s\t\t\t</DataItem>\n", indent);
        fprintf(fp, "%s\t\t</Attribute>\n", indent);

                /***** Grain Angle ****/
        fprintf(fp,
                "%s\t\t<Attribute Name=\"Grain Angle\" Center=\"Cell\" AttributeType=\"Vector\">\n",
                indent);
        fprintf(fp,
                "%s\t\t\t<DataItem ItemType=\"Function\" NumberType=\"Float\" Dimensions=\"%u %u %u 3\" Function=\"$1[$0], $2[$0], $3[$0]\">\n",
                indent, bp->gsdimz, bp->gsdimy, bp->gsdimx);
        fprintf(fp,
                "%s\t\t\t\t<DataItem Format=\"HDF\" NumberType=\"UInt\" Dimensions=\"%u %u %u\">\n",
                indent, bp->gsdimz, bp->gsdimy, bp->gsdimx);
        fprintf(fp, "%s\t\t\t\t\tSERIAL:%s:%s\n", indent,
                get_local_rank_file(i, timestamp),
                get_msubblock_path(i, "Grain"));
        fprintf(fp, "%s\t\t\t\t</DataItem>\n", indent);
        fprintf(fp,
                "%s\t\t\t\t<DataItem Reference=\"XML\" Dimensions=\"%u %u %u\">\n",
                indent, bp->gsdimz, bp->gsdimy, bp->gsdimx);
        fprintf(fp,
                "%s\t\t\t\t\t/Xdmf/Domain/DataItem[@Name=\"Grain Angle X\"]\n",
                indent);
        fprintf(fp, "%s\t\t\t\t</DataItem>\n", indent);
        fprintf(fp,
                "%s\t\t\t\t<DataItem Reference=\"XML\" Dimensions=\"%u %u %u\">\n",
                indent, bp->gsdimz, bp->gsdimy, bp->gsdimx);
        fprintf(fp,
                "%s\t\t\t\t\t/Xdmf/Domain/DataItem[@Name=\"Grain Angle Y\"]\n",
                indent);
        fprintf(fp, "%s\t\t\t\t</DataItem>\n", indent);
        fprintf(fp,
                "%s\t\t\t\t<DataItem Reference=\"XML\" Dimensions=\"%u %u %u\">\n",
                indent, bp->gsdimz, bp->gsdimy, bp->gsdimx);
        fprintf(fp,
                "%s\t\t\t\t\t/Xdmf/Domain/DataItem[@Name=\"Grain Angle Z\"]\n",
                indent);
        fprintf(fp, "%s\t\t\t\t</DataItem>\n", indent);
        fprintf(fp, "%s\t\t\t</DataItem>\n", indent);
        fprintf(fp, "%s\t\t</Attribute>\n", indent);


                /***** Fraction Solid ****/
        fprintf(fp,
                "%s\t\t<Attribute Name=\"Fraction Solid\" Center=\"Cell\">\n",
                indent);
        fprintf(fp,
                "%s\t\t\t<DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"%lu\" Dimensions=\"%u %u %u\">\n",
                indent, sizeof(double), bp->gsdimz, bp->gsdimy, bp->gsdimx);
        fprintf(fp, "%s\t\t\t\tSERIAL:%s:%s\n", indent,
                get_local_rank_file(i, timestamp),
                get_msubblock_path(i, "FracSolid"));
        fprintf(fp, "%s\t\t\t</DataItem>\n", indent);
        fprintf(fp, "%s\t\t</Attribute>\n", indent);


        fprintf(fp, "%s\t</Grid>\n", indent);
    }

    fprintf(fp, "%s</Grid>\n", indent);
}


static FILE *ts_fp = NULL;

static void
writeXDMFTimeSeriesHeader(
    )
{
    if (ts_fp == NULL)
    {
        char fn[PATH_MAX] = { 0 };
        snprintf(fn, PATH_MAX - 1, "%s_ts.xmf", bp->basefilename);
        ts_fp = fopen(fn, "w");
        if (!ts_fp)
            error("Unable to open file %s for writing.\n", fn);

        fprintf(ts_fp, "<?xml version=\"1.0\" ?>\n");
        fprintf(ts_fp, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
        fprintf(ts_fp,
                "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.0\">\n");
        fprintf(ts_fp, "\t<Domain>\n");
        fprintf(ts_fp, "\t\t<Information Name=\"Description\">\n");
        fprintf(ts_fp, "\t\t\tSimulation time series data\n");
        fprintf(ts_fp, "\t\t</Information>\n");
        fprintf(ts_fp, "\t\t<Information Name=\"Simulation Date\">\n");
        time_t t = time(NULL);
        fprintf(ts_fp, "\t\t\t%s", ctime(&t));
        fprintf(ts_fp, "\t\t</Information>\n");
        fprintf(ts_fp,
                "\t\t<Grid Name=\"Time Series\" GridType=\"Collection\" "
                "CollectionType=\"Temporal\">\n");

    }
}


static size_t *renumber_array = NULL;
static size_t renumber_array_sz = 0;

static size_t *
getRenumberArray(
    )
{
    if (bp->num_grains > renumber_array_sz)
    {
        if (renumber_array)
            xfree(renumber_array);
        xmalloc(renumber_array, size_t, bp->num_grains);
        renumber_array_sz = bp->num_grains;
        for (size_t i = 0; i < bp->num_grains; i++)
            renumber_array[i] = i;
        for (int i = 1; i < bp->num_grains; i++)
        {
            int picidx;
            do
            {
                picidx = rand() % (bp->num_grains);
            }
            while (picidx == 0);
            int tmp = renumber_array[picidx];
            renumber_array[picidx] = renumber_array[i];
            renumber_array[i] = tmp;
        }
    }
    if (!renumber_array)
    {
        xmalloc(renumber_array, size_t, 1);
        renumber_array[0] = 0;
    }
    return renumber_array;
}

static void
writeGlobalData(
    FILE * fp,
    const char *prefix)
{
    /* Grain Randomization */
    if (bp->randomize_grains)
    {
        fprintf(fp,
                "%s<DataItem Name=\"Grain Renumber\" Format=\"HDF\" Dimensions=\"%d\" NumberType=\"UInt\">\n",
                prefix, bp->num_grains);
        fprintf(fp, "%s\tSERIAL:%s:/GrainRenumber", prefix,
                get_local_global_file());
        fprintf(fp, "\n%s</DataItem>\n", prefix);
    }

    /* Grain Rotation Angles */
    const char lbl[] = { 'X', 'Y', 'Z' };
    for (int d = 0; d < 3; d++)
    {
        fprintf(fp,
                "%s<DataItem Name=\"Grain Angle %c\" Format=\"HDF\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\">\n",
                prefix, lbl[d], bp->num_grains);
        fprintf(fp, "%s\tSERIAL:%s:/GrainAngle%c\n", prefix,
                get_local_global_file(), lbl[d]);
        fprintf(fp, "%s</DataItem>\n", prefix);
    }
}


static void
writeGlobalHeavyData(
    void)
{
    hid_t file_id = H5Fcreate(get_global_file(), H5F_ACC_TRUNC,
                              H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0)
        error("Unable to open HDF5 file %s (Error code: %ld\n",
              get_global_file(), file_id);

    hsize_t ngr = bp->num_grains;
    hid_t dataspace = H5Screate_simple(1, &ngr, NULL);
    hid_t dataspace_pl = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(dataspace_pl, 1, &ngr);
    H5Pset_fill_time(dataspace_pl, H5D_FILL_TIME_NEVER);
    H5Pset_shuffle(dataspace_pl);
    H5Pset_deflate(dataspace_pl, 9);

#if __BYTE_ORDER == __LITTLE_ENDIAN
    hid_t h5type = (sizeof(size_t) == 8) ? H5T_STD_U64LE : H5T_STD_U32LE;
#else
    hid_t h5type = (sizeof(size_t) == 8) ? H5T_STD_U64BE : H5T_STD_U32BE;
#endif

    if (bp->randomize_grains)
    {
        hid_t dataset = H5Dcreate(file_id, "/GrainRenumber", H5T_NATIVE_UINT,
                                  dataspace, H5P_DEFAULT, dataspace_pl,
                                  H5P_DEFAULT);
        H5Dwrite(dataset, h5type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                 getRenumberArray());
        H5Dclose(dataset);
    }

    // Grain Rotation
    const char lbl[] = { 'X', 'Y', 'Z' };
    float *angles;
    xmalloc(angles, float,
            ngr);
    for (int d = 0; d < 3; d++)
    {
        for (size_t i = 0; i < ngr; i++)
        {
            angles[i] = (float)grain_cache[i].rotang[d];
        }

        char dsName[256] = { 0 };
        sprintf(dsName, "/GrainAngle%c", lbl[d]);

        hid_t dataset = H5Dcreate(file_id, dsName, H5T_NATIVE_FLOAT,
                                  dataspace, H5P_DEFAULT, dataspace_pl,
                                  H5P_DEFAULT);
        H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                 angles);
        H5Dclose(dataset);
    }
    xfree(angles);

    H5Fclose(file_id);

    profile(HDF_OUTPUT);
}


static void
closeXDMFTimeSeries(
    )
{
    if (ts_fp)
    {
        fclose(ts_fp);
        ts_fp = NULL;
    }
}

static void
writeXDMFTimeSeries(
    uint64_t timestamp)
{
    writeXDMFTimeSeriesHeader();

    writeXDMFGrids(ts_fp, timestamp, "\t\t\t");
    fpos_t pos;
    fgetpos(ts_fp, &pos);
    fprintf(ts_fp, "\t\t</Grid>\n");
    writeGlobalData(ts_fp, "\t\t");
    fprintf(ts_fp, "\t</Domain>\n");
    fprintf(ts_fp, "</Xdmf>\n");
    fflush(ts_fp);
    fsetpos(ts_fp, &pos);

}

#ifdef WRITE_TIMESTEPED_XMF
static void
writeXDMF(
    uint64_t timestamp)
{
    const char *filename = get_file_path(timestamp, "xmf");
    FILE *fp = fopen(filename, "w");
    if (!fp)
        error("Unable to open file %s for writing.\n", filename);


    fprintf(fp, "<?xml version=\"1.0\" ?>\n");
    fprintf(fp, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(fp,
            "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.0\">\n");
    fprintf(fp, "\t<Domain>\n");
    fprintf(fp, "\t\t<Information Name=\"Description\">\n");
    fprintf(fp, "\t\t\tSimulation timestep %lu (%g seconds)\n",
            timestamp, timestamp * bp->ts_delt);
    fprintf(fp, "\t\t</Information>\n");

    writeGlobalData(fp, "\t\t");
    writeXDMFGrids(fp, timestamp, "\t\t");

    fprintf(fp, "\t</Domain>\n");
    fprintf(fp, "</Xdmf>\n");
    fclose(fp);
}
#endif



void
xdmf_prepareIO(
    void)
{
    const char *dir = get_rank_dir(iproc);
    int ret = mkdir(dir, 0755);
    if (ret != 0 && errno != EEXIST)
        error("Failed to create directory %s:  %d\n", dir, errno);
}


void
xdmf_writeMain(
    void)
{
    assert(iproc == 0);

    // Create the XML (XMDF) file
#ifdef WRITE_TIMESTEPED_XMF
    writeXDMF(bp->timestep);
#endif
    writeXDMFTimeSeries(bp->timestep);

    writeGlobalHeavyData();
}


void
xdmf_writeSubblocks(
    void)
{
    hid_t file_id = H5Fcreate(get_rank_file(iproc, bp->timestep),
                              H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0)
        error("Unable to open HDF5 file %s (Error code: %ld\n",
              get_rank_file(iproc, bp->timestep), file_id);

    hsize_t start[3] = { 1, 1, 1 };
    hsize_t sbdims[3];
    sbdims[0] = bp->gsdimz + 2;
    sbdims[1] = bp->gsdimy + 2;
    sbdims[2] = bp->gsdimx + 2;
    hsize_t count[3];
    count[0] = bp->gsdimz;
    count[1] = bp->gsdimy;
    count[2] = bp->gsdimx;

    int totaldim = (bp->gsdimx + 2) * (bp->gsdimy + 2) * (bp->gsdimz + 2);

    // TODO:  Find best chunk size
    hsize_t chunkdims[3];
    chunkdims[0] = bp->gsdimz;
    chunkdims[1] = bp->gsdimy;
    chunkdims[2] = bp->gsdimx;

    hid_t dataspace = H5Screate_simple(3, count, NULL);
    hid_t link_pl = H5Pcreate(H5P_LINK_CREATE);
    H5Pset_create_intermediate_group(link_pl, 1 /* Positive means YES */ );

    hid_t dataspace_pl = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(dataspace_pl, 3, chunkdims);
    H5Pset_fill_time(dataspace_pl, H5D_FILL_TIME_NEVER);
    H5Pset_shuffle(dataspace_pl);
    H5Pset_deflate(dataspace_pl, 9);

    hid_t memspace = H5Screate_simple(3, sbdims, NULL);
    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, start, NULL, count, NULL);

    {
        assert(lsp != NULL);
        SB_struct *sb = lsp;

         /***** TEMPERATURE *****/
#if defined(GPU_OMP)
        float *temperature_io = sb->temperature_io;
#pragma omp target update from(temperature_io[0:totaldim])
        profile(OFFLOADING_IO);
#endif
        hid_t temp = H5Dcreate(file_id, get_subblock_path(sb, "Temperature"),
                               H5T_NATIVE_FLOAT, dataspace, link_pl,
                               dataspace_pl, H5P_DEFAULT);
        H5Dwrite(temp, H5T_NATIVE_FLOAT, memspace, H5S_ALL, H5P_DEFAULT,
                 sb->temperature_io);
        H5Dclose(temp);
        profile(HDF_OUTPUT);


        /***** SOLUTE COMPOSITION *****/
        if (bp->calc_type == DIFFUSION)
        {
#if defined(GPU_OMP)
            float *ce_io = sb->ce_io;
#pragma omp target update from(ce_io[0:totaldim])
            profile(OFFLOADING_IO);
#endif

            hid_t fce = H5Dcreate(file_id, get_subblock_path(sb, "CE"),
                                  H5T_NATIVE_FLOAT, dataspace, link_pl,
                                  dataspace_pl, H5P_DEFAULT);
            H5Dwrite(fce, H5T_NATIVE_FLOAT, memspace, H5S_ALL, H5P_DEFAULT,
                     sb->ce_io);
            H5Dclose(fce);
            profile(HDF_OUTPUT);
        }

        /***** GRAIN *****/
#if defined(GPU_OMP)
            int *gr_io = sb->gr_io;
#pragma omp target update from(gr_io[0:totaldim])
            profile(OFFLOADING_IO);
#endif
        hid_t grain = H5Dcreate(file_id, get_subblock_path(sb, "Grain"),
                                H5T_NATIVE_INT, dataspace, link_pl,
                                dataspace_pl, H5P_DEFAULT);
        H5Dwrite(grain, H5T_NATIVE_INT, memspace, H5S_ALL, H5P_DEFAULT,
                 sb->gr_io);
        H5Dclose(grain);
        profile(HDF_OUTPUT);

        /***** Fraction Solid *****/
#if defined(GPU_OMP)
            float *fs_io = sb->fs_io;
#pragma omp target update from(fs_io[0:totaldim])
            profile(OFFLOADING_IO);
#endif
        hid_t ffs = H5Dcreate(file_id, get_subblock_path(sb, "FracSolid"),
                             H5T_NATIVE_FLOAT, dataspace, link_pl,
                             dataspace_pl, H5P_DEFAULT);
        H5Dwrite(ffs, H5T_NATIVE_FLOAT, memspace, H5S_ALL, H5P_DEFAULT,
                 sb->fs_io);
        H5Dclose(ffs);
        profile(HDF_OUTPUT);
    }


    H5Fclose(file_id);

    H5Pclose(dataspace_pl);
    H5Pclose(link_pl);
    H5Sclose(dataspace);
    H5Sclose(memspace);
}

void
xdmf_closeIO(
    void)
{
    closeXDMFTimeSeries();
    if (renumber_array)
        xfree(renumber_array);
}
