/***************************************************************/
/* Copyright (c) 2023, Lang Yuan, Univeristy of South Carolina */
/* All rights reserved.                                        */
/* This file is part of muMatScale.                            */
/* See the top-level LICENSE file for details.                 */
/***************************************************************/

#include <stdio.h>
#include <inttypes.h>
#include <sys/stat.h>
#include <sys/types.h>
#define __USE_BSD   1           /* for symlink */
#include <unistd.h>
#include <errno.h>
#include "globals.h"
#include "functions.h"
#include "debug.h"
#include "ll.h"

/* Disable USE_COMPRESSION (set to 0) to write plain-text Paraview files */
#ifndef USE_COMPRESSION
#define USE_COMPRESSION 1
#endif

#if USE_COMPRESSION
#include "xmalloc.h"
#include "vtkCompress.h"
#endif

#if USE_COMPRESSION
#define VTK_FORMAT "binary"
#else
#define VTK_FORMAT "ascii"
#endif

extern SB_struct *lsp;

static void
writepvtu(
    BB_struct * b,
    int ts)
{
    char filename[256];

    sprintf(filename, "%s_%d.pvtu", bp->basefilename, ts);
    FILE *fp = fopen(filename, "w");
    if (fp == NULL)
    {
        error("%d: Unable to open pvtu file\n", iproc);
    }

    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp,
            "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    fprintf(fp, "  <PUnstructuredGrid GhostLevel=\"0\">\n");

    // Cell data
    fprintf(fp, "    <PCellData>\n");
    fprintf(fp,
            "      <PDataArray type=\"Float64\" Name=\"Temperature\"/>\n");
    fprintf(fp, "      <PDataArray type=\"Int32\" Name=\"Grain\"/>\n");
    fprintf(fp, "    </PCellData>\n");

    // Points
    fprintf(fp, "    <PPoints>\n");
#if USE_COMPRESSION
    fprintf(fp,
            "      <PDataArray type=\"Float64\" NumberOfComponents=\"3\"/>\n");
#else
    fprintf(fp,
            "      <PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>\n");
#endif
    fprintf(fp, "    </PPoints>\n");

    // Now reference the pieces coming from the tasks

    int n_sub = bp->gnsbx * bp->gnsby * bp->gnsbz;
    for (int i = 0; i < n_sub; i++)
    {
        {
            fprintf(fp, "    <Piece Source=\"%s_%04d/%d.vtu\"/>\n",
                    bp->basefilename, i, ts);
        }
    }

    // Close the xml tags
    fprintf(fp, "  </PUnstructuredGrid>\n");
    fprintf(fp, "</VTKFile>\n");

    fclose(fp);
}


static void
writeDataArrayFloat(
    FILE * fp,
    double *buf,
    const char *Name)
{
    int dimx = bp->gsdimx;
    int dimy = bp->gsdimy;
    int dimz = bp->gsdimz;
#if USE_COMPRESSION
    double *tbuf;
    vtkComp_t ptHandle;
    uint32_t buf_sz = dimy * dimz;

    xmalloc(tbuf, double,
            buf_sz);
    prepareVtkCompression(&ptHandle);
#endif

    fprintf(fp,
            "        <DataArray type=\"Float64\" Name=\"%s\" format=\"%s\">\n",
            Name, VTK_FORMAT);
    for (int i = 1; i <= dimx; i++)
    {
#if USE_COMPRESSION
        int idx = 0;
#endif
        for (int j = 1; j <= dimy; j++)
        {
            for (int k = 1; k <= dimz; k++)
            {
#if USE_COMPRESSION
                tbuf[idx++] =
                    buf[k * (dimy + 2) * (dimx + 2) + j * (dimx + 2) + i];
#else
                fprintf(fp, "%lf\n",
                        buf[k * (dimy + 2) * (dimx + 2) + j * (dimx + 2) +
                            i]);
#endif
            }
        }
#if USE_COMPRESSION
        addVtkCompressionBlock(&ptHandle, tbuf, buf_sz * sizeof(double));
#endif
    }

#if USE_COMPRESSION
    writeVtkCompressionData(&ptHandle, fp);
    xfree(tbuf);
#endif
    fprintf(fp, "\n        </DataArray>\n");

}


static void
writeDataArrayInt32(
    FILE * fp,
    int *buf,
    const char *Name)
{
    int dimx = bp->gsdimx;
    int dimy = bp->gsdimy;
    int dimz = bp->gsdimz;
#if USE_COMPRESSION
    int32_t *tbuf;
    vtkComp_t ptHandle;
    uint32_t buf_sz = dimy * dimz;

    xmalloc(tbuf, int32_t, buf_sz);
    prepareVtkCompression(&ptHandle);
#endif

    fprintf(fp,
            "        <DataArray type=\"Int32\" Name=\"%s\" format=\"%s\">\n",
            Name, VTK_FORMAT);
    for (int i = 1; i <= dimx; i++)
    {
#if USE_COMPRESSION
        int idx = 0;
#endif
        for (int j = 1; j <= dimy; j++)
        {
            for (int k = 1; k <= dimz; k++)
            {
#if USE_COMPRESSION
                tbuf[idx++] =
                    buf[k * (dimy + 2) * (dimx + 2) + j * (dimx + 2) + i];
#else
                fprintf(fp, "%d\n",
                        buf[k * (dimy + 2) * (dimx + 2) + j * (dimx + 2) +
                            i]);
#endif
            }
        }
#if USE_COMPRESSION
        addVtkCompressionBlock(&ptHandle, tbuf, buf_sz * sizeof(int32_t));
#endif
    }

#if USE_COMPRESSION
    writeVtkCompressionData(&ptHandle, fp);
    xfree(tbuf);
#endif
    fprintf(fp, "\n        </DataArray>\n");

}

static void
writevtu(
    SB_struct * sb,
    BB_struct * bb,
    int ts)
{
    int dimx = bb->gsdimx;
    int dimy = bb->gsdimy;
    int dimz = bb->gsdimz;
    char filename[256];
#if USE_COMPRESSION
    vtkComp_t ptHandle;
    int32_t *int32_buf;
    int8_t *int8_buf;
    double *float_buf;
#endif
#if 0
    alert
        ("ts %d sb %d neighbor map is now: %d, %d, %d, %d, %d, %d\n",
         ts, sb->subblockid, sb->neighbors[0][0],,
         sb->neighbors[1][0], sb->neighbors[2][0],
         sb->neighbors[3][0], sb->neighbors[4][0], sb->neighbors[5][0]);
#endif

    // The directory was created back in 'SetupInitialConditions'
    sprintf(filename, "%s_%04d/%d.vtu", bp->basefilename, sb->subblockid, ts);
    FILE *fp = fopen(filename, "w");
    if (fp == NULL)
    {
        error("%d: Unable to open vtu file\n", iproc);
    }

    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp,
            "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n");
    fprintf(fp, "  <UnstructuredGrid>\n");

    int ncells = bb->gsdimx * bb->gsdimy * bb->gsdimz;
    int npoints = (bb->gsdimx + 1) * (bb->gsdimy + 1) * (bb->gsdimz + 1);
    fprintf(fp, "      <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\" >\n",
            npoints, ncells);




    // Cell data
    fprintf(fp, "        <CellData>\n");
    writeDataArrayFloat(fp, sb->temperature, "Temperature");
    writeDataArrayInt32(fp, sb->gr, "Grain");
    fprintf(fp, "      </CellData>\n");



    // This sets up the cells

    // Now write out the point coordinates
    fprintf(fp, "      <Points>\n");
#if USE_COMPRESSION
    fprintf(fp,
            "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"binary\">\n");
    prepareVtkCompression(&ptHandle);
    uint32_t coords_sz = sizeof(double) * 3 * (dimy + 1) * (dimz + 1);
    xmalloc(float_buf, double,
            coords_sz / sizeof(double));
    for (int i = 0; i < dimx + 1; i++)
    {
        int idx = 0;
        for (int j = 0; j < dimy + 1; j++)
        {
            for (int k = 0; k < dimz + 1; k++)
            {
                double x, y, z;

                // Translate the subblock float_buf into real float_buf based on my subblock id
                scoord2realcoord(sb->coords.x, sb->coords.y, sb->coords.z, i,
                                 j, k, &x, &y, &z);
                float_buf[idx++] = x;
                float_buf[idx++] = y;
                float_buf[idx++] = z;
            }
        }

        addVtkCompressionBlock(&ptHandle, float_buf, coords_sz);
    }
    writeVtkCompressionData(&ptHandle, fp);
    xfree(float_buf);
#else
    fprintf(fp,
            "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    // write the point float_buf
    for (int i = 0; i < dimx + 1; i++)
    {
        for (int j = 0; j < dimy + 1; j++)
        {
            for (int k = 0; k < dimz + 1; k++)
            {
                float x, y, z;

                // Translate the subblock float_buf into real float_buf based on my subblock id
                scoord2realcoord(sb, i, j, k, &x, &y, &z);
                fprintf(fp, "%g %g %g\n", x, y, z);
            }
        }
    }
#endif

    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "      </Points>\n");

    fprintf(fp, "      <Cells>\n");
#if USE_COMPRESSION
    fprintf(fp,
            "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">\n");
    // Write out the extents of the cell in point indices - indexed into the above point float_buf

    prepareVtkCompression(&ptHandle);
    uint32_t cb_sz = dimy * dimz * 8;
    xmalloc(int32_buf, int32_t, cb_sz);

    int rowsize = dimz + 1;
    int planesize = (dimy + 1) * (dimz + 1);
    //int rowsize = dimx+1;
    //int planesize = (dimx+1) * (dimy+1);

    for (int i = 0; i < dimx; i++)
    {
        int idx = 0;
        for (int j = 0; j < dimy; j++)
        {
            for (int k = 0; k < dimz; k++)
            {
                int cell0 = i * planesize + j * rowsize + k;
                int32_buf[idx++] = cell0;
                int32_buf[idx++] = cell0 + 1;
                int32_buf[idx++] = cell0 + rowsize;
                int32_buf[idx++] = cell0 + rowsize + 1;
                int32_buf[idx++] = cell0 + planesize;
                int32_buf[idx++] = cell0 + planesize + 1;
                int32_buf[idx++] = cell0 + planesize + rowsize;
                int32_buf[idx++] = cell0 + planesize + rowsize + 1;
            }
        }
        addVtkCompressionBlock(&ptHandle, int32_buf, cb_sz * sizeof(int32_t));
    }
    writeVtkCompressionData(&ptHandle, fp);
    xfree(int32_buf);
#else
    fprintf(fp,
            "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
    // Write out the extents of the cell in point indices - indexed into the above point coords

    int rowsize = dimz + 1;
    int planesize = (dimy + 1) * (dimz + 1);
    //int rowsize = dimx+1;
    //int planesize = (dimx+1) * (dimy+1);

    for (int i = 0; i < dimx; i++)
    {
        for (int j = 0; j < dimy; j++)
        {
            for (int k = 0; k < dimz; k++)
            {
                //           int cell0 = i + j * rowsize + k * planesize;
                int cell0 = i * planesize + j * rowsize + k;
                fprintf(fp, "%d %d %d %d %d %d %d %d\n", cell0,
                        cell0 + 1,
                        cell0 + rowsize,
                        cell0 + rowsize + 1,
                        cell0 + planesize,
                        cell0 + planesize + 1,
                        cell0 + planesize + rowsize,
                        cell0 + planesize + rowsize + 1);
            }
        }
    }
#endif
    fprintf(fp, "        </DataArray>\n");
#if USE_COMPRESSION
    fprintf(fp,
            "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">\n");
    prepareVtkCompression(&ptHandle);
    xmalloc(int32_buf, int32_t, ncells);

    for (int n = 1; n <= ncells; n++)
    {
        int32_buf[n - 1] = n * 8;
    }
    addVtkCompressionBlock(&ptHandle, int32_buf, ncells * sizeof(int32_t));
    writeVtkCompressionData(&ptHandle, fp);
    xfree(int32_buf);

#else
    fprintf(fp,
            "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");

    for (int n = 1; n <= ncells; n++)
    {
        fprintf(fp, "%d ", n * 8);
    }
    fprintf(fp, "\n");

#endif
    fprintf(fp, "        </DataArray>\n");


    // VTU types
    // =====  =================== ============= ===
    // keys   type                n points      dim
    // =====  =================== ============= ===
    //   1   VTK_VERTEX:         1 point        2d
    //   2   VTK_POLY_VERTEX:    n points       2d
    //   3   VTK_LINE:           2 points       2d
    //   4   VTK_POLY_LINE:      n+1 points     2d
    //   5   VTK_TRIANGLE:       3 points       2d
    //   6   VTK_TRIANGLE_STRIP: n+2 points     2d
    //   7   VTK_POLYGON:        n points       2d
    //   8   VTK_PIXEL:          4 points       2d
    //   9   VTK_QUAD:           4 points       2d
    //   10  VTK_TETRA:          4 points       3d
    //   11  VTK_VOXEL:          8 points       3d
    //   12  VTK_HEXAHEDRON:     8 points       3d
    //   13  VTK_WEDGE:          6 points       3d
    //   14  VTK_PYRAMID:        5 points       3d
    // =====  =================== ============= ===

#if USE_COMPRESSION
    fprintf(fp,
            "        <DataArray type=\"UInt8\" Name=\"types\" format=\"binary\">\n");
    prepareVtkCompression(&ptHandle);
    xmalloc(int8_buf, int8_t, ncells);
    for (int n = 1; n <= ncells; n++)
    {
        int8_buf[n - 1] = 11;   // This is for VTK_HEXAHEDRON
    }
    addVtkCompressionBlock(&ptHandle, int8_buf, ncells * sizeof(int8_t));
    writeVtkCompressionData(&ptHandle, fp);
    xfree(int8_buf);
#else
    fprintf(fp,
            "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    for (int n = 1; n <= ncells; n++)
    {
        fprintf(fp, "%d ", 11);
    }
    fprintf(fp, "\n");
#endif

    fprintf(fp, "        </DataArray>\n");




    fprintf(fp, "       </Cells>\n");
    fprintf(fp, "    </Piece>\n");
    fprintf(fp, "  </UnstructuredGrid>\n");
    fprintf(fp, "</VTKFile>\n");

    fclose(fp);
}


/*************************** I/O Module functions ************************/

void
vtk_writeMain(
    )
{
    writepvtu(bp, bp->timestep);
}

void
vtk_prepareIOSubblock(
    SB_struct * sb)
{
    char filename[256];
    if (bp->data_write_freq > 0)
    {
        sprintf(filename, "%s_%04d", bp->basefilename, sb->subblockid);
        mkdir(filename, 0755);  // OK to fail, probably already created it.
    }

}

void
vtk_writeSubblocks(
    )
{
    assert(lsp != NULL);

    writevtu(lsp, bp, (int) bp->timestep);
}
