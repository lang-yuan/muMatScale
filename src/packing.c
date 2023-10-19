/***************************************************************/
/* Copyright (c) 2023, Lang Yuan, Univeristy of South Carolina */
/* All rights reserved.                                        */
/* This file is part of muMatScale.                            */
/* See the top-level LICENSE file for details.                 */
/***************************************************************/

#include "globals.h"

#include <stdlib.h>

// stride: distance between begining of two blocks of data
// bsize: number of int/double per block of data
void
pack_double(
    double *data,
    const int stride,
    const int bsize,
    const int nblocks,
    const int offset,
    double *buffer)
{
#ifdef GPU_PACK
#pragma omp target teams distribute parallel for
#endif
    for (int i = 0; i < nblocks; i++)
        for (int j = 0; j < bsize; j++)
            buffer[i * bsize + j] = data[offset + i * stride + j];

#ifdef GPU_PACK
#pragma omp target update from(buffer[0:nblocks*bsize])
#endif
}

void
pack_int(
    int *data,
    const int stride,
    const int bsize,
    const int nblocks,
    const int offset,
    int *buffer)
{
#ifdef GPU_PACK
#pragma omp target teams distribute parallel for
#endif
    for (int i = 0; i < nblocks; i++)
        for (int j = 0; j < bsize; j++)
            buffer[i * bsize + j] = data[offset + i * stride + j];

#ifdef GPU_PACK
#pragma omp target update from(buffer[0:nblocks*bsize])
#endif
}

void
pack_3double(
    double *data,
    const int stride,
    const int bsize,
    const int nblocks,
    const int offset,
    double *buffer)
{
    const int offset3 = 3 * offset;
    const int stride3 = 3 * stride;
    const int bsize3 = 3 * bsize;
#ifdef GPU_PACK
#pragma omp target teams distribute parallel for
#endif
    for (int i = 0; i < nblocks; i++)
        for (int j = 0; j < bsize3; j++)
        {
            buffer[i * bsize3 + j] = data[offset3 + i * stride3 + j];
        }

#ifdef GPU_PACK
#pragma omp target update from(buffer[0:nblocks*bsize])
#endif
}

void
pack_field(
    const size_t datasize,
    void *data,
    const int stride,
    const int bsize,
    const int nblocks,
    const int offset,
    void *buffer)
{
    switch (datasize)
    {
        case 8:
            pack_double((double *) data, stride, bsize, nblocks, offset,
                        (double *) buffer);
            break;
        case 4:
            pack_int((int *) data, stride, bsize, nblocks, offset,
                     (int *) buffer);
            break;
        case 24:
            pack_3double((double *) data, stride, bsize, nblocks, offset,
                         (double *) buffer);
            break;
        default:
            printf("error: datasize %zu not supported\n", datasize);
            break;
    }
}

void
unpack_double(
    double *data,
    const int stride,
    const int bsize,
    const int nblocks,
    const int offset,
    double *buffer)
{
#ifdef GPU_PACK
#pragma omp target update to(buffer[0:nblocks*bsize])
#pragma omp target teams distribute parallel for
#endif
    for (int i = 0; i < nblocks; i++)
        for (int j = 0; j < bsize; j++)
            data[offset + i * stride + j] = buffer[i * bsize + j];
}

void
unpack_int(
    int *data,
    const int stride,
    const int bsize,
    const int nblocks,
    const int offset,
    int *buffer)
{
#ifdef GPU_PACK
#pragma omp target update to(buffer[0:nblocks*bsize])
#pragma omp target teams distribute parallel for
#endif
    for (int i = 0; i < nblocks; i++)
        for (int j = 0; j < bsize; j++)
            data[offset + i * stride + j] = buffer[i * bsize + j];
}

void
unpack_3double(
    double *data,
    const int stride,
    const int bsize,
    const int nblocks,
    const int offset,
    double *buffer)
{
    const int offset3 = 3 * offset;
    const int stride3 = 3 * stride;
    const int bsize3 = 3 * bsize;

#ifdef GPU_PACK
#pragma omp target update to(buffer[0:3*nblocks*bsize])
#pragma omp target teams distribute parallel for
#endif
    for (int i = 0; i < nblocks; i++)
        for (int j = 0; j < bsize3; j++)
        {
            data[offset3 + i * stride3 + j] = buffer[i * bsize3 + j];
        }
}

void
unpack_field(
    const size_t datasize,
    void *data,
    const int stride,
    const int bsize,
    const int nblocks,
    const int offset,
    void *buffer)
{
    switch (datasize)
    {
        case 8:
            unpack_double((double *) data, stride, bsize, nblocks, offset,
                          (double *) buffer);
            break;
        case 4:
            unpack_int((int *) data, stride, bsize, nblocks, offset,
                       (int *) buffer);
            break;
        case 24:
            unpack_3double((double *) data, stride, bsize, nblocks, offset,
                           (double *) buffer);
            break;
        default:
            printf("error: datasize %zu not supported\n", datasize);
            break;
    }
}

void
computeHaloInfo(
    const int halo,
    int *offset,
    int *stride,
    int *bsize,
    int *nblocks)
{
    int dimx = bp->gsdimx;
    int dimy = bp->gsdimy;
    int dimz = bp->gsdimz;

    switch (halo)
    {
        case FACE_BOTTOM:
            *offset = 0;
            *stride = 1;
            *bsize = (dimx + 2) * (dimy + 2);
            *nblocks = 1;
            break;
        case FACE_TOP:
            *offset = (dimz + 1) * (dimx + 2) * (dimy + 2);
            *stride = 1;
            *bsize = (dimx + 2) * (dimy + 2);
            *nblocks = 1;
            break;
        case FACE_LEFT:
            *offset = 0;
            *stride = dimx + 2;
            *bsize = 1;
            *nblocks = (dimy + 2) * (dimz + 2);
            break;
        case FACE_RIGHT:
            *offset = dimx + 1;
            *stride = (dimx + 2);
            *bsize = 1;
            *nblocks = (dimy + 2) * (dimz + 2);
            break;
        case FACE_FRONT:
            *offset = 0;
            *stride = (dimx + 2) * (dimy + 2);
            *bsize = dimx + 2;
            *nblocks = dimz + 2;
            break;
        case FACE_BACK:
            *offset = (dimx + 2) * (dimy + 1);
            *stride = (dimx + 2) * (dimy + 2);
            *bsize = dimx + 2;
            *nblocks = dimz + 2;
            break;
    }
}

void
computeFaceInfo(
    const int face,
    int *offset,
    int *stride,
    int *bsize,
    int *nblocks)
{
    int dimx = bp->gsdimx;
    int dimy = bp->gsdimy;
    int dimz = bp->gsdimz;

    switch (face)
    {
        case FACE_BOTTOM:
            *offset = (dimx + 2) * (dimy + 2);
            *stride = 1;
            *bsize = (dimx + 2) * (dimy + 2);
            *nblocks = 1;
            break;
        case FACE_TOP:
            *offset = dimz * (dimx + 2) * (dimy + 2);
            *stride = 1;
            *bsize = (dimx + 2) * (dimy + 2);
            *nblocks = 1;
            break;
        case FACE_LEFT:
            *offset = 1;
            *stride = dimx + 2;
            *bsize = 1;
            *nblocks = (dimy + 2) * (dimz + 2);
            break;
        case FACE_RIGHT:
            *offset = dimx;
            *stride = (dimx + 2);
            *bsize = 1;
            *nblocks = (dimy + 2) * (dimz + 2);
            break;
        case FACE_FRONT:
            *offset = dimx + 2;
            *stride = (dimx + 2) * (dimy + 2);
            *bsize = dimx + 2;
            *nblocks = dimz + 2;
            break;
        case FACE_BACK:
            *offset = (dimx + 2) * dimy;
            *stride = (dimx + 2) * (dimy + 2);
            *bsize = dimx + 2;
            *nblocks = dimz + 2;
            break;
    }
}
