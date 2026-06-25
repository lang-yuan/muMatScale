/***************************************************************/
/* Copyright (c) 2023, Lang Yuan, Univeristy of South Carolina */
/* All rights reserved.                                        */
/* This file is part of muMatScale.                            */
/* See the top-level LICENSE file for details.                 */
/***************************************************************/

#include "globals.h"
#include "profiler.h"
#include "packing.h"

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
#pragma omp target teams distribute parallel for simd collapse(2)
#endif
    for (int i = 0; i < nblocks; i++)
        for (int j = 0; j < bsize; j++)
            buffer[i * bsize + j] = data[offset + i * stride + j];

    profile(PACKING);

#ifdef GPU_PACK
#ifdef CPU_MPI
#pragma omp target update from(buffer[0:nblocks*bsize])
#endif
    profile(PACKING_GPU_CPU);
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
#pragma omp target teams distribute parallel for simd collapse(2)
#endif
    for (int i = 0; i < nblocks; i++)
        for (int j = 0; j < bsize; j++)
            buffer[i * bsize + j] = data[offset + i * stride + j];

    profile(PACKING);

#ifdef GPU_PACK
#ifdef CPU_MPI
#pragma omp target update from(buffer[0:nblocks*bsize])
#endif
    profile(PACKING_GPU_CPU);
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
#pragma omp target teams distribute parallel for simd collapse(2)
#endif
    for (int i = 0; i < nblocks; i++)
        for (int j = 0; j < bsize3; j++)
        {
            buffer[i * bsize3 + j] = data[offset3 + i * stride3 + j];
        }

    profile(PACKING);

#ifdef GPU_PACK
#ifdef CPU_MPI
#pragma omp target update from(buffer[0:nblocks*bsize3])
#endif
    profile(PACKING_GPU_CPU);
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

static void
pack_faces_double(
    double *data,
    const int face_count,
    const int faces[NUM_NEIGHBORS],
    const int strides[NUM_NEIGHBORS],
    const int bsizes[NUM_NEIGHBORS],
    const int nblocks[NUM_NEIGHBORS],
    const int offsets[NUM_NEIGHBORS],
    const int buffer_slot_cells,
    double *buffer)
{
    if (face_count <= 0)
        return;

#ifdef GPU_PACK
#pragma omp target teams distribute parallel for collapse(2) schedule(static,1) \
    map(to: faces[0:face_count], strides[0:face_count], \
            bsizes[0:face_count], nblocks[0:face_count], offsets[0:face_count])
#endif
    for (int f = 0; f < face_count; f++)
    {
        for (int elem = 0; elem < buffer_slot_cells; elem++)
        {
            int elem_count = nblocks[f] * bsizes[f];
            if (elem < elem_count)
            {
                int block = elem / bsizes[f];
                int j = elem - block * bsizes[f];
                buffer[faces[f] * buffer_slot_cells + elem] =
                    data[offsets[f] + block * strides[f] + j];
            }
        }
    }

    profile(PACKING);
#ifdef GPU_PACK
    profile(PACKING_GPU_CPU);
#endif
}

static void
pack_faces_int(
    int *data,
    const int face_count,
    const int faces[NUM_NEIGHBORS],
    const int strides[NUM_NEIGHBORS],
    const int bsizes[NUM_NEIGHBORS],
    const int nblocks[NUM_NEIGHBORS],
    const int offsets[NUM_NEIGHBORS],
    const int buffer_slot_cells,
    int *buffer)
{
    if (face_count <= 0)
        return;

#ifdef GPU_PACK
#pragma omp target teams distribute parallel for collapse(2) schedule(static,1) \
    map(to: faces[0:face_count], strides[0:face_count], \
            bsizes[0:face_count], nblocks[0:face_count], offsets[0:face_count])
#endif
    for (int f = 0; f < face_count; f++)
    {
        for (int elem = 0; elem < buffer_slot_cells; elem++)
        {
            int elem_count = nblocks[f] * bsizes[f];
            if (elem < elem_count)
            {
                int block = elem / bsizes[f];
                int j = elem - block * bsizes[f];
                buffer[faces[f] * buffer_slot_cells + elem] =
                    data[offsets[f] + block * strides[f] + j];
            }
        }
    }

    profile(PACKING);
#ifdef GPU_PACK
    profile(PACKING_GPU_CPU);
#endif
}

static void
pack_faces_3double(
    double *data,
    const int face_count,
    const int faces[NUM_NEIGHBORS],
    const int strides[NUM_NEIGHBORS],
    const int bsizes[NUM_NEIGHBORS],
    const int nblocks[NUM_NEIGHBORS],
    const int offsets[NUM_NEIGHBORS],
    const int buffer_slot_cells,
    double *buffer)
{
    if (face_count <= 0)
        return;

    const int buffer_slot_elems = 3 * buffer_slot_cells;
#ifdef GPU_PACK
#pragma omp target teams distribute parallel for collapse(2) schedule(static,1) \
    map(to: faces[0:face_count], strides[0:face_count], \
            bsizes[0:face_count], nblocks[0:face_count], offsets[0:face_count])
#endif
    for (int f = 0; f < face_count; f++)
    {
        for (int elem = 0; elem < buffer_slot_elems; elem++)
        {
            int bsize3 = 3 * bsizes[f];
            int elem_count = nblocks[f] * bsize3;
            if (elem < elem_count)
            {
                int block = elem / bsize3;
                int j = elem - block * bsize3;
                buffer[faces[f] * buffer_slot_elems + elem] =
                    data[3 * offsets[f] + block * 3 * strides[f] + j];
            }
        }
    }

    profile(PACKING);
#ifdef GPU_PACK
    profile(PACKING_GPU_CPU);
#endif
}

void
pack_faces_field(
    const size_t datasize,
    void *data,
    const int face_count,
    const int faces[NUM_NEIGHBORS],
    const int strides[NUM_NEIGHBORS],
    const int bsizes[NUM_NEIGHBORS],
    const int nblocks[NUM_NEIGHBORS],
    const int offsets[NUM_NEIGHBORS],
    const int buffer_slot_cells,
    void *buffer)
{
    switch (datasize)
    {
        case 8:
            pack_faces_double((double *) data, face_count, faces, strides,
                              bsizes, nblocks, offsets, buffer_slot_cells,
                              (double *) buffer);
            break;
        case 4:
            pack_faces_int((int *) data, face_count, faces, strides, bsizes,
                           nblocks, offsets, buffer_slot_cells, (int *) buffer);
            break;
        case 24:
            pack_faces_3double((double *) data, face_count, faces, strides,
                               bsizes, nblocks, offsets, buffer_slot_cells,
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
#ifdef CPU_MPI
#pragma omp target update to(buffer[0:nblocks*bsize])
#endif

    profile(PACKING_CPU_GPU);

#pragma omp target teams distribute parallel for simd collapse(2)
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
#ifdef CPU_MPI
#pragma omp target update to(buffer[0:nblocks*bsize])
#endif

    profile(PACKING_CPU_GPU);

#pragma omp target teams distribute parallel for simd collapse(2)
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
#ifdef CPU_MPI
#pragma omp target update to(buffer[0:3*nblocks*bsize])
#endif

    profile(PACKING_CPU_GPU);

#pragma omp target teams distribute parallel for simd collapse(2)
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
    profile(UNPACKING);
}

static void
unpack_faces_double(
    double *data,
    const int face_count,
    const int faces[NUM_NEIGHBORS],
    const int strides[NUM_NEIGHBORS],
    const int bsizes[NUM_NEIGHBORS],
    const int nblocks[NUM_NEIGHBORS],
    const int offsets[NUM_NEIGHBORS],
    const int buffer_slot_cells,
    double *buffer)
{
    if (face_count <= 0)
        return;

#ifdef GPU_PACK
    profile(PACKING_CPU_GPU);
#pragma omp target teams distribute parallel for collapse(2) schedule(static,1) \
    map(to: faces[0:face_count], strides[0:face_count], \
            bsizes[0:face_count], nblocks[0:face_count], offsets[0:face_count])
#endif
    for (int f = 0; f < face_count; f++)
    {
        for (int elem = 0; elem < buffer_slot_cells; elem++)
        {
            int elem_count = nblocks[f] * bsizes[f];
            if (elem < elem_count)
            {
                int block = elem / bsizes[f];
                int j = elem - block * bsizes[f];
                data[offsets[f] + block * strides[f] + j] =
                    buffer[faces[f] * buffer_slot_cells + elem];
            }
        }
    }
}

static void
unpack_faces_int(
    int *data,
    const int face_count,
    const int faces[NUM_NEIGHBORS],
    const int strides[NUM_NEIGHBORS],
    const int bsizes[NUM_NEIGHBORS],
    const int nblocks[NUM_NEIGHBORS],
    const int offsets[NUM_NEIGHBORS],
    const int buffer_slot_cells,
    int *buffer)
{
    if (face_count <= 0)
        return;

#ifdef GPU_PACK
    profile(PACKING_CPU_GPU);
#pragma omp target teams distribute parallel for collapse(2) schedule(static,1) \
    map(to: faces[0:face_count], strides[0:face_count], \
            bsizes[0:face_count], nblocks[0:face_count], offsets[0:face_count])
#endif
    for (int f = 0; f < face_count; f++)
    {
        for (int elem = 0; elem < buffer_slot_cells; elem++)
        {
            int elem_count = nblocks[f] * bsizes[f];
            if (elem < elem_count)
            {
                int block = elem / bsizes[f];
                int j = elem - block * bsizes[f];
                data[offsets[f] + block * strides[f] + j] =
                    buffer[faces[f] * buffer_slot_cells + elem];
            }
        }
    }
}

static void
unpack_faces_3double(
    double *data,
    const int face_count,
    const int faces[NUM_NEIGHBORS],
    const int strides[NUM_NEIGHBORS],
    const int bsizes[NUM_NEIGHBORS],
    const int nblocks[NUM_NEIGHBORS],
    const int offsets[NUM_NEIGHBORS],
    const int buffer_slot_cells,
    double *buffer)
{
    if (face_count <= 0)
        return;

    const int buffer_slot_elems = 3 * buffer_slot_cells;
#ifdef GPU_PACK
    profile(PACKING_CPU_GPU);
#pragma omp target teams distribute parallel for collapse(2) schedule(static,1) \
    map(to: faces[0:face_count], strides[0:face_count], \
            bsizes[0:face_count], nblocks[0:face_count], offsets[0:face_count])
#endif
    for (int f = 0; f < face_count; f++)
    {
        for (int elem = 0; elem < buffer_slot_elems; elem++)
        {
            int bsize3 = 3 * bsizes[f];
            int elem_count = nblocks[f] * bsize3;
            if (elem < elem_count)
            {
                int block = elem / bsize3;
                int j = elem - block * bsize3;
                data[3 * offsets[f] + block * 3 * strides[f] + j] =
                    buffer[faces[f] * buffer_slot_elems + elem];
            }
        }
    }
}

void
unpack_faces_field(
    const size_t datasize,
    void *data,
    const int face_count,
    const int faces[NUM_NEIGHBORS],
    const int strides[NUM_NEIGHBORS],
    const int bsizes[NUM_NEIGHBORS],
    const int nblocks[NUM_NEIGHBORS],
    const int offsets[NUM_NEIGHBORS],
    const int buffer_slot_cells,
    void *buffer)
{
    switch (datasize)
    {
        case 8:
            unpack_faces_double((double *) data, face_count, faces, strides,
                                bsizes, nblocks, offsets, buffer_slot_cells,
                                (double *) buffer);
            break;
        case 4:
            unpack_faces_int((int *) data, face_count, faces, strides, bsizes,
                             nblocks, offsets, buffer_slot_cells,
                             (int *) buffer);
            break;
        case 24:
            unpack_faces_3double((double *) data, face_count, faces, strides,
                                 bsizes, nblocks, offsets, buffer_slot_cells,
                                 (double *) buffer);
            break;
        default:
            printf("error: datasize %zu not supported\n", datasize);
            break;
    }
    profile(UNPACKING);
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

int
face_is_contiguous_plane(
    const int face)
{
    return face == FACE_BOTTOM || face == FACE_TOP;
}
