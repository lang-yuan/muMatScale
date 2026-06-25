/***************************************************************/
/* Copyright (c) 2023, Lang Yuan, Univeristy of South Carolina */
/* All rights reserved.                                        */
/* This file is part of muMatScale.                            */
/* See the top-level LICENSE file for details.                 */
/***************************************************************/

#ifndef PACKING_H_
#define PACKING_H_

#include <stddef.h>

#include "globals.h"

void pack_field(
    const size_t datasize,
    void *data,
    const int stride,
    const int bsize,
    const int nblocks,
    const int offset,
    void *buffer);

void pack_faces_field(
    const size_t datasize,
    void *data,
    const int face_count,
    const int faces[NUM_NEIGHBORS],
    const int strides[NUM_NEIGHBORS],
    const int bsizes[NUM_NEIGHBORS],
    const int nblocks[NUM_NEIGHBORS],
    const int offsets[NUM_NEIGHBORS],
    const int buffer_slot_cells,
    void *buffer);

void unpack_field(
    const size_t datasize,
    void *data,
    const int stride,
    const int bsize,
    const int nblocks,
    const int offset,
    void *buffer);

void unpack_faces_field(
    const size_t datasize,
    void *data,
    const int face_count,
    const int faces[NUM_NEIGHBORS],
    const int strides[NUM_NEIGHBORS],
    const int bsizes[NUM_NEIGHBORS],
    const int nblocks[NUM_NEIGHBORS],
    const int offsets[NUM_NEIGHBORS],
    const int buffer_slot_cells,
    void *buffer);

void computeHaloInfo(
    const int halo,
    int *offset,
    int *stride,
    int *bsize,
    int *nblocks);

void computeFaceInfo(
    const int face,
    int *offset,
    int *stride,
    int *bsize,
    int *nblocks);

int face_is_contiguous_plane(
    const int face);

#endif /* PACKING_H_ */
