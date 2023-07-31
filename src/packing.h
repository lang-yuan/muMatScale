/***************************************************************/
/* Copyright (c) 2023, Lang Yuan, Univeristy of South Carolina */
/* All rights reserved.                                        */
/* This file is part of muMatScale.                            */
/* See the top-level LICENSE file for details.                 */
/***************************************************************/

void pack_field(
    const size_t datasize,
    void *data,
    const int stride,
    const int bsize,
    const int nblocks,
    const int offset,
    void *buffer);

void unpack_field(
    const size_t datasize,
    void *data,
    const int stride,
    const int bsize,
    const int nblocks,
    const int offset,
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
