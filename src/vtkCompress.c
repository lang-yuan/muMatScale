/***************************************************************/
/* Copyright (c) 2023, Lang Yuan, Univeristy of South Carolina */
/* All rights reserved.                                        */
/* This file is part of muMatScale.                            */
/* See the top-level LICENSE file for details.                 */
/***************************************************************/

#include <stdio.h>
#include <sys/types.h>
#include <inttypes.h>
#include <zlib.h>
#include "vtkCompress.h"
#include "xmalloc.h"
#include "debug.h"

static const unsigned char encodeTable[] =
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ" "abcdefghijklmnopqrstuvwxyz" "0123456789+/";


static void
encodeTriplet(
    uint8_t * d,
    size_t len,
    char *c0,
    char *c1,
    char *c2,
    char *c3)
{
    if (len >= 3)
    {
        *c0 = encodeTable[(d[0] >> 2) & 0x3F];
        *c1 = encodeTable[((d[0] << 4) & 0x30) | ((d[1] >> 4) & 0x0F)];
        *c2 = encodeTable[((d[1] << 2) & 0x3C) | ((d[2] >> 6) & 0x03)];
        *c3 = encodeTable[d[2] & 0x3F];
    }
    else if (len == 2)
    {
        *c0 = encodeTable[(d[0] >> 2) & 0x3F];
        *c1 = encodeTable[((d[0] << 4) & 0x30) | ((d[1] >> 4) & 0x0F)];
        *c2 = encodeTable[(d[1] << 2) & 0x3C];
        *c3 = '=';
    }
    else if (len == 1)
    {
        *c0 = encodeTable[(d[0] >> 2) & 0x3F];
        *c1 = encodeTable[((d[0] << 4) & 0x30)];
        *c2 = '=';
        *c3 = '=';
    }
    else
    {
        error("Invalid Length\n");
    }
}


void
writeBase64(
    FILE * stream,
    void *data,
    size_t len)
{
    uint8_t *d = (uint8_t *) data;
    while (len != 0)
    {
        char o0, o1, o2, o3;
        encodeTriplet(d, len, &o0, &o1, &o2, &o3);
        len = (len < 3) ? 0 : (len - 3);
        d += 3;
        fprintf(stream, "%c%c%c%c", o0, o1, o2, o3);
    }
}



void
prepareVtkCompression(
    vtkComp_t * handle)
{
    memset(handle, 0, sizeof(vtkComp_t));
    xmalloc(handle->hdr, uint32_t, 3);
}

#define NBLOCKS(h)  h->hdr[0]
#define BLK_SZ(h)   h->hdr[1]
#define COMP_SZ(h, idx) h->hdr[3+idx]
void
addVtkCompressionBlock(
    vtkComp_t * h,
    void *data,
    uint32_t len)
{
    uint32_t blk_idx = NBLOCKS(h);
    NBLOCKS(h)++;               // # blocks
    if (BLK_SZ(h) == 0)
    {
        BLK_SZ(h) = len;
        h->cmp_len = compressBound(len);
        xmalloc(h->cmp_buf, uint8_t, h->cmp_len);
    }
    else
    {
        assert(BLK_SZ(h) == len);
    }
    xrealloc(h->hdr, uint32_t, 3 + h->hdr[0]);

    uLongf comp_len = h->cmp_len;
    int zlibret = compress(h->cmp_buf, &comp_len, data, len);
    assert(zlibret == Z_OK);

    COMP_SZ(h, blk_idx) = comp_len;
    xrealloc(h->bigbuf, uint8_t, h->buf_sz + comp_len);
    memcpy(h->bigbuf + h->buf_sz, h->cmp_buf, comp_len);
    h->buf_sz += comp_len;
}


void
writeVtkCompressionData(
    vtkComp_t * h,
    FILE * fp)
{
    writeBase64(fp, h->hdr, (3 + NBLOCKS(h)) * sizeof(uint32_t));
    writeBase64(fp, h->bigbuf, h->buf_sz);

    xfree(h->hdr);
    xfree(h->bigbuf);
    xfree(h->cmp_buf);
}
