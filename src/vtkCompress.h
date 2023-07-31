/***************************************************************/
/* Copyright (c) 2023, Lang Yuan, Univeristy of South Carolina */
/* All rights reserved.                                        */
/* This file is part of muMatScale.                            */
/* See the top-level LICENSE file for details.                 */
/***************************************************************/

#ifndef VTKCOMPRESS_H_
#define VTKCOMPRESS_H_

#include <stdio.h>
#include <sys/types.h>


void writeBase64(
    FILE * stream,
    void *data,
    size_t len);

typedef struct vtkComp
{
    uint32_t *hdr;
    uint8_t *bigbuf;
    uint32_t buf_sz;
    uint8_t *cmp_buf;
    uint32_t cmp_len;

} vtkComp_t;

void prepareVtkCompression(
    vtkComp_t * handle);
void addVtkCompressionBlock(
    vtkComp_t * handle,
    void *data,
    uint32_t len);
void writeVtkCompressionData(
    vtkComp_t * handle,
    FILE * fp);

#endif
