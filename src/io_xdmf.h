/***************************************************************/
/* Copyright (c) 2023, Lang Yuan, Univeristy of South Carolina */
/* All rights reserved.                                        */
/* This file is part of muMatScale.                            */
/* See the top-level LICENSE file for details.                 */
/***************************************************************/

#ifndef __IO_XDMF_H_
#define __IO_XDMF_H_

#include <inttypes.h>

void xdmf_prepareIO(
    void);
void xdmf_writeMain(
    int timestep);
void xdmf_writeSubblocks(
    int timestep);
void xdmf_closeIO(
    void);

#endif
