/***************************************************************/
/* Copyright (c) 2023, Lang Yuan, Univeristy of South Carolina */
/* All rights reserved.                                        */
/* This file is part of muMatScale.                            */
/* See the top-level LICENSE file for details.                 */
/***************************************************************/

#ifndef __IO_VTK_H__
#define __IO_VTK_H__

#include <inttypes.h>
#include "globals.h"

void vtk_writeMain(
    );
void vtk_prepareIOSubblock(
    SB_struct * sb);
void vtk_writeSubblocks(
    );


#endif
