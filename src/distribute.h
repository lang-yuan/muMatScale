/***************************************************************/
/* Copyright (c) 2023, Lang Yuan, Univeristy of South Carolina */
/* All rights reserved.                                        */
/* This file is part of muMatScale.                            */
/* See the top-level LICENSE file for details.                 */
/***************************************************************/

#ifndef __MAIN_H__
#define __MAIN_H__

#include <inttypes.h>
#include "globals.h"

typedef struct task_s
{
    int rank;
    int assigned_sb;
} task_t;

typedef struct tasks_s
{
    task_t *task;
    nbr_info_t *nbrbuf;
    int *nbrcnts;
    int *nbrdisp;
} tasks_t;


void outputvis(
    );
void init_gmsp(
    );
void sendSubblockInfo(
    );
void completeSendSubblockInfo(
    );
void clearSubblocks(
    );
void recvSubblockInfo(
    );

extern int g_interactive_mode;
extern tasks_t gtasks;

#endif
