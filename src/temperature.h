/***************************************************************/
/* Copyright (c) 2023, Lang Yuan, Univeristy of South Carolina */
/* All rights reserved.                                        */
/* This file is part of muMatScale.                            */
/* See the top-level LICENSE file for details.                 */
/***************************************************************/

#ifndef __TEMPERATURE_H__
#define __TEMPERATURE_H__

#include <sys/types.h>
#include <stdbool.h>
#include "globals.h"

#define TS_NEVER (~0ULL)
uint64_t timeToReachTemp(
    volume_t * volume,
    double target);
uint64_t timeToReachTempSBID(
    size_t sbid,
    double target);

void tempPrepareSB(
    SB_struct * sb);
void tempUpdate(
    bool updateAll);

double totalNonMoldVolume(
    void);


#endif
