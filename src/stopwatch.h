/***************************************************************/
/* Copyright (c) 2023, Lang Yuan, Univeristy of South Carolina */
/* All rights reserved.                                        */
/* This file is part of muMatScale.                            */
/* See the top-level LICENSE file for details.                 */
/***************************************************************/

#ifndef STOPWATCH_H__
#define STOPWATCH_H__
#include <sys/time.h>           //System time

typedef struct stopwatch
{
    struct timeval start_tv;
    struct timeval stop_tv;
    struct timeval temp_tv;
    struct timeval temp2_tv;
    struct timeval lap_tv;
    int started;
    int stopped;
    int lap;
} Stopwatch;

double timer_diff(
    struct timeval *start,
    struct timeval *stop);

/* Stop watch operations */
void timer_init(
    );
void timer_start(
    );
void timer_stop(
    );
void timer_reset(
    );

/* Always get elapsed time since the previous start time
 * Interval since last start time (if the watch is not stopped) - |now - start|
 * OR between last start & stop times (if the watch is stopped) - |stop - start|
 */
double timer_elapsed(
    );

/* Get elapsed time since last call w/o stopping the watch
 * Interval since last start (first time call) - |now - start|
 * OR since last call to lapTime (following calls) - |now - previous|
 */
double timer_lapTime(
    );

extern Stopwatch timer;

#endif
