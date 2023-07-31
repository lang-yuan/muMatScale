/***************************************************************/
/* Copyright (c) 2023, Lang Yuan, Univeristy of South Carolina */
/* All rights reserved.                                        */
/* This file is part of muMatScale.                            */
/* See the top-level LICENSE file for details.                 */
/***************************************************************/

#include <string.h>             //NULL definition
#include "stopwatch.h"

Stopwatch timer;

double
timer_diff(
    struct timeval *start,
    struct timeval *stop)
{
    double d1, d2;
    d1 = ((double) start->tv_sec) * 1e6 + (double) start->tv_usec;      //microseconds
    d2 = ((double) stop->tv_sec) * 1e6 + (double) stop->tv_usec;
    return ((d2 - d1) * 1e-6);  //seconds
}

void
timer_init(
    )
{
    timer.started = 0;
    timer.stopped = 0;
    timer.lap = 0;
}

void
timer_start(
    )
{
    timer.started = 1;
    timer.stopped = 0;
    timer.lap = 0;
    gettimeofday(&timer.start_tv, NULL);
}

void
timer_stop(
    )
{
    if (timer.started)
    {
        gettimeofday(&timer.stop_tv, NULL);
        timer.stopped = 1;
        timer.started = 0;
        timer.lap = 0;
    }
}

void
timer_reset(
    )
{
    timer.stopped = 1;
    timer.started = 0;
    timer.lap = 0;
    timer.start_tv = timer.stop_tv;
}

double
timer_elapsed(
    )
{
    if (timer.stopped)
    {
        return timer_diff(&(timer.start_tv), &(timer.stop_tv));
    }
    else if (timer.started)
    {
        gettimeofday(&(timer.temp_tv), NULL);
        return timer_diff(&(timer.start_tv), &(timer.temp_tv));
    }
    else
    {
        return (double) 0;      //err: watch has not started yet
    }
}

double
timer_lapTime(
    )
{
    if (timer.started)
    {
        gettimeofday(&(timer.temp_tv), NULL);
        timer.temp2_tv = timer.lap_tv;
        timer.lap_tv = timer.temp_tv;

        if (timer.lap)
        {
            return timer_diff(&(timer.temp2_tv), &(timer.temp_tv));
        }
        else
        {
            timer.lap = 1;
            return timer_diff(&(timer.start_tv), &(timer.temp_tv));
        }
    }
    else
    {
        return (double) 0;      //err: watch has not started yet
    }
}
