/***************************************************************/
/* Copyright (c) 2023, Lang Yuan, Univeristy of South Carolina */
/* All rights reserved.                                        */
/* This file is part of muMatScale.                            */
/* See the top-level LICENSE file for details.                 */
/***************************************************************/

#include <mpi.h>
#include "debug.h"
#include "xmalloc.h"
#include "profiler.h"
#include "globals.h"
#include "functions.h"

static char *bucket_names[] = {
    "Setup",
    "Temp Prep",
    "Idle",
    "Chkpt",
    "Distribution",
    "Output",
    "Face Xchg (local)",
    "Face Xchg SEND",
    "Face Xchg WAIT",
    "Temp Update",
    "Nucleation",
    "FS Change",
    "Grow Octahedra",
    "Cell Index",
    "Capture Octahedra",
    "Diffuse Alloy",
    "Grain Activation",
    "Grain Sync 1",
    "Grain Sync 2",
    "Offloading CPU-GPU",
    "Offloading GPU-CPU",
    "Timer Temp1",
    "Timer Temp2",
    "unPacking",
    "Packing",
    "Initialization",
    "Computation",
    "Communication",
    "Local Exchange",
    "Synchronization",
    "File IO",
};

static double *bucket_times = NULL;
static double init_time = 0.0;
static double last_recorded = 0.0;

void
profiler_init(
    )
{
    /* Verify that our sizes are the same */
    assert(sizeof(bucket_names) / sizeof(char *) == NUM_BUCKETS);

    /* Each process stores timing info into the bucket array entries */
    xmalloc(bucket_times, double,
            NUM_BUCKETS);

    init_time = MPI_Wtime();
    last_recorded = init_time;
}


void
profile(
    bucket_tag tag)
{
    double time_now = MPI_Wtime();
    bucket_times[tag] += (time_now - last_recorded);
    last_recorded = time_now;
}

/* 	[__PROFILE]
	Accumulate the given elapsed time in one of the bucket
	- Initialization, Computation, Communication, and fileIO
*/
void
timing(
    bucket_tag tag,
    double t)
{
    timer_stop();
    bucket_times[tag] += t;
    timer_start();
}

void
profiler_collate(
    )
{
#ifndef PATH_MAX
#define PATH_MAX 512
#endif
    double (
    *allbuckets)[NUM_BUCKETS] = NULL;

    if (iproc == 0)
    {
        allbuckets =
            (double (*)[NUM_BUCKETS]) malloc(sizeof(double[NUM_BUCKETS]) *
                                             nproc);
    }
    MPI_Gather(bucket_times, NUM_BUCKETS, MPI_DOUBLE,
               allbuckets, NUM_BUCKETS, MPI_DOUBLE, 0, mpi_comm_new);

    if (iproc == 0)
    {
        char profiler_file[PATH_MAX] = { 0 };
        snprintf(profiler_file, PATH_MAX - 1, "%s_profile_%d.csv",
                 bp->basefilename, bp->timestep);
        FILE *fp = fopen(profiler_file, "w");

        int j = 0;
        fprintf(fp, "Rank,Num SB,");
        for (; j < (NUM_BUCKETS - SYS_PROFILE); j++)
        {
            fprintf(fp, "%s,", bucket_names[j]);
        }
        fprintf(fp, "Total,");

        for (; j < NUM_BUCKETS; j++)
        {                       //[__PROFILE]
            fprintf(fp, "%s,", bucket_names[j]);
        }
        fprintf(fp, "System Total\n");


        for (int i = 0; i < nproc; i++)
        {
            j = 0;
            fprintf(fp, "%d,%d,", i, 1);
            double ranksum = 0.0;
            for (; j < (NUM_BUCKETS - SYS_PROFILE); j++)
            {
                ranksum += allbuckets[i][j];
                fprintf(fp, "%lg,", allbuckets[i][j]);
            }
            fprintf(fp, "%lg,", ranksum);
            ranksum = 0.0;
            for (; j < NUM_BUCKETS; j++)
            {                   //[__PROFILE]
                ranksum += allbuckets[i][j];
                fprintf(fp, "%lg,", allbuckets[i][j]);
            }
            fprintf(fp, "%lg\n", ranksum);
        }
        fclose(fp);

        xfree(allbuckets);
    }
}


void
profiler_write_stats(
    )
{
    double accum_time = 0.0;
    double time_now = MPI_Wtime();
    printf("Timing Statistics for rank %d\n", iproc);
    for (int i = 0; i < NUM_BUCKETS; i++)
    {
        accum_time += bucket_times[i];
        printf("%40s     %lg\n", bucket_names[i], bucket_times[i]);
    }
    double total_time = time_now - init_time;
    printf("Total time:  %lg\n", total_time);
    printf("Accumulated time:  %lg\n", accum_time);
    printf("Unaccounted for time:  %lg\n", total_time - accum_time);
}
