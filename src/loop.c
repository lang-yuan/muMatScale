/***************************************************************/
/* Copyright (c) 2023, Lang Yuan, Univeristy of South Carolina */
/* All rights reserved.                                        */
/* This file is part of muMatScale.                            */
/* See the top-level LICENSE file for details.                 */
/***************************************************************/

#include "debug.h"
#include "profiler.h"
#include "globals.h"
#include "distribute.h"
#include "setup.h"
#include "grain.h"
#include "xmalloc.h"
#include "temperature.h"
#include "calculate.h"
#include "file_io.h"
#include "checkpoint.h"
#include "growth.h"

#include <time.h>
#include <omp.h>

extern SB_struct *lsp;

static char *output_leader = "";
static char *output_tailer = "\n";
int g_interactive_mode = 0;

/**
 * Writes the header of the status information to the screen.
 * Also sets, based off the interactive mode flag, the line terminator
 * characters.
 */
static void
outputscreenHeader(
    )
{
    printf("\n%10s  %10s   %8s  %8s\n",
           "Wall Time", "Sim Time", "Fraction", "Grain NO");
    printf("%10s  %10s   %8s\n", "(sec)  ", "(sec)  ", "Solid");

    if (g_interactive_mode)
    {
        output_leader = "        \r";
        output_tailer = "";
    }
    else
    {
        output_leader = "";
        output_tailer = "\n";
    }
}

/**
 * Writes the current status to the screen
 */
void
outputscreen(
    double gStartTime, double fs)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    double t_end = tv.tv_sec + ((double) tv.tv_usec * 1.0e-6);

    printf("%s%10g  %10g  %8.4lf%% %8d  %s \n", output_leader, t_end - gStartTime,      //elapsed wallclock time
           bp->timestep * bp->ts_delt,  //elapsed virtual time
           100.0 * fs,      //solid fraction attained
           bp->num_grains - 1, output_tailer);
    fflush(stdout);
}

// copy data into (lower precision) arrays for IO
void
cache_io_data()
{
    const int dimx = bp->gsdimx;
    const int dimy = bp->gsdimy;
    const int dimz = bp->gsdimz;

    int totaldim = (dimx + 2) * (dimy + 2) * (dimz + 2);

    double* fs = lsp->fs;
    float* fs_io = lsp->fs_io;
#if defined(GPU_OMP)
#pragma omp target teams distribute parallel for schedule(static,1)
#endif
    for (int i = 0; i < totaldim; i++)
        fs_io[i] = (float)fs[i];

    double* ce = lsp->ce;
    float* ce_io = lsp->ce_io;
#if defined(GPU_OMP)
#pragma omp target teams distribute parallel for schedule(static,1)
#endif
    for (int i = 0; i < totaldim; i++)
        ce_io[i] = (float)ce[i];

    double* temperature = lsp->temperature;
    float* temperature_io = lsp->temperature_io;
#if defined(GPU_OMP)
#pragma omp target teams distribute parallel for schedule(static,1)
#endif
    for (int i = 0; i < totaldim; i++)
        temperature_io[i] = (float)temperature[i];

    int *gr = lsp->gr;
    int *gr_io = lsp->gr_io;
#if defined(GPU_OMP)
#pragma omp target teams distribute parallel for schedule(static,1)
#endif
    for (int i = 0; i < totaldim; i++)
        gr_io[i] = gr[i];
}

void
writeData(int timestep)
{
    // Write out the visualization files for each subblock
    if (iproc == 0)
        writeMain(timestep);
    writeSubblocks(timestep);

    MPI_Barrier(mpi_comm_new);
    profile(PROF_OUTPUT);
    timing(FILEIO, timer_elapsed());
}

void
writeCheckpoint()
{
    printf("Checkpoint at timestep %lu \n", bp->timestep);
    if (iproc == 0)
    {
        writeMainCheckpoint();
    }
    writeTaskCheckpoint(bp->timestep);

    MPI_Barrier(mpi_comm_new);
    clean_checkpoint();
    MPI_Barrier(mpi_comm_new);
    profile(PROF_CHECKPOINT);
    timing(FILEIO, timer_elapsed());
}

int
checkDone(const double fs)
{
    int check = 0;

    /* Check for exit conditions */
    if (iproc == 0)
    {
        if (bp->finish_time
            && ((bp->timestep * bp->ts_delt) >= bp->finish_time))
        {
            check = 1;
            printf("\nAuto-Finish time met.\n");
        }
        if (fs >= bp->fs_finish)
        {
            printf("\nSolid Fraction termination condition met.\n");
            check = 1;
        }
    }
    MPI_Bcast(&check, 1, MPI_INT, 0, mpi_comm_new);
    return check;
}

double
computeFS()
{
    double gsolid_volume = solid_volume(lsp);
    double svol = -1.;
    MPI_Reduce(&gsolid_volume, &svol, 1, MPI_DOUBLE,
               MPI_SUM, 0, mpi_comm_new);

    profile(REDUCE_FS);

    return svol / totalNonMoldVolume();
}

// main work loop
void
loop(
    uint64_t restart,
    double gStartTime)
{
    if (iproc == 0)
    {
        printf("Initial Subblocks assigned.  Starting calculations.\n");
        time_t now = time(NULL);
        printf("Current time: %s\n", ctime(&now));
        fflush(stdout);
    }
    timing(COMPUTATION, timer_elapsed());

    int write_step = 0;
    if (bp->data_write_freq > 0 && !restart)
    {
        cache_io_data();
        writeData(bp->timestep);
        write_step = bp->timestep;
    }

    if (bp->screenpfreq > 0 && iproc == 0)
    {
        outputscreenHeader();
        outputscreen(gStartTime,0.);
    }

    int done = 0;
    double fs = 0.;
    int io_time_step = bp->timestep;

#pragma omp parallel shared(done,fs,io_time_step,write_step)
{
    #pragma omp single
    while (!done)
    {
        //printf("bp->timestep = %d\n",bp->timestep);
        if (iproc == 0)
            dwrite(DEBUG_MAIN_CTRL,
                   "-------------------------------------------------\n");
        if (bp->timestep >= bp->data_write_start)
        if (bp->data_write_freq > 0)
        if (bp->timestep % bp->data_write_freq == 0)
        {
            //printf("cache_io_data...\n");
            cache_io_data();
            io_time_step = bp->timestep;
        }

        // compute task
        #pragma omp task
        {
            //printf("thread % d compute...\n",omp_get_thread_num());
            for(int i = 0; i<bp->screenpfreq;i++)
            {
                // Advance to computing the next timestep
                bp->timestep++;

                doiteration();
            }
        }

        // IO task
        // Warning: may involve some collective MPI calls that may interfere
        //          with MPI calls from other threads
        #pragma omp task
        if (io_time_step>write_step)
        {
            //printf("thread %d writeData...\n",omp_get_thread_num());
            writeData(io_time_step);
            write_step = io_time_step;
        }

        #pragma omp taskwait

        if (bp->checkpointfreq > 0 )
        if (bp->timestep % bp->checkpointfreq == 0)
        {
            writeCheckpoint();
        }

        // check for termination should be done by one thread only
        {
            fs = computeFS();
            if( iproc == 0)
                outputscreen(gStartTime,fs);

            //rintf("thread %d, fs = %le\n", omp_get_thread_num(), fs);
            done = checkDone(fs);
        }
    }                           // end of simulation main loop

} // end omp parallel

    // Output final solution (only if we haven't already)
    if (bp->data_write_freq > 0 && bp->timestep > write_step)
    {
        cache_io_data();
        writeData(bp->timestep);
    }

    output_grains(lsp);
    grainShutdown();

    clearSubblocks();

    timing(COMPUTATION, timer_elapsed());

    closeIO();
    timing(FILEIO, timer_elapsed());

    if (iproc == 0)
    {
        printf("\nSimulation complete.  Shutting down.\n");
    }
}
