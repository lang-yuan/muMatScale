/***************************************************************/
/* Copyright (c) 2023, Lang Yuan, Univeristy of South Carolina */
/* All rights reserved.                                        */
/* This file is part of muMatScale.                            */
/* See the top-level LICENSE file for details.                 */
/***************************************************************/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <signal.h>

#include "globals.h"
#include "functions.h"
#include "distribute.h"
#include "xmalloc.h"
#include "setup.h"
#include "profiler.h"
#include "read_ctrl.h"
#include "file_io.h"
#include "grain.h"
#include "ll.h"
#include "checkpoint.h"
#include "loop.h"

#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef USE_MPE
#include <mpe.h>
int gLogStateBeginID;
int gLogStateEndID;
#endif

extern char *optarg;

// list of subblocks for task 0
extern ll_t *gsubblock_list;
extern SB_struct *lsp;

void
print_usage(
    const char *app_name)
{
    fprintf(stderr, "Usage:\n\t%s [-i] [-h] -c <file> [-r <time>]\n\n",
            app_name);
    fprintf(stderr, "\t-i Turns on \"interactive\" mode for screen stats.\n");
    fprintf(stderr, "\t-c <control file>\n");
    fprintf(stderr, "\t-h Show this help.\n");
}

void
signal_handler(
    int signum)
{
    switch (signum)
    {
        case SIGINT:
            alert("\n\n*** Ctrl-C Caught! *** Shutting Down ***\n\n");
            break;
        case SIGUSR2:
            alert
                ("\n\n*** Wall Time Limit Reached!  *** Shutting Down ***\n\n");
            break;
        default:
            alert
                ("\n\n*** Unexpected signal %d Received!   *** Shutting Down ***\n\n",
                 signum);
    }
}

int
main(
    int argc,
    char *argv[])
{
    timer_init();
    timer_start();

    signal(SIGINT, signal_handler);
    signal(SIGUSR2, signal_handler);

    int retval = 0;
    uint64_t restart = 0;
    int flag;
    char *ctrl_fname = NULL;

    xmalloc(bp, BB_struct, 1);

    gMPIThreadLevel = MPI_THREAD_SINGLE;
#ifdef _OPENMP
    // Disable OpenMP if we've not explicitly asked for threads
    if (getenv("OMP_NUM_THREADS") == NULL)
    {
        omp_set_num_threads(1);
    }
    else
    {
        // Enable threads to MPI (potentially)
        gMPIThreadLevel = MPI_THREAD_MULTIPLE;
    }
#endif

    MPI_Init_thread(&argc, &argv, gMPIThreadLevel, &gMPIThreadLevel);

    // Switch error handlers for debugging, from default (abort) to return
    //MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

#ifdef USE_MPE
    MPE_Log_get_state_eventIDs(&gLogStateBeginID, &gLogStateEndID);
    MPE_Describe_state(gLogStateBeginID, gLogStateEndID, "Temperature Update",
                       "orange");
#endif


    int flags_ok = true;
    while ((flag = getopt(argc, argv, "ic:n:r:h")) != -1)
    {
        switch (flag)
        {
            case 'i':
                g_interactive_mode = 1;
                break;
            case 'c':
                ctrl_fname = optarg;
                break;
            case 'r':
                restart = atoi(optarg);
                break;
            case 'h':
                break;
            default:
                flags_ok = false;
                break;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
    fprintf(stderr, "Communicator: iproc=%d, nproc=%d \n", iproc, nproc);

    if (iproc == 0)
        printf("Start the program... \n");

    if (!flags_ok)
    {
        if (iproc == 0)
            print_usage(argv[0]);
        retval = 1;
        goto out;
    }

    if (ctrl_fname == NULL)
    {
        if (iproc == 0)
        {
            fprintf(stderr, "No Control or Restart file specified.\n");
            print_usage(argv[0]);
        }
        retval = 1;
        goto out;
    }

    setup_mpi_datatypes();

    if (iproc == 0)
    {
        printf("Start to run Microstructure Solidification Code $\n");
#ifdef _OPENMP
        dprintf("Running on %d MPI Ranks, %d threads / rank\n", nproc,
                omp_get_max_threads());
#else
        dprintf("Running on %d MPI Ranks\n", nproc);
#endif
    }

    {
        char name[256];
        int nl = 255;
        MPI_Get_processor_name(name, &nl);
        dprintf("Rank %d on Host %s\n", iproc, name);
    }

    profiler_init();

    set_defaults();
    if (iproc == 0)
    {
        read_config(ctrl_fname);
        print_config(stdout);
    }
    /* Send configuration parameters to the tasks */
    MPI_Bcast(bp, sizeof(BB_struct), MPI_BYTE, 0, MPI_COMM_WORLD);

    int dims[3] = { bp->gnsbz, bp->gnsby, bp->gnsbx };
    int periods[3] = { 1, 1, 1 };
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &mpi_comm_new);

    MPI_Barrier(mpi_comm_new);
    MPI_Comm_size(mpi_comm_new, &nproc);
    MPI_Comm_rank(mpi_comm_new, &iproc);

    //MPI_Cart_coords not used for now (not compatible with other parts of code)
    //int coords[3];
    //MPI_Cart_coords(mpi_comm_new, iproc, nproc, coords);
    //printf("My MPI rank is %d, MPI coords are %d,%d,%d \n",
    //       iproc, coords[2], coords[1], coords[0]);

    // Set our random seed
    srand(iproc + 1 + bp->base_random_seed);

    // Setup window for Transferring Grain information
    grainSetup();

    prepareIO();

    if (iproc == 0)
    {
        gsubblock_list = ll_init(NULL, NULL);
        init_gmsp(restart);
    }

    if (restart)
    {
        if (iproc == 0)
            mainRestart(restart);
        taskRestart(restart);
    }
    else
        findStartTime();        // Now, as a group, find start time

    if (iproc == 0)
        printf
            ("Liquidus temperature will first be reached at %lg seconds\n",
             bp->timestep * bp->ts_delt);

    fflush(stdout);

    if (iproc == 0 && !restart)
        init_nucleation_sites();

    MPI_Barrier(MPI_COMM_WORLD);
    if (iproc == 0)
    {
        computeSubblocks();

        sendSubblockInfo();
    }
    recvSubblockInfo();

    if (iproc == 0)
        completeSendSubblockInfo();

    MPI_Barrier(MPI_COMM_WORLD);
    //printf("proc: %d -> %d,%d,%d\n",iproc,lsp->coords,lsp->coords.y, lsp->coords.z);
    //lsp->coords.x = coords[2];
    //lsp->coords.y = coords[1];
    //lsp->coords.z = coords[0];

    // Figure out who my neighbors are for halos exchange
    int neighbors[NUM_NEIGHBORS];
    determine_3dneighbors(iproc, neighbors);
    for (int i = 0; i < NUM_NEIGHBORS; i++)
    {
        lsp->neighbors[i][0] = neighbors[i];
    }

    allocateFields();

    profile(SETUP);
    timing(INITIALIZATION, timer_elapsed());

    struct timeval tv;
    gettimeofday(&tv, NULL);
    double gStartTime = tv.tv_sec + ((double) tv.tv_usec * 1.0e-6);

    loop(restart, gStartTime);

    if (iproc == 0)
    {
        xfree(gmsp);
    }

    timer_stop();

    gettimeofday(&tv, NULL);
    double t_end = tv.tv_sec + ((double) tv.tv_usec * 1.0e-6);
    if (iproc == 0)
        printf("Elapsed time:  %lu iterations, %lg sec\n",
               bp->timestep, t_end - gStartTime);

    if (iproc == 0)
        ll_destroy(gsubblock_list);

    profiler_collate();

  out:
    MPI_Finalize();

    xfree(bp);

    return retval;;
}
