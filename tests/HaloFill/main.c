#include "face_util.h"
#include "globals.h"
#include "xmalloc.h"
#include "read_ctrl.h"
#include "packing.h"
#include "distribute.h"
#include "functions.h"
#include "profiler.h"

#include <math.h>

int neighbors[NUM_NEIGHBORS];
int dim2;


int
check_values(
    double *field,
    const double value,
    int i0,
    int i1,
    int j0,
    int j1,
    int k0,
    int k1)
{
    int ret = 0;
    // check field with halo values
    int count = 0;
    double tol = 1.e-8;
    int dim3 = (bp->gsdimx + 2) * (bp->gsdimy + 2) * (bp->gsdimz + 2);

#ifdef GPU_PACK
#pragma omp target update from(field[:dim3])
#endif

    for (int k = k0; k <= k1; k++)
    {
        for (int j = j0; j <= j1; j++)
        {
            for (int i = i0; i <= i1; i++)
            {
                int idx =
                    k * (bp->gsdimy + 2) * (bp->gsdimx + 2) +
                    j * (bp->gsdimx + 2) + i;
                count++;
                if (fabs(field[idx] - value) > tol)
                {
                    printf("%d %d %d\n", i, j, k);
                    printf
                        ("iproc = %d, Expected value: %le. Value found: %le\n",
                         iproc, value, field[idx]);
                    ret = 1;
                }
            }
        }
    }
    printf("Number of values tested: %d\n", count);

    return ret;
}

// initialize field (no halo)
void
init_field(
    double *field,
    const double value,
    const int dimx, const int dimy, const int dimz)
{
#ifdef GPU_PACK
#pragma omp target teams distribute parallel for
#endif
    for (int k = 1; k <= dimz; k++)
    {
        for (int j = 1; j <= dimy; j++)
        {
            for (int i = 1; i <= dimx; i++)
            {
                int idx =
                    k * (dimy + 2) * (dimx + 2) +
                    j * (dimx + 2) + i;
                field[idx] = value;
            }
        }
    }
}

void
exchange_data(
    int face,
    int halo,
    double *field,
    void *send_buffer[6],
    void *recv_buffer[6],
    MPI_Comm comm)
{
    int dest = neighbors[face];
    assert(dest < nproc);
    assert(dest >= 0);
    int src = neighbors[halo];
    assert(src < nproc);
    assert(src >= 0);

    MPI_Request req[2];

    printf("Send_Plane...\n");
    Send_Plane(field, sizeof(double), face, dest, 0, send_buffer, &req[0]);

    Recv_Plane(field, sizeof(double), halo, src, 0, recv_buffer, &req[1]);

    MPI_Status mpi_status;
    MPI_Waitall(2, req, MPI_STATUSES_IGNORE);

    // recv data in halo cells
    printf("unpack plane...\n");
    unpack_plane(field, sizeof(double), halo, recv_buffer);
}

int
main(
    int argc,
    char *argv[])
{
    char *ctrl_fname = argv[1];
    int ret = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    // set muMatScale communicator to test communicator
    mpi_comm_new = comm;
    MPI_Comm_size(comm, &nproc);
    MPI_Comm_rank(comm, &iproc);

    if (iproc == 0)
        printf("Test packing functions...\n");

    xmalloc(bp, BB_struct, 1);

    set_defaults();
    if (iproc == 0)
    {
        read_config(ctrl_fname);
        print_config(stdout);
    }

    MPI_Bcast(bp, sizeof(BB_struct), MPI_BYTE, 0, comm);
    int dim3 = (bp->gsdimx + 2) * (bp->gsdimy + 2) * (bp->gsdimz + 2);
    if (iproc == 0)
        printf("Local array size = %d\n", dim3);

    profiler_init();

    MPI_Barrier(comm);

    double *field = malloc(dim3 * sizeof(double));
    int dimxy = (bp->gsdimx + 2) * (bp->gsdimy + 2);
    int dimxz = (bp->gsdimx + 2) * (bp->gsdimz + 2);
    int dimyz = (bp->gsdimy + 2) * (bp->gsdimz + 2);

    // create buffers large enough
    dim2 = (dimxy > dimxz) ? dimxy : dimxz;
    dim2 = (dim2 > dimyz) ? dim2 : dimyz;
    void *send_buffer[NUM_NEIGHBORS];
    void *recv_buffer[NUM_NEIGHBORS];

    for (int face = 0; face < NUM_NEIGHBORS; face++)
    {
        send_buffer[face] = malloc(dim2 * sizeof(double));
        recv_buffer[face] = malloc(dim2 * sizeof(double));
#ifdef GPU_PACK
        double* dbuf;
#pragma omp target enter data map(alloc:field[:dim3])
        dbuf = (double*)send_buffer[face];
#pragma omp target enter data map(alloc:dbuf[:dim2])
        dbuf = (double*)recv_buffer[face];
#pragma omp target enter data map(alloc:dbuf[:dim2])
#endif
    }

    if (iproc == 0)
        printf("Determine neighbors...\n");

    determine_3dneighbors(iproc, neighbors);

    // initialize field (no halo)
    double value = 3.33;

    int dimx = bp->gsdimx;
    int dimy = bp->gsdimy;
    int dimz = bp->gsdimz;

    // check halo fill one direction at a time
    {
        if (iproc == 0)
           printf("init_field...\n");
        init_field(field, value, dimx, dimy, dimz);
       if (iproc == 0)
            printf("Check FACE_TOP -> FACE_BOTTOM\n");
        // send data from face cells
        int face = FACE_TOP;
        int halo = FACE_BOTTOM;

        exchange_data(face, halo, field, send_buffer, recv_buffer, comm);

        // check field with halo values
        if (iproc == 0)
           printf("check_values()...\n");
        ret = check_values(field, value, 1, bp->gsdimx, 1, bp->gsdimy, 0, 0);
    }

    value += 1.;
    {
        init_field(field, value, dimx, dimy, dimz);

        if (iproc == 0)
            printf("Check FACE_BOTTOM -> FACE_TOP\n");
        // send data from face cells
        int face = FACE_BOTTOM;
        int halo = FACE_TOP;

        exchange_data(face, halo, field, send_buffer, recv_buffer, comm);

        // check field with halo values
        ret +=
            check_values(field, value, 1, bp->gsdimx, 1, bp->gsdimy,
                         bp->gsdimz + 1, bp->gsdimz + 1);
    }

    value += 1.;
    {
        init_field(field, value, dimx, dimy, dimz);
        if (iproc == 0)
            printf("Check FACE_LEFT -> FACE_RIGHT\n");
        // send data from face cells
        int face = FACE_LEFT;
        int halo = FACE_RIGHT;

        exchange_data(face, halo, field, send_buffer, recv_buffer, comm);

        // check field with halo values
        ret +=
            check_values(field, value, bp->gsdimx + 1, bp->gsdimx + 1, 1,
                         bp->gsdimy, 1, bp->gsdimz);
    }

    value += 1.;
    {
        init_field(field, value, dimx, dimy, dimz);
        if (iproc == 0)
            printf("Check FACE_RIGHT -> FACE_LEFT\n");
        // send data from face cells
        int face = FACE_RIGHT;
        int halo = FACE_LEFT;

        exchange_data(face, halo, field, send_buffer, recv_buffer, comm);

        // check field with halo values
        ret += check_values(field, value, 0, 0, 1, bp->gsdimy, 1, bp->gsdimz);
    }

    value += 1.;
    {
        init_field(field, value, dimx, dimy, dimz);
        if (iproc == 0)
            printf("Check FACE_FRONT -> FACE_BACK\n");
        // send data from face cells
        int face = FACE_FRONT;
        int halo = FACE_BACK;

        exchange_data(face, halo, field, send_buffer, recv_buffer, comm);

        // check field with halo values
        ret =
            check_values(field, value, 1, bp->gsdimx, bp->gsdimy + 1,
                         bp->gsdimy + 1, 1, bp->gsdimz);
    }

    value += 1.;
    {
        init_field(field, value, dimx, dimy, dimz);
        if (iproc == 0)
            printf("Check FACE_BACK -> FACE_FRONT\n");
        // send data from face cells
        int face = FACE_BACK;
        int halo = FACE_FRONT;

        exchange_data(face, halo, field, send_buffer, recv_buffer, comm);

        // check field with halo values
        ret = check_values(field, value, 1, bp->gsdimx, 0, 0, 1, bp->gsdimz);
    }

    for (int face = 0; face < NUM_NEIGHBORS; face++)
    {
#ifdef GPU_PACK
        double* dbuf;
#pragma omp target exit data map(delete:field[:dim3])
        dbuf = (double*)send_buffer[face];
#pragma omp target exit data map(delete:dbuf[:dim2])
        dbuf = (double*)recv_buffer[face];
#pragma omp target exit data map(delete:dbuf[:dim2])
#endif
        free(recv_buffer[face]);
        free(send_buffer[face]);
    }
    free(field);

    MPI_Finalize();

    return ret;
}
