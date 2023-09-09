#include "globals.h"
#include "xmalloc.h"
#include "read_ctrl.h"
#include "packing.h"
#include "distribute.h"
#include "functions.h"

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
#ifdef GPU_PACK
#pragma omp declare target
#endif
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
#ifdef GPU_PACK
#pragma omp end declare target
#endif

void
exchange_data(
    int face,
    int halo,
    double *field,
    double *send_buffer,
    double *recv_buffer)
{
    int dest = neighbors[face];
    assert(dest < nproc);
    assert(dest >= 0);
    int src = neighbors[halo];
    assert(src < nproc);
    assert(src >= 0);

    int offset;
    int stride;
    int bsize;
    int nblocks;
    computeFaceInfo(face, &offset, &stride, &bsize, &nblocks);

    printf("pack_field...\n");
#ifdef GPU_PACK
#pragma omp target
#endif
{
    pack_field(sizeof(double), field, stride, bsize, nblocks, offset,
               send_buffer);
}

#ifdef GPU_PACK
#pragma omp target update from(send_buffer[0:nblocks*bsize])
#endif

    MPI_Request req;
    MPI_Irecv(recv_buffer, dim2, MPI_DOUBLE, src, 0, MPI_COMM_WORLD, &req);
    MPI_Send(send_buffer, dim2, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
    MPI_Status mpi_status;
    MPI_Wait(&req, &mpi_status);

    // recv data in halo cells
    computeHaloInfo(halo, &offset, &stride, &bsize, &nblocks);

    printf("unpack_field...\n");
#ifdef GPU_PACK
#pragma omp target update to(recv_buffer[0:nblocks*bsize])
#pragma omp target
#endif
{
    unpack_field(sizeof(double), field, stride, bsize, nblocks, offset,
                 recv_buffer);
}

}

int
main(
    int argc,
    char *argv[])
{
    char *ctrl_fname = argv[1];

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);

    xmalloc(bp, BB_struct, 1);

    set_defaults();
    if (iproc == 0)
    {
        read_config(ctrl_fname);
        print_config(stdout);
    }

    MPI_Bcast(bp, sizeof(BB_struct), MPI_BYTE, 0, MPI_COMM_WORLD);
    int dim3 = (bp->gsdimx + 2) * (bp->gsdimy + 2) * (bp->gsdimz + 2);
    if (iproc == 0)
        printf("Local array size = %d\n", dim3);

    double *field = malloc(dim3 * sizeof(double));
    int dimxy = (bp->gsdimx + 2) * (bp->gsdimy + 2);
    int dimxz = (bp->gsdimx + 2) * (bp->gsdimz + 2);
    int dimyz = (bp->gsdimy + 2) * (bp->gsdimz + 2);

    // create buffers large enough
    dim2 = (dimxy > dimxz) ? dimxy : dimxz;
    dim2 = (dim2 > dimyz) ? dim2 : dimyz;
    double *send_buffer = malloc(dim2 * sizeof(double));
    double *recv_buffer = malloc(dim2 * sizeof(double));

#ifdef GPU_PACK
#pragma omp target enter data map(alloc:field[:dim3])
#pragma omp target enter data map(alloc:send_buffer[:dim2])
#pragma omp target enter data map(alloc:recv_buffer[:dim2])
#endif

    determine_3dneighbors(iproc, neighbors);

    // initialize field (no halo)
    double value = 3.33;

    for (int i = 0; i < dim2; i++)
        send_buffer[i] = 1.11;
    for (int i = 0; i < dim2; i++)
        recv_buffer[i] = 2.22;

    int ret = 0;

    int dimx = bp->gsdimx;
    int dimy = bp->gsdimy;
    int dimz = bp->gsdimz;

    // check halo fill one direction at a time
    {
        if (iproc == 0)
           printf("init_field...\n");
#ifdef GPU_PACK
#pragma omp target
#endif
{
        init_field(field, value, dimx, dimy, dimz);
}
        if (iproc == 0)
            printf("Check FACE_TOP -> FACE_BOTTOM\n");
        // send data from face cells
        int face = FACE_TOP;
        int halo = FACE_BOTTOM;

        exchange_data(face, halo, field, send_buffer, recv_buffer);

        // check field with halo values
#ifdef GPU_PACK
#pragma omp target update from(field[:dim3])
#endif

        if (iproc == 0)
           printf("check_values()...\n");
        ret = check_values(field, value, 1, bp->gsdimx, 1, bp->gsdimy, 0, 0);
    }

    value += 1.;
    {
#ifdef GPU_PACK
#pragma omp target
#endif
{
        init_field(field, value, dimx, dimy, dimz);
}
        if (iproc == 0)
            printf("Check FACE_BOTTOM -> FACE_TOP\n");
        // send data from face cells
        int face = FACE_BOTTOM;
        int halo = FACE_TOP;

        exchange_data(face, halo, field, send_buffer, recv_buffer);

        // check field with halo values
#ifdef GPU_PACK
#pragma omp target update from(field[:dim3])
#endif
        ret +=
            check_values(field, value, 1, bp->gsdimx, 1, bp->gsdimy,
                         bp->gsdimz + 1, bp->gsdimz + 1);
    }

    value += 1.;
    {
#ifdef GPU_PACK
#pragma omp target
#endif
{
        init_field(field, value, dimx, dimy, dimz);
}
        if (iproc == 0)
            printf("Check FACE_LEFT -> FACE_RIGHT\n");
        // send data from face cells
        int face = FACE_LEFT;
        int halo = FACE_RIGHT;

        exchange_data(face, halo, field, send_buffer, recv_buffer);

        // check field with halo values
#ifdef GPU_PACK
#pragma omp target update from(field[:dim3])
#endif
        ret +=
            check_values(field, value, bp->gsdimx + 1, bp->gsdimx + 1, 1,
                         bp->gsdimy, 1, bp->gsdimz);
    }

    value += 1.;
    {
#ifdef GPU_PACK
#pragma omp target
#endif
{
        init_field(field, value, dimx, dimy, dimz);
}
        if (iproc == 0)
            printf("Check FACE_RIGHT -> FACE_LEFT\n");
        // send data from face cells
        int face = FACE_RIGHT;
        int halo = FACE_LEFT;

        exchange_data(face, halo, field, send_buffer, recv_buffer);

        // check field with halo values
#ifdef GPU_PACK
#pragma omp target update from(field[:dim3])
#endif
        ret += check_values(field, value, 0, 0, 1, bp->gsdimy, 1, bp->gsdimz);
    }

    value += 1.;
    {
#ifdef GPU_PACK
#pragma omp target
#endif
{
        init_field(field, value, dimx, dimy, dimz);
}
        if (iproc == 0)
            printf("Check FACE_FRONT -> FACE_BACK\n");
        // send data from face cells
        int face = FACE_FRONT;
        int halo = FACE_BACK;

        exchange_data(face, halo, field, send_buffer, recv_buffer);

        // check field with halo values
#ifdef GPU_PACK
#pragma omp target update from(field[:dim3])
#endif
        ret =
            check_values(field, value, 1, bp->gsdimx, bp->gsdimy + 1,
                         bp->gsdimy + 1, 1, bp->gsdimz);
    }

    value += 1.;
    {
#ifdef GPU_PACK
#pragma omp target
#endif
{
        init_field(field, value, dimx, dimy, dimz);
}
        if (iproc == 0)
            printf("Check FACE_BACK -> FACE_FRONT\n");
        // send data from face cells
        int face = FACE_BACK;
        int halo = FACE_FRONT;

        exchange_data(face, halo, field, send_buffer, recv_buffer);

        // check field with halo values
#ifdef GPU_PACK
#pragma omp target update from(field[:dim3])
#endif
        ret = check_values(field, value, 1, bp->gsdimx, 0, 0, 1, bp->gsdimz);
    }

#ifdef GPU_PACK
#pragma omp target exit data map(delete:field[:dim3])
#pragma omp target exit data map(delete:send_buffer[:dim2])
#pragma omp target exit data map(delete:recv_buffer[:dim2])
#endif
    free(recv_buffer);
    free(send_buffer);
    free(field);

    MPI_Finalize();

    return ret;
}
