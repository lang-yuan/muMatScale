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

// arbitrary, non-uniform periodic function
double f(double x, double y, double z)
{
    double fx = 0.2*sin(2.*M_PI*x/(bp->gdimx*bp->cellSize));
    double fy = 0.3*sin(2.*M_PI*y/(bp->gdimy*bp->cellSize));
    double fz = 0.4*sin(2.*M_PI*z/(bp->gdimz*bp->cellSize));

    return fx + fy +fz;
}

// check range of indices for function matching
int
check_values(
    double *field,
    const int dimx, const int dimy, const int dimz,
    uint32_t sbx, uint32_t sby, uint32_t sbz,
    int i0, int i1, int j0, int j1, int k0, int k1)
{
    int ret = 0;
    // check field with halo values
    int count = 0;
    double tol = 1.e-8;

    int dim3 = (dimx + 2) * (dimy + 2) * (dimz + 2);
#ifdef GPU_PACK
#pragma omp target update from(field[:dim3])
#endif

    // loop over single subarray specified by i0, i1, j0, j1, k0, k1
    for (int k = k0; k <= k1; k++)
    {
        for (int j = j0; j <= j1; j++)
        {
            for (int i = i0; i <= i1; i++)
            {
                int idx =
                    k * (dimy + 2) * (dimx + 2) +
                    j * (dimx + 2) + i;
                count++;

                double rx, ry, rz;
                scoord2realcoord(sbx, sby, sbz, i-1, j-1, k-1, &rx, &ry, &rz);
                double value = f(rx,ry,rz);

                if (fabs(field[idx] - value) > tol)
                {
                    printf("%d %d %d %le %le %le\n", i, j, k, rx, ry, rz);
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

// initialize field (no halo) as function of x,y,z
void
init_field(
    double *field,
    uint32_t sbx, uint32_t sby, uint32_t sbz,
    const int dimx, const int dimy, const int dimz)
{
    // loop over "interior" values
    for (int k = 1; k <= dimz; k++)
    {
        for (int j = 1; j <= dimy; j++)
        {
            for (int i = 1; i <= dimx; i++)
            {
                double rx, ry, rz;
                scoord2realcoord(sbx, sby, sbz, i-1, j-1, k-1, &rx, &ry, &rz);
                double value = f(rx,ry,rz);
                int idx =
                    k * (dimy + 2) * (dimx + 2) +
                    j * (dimx + 2) + i;
                field[idx] = value;
            }
        }
    }

    int dim3 = (dimx + 2) * (dimy + 2) * (dimz + 2);
#ifdef GPU_PACK
#pragma omp target update to(field[:dim3])
#endif
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

    int dims[3] = { bp->gnsbz, bp->gnsby, bp->gnsbx };
    int periods[3] = { 1, 1, 1 };
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &mpi_comm_new);

    MPI_Barrier(mpi_comm_new);
    MPI_Comm_size(mpi_comm_new, &nproc);
    MPI_Comm_rank(mpi_comm_new, &iproc);

    // get coordinates of local subblock in MPI grid
    int sb_coords[3];
    MPI_Cart_coords(mpi_comm_new, iproc, 3, sb_coords);
    // swap x and z value to be consistent with CA code
    int tmp = sb_coords[0];
    sb_coords[0] = sb_coords[2];
    sb_coords[2] = tmp;

    profiler_init();

    MPI_Barrier(mpi_comm_new);

    double *field = malloc(dim3 * sizeof(double));
    int dimxy = (bp->gsdimx + 2) * (bp->gsdimy + 2);
    int dimxz = (bp->gsdimx + 2) * (bp->gsdimz + 2);
    int dimyz = (bp->gsdimy + 2) * (bp->gsdimz + 2);

    // create buffers large enough for all faces
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

    // local block size (without halo)
    int dimx = bp->gsdimx;
    int dimy = bp->gsdimy;
    int dimz = bp->gsdimz;

    // initialize values of field in block (without halos)
    init_field(field, sb_coords[0], sb_coords[1], sb_coords[2], dimx, dimy, dimz);

    // halo fill one face at a time
    {
        if (iproc == 0)
            printf("Check FACE_TOP -> FACE_BOTTOM\n");
        // send data from face cells
        int face = FACE_TOP;
        int halo = FACE_BOTTOM;

        exchange_data(face, halo, field, send_buffer, recv_buffer, mpi_comm_new);

    }

    {
        if (iproc == 0)
            printf("Check FACE_BOTTOM -> FACE_TOP\n");
        // send data from face cells
        int face = FACE_BOTTOM;
        int halo = FACE_TOP;

        exchange_data(face, halo, field, send_buffer, recv_buffer, mpi_comm_new);
    }

    {
        if (iproc == 0)
            printf("Check FACE_LEFT -> FACE_RIGHT\n");
        // send data from face cells
        int face = FACE_LEFT;
        int halo = FACE_RIGHT;

        exchange_data(face, halo, field, send_buffer, recv_buffer, mpi_comm_new);
    }

    {
        if (iproc == 0)
            printf("Check FACE_RIGHT -> FACE_LEFT\n");
        // send data from face cells
        int face = FACE_RIGHT;
        int halo = FACE_LEFT;

        exchange_data(face, halo, field, send_buffer, recv_buffer, mpi_comm_new);
    }

    {
        if (iproc == 0)
            printf("Check FACE_FRONT -> FACE_BACK\n");
        // send data from face cells
        int face = FACE_FRONT;
        int halo = FACE_BACK;

        exchange_data(face, halo, field, send_buffer, recv_buffer, mpi_comm_new);
    }

    {
        if (iproc == 0)
            printf("Check FACE_BACK -> FACE_FRONT\n");
        // send data from face cells
        int face = FACE_BACK;
        int halo = FACE_FRONT;

        exchange_data(face, halo, field, send_buffer, recv_buffer, mpi_comm_new);
    }

    // check values in all 6 halos
    ret = check_values(field, dimx, dimy, dimz,
                       sb_coords[0], sb_coords[1], sb_coords[2],
                       1, dimx, 1, dimy, 0, 0);
    ret += check_values(field, dimx, dimy, dimz,
                        sb_coords[0], sb_coords[1], sb_coords[2],
                        1, dimx, 1, dimy, dimz + 1, dimz + 1);
    ret += check_values(field, dimx, dimy, dimz,
                        sb_coords[0], sb_coords[1], sb_coords[2],
                        dimx + 1, dimx + 1, 1, dimy, 1, dimz);
    ret += check_values(field, dimx, dimy, dimz,
                        sb_coords[0], sb_coords[1], sb_coords[2],
                        0, 0, 1, dimy, 1, dimz);
    ret += check_values(field, dimx, dimy, dimz,
                        sb_coords[0], sb_coords[1], sb_coords[2],
                        1, dimx, dimy + 1, dimy + 1, 1, dimz);
    ret += check_values(field, dimx, dimy, dimz,
                        sb_coords[0], sb_coords[1], sb_coords[2],
                        1, dimx, 0, 0, 1, dimz);

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
