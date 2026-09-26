#include "face_util.h"
#include "globals.h"
#include "xmalloc.h"
#include "read_ctrl.h"
#include "packing.h"
#include "distribute.h"
#include "functions.h"
#include "profiler.h"
#include "variables.h"
#include "ll.h"

#include <math.h>

extern ll_t *gsubblock_list;
extern SB_struct *lsp;

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
    double *field, const int depth,
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
#pragma omp target update from(field[:depth*dim3])
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

                for (int d = 0; d < depth; d++)
                if (fabs(field[depth*idx+d] - value) > tol)
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
    double *field, const int depth,
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
                for (int d = 0; d < depth; d++)
                   field[depth*idx+d] = value;
            }
        }
    }

    int dim3 = (dimx + 2) * (dimy + 2) * (dimz + 2);
#ifdef GPU_PACK
#pragma omp target update to(field[:depth*dim3])
#endif
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

    profiler_init();

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

    MPI_Barrier(MPI_COMM_WORLD);
    if (iproc == 0)
    {
        gsubblock_list = ll_init(NULL, NULL);
        init_gmsp();
        sendSubblockInfo();
    }
    recvSubblockInfo();

    if (iproc == 0)
        completeSendSubblockInfo();

    MPI_Barrier(MPI_COMM_WORLD);
    lsp->coords.x = sb_coords[0];
    lsp->coords.y = sb_coords[1];
    lsp->coords.z = sb_coords[2];
    lsp->subblockid = iproc;

    profiler_init();

    MPI_Barrier(mpi_comm_new);

    int dimxy = (bp->gsdimx + 2) * (bp->gsdimy + 2);
    int dimxz = (bp->gsdimx + 2) * (bp->gsdimz + 2);
    int dimyz = (bp->gsdimy + 2) * (bp->gsdimz + 2);

    if (iproc == 0)
        printf("Determine neighbors...\n");

    int neighbors[NUM_NEIGHBORS];
    determine_3dneighbors(iproc, neighbors);
    for (int i = 0; i < NUM_NEIGHBORS; i++)
    {
        lsp->neighbors[i][0] = neighbors[i];
    }

    // local block size (without halo)
    int dimx = bp->gsdimx;
    int dimy = bp->gsdimy;
    int dimz = bp->gsdimz;

for(int depth = 1; depth < 4; depth+=2)
{
    printf("Test depth %d\n", depth);

    double *field = malloc(dim3 * depth * sizeof(double));
#ifdef GPU_PACK
#pragma omp target enter data map(alloc:field[:depth*dim3])
#endif

    // initialize values of field in block (without halos)
    if (iproc == 0)
        printf("Initialize field...\n");
    init_field(field, depth, sb_coords[0], sb_coords[1], sb_coords[2], dimx, dimy, dimz);

    if (iproc == 0)
        printf("Exchange data...\n");

    MPI_Barrier(mpi_comm_new);

    int var = registerCommInfo(depth*sizeof(double));
    ExchangeFacesForVar(var, field);

    FinishExchangeForVar(var, field);

    // check values in all 6 halos
    MPI_Barrier(MPI_COMM_WORLD);
    if (iproc == 0)
        printf("Check 1st halo...\n");
    ret += check_values(field, depth, dimx, dimy, dimz,
                       sb_coords[0], sb_coords[1], sb_coords[2],
                       1, dimx, 1, dimy, 0, 0);

    MPI_Barrier(MPI_COMM_WORLD);
    if (iproc == 0)
        printf("Check 2nd halo...\n");
    ret += check_values(field, depth, dimx, dimy, dimz,
                        sb_coords[0], sb_coords[1], sb_coords[2],
                        1, dimx, 1, dimy, dimz + 1, dimz + 1);

    MPI_Barrier(MPI_COMM_WORLD);
    if (iproc == 0)
        printf("Check 3rd halo...\n");
    ret += check_values(field, depth, dimx, dimy, dimz,
                        sb_coords[0], sb_coords[1], sb_coords[2],
                        dimx + 1, dimx + 1, 1, dimy, 1, dimz);
    MPI_Barrier(MPI_COMM_WORLD);
    if (iproc == 0)
        printf("Check 4th halo...\n");
    ret += check_values(field, depth, dimx, dimy, dimz,
                        sb_coords[0], sb_coords[1], sb_coords[2],
                        0, 0, 1, dimy, 1, dimz);

    MPI_Barrier(MPI_COMM_WORLD);
    if (iproc == 0)
        printf("Check 5th halo...\n");
    ret += check_values(field, depth, dimx, dimy, dimz,
                        sb_coords[0], sb_coords[1], sb_coords[2],
                        1, dimx, dimy + 1, dimy + 1, 1, dimz);
    MPI_Barrier(MPI_COMM_WORLD);
    if (iproc == 0)
        printf("Check 6th halo...\n");
    ret += check_values(field, depth, dimx, dimy, dimz,
                        sb_coords[0], sb_coords[1], sb_coords[2],
                        1, dimx, 0, 0, 1, dimz);

#ifdef GPU_PACK
#pragma omp target exit data map(delete:field[:dim3])
#endif
    free(field);
}

    if (iproc == 0){
        ll_destroy(gsubblock_list);
        xfree(gmsp);
    }

    MPI_Finalize();

    xfree(bp);

    return ret;
}
