/***************************************************************/
/* Copyright (c) 2023, Lang Yuan, Univeristy of South Carolina */
/* All rights reserved.                                        */
/* This file is part of muMatScale.                            */
/* See the top-level LICENSE file for details.                 */
/***************************************************************/

#include <unistd.h>
#include <stddef.h>

#include <math.h>
#include "grain.h"
#include "globals.h"
#include "functions.h"
#include "xmalloc.h"
#include "debug.h"
#include "ll.h"
#include "profiler.h"

/* Constants taken from uMatIC */
#define G_NUC_MIN_UND 1.0e-5
#define G_NUC_SIG_MULT 6.0

typedef struct nuc_volume
{
    double density;
    double minx;
    double miny;
    double minz;
    double maxx;
    double maxy;
    double maxz;
    double threshold;
    size_t minCellX;
    size_t minCellY;
    size_t minCellZ;
    size_t maxCellX;
    size_t maxCellY;
    size_t maxCellZ;
} nuc_volume_t;

/* Temporaries for startup, reading config, etc */
static ll_t *fixed_nucs_ll = NULL;
static ll_t *nuc_volumes_ll = NULL;

static MPI_Datatype grainCommType;
extern SB_struct *lsp;

grain_t *grain_cache;

#ifdef GPU_OMP
#pragma omp declare target
#endif

size_t numNewGrains = 0;
#ifndef GPU_OMP_NUC
static size_t newGrainSize = 0;
#endif
loc_t *newGrainLocs = NULL;

#ifdef GPU_OMP
#pragma omp end declare target
#endif

static nucleation_t *g_nuc_pt;  // initial Nucleation Points
//static size_t g_nnuc;
static size_t g_nnuc = 0;

struct select_if_not_fixed_ud
{
    SB_struct *sb;
    size_t *cells;
    double gaussCenter;
    double gaussSigma;
};


/**
 * Adds a new cell location to this list of potential nucleations
 * to be handled at the end of the timestep.
 *
 * \param[in] x  global cell index
 * \param[in] y  global cell index
 * \param[in] z  global cell index
 */

#ifndef GPU_OMP_NUC
static void
addNewGrain(
    uint32_t x,
    uint32_t y,
    uint32_t z)
{
    {
        if (numNewGrains + 1 > newGrainSize)
        {
            newGrainSize = MAX(8, 2 * newGrainSize);
            xrealloc(newGrainLocs, loc_t, newGrainSize);
        }
        // Can't nucleate any more.  We've maxed out!
        if ((bp->num_grains + numNewGrains) < (bp->maxTotalGrains + 1)) // +1 for Grain 0
        {
            newGrainLocs[numNewGrains].x = x;
            newGrainLocs[numNewGrains].y = y;
            newGrainLocs[numNewGrains].z = z;
            numNewGrains++;
        }
    }
}
#endif

/**
 * Sets the orientation for a nucleation's growth.
 *
 * Either oriented vertically, or random if \c GNOriented is set
 */
static void
setNucOrientation(
    nucleation_t * nuc)
{
    nuc->gx = nuc->gy = nuc->gz = 0.0;
    if (bp->gnUseOrientation)
    {
        // gx/gy/gz should be in radians
        float acon = 4.0;
        nuc->gx = M_PI / acon - 2 * (M_PI / acon * getRandScale());
        nuc->gy = M_PI / acon - 2 * (M_PI / acon * getRandScale());
        nuc->gz = M_PI / acon - 2 * (M_PI / acon * getRandScale());
    }
}


/**
 * Pict a random threshold for undercooling with gaussian distribution
 *
 * \param[in] gaussCenter  Center of the distribution
 * \param[in] gaussSigma   Standard Deviation of the distribution
 * \return a undercooling threshold
 */
static float
genRandomThreshold(
    double gaussCenter,
    double gaussSigma)
{
    float t = 0.0;
    do
    {
        t = (float) (gaussSigma * getRandomGaussian() + gaussCenter);
    }
    while ((t < G_NUC_MIN_UND) ||
           (fabs(t - gaussCenter) > (G_NUC_SIG_MULT * gaussSigma)));

    return t;
}


typedef void (
    *select_func) (
    size_t index,
    void *userdata);


/**
 * Selects a site as a nucleation point.
 *
 * \param[in] index    index into an array to become a nucleation point
 * \param[in] userdata \c void* cast \c nuc_volume that describes the
 *                     volume in which the nucleation sites are being
 *                     distributed
 *
 * \see distribute_nucleations
 */
static void
select_nuc_volume(
    size_t index,
    void *userdata)
{
    nuc_volume_t *v = (nuc_volume_t *) userdata;

    // Index into the virtual "list" of available cells
    // nx, ny, nz = num of cells in x, y, z
    size_t nx = v->maxCellX - v->minCellX;
    size_t ny = v->maxCellY - v->minCellY;

    // vx, vy, vz = volume indicies, x, y, z
    size_t vx = index % nx;
    size_t vy = (index % (nx * ny)) / (nx);
    size_t vz = index / (nx * ny);

    nucleation_t nuc;
    // cx, cy, cz = actual coordinates of site
    nuc.cx = v->minCellX + vx;
    nuc.cy = v->minCellY + vy;
    nuc.cz = v->minCellZ + vz;
    nuc.threshold = v->threshold;

    setNucOrientation(&nuc);

    nuc.x = CELLCOORD2REALCOORD0(nuc.cx);
    nuc.y = CELLCOORD2REALCOORD1(nuc.cy);
    nuc.z = CELLCOORD2REALCOORD2(nuc.cz);

    /* We could do a uniqueness test, but that would be  O(n^2)
     * Instead, just set it, and later, when we're looking for nuc
     * points, we'll find the earlier one first anyway.
     */
    g_nuc_pt[g_nnuc++] = nuc;


}


/**
 * Sets a random nucleation threshold at a cell if one doesn't
 * already exist.
 *
 * Designed to be used with \c distribute_nucleation
 *
 * \param[in] index    index into array of the site
 * \param[in] userdata \c void* casted \c select_if_not_fixed_ud which
 *                     contains the subblock and the cell index lookup map
 *
 * \see distribute_nucleations
 */
static void
select_if_not_fixed(
    size_t index,
    void *userdata)
{
    // If nobody else is already here (fixed grain, etc), add threshold
    struct select_if_not_fixed_ud *data =
        (struct select_if_not_fixed_ud *) userdata;
    size_t myindex = data->cells[index];
    if (data->sb->nuc_threshold[myindex] == INFINITY)
        data->sb->nuc_threshold[myindex] =
            genRandomThreshold(data->gaussCenter, data->gaussSigma);
}


/**
 * Selects random cells from a list of cells to be used for nucleation
 *
 * \param[in] nc                 number of cells in array
 * \param[in] num_to_select      number of cells to pick
 * \param[in] selection_function Function Pointer to be called on
 *                               selected cells
 * \param[in] userdata           caller-specified data to be passed
 *                               to the selection function
 */
static void
distribute_nucleations(
    size_t nc,
    size_t num_to_select,
    select_func selection_function,
    void *userdata)
{
    dprintf
        ("Distributing %lu Nucleation points over %lu cells with func %p\n",
         num_to_select, nc, selection_function);


    /* See Knuth Volume 2, Algorithm S (pg 122)
     * Random Sampling
     */
    size_t n = MIN(num_to_select, nc);  // Guard against number too high
    size_t N = nc;
    size_t t = 0, m = 0;
    size_t idx = 0;

    // short-circuit
    if (n == 0)
        return;

    do
    {
        double U = getRandScale();
        if ((N - t) * U >= (n - m))
        {
            t++;
        }
        else
        {
            // Select
            selection_function(idx, userdata);
            t++;
            m++;
        }
        idx++;
    }
    while (m < n);
    if (!(idx <= nc))
    {
        error("%lu %lu %lu %lu %lu\n", t, m, n, N, idx);
    }

}


/* Task */
static void
nucleate_grain(
    nucleation_t * nuc,
    grain_t * grain)
{

    double c0, s0, c1, s1, c2, s2;
    double a0, a1, a2;

    // gx/gy/gz should be in radians
    a0 = nuc->gx;
    a1 = nuc->gy;
    a2 = nuc->gz;

    c0 = cos(a0);
    s0 = sin(a0);
    c1 = cos(a1);
    s1 = sin(a1);
    c2 = cos(a2);
    s2 = sin(a2);

    /* rotation matrix */
    grain->rotang[0] = nuc->gx;
    grain->rotang[1] = nuc->gy;
    grain->rotang[2] = nuc->gz;

    grain->rotmat[0][0] = c2 * c0 - s2 * s1 * s0;
    grain->rotmat[1][0] = c2 * s0 + s2 * c1 * c0;
    grain->rotmat[2][0] = s2 * s1;

    grain->rotmat[0][1] = -1. * s2 * c0 - c2 * c1 * s0;
    grain->rotmat[1][1] = -1. * s2 * s0 + c2 * c1 * c0;
    grain->rotmat[2][1] = c2 * s1;

    grain->rotmat[0][2] = s1 * s0;
    grain->rotmat[1][2] = -1. * s1 * c0;
    grain->rotmat[2][2] = c1;

    grain->nuc_x = nuc->x;
    grain->nuc_y = nuc->y;
    grain->nuc_z = nuc->z;
    grain->nuc_timestep = bp->timestep;
    grain->nuc_threshold = nuc->threshold;

}


static void
createGrains(
    size_t count,
    size_t base_gr_num)
{
    for (size_t g = 0; g < count; g++)
    {
        int gr_num = base_gr_num + g;
        /* Set the grain number in the array */
        if (gr_num > bp->maxTotalGrains)
        {
            error("This shouldn't happed due to previous conditional!\n");
        }

        int z = newGrainLocs[g].z;
        int x = newGrainLocs[g].x;
        int y = newGrainLocs[g].y;


        /* Step 1 */

        // Adjust to subblock-internal coordinates
        int sx = 1 + (x % bp->gsdimx);
        int sy = 1 + (y % bp->gsdimy);
        int sz = 1 + (z % bp->gsdimz);

        int idx =
            sx + sy * (bp->gsdimx + 2) + sz * (bp->gsdimx + 2) * (bp->gsdimy +
                                                                  2);
        SB_struct *sb = lsp;
        sb->gr[idx] = gr_num;

        /* Step 2 */
        int found = 0;
        nucleation_t nuc;

        // if the nucleation is a 'fixed nuc', get its data
        // otherwise, we do some randomization
        for (int j = 0; j < sb->nnuc; j++)
        {
            if ((sb->nuc_pt[j].cx == x) &&
                (sb->nuc_pt[j].cy == y) && (sb->nuc_pt[j].cz == z))
            {
                nuc = sb->nuc_pt[j];
                found = 1;
                break;
            }
        }
        if (!found)
        {
            // Generate a new nucleation, including orientation
            nuc.cx = x;
            nuc.cy = y;
            nuc.cz = z;
            nuc.threshold = sb->nuc_threshold[idx];
            nuc.x = CELLCOORD2REALCOORD0(x);
            nuc.y = CELLCOORD2REALCOORD1(y);
            nuc.z = CELLCOORD2REALCOORD2(z);

            setNucOrientation(&nuc);
        }

        nucleate_grain(&nuc, &grain_cache[gr_num]);
        dprintf("putting a grain into offset %d\n", gr_num);
    }
}



/***** PUBLIC METHODS *****/

/**
 * Initialize the grain module.
 * Called by both Main & Task during initialization
 */
void
grainSetup(
    void)
{
    // Build grainCommType
    const int count = 3;
    int blocklens[] = { 3, 1, 13 };
    MPI_Aint indicies[] = {
        0,
        offsetof(grain_t, nuc_timestep),
        offsetof(grain_t, nuc_threshold)
    };
    MPI_Datatype componentTypes[] = { MPI_DOUBLE,
        MPI_UNSIGNED_LONG_LONG,
        MPI_DOUBLE
    };

    //MPI_Type_struct(count, blocklens, indicies, componentTypes, &grainCommType);
    MPI_Type_create_struct(count, blocklens, indicies, componentTypes,
                           &grainCommType);
    MPI_Type_commit(&grainCommType);


    // pre-calculate the potential grain numbers... ly

    if (bp->doSurfaceNuc && (bp->maxGrainDensitySurf > 0.0))
    {
        bp->pre_num_grains +=
            ceil(bp->maxGrainDensitySurf *
                 (bp->gdimx * bp->gdimy + bp->gdimx * bp->gdimz +
                  bp->gdimz * bp->gdimy) * bp->cellSize * bp->cellSize * 2.0);
    }

    if (bp->temp_type == INTERNAL)
    {
        int nuc_nuc = 0;
        nuc_nuc = ceil(bp->maxGrainDensity * bp->gdimx * bp->gdimy * bp->gdimz
                       * bp->cellSize * bp->cellSize * bp->cellSize);
        bp->pre_num_grains += nuc_nuc;
    }

    // minimize the memory allocation for bp->grain_cache
    //bp->pre_num_grains = (bp->pre_num_grains > bp->maxTotalGrains)? bp->maxTotalGrains : bp->pre_num_grains ;

    xmalloc(grain_cache, grain_t, bp->maxTotalGrains + 1);
#pragma omp target enter data map(to:grain_cache[0:bp->maxTotalGrains + 1])
    //xmalloc(bp->grain_cache, grain_t, bp->pre_num_grains +1);
}


/**
 * Deinitialize the grain module.
 * Called by both Main & Task during shutdown
 */
void
grainShutdown(
    void)
{
    MPI_Type_free(&grainCommType);
    xfree(grain_cache);
}


/**
 * Takes a coordinate and orientation angle to specify a nucleation point.
 * Adds this information to a list to later be created as a nucleation point.
 * Called by the main during configuration reading.
 *
 * \see init_nucleation_sites
 *
 * \param[in] x  global coordinate
 * \param[in] y  global coordinate
 * \param[in] z  global coordinate
 * \param[in] ax Rotation Angle
 * \param[in] ay Rotation Angle
 * \param[in] az Rotation Angle
 */
void
add_fixed_nuc(
    double x,
    double y,
    double z,
    double ax,
    double ay,
    double az,
    double tsh)
{
    if (fixed_nucs_ll == NULL)
        fixed_nucs_ll = ll_init(NULL, NULL);

    nucleation_t *n;
    xmalloc(n, nucleation_t, 1);
    n->x = x;
    n->y = y;
    n->z = z;
    // gx/gy/gz should be in radians, ax/ay/xz are in degrees
    n->gx = (M_PI / 180.0) * ax;
    n->gy = (M_PI / 180.0) * ay;
    n->gz = (M_PI / 180.0) * az;
    n->threshold = tsh;
    ll_insert(fixed_nucs_ll, n);
    fprintf(stderr,
            "Fixed grains: location (%g, %g, %g) and angle(%g, %g, %g)\n", x,
            y, z, ax, ay, az);

}


/**
 * Specifies a volume in which to seed nucleation points.
 * Stores this information for later processing.
 * \see init_nucleation_sites
 * \param[in] density  The density (as a percentage [0 to 1]) of nucleation sites in the volume
 * \param[in] minx     Minimum X coordinate for the volume
 * \param[in] miny     Minimum Y coordinate for the volume
 * \param[in] minz     Minimum Z coordinate for the volume
 * \param[in] maxx     Maximum X coordinate for the volume
 * \param[in] maxy     Maximum Y coordinate for the volume
 * \param[in] maxz     Maximum Z coordinate for the volume
 */
void
add_nuc_volume(
    double density,
    double minx,
    double miny,
    double minz,
    double maxx,
    double maxy,
    double maxz,
    double tsh)
{
    if (nuc_volumes_ll == NULL)
        nuc_volumes_ll = ll_init(NULL, NULL);

    nuc_volume_t *v;
    xmalloc(v, nuc_volume_t, 1);
    v->density = density;
    v->minx = minx;
    v->miny = miny;
    v->minz = minz;
    v->maxx = maxx;
    v->maxy = maxy;
    v->maxz = maxz;
    v->threshold = tsh;
    ll_insert(nuc_volumes_ll, v);
}


/**
 * Initializes all the nucleation sites.
 *
 * There are 3 ways to specify nucleation sites:
 *  \li FixedNuc: Fixed Nucleation Sites
 *  \li NucVolume: Random fill inside specific volume
 *  \li "Purely" Random 'GNGaussian' in the body or on the surface
 *
 * Each potential nucleation site will have an undercooling threshold, that
 * when reached, will cause a grain to form at that site (unless one has
 * already grown into the cell).
 *
 * This function takes all the specified sites and generates their info.
 *
 * \see add_fixed_nuc
 * \see add_nuc_volume
 */
void
init_nucleation_sites(
    void)
{
    /*
     * 3 Types of Nucleations:
     *   * FixedNuc: Fixed Nucleation Sites
     *   * NucVolume: Random fill inside specific volume
     *   * "Purely" Random 'GNGaussian' from uMatIC
     *
     * We'll init the Fixed Nuc's first, then NucVolume, then Gaussian
     */

    /* Step 1:  Fixed Nucleation Points! */
    if (fixed_nucs_ll != NULL)
    {
        size_t n_fnuc = ll_count(fixed_nucs_ll);
        printf("Initializing %lu fixed nucleation points\n", n_fnuc);
        xrealloc(g_nuc_pt, nucleation_t, n_fnuc);
        lli_t *llp = fixed_nucs_ll->head;
        while (llp != NULL)
        {
            nucleation_t *n = (nucleation_t *) llp->data;
            realcoord2scoord(n->x, n->y, n->z, &n->cx, &n->cy, &n->cz);
            llp = llp->next;

            // Copy over the data;
            g_nuc_pt[g_nnuc++] = *n;

            xfree(n);
        }
        ll_destroy(fixed_nucs_ll);
        fixed_nucs_ll = NULL;
    }

    /* Step 2:  NucVolumes */
    if (nuc_volumes_ll != NULL)
    {
        size_t n_fnuc = ll_count(nuc_volumes_ll);
        printf("Initializing %lu random volumes of nucleation points\n",
               n_fnuc);
        lli_t *llp = nuc_volumes_ll->head;
        while (llp != NULL)
        {
            nuc_volume_t *v = (nuc_volume_t *) llp->data;
            v->minx = MAX(v->minx, bp->origin_offset[0]);
            v->maxx =
                MIN(v->maxx,
                    bp->origin_offset[0] + (bp->gdimx * bp->cellSize));
            v->miny = MAX(v->miny, bp->origin_offset[1]);
            v->maxy =
                MIN(v->maxy,
                    bp->origin_offset[1] + (bp->gdimy * bp->cellSize));
            v->minz = MAX(v->minz, bp->origin_offset[2]);
            v->maxz =
                MIN(v->maxz,
                    bp->origin_offset[2] + (bp->gdimz * bp->cellSize));
            realcoord2scoord(v->minx, v->miny, v->minz, &v->minCellX,
                             &v->minCellY, &v->minCellZ);
            realcoord2scoord(v->maxx, v->maxy, v->maxz, &v->maxCellX,
                             &v->maxCellY, &v->maxCellZ);

            printf
                ("Volume is from [%g, %g, %g] to [%g, %g, %g] with %g%% density\n",
                 v->minx, v->miny, v->minz, v->maxx, v->maxy, v->maxz,
                 v->density * 100.0);
            llp = llp->next;
            double volume_size =
                (v->maxx - v->minx) * (v->maxy - v->miny) * (v->maxz -
                                                             v->minz);
            size_t volume_count =
                (v->maxCellX - v->minCellX) * (v->maxCellY -
                                               v->minCellY) * (v->maxCellZ -
                                                               v->minCellZ);
            size_t max_nucs = volume_count * v->density;
            printf("Volume size: %lg, count: %lu, max nucs: %lu\n",
                   volume_size, volume_count, max_nucs);
            xrealloc(g_nuc_pt, nucleation_t, g_nnuc + max_nucs);

            distribute_nucleations(volume_count, max_nucs, select_nuc_volume,
                                   v);

            xfree(v);
        }
        ll_destroy(nuc_volumes_ll);
        nuc_volumes_ll = NULL;
    }


    for (size_t i = 0; i < g_nnuc; i++)
    {
        double x = g_nuc_pt[i].x;
        double y = g_nuc_pt[i].y;
        double z = g_nuc_pt[i].z;
        int sbid = sbID_from_realcoords(x, y, z);
        MSB_struct *sb = &gmsp[sbid];
        if (sb->nnuc + 1 > sb->nuc_pt_len)
        {
            // alloc new space
            size_t newlen =
                MAX(16, MIN(sb->nuc_pt_len * 2, sb->nuc_pt_len + 256));
            dprintf("reallocing %d->nuc_pt from %zu to %zu\n", sbid,
                    sb->nuc_pt_len, newlen);
            xrealloc(sb->nuc_pt, nucleation_t, newlen);
            sb->nuc_pt_len = newlen;
        }
        sb->nuc_pt[sb->nnuc] = g_nuc_pt[i];
        sb->nnuc++;
    }

    g_nnuc = 0;
}


static void
calculateGrainSizes(
    SB_struct * sb,
    size_t * sizes)
{
    int dimx = bp->gsdimx;
    int dimy = bp->gsdimy;
    int dimz = bp->gsdimz;
    for (int z = 1; z <= dimz; z++)
    {
        for (int y = 1; y <= dimy; y++)
        {
            for (int x = 1; x <= dimx; x++)
            {
                if (sb->gr[SBIDX(x, y, z)] > 0)
                {
                    sizes[sb->gr[SBIDX(x, y, z)]]++;
                }
            }
        }
    }
}


/**
 * Creates a CSV file with information about every grain.
 * Called by ALL
 */
void
output_grains(
    SB_struct * sb)
{
#ifndef PATH_MAX
#define PATH_MAX 512
#endif
    size_t *grainsizes;
    xmalloc(grainsizes, size_t, bp->num_grains);
    // First, collate grain sizes
    calculateGrainSizes(sb, grainsizes);

    void *sendbuf, *recvbuf;
    if (iproc == 0)
    {
        sendbuf = MPI_IN_PLACE;
        recvbuf = grainsizes;
    }
    else
    {
        sendbuf = grainsizes;
        recvbuf = NULL;
    }
    timing(COMPUTATION, timer_elapsed());

    MPI_Reduce(sendbuf, recvbuf, bp->num_grains, MPI_UNSIGNED_LONG, MPI_SUM,
               0, mpi_comm_new);
    timing(COMMUNICATION, timer_elapsed());

    // Now, output data
    if (iproc == 0)
    {
        char fname[PATH_MAX] = { 0 };
        snprintf(fname, PATH_MAX - 1, "%s_grains.csv", bp->basefilename);
        FILE *fp = fopen(fname, "w");
        if (!fp)
        {
            alert("Unable to open grain file %s for writing!\n", fname);
            return;
        }

        fprintf(fp,
                "Number,Nuc Timestep,Nuc Time,Nuc Threshold,X,Y,Z,CellX,CellY,CellZ,Rot1,Rot2,Rot3,Rotation Matrix,Size mm^3\n");
        size_t minSize = grainsizes[1];
        size_t maxSize = grainsizes[1];
        size_t avgSize = 0;
        size_t zeroCount = 0;
        double cellVolume = pow(bp->cellSize * 1000, 3.0);      // Volume in mm^3
        for (int i = 1; i < bp->num_grains; i++)
        {
            size_t cx, cy, cz;
            grain_t *g = &grain_cache[i];
            realcoord2scoord(g->nuc_x, g->nuc_y, g->nuc_z, &cx, &cy, &cz);
            fprintf(fp,
                    "%d,%lu,%lg,%lg,%lg,%lg,%lg,%lu,%lu,%lu,%lg,%lg,%lg,[[%lg %lg %lg][%lg %lg %lg][%lg %lg %lg]],%lg\n",
                    i, g->nuc_timestep, g->nuc_timestep * bp->ts_delt,
                    g->nuc_threshold, g->nuc_x, g->nuc_y, g->nuc_z, cx, cy,
                    cz, g->rotang[0], g->rotang[1], g->rotang[2],
                    g->rotmat[0][0], g->rotmat[0][1], g->rotmat[0][2],
                    g->rotmat[1][0], g->rotmat[1][1], g->rotmat[1][2],
                    g->rotmat[2][0], g->rotmat[2][1], g->rotmat[2][2],
                    grainsizes[i] * cellVolume);
            if (grainsizes[i] < minSize)
                minSize = grainsizes[i];
            if (grainsizes[i] > maxSize)
                maxSize = grainsizes[i];
            avgSize += grainsizes[i];
            if (grainsizes[i] == 0)
                zeroCount++;
        }
        printf
            ("%d grains.  Min: %lg mm^3   Max: %lg mm^3   Average:  %lg mm^3\n",
             bp->num_grains - 1, minSize * cellVolume, maxSize * cellVolume,
             ((double) avgSize /
              ((double) ((bp->num_grains - 1) - zeroCount)) * cellVolume));

        fclose(fp);
        timing(FILEIO, timer_elapsed());
    }
    xfree(grainsizes);
    timing(COMPUTATION, timer_elapsed());
}


/**
 * Organizes and activates any new grains, communicating between
 * all hosts  (Main & Tasks)
 */
void
activateNewGrains(
    void)
{
    static int *newGrainArray = NULL;
    static int *displacements = NULL;
    if (newGrainArray == NULL)
    {
        xmalloc(newGrainArray, int,
                nproc);
        xmalloc(displacements, int,
                nproc);
    }

    int grains = numNewGrains;
    timing(COMPUTATION, timer_elapsed());

    MPI_Allgather(&grains, 1, MPI_INT, newGrainArray, 1, MPI_INT,
                  mpi_comm_new);
    timing(COMMUNICATION, timer_elapsed());
    profile(GRAIN_ACTIV_SYNC1);

    int new_activations = 0;
    for (int i = 0; i < nproc; i++)
    {
        new_activations += newGrainArray[i];
    }
    static int warned = 0;
    if (bp->num_grains + new_activations > (bp->maxTotalGrains + 1))    // +1 to count for grain 0
    {
        if (!warned && (iproc == 0))
        {
            printf
                ("\nReaced Maximum Total Grains (%lu).  No more nucleations allowed.\n",
                 bp->maxTotalGrains);
            warned = 1;
        }
        // Cap our activations at the remaining slots
        new_activations = (bp->maxTotalGrains + 1) - bp->num_grains;
        size_t gr_sum = 0;
        for (int i = 0; i < nproc; i++)
        {
            if (gr_sum + newGrainArray[i] > bp->maxTotalGrains)
            {
                newGrainArray[i] = bp->maxTotalGrains - gr_sum;
            }
            gr_sum += newGrainArray[i];
        }
    }

    if (new_activations)
    {
        // Calculate offsets
        displacements[0] = bp->num_grains;
        for (int i = 1; i < nproc; i++)
        {
            displacements[i] = displacements[i - 1] + newGrainArray[i - 1];
        }

        // Create the grains
        createGrains(newGrainArray[iproc], displacements[iproc]);
        profile(GRAIN_ACTIVATION);
        timing(COMPUTATION, timer_elapsed());

        MPI_Allgatherv(MPI_IN_PLACE, newGrainArray[iproc], grainCommType,
                       grain_cache, newGrainArray, displacements,
                       grainCommType, mpi_comm_new);
        timing(COMMUNICATION, timer_elapsed());
        profile(GRAIN_ACTIV_SYNC2);
#pragma omp target update to(grain_cache[0:bp->maxTotalGrains + 1])
    }

    bp->num_grains += new_activations;
    numNewGrains = 0;
}



/**
 * Called when a new subblock is being activated, this function fills in
 * the potential nucleation sites that reside within this subblock's
 * boundaries.
 *
 * Will assign both fixed nuclei and gaussian-distributed surface and
 * volume nucleation points.
 *
 * \param[in] sb  The subblock
 */
void
init_sb_nucleation(
    SB_struct * sb)
{
    int dimx = bp->gsdimx;
    int dimy = bp->gsdimy;
    // Set any potential nucleation points
    for (int nidx = 0; nidx < sb->nnuc; nidx++)
    {
        nucleation_t *n = &sb->nuc_pt[nidx];

        // Convert to local indicies
        size_t x = 1 + n->cx % bp->gsdimx;
        size_t y = 1 + n->cy % bp->gsdimy;
        size_t z = 1 + n->cz % bp->gsdimz;

        size_t index = z * (dimy + 2) * (dimx + 2) + y * (dimx + 2) + x;

        if (!sb->mold[index])   // Don't set for grains in the mold
        {
            // Don't activate yet, will do so in calculation
            // depending on threshold
            sb->nuc_threshold[index] = n->threshold;
        }
    }

    /* Perform Gaussian SB Creationr
     *
     * Two types: Parameters for surface cells are different than volume cells
     */
    size_t *suf_cells;
    size_t *vol_cells;
    xmalloc(suf_cells, size_t, bp->gsdimx * bp->gsdimy * bp->gsdimz);
    xmalloc(vol_cells, size_t, bp->gsdimx * bp->gsdimy * bp->gsdimz);
    size_t s = 0, v = 0;
    for (int z = 1; z <= bp->gsdimz; z++)
        for (int y = 1; y <= bp->gsdimy; y++)
            for (int x = 1; x <= bp->gsdimx; x++)
            {
                size_t i = SBIDX(x, y, z);
                if (!sb->mold[i])
                {
                    vol_cells[v++] = i;
                    // Check to see if we're on the surface, too
                    int surf = 0;
                    surf += sb->mold[SBIDX(x - 1, y, z)];
                    surf += sb->mold[SBIDX(x + 1, y, z)];
                    surf += sb->mold[SBIDX(x, y - 1, z)];
                    surf += sb->mold[SBIDX(x, y + 1, z)];
                    surf += sb->mold[SBIDX(x, y, z - 1)];
                    surf += sb->mold[SBIDX(x, y, z + 1)];

                    if (surf)
                        suf_cells[s++] = i;
                }
            }

    struct select_if_not_fixed_ud ud;

    ud.sb = sb;
    ud.cells = vol_cells;
    ud.gaussCenter = bp->gnGaussCenter;
    ud.gaussSigma = bp->gnGaussSigma;
    size_t gr_count = bp->maxGrainDensity * (v * pow(bp->cellSize, 3.0));
    dprintf("Distributing %zu grains in the volume\n", gr_count);
    fprintf(stderr, "Distributing %lu grains in the volume\n", gr_count);
    distribute_nucleations(v, gr_count, select_if_not_fixed, &ud);

    if (bp->doSurfaceNuc && (bp->maxGrainDensitySurf > 0.0))
    {
        ud.cells = suf_cells;
        ud.gaussCenter = bp->gnGaussCenterSurf;
        ud.gaussSigma = bp->gnGaussSigmaSurf;
        gr_count = bp->maxGrainDensitySurf * (s * pow(bp->cellSize, 2.0));
        dprintf("Distributing %zu grains on the surface\n", gr_count);
        distribute_nucleations(s, gr_count, select_if_not_fixed, &ud);
    }

    xfree(suf_cells);
    xfree(vol_cells);
}


/**
 * Given a grain number, retrieve the grain information.
 *
 * \param[in] grnum  Grain Identifying number
 * \return pointer to a grain structure
 */

grain_t *
getGrain2(
    int grnum,
    BB_struct * bp)
{
    grain_t *gr = NULL;
    if (grnum > 0 && grnum <= bp->maxTotalGrains)
    {
        gr = &grain_cache[grnum];
    }
    return gr;
}


/**
 * Determines if a cell is past its nucleation undercooling limit.
 * If so, will call for a new grain to be created.
 *
 * Designed to be called from a (potentially parallel) \c ll_walk()
 *
 * \param[in] vsb  Void pointer to a Subblock
 */
void
cell_nucleation(
    SB_struct * sb,
    void * __attribute__ ((__unused__)) __unused)
{
    int dimx = bp->gsdimx;
    int dimy = bp->gsdimy;
    int dimz = bp->gsdimz;

    double m_solute0_a = bp->m_solute0_a;
    double m_solute0_b = bp->m_solute0_b;
    double m_solute0_c = bp->m_solute0_c;
    double cinit = bp->Cinit;
    int nindex2 = 0;
    double sim_time = bp->timestep * bp->ts_delt;

#ifdef GPU_OMP_NUC
    int numGrainPerSub = (int) bp->maxTotalGrains / nproc;
    if (newGrainLocs == NULL)
    {
        xrealloc(newGrainLocs, loc_t, numGrainPerSub);
#pragma omp target enter data map(to:newGrainLocs[0:numGrainPerSub])
    }
#pragma omp target update to(numNewGrains)
#endif


#ifndef NUC_SEP

    int nng = 0;
    double *temperature = sb->temperature;
    float *nuc_threshold = sb->nuc_threshold;
    int *lsindex = sb->lsindex;
#ifdef GPU_OMP_NUC
#pragma omp target map(from:numNewGrains)
#pragma omp teams distribute
#endif
    for (int k = 1; k <= dimz; k++)
    {
#ifdef GPU_OMP_NUC
#pragma omp parallel for collapse(2) schedule(static,1)
#endif
        for (int j = 1; j <= dimy; j++)
        {
            for (int i = 1; i <= dimx; i++)
            {
                int idx = k * (dimy + 2) * (dimx + 2) + j * (dimx + 2) + i;
                double Ttemp = temperature[idx];

                int layerno = 1;
                int initnuc = 1;

                if (bp->am)
                {
                    if (bp->amlayer == 1)
                        layerno = 1;
                    else if (bp->amlayer > 1 && lsindex[idx] == 1)
                        layerno = 1;
                    else
                        layerno = 0;

                    if (sim_time <= 1e-4)
                        initnuc = 1;
                    else if (Ttemp >= bp->solidusTemp)
                        initnuc = 1;
                    else
                        initnuc = 0;
                }

                if ((!sb->mold[idx]) && (sb->gr[idx] <= 0) && lyaerno
                    && initnuc)
                    // Liquid, and a nuc site!
                {
                    double Tunder = bp->liquidusTemp - Ttemp;
                    double m_solute0 =
                        m_solute0_a * Ttemp * Ttemp +
                        m_solute0_b * Ttemp + m_solute0_c;
                    Tunder += m_solute0 * (sb->cl[idx] - bp->Cinit);
                    // If this isn't a nucleation site, nuc_threshold will be
                    // INFINITY, so this give us a false positive
                    if (Tunder > nuc_threshold[idx])
                    {
                        // This nuc site's threshold has been exceeded.  Nucleate!
                        uint32_t gx = (i - 1) + sb->coords.x * bp->gsdimx;
                        uint32_t gy = (j - 1) + sb->coords.y * bp->gsdimy;
                        uint32_t gz = (k - 1) + sb->coords.z * bp->gsdimz;
#ifndef GPU_OMP_NUC
                        addNewGrain(gx, gy, gz);
#else
                        newGrainLocs[numNewGrains].x = gx;
                        newGrainLocs[numNewGrains].y = gy;
                        newGrainLocs[numNewGrains].z = gz;
#pragma omp atomic capture
                        {
                            nng = nng + 1;
                            numNewGrains = nng;
                        }
#endif // GPU_OMP_NUC
                    }
                }
            }
        }
    }

#else // i.e. NUC_SEP defined

    profile(CALC_NUCLEATION);

#ifdef NUC_THREELIST

    int nindex = 0;
    int *nuc_id = sb->nuc_id;

#ifndef NUC_PRELIST

#ifdef GPU_OMP_NUC
#pragma omp target map(tofrom:nindex)
#pragma omp teams distribute
#endif
    for (int k = 1; k <= dimz; k++)
    {
#ifdef GPU_OMP_NUC
#pragma omp parallel for collapse(2) schedule(static,1)
#endif

        for (int j = 1; j <= dimy; j++)
        {
            for (int i = 1; i <= dimx; i++)
            {
                int idx = k * (dimy + 2) * (dimx + 2) + j * (dimx + 2) + i;
                if (sb->gr[idx] <= 0)
                {
                    int nn = 0;
#pragma omp atomic capture
                    {
                        nindex = nindex + 1;
                        nn = nindex;
                    }
                    nuc_id[nn] = idx;
                }
            }
        }
    }
#else // i.e. NUC_PRELIST defined
    nindex = sb->nindex;
#endif // NUC_PRELIST

//profile(TIMER_1);
    int *nuc_id2 = sb->nuc_id2;
#ifdef GPU_OMP_NUC
#pragma omp target map(tofrom:nindex2)
#pragma omp teams distribute parallel for
#endif
    for (int i = 1; i <= nindex; i++)
    {
        int idx = nuc_id[i];
        double Ttemp = temperature[idx];
        double Tunder = bp->liquidusTemp - Ttemp;
        double m_solute0 =
            m_solute0_a * Ttemp * Ttemp + m_solute0_b * Ttemp + m_solute0_c;
        Tunder += m_solute0 * (sb->cl[idx] - cinit);

        if (Tunder > nuc_threshold[idx])
        {
            {
                int nn = 0;
#pragma omp atomic capture
                {
                    nindex2 = nindex2 + 1;
                    nn = nindex2;
                }
                nuc_id2[nn] = idx;
            }
        }
    }

//profile(TIMER_2);
#else // i.e. NUC_THREELIST not defined

    double *temperature = sb->temperature;
    float *nuc_threshold = sb->nuc_threshold;
    int *nuc_id2 = sb->nuc_id2;
    int *lsindex = sb->lsindex;
int *gr = sb->gr;
double *cl = sb->cl;
#ifdef GPU_OMP_NUC
#pragma omp target map(tofrom:nindex2)
#pragma omp teams distribute
#endif
    for (int k = 1; k <= dimz; k++)
    {
//        int *gr = sb->gr;
//        double *cl = sb->cl;
#ifdef GPU_OMP_NUC
#pragma omp parallel for collapse(2) schedule(static,1)
#endif
        for (int j = 1; j <= dimy; j++)
        {
            for (int i = 1; i <= dimx; i++)
            {
                int idx = k * (dimy + 2) * (dimx + 2) + j * (dimx + 2) + i;
                double Ttemp = temperature[idx];

                int layerno = 1;
                int initnuc = 1;

                if (bp->am)
                {
                    if (bp->amlayer == 1)
                        layerno = 1;
                    else if (bp->amlayer > 1 && lsindex[idx] == 1)
                        layerno = 1;
                    else
                        layerno = 0;

                    if (sim_time <= 1e-4)
                        initnuc = 1;
                    else if (Ttemp >= bp->solidusTemp)
                        initnuc = 1;
                    else
                        initnuc = 0;
                }

                if (gr[idx] <= 0 && layerno && initnuc)
                {
                    double Tunder = bp->liquidusTemp - Ttemp;
                    double m_solute0 =
                        m_solute0_a * Ttemp * Ttemp +
                        m_solute0_b * Ttemp + m_solute0_c;
                    Tunder += m_solute0 * (cl[idx] - cinit);

                    if (Tunder > nuc_threshold[idx])
                    {
                        {
                            int nn = 0;
#pragma omp atomic capture
                            {
                                nindex2 = nindex2 + 1;
                                nn = nindex2;
                            }
                            nuc_id2[nn] = idx;
                        }
                    }
                }
            }
        }
    }

#endif // NUC_THREELIST


    const int coordx = sb->coords.x * bp->gsdimx;
    const int coordy = sb->coords.y * bp->gsdimy;
    const int coordz = sb->coords.z * bp->gsdimz;

#ifdef GPU_OMP_NUC
#pragma omp target teams distribute parallel for
#endif
    for (int i = 1; i <= nindex2; i++)
    {
        int idx = nuc_id2[i];
        int dimxy = (dimx + 2) * (dimy + 2);
        int kk = (int) (idx / dimxy);
        int jj = (int) ((idx - kk * dimxy) / (dimx + 2));
        int ii = (int) (idx - kk * dimxy - jj * (dimx + 2));

        uint32_t gx = (ii - 1) + coordx;
        uint32_t gy = (jj - 1) + coordy;
        uint32_t gz = (kk - 1) + coordz;
#ifndef GPU_OMP_NUC
        addNewGrain(gx, gy, gz);
#else
        newGrainLocs[i - 1].x = gx;
        newGrainLocs[i - 1].y = gy;
        newGrainLocs[i - 1].z = gz;
#endif
    }

//profile(TIMER_3);
    numNewGrains = nindex2;
#endif

#ifdef GPU_OMP_NUC
#pragma omp target update from(newGrainLocs[0:nindex2])
#endif
}
