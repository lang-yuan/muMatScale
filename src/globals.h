/***************************************************************/
/* Copyright (c) 2023, Lang Yuan, Univeristy of South Carolina */
/* All rights reserved.                                        */
/* This file is part of muMatScale.                            */
/* See the top-level LICENSE file for details.                 */
/***************************************************************/

#ifndef __GLOBALS_H__
#define __GLOBALS_H__
#include <inttypes.h>
#include <stdio.h>
#include <mpi.h>

//#include "ll.h"

//#define OMP_RELEASE

#define MAX(x,y)    (((x) > (y))? (x) : (y))
#define MIN(x,y)    (((x) < (y))? (x) : (y))
#define MINMAX(min,x,max) \
	do { \
		__typeof__(x) ___x = (x); \
		if ( ___x < (min) ) \
			min = ___x; \
		else if ( ___x > (max) ) \
			max = ___x; \
	} while (0)

//#define ABS(x)      (((x) >= 0)? (x): (-(x)))
#define ABS(x)      fabs(x)
#define IDX(z,y,x)  ((x) + ((y)*NX) + ((z)*NX*NY))

//ADD_CURV_LY
#define SGN(x)  	((x > 0) ? 1 : ((x < 0) ? -1 : 0))

#ifndef __cplusplus
#define true 1
#define false 0
#endif

// Message types
#define SYNCTAG 0
#define ACTIVATETAG 1
#define GRAINTAG 2

#define FACE_TOP     0
#define FACE_BOTTOM  1
#define FACE_LEFT    2
#define FACE_RIGHT   3
#define FACE_FRONT   4
#define FACE_BACK    5

#define MAX_STRING_LEN 256

#define NUM_NEIGHBORS   6
#define NBOR_INFO       2
                         /* NBOR_INFO = Holds information for each neighbor.
                          * Currently, NBOR_INFO = 0 ==> MPI Rank where neighbor lives
                          *            NBOR_INFO = 1 ==> unused
                          */


#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif
#define FLAT_GROWTH    7e-5
#define LEAST_GROWTH_STEP  5    //each time step, at least three time steps to finish...
#define COURANT_LIMIT      0.2  //diffusion stability

// Subblock IDX (index into a SB array)
#define SBIDX(x, y, z) \
	((x) + \
	 ((y)*(bp->gsdimx+2)) + \
	 ((z)*(bp->gsdimy+2)*(bp->gsdimx+2)))
// Subblock IDX (index into a SB array that has No Halo)
// Useful for when iterating over a halo'd array, and wanting access to
// the "same" spot in a non-halo'd array.
#define SBIDXNH(x, y, z) \
		(((x)-1) + \
		 (((y)-1)*(bp->gsdimx)) + \
		 (((z)-1)*(bp->gsdimy)*(bp->gsdimx)))


typedef int neighbors_t[NUM_NEIGHBORS][NBOR_INFO];

// If you change this order, make sure you update file_io.c
typedef enum vis_format
{
    NONE = 0,
    VTK,                        /* Paraview default */
    XDMF,
} vis_format_t;

typedef enum calc_type
{
    DIFFUSION = 0,
} calc_type_t;

typedef enum temp_type
{
    INTERNAL = 0,
} temp_type_t;

typedef struct commTypes
{
    size_t dsize;
} commTypes_t;

typedef struct
{
    double x;
    double y;
    double z;
} dloc_t;


#ifdef GPU_OMP
#pragma omp declare target
#endif
typedef struct
{
    uint32_t x;
    uint32_t y;
    uint32_t z;
} loc_t;

#ifdef GPU_OMP
#pragma omp end declare target
#endif

typedef struct volume_s
{
    double minX;
    double minY;
    double minZ;
    double maxX;
    double maxY;
    double maxZ;
} volume_t;

typedef struct decentered_s
{
    double x, y, z;
} __attribute__ ((packed)) decentered_t;


typedef struct
{
    double x;                   // Nuc Location
    double y;
    double z;

    size_t cx, cy, cz;          // Global Cell Index

    double gx;                  // Growth Direction
    double gy;
    double gz;

    double threshold;
} nucleation_t;

typedef enum
{
    INACTIVE = 0,
    ACTIVE,
} activestate_t;

// subblockid is i + j * bp->gnsbx + k * bp->gnsbx * bp->gnsby
/* Note, if you update this, make sure to update globals.c: CopyMSBtoSB() */
#define COMMON_SUBBLOCK_DEFS \
	neighbors_t neighbors;   \
    loc_t coords;            \
    int subblockid;          \
    int procid;              \
	int nnuc;				 \
	nucleation_t *nuc_pt


typedef struct mainsubblock
{
    COMMON_SUBBLOCK_DEFS;
    size_t nuc_pt_len;
    uint64_t activate_ts;
} MSB_struct;

typedef struct subblock
{
    // Reminder: These should be replicated in the MSB_struct
    // Each variable must be copied from the MSB_struct to the SB_struct
    // Copying is done in the CopyMSBtoSB() in globals.c
    COMMON_SUBBLOCK_DEFS;
    int last_active_ts;

    // Reminder: These are not in the mainsubblock
    void *temp_info;
    double *temperature;
    int8_t *mold;
    size_t mold_size;


    double *fs;
    double *ce, *oce;
    double *cl;

    double *curv;

    int *diff_id;
    int *fs_id;
    int fsindex;
    int gindex;
    int nindex;
    int growindex;
    int *nuc_id;
    int *nuc_id2;
    int totaldim;
    //int nindex;

    // Decentered octahedron variables
    decentered_t *dc;
    double *d;                  // diagonal of decentered octahedron
    int *gr;                    // grain number

    int *ogr;                   /// added to make sure it won't double use the gr number... ly

    float *nuc_threshold;

    // Fluid flow
    double *cell_u;
    double *cell_v;
    double *cell_w;

    //AM modle
    int *lsindex;

} SB_struct;


// Note, if you change this, update 'grainSetup' in grain.c
typedef struct
{
    double nuc_x;
    double nuc_y;
    double nuc_z;
    uint64_t nuc_timestep;
    double nuc_threshold;
    // Rotation matrix
    double rotmat[3][3];
    double rotang[3];
} grain_t;

typedef struct internal_temp_ctrl
{
    double initialTemperature;
    double velo_coef;
    double grad_coef;
    double grad_slope;
    double slope_coef;
    double gradient;
    double velocity;
    double iso_coef1;
    double iso_coef2;
    double coef_iso2;
    double cool_rate;
} internal_temp_ctrl_t;


typedef struct bbstruct
{
    // Global information
    int gdimx, gdimy, gdimz;    // global dimensions of simulation in cells
    int gsdimx, gsdimy, gsdimz; // global subblock dimensions in cells
    int gnsbx, gnsby, gnsbz;    // global dimensions in number of subblocks
    int padBndyX, padBndyY, padBndyZ;   // 1 = Pad, 0 = wrap
    double origin_offset[3];    // In meters, offset from 0,0,0 to minx, miny, minz of the array
    double cellSize;

    unsigned int base_random_seed;

    uint64_t timestep;

    // Communication Types
    MPI_Datatype MPI_DECENTERED;

    // Constants
    calc_type_t calc_type;
    temp_type_t temp_type;
    char ext_temp_filename[256];
    int pc_material_group;
    internal_temp_ctrl_t temp_ctrl;

    int pc_region;              //calculate regions
    volume_t pc_region_vol;

    //ADD_CURV_LY
    //----------
    int tip_curv;
    double gibbs_coef;
    //-----------

    double fsgrow;
    double rev_fsgrow;
    double liquidusTemp;
    double solidusTemp;
    double temp_pure;
    double ts_delt;

    int num_grains;             // Count includes grain 0 -> This is the extent of the grain_cache array usage
    grain_t *grain_cache;
    uint64_t maxTotalGrains;
    int pre_num_grains;
    int randomize_grains;
    int doSurfaceNuc;           // boolean - Use Surface/Mold based nucleation

    // Gaussian Grain Info
    double maxGrainDensity;
    double gnGaussCenter;
    double gnGaussSigma;
    int gnUseOrientation;

    double maxGrainDensitySurf;
    double gnGaussCenterSurf;
    double gnGaussSigmaSurf;

    // Alloy(0) Properties  Should be an array of these?
    double Cinit;
    double m_solute0;
    double m_solute0_a;
    double m_solute0_b;
    double m_solute0_c;

    double Dsol;
    double Dliq;
//    double Tliq;
    double part_coef;
    double part_coef_a;
    double part_coef_b;
    double part_coef_c;
    double melt_dt;

    // Growth variables
    double gg_const;
    double gg_cub;

    // Termination Conditions
    double fs_finish;           // Fraction of Solid to finish at
    double finish_time;

    // Current status
    double fs;                  // Solid Fraction

    // Output Control
    int screenpfreq;
    char basefilename[256];
    vis_format_t vis_format;
    int data_write_freq;
    int data_write_start;
    int checkpointfreq;
    int profile_write_freq;

    // Fluid Flow
    int fluidflow;
    int ffstepgap;
    double ffstarttime;
    double ffdelt;

    double rho;                 /* density [g/mm^3]                     */
    double beta_T;
    double beta_c;
    double ref_T;
    double ref_c;
    double viscosity;
    double gravity_x;
    double gravity_y;
    double gravity_z;
    double initvelo[3];

    // laser analytical thermal solution
    int thermal_ana;
    int am;
    double lpow;
    double lv;
    double lsd;
    double lab;
    double lhs;
    double liniT;
    double lthcon;
    double lcp;
    double x0, y0, z0;

    int amlayer;                // number of layers
    int amtracks;               // number of tracks

    //Goldak input
    double lga, lgb, lgcr, lgcf;

} BB_struct;


// If you update this, update 'functions.c:setup_mpi_datatypes'
typedef struct nbr_info_s
{
    int id;
    neighbors_t nbors;
} nbr_info_t;


extern int iproc;               // My process rank
extern int nproc;               // Number of processors

extern MSB_struct *gmsp;

#ifdef GPU_OMP
#pragma omp declare target
#endif
extern BB_struct *bp;
#ifdef GPU_OMP
#pragma omp end declare target
#endif

extern const int gface_reverse[NUM_NEIGHBORS];

extern int gMPIThreadLevel;
extern MPI_Datatype MPI_Type_Neighbor;

extern MPI_Comm mpi_comm_new;

void CopyMSBtoSB(
    MSB_struct * m,
    SB_struct * s);

#endif
