/***************************************************************/
/* Copyright (c) 2023, Lang Yuan, Univeristy of South Carolina */
/* All rights reserved.                                        */
/* This file is part of muMatScale.                            */
/* See the top-level LICENSE file for details.                 */
/***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define __USE_BSD 1             /* For strdup */
#include <string.h>
#include <strings.h>
#include <mpi.h>
#include <time.h>
#include <limits.h>
#include <float.h>
#include "globals.h"
#include "file_io.h"
#include "debug.h"
#include "functions.h"
#include "grain.h"
#include "face_util.h"

/* Because there's really no fundamental reason to split config files, just
 * use recursion for Geometry File and Material Files
 */

void
set_defaults(
    )
{
    sprintf(bp->basefilename, "output");
    bp->origin_offset[0] = bp->origin_offset[1] = bp->origin_offset[2] = 0.0;
    bp->ext_temp_filename[0] = '\0';
    bp->padBndyX = bp->padBndyY = bp->padBndyZ = 1;
    bp->cellSize = 0.01f;
    bp->temp_type = INTERNAL;
    bp->temp_ctrl.initialTemperature = 100.0f;
    bp->temp_ctrl.cool_rate = 0.0f;
    bp->liquidusTemp = 1180.0f;
    bp->solidusTemp = 1120.0f;
    bp->temp_pure = 0.0f;
    bp->pc_material_group = -1.0;
    bp->pc_region = 0;
    //ADD_CURV_LY
    bp->tip_curv = 0;
    bp->gibbs_coef = 0.0;

    bp->num_grains = 1;         // The "zero" grain

    //bp->Tsol = 1120.0f; // Is this the same as liquidusTemp?

    bp->melt_dt = 0.0f;         // This one is never read in in umatic - always seems to be 0

    bp->base_random_seed = time(NULL);  // Default to random

    bp->fs_finish = 1.0f;
    bp->finish_time = 0.0f;
    bp->ts_delt = 0.1f;
    bp->screenpfreq = 1;
    bp->vis_format = NONE;
    bp->data_write_freq = 0;
    bp->data_write_start = 0;
    bp->checkpointfreq = 0;
    bp->gg_const = 1e-5;
    bp->gg_cub = 0.0;
    bp->profile_write_freq = 0;

    // Used in grain.c
    bp->randomize_grains = 1;
    bp->maxGrainDensity = 0.0;
    bp->gnGaussCenter = 0.0;
    bp->gnGaussSigma = 0.0;
    bp->gnUseOrientation = 0;
    bp->maxGrainDensitySurf = 0.0;
    bp->gnGaussCenterSurf = 0.0;
    bp->gnGaussSigmaSurf = 0.0;
    bp->maxTotalGrains = 1e6;   // Default to a million grains
    bp->doSurfaceNuc = 1;
    bp->pre_num_grains = 0;     //pre-count potential grain numbers

    bp->calc_type = DIFFUSION;

    // Alloy Properties
    bp->m_solute0 = -2.0;
    bp->m_solute0_a = 0.0;
    bp->m_solute0_b = 0.0;
    bp->m_solute0_c = -2.0;

    bp->Cinit = 5.0;
    //bp->Tliq = 1180.0f;
    bp->Dsol = 1e-13;
    bp->Dliq = 1e-9;
    bp->part_coef = 0.3;
    bp->part_coef_a = 0.0;
    bp->part_coef_b = 0.0;
    bp->part_coef_c = 0.3;

    bp->fluidflow = 0;

    //AM
    bp->am = 0;
    bp->thermal_ana = 0;
    bp->amlayer = 1;
    bp->amtracks = 1;
}


void
print_config(
    FILE * fp)
{
    fprintf(fp, "Global Configuration:\n");
    fprintf(fp, "Base filename:  %s\n", bp->basefilename);
    fprintf(fp, "Dimensions (Cells):  %d, %d, %d\n",
            bp->gdimx, bp->gdimy, bp->gdimz);
    fprintf(fp, "Dimensions (Meters):  %g, %g, %g\n",
            bp->gdimx * bp->cellSize,
            bp->gdimy * bp->cellSize, bp->gdimz * bp->cellSize);
    fprintf(fp, "Simulated Volume:  [%lg %lg] x [%lg %lg] x [%lg %lg]\n",
            bp->origin_offset[0],
            (bp->gdimx * bp->cellSize) + bp->origin_offset[0],
            bp->origin_offset[1],
            (bp->gdimy * bp->cellSize) + bp->origin_offset[1],
            bp->origin_offset[2],
            (bp->gdimz * bp->cellSize) + bp->origin_offset[2]);
    fprintf(fp, "Subblock Size (Cells):  %d, %d, %d\n", bp->gsdimx,
            bp->gsdimy, bp->gsdimz);
    fprintf(fp, "Cell Size:  %g m\n", bp->cellSize);
    fprintf(fp, "Num Subblocks:  %d, %d, %d\n",
            bp->gnsbx, bp->gnsby, bp->gnsbz);
    fprintf(fp, "Boundry Conditions:  X: %s, Y: %s, Z: %s\n",
            bp->padBndyX ? "Pad" : "Wrap",
            bp->padBndyY ? "Pad" : "Wrap", bp->padBndyZ ? "Pad" : "Wrap");
    //ADD_CURV_LY
    fprintf(fp, "Dendrite tip curvature : %d \n", bp->tip_curv);
    if (bp->tip_curv)
        fprintf(fp, "Gibbs Thomson Coefficient : %e \n", bp->gibbs_coef);

    if (bp->fluidflow)
    {
        fprintf(fp, "Fluid flow mode is ON \n");
        fprintf(fp, "\tInitial velosity: %e %e %e \n", bp->initvelo[0],
                bp->initvelo[1], bp->initvelo[2]);
    }
    fprintf(fp, "Initial Temperature: %g\n",
            bp->temp_ctrl.initialTemperature);
    fprintf(fp, "Temperature gradient: %g\n", bp->temp_ctrl.gradient);
    fprintf(fp, "Casting speed: %g\n", bp->temp_ctrl.velocity);
    fprintf(fp, "Liquidus Temperature: %g\n", bp->liquidusTemp);
    fprintf(fp, "Solidus Temperature: %g\n", bp->solidusTemp);
    fprintf(fp, "Liquidus slope: %g\n", bp->m_solute0);
    fprintf(fp, "Partition Coef: %g\n", bp->part_coef);
    fprintf(fp, "Time Step: %g sec/iter\n", bp->ts_delt);
    fprintf(fp, "Volume Grain Info:  Density: %g Center: %g  Sigma: %g\n",
            bp->maxGrainDensity, bp->gnGaussCenter, bp->gnGaussSigma);
    fprintf(fp, "Surface Grain Info:  Density: %g Center: %g  Sigma: %g\n",
            bp->maxGrainDensitySurf, bp->gnGaussCenterSurf,
            bp->gnGaussSigmaSurf);
    fprintf(fp, "Beginning simulation at time %g (sec)\n",
            bp->timestep * bp->ts_delt);
    fprintf(fp, "Auto-Terminate time:  %g (sec) (0=disable)\n",
            bp->finish_time);
    fprintf(fp, "Solid Fraction Termination: %g%%\n", 100 * bp->fs_finish);
    fprintf(fp, "Frequencey to Write to Screen:  %d (timesteps)\n",
            bp->screenpfreq);
    fprintf(fp, "Datafile format: %s\n", getIOTypeName());
    fprintf(fp, "Frequencey to Write Datafile:  %d (timesteps)\n",
            bp->data_write_freq);
    fprintf(fp, "Frequencey to Write Profile:  %d (timesteps)\n",
            bp->profile_write_freq);
    fprintf(fp, "Timestep to begin Writing Datafile:  %g (timesteps)\n",
            bp->data_write_start * bp->ts_delt);
    fprintf(fp, "Frequencey to Checkpoint:  %d (timesteps)\n",
            bp->checkpointfreq);
    fprintf(fp, "GG_const:  %g\n", bp->gg_const);
    fprintf(fp, "GG_cub:  %g\n", bp->gg_cub);

    fprintf(fp, "T_pure: %g\n", bp->temp_pure);
    //fprintf(fp, "Tliq: %g\n", bp->Tliq);
    fprintf(fp, "Dliq: %g\n", bp->Dliq);
    fprintf(fp, "Dsol: %g\n", bp->Dsol);
    fprintf(fp, "Cinit: %g\n", bp->Cinit);

    if (bp->am)
    {
        if (bp->temp_type != INTERNAL)
        {
            fprintf(fp, "Temp_type has to be INTERNAL (enforced now) \n");
            bp->temp_type = INTERNAL;
        }
        if (bp->thermal_ana == 1)
        {
            fprintf(fp,
                    "Analytical Thermal is ON, based on Rosenthal solution \n");
        }
        else if (bp->thermal_ana == 2)
        {
            fprintf(fp,
                    "Analytical Thermal is ON, based on 3D double ellipsoidal heat source by Goldak \n");
            fprintf(fp, "\t Golak semi-axis for side: %g \n", bp->lga);
            fprintf(fp, "\t Golak semi-axis for depth: %g \n", bp->lgb);
            fprintf(fp, "\t Golak semi-axis for rear: %g \n", bp->lgcr);
            fprintf(fp, "\t Golak semi-axis for front: %g \n", bp->lgcf);
        }

        fprintf(fp, "\t Laser Power: %g W\n", bp->lpow);
        fprintf(fp, "\t Laser Speed: %g m/s\n", bp->lv);
        fprintf(fp, "\t Laser Spot Size: %g m\n", bp->lsd);
        fprintf(fp, "\t Laser Hatch Spacing: %g m\n", bp->lhs);
        fprintf(fp, "\t Laser Absorption Coef.: %g \n", bp->lab);

    }
}

/*
void SetCalculatedValues(void)
{
    bp->Tliq = bp->temp_pure + (bp->m_solute0 * bp->Cinit);
}
*/

void
verify_config(
    void)
{
    if (bp->temp_pure == 0.0 && bp->m_solute0_a == 0 && bp->m_solute0_b == 0)
        bp->temp_pure = bp->liquidusTemp - (bp->m_solute0_c * bp->Cinit);
}


void
read_config(
    const char *ctrl_fname)
{
    char line[MAX_STRING_LEN + 2];
    char token[MAX_STRING_LEN + 2];
    char strarg[MAX_STRING_LEN + 2];

    FILE *fp = fopen(ctrl_fname, "r");
    if (!fp)
    {
        error("Unable to open control file: %s\n", ctrl_fname);
    }

    while (fgets(line, MAX_STRING_LEN, fp) != NULL)
    {
        int idx;
        char *lp;
        int c;
        lp = line;
        c = sscanf(lp, " %s%n", token, &idx);
        lp += idx;
        dwrite(DEBUG_PARSING, "Parsing line:  '%s':  '%s'\n", line, token);
        if ((c <= 0) || (token[0] == '#' || token[0] == '%'))
        {
            continue;
        }
        else if (strcasecmp(token, "BaseFileName") == 0)
        {
            c = sscanf(lp, "%s%n", bp->basefilename, &idx);
            assert(c == 1);
            lp += idx;
        }
        else if (strcasecmp(token, "NSubblocks") == 0)
        {
            c = sscanf(lp, "%d %d %d", &bp->gnsbx, &bp->gnsby, &bp->gnsbz);
            assert(c == 3);
        }
        else if (strcasecmp(token, "NCellsPerSB") == 0)
        {
            c = sscanf(lp, "%d %d %d", &bp->gsdimx, &bp->gsdimy, &bp->gsdimz);
            assert(c == 3);
        }
        else if (strcasecmp(token, "CellSize") == 0)
        {
            c = sscanf(lp, "%lf", &bp->cellSize);
            assert(c == 1);
        }
        else if (strcasecmp(token, "FaceCtrl") == 0)
        {
            int x, y, z;
            c = sscanf(lp, "%d %d %d %d %d %d",
                       &bp->padBndyX, &x, &bp->padBndyY, &y, &bp->padBndyZ,
                       &z);
            assert(c == 6);
            assert(bp->padBndyX == x);
            assert(bp->padBndyY == y);
            assert(bp->padBndyZ == z);
        }
        else if (strcasecmp(token, "FixedNuc") == 0)
        {
            double x, y, z;
            double ax, ay, az;
            double tsh;         //threshold
            // There are six values for each fixed nucleation site
            // location in x,y,z and preferred growth angle in x,y,z
            c = sscanf(lp, "%lf %lf %lf %lf %lf %lf %lf",
                       &x, &y, &z, &ax, &ay, &az, &tsh);
            assert(c == 7);
            bp->pre_num_grains++;
            add_fixed_nuc(x, y, z, ax, ay, az, tsh);
        }
        else if (strcasecmp(token, "NucVolume") == 0)
        {
            double density;
            double minx, miny, minz, maxx, maxy, maxz;
            double tsh;
            // # nuc points
            // minx, miny, minz, maxx, maxy, maxz
            c = sscanf(lp, "%lg %lg %lg %lg %lg %lg %lg %lg\n",
                       &density, &minx, &miny, &minz, &maxx, &maxy, &maxz,
                       &tsh);
            assert(c == 8);
            int nuc_num = 0;
            nuc_num =
                ceil((maxx - minx) * (maxy - miny) * (maxz -
                                                      minz) / bp->cellSize /
                     bp->cellSize / bp->cellSize * density / 100);
            bp->pre_num_grains += nuc_num;
            add_nuc_volume(density / 100.0, minx, miny, minz, maxx, maxy,
                           maxz, tsh);
        }
        else if (strcasecmp(token, "MaxGrainDensity") == 0)
        {                       /* Note, this is max density per m^3 */
            c = sscanf(lp, "%lf", &bp->maxGrainDensity);
            assert(c == 1);
        }
        else if (strcasecmp(token, "MaxGrainDensitySurf") == 0)
        {                       /* Note, this is max density per m^2 */
            c = sscanf(lp, "%lf", &bp->maxGrainDensitySurf);
            assert(c == 1);
        }
        else if (strcasecmp(token, "MaxTotGrains") == 0)
        {
            c = sscanf(lp, "%lu", &bp->maxTotalGrains);
            assert(c == 1);
        }
        else if (strcasecmp(token, "RandomizeGrains") == 0)
        {
            c = sscanf(lp, "%d", &bp->randomize_grains);
            assert(c == 1);
        }
        else if (strcasecmp(token, "MouldNuc") == 0
                 || strcasecmp(token, "MoldNuc") == 0)
        {
            c = sscanf(lp, "%d", &bp->doSurfaceNuc);
            assert(c == 1);
        }
        else if (strcasecmp(token, "GNGaussCentre") == 0
                 || strcasecmp(token, "GNGaussCenter") == 0)
        {
            c = sscanf(lp, "%lf", &bp->gnGaussCenter);
            assert(c == 1);
        }
        else if (strcasecmp(token, "GNGaussSigma") == 0)
        {
            c = sscanf(lp, "%lf", &bp->gnGaussSigma);
            assert(c == 1);
        }
        else if (strcasecmp(token, "GNGaussCentreSurf") == 0
                 || strcasecmp(token, "GNGaussCenterSurf") == 0)
        {
            c = sscanf(lp, "%lf", &bp->gnGaussCenterSurf);
            assert(c == 1);
        }
        else if (strcasecmp(token, "GNGaussSigmaSurf") == 0)
        {
            c = sscanf(lp, "%lf", &bp->gnGaussSigmaSurf);
            assert(c == 1);
        }
        else if (strcasecmp(token, "GNOriented") == 0)
        {
            c = sscanf(lp, "%d", &bp->gnUseOrientation);
            assert(c == 1);
        }
        else if (strcasecmp(token, "LiquidusTemp") == 0)
        {
            c = sscanf(lp, "%lf", &bp->liquidusTemp);
            assert(c == 1);
        }
        else if (strcasecmp(token, "SolidusTemp") == 0)
        {
            c = sscanf(lp, "%lf", &bp->solidusTemp);
            assert(c == 1);
        }
        else if (strcasecmp(token, "T_pure") == 0)
        {
            c = sscanf(lp, "%lf", &bp->temp_pure);
            assert(c == 1);
        }
        else if (strcasecmp(token, "Cinit") == 0)
        {
            c = sscanf(lp, "%lf", &bp->Cinit);
            assert(c == 1);
        }
        else if (strcasecmp(token, "m_solute0") == 0)
        {
            c = sscanf(lp, "%lf %lf %lf", &bp->m_solute0_a, &bp->m_solute0_b,
                       &bp->m_solute0_c);
            bp->m_solute0 = bp->m_solute0_c;
            assert(c == 3);
        }
        else if (strcasecmp(token, "part_coef0") == 0)
        {
            c = sscanf(lp, "%lf %lf %lf", &bp->part_coef_a, &bp->part_coef_b,
                       &bp->part_coef_c);
            bp->part_coef = bp->part_coef_c;
            assert(c == 3);
        }
        else if (strcasecmp(token, "Dsol") == 0)
        {
            c = sscanf(lp, "%lf", &bp->Dsol);
            assert(c == 1);
        }
        else if (strcasecmp(token, "Dliq") == 0)
        {
            c = sscanf(lp, "%lf", &bp->Dliq);
            assert(c == 1);
        }
        else if (strcasecmp(token, "CalculationModel") == 0)
        {
            c = sscanf(lp, "%s%n", strarg, &idx);
            assert(c == 1);
            lp += idx;
            bp->calc_type = DIFFUSION;
            {
                error("Unknown Calculation Model '%s'.\n"
                      "Should be: DIFFUSION\n", strarg);
            }
        }
        else if (strcasecmp(token, "TemperatureModel") == 0)
        {
            c = sscanf(lp, "%s%n", strarg, &idx);
            assert(c == 1);
            lp += idx;
            if (strcasecmp(strarg, "INTERNAL") == 0)
                bp->temp_type = INTERNAL;
            else
            {
                error("Unknown Temperature Model '%s'.\n"
                      "Should be :  INTERNAL\n", strarg);
            }
        }
        else if (strcasecmp(token, "ExternalTempFile") == 0)
        {
            c = sscanf(lp, "%s%n", bp->ext_temp_filename, &idx);
            assert(c == 1);
            lp += idx;
        }
        else if (strcasecmp(token, "ProcastMaterialGroup") == 0)
        {
            c = sscanf(lp, "%d", &bp->pc_material_group);
            assert(c == 1);
        }
        else if (strcasecmp(token, "CalculateProcastRegion") == 0)
        {
            c = sscanf(lp, "%d", &bp->pc_region);
            assert(c == 1);
        }
        else if (strcasecmp(token, "ProcastRegion") == 0)
        {
            c = sscanf(lp, "%lg %lg %lg %lg %lg %lg\n",
                       &bp->pc_region_vol.minX,
                       &bp->pc_region_vol.minY,
                       &bp->pc_region_vol.minZ,
                       &bp->pc_region_vol.maxX,
                       &bp->pc_region_vol.maxY, &bp->pc_region_vol.maxZ);
            assert(c == 6);
        }
        //ADD_CURV_LY
        //-------------------------------------------
        else if (strcasecmp(token, "TipCurvature") == 0)
        {
            c = sscanf(lp, "%d", &bp->tip_curv);
            assert(c == 1);
        }
        else if (strcasecmp(token, "GibbsThomsonCoef") == 0)
        {
            c = sscanf(lp, "%lg", &bp->gibbs_coef);
            assert(c == 1);
        }
        //-------------------------------------------

        else if (strcasecmp(token, "InitialTemperature") == 0)
        {
            c = sscanf(lp, "%lf", &bp->temp_ctrl.initialTemperature);
            assert(c == 1);
        }
        else if (strcasecmp(token, "Velo_Coef") == 0)
        {
            c = sscanf(lp, "%lf", &bp->temp_ctrl.velo_coef);
            assert(c == 1);
        }
        else if (strcasecmp(token, "Grad_Coef") == 0)
        {
            c = sscanf(lp, "%lf", &bp->temp_ctrl.grad_coef);
            assert(c == 1);
        }
        else if (strcasecmp(token, "GradSlope") == 0)
        {
            c = sscanf(lp, "%lf", &bp->temp_ctrl.grad_slope);
            assert(c == 1);
        }
        else if (strcasecmp(token, "Slope_coef") == 0)
        {
            c = sscanf(lp, "%lf", &bp->temp_ctrl.slope_coef);
            assert(c == 1);
        }
        else if (strcasecmp(token, "Gradient") == 0)
        {
            c = sscanf(lp, "%lf", &bp->temp_ctrl.gradient);
            assert(c == 1);
        }
        else if (strcasecmp(token, "Velocity") == 0)
        {
            c = sscanf(lp, "%lf", &bp->temp_ctrl.velocity);
            assert(c == 1);
        }
        else if (strcasecmp(token, "Iso_Coef_One") == 0)
        {
            c = sscanf(lp, "%lf", &bp->temp_ctrl.iso_coef1);
            assert(c == 1);
        }
        else if (strcasecmp(token, "Iso_coef_Two") == 0)
        {
            c = sscanf(lp, "%lf", &bp->temp_ctrl.iso_coef2);
            assert(c == 1);
        }
        else if (strcasecmp(token, "Coef_Iso2") == 0)
        {
            c = sscanf(lp, "%lf", &bp->temp_ctrl.coef_iso2);
            assert(c == 1);
        }
        else if (strcasecmp(token, "CoolingRate") == 0)
        {
            c = sscanf(lp, "%lf", &bp->temp_ctrl.cool_rate);
            assert(c == 1);
        }

        else if (strcasecmp(token, "FSGrow") == 0)
        {
            c = sscanf(lp, "%lf", &bp->fsgrow);
            bp->rev_fsgrow = 1.0 / bp->fsgrow;
            assert(c == 1);
        }
        else if (strcasecmp(token, "GG_Constant") == 0)
        {
            c = sscanf(lp, "%lf", &bp->gg_const);
            assert(c == 1);
        }
        else if (strcasecmp(token, "GG_Cub") == 0)
        {
            c = sscanf(lp, "%lf", &bp->gg_cub);
            assert(c == 1);
        }
        else if (strcasecmp(token, "TimeStep") == 0)
        {
            c = sscanf(lp, "%lf", &bp->ts_delt);
            assert(c == 1);
        }
        else if (strcasecmp(token, "FinishTime") == 0)
        {
            c = sscanf(lp, "%lf", &bp->finish_time);
            assert(c == 1);
        }
        else if (strcasecmp(token, "FsFinish") == 0)
        {
            c = sscanf(lp, "%lf", &bp->fs_finish);
            assert(c == 1);
        }
        else if (strcasecmp(token, "ScreenPFreq") == 0)
        {
            c = sscanf(lp, "%d", &bp->screenpfreq);
            assert(c == 1);
        }
        else if (strcasecmp(token, "DataWriteFreq") == 0)
        {
            c = sscanf(lp, "%d", &bp->data_write_freq);
            assert(c == 1);
        }
        else if (strcasecmp(token, "ProfileWriteFreq") == 0)
        {
            c = sscanf(lp, "%d", &bp->profile_write_freq);
            assert(c == 1);
        }
        else if (strcasecmp(token, "DataWriteStart") == 0)
        {
            c = sscanf(lp, "%d", &bp->data_write_start);
            assert(c == 1);
            c = sscanf(lp, "%d", &bp->checkpointfreq);
            assert(c == 1);
        }
        else if (strcasecmp(token, "DataFormat") == 0)
        {
            c = sscanf(lp, "%s%n", strarg, &idx);
            assert(c == 1);
            lp += idx;
            if (strcasecmp(strarg, "VTK") == 0)
            {
                bp->vis_format = VTK;
            }
            else if (strcasecmp(strarg, "XDMF") == 0)
            {
                bp->vis_format = XDMF;
            }
            else
            {
                fprintf(stderr, "Supported File types:\n"
                        "\tNONE    (Do not write data files)\n"
                        "\tVTK     (Parallel VTK XML format)\n"
                        "\tXDMF    (XDMF XML+HDF5 files)\n");
                error("Unsupported data file type: %s\n\n", strarg);
            }
        }
        else if ((strcasecmp(token, "GeoFileName") == 0) ||
                 (strcasecmp(token, "MatFileName") == 0) ||
                 (strcasecmp(token, "AlloyPropsFile0") == 0))
        {
            c = sscanf(lp, "%s%n", strarg, &idx);
            assert(c == 1);
            lp += idx;
            read_config(strarg);
        }
        else if (strcasecmp(token, "RandSeedVal") == 0)
        {
            c = sscanf(lp, "%u", &bp->base_random_seed);
            assert(c == 1);
        }
//
//      Parameters for fluid flow module
//
        else if (strcasecmp(token, "FluidFlowMode") == 0)
        {
            c = sscanf(lp, "%d", &bp->fluidflow);
            assert(c == 1);
        }
        else if (strcasecmp(token, "FluidFlowTimeStep") == 0)
        {
            c = sscanf(lp, "%lf", &bp->ffdelt);
            assert(c == 1);
        }
        else if (strcasecmp(token, "FFStartTime") == 0)
        {
            c = sscanf(lp, "%lf", &bp->ffstarttime);
            assert(c == 1);
        }
        else if (strcasecmp(token, "FFStepGap") == 0)
        {
            c = sscanf(lp, "%d", &bp->ffstepgap);
            assert(c == 1);
        }
        else if (strcasecmp(token, "Density") == 0)
        {
            c = sscanf(lp, "%lf", &bp->rho);
            assert(c == 1);
        }
        else if (strcasecmp(token, "Viscosity") == 0)
        {
            c = sscanf(lp, "%lf", &bp->viscosity);
            assert(c == 1);
        }
        else if (strcasecmp(token, "ThermExpCoe") == 0)
        {
            c = sscanf(lp, "%lf", &bp->beta_T);
            assert(c == 1);
        }
        else if (strcasecmp(token, "SoluteExpCoe") == 0)
        {
            c = sscanf(lp, "%lf", &bp->beta_c);
            assert(c == 1);
        }
        else if (strcasecmp(token, "ReferenceTemp") == 0)
        {
            c = sscanf(lp, "%lf", &bp->ref_T);
            assert(c == 1);
        }
        else if (strcasecmp(token, "ReferSolute") == 0)
        {
            c = sscanf(lp, "%lf", &bp->ref_c);
            assert(c == 1);
        }
        else if (strcasecmp(token, "Gravity_x") == 0)
        {
            c = sscanf(lp, "%lf", &bp->gravity_x);
            assert(c == 1);
        }
        else if (strcasecmp(token, "Gravity_y") == 0)
        {
            c = sscanf(lp, "%lf", &bp->gravity_x);
            assert(c == 1);
        }
        else if (strcasecmp(token, "Gravity_z") == 0)
        {
            c = sscanf(lp, "%lf", &bp->gravity_x);
            assert(c == 1);
        }
        else if (strcasecmp(token, "InitialVelocity") == 0)
        {
            c = sscanf(lp, "%lf %lf %lf", &bp->initvelo[0], &bp->initvelo[1],
                       &bp->initvelo[2]);
            assert(c == 3);
        }
//Parameters for laser thermal analytical thermal
        else if (strcasecmp(token, "AnalyticalThermal") == 0)
        {
            c = sscanf(lp, "%d", &bp->thermal_ana);
            assert(c == 1);
            if (bp->thermal_ana)
                bp->am = 1;
        }
        else if (strcasecmp(token, "LaserPower") == 0)
        {
            c = sscanf(lp, "%lf", &bp->lpow);
            assert(c == 1);
        }
        else if (strcasecmp(token, "LaserSpeed") == 0)
        {
            c = sscanf(lp, "%lf", &bp->lv);
            assert(c == 1);
        }
        else if (strcasecmp(token, "LaserSpotSize") == 0)
        {
            c = sscanf(lp, "%lf", &bp->lsd);
            assert(c == 1);
        }
        else if (strcasecmp(token, "LaserAbsorption") == 0)
        {
            c = sscanf(lp, "%lf", &bp->lab);
            assert(c == 1);
        }
        else if (strcasecmp(token, "LaserHatchSpacing") == 0)
        {
            c = sscanf(lp, "%lf", &bp->lhs);
            assert(c == 1);
        }
        else if (strcasecmp(token, "SubInitialTemp") == 0)
        {
            c = sscanf(lp, "%lf", &bp->liniT);
            assert(c == 1);
        }
        else if (strcasecmp(token, "ThermalConductivity") == 0)
        {
            c = sscanf(lp, "%lf", &bp->lthcon);
            assert(c == 1);
        }
        else if (strcasecmp(token, "SpecificHeat") == 0)
        {
            c = sscanf(lp, "%lf", &bp->lcp);
            assert(c == 1);
        }
        else if (strcasecmp(token, "LaserStartLocation") == 0)
        {
            c = sscanf(lp, "%lf %lf %lf", &bp->x0, &bp->y0, &bp->z0);
            assert(c == 3);
        }
        else if (strcasecmp(token, "NoOfLayers") == 0)
        {
            c = sscanf(lp, "%d", &bp->amlayer);
            assert(c == 1);
        }
        else if (strcasecmp(token, "NoOfTracks") == 0)
        {
            c = sscanf(lp, "%d", &bp->amtracks);
            assert(c == 1);
        }
        else if (strcasecmp(token, "SemiaxisSide") == 0)
        {
            c = sscanf(lp, "%lf", &bp->lga);
            assert(c == 1);
        }
        else if (strcasecmp(token, "SemiaxisDepth") == 0)
        {
            c = sscanf(lp, "%lf", &bp->lgb);
            assert(c == 1);
        }
        else if (strcasecmp(token, "SemiaxisRear") == 0)
        {
            c = sscanf(lp, "%lf", &bp->lgcr);
            assert(c == 1);
        }
        else if (strcasecmp(token, "SemiaxisFront") == 0)
        {
            c = sscanf(lp, "%lf", &bp->lgcf);
            assert(c == 1);
        }
#if 0
        else
        {
            fprintf(stderr, "Warning:  Unknown Configuration Option %s\n",
                    token);
        }
#endif
    }

    // Determine the number of subblocks in each dimension
    bp->gdimx = bp->gnsbx * bp->gsdimx;
    bp->gdimy = bp->gnsby * bp->gsdimy;
    bp->gdimz = bp->gnsbz * bp->gsdimz;

    verify_config();
//    SetCalculatedValues();

    fclose(fp);
}
