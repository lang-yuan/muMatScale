/***************************************************************/
/***************************************************************/
/* Yuan, Univeristy of South Carolina		               */
/* 					                       */
/* All rights reserved.                                        */
/***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#define __USE_BSD 1             /* For strdup */
#include <string.h>
#include <strings.h>
#include <time.h>
#include <limits.h>
#include "globals.h"
#include "debug.h"

void
set_defaults(
    )
{
    sprintf(bp->basefilename, "output");
    bp->origin_offset[0] = bp->origin_offset[1] = bp->origin_offset[2] = 0.0;
    bp->ext_temp_filename[0] = '\0';
    bp->padBndyX = bp->padBndyY = bp->padBndyZ = 1;
    bp->cellSize = 0.01f;
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
        else
        {
            fprintf(stderr, "Warning:  Unknown Configuration Option %s\n",
                    token);
        }
    }

    // Determine the number of subblocks in each dimension
    bp->gdimx = bp->gnsbx * bp->gsdimx;
    bp->gdimy = bp->gnsby * bp->gsdimy;
    bp->gdimz = bp->gnsbz * bp->gsdimz;

    fclose(fp);
}
