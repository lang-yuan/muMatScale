/***************************************************************/
/* Copyright (c) 2023, Lang Yuan, Univeristy of South Carolina */
/* All rights reserved.                                        */
/* This file is part of muMatScale.                            */
/* See the top-level LICENSE file for details.                 */
/***************************************************************/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>

#include "globals.h"
#include "functions.h"
#include "grain.h"
#include "debug.h"
#include "temperature.h"


/**
 * ADD_CURV_LY
 * Calculate curvature based on the local solid fraction
 * Based on Parker and Youngs method to determine the interface normal
 * Use Height function to calucate curvature... second-order accuracy
 *
*/


void
sb_curvature(
    void *vlsp,
    void * __attribute__ ((__unused__)) __unused)
{
    SB_struct *lsp = (SB_struct *) vlsp;

    int dimx = bp->gsdimx;
    int dimy = bp->gsdimy;
    int dimz = bp->gsdimz;

    int idx;
    int idl, idr, idf, idb, idu, idd;
    int idur, idul, iddr, iddl, idfr, idfl, idbr, idbl, idub, iduf, iddf,
        iddb;
    int idufr, idufl, idubr, idubl, iddbr, iddbl, iddfl, iddfr;

    double fsl, fsr, fsf, fsb, fsu, fsd;

    double normal_x, normal_y, normal_z;

    // smoother parameters for better accurary
    // Suggest not to change... LY
    int alpha = 2;
    int beta = 4;
    double gamma;

    int n, i, j, k;

    int nn[6][3] = { {1, 0, 0},
    {0, 1, 0},
    {-1, 0, 0},
    {0, -1, 0},
    {0, 0, 1},
    {0, 0, -1}
    };
    int nn_num = 6;

    int jgap = dimx + 2;
    int kgap = (dimy + 2) * (dimx + 2);

    double hl, hfl, hbl, hfr, hr, hbr, h;
    double hu, hd, hul, hur, hdr, hdl;
    double hb, hf, huf, hub, hdb, hdf;
    double hx, hy, hxx, hyy, hxy;
    int t;
    int nor_dir = 0;
    double normal_all = 0.0;

    hx = hy = hxx = hyy = hxy = 0.0;

    for (k = 2; k <= dimz - 1; k++)
    {
        for (j = 2; j <= dimy - 1; j++)
        {
            for (i = 2; i <= dimx - 1; i++)
            {
                idx = k * kgap + j * jgap + i;
                lsp->curv[idx] = 0.0;
                if ((lsp->mold[idx]) || (lsp->gr[idx] <= 0)
                    || (fs[idx] == 0))
                {
                    continue;
                }
                int nbindex = 0;
                for (n = 0; n < nn_num; n++)
                {
                    int ndx = (k + nn[n][2]) * kgap +
                        (j + nn[n][1]) * jgap + (i + nn[n][0]);
                    nbindex += lsp->gr[ndx];
                }
                //Ensure the curvautre is only calcuated on the interface
                if (nbindex ? 1 : 0)
                {

                    //Define the index of surrounding Nbs
                    idl = idx + 1;
                    idr = idx - 1;
                    idf = idx + jgap;
                    idb = idx - jgap;
                    idu = idx + kgap;
                    idd = idx - kgap;

                    idur = idr + kgap;
                    idul = idl + kgap;
                    iddr = idr - kgap;
                    iddl = idl - kgap;

                    idfr = idf - 1;
                    idfl = idf + 1;
                    idbr = idb - 1;
                    idbl = idb + 1;

                    idub = idb + kgap;
                    iduf = idf + kgap;
                    iddf = idf - kgap;
                    iddb = idb - kgap;

                    idufr = iduf - 1;
                    idufl = iduf + 1;
                    idubr = idub - 1;
                    idubl = idub + 1;

                    iddbr = iddb - 1;
                    iddbl = iddb + 1;
                    iddfr = iddf - 1;
                    iddfl = iddf + 1;

                    fsl =
                        (fs[idufl] + fs[idubl] + fs[iddbl] +
                         fs[iddfl] + alpha * (fs[idul] +
                                                   fs[iddl] +
                                                   fs[idfl] +
                                                   fs[idbl]) +
                         beta * fs[idl]) / (1 + 4 * alpha + beta);

                    fsr =
                        (fs[idufr] + fs[idubr] + fs[iddbr] +
                         fs[iddfr] + alpha * (fs[idur] +
                                                   fs[iddr] +
                                                   fs[idfr] +
                                                   fs[idbr]) +
                         beta * fs[idr]) / (1 + 4 * alpha + beta);

                    fsf =
                        (fs[idufl] + fs[idufr] + fs[iddfr] +
                         fs[iddfl] + alpha * (fs[idfl] +
                                                   fs[idfr] +
                                                   fs[iduf] +
                                                   fs[iddf]) +
                         beta * fs[idf]) / (1 + 4 * alpha + beta);

                    fsb =
                        (fs[idubr] + fs[idubl] + fs[iddbl] +
                         fs[iddbr] + alpha * (fs[idbl] +
                                                   fs[idbr] +
                                                   fs[idub] +
                                                   fs[iddb]) +
                         beta * fs[idb]) / (1 + 4 * alpha + beta);

                    fsu =
                        (fs[idufl] + fs[idufr] + fs[idubr] +
                         fs[idubl] + alpha * (fs[idul] +
                                                   fs[idur] +
                                                   fs[idub] +
                                                   fs[iduf]) +
                         beta * fs[idl]) / (1 + 4 * alpha + beta);

                    fsd =
                        (fs[iddbr] + fs[iddbl] + fs[iddfl] +
                         fs[iddfr] + alpha * (fs[iddr] +
                                                   fs[iddl] +
                                                   fs[iddf] +
                                                   fs[iddb]) +
                         beta * fs[idd]) / (1 + 4 * alpha + beta);

                    normal_x = (fsl - fsr) / 2.0;
                    normal_y = (fsf - fsb) / 2.0;
                    normal_z = (fsu - fsd) / 2.0;

                    normal_all =
                        sqrt(normal_x * normal_x + normal_y * normal_y +
                             normal_z * normal_z);

                    if (!(normal_all))
                    {
                        normal_x = normal_x / normal_all;
                        normal_y = normal_y / normal_all;
                        normal_z = normal_z / normal_all;
                    }

                    // determine the height function direction
                    nor_dir = 0;
                    if (!(normal_x) || !(normal_y) || !(normal_z))
                    {
                        if ((ABS(normal_x) >= ABS(normal_y))
                            && (ABS(normal_x) >= ABS(normal_z)))
                            nor_dir = 1;
                        else if ((ABS(normal_y) >= ABS(normal_z))
                                 && (ABS(normal_y) >= ABS(normal_x)))
                            nor_dir = 2;
                        else if ((ABS(normal_z) >= ABS(normal_y))
                                 && (ABS(normal_z) >= ABS(normal_x)))
                            nor_dir = 3;
                    }



                    switch (nor_dir)
                    {
                            //z direction
                        case 3:

                            //determine the gamma value to smooth the interface
                            if (acos(ABS(normal_z)) < 0.8)
                                gamma = 0.0;
                            else
                                gamma = 0.2;

                            if (((k + 2) <= dimz) && (k - 2) >= 1)
                                t = 3;
                            else
                                t = 1;

                            h = hf = hb = hr = hl = hfl = hfr = hbr = hbl =
                                0.0;

                            for (n = -t; n <= t; n++)
                            {
                                h += fs[idx + kgap * n];
                                hf += fs[idf + kgap * n];
                                hb += fs[idb + kgap * n];
                                hl += fs[idl + kgap * n];
                                hr += fs[idr + kgap * n];
                                hfl += fs[idfl + kgap * n];
                                hfr += fs[idfr + kgap * n];
                                hbr += fs[idbr + kgap * n];
                                hbl += fs[idbl + kgap * n];
                            }

                            hx = (gamma * (hfl - hfr) + hl - hr +
                                  gamma * (hbl - hbr)) / (1 +
                                                          2 * gamma) / 2.0;
                            hy = (gamma * (hfl - hbl) + hf - hb +
                                  gamma * (hfr - hbr)) / (1 +
                                                          2 * gamma) / 2.0;
                            hxx =
                                (gamma * (hfl - 2 * hf + hfr) + hl - 2 * h +
                                 hr + gamma * (hbl - 2 * hb +
                                               hbr)) / (bp->cellSize) / (1 +
                                                                         2 *
                                                                         gamma);
                            hyy =
                                (gamma * (hfl - 2 * hl + hbl) + hf - 2 * h +
                                 hb + gamma * (hfr - 2 * hr +
                                               hbr)) / (bp->cellSize) / (1 +
                                                                         2 *
                                                                         gamma);
                            hxy =
                                (hfl - hbl - hfr + hbr) / bp->cellSize / 4.0;

                            //curvature for idx cell
                            lsp->curv[idx] =
                                (hxx + hyy + hxx * hy * hy + hyy * hx * hx -
                                 2 * hxy * hx * hy) / pow((1 + hx * hx +
                                                           hy * hy), 1.5);

                            break;
                            // y direction
                        case 2:
                            //determine the gamma value to smooth the interface
                            if (acos(ABS(normal_y)) < 0.8)
                                gamma = 0.0;
                            else
                                gamma = 0.2;

                            if (((j + 2) <= dimy) && (j - 2) >= 1)
                                t = 3;
                            else
                                t = 1;

                            h = hu = hd = hr = hl = hul = hur = hdr = hdl =
                                0.0;

                            for (n = -t; n <= t; n++)
                            {
                                h += fs[idx + jgap * n];
                                hu += fs[idu + jgap * n];
                                hd += fs[idd + jgap * n];
                                hl += fs[idl + jgap * n];
                                hr += fs[idr + jgap * n];
                                hul += fs[idul + jgap * n];
                                hur += fs[idur + jgap * n];
                                hdr += fs[iddr + jgap * n];
                                hdl += fs[iddl + jgap * n];
                            }

                            hx = (gamma * (hul - hur) + hl - hr +
                                  gamma * (hdl - hdr)) / (1 +
                                                          2 * gamma) / 2.0;
                            hy = (gamma * (hul - hdl) + hu - hd +
                                  gamma * (hur - hdr)) / (1 +
                                                          2 * gamma) / 2.0;
                            hxx =
                                (gamma * (hul - 2 * hu + hur) + hl - 2 * h +
                                 hr + gamma * (hdl - 2 * hd +
                                               hdr)) / (bp->cellSize) / (1 +
                                                                         2 *
                                                                         gamma);
                            hyy =
                                (gamma * (hul - 2 * hl + hdl) + hu - 2 * h +
                                 hd + gamma * (hur - 2 * hr +
                                               hdr)) / (bp->cellSize) / (1 +
                                                                         2 *
                                                                         gamma);
                            hxy =
                                (hul - hdl - hur + hdr) / bp->cellSize / 4.0;

                            //curvature for idx cell
                            lsp->curv[idx] =
                                (hxx + hyy + hxx * hy * hy + hyy * hx * hx -
                                 2 * hxy * hx * hy) / pow((1 + hx * hx +
                                                           hy * hy), 1.5);

                            break;


                            // x direction
                        case 1:
                            //determine the gamma value to smooth the interface
                            if (acos(ABS(normal_x)) < 0.8)
                                gamma = 0.0;
                            else
                                gamma = 0.2;

                            if (((i + 2) <= dimx) && (i - 2) >= 1)
                                t = 3;
                            else
                                t = 1;

                            h = hu = hd = hb = hf = huf = hub = hdb = hdf =
                                0.0;

                            for (n = -t; n <= t; n++)
                            {
                                h += fs[idx + n];
                                hu += fs[idu + n];
                                hd += fs[idd + n];
                                hf += fs[idf + n];
                                hb += fs[idb + n];
                                huf += fs[iduf + n];
                                hub += fs[idub + n];
                                hdb += fs[iddb + n];
                                hdf += fs[iddf + n];
                            }

                            hx = (gamma * (huf - hub) + hf - hb +
                                  gamma * (hdf - hdb)) / (1 +
                                                          2 * gamma) / 2.0;
                            hy = (gamma * (huf - hdf) + hu - hd +
                                  gamma * (hub - hdb)) / (1 +
                                                          2 * gamma) / 2.0;
                            hxx =
                                (gamma * (huf - 2 * hu + hub) + hf - 2 * h +
                                 hb + gamma * (hdf - 2 * hd +
                                               hdb)) / (bp->cellSize) / (1 +
                                                                         2 *
                                                                         gamma);
                            hyy =
                                (gamma * (huf - 2 * hf + hdf) + hu - 2 * h +
                                 hd + gamma * (hub - 2 * hb +
                                               hdb)) / (bp->cellSize) / (1 +
                                                                         2 *
                                                                         gamma);
                            hxy =
                                (huf - hdf - hub + hdb) / bp->cellSize / 4.0;

                            //curvature for idx cell
                            lsp->curv[idx] =
                                (hxx + hyy + hxx * hy * hy + hyy * hx * hx -
                                 2 * hxy * hx * hy) / pow((1 + hx * hx +
                                                           hy * hy), 1.5);

                            break;

                        default:
                            //      printf("Curvature: No value\n");
                            break;
                    }
                }
            }
        }
    }


}
