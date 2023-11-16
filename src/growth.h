/***************************************************************/
/* Copyright (c) 2023, Lang Yuan, Univeristy of South Carolina */
/* All rights reserved.                                        */
/* This file is part of muMatScale.                            */
/* See the top-level LICENSE file for details.                 */
/***************************************************************/

#ifndef CALC_DIFFUSION_H_
#define CALC_DIFFUSION_H_

void sb_diffuse_alloy_decentered(
    SB_struct * lsp,
    void * __attribute__ ((__unused__)) __unused);
void fs_change_diffuse(
    SB_struct * lsp,
    void * __attribute__ ((__unused__)) __unused);
void capture_octahedra_diffuse(
    SB_struct * lsp,
    void *vsolid_volume);
void grow_octahedra(
    SB_struct * lsp,
    void * __attribute__ ((__unused__)) __unused);
void grow_cell_reduction(
    SB_struct * lsp,
    void * __attribute__ ((__unused__)) __unused);
void pre_cell_reduction(
    SB_struct * lsp,
    void * __attribute__ ((__unused__)) __unused);

void gr_dataexchange_to(
    SB_struct * lsp,
    void * __attribute__ ((__unused__)) __unused);
void gr_dataexchange_from(
    SB_struct * lsp,
    void * __attribute__ ((__unused__)) __unused);
#endif
