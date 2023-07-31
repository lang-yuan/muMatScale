/***************************************************************/
/* Copyright (c) 2023, Lang Yuan, Univeristy of South Carolina */
/* All rights reserved.                                        */
/* This file is part of muMatScale.                            */
/* See the top-level LICENSE file for details.                 */
/***************************************************************/


#ifndef __READ_CTRL_H__
#define __READ_CTRL_H__


void set_defaults(
    void);
void verify_config(
    void);
void read_config(
    const char *ctrl_fname);
void print_config(
    FILE * fp);



#endif
