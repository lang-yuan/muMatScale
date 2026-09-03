/***************************************************************/
/* Copyright (c) 2023, Lang Yuan, Univeristy of South Carolina */
/* All rights reserved.                                        */
/* This file is part of muMatScale.                            */
/* See the top-level LICENSE file for details.                 */
/***************************************************************/

#ifndef __VARIABLES_H__
#define __VARIABLES_H__

#include <mpi.h>

int registerCommInfo(size_t datasize);
void FinishExchangeForVar(int variable_key, void* data);
void ExchangeFacesForVar(int variable_key, void* d);

#endif
