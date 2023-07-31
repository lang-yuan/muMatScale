/***************************************************************/
/* Copyright (c) 2023, Lang Yuan, Univeristy of South Carolina */
/* All rights reserved.                                        */
/* This file is part of muMatScale.                            */
/* See the top-level LICENSE file for details.                 */
/***************************************************************/
/*
 * Debugging routines
 *  - Defines:  dprintf(), alert() and error()
 */

#ifndef DEBUG_H_
#define DEBUG_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "globals.h"

#ifndef EBUG
#define EBUG 0
#endif

#define DEBUG_CHECKPOINT 0

#define DEBUG_MPI         (1<<0)
#define DEBUG_MAIN_CTRL (1<<1)
#define DEBUG_TASK_CTRL (1<<2)
#define DEBUG_PARSING     (1<<3)
#define DEBUG_VERIFY_CORRECT     (1<<4)

#ifdef MPI_VERSION
#define ABORT MPI_Abort(mpi_comm_new, 1); abort()
#else
#define ABORT abort()
#endif


#define __display(type, fmt, ...) \
		fprintf(stderr,  "%3d %s: [{%s} %s:%d] " fmt, \
				iproc, type, __func__, __FILE__, __LINE__, ##__VA_ARGS__ ); \
		fflush(stderr);


#define alert(fmt, ...) \
	do { \
		__display("ALERT", fmt, ##__VA_ARGS__); \
	} while (0)

#define derror(str) error(str ": %s\n", strerror(errno))
#define werror(str) alert(str ": %s\n", strerror(errno))

#define error(fmt, ...) \
	do { \
		__display("ERROR", fmt, ##__VA_ARGS__); \
		ABORT; \
	} while (0)

#define assert(a) if (!(a)) error("Assertion (" #a ") failed!\n")

#if !(EBUG)
#define dprintf(fmt, ...) do { } while(0)
#else
#define dprintf(fmt, ...) \
	do { \
		__display("DEBUG", fmt, ##__VA_ARGS__); \
	} while(0)
#endif

#ifndef NODEBUG
#define dwrite(test, ...) \
	do { \
		if ( test & EBUG ) { \
			__display("DEBUG", ##__VA_ARGS__); \
		} \
	} while(0)
#else
#define dwrite(test, ...) do { } while(0)
#endif

#endif /* DEBUG_H_ */
