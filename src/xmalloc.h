/***************************************************************/
/* Copyright (c) 2023, Lang Yuan, Univeristy of South Carolina */
/* All rights reserved.                                        */
/* This file is part of muMatScale.                            */
/* See the top-level LICENSE file for details.                 */
/***************************************************************/
#ifndef XMALLOC_H_
#define XMALLOC_H_


#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include "debug.h"



/**
 * A "safe" \c malloc that will give an error and exit if \c malloc fails.
 * Will also set the newly allocated buffer to zeros.
 * \param[in] var A pointer variable to which the buffer should be assigned.
 * \param[in] type The type of each element of the buffer.
 * \param[in] count How many elements should be allocated in the buffer.
 */
#define xmalloc(var, type, count) \
	do { \
		size_t _size = (count) * sizeof(type); \
		if ( ( var = (type *)malloc( _size ) ) == NULL ) { \
			derror("malloc"); \
		} else { \
			memset( var, 0, _size ); \
		} \
	} while ( 0 )

/**
 * A "safe" \c realloc that will give an error if it fails
 * \param[in] var A pointer variable to which the buffer should be assigned.
 * \param[in] type The type of each element of the buffer.
 * \param[in] count How many elements should be reallocated in the buffer.
 */
#define xrealloc(var, type, count) \
	do { \
		type *_tmp = NULL; \
		size_t _size = (count) * sizeof(type);  \
		_tmp = (type *)realloc( var, _size ); \
		if ( _tmp == NULL ) { \
			derror("realloc"); \
		} else { \
			var = _tmp; \
		} \
	} while ( 0 )
/**
 * A "safe" \c free call that will only free if the variable is not \c NULL.
 * Will set the variable to \c NULL.
 */
#define xfree( var ) \
	do { \
		if ( var != NULL ) { \
			free( var ); \
			var = NULL; \
		} else { \
			alert("Double Free on " #var "!\n"); \
		} \
	} while ( 0 )

#endif /* XMALLOC_H_ */
