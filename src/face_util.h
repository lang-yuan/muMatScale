/***************************************************************/
/* Copyright (c) 2023, Lang Yuan, Univeristy of South Carolina */
/* All rights reserved.                                        */
/* This file is part of muMatScale.                            */
/* See the top-level LICENSE file for details.                 */
/***************************************************************/

#ifndef FACE_UTIL_H_
#define FACE_UTIL_H_
#include <mpi.h>
#include "globals.h"


#define TAG_FACE_SHIFT 3        /* Support 6 faces, need 3 bits */
#define TAG_DATA_KEY_SHIFT 3    /* Support 8 data keys (additional tag info) */


/**
 * Initializes the MPI Comm datatypes needed for communication
 * for an array of a particular base type.
 *
 * \param[in] origType The Base MPI Datatype for this type.
 * \param[in,out] types A pointer to an existing \c commTypes_t datastructure.
 *                      that will be filled out.
 * \param[in] nz Size of the array in Z
 * \param[in] ny Size of the array in Y
 * \param[in] nx Size of the array in X
 *
 * \return 0 on success, otherwise will error
 */
int createCommTypes(
    MPI_Datatype origType,
    commTypes_t * types);
/**
 * Destroys / Deallocates the MPI datatypes used for face communication
 * of the passed variable
 *
 * \param[in] types pointer to structure with types to be deallocated
 */
void destroyCommTypes(
    commTypes_t * types);



/**
 * Send/Recv the halo information for an array in a Non-Blocking fashion
 *
 * \param[in] data ter to array to be updated
 * \param[in] data_key an identifier for this array. {User specified}
 * \param[in] t The communication datatypes needed for this array
 * \param[in] myID  ID of the subblock being updated
 * \param[in] connMap Neighbor connection map
 * \param[in] sbuf Buffer for sending data
 * \param[in] rbuf Buffer to receive data
 * \param[in,out] reqs pointer to array of MPI_Requests that will be filled.
 *
 * \return number of MPI Requests in the request buffer
 */
int SendRecvHalosNB(
    void *data,
    int data_key,
    size_t datasize,
    int connMap[][2],
    void *sbuf[6],
    void *rbuf[6],
    MPI_Request * reqs);

void unpack_plane(
    void *data,
    size_t datasize,
    int halo,
    void *rbuf[6]);

#endif /* FACE_UTIL.H */
