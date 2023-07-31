/***************************************************************/
/* Copyright (c) 2023, Lang Yuan, Univeristy of South Carolina */
/* All rights reserved.                                        */
/* This file is part of muMatScale.                            */
/* See the top-level LICENSE file for details.                 */
/***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <stddef.h>
#include <omp.h>
#include "face_util.h"
#include "globals.h"
#include "functions.h"
#include "xmalloc.h"
#include "debug.h"
#include "profiler.h"
#include "packing.h"

/**
 * Initializes the MPI Comm datatypes needed for communication
 * for an array of a particular base type.
 *
 * \param[in] origType The Base MPI Datatype for this type.
 * \param[in,out] types A pointer to an existing \c commTypes_t datastructure.
 *                      that will be filled out.
 *
 * \return 0 on success, otherwise will error
 */
int
createCommTypes(
    MPI_Datatype origType,
    commTypes_t * types)
{
    if (origType == MPI_DOUBLE)
        types->dsize = sizeof(double);
    if (origType == MPI_INT)
        types->dsize = sizeof(int);
    if (origType == bp->MPI_DECENTERED)
        types->dsize = sizeof(decentered_t);

    return 0;
}

#define MKTAG(sbid, data_key, face) \
	(((sbid)<<(TAG_FACE_SHIFT+TAG_DATA_KEY_SHIFT)) | \
	 ((data_key)<<(TAG_FACE_SHIFT)) | \
	 (face))

/**
 * Receives halo information
 *
 * \param[in]  data    The data buffer to receive into
 * \param[in]  t       The types datastructure
 * \param[in]  halo    Which halo to recive into
 * \param[in]  from    MPI Rank of sender
 * \param[in]  tag     MPI Tag
 * \param[in]  rbuffer MPI buffer
 * \param[out] req     MPI_Request to fill
 * \return MPI_SUCCESS on no error, otherwise, an error code
 */
static int
Recv_Plane(
    void *data,
    commTypes_t * t,
    int halo,
    int from,
    int tag,
    void *rbuffer[6],
    MPI_Request * req)
{
    assert(rbuffer[halo] != NULL);
    int offset;
    int stride;
    int bsize;
    int nblocks;
    computeHaloInfo(halo, &offset, &stride, &bsize, &nblocks);

    int sizeb = bsize * nblocks * t->dsize;
    return MPI_Irecv(rbuffer[halo], sizeb, MPI_BYTE, from, tag, mpi_comm_new,
                     req);
}

void
unpack_plane(
    void *data,
    commTypes_t * t,
    int halo,
    void *rbuffer[6])
{
    assert(rbuffer[halo] != NULL);

    int offset;
    int stride;
    int bsize;
    int nblocks;
    computeHaloInfo(halo, &offset, &stride, &bsize, &nblocks);

    unpack_field(t->dsize, data, stride, bsize, nblocks, offset,
                 rbuffer[halo]);
}

/**
 * Receives halo information, Non-Blocking
 *
 * \param[in] data     The data buffer to receive into
 * \param[in] data_key user-supplied "key" to mix into tag
 * \param[in] t        Datatypes structure to be used
 * \param[in] connMap  Connection map
 * \param[in] rbuf     Recv buffer
 * \param[in,out] reqs MPI Requests array, will have requests added to it
 * \return number of MPI_Requests added to array
 */
static int
RecvHalosNB(
    void *data,
    int data_key,
    commTypes_t * t,
    int connMap[][2],
    void *rbuf[6],
    MPI_Request * reqs)
{
    int n_req = 0;
    int err = 0;

    for (int face = 0; face < NUM_NEIGHBORS; face++)
    {
        int rank = connMap[face][0];
        // A rank of less than 0 means that it isn't assigned
        if (rank >= 0 && rank != iproc)
        {
            dwrite(DEBUG_MPI,
                   "Posting recv from %d with tag %d from face %d\n", rank,
                   MKTAG(0, data_key, gface_reverse[face]), face);
            timing(COMPUTATION, timer_elapsed());

            err = Recv_Plane(data, t, face, rank,
                             MKTAG(0, data_key,
                                   gface_reverse[face]), rbuf,
                             &reqs[n_req++]);
            timing(COMMUNICATION, timer_elapsed());
            if (err)
                goto err;
        }
    }

    return n_req;
  err:
    error("Received Recv Error: %d\n", err);
}

/**
 * Sends a face of a 3D array to a target MPI task
 *
 * \param[in]  data   The data to send
 * \param[in]  t      The types datastructure
 * \param[in]  face   Which face to send from
 * \param[in]  to     MPI Rank of receiver
 * \param[in]  tag    MPI Tag
 * \param[in]  sbuffer The buffers to send data
 * \param[out] req    MPI_Request to fill
 * \return MPI_SUCCESS on no error, otherwise, an error code
 */
static int
Send_Plane(
    void *data,
    commTypes_t * t,
    int face,
    int to,
    int tag,
    void *sbuffer[6],
    MPI_Request * req)
{
    assert(sbuffer[face] != NULL);
    int offset;
    int stride;
    int bsize;
    int nblocks;
    computeFaceInfo(face, &offset, &stride, &bsize, &nblocks);

    pack_field(t->dsize, data, stride, bsize, nblocks, offset, sbuffer[face]);

    int sizeb = bsize * nblocks * t->dsize;
    return MPI_Isend(sbuffer[face], sizeb, MPI_BYTE, to, tag, mpi_comm_new,
                     req);
}

/**
 * Sends faces to a halo, Non-Blocking
 *
 * \param[in] data     Array to send
 * \param[in] data_key user-supplied "key" to mix into tag
 * \param[in] t        Datatypes structure to be used
 * \param[in] connMap  Connection map
 * \param[in]  sbuf    The buffer to send data
 * \param[in,out] reqs MPI Requests array, will have requests added to it
 * \return number of MPI_Requests added to array
 */
int
SendFacesNB(
    void *data,
    int data_key,
    commTypes_t * t,
    int connMap[][2],
    void *sbuf[6],
    MPI_Request * reqs)
{
    int n_req = 0;
    int err = 0;

    for (int face = 0; face < NUM_NEIGHBORS; face++)
    {
        int rank = connMap[face][0];
        // A rank of less than 0 means that it isn't assigned
        if (rank >= 0 && rank != iproc)
        {
            dwrite(DEBUG_MPI, "Posting send to %d with tag %d of face %d\n",
                   rank, MKTAG(0, data_key, face), face);
            timing(COMPUTATION, timer_elapsed());

            err = Send_Plane(data, t, face, rank,
                             MKTAG(0, data_key, face), sbuf, &reqs[n_req++]);
            timing(COMMUNICATION, timer_elapsed());
            if (err)
                goto err;
        }
    }

    return n_req;
  err:
    error("Received Send Error: %d\n", err);
}

/**
 * Send/Recv the halo information for an array in a Non-Blocking fashion
 *
 * \param[in] data ter to array to be updated
 * \param[in] data_key an identifier for this array. {User specified}
 * \param[in] t The communication datatypes needed for this array
 * \param[in] connMap Neighbor connection map
 * \param[in] sbuf Buffer for sending data
 * \param[in] rbuf Buffer to receive data
 * \param[in,out] reqs pointer to array of MPI_Requests that will be filled.
 *
 * \return number of MPI Requests in the request buffer
 */
int
SendRecvHalosNB(
    void *data,
    int data_key,
    commTypes_t * t,
    int connMap[][2],
    void *sbuf[6],
    void *rbuf[6],
    MPI_Request * reqs)
{
    int n_req = 0;
    n_req += RecvHalosNB(data, data_key, t, connMap, rbuf, &reqs[0]);
    n_req += SendFacesNB(data, data_key, t, connMap, sbuf, &reqs[n_req]);

    return n_req;
}
