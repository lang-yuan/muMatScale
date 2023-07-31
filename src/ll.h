/***************************************************************/
/* Copyright (c) 2023, Lang Yuan, Univeristy of South Carolina */
/* All rights reserved.                                        */
/* This file is part of muMatScale.                            */
/* See the top-level LICENSE file for details.                 */
/***************************************************************/
#ifndef LL_H_
#define LL_H_

#include "globals.h"

struct lli_s
{
    struct lli_s *next;
    struct lli_s *prev;
    void *data;
};
typedef struct lli_s lli_t;

struct ll_s
{
    lli_t *head;
    lli_t *tail;
    int count;
    int (
    *cmp) (
    const void *v1,
    const void *v2,
    const void *v3);
    const void *ud;
};
typedef struct ll_s ll_t;


ll_t *ll_init(
    int (*cmp) (const void *v1,
                const void *v2,
                const void *v3),
    const void *ud);            /* Initialize LL  */
int ll_count(
    ll_t * ll);
int ll_insert(
    ll_t * ll,
    void *data);                /* Insert at Head */
int ll_append(
    ll_t * ll,
    void *data);                /* Insert at Tail */
int ll_delete(
    ll_t * ll,
    void *data);                /* Delete Element */
int lli_delete(
    ll_t * ll,
    lli_t * lli);               /* Delete Element */
void *ll_head(
    ll_t * ll);                 /* returns first element */
void *ll_tail(
    ll_t * ll);                 /* returns last element */
void *ll_find(
    ll_t * ll,
    void *data);                /* find an element in the list */
void *ll_pop(
    ll_t * ll);                 /* Pops the first element off */
void ll_walk(
    ll_t * ll,
    void (*func) (SB_struct * data,
                  void *ud),
    void *ud);
void ll_destroy(
    ll_t * ll);
void ll_print(
    ll_t * ll);

#endif /* LL_H_ */
