/***************************************************************/
/* Copyright (c) 2023, Lang Yuan, Univeristy of South Carolina */
/* All rights reserved.                                        */
/* This file is part of muMatScale.                            */
/* See the top-level LICENSE file for details.                 */
/***************************************************************/

#include <errno.h>
#include "xmalloc.h"
#include "ll.h"

static int ll_cmp(
    const void *v1,
    const void *v2,
    const void *v3);
static lli_t *lli_find(
    ll_t * ll,
    void *data);
static void lli_del(
    lli_t * lli);

ll_t *
ll_init(
    int (*cmp) (const void *v1,
                const void *v2,
                const void *v3),
    const void *ud)
{
    ll_t *ll;

    xmalloc(ll, ll_t, 1);

    ll->head = ll->tail = NULL;
    ll->count = 0;
    ll->cmp = cmp ? cmp : ll_cmp;
    ll->ud = ud;

    return ll;
}


int
ll_count(
    ll_t * ll)
{
    if (!ll)
    {
        errno = EINVAL;
        return -1;
    }
    else
    {
        return ll->count;
    }
}


int
ll_insert(
    ll_t * ll,
    void *data)
{
    lli_t *lli = NULL;
    if (!ll)
    {
        errno = EINVAL;
        return 1;
    }

    xmalloc(lli, lli_t, 1);

    lli->prev = NULL;
    lli->next = ll->head;
    lli->data = data;
    if (lli->next)
    {
        lli->next->prev = lli;
    }
    else
    {
        /* ll->head was NULL */
        ll->tail = lli;
    }
    ll->head = lli;
    ll->count++;

    return 0;
}

int
ll_append(
    ll_t * ll,
    void *data)
{
    lli_t *lli = NULL;
    if (!ll)
    {
        errno = EINVAL;
        return EINVAL;
    }

    xmalloc(lli, lli_t, 1);

    lli->next = NULL;
    lli->prev = ll->tail;
    lli->data = data;
    if (lli->prev)
    {
        lli->prev->next = lli;
    }
    else
    {
        /* ll->tail was null */
        ll->head = lli;
    }
    ll->tail = lli;

    ll->count++;

    return 0;
}

int
ll_delete(
    ll_t * ll,
    void *data)
{
    lli_t *lli = NULL;
    if (!ll)
    {
        errno = EINVAL;
        return 1;
    }

    lli = lli_find(ll, data);
    if (!lli)
    {
        errno = ENOENT;
        return 1;
    }

    return lli_delete(ll, lli);
}


int
lli_delete(
    ll_t * ll,
    lli_t * lli)
{
    if (!ll || !lli)
    {
        errno = EINVAL;
        return 1;
    }

    if (lli->prev)
    {
        lli->prev->next = lli->next;
    }
    else
    {                           /* Head Node being removed */
        ll->head = lli->next;
    }
    if (lli->next)
    {
        lli->next->prev = lli->prev;
    }
    else
    {                           /* Tail Node being removed */
        ll->tail = lli->prev;
    }

    xfree(lli);

    ll->count--;

    return 0;
}


void *
ll_head(
    ll_t * ll)
{
    if (!ll)
    {
        errno = EINVAL;
        return NULL;
    }

    if (ll->head)
    {
        return ll->head->data;
    }
    else
    {
        errno = 0;
        return NULL;
    }
}


void *
ll_tail(
    ll_t * ll)
{
    if (!ll)
    {
        errno = EINVAL;
        return NULL;
    }

    if (ll->tail)
    {
        return ll->tail->data;
    }
    else
    {
        return NULL;
    }
}


void *
ll_find(
    ll_t * ll,
    void *data)
{
    lli_t *lli = lli_find(ll, data);
    return lli ? lli->data : NULL;
}

void *
ll_pop(
    ll_t * ll)
{
    void *data = NULL;

    if (!ll)
    {
        errno = EINVAL;
        return NULL;
    }

    data = ll_head(ll);
    ll_delete(ll, data);
    return data;
}


void
ll_walk(
    ll_t * ll,
    void (*func) (SB_struct * data,
                  void *ud),
    void *ud)
{
    lli_t *lli = NULL;
    if (!ll || !func)
    {
        errno = EINVAL;
        return;
    }

    lli = ll->head;
    while (lli)
    {
        SB_struct *lsp = (SB_struct *) (lli->data);
        (*func) (lsp, ud);

        lli = lli->next;
    }
}


void
ll_destroy(
    ll_t * ll)
{
    if (!ll)
    {
        errno = EINVAL;
        return;
    }

    lli_del(ll->head);
    ll->head = ll->tail = NULL;

    xfree(ll);
}


void
ll_print(
    ll_t * ll)
{
    lli_t *lli;
    if (!ll)
    {
        errno = EINVAL;
        return;
    }

    lli = ll->head;
    //fprintf(stderr,"\nPrinting List: %p\n", (void*)ll );
    while (lli)
    {
        //fprintf(stderr,"%p (%p) ->\n", (void*)lli, lli->data);
        lli = lli->next;
    }
    //fprintf(stderr, "(0x0) [tail]\n"); fflush(stderr);
}


static int
ll_cmp(
    const void *v1,
    const void *v2,
    const void *v3)
{
    return !(v1 == v2);
}

static lli_t *
lli_find(
    ll_t * ll,
    void *data)
{
    lli_t *lli = NULL;

    if (!ll)
    {
        errno = EINVAL;
        return NULL;
    }

    lli = ll->head;

    while (lli && (*(ll->cmp)) (lli->data, data, ll->ud))
    {
        lli = lli->next;
    }

    return lli;
}

static void
lli_del(
    lli_t * lli)
{
    if (!lli)
    {
        errno = EINVAL;
        return;
    }

    lli_del(lli->next);
    xfree(lli);
}
