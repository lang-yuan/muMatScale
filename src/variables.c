#include "variables.h"
#include "face_util.h"

#include <assert.h>

int var_count = 0;
variable_registration var_regs[1 << TAG_DATA_KEY_SHIFT];

int
registerCommInfo(
    size_t datasize)
{
    assert(var_count < (1 << TAG_DATA_KEY_SHIFT));
    int varnum = var_count++;

    variable_registration *v = &var_regs[varnum];
    v->reqs = NULL;
    v->nreq = 0;
    v->datasize = datasize;
    v->rbuf_base = NULL;
    v->sbuf_base = NULL;
    v->buffer_slot_cells = 0;
    v->buffer_slot_bytes = 0;
    for (int i = 0; i < 6; i++)
    {
        v->rbuf[i] = NULL;
        v->sbuf[i] = NULL;
    }
    return varnum;
}

