#ifndef SSH3_SEEDING_H
#define SSH3_SEEDING_H

#include "ssh3_deletions.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Pivot-based multiple-seed interface. Tables contains S independent sets;
   R leading pivots are used from EACH set. Each callback receives R*S seeds
   in set-major, then pivot-major order (index = set*R+pivot).
   OUTWARD=0 is the default traversal; INWARD selects the alternate order.
   Positions are relative to the current window. Callback results are borrowed
   until it returns. Return 0 on success, -1 for invalid arguments, allocation/
   seeding failure or nonzero callback return. Input symbols must be 0..3;
   length>=window, R>0, S>0. Tables stay immutable throughout the call.
   Passing length==window produces one batch for a single window.
   Outward reuses anchor DPs across overlapping windows for every set. */
typedef SubseqHash3RepeatsWindowCallback SubseqHash3SubSeedsCallback;
typedef SubseqHash3DeletionWindowCallback SubseqHash3IndSeedsCallback;

/* One deletion, exactly E substitutions. 2<=W<=257, E<=W-1, R<=W-1. */
int subseqhash3_seeds_sub(const uint8_t *sequence,size_t length,size_t window,
    size_t substitutions,const SubseqHash3PivotTables *tables,
    size_t pivots_per_set,size_t sets,SubseqHash3Traversal traversal,
    SubseqHash3SubSeedsCallback callback,void *context);

/* Exactly E deletions, no substitutions. 1<=W<=256, E<W, R<=W-E.
   Accepts the SAME pivot tables as seeds_sub; conversion is once per call. */
int subseqhash3_seeds_ind(const uint8_t *sequence,size_t length,size_t window,
    size_t deletions,const SubseqHash3PivotTables *tables,
    size_t pivots_per_set,size_t sets,SubseqHash3Traversal traversal,
    SubseqHash3IndSeedsCallback callback,void *context);

#ifdef __cplusplus
}
#endif
#endif
