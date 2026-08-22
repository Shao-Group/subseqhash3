#include "subseqhash_edits_i8.h"

#include <assert.h>
#include <stdio.h>
#include <string.h>

int main(void) {
    size_t trial;
    for (trial=0;trial<100;++trial) {
        const size_t n=3+trial%5,edits=trial%3,d=3+trial%13;
        uint8_t s[7]; size_t i,split;
        SubseqHashI8Tables compact;
        SubseqHashTables reference;
        SubseqHashEditResult actual[3],legacy[3],twopass[3],expected[3];
        for (i=0;i<n;++i) s[i]=(uint8_t)((trial*11+i*3+i*i)%4);
        assert(!subseqhash_i8_tables_init(&compact,n,d,12345+trial));
        assert(!subseqhash_i8_tables_make_double(&compact,&reference));
        assert(!subseqhash_compute_all_edit_splits_i8(
            s,n,edits,&compact,actual));
        assert(!subseqhash_compute_all_edit_splits_i8_legacy(
            s,n,edits,&compact,legacy));
        assert(!subseqhash_compute_all_edit_splits_i8_twopass(
            s,n,edits,&compact,twopass));
        assert(!subseqhash_compute_all_edit_splits(
            s,n,edits,&reference,expected));
        for (split=0;split<=edits;++split) {
            const size_t k=actual[split].length;
            assert(actual[split].deletions==expected[split].deletions);
            assert(actual[split].substitutions==expected[split].substitutions);
            assert(actual[split].psi==expected[split].psi);
            assert(actual[split].omega==expected[split].omega);
            assert(!memcmp(actual[split].seed,expected[split].seed,k));
            assert(!memcmp(actual[split].positions,expected[split].positions,
                           k*sizeof(*actual[split].positions)));
            assert(!memcmp(actual[split].substituted,
                           expected[split].substituted,k));
            assert(actual[split].psi==legacy[split].psi);
            assert(actual[split].omega==legacy[split].omega);
            assert(!memcmp(actual[split].seed,legacy[split].seed,k));
            assert(!memcmp(actual[split].positions,legacy[split].positions,
                           k*sizeof(*actual[split].positions)));
            assert(!memcmp(actual[split].substituted,
                           legacy[split].substituted,k));
            assert(actual[split].psi==twopass[split].psi);
            assert(actual[split].omega==twopass[split].omega);
            assert(!memcmp(actual[split].seed,twopass[split].seed,k));
            assert(!memcmp(actual[split].positions,twopass[split].positions,
                           k*sizeof(*actual[split].positions)));
            assert(!memcmp(actual[split].substituted,
                           twopass[split].substituted,k));
        }
        subseqhash_edit_results_destroy(expected,edits+1);
        subseqhash_edit_results_destroy(twopass,edits+1);
        subseqhash_edit_results_destroy(legacy,edits+1);
        subseqhash_edit_results_destroy(actual,edits+1);
        subseqhash_tables_destroy(&reference);
        subseqhash_i8_tables_destroy(&compact);
    }
    printf("all i8 edit tests passed (%s)\n",
           subseqhash_edits_i8_backend());
    return 0;
}
