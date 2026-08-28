#include "subseqhash_edits.h"

#include <stdlib.h>
#include <string.h>

void subseqhash_edit_results_destroy(SubseqHashEditResult *results,
                                     size_t count) {
    size_t i;
    if (!results) return;
    for (i=0;i<count;++i) {
        free(results[i].substituted);
        free(results[i].positions);
        free(results[i].seed);
        memset(results+i,0,sizeof(*results));
    }
}
