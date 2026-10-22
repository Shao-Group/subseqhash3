#ifndef SSH3_PIVOT_TABLE_FILE_H
#define SSH3_PIVOT_TABLE_FILE_H

#include "ssh3_pivots.h"
#ifdef __cplusplus
extern "C" {
#endif

/* Shared on-disk format for pivot seeding. n-1 is row capacity;
   seeds-sub needs n>=W, seeds-ind needs n-1>=W-E. Both use S sets.
   Owned, writable arrays, unlike the legacy mmap-backed table-file API.
   count is the number of independent families, not the number of seeds. */
typedef struct {
    size_t n,d,count;
    SubseqHash3PivotTables *tables;
} SubseqHash3PivotTableFile;

/* Preferred interface: W sizes one set for both seed modes, including E=0
   indels. Stores min(W,256) rows per set; 1<=W<=257, D>0, sets>0. */
int subseqhash3_seed_tables_create(const char *path,size_t window,size_t residues,
    size_t sets,uint64_t seed);

/* Low-level row-capacity interface retained for existing C callers. */
int subseqhash3_pivot_table_file_create(const char *path,size_t n,size_t d,
    size_t count,uint64_t seed);
int subseqhash3_pivot_table_file_open(SubseqHash3PivotTableFile *file,const char *path);
void subseqhash3_pivot_table_file_close(SubseqHash3PivotTableFile *file);

#ifdef __cplusplus
}
#endif
#endif
