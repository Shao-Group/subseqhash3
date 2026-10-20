#ifndef SSH3_PIVOT_TABLE_FILE_H
#define SSH3_PIVOT_TABLE_FILE_H

#include "ssh3_pivots.h"
#ifdef __cplusplus
extern "C" {
#endif

/* Shared on-disk format for repeats and deletions. n-1 is row capacity;
   repeats needs n>=W, deletions needs n-1>=W-E (and uses family zero).
   Owned, writable arrays, unlike the legacy mmap-backed table-file API.
   count is the number of independent families, not the number of seeds. */
typedef struct {
    size_t n,d,count;
    SubseqHash3PivotTables *tables;
} SubseqHash3PivotTableFile;

int subseqhash3_pivot_table_file_create(const char *path,size_t n,size_t d,
    size_t count,uint64_t seed);
int subseqhash3_pivot_table_file_open(SubseqHash3PivotTableFile *file,const char *path);
void subseqhash3_pivot_table_file_close(SubseqHash3PivotTableFile *file);

#ifdef __cplusplus
}
#endif
#endif
