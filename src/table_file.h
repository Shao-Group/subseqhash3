#ifndef SUBSEQHASH3_TABLE_FILE_H
#define SUBSEQHASH3_TABLE_FILE_H

#include "subseqhash3.h"

typedef struct {
    void *mapping;
    size_t mapping_size,n,d,count;
    SubseqHash3Tables *tables;
} SubseqHash3TableFile;

int subseqhash3_table_file_create(const char *path,size_t n,size_t d,
                                  size_t count,uint64_t seed);
int subseqhash3_table_file_open(SubseqHash3TableFile *file,const char *path);
void subseqhash3_table_file_close(SubseqHash3TableFile *file);

#endif
