#ifndef SSH3_CLI_H
#define SSH3_CLI_H

#include "pivot_table_file.h"
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

static int ssh3_cli_number(const char *text,uint64_t *value) {
    char *end; unsigned long long x;
    if (!text || *text<'0' || *text>'9') return -1;
    errno=0; x=strtoull(text,&end,10);
    if (errno || *end) return -1;
    *value=(uint64_t)x; return 0;
}
static int ssh3_cli_alias(const char *a,const char *b) {
    struct stat x,y;
    return !strcmp(a,b) || (!stat(a,&x) && !stat(b,&y) &&
        x.st_dev==y.st_dev && x.st_ino==y.st_ino);
}
/* The caller owns generated and the opened table file, including on failure. */
static int ssh3_cli_tables(SubseqHash3PivotTableFile *file,const char *path,
    char **generated,const char *base,const char *input,size_t w,size_t d,
    size_t sets,size_t required_rows,uint64_t seed) {
    if (!path) {
        const size_t len=strlen(base);
        if (len>SIZE_MAX-5) return -1;
        *generated=malloc(len+5); if (!*generated) return -1;
        memcpy(*generated,base,len); memcpy(*generated+len,".piv",5);
        if (ssh3_cli_alias(input,*generated) || ssh3_cli_alias(base,*generated)) {
            fputs("generated table path aliases INPUT or OUTPUT\n",stderr); return -1;
        }
        if (subseqhash3_seed_tables_create(*generated,w,d,sets,seed)) return -1;
        path=*generated;
        fprintf(stderr,"pivot table seed: %llu\n",(unsigned long long)seed);
    }
    if (subseqhash3_pivot_table_file_open(file,path)) {
        fputs("expected a pivot table file from 'tables W D S OUTPUT'\n",stderr); return -1;
    }
    if (file->d!=d || file->count<sets || file->n-1<required_rows) {
        fputs("table file must have matching D, at least S sets and enough rows for the seed length\n",stderr);
        return -1;
    }
    return 0;
}
#endif
