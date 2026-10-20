#define _POSIX_C_SOURCE 200809L
#include "pivot_table_file.h"

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

/* Version 1: 64-byte little-endian header, followed by actual table data.
   No pointers, native-size integers, cached masks, or struct padding on disk. */
#define PIVOT_HEADER 64u
#define PIVOT_VERSION 1u

static uint64_t get_u64(const unsigned char *p) {
    uint64_t value=0;
    unsigned i;
    for (i=0;i<8;++i) value|=(uint64_t)p[i]<<(8*i);
    return value;
}

static void put_u64(unsigned char *p,uint64_t value) {
    unsigned i;
    for (i=0;i<8;++i) p[i]=(unsigned char)(value>>(8*i));
}

static int geometry(size_t n,size_t d,size_t count,size_t *bytes) {
    size_t rows,cells,stride;
    if (n<2 || n>257 || !d || !count) return -1;
    rows=4*(n-1);
    if (d>SIZE_MAX/rows) return -1;
    cells=rows*d;
    /* Two sides: A/B1/B2 and 64-bit C; pivot A/signs/64-bit C. */
    if (cells>(SIZE_MAX-27*rows)/6) return -1;
    stride=6*cells+27*rows;
    if (count>(SIZE_MAX-PIVOT_HEADER)/stride ||
        count>SIZE_MAX/sizeof(SubseqHash3PivotTables)) return -1;
    *bytes=PIVOT_HEADER+count*stride;
    if ((uintmax_t)*bytes>(uintmax_t)INT64_MAX) return -1;
    return 0;
}

static int write_residues(FILE *out,const size_t *values,size_t count) {
    size_t i;
    unsigned char bytes[8];
    for (i=0;i<count;++i) {
        put_u64(bytes,(uint64_t)values[i]);
        if (fwrite(bytes,1,8,out)!=8) return -1;
    }
    return 0;
}

static int read_residues(FILE *in,size_t *values,size_t count,size_t d) {
    size_t i;
    unsigned char bytes[8];
    for (i=0;i<count;++i) {
        uint64_t value;
        if (fread(bytes,1,8,in)!=8) return -1;
        value=get_u64(bytes);
        if (value>=d) return -1;
        values[i]=(size_t)value;
    }
    return 0;
}

static int write_family(FILE *out,const SubseqHash3PivotTables *t) {
    const size_t rows=4*t->k,cells=rows*t->d;
    unsigned side;
    for (side=0;side<2;++side) {
        const SubseqHash3Tables *s=side?&t->right:&t->left;
        if (fwrite(s->A,1,cells,out)!=cells ||
            fwrite(s->B1,1,cells,out)!=cells ||
            fwrite(s->B2,1,cells,out)!=cells ||
            write_residues(out,s->C,rows)) return -1;
    }
    if (fwrite(t->pivot_A,1,rows,out)!=rows ||
        fwrite(t->pivot_left_sign,1,rows,out)!=rows ||
        fwrite(t->pivot_right_sign,1,rows,out)!=rows ||
        write_residues(out,t->pivot_C,rows)) return -1;
    return 0;
}

int subseqhash3_pivot_table_file_create(const char *path,size_t n,size_t d,
    size_t count,uint64_t seed) {
    unsigned char header[PIVOT_HEADER]={0};
    char *temporary=NULL;
    FILE *out=NULL;
    size_t bytes,family,length;
    int fd=-1,status=-1;
    if (!path || geometry(n,d,count,&bytes)) {
        fputs("invalid or oversized pivot table dimensions\n",stderr); return -1;
    }
    length=strlen(path);
    if (length>SIZE_MAX-12) return -1;
    temporary=malloc(length+12);
    if (!temporary) return -1;
    memcpy(temporary,path,length);
    memcpy(temporary+length,".tmp.XXXXXX",12);
    fd=mkstemp(temporary);
    if (fd<0) { perror(path); goto done; }
    out=fdopen(fd,"wb");
    if (!out) { close(fd); goto done; }
    memcpy(header,"SQH3PIV\0",8);
    put_u64(header+8,PIVOT_VERSION);
    put_u64(header+16,n); put_u64(header+24,d); put_u64(header+32,count);
    put_u64(header+40,bytes); put_u64(header+48,seed);
    if (fwrite(header,1,sizeof(header),out)!=sizeof(header)) goto done;
    for (family=0;family<count;++family) {
        SubseqHash3PivotTables t={0};
        int failed;
        if (subseqhash3_pivot_tables_init(&t,n,d,
                seed+family*UINT64_C(0x9e3779b97f4a7c15))) goto done;
        failed=write_family(out,&t);
        subseqhash3_pivot_tables_destroy(&t);
        if (failed) goto done;
    }
    status=-1;
    if (fflush(out) || fsync(fileno(out))) goto done;
    if (fclose(out)) { out=NULL; goto done; }
    out=NULL;
    if (rename(temporary,path)) { perror(path); goto done; }
    status=0;
done:
    if (out) fclose(out);
    if (status) {
        fprintf(stderr,"could not create pivot table file: %s\n",path);
        if (fd>=0) unlink(temporary);
    }
    free(temporary);
    return status;
}

static void *aligned_bytes(size_t bytes) {
    void *p=NULL;
    return posix_memalign(&p,64,bytes)?NULL:p;
}

/* Allocate without random-table generation: loading is linear in file size. */
static int allocate_family(SubseqHash3PivotTables *t,size_t n,size_t d) {
    const size_t k=n-1,rows=4*k,cells=rows*d,blocks=(d+31)/32;
    unsigned side;
    t->k=k; t->d=d;
    for (side=0;side<2;++side) {
        SubseqHash3Tables *s=side?&t->right:&t->left;
        size_t tail,lane;
        s->k=k; s->d=d; s->blocks=blocks;
        s->A=malloc(cells); s->B1=malloc(cells); s->B2=malloc(cells);
        s->signed_A=aligned_bytes(cells);
        s->B1_masks=calloc(rows*blocks,sizeof(uint32_t));
        s->C=malloc(rows*sizeof(size_t));
        s->wrap_indices=malloc(1024*sizeof(uint16_t));
        if (!s->A || !s->B1 || !s->B2 || !s->signed_A ||
            !s->B1_masks || !s->C || !s->wrap_indices) return -1;
        for (tail=0;tail<32;++tail) for (lane=0;lane<32;++lane)
            s->wrap_indices[32*tail+lane]=(uint16_t)(
                lane<tail?lane:32+lane-tail);
    }
    t->pivot_A=malloc(rows); t->pivot_left_sign=malloc(rows);
    t->pivot_right_sign=malloc(rows); t->pivot_C=malloc(rows*sizeof(size_t));
    t->right_reversed_A=aligned_bytes(cells);
    t->right_reversed_B1_masks=calloc(rows*blocks,sizeof(uint32_t));
    return !t->pivot_A || !t->pivot_left_sign || !t->pivot_right_sign ||
        !t->pivot_C || !t->right_reversed_A || !t->right_reversed_B1_masks?-1:0;
}

static int read_family(FILE *in,SubseqHash3PivotTables *t) {
    const size_t rows=4*t->k,d=t->d,cells=rows*d;
    size_t i,row,r;
    unsigned side;
    for (side=0;side<2;++side) {
        SubseqHash3Tables *s=side?&t->right:&t->left;
        if (fread(s->A,1,cells,in)!=cells || fread(s->B1,1,cells,in)!=cells ||
            fread(s->B2,1,cells,in)!=cells || read_residues(in,s->C,rows,d))
            return -1;
        for (i=0;i<cells;++i) {
            if (s->A[i]<1 || (s->B1[i]!=1 && s->B1[i]!=-1) ||
                (s->B2[i]!=1 && s->B2[i]!=-1)) return -1;
            s->signed_A[i]=(int8_t)(s->A[i]*s->B2[i]);
            if (s->B1[i]>0)
                s->B1_masks[(i/d)*s->blocks+(i%d)/32]|=
                    UINT32_C(1)<<((i%d)%32);
        }
    }
    if (fread(t->pivot_A,1,rows,in)!=rows ||
        fread(t->pivot_left_sign,1,rows,in)!=rows ||
        fread(t->pivot_right_sign,1,rows,in)!=rows ||
        read_residues(in,t->pivot_C,rows,d)) return -1;
    for (i=0;i<rows;++i)
        if (!t->pivot_A[i] || t->pivot_A[i]==INT8_MIN ||
            (t->pivot_left_sign[i]!=1 && t->pivot_left_sign[i]!=-1) ||
            (t->pivot_right_sign[i]!=1 && t->pivot_right_sign[i]!=-1)) return -1;
    for (row=0;row<rows;++row) for (r=0;r<d;++r) {
        const size_t src=row*d+(r?d-r:0);
        t->right_reversed_A[row*d+r]=t->right.signed_A[src];
        if (t->right.B1[src]>0)
            t->right_reversed_B1_masks[row*t->right.blocks+r/32]|=
                UINT32_C(1)<<(r%32);
    }
    return 0;
}

void subseqhash3_pivot_table_file_close(SubseqHash3PivotTableFile *file) {
    size_t i;
    if (!file) return;
    if (file->tables) for (i=0;i<file->count;++i)
        subseqhash3_pivot_tables_destroy(file->tables+i);
    free(file->tables);
    memset(file,0,sizeof(*file));
}

int subseqhash3_pivot_table_file_open(SubseqHash3PivotTableFile *file,const char *path) {
    unsigned char header[PIVOT_HEADER];
    struct stat info;
    uint64_t n,d,count;
    size_t bytes,family;
    FILE *in;
    int status=-1;
    if (!file || !path) return -1;
    memset(file,0,sizeof(*file));
    in=fopen(path,"rb");
    if (!in) { perror(path); return -1; }
    if (fstat(fileno(in),&info) || fread(header,1,sizeof(header),in)!=sizeof(header) ||
        memcmp(header,"SQH3PIV\0",8) || get_u64(header+8)!=PIVOT_VERSION ||
        get_u64(header+56)) goto done;
    n=get_u64(header+16); d=get_u64(header+24); count=get_u64(header+32);
    if (n>SIZE_MAX || d>SIZE_MAX || count>SIZE_MAX ||
        geometry((size_t)n,(size_t)d,(size_t)count,&bytes) ||
        info.st_size<0 || (uintmax_t)info.st_size!=(uintmax_t)bytes ||
        get_u64(header+40)!=bytes) goto done;
    file->n=(size_t)n; file->d=(size_t)d; file->count=(size_t)count;
    file->tables=calloc(file->count,sizeof(*file->tables));
    if (!file->tables) goto done;
    for (family=0;family<file->count;++family)
        if (allocate_family(file->tables+family,file->n,file->d) ||
            read_family(in,file->tables+family)) goto done;
    if (fgetc(in)!=EOF || ferror(in)) goto done;
    status=0;
done:
    if (fclose(in)) status=-1;
    if (status) {
        fprintf(stderr,"invalid or incompatible pivot table file: %s "
            "(use 'tables --pivots')\n",path);
        subseqhash3_pivot_table_file_close(file);
    }
    return status;
}
