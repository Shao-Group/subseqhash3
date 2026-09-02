#define _POSIX_C_SOURCE 200809L
#include "table_file.h"

#include <errno.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#define TABLE_VERSION UINT32_C(2)
#define TABLE_ENDIAN UINT32_C(0x01020304)

typedef struct {
    char magic[8];
    uint32_t version,endian,word_size,reserved;
    uint64_t n,d,count,file_size,set_base,set_stride,wrap_offset;
    uint64_t a_offset,b1_offset,b2_offset,signed_a_offset,masks_offset,c_offset;
    uint64_t f_offset,f_stride;
    uint64_t generation_seed;
    uint8_t padding[128];
} TableHeader;

static size_t align64(size_t value) {
    return (value+63u)&~(size_t)63u;
}

static int layout(size_t n,size_t d,size_t count,TableHeader *header) {
    size_t cells,blocks,mask_cells,f_cells,offset,stride,total;
    if (!n || !d || !count || n>SIZE_MAX/d || n*d>SIZE_MAX/4) return -1;
    cells=n*d*4; blocks=(d+31)/32;
    if (n>SIZE_MAX/128 || 128*n==SIZE_MAX || d>SIZE_MAX/(128*n+1)) return -1;
    f_cells=d*(128*n+1);
    if (n>SIZE_MAX/4 || n*4>SIZE_MAX/blocks) return -1;
    mask_cells=n*4*blocks;
    memset(header,0,sizeof(*header));
    memcpy(header->magic,"SQH3TAB\0",8);
    header->version=TABLE_VERSION; header->endian=TABLE_ENDIAN;
    header->word_size=(uint32_t)sizeof(size_t);
    header->n=n; header->d=d; header->count=count;
    header->wrap_offset=align64(sizeof(*header));
    header->set_base=align64((size_t)header->wrap_offset+1024*sizeof(uint16_t));
    offset=0;
    header->a_offset=offset; offset=align64(offset+cells);
    header->b1_offset=offset; offset=align64(offset+cells);
    header->b2_offset=offset; offset=align64(offset+cells);
    header->signed_a_offset=offset; offset=align64(offset+cells);
    header->masks_offset=offset;
    if (mask_cells>SIZE_MAX/sizeof(uint32_t)) return -1;
    offset=align64(offset+mask_cells*sizeof(uint32_t));
    header->c_offset=offset;
    if (n*4>SIZE_MAX/sizeof(size_t)) return -1;
    offset=align64(offset+n*4*sizeof(size_t));
    header->f_offset=offset; header->f_stride=128*n+1;
    if (f_cells>(SIZE_MAX-offset)/sizeof(float)) return -1;
    stride=align64(offset+f_cells*sizeof(float));
    if (count>(SIZE_MAX-(size_t)header->set_base)/stride) return -1;
    total=(size_t)header->set_base+count*stride;
    header->set_stride=stride; header->file_size=total;
    return 0;
}

static void copy_set(uint8_t *base,const TableHeader *h,
                     const SubseqHash3Tables *table) {
    size_t cells=(size_t)h->n*(size_t)h->d*4;
    size_t masks=(size_t)h->n*4*(((size_t)h->d+31)/32);
    memcpy(base+h->a_offset,table->A,cells);
    memcpy(base+h->b1_offset,table->B1,cells);
    memcpy(base+h->b2_offset,table->B2,cells);
    memcpy(base+h->signed_a_offset,table->signed_A,cells);
    memcpy(base+h->masks_offset,table->B1_masks,masks*sizeof(uint32_t));
    memcpy(base+h->c_offset,table->C,(size_t)h->n*4*sizeof(size_t));
    memcpy(base+h->f_offset,table->F,
           (size_t)h->d*(size_t)h->f_stride*sizeof(float));
}

int subseqhash3_table_file_create(const char *path,size_t n,size_t d,
                                  size_t count,uint64_t seed) {
    TableHeader header;
    uint8_t *mapping;
    int fd=-1,status=-1;
    size_t set;
    if (layout(n,d,count,&header)) {
        fputs("invalid or oversized table dimensions\n",stderr); return -1;
    }
    header.generation_seed=seed;
    fd=open(path,O_RDWR|O_CREAT|O_TRUNC,0666);
    if (fd<0) { perror(path); return -1; }
    if (ftruncate(fd,(off_t)header.file_size)) { perror(path); goto done; }
    mapping=mmap(NULL,(size_t)header.file_size,PROT_READ|PROT_WRITE,
                 MAP_SHARED,fd,0);
    if (mapping==MAP_FAILED) { perror("mmap"); goto done; }
    memcpy(mapping,&header,sizeof(header));
    for (set=0;set<count;++set) {
        SubseqHash3Tables table={0};
        if (subseqhash3_tables_init(&table,n,d,
                seed+set*UINT64_C(0x9e3779b97f4a7c15))) {
            subseqhash3_tables_destroy(&table); goto unmap;
        }
        if (!set)
            memcpy(mapping+header.wrap_offset,table.wrap_indices,
                   1024*sizeof(uint16_t));
        copy_set(mapping+header.set_base+set*header.set_stride,&header,&table);
        subseqhash3_tables_destroy(&table);
    }
    if (msync(mapping,(size_t)header.file_size,MS_SYNC)) { perror("msync"); goto unmap; }
    status=0;
unmap:
    if (munmap(mapping,(size_t)header.file_size) && !status) status=-1;
done:
    if (close(fd) && !status) status=-1;
    return status;
}

int subseqhash3_table_file_open(SubseqHash3TableFile *file,const char *path) {
    struct stat info;
    TableHeader *h,expected;
    uint8_t *mapping;
    size_t set;
    int fd;
    memset(file,0,sizeof(*file));
    fd=open(path,O_RDONLY);
    if (fd<0) { perror(path); return -1; }
    if (fstat(fd,&info) || info.st_size<(off_t)sizeof(*h)) {
        fprintf(stderr,"invalid table file: %s\n",path); close(fd); return -1;
    }
    mapping=mmap(NULL,(size_t)info.st_size,PROT_READ,MAP_PRIVATE,fd,0);
    close(fd);
    if (mapping==MAP_FAILED) { perror("mmap"); return -1; }
    h=(TableHeader *)mapping;
    if (memcmp(h->magic,"SQH3TAB\0",8) || h->version!=TABLE_VERSION ||
        h->endian!=TABLE_ENDIAN || h->word_size!=sizeof(size_t) ||
        h->file_size!=(uint64_t)info.st_size || !h->n || !h->d || !h->count ||
        h->n>SIZE_MAX || h->d>SIZE_MAX || h->count>SIZE_MAX ||
        layout((size_t)h->n,(size_t)h->d,(size_t)h->count,&expected) ||
        h->file_size!=expected.file_size || h->set_base!=expected.set_base ||
        h->set_stride!=expected.set_stride ||
        h->wrap_offset!=expected.wrap_offset || h->a_offset!=expected.a_offset ||
        h->b1_offset!=expected.b1_offset || h->b2_offset!=expected.b2_offset ||
        h->signed_a_offset!=expected.signed_a_offset ||
        h->masks_offset!=expected.masks_offset || h->c_offset!=expected.c_offset ||
        h->f_offset!=expected.f_offset || h->f_stride!=expected.f_stride) {
        fprintf(stderr,"invalid or incompatible table file: %s\n",path);
        munmap(mapping,(size_t)info.st_size); return -1;
    }
    file->tables=calloc((size_t)h->count,sizeof(*file->tables));
    if (!file->tables) { munmap(mapping,(size_t)info.st_size); return -1; }
    file->mapping=mapping; file->mapping_size=(size_t)info.st_size;
    file->n=(size_t)h->n; file->d=(size_t)h->d; file->count=(size_t)h->count;
    for (set=0;set<file->count;++set) {
        uint8_t *base=mapping+h->set_base+set*h->set_stride;
        SubseqHash3Tables *t=file->tables+set;
        t->k=file->n; t->d=file->d; t->blocks=(file->d+31)/32;
        t->A=(int8_t *)(base+h->a_offset);
        t->B1=(int8_t *)(base+h->b1_offset);
        t->B2=(int8_t *)(base+h->b2_offset);
        t->signed_A=(int8_t *)(base+h->signed_a_offset);
        t->B1_masks=(uint32_t *)(base+h->masks_offset);
        t->C=(size_t *)(base+h->c_offset);
        t->F=(float *)(base+h->f_offset);
        t->f_stride=(size_t)h->f_stride;
        t->wrap_indices=(uint16_t *)(mapping+h->wrap_offset);
    }
    return 0;
}

void subseqhash3_table_file_close(SubseqHash3TableFile *file) {
    if (!file) return;
    free(file->tables);
    if (file->mapping) munmap(file->mapping,file->mapping_size);
    memset(file,0,sizeof(*file));
}
