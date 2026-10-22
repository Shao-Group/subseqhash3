#define _POSIX_C_SOURCE 200809L
#include "pivot_table_file.h"
#include "table_file.h"
#include "ssh3_commands.h"

#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

static void same_side(const SubseqHash3Tables *a,const SubseqHash3Tables *b) {
    const size_t cells=a->k*4*a->d,rows=a->k*4;
    assert(a->k==b->k && a->d==b->d && a->blocks==b->blocks);
    assert(!memcmp(a->A,b->A,cells));
    assert(!memcmp(a->B1,b->B1,cells));
    assert(!memcmp(a->B2,b->B2,cells));
    assert(!memcmp(a->signed_A,b->signed_A,cells));
    assert(!memcmp(a->B1_masks,b->B1_masks,rows*a->blocks*sizeof(uint32_t)));
    assert(!memcmp(a->C,b->C,rows*sizeof(size_t)));
    assert(!memcmp(a->wrap_indices,b->wrap_indices,1024*sizeof(uint16_t)));
}

static void same_family(const SubseqHash3PivotTables *a,const SubseqHash3PivotTables *b) {
    const size_t rows=a->k*4;
    assert(a->k==b->k && a->d==b->d);
    same_side(&a->left,&b->left); same_side(&a->right,&b->right);
    assert(!memcmp(a->pivot_A,b->pivot_A,rows));
    assert(!memcmp(a->pivot_left_sign,b->pivot_left_sign,rows));
    assert(!memcmp(a->pivot_right_sign,b->pivot_right_sign,rows));
    assert(!memcmp(a->pivot_C,b->pivot_C,rows*sizeof(size_t)));
    assert(!memcmp(a->right_reversed_A,b->right_reversed_A,rows*a->d));
    assert(!memcmp(a->right_reversed_B1_masks,b->right_reversed_B1_masks,
        rows*a->right.blocks*sizeof(uint32_t)));
}

static void same_result(const SubseqHash3EditResult *a,const SubseqHash3EditResult *b) {
    assert(a->length==b->length && a->psi==b->psi && a->omega==b->omega);
    assert(a->deletions==b->deletions && a->substitutions==b->substitutions);
    assert(!memcmp(a->seed,b->seed,a->length));
    assert(!memcmp(a->positions,b->positions,a->length*sizeof(size_t)));
    assert(!memcmp(a->substituted,b->substituted,a->length));
}

static void round_trip(const char *path,size_t n,size_t d) {
    SubseqHash3PivotTables expected[3]={{0}};
    SubseqHash3PivotTableFile actual={0};
    SubseqHash3PivotWorkspace *workspace=subseqhash3_pivot_workspace_create();
    assert(workspace);
    const uint64_t seed=12345;
    size_t i,length;
    assert(!subseqhash3_pivot_table_file_create(path,n,d,3,seed));
    assert(!subseqhash3_pivot_table_file_open(&actual,path));
    assert(actual.n==n && actual.d==d && actual.count==3);
    for (i=0;i<3;++i) {
        assert(!subseqhash3_pivot_tables_init(expected+i,n,d,
            seed+i*UINT64_C(0x9e3779b97f4a7c15)));
        same_family(actual.tables+i,expected+i);
    }
    for (length=2;length<=n;++length) {
        uint8_t sequence[8];
        SubseqHash3EditResult a[16]={{0}},b[16]={{0}};
        const size_t repeats=2*(length-1)+1,e=length/2;
        for (i=0;i<length;++i) sequence[i]=(uint8_t)((i*3+length)%4);
        assert(!subseqhash3_compute_pivot_repeats(sequence,length,e,
            actual.tables,3,repeats,a));
        assert(!subseqhash3_compute_pivot_repeats(sequence,length,e,
            expected,3,repeats,b));
        for (i=0;i<repeats;++i) same_result(a+i,b+i);
        subseqhash3_edit_results_destroy(b,repeats);
        SubseqHash3EditResult c[16]={{0}};
        assert(!subseqhash3_compute_outward_pivot_repeats_with_workspace(sequence,length,e,
            actual.tables,3,repeats,b,workspace));
        assert(!subseqhash3_compute_outward_pivot_repeats(sequence,length,e,
            expected,3,repeats,c));
        for (i=0;i<repeats;++i) same_result(b+i,c+i);
        subseqhash3_edit_results_destroy(c,repeats);
        assert(!subseqhash3_compute_pivot_repeats_with_workspace(sequence,length,e,
            actual.tables,3,repeats,c,workspace));
        for (i=0;i<repeats;++i) same_result(a+i,c+i);
        subseqhash3_edit_results_destroy(c,repeats);
        assert(!subseqhash3_compute_outward_pivot_repeats_with_workspace(sequence,length,e,
            actual.tables,3,repeats,c,workspace));
        for (i=0;i<repeats;++i) same_result(b+i,c+i);
        subseqhash3_edit_results_destroy(a,repeats);
        subseqhash3_edit_results_destroy(b,repeats);
        subseqhash3_edit_results_destroy(c,repeats);
    }
    subseqhash3_pivot_workspace_destroy(workspace);
    for (i=0;i<3;++i) subseqhash3_pivot_tables_destroy(expected+i);
    subseqhash3_pivot_table_file_close(&actual);
}

static void bad_files(const char *path) {
    const size_t n=5,d=31,cells=(n-1)*4*d;
    size_t test;
    for (test=0;test<13;++test) {
        SubseqHash3PivotTableFile file={0};
        FILE *out;
        assert(!subseqhash3_pivot_table_file_create(path,n,d,1,7));
        out=fopen(path,"r+b"); assert(out);
        if (test==0) { assert(!fseek(out,0,SEEK_SET)); assert(fputc('X',out)!=EOF); }
        if (test==1) { assert(!fseek(out,8,SEEK_SET)); assert(fputc(2,out)!=EOF); }
        if (test==2) { assert(!fseek(out,16,SEEK_SET)); assert(fputc(1,out)!=EOF); }
        if (test==3) { assert(!fseek(out,24,SEEK_SET)); assert(fputc(0,out)!=EOF); }
        if (test==4) { assert(!fseek(out,32,SEEK_SET)); assert(fputc(0,out)!=EOF); }
        if (test==5) { assert(!fseek(out,40,SEEK_SET)); assert(fputc(0,out)!=EOF); }
        if (test==6) { assert(!fseek(out,56,SEEK_SET)); assert(fputc(1,out)!=EOF); }
        if (test==7) { assert(!fseek(out,64,SEEK_SET)); assert(fputc(0,out)!=EOF); }
        if (test==8) { assert(!fseek(out,(long)(64+cells),SEEK_SET)); assert(fputc(2,out)!=EOF); }
        if (test==9) { assert(!fseek(out,(long)(64+3*cells),SEEK_SET)); assert(fputc(31,out)!=EOF); }
        if (test==10) assert(!ftruncate(fileno(out),63));
        if (test==11) {
            struct stat info;
            assert(!fstat(fileno(out),&info));
            assert(!ftruncate(fileno(out),info.st_size-1));
        }
        if (test==12) { assert(!fseek(out,0,SEEK_END)); assert(fputc(0,out)!=EOF); }
        assert(!fclose(out));
        assert(subseqhash3_pivot_table_file_open(&file,path)==-1);
        assert(!file.tables && !file.n && !file.count);
        subseqhash3_pivot_table_file_close(&file);
    }
}

int main(void) {
    char directory[]="/tmp/ssh3-pivot-files-XXXXXX",tables[256];
    const size_t ds[]={1,31,32,33,1025};
    assert(mkdtemp(directory));
    snprintf(tables,sizeof(tables),"%s/tables.piv",directory);
    for (size_t i=0;i<sizeof(ds)/sizeof(ds[0]);++i) round_trip(tables,7,ds[i]);
    round_trip(tables,2,1); bad_files(tables);
    for (size_t w=1;w<=257;w+=256) {
        SubseqHash3PivotTableFile file={0};
        assert(!subseqhash3_seed_tables_create(tables,w,1,2,42));
        assert(!subseqhash3_pivot_table_file_open(&file,tables));
        assert(file.n==(w==257?257:w+1) && file.count==2);
        subseqhash3_pivot_table_file_close(&file);
    }
    assert(subseqhash3_seed_tables_create(tables,0,1,1,42));
    assert(subseqhash3_seed_tables_create(tables,258,1,1,42));
    assert(!unlink(tables)); assert(!rmdir(directory));
    puts("pivot table round trips, independent sets, malformed files and W-sized generation passed");
    return 0;
}
