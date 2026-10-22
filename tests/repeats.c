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

static void put_text(const char *path,const char *text) {
    FILE *out=fopen(path,"wb");
    assert(out);
    assert(fwrite(text,1,strlen(text),out)==strlen(text));
    assert(!fclose(out));
}

static void same_file(const char *a,const char *b) {
    FILE *x=fopen(a,"rb"),*y=fopen(b,"rb");
    int p,q;
    assert(x && y);
    do { p=fgetc(x); q=fgetc(y); assert(p==q); } while (p!=EOF);
    assert(!ferror(x) && !ferror(y));
    assert(!fclose(x) && !fclose(y));
}

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

static void check_csv(const char *path,const char *tables_path,size_t w,size_t e,size_t r,int outward) {
    const uint8_t inputs[2][6]={{0,1,2,3,0,1},{3,3,3,3,0,1}};
    SubseqHash3PivotTableFile tables={0};
    SubseqHash3EditResult results[10]={{0}};
    FILE *in=fopen(path,"r");
    char actual[512],expected[512],seed[8];
    size_t sequence,start,repeat,j;
    assert(r<=10 && in);
    assert(!subseqhash3_pivot_table_file_open(&tables,tables_path));
    assert(fgets(actual,sizeof(actual),in));
    assert(!strcmp(actual,"repeat,family,pivot,sequence,window_start,substitutions,"
        "deletions,seed,psi,omega\n"));
    for (sequence=0;sequence<2;++sequence) for (start=0;start<=6-w;++start) {
        assert(!(outward?subseqhash3_compute_outward_pivot_repeats:
            subseqhash3_compute_pivot_repeats)(inputs[sequence]+start,w,e,
            tables.tables,tables.count,r,results));
        for (repeat=0;repeat<r;++repeat) {
            for (j=0;j<w-1;++j) seed[j]="ACGT"[results[repeat].seed[j]];
            seed[w-1]=0;
            snprintf(expected,sizeof(expected),"%zu,%zu,%zu,%zu,%zu,%zu,1,%s,%zu,%.0f\n",
                repeat,repeat/(w-1),repeat%(w-1),sequence,start,e,seed,
                results[repeat].psi,results[repeat].omega);
            assert(fgets(actual,sizeof(actual),in));
            assert(!strcmp(actual,expected));
        }
        subseqhash3_edit_results_destroy(results,r);
    }
    assert(fgetc(in)==EOF && !ferror(in));
    assert(!fclose(in));
    subseqhash3_pivot_table_file_close(&tables);
}

int main(void) {
    char directory[]="/tmp/ssh3-repeats-tests-XXXXXX";
    char input[256],tables[256],legacy[256],csv[256],other[256],generated[260];
    const size_t ds[]={1,31,32,33,1025};
    size_t i;
    assert(mkdtemp(directory));
#define PATH(NAME,SUFFIX) assert(snprintf(NAME,sizeof(NAME),"%s/%s",directory,SUFFIX)>0)
    PATH(input,"input.txt"); PATH(tables,"tables.piv"); PATH(legacy,"legacy.abc");
    PATH(csv,"output.csv"); PATH(other,"other.csv");
#undef PATH
    for (i=0;i<sizeof(ds)/sizeof(ds[0]);++i) round_trip(tables,7,ds[i]);
    round_trip(tables,2,1);
    bad_files(tables);
    {
        char *make[]={"tables","--pivots","4","31","10",tables,"--seed","42"};
        char *make_again[]={"tables","4","31","10",other,"--seed","42","--pivots"};
        char *run[]={"repeats",input,"4","2","31","10",csv,tables};
        SubseqHash3PivotTableFile file={0};
        put_text(input,"acgTaC\r\nTTTTAC\n");
        assert(!ssh3_tables_main(8,make));
        assert(!ssh3_tables_main(8,make_again));
        same_file(tables,other);
        assert(!subseqhash3_pivot_table_file_open(&file,tables));
        assert(file.n==4 && file.d==31 && file.count==4);
        subseqhash3_pivot_table_file_close(&file);
        assert(!ssh3_repeats_main(8,run));
        check_csv(csv,tables,4,2,10,1);
        {
            char *outward[]={"repeats","--outward",input,"4","2","31","10",csv,tables};
            assert(!ssh3_repeats_main(9,outward));
            check_csv(csv,tables,4,2,10,1);
            char *inward[]={"repeats","--inward",input,"4","2","31","10",csv,tables};
            assert(!ssh3_repeats_main(9,inward));
            check_csv(csv,tables,4,2,10,0);
            inward[1]=input; inward[2]="--inward";
            assert(!ssh3_repeats_main(9,inward));
            check_csv(csv,tables,4,2,10,0);
            char *duplicate[]={"repeats",input,"4","2","31","10",csv,tables,"--outward","--outward"};
            assert(ssh3_repeats_main(10,duplicate)!=0);
            duplicate[8]="--inward"; duplicate[9]="--inward";
            assert(ssh3_repeats_main(10,duplicate)!=0);
            duplicate[9]="--outward";
            assert(ssh3_repeats_main(10,duplicate)!=0);
            duplicate[8]="--outward"; duplicate[9]="--inward";
            assert(ssh3_repeats_main(10,duplicate)!=0);
        }
        for (i=0;i<2;++i) {
            run[3]=i?"3":"0";
            assert(!ssh3_repeats_main(8,run));
            check_csv(csv,tables,4,i?3:0,10,1);
        }
        run[2]="3"; run[3]="0"; run[5]="6";
        assert(!ssh3_repeats_main(8,run));
        check_csv(csv,tables,3,0,6,1);
        run[5]="10";
        assert(ssh3_repeats_main(8,run)!=0); /* Shorter windows need more families. */
        run[2]="4"; run[3]="2";
        assert(!ssh3_repeats_main(7,run)); /* Auto-created tables. */
        snprintf(generated,sizeof(generated),"%s.piv",csv);
        run[6]=other; run[7]=generated;
        assert(!ssh3_repeats_main(8,run));
        same_file(csv,other);
        check_csv(other,generated,4,2,10,1);
        assert(!unlink(generated));
        run[6]=csv; run[7]=tables;
        {
            const char *bad_w[]={"0","1","258","-1","18446744073709551616"};
            for (i=0;i<5;++i) {
                run[2]=(char *)bad_w[i]; assert(ssh3_repeats_main(8,run)!=0);
            }
        }
        run[2]="4"; run[3]="4"; assert(ssh3_repeats_main(8,run)!=0);
        run[3]="-1"; assert(ssh3_repeats_main(8,run)!=0);
        run[3]="2"; run[4]="0"; assert(ssh3_repeats_main(8,run)!=0);
        run[4]="31"; run[5]="0"; assert(ssh3_repeats_main(8,run)!=0);
        run[5]="10"; run[6]=input; assert(ssh3_repeats_main(8,run)!=0);
        run[6]=csv; run[7]=csv; assert(ssh3_repeats_main(8,run)!=0);
        run[7]=tables;
        put_text(input,"ACNTAC\n"); assert(ssh3_repeats_main(8,run)!=0);
        put_text(input,"\n"); assert(ssh3_repeats_main(8,run)!=0);
        put_text(input,"AC\n"); assert(ssh3_repeats_main(8,run)!=0);
        put_text(input,"ACGTAC\n");
        assert(!subseqhash3_table_file_create(legacy,4,31,2,42));
        run[7]=legacy;
        put_text(csv,"KEEP\n"); put_text(other,"KEEP\n");
        assert(ssh3_repeats_main(8,run)!=0);
        same_file(csv,other); /* Bad tables do not truncate an existing output. */
        {
            SubseqHash3TableFile old={0};
            assert(subseqhash3_table_file_open(&old,tables)!=0);
            subseqhash3_table_file_close(&old);
        }
    }
    {
        char *make[]={"tables","4","31","2",legacy,"--seed","42"};
        char *run[]={"seeds",input,"4","2","31","2",other,legacy};
        char *residues[]={"32","1024","1025"};
        assert(!ssh3_tables_main(7,make));
        assert(!ssh3_seeds_main(8,run));
        for (i=0;i<3;++i) {
            make[2]=residues[i]; run[4]=residues[i];
            assert(!ssh3_tables_main(7,make));
            assert(!ssh3_seeds_main(8,run));
        }
        make[1]="-1"; assert(ssh3_tables_main(7,make)!=0);
        make[1]="4"; make[6]="-1"; assert(ssh3_tables_main(7,make)!=0);
    }
    {
        char *make[]={"tables","--pivots","2","1","5",tables,"--seed","77"};
        char *run[]={"repeats",input,"2","1","1","5",csv,tables};
        SubseqHash3PivotTableFile loaded={0};
        put_text(input,"AC\n");
        assert(!ssh3_tables_main(8,make));
        assert(!subseqhash3_pivot_table_file_open(&loaded,tables));
        assert(loaded.n==2 && loaded.count==5);
        subseqhash3_pivot_table_file_close(&loaded);
        assert(!ssh3_repeats_main(8,run));
        make[2]="257"; make[4]="1";
        run[2]="257"; run[3]="256"; run[5]="1";
        {
            char sequence[259],row[1024];
            const char prefix[]="0,0,0,0,0,256,1,";
            FILE *out;
            memset(sequence,'A',257); sequence[257]='\n'; sequence[258]=0;
            put_text(input,sequence);
            assert(!ssh3_tables_main(8,make));
            assert(!ssh3_repeats_main(8,run));
            char *inward_args[9];
            for (size_t arg=0;arg<8;++arg) inward_args[arg]=run[arg];
            inward_args[8]="--inward";
            assert(!ssh3_repeats_main(9,inward_args));
            out=fopen(csv,"r"); assert(out);
            assert(fgets(row,sizeof(row),out));
            assert(fgets(row,sizeof(row),out));
            assert(!strncmp(row,prefix,strlen(prefix)));
            assert(strspn(row+strlen(prefix),"CGT")==256);
            assert(fgetc(out)==EOF && !ferror(out));
            assert(!fclose(out));
        }
        make[2]="258"; assert(ssh3_tables_main(8,make)!=0);
        make[2]="1"; assert(ssh3_tables_main(8,make)!=0);
        /* Same-inode aliases are rejected without truncating the input. */
        assert(!unlink(other)); assert(!link(input,other));
        run[6]=other; assert(ssh3_repeats_main(8,run)!=0);
        same_file(input,other);
        assert(!unlink(other)); put_text(other,"KEEP\n");
        /* Automatic table generation must not overwrite its input. */
        snprintf(generated,sizeof(generated),"%s.piv",csv);
        put_text(generated,"ACGT\n");
        run[1]=generated; run[2]="4"; run[3]="0"; run[5]="1"; run[6]=csv;
        assert(ssh3_repeats_main(7,run)!=0);
        assert(!unlink(generated));
        /* Invalid table dimensions do not replace an existing destination. */
        put_text(csv,"KEEP\n");
        assert(subseqhash3_pivot_table_file_create(csv,1,31,1,42)!=0);
        same_file(csv,other);
    }
    assert(!unlink(input)); assert(!unlink(tables)); assert(!unlink(legacy));
    assert(!unlink(csv)); assert(!unlink(other)); assert(!rmdir(directory));
    puts("pivot table round trips, malformed files, CSV/API equivalence, auto/reused "
         "tables, edit limits, and legacy CLI compatibility passed");
    return 0;
}
