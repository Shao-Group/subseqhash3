#define _POSIX_C_SOURCE 200809L
#include "subseqhash_edits_i8.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static double seconds(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC,&t);
    return (double)t.tv_sec+(double)t.tv_nsec*1e-9;
}
static int parse_size(const char *text,size_t *value) {
    char *end; unsigned long long parsed;
    errno=0; parsed=strtoull(text,&end,10);
    if (errno || !*text || *end || parsed>(unsigned long long)SIZE_MAX)
        return -1;
    *value=(size_t)parsed; return 0;
}
static void trim(char *line,size_t *length) {
    while (*length && (line[*length-1]=='\n' || line[*length-1]=='\r'))
        line[--*length]='\0';
}
static int encode(const char *text,size_t n,uint8_t *s) {
    size_t i;
    for (i=0;i<n;++i) {
        switch (text[i]) {
            case 'A': s[i]=0; break;
            case 'C': s[i]=1; break;
            case 'G': s[i]=2; break;
            case 'T': s[i]=3; break;
            default: return -1;
        }
    }
    return 0;
}

static char *make_seed_filename(const char *input_filename) {
    static const char suffix[]=".i8_twopass_seeds.csv";
    size_t length=strlen(input_filename);
    char *name;
    if (length>SIZE_MAX-sizeof(suffix)) return NULL;
    name=malloc(length+sizeof(suffix));
    if (!name) return NULL;
    memcpy(name,input_filename,length);
    memcpy(name+length,suffix,sizeof(suffix));
    return name;
}

static int write_seed_records(FILE *output,size_t repeat,size_t pair,
                              unsigned side,const char *input,
                              const SubseqHashEditResult *results,
                              size_t edits) {
    static const char alphabet[4]={'A','C','G','T'};
    size_t split,i;
    for (split=0;split<=edits;++split) {
        const SubseqHashEditResult *result=results+split;
        if (fprintf(output,"%zu,%zu,%u,%s,%zu,%zu,%zu,%.0f,",
                    repeat,pair,side,input,split,edits-split,
                    result->psi,result->omega)<0) return -1;
        for (i=0;i<result->length;++i) {
            if (result->seed[i]>3 ||
                fputc(alphabet[result->seed[i]],output)==EOF) return -1;
        }
        if (fputc('\n',output)==EOF) return -1;
    }
    return 0;
}

static int benchmark_file(const char *filename,size_t n,size_t k,size_t d,
                          size_t edits,size_t repeats,uint64_t seed_base) {
    const size_t baseline_edits=n-k;
    FILE *input=fopen(filename,"r");
    FILE *seed_output=NULL;
    char *seed_filename=NULL;
    char *line1=NULL,*line2=NULL;
    size_t cap1=0,cap2=0,repeat,split;
    uint8_t *s1=NULL,*s2=NULL;
    double pairs_sum=0,baseline_time_sum=0,baseline_rate_sum=0;
    double baseline_collisions_sum=0,edit_time_sum=0,edit_rate_sum=0;
    double any_sum=0,all_sum=0,split_sum[32]={0};
    int status=-1;
    if (!input) { perror(filename); return -1; }
    seed_filename=make_seed_filename(filename);
    if (!seed_filename || !(seed_output=fopen(seed_filename,"w"))) {
        if (seed_filename) perror(seed_filename);
        else fputs("seed filename allocation failed\n",stderr);
        goto done;
    }
    if (fputs("repeat,pair,side,input,Es,Ed,psi,omega,seed\n",
              seed_output)==EOF) goto done;
    s1=malloc(n?n:1); s2=malloc(n?n:1);
    if (!s1 || !s2) goto done;

    for (repeat=0;repeat<repeats;++repeat) {
        SubseqHashI8Tables tables={0};
        size_t pairs=0,baseline_collisions=0,any_collisions=0;
        size_t all_collisions=0,split_collisions[32]={0};
        double baseline_time=0,edit_time=0;
        int run_status=-1;
        if (subseqhash_i8_tables_init(
                &tables,n,d,
                seed_base+repeat*UINT64_C(0x9e3779b97f4a7c15)))
            goto done;
        rewind(input); clearerr(input);
        for (;;) {
            ssize_t read1=getline(&line1,&cap1,input),read2;
            size_t length1,length2;
            SubseqHashEditResult baseline1[32],baseline2[32];
            SubseqHashEditResult left[32],right[32];
            double t0; int any=0,all=1;
            if (read1<0) { run_status=ferror(input)?-1:0; break; }
            read2=getline(&line2,&cap2,input);
            if (read2<0) {
                fprintf(stderr,"%s has an odd line count\n",filename); break;
            }
            length1=(size_t)read1; length2=(size_t)read2;
            trim(line1,&length1); trim(line2,&length2);
            if (length1!=n || length2!=n ||
                encode(line1,n,s1) || encode(line2,n,s2)) {
                fprintf(stderr,"invalid pair %zu in %s\n",pairs+1,filename);
                break;
            }

            t0=seconds();
            if (subseqhash_compute_all_edit_splits_i8_twopass(
                    s1,n,baseline_edits,&tables,baseline1) ||
                subseqhash_compute_all_edit_splits_i8_twopass(
                    s2,n,baseline_edits,&tables,baseline2)) {
                subseqhash_edit_results_destroy(baseline1,baseline_edits+1);
                break;
            }
            baseline_time+=seconds()-t0;
            baseline_collisions+=(size_t)!memcmp(
                baseline1[0].seed,baseline2[0].seed,k);
            subseqhash_edit_results_destroy(baseline2,baseline_edits+1);
            subseqhash_edit_results_destroy(baseline1,baseline_edits+1);

            t0=seconds();
            if (subseqhash_compute_all_edit_splits_i8_twopass(
                    s1,n,edits,&tables,left) ||
                subseqhash_compute_all_edit_splits_i8_twopass(
                    s2,n,edits,&tables,right)) {
                subseqhash_edit_results_destroy(left,edits+1);
                break;
            }
            edit_time+=seconds()-t0;
            if (write_seed_records(seed_output,repeat,pairs+1,1,line1,
                                   left,edits) ||
                write_seed_records(seed_output,repeat,pairs+1,2,line2,
                                   right,edits)) {
                fputs("failed writing seed file\n",stderr);
                subseqhash_edit_results_destroy(right,edits+1);
                subseqhash_edit_results_destroy(left,edits+1);
                run_status=-1; break;
            }
            for (split=0;split<=edits;++split) {
                const int equal=left[split].length==right[split].length &&
                    !memcmp(left[split].seed,right[split].seed,
                            left[split].length);
                any|=equal; all&=equal;
                split_collisions[split]+=(size_t)equal;
            }
            subseqhash_edit_results_destroy(right,edits+1);
            subseqhash_edit_results_destroy(left,edits+1);
            any_collisions+=(size_t)any;
            all_collisions+=(size_t)all;
            ++pairs;
        }
        subseqhash_i8_tables_destroy(&tables);
        if (run_status) goto done;
        pairs_sum+=pairs; baseline_time_sum+=baseline_time;
        baseline_rate_sum+=baseline_time?(double)pairs/baseline_time:0;
        baseline_collisions_sum+=baseline_collisions;
        edit_time_sum+=edit_time;
        edit_rate_sum+=edit_time?(double)pairs/edit_time:0;
        any_sum+=any_collisions; all_sum+=all_collisions;
        for (split=0;split<=edits;++split)
            split_sum[split]+=split_collisions[split];
    }
    if (fflush(seed_output) || ferror(seed_output)) goto done;
    fprintf(stderr,"seed file: %s\n",seed_filename);
    status=0;
    printf("%s,%zu,%.3f,%.6f,%.3f,%.3f,%.6f,%.3f,%.3f,%.3f",
           filename,repeats,pairs_sum/repeats,baseline_time_sum/repeats,
           baseline_rate_sum/repeats,baseline_collisions_sum/repeats,
           edit_time_sum/repeats,edit_rate_sum/repeats,
           any_sum/repeats,all_sum/repeats);
    for (split=0;split<=edits;++split)
        printf(",%.3f",split_sum[split]/repeats);
    putchar('\n');
done:
    if (seed_output && fclose(seed_output) && !status) status=-1;
    free(seed_filename);
    free(s2); free(s1); free(line2); free(line1); fclose(input);
    return status;
}

int main(int argc,char **argv) {
    size_t n,k,d,edits,repeats,i;
    uint64_t seed_base=(uint64_t)time(NULL)^((uint64_t)clock()<<32);
    int status=0;
    if (argc<7 || parse_size(argv[1],&n) || parse_size(argv[2],&k) ||
        parse_size(argv[3],&d) || parse_size(argv[4],&edits) ||
        parse_size(argv[5],&repeats) || !d || (d&15u) || !repeats ||
        k>n || edits>n || edits>=32 || n-k>=32) {
        fprintf(stderr,
                "usage: %s n k D edits repeats FILE...  (D %% 16 == 0)\n",
                argv[0]);
        return EXIT_FAILURE;
    }
    fprintf(stderr,"backend: %s\n",subseqhash_edits_i8_backend());
    printf("file,repeats,average_pairs,average_i8_twopass_k%zu_seconds,"
           "average_i8_twopass_k%zu_pairs_per_second,"
           "average_i8_twopass_k%zu_collisions,average_i8_twopass_edits_seconds,"
           "average_i8_twopass_edits_pairs_per_second,"
           "average_edit_any_seed_collisions,"
           "average_edit_all_seeds_collisions",k,k,k);
    for (i=0;i<=edits;++i)
        printf(",average_edit_Es%zu_Ed%zu_collisions",i,edits-i);
    putchar('\n');
    for (i=6;i<(size_t)argc;++i)
        if (benchmark_file(argv[i],n,k,d,edits,repeats,seed_base)) {
            status=1; break;
        }
    return status?EXIT_FAILURE:EXIT_SUCCESS;
}
