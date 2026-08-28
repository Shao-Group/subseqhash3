#define _POSIX_C_SOURCE 200809L
#include "table_file.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static double seconds(clockid_t clock_id) {
    struct timespec value;
    if (clock_gettime(clock_id,&value)) return 0.0;
    return (double)value.tv_sec+(double)value.tv_nsec*1e-9;
}

static int parse_size(const char *text,size_t *value) {
    char *end; unsigned long long parsed;
    errno=0; parsed=strtoull(text,&end,10);
    if (errno || !*text || *end || parsed>(unsigned long long)SIZE_MAX)
        return -1;
    *value=(size_t)parsed; return 0;
}

static void help(FILE *stream) {
    fputs("Estimate SubseqHash3 seed-collision probabilities from DNA pairs.\n\n"
          "Usage:\n  subseqhash3 collisions PAIRS_FILE E D R [TABLE_FILE]\n\n"
          "Arguments:\n"
          "  PAIRS_FILE  CSV input with one STRING1,STRING2 pair per line.\n"
          "               Strings in a pair may have different lengths.\n"
          "  E           Edit budget applied independently to each string.\n"
          "  D           Number of hash residues.\n"
          "  R           Number of table sets and repeated trials.\n"
          "  TABLE_FILE  Optional reusable table file matching D and R.\n\n"
          "Only left/right edit splits producing equal-length seeds are compared.\n"
          "If TABLE_FILE is omitted, new tables are written to PAIRS_FILE.ABC.\n",
          stream);
}

static void trim(char *line,size_t *length) {
    while (*length && (line[*length-1]=='\n' || line[*length-1]=='\r'))
        line[--*length]='\0';
}

static int valid_dna(const char *text,size_t length) {
    size_t i;
    for (i=0;i<length;++i)
        if (text[i]!='A' && text[i]!='a' && text[i]!='C' && text[i]!='c' &&
            text[i]!='G' && text[i]!='g' && text[i]!='T' && text[i]!='t')
            return 0;
    return 1;
}

static int parse_pair(char *line,size_t length,char **left,size_t *left_length,
                      char **right,size_t *right_length) {
    char *comma;
    trim(line,&length);
    comma=strchr(line,',');
    if (!comma || strchr(comma+1,',')) return -1;
    *left=line; *left_length=(size_t)(comma-line);
    *right=comma+1; *right_length=length-*left_length-1;
    *comma='\0';
    if (!*left_length || !*right_length ||
        !valid_dna(*left,*left_length) || !valid_dna(*right,*right_length))
        return -1;
    return 0;
}

static int scan_input(const char *filename,size_t edits,size_t *maximum_length,
                      size_t *pair_count) {
    FILE *input=fopen(filename,"r");
    char *line=NULL; size_t capacity=0,line_number=0;
    ssize_t read; int status=-1;
    if (!input) { perror(filename); return -1; }
    *maximum_length=0; *pair_count=0;
    while ((read=getline(&line,&capacity,input))>=0) {
        char *left,*right; size_t n1,n2,difference;
        ++line_number;
        if (parse_pair(line,(size_t)read,&left,&n1,&right,&n2)) {
            fprintf(stderr,"invalid pair on line %zu of %s\n",line_number,filename);
            goto done;
        }
        difference=n1>n2?n1-n2:n2-n1;
        if (n1<edits || n2<edits) {
            fprintf(stderr,"line %zu of %s contains a string shorter than E=%zu\n",
                    line_number,filename,edits);
            goto done;
        }
        if (difference>edits) {
            fprintf(stderr,"line %zu of %s has length difference %zu, greater than E=%zu\n",
                    line_number,filename,difference,edits);
            goto done;
        }
        if (n1>*maximum_length) *maximum_length=n1;
        if (n2>*maximum_length) *maximum_length=n2;
        ++*pair_count;
    }
    if (ferror(input)) { perror(filename); goto done; }
    if (!*pair_count) { fprintf(stderr,"%s contains no pairs\n",filename); goto done; }
    status=0;
done:
    free(line); fclose(input); return status;
}

static int encode(const char *text,size_t length,uint8_t *sequence) {
    size_t i;
    for (i=0;i<length;++i) {
        switch (text[i]) {
        case 'A': case 'a': sequence[i]=0; break;
        case 'C': case 'c': sequence[i]=1; break;
        case 'G': case 'g': sequence[i]=2; break;
        default: sequence[i]=3; break;
        }
    }
    return 0;
}

static int estimate(const char *filename,size_t edits,size_t d,size_t repeats,
                    size_t maximum_length,size_t expected_pairs,
                    const SubseqHash3TableFile *table_file) {
    const size_t per=edits+1;
    FILE *input=fopen(filename,"r");
    char *line=NULL; size_t capacity=0,repeat,left_split,right_split;
    double pairs_sum=0,any_sum=0,all_sum=0,wall_sum=0,cpu_sum=0;
    double *collisions=NULL,*denominators=NULL;
    int status=-1;
    if (!input) { perror(filename); return -1; }
    if (per>SIZE_MAX/per || per*per>SIZE_MAX/sizeof(*collisions)) goto done;
    collisions=calloc(per*per,sizeof(*collisions));
    denominators=calloc(per*per,sizeof(*denominators));
    if (!collisions || !denominators) goto done;

    for (repeat=0;repeat<repeats;++repeat) {
        const SubseqHash3Tables *tables=table_file->tables+repeat;
        size_t pairs=0,any_collisions=0,all_collisions=0;
        double wall_start=seconds(CLOCK_MONOTONIC);
        double cpu_start=seconds(CLOCK_PROCESS_CPUTIME_ID);
        ssize_t read;
        rewind(input); clearerr(input);
        while ((read=getline(&line,&capacity,input))>=0) {
            char *left_text,*right_text; size_t n1,n2,compatible=0;
            uint8_t *s1=NULL,*s2=NULL;
            SubseqHash3EditResult *left=NULL,*right=NULL;
            int any=0,all=1;
            if (parse_pair(line,(size_t)read,&left_text,&n1,&right_text,&n2))
                goto repeat_failed;
            s1=malloc(n1); s2=malloc(n2);
            left=calloc(per,sizeof(*left)); right=calloc(per,sizeof(*right));
            if (!s1 || !s2 || !left || !right) goto pair_failed;
            encode(left_text,n1,s1); encode(right_text,n2,s2);
            if (subseqhash3_compute_all_edit_splits(s1,n1,edits,tables,left) ||
                subseqhash3_compute_all_edit_splits(s2,n2,edits,tables,right))
                goto pair_failed;
            for (left_split=0;left_split<per;++left_split)
                for (right_split=0;right_split<per;++right_split) {
                    size_t index=left_split*per+right_split;
                    int equal;
                    if (left[left_split].length!=right[right_split].length)
                        continue;
                    ++compatible; denominators[index]+=1.0;
                    equal=!memcmp(left[left_split].seed,right[right_split].seed,
                                  left[left_split].length);
                    collisions[index]+=(double)equal;
                    any|=equal; all&=equal;
                }
            if (!compatible) goto pair_failed;
            subseqhash3_edit_results_destroy(left,per);
            subseqhash3_edit_results_destroy(right,per);
            free(left); free(right); free(s1); free(s2);
            any_collisions+=(size_t)any; all_collisions+=(size_t)all; ++pairs;
            continue;
pair_failed:
            subseqhash3_edit_results_destroy(left,per);
            subseqhash3_edit_results_destroy(right,per);
            free(left); free(right); free(s1); free(s2);
            goto repeat_failed;
        }
        if (ferror(input) || pairs!=expected_pairs) goto repeat_failed;
        wall_sum+=seconds(CLOCK_MONOTONIC)-wall_start;
        cpu_sum+=seconds(CLOCK_PROCESS_CPUTIME_ID)-cpu_start;
        pairs_sum+=pairs; any_sum+=any_collisions; all_sum+=all_collisions;
        continue;
repeat_failed:
        goto done;
    }

    printf("file %s\n",filename);
    printf("maximum_input_length %zu\n",maximum_length);
    printf("edits %zu\nresidues %zu\nrepeats %zu\n",edits,d,repeats);
    printf("average_pairs %.3f\n",pairs_sum/repeats);
    printf("average_wall_seconds %.6f\n",wall_sum/repeats);
    printf("average_cpu_seconds %.6f\n",cpu_sum/repeats);
    printf("any_seed_collision_probability %.9f\n",
           pairs_sum?any_sum/pairs_sum:0.0);
    printf("all_compatible_seeds_collision_probability %.9f\n",
           pairs_sum?all_sum/pairs_sum:0.0);
    for (left_split=0;left_split<per;++left_split)
        for (right_split=0;right_split<per;++right_split) {
            size_t index=left_split*per+right_split;
            if (denominators[index])
                printf("left_Es%zu_Ed%zu_right_Es%zu_Ed%zu_collision_probability %.9f\n",
                       left_split,edits-left_split,right_split,
                       edits-right_split,collisions[index]/denominators[index]);
        }
    status=0;
done:
    free(denominators); free(collisions); free(line); fclose(input);
    return status;
}

int ssh3_collisions_main(int argc,char **argv) {
    size_t edits,d,repeats,maximum_length,pairs;
    uint64_t seed=(uint64_t)time(NULL)^((uint64_t)clock()<<32);
    SubseqHash3TableFile tables={0};
    char *generated_path=NULL;
    int result=EXIT_FAILURE;
    if (argc==2 && (!strcmp(argv[1],"--help") || !strcmp(argv[1],"-h"))) {
        help(stdout); return EXIT_SUCCESS;
    }
    if ((argc!=5 && argc!=6) || parse_size(argv[2],&edits) ||
        parse_size(argv[3],&d) || parse_size(argv[4],&repeats) ||
        !d || !repeats || edits==SIZE_MAX) {
        help(stderr); return EXIT_FAILURE;
    }
    if (scan_input(argv[1],edits,&maximum_length,&pairs)) goto done;
    if (argc==5) {
        size_t length=strlen(argv[1]);
        generated_path=malloc(length+5);
        if (!generated_path) goto done;
        memcpy(generated_path,argv[1],length);
        memcpy(generated_path+length,".ABC",5);
        if (subseqhash3_table_file_create(
                generated_path,maximum_length,d,repeats,seed)) goto done;
    }
    if (subseqhash3_table_file_open(&tables,argc==6?argv[5]:generated_path))
        goto done;
    if (tables.n<maximum_length || tables.d!=d || tables.count!=repeats) {
        fputs("table file dimensions do not match the input, D, and R\n",stderr);
        goto done;
    }
    fprintf(stderr,"backend: %s\n",subseqhash3_backend());
    result=estimate(argv[1],edits,d,repeats,maximum_length,pairs,&tables)?
        EXIT_FAILURE:EXIT_SUCCESS;
done:
    free(generated_path); subseqhash3_table_file_close(&tables);
    return result;
}
