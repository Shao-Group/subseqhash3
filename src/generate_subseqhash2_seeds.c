#define _POSIX_C_SOURCE 200809L
#include "subseqhash2_edits_i8.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int parse_size(const char *text,size_t *value) {
    char *end; unsigned long long parsed;
    errno=0; parsed=strtoull(text,&end,10);
    if (errno || !*text || *end || parsed>(unsigned long long)SIZE_MAX)
        return -1;
    *value=(size_t)parsed; return 0;
}

static int encode_sequence(const char *line,size_t length,uint8_t **encoded) {
    uint8_t *s=malloc(length?length:1); size_t i;
    if (!s) return -1;
    for (i=0;i<length;++i) {
        switch (line[i]) {
        case 'A': case 'a': s[i]=0; break;
        case 'C': case 'c': s[i]=1; break;
        case 'G': case 'g': s[i]=2; break;
        case 'T': case 't': s[i]=3; break;
        default:
            free(s); return -1;
        }
    }
    *encoded=s; return 0;
}

static int write_result(FILE *out,size_t sequence_index,size_t window_start,
                        size_t order,const SubseqHashEditResult *r) {
    static const char alphabet[]="ACGT"; size_t i;
    if (fprintf(out,"%zu\t%zu\t%zu\t%zu\t%zu\t%zu\t",
                sequence_index,window_start,r->substitutions,r->deletions,
                order,r->length)<0) return -1;
    for (i=0;i<r->length;++i)
        if (fputc(alphabet[r->seed[i]],out)==EOF) return -1;
    return fprintf(out,"\t%zu\t%.0f\n",r->psi,r->omega)<0?-1:0;
}

static int process_sequence(FILE *out,const char *line,size_t length,
                            size_t sequence_index,size_t W,size_t E,
                            size_t mode,const SubseqHash2I8Tables *tables) {
    uint8_t *encoded=NULL; SubseqHashEditResult *results=NULL;
    size_t windows,per,total,w,within=0,s,order,index=0; int status=-1;
    if (length<W) {
        fprintf(stderr,"sequence %zu has length %zu, shorter than W=%zu\n",
                sequence_index,length,W);
        return -1;
    }
    if (encode_sequence(line,length,&encoded)) {
        fprintf(stderr,"sequence %zu contains a character outside A/C/G/T\n",
                sequence_index);
        return -1;
    }
    windows=length-W+1; per=subseqhash2_i8_results_per_window(W,E);
    if (per==SIZE_MAX || windows>SIZE_MAX/per ||
        windows*per>SIZE_MAX/sizeof(*results)) goto done;
    total=windows*per;
    results=calloc(total?total:1,sizeof(*results));
    if (!results) goto done;
    if (subseqhash2_compute_long_edits_i8(
            encoded,length,W,E,tables,results)) goto done;
    for (w=0;w<windows;++w) {
        within=0;
        for (s=0;s<=E;++s) {
            const size_t K=W-E+s;
            if (mode==1) {
                for (order=0;order<K;++order,++within,++index)
                    if (write_result(out,sequence_index,w,order,
                                     results+index)) goto done;
            } else {
                size_t best=0;
                for (order=1;order<K;++order)
                    if (results[index+order].omega>results[index+best].omega)
                        best=order;
                if (write_result(out,sequence_index,w,best,
                                 results+index+best)) goto done;
                index+=K; within+=K;
            }
        }
        if (within!=per) goto done;
    }
    status=0;
done:
    if (results) {
        subseqhash_edit_results_destroy(results,total);
        free(results);
    }
    free(encoded); return status;
}

int main(int argc,char **argv) {
    FILE *input=NULL,*output=NULL; char *line=NULL; size_t capacity=0;
    ssize_t got; size_t W,E,D,mode,sequence_index=0; int status=EXIT_FAILURE;
    SubseqHash2I8Tables tables;
    memset(&tables,0,sizeof(tables));
    if (argc!=7 || parse_size(argv[2],&W) || parse_size(argv[3],&E) ||
        parse_size(argv[4],&D) || parse_size(argv[5],&mode) || !W || !D ||
        E>W || (mode!=1 && mode!=2)) {
        fprintf(stderr,"usage: %s INPUT W E D MODE OUTPUT\n",argv[0]);
        return EXIT_FAILURE;
    }
    input=fopen(argv[1],"r");
    if (!input) { perror(argv[1]); goto done; }
    output=fopen(argv[6],"w");
    if (!output) { perror(argv[6]); goto done; }
    if (subseqhash2_i8_tables_init(
            &tables,W,D,UINT64_C(0x7375627365716832))) {
        fputs("failed to initialize SubseqHash2 tables\n",stderr); goto done;
    }
    if (fputs("sequence\twindow_start\tsubstitutions\tdeletions\torder"
              "\tseed_length\tseed\tpsi\tomega\n",output)==EOF)
        goto done;
    while ((got=getline(&line,&capacity,input))>=0) {
        size_t length=(size_t)got;
        while (length && (line[length-1]=='\n'||line[length-1]=='\r')) --length;
        if (!length) {
            fprintf(stderr,"sequence %zu is empty\n",sequence_index); goto done;
        }
        if (process_sequence(
                output,line,length,sequence_index,W,E,mode,&tables))
            goto done;
        ++sequence_index;
    }
    if (ferror(input)) { perror(argv[1]); goto done; }
    if (fflush(output)) { perror(argv[6]); goto done; }
    status=EXIT_SUCCESS;
done:
    free(line); subseqhash2_i8_tables_destroy(&tables);
    if (input && fclose(input)) status=EXIT_FAILURE;
    if (output && fclose(output)) status=EXIT_FAILURE;
    return status;
}
