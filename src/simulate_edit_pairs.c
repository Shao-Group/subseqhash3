#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static uint64_t random_word(uint64_t *state) {
    uint64_t z=(*state+=UINT64_C(0x9e3779b97f4a7c15));
    z=(z^(z>>30))*UINT64_C(0xbf58476d1ce4e5b9);
    z=(z^(z>>27))*UINT64_C(0x94d049bb133111eb);
    return z^(z>>31);
}

static int parse_size(const char *text,size_t *value) {
    char *end; unsigned long long parsed;
    errno=0; parsed=strtoull(text,&end,10);
    if (errno || !*text || *end || parsed>(unsigned long long)SIZE_MAX)
        return -1;
    *value=(size_t)parsed;
    return 0;
}

static size_t minimum3(size_t a,size_t b,size_t c) {
    size_t m=a<b?a:b;
    return m<c?m:c;
}

static size_t edit_distance(const char *a,size_t n,const char *b,size_t m,
                            size_t *previous,size_t *current) {
    size_t i,j,*tmp;
    for (j=0;j<=m;++j) previous[j]=j;
    for (i=1;i<=n;++i) {
        current[0]=i;
        for (j=1;j<=m;++j)
            current[j]=minimum3(previous[j]+1,current[j-1]+1,
                                previous[j-1]+(a[i-1]!=b[j-1]));
        tmp=previous; previous=current; current=tmp;
    }
    return previous[m];
}

static void shuffle(size_t *values,size_t count,uint64_t *rng) {
    size_t i;
    for (i=count;i>1;--i) {
        const size_t j=(size_t)(random_word(rng)%i);
        const size_t tmp=values[i-1];
        values[i-1]=values[j]; values[j]=tmp;
    }
}

static int simulate(FILE *output,size_t pairs,size_t n,size_t substitutions,
                    size_t deletions,uint64_t seed) {
    static const char alphabet[]="ACGT";
    const size_t core_length=n-deletions;
    const size_t target_distance=substitutions+2*deletions;
    char *source=NULL,*core=NULL,*target=NULL;
    size_t *indices=NULL,*previous=NULL,*current=NULL;
    uint8_t *deleted=NULL,*inserted=NULL;
    size_t produced=0,attempts=0,i,j;
    int result=-1;

    source=malloc(n?n:1); core=malloc(core_length?core_length:1);
    target=malloc(n?n:1);
    indices=malloc((n?n:1)*sizeof(*indices));
    deleted=malloc(n?n:1); inserted=malloc(n?n:1);
    previous=malloc((n+1)*sizeof(*previous));
    current=malloc((n+1)*sizeof(*current));
    if (!source || !core || !target || !indices || !deleted || !inserted ||
        !previous || !current)
        goto done;

    while (produced<pairs) {
        if (++attempts==SIZE_MAX) goto done;
        for (i=0;i<n;++i) {
            source[i]=alphabet[random_word(&seed)&3u];
            indices[i]=i;
            deleted[i]=0;
            inserted[i]=0;
        }
        shuffle(indices,n,&seed);
        for (i=0;i<deletions;++i) deleted[indices[i]]=1;
        for (i=0,j=0;i<n;++i)
            if (!deleted[i]) core[j++]=source[i];

        for (i=0;i<core_length;++i) indices[i]=i;
        shuffle(indices,core_length,&seed);
        for (i=0;i<substitutions;++i) {
            const size_t p=indices[i];
            unsigned replacement=(unsigned)(random_word(&seed)%3u);
            const unsigned original=(unsigned)(strchr(alphabet,core[p])-alphabet);
            if (replacement>=original) ++replacement;
            core[p]=alphabet[replacement];
        }

        for (i=0;i<n;++i) indices[i]=i;
        shuffle(indices,n,&seed);
        for (i=0;i<deletions;++i) inserted[indices[i]]=1;
        for (i=0,j=0;i<n;++i)
            target[i]=inserted[i]?alphabet[random_word(&seed)&3u]:core[j++];

        if (edit_distance(source,n,target,n,previous,current)!=target_distance)
            continue;
        if (fwrite(source,1,n,output)!=n || fputc('\n',output)==EOF ||
            fwrite(target,1,n,output)!=n || fputc('\n',output)==EOF)
            goto done;
        ++produced;
    }
    result=0;
done:
    free(current); free(previous); free(inserted); free(deleted); free(indices);
    free(target); free(core); free(source);
    if (result)
        fprintf(stderr,"generation failed after %zu of %zu pairs\n",
                produced,pairs);
    else
        fprintf(stderr,"generated %zu pairs in %zu attempts\n",
                produced,attempts);
    return result;
}

int main(int argc,char **argv) {
    size_t pairs,n,substitutions,deletions;
    FILE *output=stdout;
    uint64_t seed;
    int result;
    if ((argc!=5 && argc!=6) || parse_size(argv[1],&pairs) ||
        parse_size(argv[2],&n) || parse_size(argv[3],&substitutions) ||
        parse_size(argv[4],&deletions) || deletions>n ||
        substitutions>n-deletions ||
        deletions>(SIZE_MAX-substitutions)/2 ||
        substitutions+2*deletions>n) {
        fprintf(stderr,"usage: %s N n Es Ed [OUTPUT]\n",argv[0]);
        return EXIT_FAILURE;
    }
    if (argc==6) {
        output=fopen(argv[5],"w");
        if (!output) { perror(argv[5]); return EXIT_FAILURE; }
    }
    seed=UINT64_C(0x73696d756c617465) ^
         ((uint64_t)pairs<<32) ^ ((uint64_t)n<<16) ^
         ((uint64_t)substitutions<<8) ^ (uint64_t)deletions;
    result=simulate(output,pairs,n,substitutions,deletions,seed);
    if (output!=stdout && fclose(output)) result=-1;
    return result?EXIT_FAILURE:EXIT_SUCCESS;
}
