#define _POSIX_C_SOURCE 200809L
#include "ssh3_seeding.h"
#include "ssh3_cli.h"
#include "ssh3_collisions.h"

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef struct {
    uint64_t *bits;
    long double paths;
} EditDfaState;

typedef struct {
    const char *text;
    size_t n,es,ed,cells,words;
} EditDfa;

static size_t edit_cell(const EditDfa *a,size_t p,size_t d,size_t s) {
    return (p*(a->ed+1)+d)*(a->es+1)+s;
}

static void bit_set(uint64_t *bits,size_t i) { bits[i>>6]|=UINT64_C(1)<<(i&63); }
static int bit_get(const uint64_t *bits,size_t i) { return (int)((bits[i>>6]>>(i&63))&1); }

/* Add all states reachable by deleting bases. */
static void deletion_closure(const EditDfa *a,uint64_t *bits) {
    size_t p,d,s;
    for (p=0;p<a->n;++p) for (d=0;d<a->ed;++d) for (s=0;s<=a->es;++s)
        if (bit_get(bits,edit_cell(a,p,d,s)))
            bit_set(bits,edit_cell(a,p+1,d+1,s));
}

static int dfa_step(const EditDfa *a,const uint64_t *from,unsigned symbol,
                    uint64_t *to) {
    size_t p,d,s; int nonempty=0;
    memset(to,0,a->words*sizeof(*to));
    for (p=0;p<a->n;++p) for (d=0;d<=a->ed;++d) for (s=0;s<=a->es;++s) {
        size_t next_s;
        unsigned original;
        if (!bit_get(from,edit_cell(a,p,d,s))) continue;
        original=a->text[p]=='A'||a->text[p]=='a'?0:
                 a->text[p]=='C'||a->text[p]=='c'?1:
                 a->text[p]=='G'||a->text[p]=='g'?2:3;
        next_s=s+(symbol!=original);
        if (next_s<=a->es) {
            bit_set(to,edit_cell(a,p+1,d,next_s)); nonempty=1;
        }
    }
    if (nonempty) deletion_closure(a,to);
    return nonempty;
}

static int add_state(EditDfaState **states,size_t *count,size_t *capacity,
                     const uint64_t *bits,size_t words,long double paths) {
    size_t i;
    for (i=0;i<*count;++i)
        if (!memcmp((*states)[i].bits,bits,words*sizeof(*bits))) {
            (*states)[i].paths+=paths; return 0;
        }
    if (*count==*capacity) {
        size_t next=*capacity?*capacity*2:8;
        EditDfaState *grown=realloc(*states,next*sizeof(**states));
        if (!grown) return -1;
        *states=grown; *capacity=next;
    }
    (*states)[*count].bits=malloc(words*sizeof(*bits));
    if (!(*states)[*count].bits) return -1;
    memcpy((*states)[*count].bits,bits,words*sizeof(*bits));
    (*states)[*count].paths=paths; ++*count; return 0;
}

static void free_states(EditDfaState *states,size_t count) {
    size_t i; for (i=0;i<count;++i) free(states[i].bits); free(states);
}

/* Count the language of a, or its intersection with b when b is non-NULL. */
static int count_language(const EditDfa *a,const EditDfa *b,long double *answer) {
    EditDfaState *current=NULL,*next=NULL; size_t nc=0,cc=0,nn,cn,i,k,length=a->n-a->ed;
    size_t pair_words=a->words+(b?b->words:0),step_words=pair_words;
    uint64_t *initial=calloc(pair_words,sizeof(*initial)),*stepped=calloc(step_words,sizeof(*stepped));
    if (!initial || !stepped) goto fail;
    bit_set(initial,edit_cell(a,0,0,0)); deletion_closure(a,initial);
    if (b) { bit_set(initial+a->words,edit_cell(b,0,0,0)); deletion_closure(b,initial+a->words); }
    if (add_state(&current,&cc,&nc,initial,pair_words,1.0L)) goto fail;
    for (k=0;k<length;++k) {
        nn=cn=0; next=NULL;
        for (i=0;i<cc;++i) for (unsigned x=0;x<4;++x) {
            if (!dfa_step(a,current[i].bits,x,stepped)) continue;
            if (b && !dfa_step(b,current[i].bits+a->words,x,stepped+a->words)) continue;
            if (add_state(&next,&nn,&cn,stepped,pair_words,current[i].paths)) { free_states(next,nn); goto fail; }
        }
        free_states(current,cc); current=next; cc=nn; nc=cn;
    }
    *answer=0;
    for (i=0;i<cc;++i) {
        int accept=bit_get(current[i].bits,edit_cell(a,a->n,a->ed,a->es));
        if (b) accept&=bit_get(current[i].bits+a->words,edit_cell(b,b->n,b->ed,b->es));
        if (accept) *answer+=current[i].paths;
    }
    free(initial); free(stepped); free_states(current,cc); return 0;
fail:
    free(initial); free(stepped); free_states(current,cc); return -1;
}

static int init_dfa(EditDfa *a,const char *text,size_t es,size_t ed) {
    size_t n=strlen(text),i;
    if (ed>n || es>n-ed) return -1;
    for (i=0;i<n;++i) if (!strchr("ACGTacgt",text[i])) return -1;
    if (n+1>SIZE_MAX/(ed+1) || (n+1)*(ed+1)>SIZE_MAX/(es+1)) return -1;
    a->text=text; a->n=n; a->es=es; a->ed=ed;
    a->cells=(n+1)*(ed+1)*(es+1); a->words=(a->cells+63)/64; return 0;
}

double ssh3_edit_neighborhood_jaccard(const char *s,const char *t,size_t es,size_t ed) {
    EditDfa a,b; long double na,nb,intersection,union_size;
    if (!s || !t || init_dfa(&a,s,es,ed) || init_dfa(&b,t,es,ed)) return NAN;
    if (count_language(&a,NULL,&na) || count_language(&b,NULL,&nb)) return NAN;
    if (a.n-a.ed!=b.n-b.ed) intersection=0;
    else if (count_language(&a,&b,&intersection)) return NAN;
    union_size=na+nb-intersection;
    return union_size?(double)(intersection/union_size):1.0;
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


typedef struct {
    size_t count;
    uint8_t **seed;
    size_t *length;
} SeedBatch;
static void batch_free(SeedBatch *a) {
    if (a->seed) for (size_t i=0;i<a->count;++i) free(a->seed[i]);
    free(a->seed); free(a->length); memset(a,0,sizeof(*a));
}
static int batch_copy(SeedBatch *a,size_t i,const uint8_t *seed,size_t length) {
    a->seed[i]=malloc(length); if (!a->seed[i]) return -1;
    memcpy(a->seed[i],seed,length); a->length[i]=length; return 0;
}
static int collect_sub(size_t start,const SubseqHash3EditResult *seeds,size_t count,void *ctx) {
    SeedBatch *a=ctx; (void)start;
    for (size_t i=0;i<count;++i) if (batch_copy(a,i,seeds[i].seed,seeds[i].length)) return -1;
    return 0;
}
static int collect_ind(size_t start,const SubseqHash3DeletionResult *seeds,size_t count,void *ctx) {
    SeedBatch *a=ctx; (void)start;
    for (size_t i=0;i<count;++i) if (batch_copy(a,i,seeds[i].seed,seeds[i].length)) return -1;
    return 0;
}
static int batch_seed(SeedBatch *a,const char *dna,size_t length,size_t e,
    const SubseqHash3PivotTables *tables,size_t r,size_t sets,int ind,SubseqHash3Traversal traversal) {
    a->count=r*sets;
    a->seed=calloc(a->count,sizeof(*a->seed)); a->length=calloc(a->count,sizeof(*a->length));
    uint8_t *encoded=malloc(length);
    if (!a->seed || !a->length || !encoded) { free(encoded); return -1; }
    static const char alphabet[]="ACGTacgt";
    for (size_t i=0;i<length;++i) encoded[i]=(uint8_t)((strchr(alphabet,dna[i])-alphabet)%4);
    int status=ind?subseqhash3_seeds_ind(encoded,length,length,e,tables,r,sets,traversal,collect_ind,a):
        subseqhash3_seeds_sub(encoded,length,length,e,tables,r,sets,traversal,collect_sub,a);
    free(encoded); return status;
}
static int equal_seed(const SeedBatch *a,size_t i,const SeedBatch *b,size_t j) {
    return a->length[i]==b->length[j] && !memcmp(a->seed[i],b->seed[j],a->length[i]);
}
static void help(FILE *stream) {
    fputs("Evaluate pivot seeds using the same R pivots and S sets as seeding.\n\n"
        "Usage:\n  subseqhash3 collisions sub|ind PAIRS_FILE E D R S [TABLE_FILE]\n"
        "                         [--inward] [--seed SEED] [--jaccard]\n\n"
        "  sub        One deletion and exactly E substitutions per string.\n"
        "  ind        Exactly E deletions and no substitutions per string.\n"
        "  PAIRS_FILE One ACGT_STRING1,ACGT_STRING2 pair per line; no header.\n"
        "  R          Leading pivots per set, valid for BOTH strings in every pair.\n"
        "  S          Independent table sets; each string produces R*S seeds.\n"
        "  TABLE_FILE Optional pivot tables with matching D and at least S sets.\n"
        "  --inward   Alternative traversal; default is outward.\n"
        "  --seed     Reproducible automatic tables, stored in PAIRS_FILE.piv.\n"
        "  --jaccard  Also calculate exact edit-neighborhood Jaccard (can be slow).\n\n"
        "Whole strings are seeded as single windows. Different seed lengths do not\n"
        "match. Report any seed-pair intersection and aligned set/pivot matches.\n",stream);
}

int ssh3_collisions_main(int argc,char **argv) {
    const char *args[7]; size_t nargs=0,max_length=0,pairs=0,capacity=0;
    uint64_t v[4],seed=(uint64_t)time(NULL)^((uint64_t)clock()<<32);
    int inward=0,jaccard=0,seed_set=0,status=EXIT_FAILURE;
    FILE *input=NULL; char *line=NULL,*generated=NULL;
    size_t *aligned=NULL,any=0,any_aligned=0,all_aligned=0;
    double jac_sum=0,jac_sq=0;
    SubseqHash3PivotTableFile file={0};
    if (argc==2 && (!strcmp(argv[1],"--help") || !strcmp(argv[1],"-h"))) { help(stdout); return 0; }
    for (int i=1;i<argc;++i) {
        if (!strcmp(argv[i],"--inward")) { if (inward) goto invalid; inward=1; }
        else if (!strcmp(argv[i],"--jaccard")) { if (jaccard) goto invalid; jaccard=1; }
        else if (!strcmp(argv[i],"--seed")) {
            if (seed_set || i+1>=argc || ssh3_cli_number(argv[++i],&seed)) goto invalid;
            seed_set=1;
        } else {
            if (!strncmp(argv[i],"--",2) || nargs==7) goto invalid;
            args[nargs++]=argv[i];
        }
    }
    if (nargs<6 || (strcmp(args[0],"sub") && strcmp(args[0],"ind")) || (nargs==7 && seed_set)) goto invalid;
    for (size_t i=0;i<4;++i) if (ssh3_cli_number(args[i+2],v+i) || v[i]>SIZE_MAX) goto invalid;
    const int ind=!strcmp(args[0],"ind");
    const size_t e=(size_t)v[0],d=(size_t)v[1],r=(size_t)v[2],sets=(size_t)v[3];
    if (!d || !r || !sets || sets>SIZE_MAX/r || r*sets>SIZE_MAX/sizeof(*aligned)) goto invalid;
    input=fopen(args[1],"r"); if (!input) { perror(args[1]); goto done; }
    ssize_t got;
    /* Validate every pair and dimension before generating tables or emitting statistics. */
    while ((got=getline(&line,&capacity,input))>=0) {
        char *left,*right; size_t n1,n2;
        if (parse_pair(line,(size_t)got,&left,&n1,&right,&n2)) goto bad_pair;
        if (n1<(ind?1u:2u) || n2<(ind?1u:2u) || n1>(ind?256u:257u) || n2>(ind?256u:257u) ||
            e>=n1 || e>=n2 || r>(ind?n1-e:n1-1) || r>(ind?n2-e:n2-1)) goto bad_pair;
        if (n1>max_length) max_length=n1;
        if (n2>max_length) max_length=n2;
        ++pairs; continue;
 bad_pair:
        fprintf(stderr,"invalid DNA, edit budget or R on pair line %zu\n",pairs+1); goto done;
    }
    if (ferror(input) || !pairs) { fputs("pair file is empty or unreadable\n",stderr); goto done; }
    if (ssh3_cli_tables(&file,nargs==7?args[6]:NULL,&generated,args[1],args[1],max_length,d,sets,
        ind?max_length-e:max_length-1,seed)) goto done;
    aligned=calloc(r*sets,sizeof(*aligned)); if (!aligned) goto done;
    rewind(input); clearerr(input);
    size_t seen=0;
    while ((got=getline(&line,&capacity,input))>=0) {
        char *left,*right; size_t n1,n2; SeedBatch a={0},b={0};
        if (parse_pair(line,(size_t)got,&left,&n1,&right,&n2)) goto done;
        if (batch_seed(&a,left,n1,e,file.tables,r,sets,ind,inward?SUBSEQHASH3_INWARD:SUBSEQHASH3_OUTWARD) ||
            batch_seed(&b,right,n2,e,file.tables,r,sets,ind,inward?SUBSEQHASH3_INWARD:SUBSEQHASH3_OUTWARD)) {
            batch_free(&a); batch_free(&b); goto done;
        }
        int some=0,all=1,intersection=0;
        for (size_t i=0;i<r*sets;++i) {
            const int equal=equal_seed(&a,i,&b,i);
            aligned[i]+=(size_t)equal; some|=equal; all&=equal;
            for (size_t j=0;!intersection && j<r*sets;++j) intersection|=equal_seed(&a,i,&b,j);
        }
        any+=(size_t)intersection; any_aligned+=(size_t)some; all_aligned+=(size_t)all;
        batch_free(&a); batch_free(&b);
        if (jaccard) {
            const double jac=ssh3_edit_neighborhood_jaccard(left,right,ind?0:e,ind?e:1);
            if (!isfinite(jac)) goto done;
            jac_sum+=jac; jac_sq+=jac*jac;
        }
        ++seen;
    }
    if (ferror(input) || seen!=pairs) goto done;
    printf("mode %s\ntraversal %s\npairs %zu\nresidues %zu\npivots_per_set %zu\nsets %zu\nseeds_per_string %zu\n",
        args[0],inward?"inward":"outward",pairs,d,r,sets,r*sets);
    printf("any_seed_collision_probability %.9f\n",(double)any/pairs);
    printf("any_same_set_pivot_collision_probability %.9f\n",(double)any_aligned/pairs);
    printf("all_same_set_pivot_collision_probability %.9f\n",(double)all_aligned/pairs);
    for (size_t i=0;i<r*sets;++i)
        printf("set_%zu_pivot_%zu_collision_probability %.9f\n",i/r,i%r,(double)aligned[i]/pairs);
    if (jaccard) {
        const double mean=jac_sum/pairs;
        printf("average_jaccard %.9f\njaccard_std %.9f\n",mean,sqrt(fmax(0,jac_sq/pairs-mean*mean)));
    }
    status=ferror(stdout)?EXIT_FAILURE:EXIT_SUCCESS;
done:
    if (input && fclose(input)) status=EXIT_FAILURE;
    free(line); free(generated); free(aligned); subseqhash3_pivot_table_file_close(&file); return status;
invalid:
    help(stderr); return EXIT_FAILURE;
}
