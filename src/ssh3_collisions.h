#ifndef SSH3_COLLISIONS_H
#define SSH3_COLLISIONS_H

#include <stddef.h>

/* Return the Jaccard index of the distinct DNA strings obtainable from s and t
 * with exactly substitutions substitutions and deletions deletions.  A
 * substitution changes a retained base to one of the other three DNA bases.
 * NAN is returned for invalid DNA or impossible edit counts. */
double ssh3_edit_neighborhood_jaccard(const char *s,const char *t,
                                      size_t substitutions,size_t deletions);

#endif
