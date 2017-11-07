#ifndef RMODEL_H_INCLUDED
#define RMODEL_H_INCLUDED

#include "defs.h"
#include "buffer.h"
#include "hash.h"
#include "pos.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// REPEAT MODELS TO HANDLE LONG SEGMENTS. DATA SUBSTITUTIONS DO NOT AFFECT THE
// PERFORMANCE SO MUCH AS IN CONTEXT MODELS.
//

typedef struct{
  uint32_t kmer;     // CONTEXT TEMPLATE SIZE FOR REPEAT MODEL
  uint32_t minSize;  // MINIMUM BLOCK SIZE
  uint64_t mult;     // INDEX MULTIPLIER
  uint64_t idx;      // CURRENT CONTEXT INDEX
  uint64_t pos;      // POSITION ALONG THE SEQUENCE
  uint64_t nBases;   // NUMBER OF ACGTN's
  }
RMODEL;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint64_t   CalcMult         (uint32_t);
RMODEL     *CreateRModel    (uint32_t, uint32_t);
uint64_t   GetIdxRM         (uint8_t *, RMODEL *);
void       UpdateRM         (RMODEL *, HASH *, uint64_t);
void       ResetIdxRM       (RMODEL *);
void       RemoveRModel     (RMODEL *);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif
