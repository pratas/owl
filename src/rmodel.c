#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include "defs.h"
#include "mem.h"
#include "common.h"
#include "rmodel.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// CALCULATION OF CONTEXT MULTIPLICATOR FOR INDEX FUNCTION USAGE
//
uint64_t CalcMult(uint32_t c){
  uint32_t n;
  uint64_t x[c], p = 1;
  for(n = 0 ; n < c ; ++n){
    x[n] = p;
    p <<= 2;
    }
  return x[c-1];
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// CREATES THE RMODEL BASIC STRUCTURE 
//
RMODEL *CreateRModel(uint32_t min, uint32_t k){

  RMODEL *R  = (RMODEL *)  Calloc(1,   sizeof(RMODEL));
  R->idx      = 0;
  R->kmer     = k;
  R->nBases   = 0;
  R->minSize  = min;
  R->pos      = 0;
  R->mult     = CalcMult(k);

  return R;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// INDEX CALC BASED ON PAST SYMBOLS FOR REPEAT MODEL
//
uint64_t GetIdxRM(uint8_t *p, RMODEL *R){
  return (R->idx = ((R->idx-*(p-R->kmer)*R->mult)<<2)+*p);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// RESET INDEX IN REPEAT MODEL
//
void ResetIdxRM(RMODEL *R){
  R->idx = 0;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// UPDATE REPEAT MODEL
//
void UpdateRM(RMODEL *R, HASH *H, uint64_t pos){
  InsertKmerPos(H, R->idx, pos);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// GET THE POSITION IN THE HASH
//
int64_t GetPositionRM(RMODEL *R, HASH *H){
  return GetHashPosition(H, R->idx);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// REMOVE PERMANENTLY RMODEL
//
void RemoveRModel(RMODEL *R){
  Free(R);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
