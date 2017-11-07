#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "lines.h"
#include "defs.h"
#include "mem.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

LCACHE *CreateLCache(size_t size){
  LCACHE *LC = (LCACHE *) Calloc(1, sizeof(LCACHE));
  LC->max    = size;
  LC->idx    = 0;
  LC->Lines  = (LINE *) Calloc(LC->max+1, sizeof(LINE));
  return LC;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void PrintLine(LCACHE *LC, FILE *F, uint32_t idx){
  fprintf(F, "%s\t%"PRIi64"\t%"PRIi64"\t%"PRIi64"\t%"PRIi64"\t%s\t"
                 "%"PRIi64"\t%"PRIi64"\t%"PRIi64"\t%"PRIi64"\n",
          LC->Lines[idx].contigs_name, 
          LC->Lines[idx].contigs_relative_init_pos,
          LC->Lines[idx].contigs_relative_end_pos,
          LC->Lines[idx].contigs_absolute_init_pos,
          LC->Lines[idx].contigs_absolute_end_pos,
          LC->Lines[idx].reference_name, 
          LC->Lines[idx].reference_relative_init_pos,
          LC->Lines[idx].reference_relative_end_pos, 
          LC->Lines[idx].reference_absolute_init_pos,
          LC->Lines[idx].reference_absolute_end_pos);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void ResetChar2Bar0(char *str){
  uint32_t n;
  for(n = 0 ; n < MAX_STR ; ++n)
    str[n] = '\0';
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void UpdateLCacheIdxInit(LCACHE *LC){
  if(++LC->idx == LC->max){
    LINE *LI = &LC->Lines[0], *LE = &LC->Lines[LC->max];
    strcpy(LI->contigs_name, LE->contigs_name);
    LI->contigs_relative_init_pos   = LE->contigs_relative_init_pos;
    LI->contigs_relative_end_pos    = LE->contigs_relative_end_pos;
    LI->contigs_absolute_init_pos   = LE->contigs_absolute_init_pos;
    LI->contigs_absolute_end_pos    = LE->contigs_absolute_end_pos;
    strcpy(LI->reference_name, LE->reference_name);   
    LI->reference_relative_init_pos = LE->reference_relative_init_pos;
    LI->reference_relative_end_pos  = LE->reference_relative_end_pos;
    LI->reference_absolute_init_pos = LE->reference_absolute_init_pos;
    LI->reference_absolute_end_pos  = LE->reference_absolute_end_pos;
    LC->idx = 1; // SET IDX TO NEXT
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void UpdateLCacheIdx(LCACHE *LC){
  if(++LC->idx == LC->max)
    LC->idx = 1; // SET IDX TO NEXT
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void RemoveLCache(LCACHE *LC){
  uint32_t n;
  for(n = 0 ; n < LC->max ; ++n)
    Free(LC->Lines);
  Free(LC);
  } 

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
