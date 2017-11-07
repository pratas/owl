#ifndef LINES_H_INCLUDED
#define LINES_H_INCLUDED

#include "defs.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

typedef struct{
  // CONTIGS
  char     contigs_name[MAX_STR];
  int64_t  contigs_relative_init_pos;
  int64_t  contigs_relative_end_pos;
  int64_t  contigs_absolute_init_pos;
  int64_t  contigs_absolute_end_pos;
  // REFERENCE
  char     reference_name[MAX_STR];
  int64_t  reference_relative_init_pos;
  int64_t  reference_relative_end_pos;
  int64_t  reference_absolute_init_pos;
  int64_t  reference_absolute_end_pos;
  }
LINE;

typedef struct{
  LINE     *Lines;
  uint64_t idx;
  uint64_t max;
  }
LCACHE;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

LCACHE        *CreateLCache        (size_t);
void          PrintLine            (LCACHE *, FILE *, uint32_t);
void          ResetChar2Bar0       (char *);
void          UpdateLCacheIdxInit  (LCACHE *);
void          UpdateLCacheIdx      (LCACHE *);
void          RemoveLCache         (LCACHE *);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif
