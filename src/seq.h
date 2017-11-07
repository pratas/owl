#ifndef SEQ_H_INCLUDED
#define SEQ_H_INCLUDED

#include "defs.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

typedef struct{
  uint8_t  *buf;
  uint64_t size;
  uint64_t init;
  uint64_t idx;
  }
SEQ;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

SEQ *CreateSeq(uint32_t);
void UpdateSeq(SEQ *, uint8_t);
void RemoveSeq(SEQ *);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif
