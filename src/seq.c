#include "seq.h"
#include "mem.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

SEQ *CreateSeq(uint32_t size){
  SEQ *Sequence  = (SEQ *) Calloc(1, sizeof(SEQ));
  Sequence->size = size;
  Sequence->init = size;
  Sequence->idx  = 0;
  Sequence->buf  = (uint8_t *) Calloc(size, sizeof(uint8_t));
  return Sequence;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void UpdateSeq(SEQ *Sequence, uint8_t sym){
  if(Sequence->idx == Sequence->size)
    Sequence->buf = (uint8_t *) Realloc(Sequence->buf, (Sequence->size += 
    Sequence->init) * sizeof(uint8_t));
  Sequence->buf[Sequence->idx++] = sym;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void RemoveSeq(SEQ *Sequence){
  Free(Sequence->buf);
  Free(Sequence);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
