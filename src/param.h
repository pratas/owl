#ifndef PARAM_H_INCLUDED
#define PARAM_H_INCLUDED

#include "defs.h"
#include "rmodel.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

typedef struct{
  U8       help;
  U8       verbose;
  U8       force;
  U8       order;
  U8       link;
  U32      kmer;
  U32      minimum;
  char     *reference;
  }
Parameters;

Parameters *P;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif
