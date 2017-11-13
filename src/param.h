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
  U8       delete;
  U8       header;
  U32      kmer;
  U32      minimum;
  char     *reference;
  char     *label;
  }
Parameters;

Parameters *P;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif
