#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <ctype.h>
#include <time.h>
#include <pthread.h>
#include <unistd.h>
#include "mem.h"
#include "seq.h"
#include "pos.h"
#include "time.h"
#include "defs.h"
#include "param.h"
#include "msg.h"
#include "parser.h"
#include "reads.h"
#include "buffer.h"
#include "common.h"
#include "rmodel.h"

HASH *Hash; // HASH MEMORY IS PUBLIC

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - S O R T I N G - - - - - - - - - - - - - - - -
void MaxHeapify(int64_t a[], int i, int size){
  int64_t tmp, largest;
  int64_t l = (2 * i) + 1;
  int64_t r = (2 * i) + 2;
  if((l <= size) && (a[l] > a[i])) largest = l;
  else                             largest = i;
  if((r <= size) && (a[r] > a[largest])) largest = r;
  if(largest != i){
    tmp        = a[i];
    a[i]       = a[largest];
    a[largest] = tmp;
    MaxHeapify(a, largest, size);
    }
  }

void BuildMaxHeap(int64_t a[], int size){
  int i;
  for(i = size/2; i >= 0; i--)
    MaxHeapify(a, i, size);
  }

void HeapSort(int64_t a[], int size){
  int64_t tmp;
  int i;
  BuildMaxHeap(a, size);
  for(i = size; i > 0; i--){
    tmp  = a[i];
    a[i] = a[0];
    a[0] = tmp;
    --size;
    MaxHeapify(a, 0, size);
    }
  }

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - C O M P R E S S I N G - - - - - - - - - - - - - 
void MapTarget(void){
  uint64_t    nBase = 0, r = 0, nSymbol, initNSymbol;
  uint32_t    n, k, idxPos;
  PARSER      *PA = CreateParser();
  CBUF        *symBuf = CreateCBuffer(BUFFER_SIZE, BGUARD);
  uint8_t     sym, *pos;
  RMODEL      *RM = CreateRModel(P->minimum, P->kmer);

  FileType(PA, stdin);
  if(PA->type != 2){
    fprintf(stderr, "Error: input file must be in FASTQ format!\n");
    exit(1);
    }

  srand(0);

  Read *Read = CreateRead(10000, 40000);
  while((Read = GetRead(stdin, Read)) != NULL){

    nBase = strlen(Read->bases) - 1; // IT ALSO LOADS '\n' AT THE END
    int64_t positions[nBase+1], idx = 0;
    uint64_t base = 0;

    for(idxPos = 0 ; idxPos < nBase ; ++idxPos){

      sym = Read->bases[idxPos];
      if(sym == 'N') sym = rand() % 4; 
      else           sym = DNASymToNum(sym);
      symBuf->buf[symBuf->idx] = sym;

      n = 0;
      pos = &symBuf->buf[symBuf->idx-1];
      GetIdxRM(pos, RM);
      if(++base > P->kmer)
        positions[idx++] = GetPositionRM(RM, Hash);
      UpdateCBuffer(symBuf);
      }

    // TODO: EVALUATE BASES -> FUZZY MODE

    HeapSort(positions, idx-1);

    int x;
    for (x = 0; x < idx; ++x)
      printf("%li,", positions[x]);
    printf("\n");
    
    PA->nRead++;
    }

  if(P->verbose)
    ; //fprintf(stderr, "Processed %"PRI64" reads.\n");

  RemoveRModel(RM);
  RemoveCBuffer(symBuf);
  RemoveParser(PA);
  }

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - R E F E R E N C E - - - - - - - - - - - - - -
void LoadReference(void){
  FILE      *Reader = Fopen(P->reference, "r");
  uint32_t  n;
  uint64_t  idx = 0;
  uint64_t  k, idxPos, pos = 0, i = 0;
  PARSER    *PA = CreateParser();
  CBUF      *symBuf  = CreateCBuffer(BUFFER_SIZE, BGUARD);
  uint8_t   *readBuf = Calloc(BUFFER_SIZE, sizeof(uint8_t));
  uint8_t   sym;
  FileType(PA, Reader);
  uint64_t  file_length = NBytesInFile(Reader);
  RMODEL    *RM = CreateRModel(P->minimum, P->kmer);

  while((k = fread(readBuf, 1, BUFFER_SIZE, Reader)))
    for(idxPos = 0 ; idxPos < k ; ++idxPos){

      if(ParseSym(PA, (sym = readBuf[idxPos])) == -1){
        idx = 0;
        continue;
        }
 
      symBuf->buf[symBuf->idx] = sym = DNASymToNum(sym);
  
      if(sym != 4){
        symBuf->buf[symBuf->idx] = sym;
        GetIdxRM(symBuf->buf+symBuf->idx-1, RM);
        UpdateRM(RM, Hash, pos++);
        UpdateCBuffer(symBuf);
        }

      CalcProgress(file_length, ++i);
      }

  RemoveRModel(RM);
  RemoveCBuffer(symBuf);
  Free(readBuf);
  RemoveParser(PA);
  fclose(Reader);
  }

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - - M A P P E R - - - - - - - - - - - - - - - -
void MapAction(){
  uint32_t n;
  pthread_t t[P->nThreads];
  Threads  *T = (Threads *) Calloc(P->nThreads, sizeof(Threads));
  for(n = 0 ; n < P->nThreads ; ++n) T[n].id = n; 

  fprintf(stderr, "  [+] Building hash ...\n");
  Hash = CreateHash();
  fprintf(stderr, "      Done!                \n");

  fprintf(stderr, "  [+] Loading reference ...\n");
  LoadReference();
  fprintf(stderr, "      Done!                \n");

  fprintf(stderr, "  [+] Map contigs ... \n");
  MapTarget();
  fprintf(stderr, "\r      Done!                   \n");
  }

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - M A I N - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
int32_t main(int argc, char *argv[]){
  char **p = *&argv;

  P = (Parameters *) Malloc(1 * sizeof(Parameters));
  if((P->help = ArgsState(DEF_HELP, p, argc, "-h")) == 1 || argc < 2){
    PrintMenu();
    return EXIT_SUCCESS;
    }

  if(ArgsState(DEF_VERSION, p, argc, "-V")){
    PrintVersion();
    return EXIT_SUCCESS;
    }

  P->verbose    = ArgsState (DEF_VERBOSE, p, argc, "-v" );
  P->force      = ArgsState (DEF_FORCE,   p, argc, "-F" );
  P->kmer       = ArgsNum   (DEF_KMER,    p, argc, "-k", MIN_KMER, MAX_KMER);
  P->minimum    = ArgsNum   (DEF_MINI,    p, argc, "-m", MIN_MINI, MAX_MINI);
  P->reference  = argv[argc-1];

  fprintf(stderr, "\n");
  if(P->verbose) PrintArgs(P);

  fprintf(stderr, "==[ PROCESSING ]====================\n");
  TIME *Time = CreateClock(clock());
  MapAction();
  StopTimeNDRM(Time, clock());
  fprintf(stderr, "\n");

  fprintf(stderr, "==[ STATISTICS ]====================\n");
  StopCalcAll(Time, clock());
  fprintf(stderr, "\n");
  RemoveClock(Time);

  return EXIT_SUCCESS;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
