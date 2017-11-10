#include "msg.h"
#include <stdio.h>
#include <stdlib.h>

void PrintMenu(void){
  fprintf(stderr,
  "Usage: OWL [OPTIONS]... [FILE] [FILE]                                    \n"
  "A tool to sort FASTQ reads using cluster mapping.                        \n"
  "                                                                         \n"
  "Non-mandatory arguments:                                                 \n"
  "                                                                         \n"
  "  -h                         give this help,                             \n"
  "  -V                         display version number,                     \n"
  "  -v                         verbose mode (more information),            \n"
  "  -k <k-mer>                 k-mer size [1;20],                          \n"
  "  -m <minimum>               minimum block size,                         \n"
  "                                                                         \n"
  "Mandatory arguments:                                                     \n"
  "                                                                         \n"
  "  <FILE>                     reference file,                             \n"
  "                                                                         \n"
  "  <  <FILE>                  stdin input FASTQ file,                     \n"
  "  >  <FILE>                  stdout output sorted FASTQ file,            \n"
  "                                                                         \n"
  "Example:                                                                 \n"
  "                                                                         \n"
  "  ./OWL -v -k 16 -m 40 reference.fa < ex1.fq > ex1-sort.fq               \n"
  "                                                                         \n"
  "Report bugs to <pratas@ua.pt>.                                           \n");
  }

void PrintVersion(void){
  fprintf(stderr,
  "                                                                         \n"
  "                         =====================                           \n"
  "                         |      OWL %u.%u      |                         \n"
  "                         =====================                           \n"
  "                                                                         \n"
  "            A tool to sort FASTQ reads using cluster mapping             \n"
  "                                                                         \n"
  "Copyright (C) 2017-2018 University of Aveiro. This is a Free software.   \n"
  "You may redistribute copies of it under the terms of the GNU - General   \n"
  "Public License v3 <http://www.gnu.org/licenses/gpl.html>. There is not   \n"
  "ant warranty, to the extent permitted by law.                            \n" 
  "                                                                         \n"
  "           Developed and written by D. Pratas <pratas@ua.pt>.          \n\n", 
  VERSION, RELEASE);
  }

