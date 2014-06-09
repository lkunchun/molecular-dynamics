#ifndef PARSER_H
#define PARSER_H

#include <stdio.h>
#include <stdlib.h>

#define MAX_BUFFER 80  //do not change this number
#define MAX_DATA 12

enum INPUTS {
   PARSER_ERROR = -2,
   QUITTING = -1,
   RUN = 0,
   ATOM_INPUT,
   TEMPERATURE_COEFF,
   BOX_COEFF,
   TIME_COEFF,
   ANDERSEN_COEFF,
   VELOCITY_SCALE_COEFF,
   CONFIG_OUTFILE,
   COMMENT,
   SET_SEED,
   RESTART_OUTFILE,
   THERMO_OUTFILE,
   GROW_CRYSTAL
};


extern char* fortran_buffer; 

#endif
