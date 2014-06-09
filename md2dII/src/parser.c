#include "parser.h"

int lexscan(char*, double*, int*); 

copyparse_(inst, str, data, flag) 
char* str;
double* data;
int* inst;
int* flag;
{
   int i;
   fortran_buffer=str;

   *flag = 0;
   (*inst) = lexscan(str,data,flag);

}
