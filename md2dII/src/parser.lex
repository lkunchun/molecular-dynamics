%option noyywrap

/* 
 * Lexer for the md2d program
 * text scanner
 *
 */

%{

//#undef YY_MAIN 
//#include "y.tab.h"
#include "parser.h"

char* fortran_buffer; 
//int parser_counter = 0;
static int counter=0;

#undef YY_DECL
#define YY_DECL int lexscan(str,data,flag) char* str; double* data; int* flag; 

#undef YY_INPUT
#define YY_INPUT(buffer, result, max) (result = myinput(buffer, max))

%}

%s thermo restart seed configout atom box time temperature andersen velscale
%%

#.*(\n)? {return COMMENT;}

[qQ](uit)? {fprintf(stderr,"quitting program . . .\n");return QUITTING;}

@?r(un)? |
<<EOF>> {return RUN;}


(@*)?([ ]*)?s(eed)? {printf("@%s\n", fortran_buffer);BEGIN(seed);}
(@*)?([ ]*)?a(tom)?(s)?  {printf("@%s\n", fortran_buffer);BEGIN(atom);}
i(atom)?  {BEGIN(atom);}
(@*)?([ ]*)?b(ox)? {printf("@%s\n", fortran_buffer);BEGIN(box);}
(@*)?([ ]*)?config(out(file)?)? {printf("@%s\n", fortran_buffer);BEGIN(configout);}
(@*)?([ ]*)?thermo(out(file)?)? {printf("@%s\n", fortran_buffer);BEGIN(thermo);}
(@*)?([ ]*)?restart(out(file)?)? {printf("@%s\n", fortran_buffer);BEGIN(restart);}
(@*)?([ ]*)?time {printf("@%s\n", fortran_buffer);BEGIN(time);}
(@*)?([ ]*)?temperature {printf("@%s\n", fortran_buffer);BEGIN(temperature);}
(@*)?([ ]*)?andersen {printf("@%s\n", fortran_buffer);BEGIN(andersen);}
(@*)?([ ]*)?velscale {printf("@%s\n", fortran_buffer);BEGIN(velscale);}

<*>[\t ] /* do nothing*/ ;

<*>-?([0-9]*)(\.([0-9]+)?)? { data[*flag]=atof(yytext); (*flag)++;}

<configout,thermo,restart>[a-zA-Z][a-zA-Z_0-9]* {

         for(counter=0;counter<MAX_BUFFER;counter++){
            if(counter<yyleng)
               str[counter]=yytext[counter];
            else
               str[counter]= ' ';
         }
         (*flag)++;
}

<atom>\n {if (*flag>0) return ATOM_INPUT; else return PARSER_ERROR;}

<seed>\n {if (*flag>0) return SET_SEED; else return PARSER_ERROR;}

<box>\n {if (*flag>0)  return BOX_COEFF; else return PARSER_ERROR;}

<time>\n {if (*flag>0)  return TIME_COEFF; else return PARSER_ERROR;}

<temperature>\n {if (*flag>0) return TEMPERATURE_COEFF; else return PARSER_ERROR;}

<andersen>\n {if (*flag>0) return ANDERSEN_COEFF; else return PARSER_ERROR;}

<velscale>\n {if (*flag>0) return VELOCITY_SCALE_COEFF; else return PARSER_ERROR;}

<configout>\n {if (*flag>0) return CONFIG_OUTFILE; else return PARSER_ERROR;}

<thermo>\n {if (*flag>0) return THERMO_OUTFILE; else return PARSER_ERROR;}

<restart>\n {if (*flag>0) return RESTART_OUTFILE; else return PARSER_ERROR;}

<INITIAL>\n {BEGIN(INITIAL); return COMMENT; //assume it is a comment
                  } 

. {fprintf(stderr,"unknow command: %s\n", yytext); return PARSER_ERROR;}

%%


int myinput(char* buffer, int max) {

   int i;
   if(MAX_BUFFER>max) {
      fprintf(stderr,"error has occured\n");
      exit(0);
   } else {

      for(i = 0; i<MAX_BUFFER;i++) {
         buffer[i] = fortran_buffer[i];
      }

      buffer[MAX_BUFFER - 1]='\n';
      return MAX_BUFFER;

   }

}

