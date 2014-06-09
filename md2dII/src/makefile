HEADER = makefile md2d.h parser.h parser_f.h
OBJ = parser_f.o parser.o lex.yy.o print.o random.o neighbor.o move.o init.o thermostats.o dforce.o stringlib_unix.o 
XOBJ = viewbox.o
F77 = gfortran 
CC = gcc
OPT = -O 
CFLAG = $(OPT) -c
XLIB_FLAG = -L/usr/X11R6/lib/ -lstdc++ -lsupc++ -lX11 -lXext -lm 
FLTK_LIB = -Wl,-Bsymbolic-functions -Wl,-z,relro /usr/lib/i386-linux-gnu/libfltk_gl.a -lGL \
/usr/lib/i386-linux-gnu/libfltk.a -lXext -lXft -lfontconfig -lfontconfig \
-lXinerama -ldl -lm -lX11


run.exe:	$(OBJ) $(XOBJ) $(HEADER)
			$(F77) $(OPT) -o sepme seps.f $(OBJ)
			$(F77) $(OPT) -o genme lattice.f $(OBJ)
		   $(F77) $(OPT) -o runme md2d.f $(OBJ) $(XOBJ) $(FLTK_LIB) $(XLIB_FLAG)
			mv $(OBJ) $(XOBJ) obj

.cpp.o:	$(HEADER)
		$(CC) $(OPT) $(CFLAG) $< 

.c.o:	$(HEADER)
		$(CC) $(CFLAG) $<

.f.o:	$(HEADER)
		$(F77) $(CFLAG) $<

lex.yy.o:	$(HEADER) parser.lex
#			lex parser.lex
			$(CC) $(CFLAG) lex.yy.c

clean:
	rm -f *.o
