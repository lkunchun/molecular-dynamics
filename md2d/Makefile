FFLAGS = -O 

OBJECTS = advance.o dforce.o dqint.o energy.o md2d.o vave.o vint.o

F77 = gfortran $(FFLAGS)

md2d:   $(OBJECTS)
	$(F77) -o md2d $(OBJECTS)
	
.f.o:	
	$(F77) -c $<
