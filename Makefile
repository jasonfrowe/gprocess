#Name of C compiler
CCMPL = gcc
#Name of Fortran compiler
F90 = ifort
F77 = ifort
#compiling object file flags
OPT1 = -O3 -fast -parallel -ipo
#OPT1 = -O0 -g -CB
#OPT1 = -O3
#OPT2 = -heap-arrays 0
OPT2 = 
FFLAGS = -c $(OPT1) $(OPT2)
#FFLAGS = -c -O0 -g -CB -warn $(OPT2)
#linking flags
LFLAGS = $(OPT1) $(OPT2)
#LFLAGS = -O0 -g -CB -warn $(OPT2)
#testing flags
#LFLAGS = -O0 -g -CB
#fitsio libraries
FITSIODIR = /usr/local/lib
#Pgplot plot libraries
PGPLOTDIR = /usr/local/lib
#X11 libraries
X11DIR = /usr/X11R6/lib
# libraries for linking PGPLOT
PLPLOTDIR = -I/usr/local/Cellar/plplot/5.9.11/lib/fortran/modules/plplot -I/usr/local/Cellar/plplot/5.9.11/include/plplot -L/usr/local/Cellar/plplot/5.9.11/lib
#LIBS = $(PLPLOTDIR) -lplplotf95d -lplplotf95cd
LIBS = -L$(PGPLOTDIR) -L$(X11DIR) -lX11 -lpgplot -lpng
# libraries for linking CFITSIO
LIBS2 = -L$(PGPLOTDIR) -L$(X11DIR) -L$(FITSIODIR) -lX11 -lpgplot -lcfitsio -lpng
#Directory where executable are placed
BIN = /Users/rowe/Documents/gprocess/
#utils source directory
UTILS = utils/

#Listing of programs to create.
all: gptest

gptestincl = precision.o getdata.o plotdata.o
gptest: gptest.f90 $(gptestincl)
	$(F90) $(LFLAGS) -o $(BIN)$@ gptest.f90 $(gptestincl) $(LIBS)

#building object libraries
precision.o: $(UTILS)precision.f90
	$(F90) $(FFLAGS) $(UTILS)precision.f90
getdata.o: $(UTILS)getdata.f90
	$(F90) $(FFLAGS) $(UTILS)getdata.f90
plotdata.o: $(UTILS)plotdata.f90
	$(F90) $(FFLAGS) $(UTILS)plotdata.f90

# Removing object files
.PHONY : clean
clean :
	rm *.o
	rm *.mod
