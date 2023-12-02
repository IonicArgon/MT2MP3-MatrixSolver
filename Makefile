
CC=gcc
CFLAGS=-lm -Ofast
DFLAGS=-Wall -Wextra -W -g -O0 -lm


FLAGSGaussSeidel= -DGAUSS_SEIDEL -DITERATIONS=100 -DPRECONDITIONING
FLAGSJacobi= -DJACOBI -DITERATIONS=100 -DPRECONDITIONING
FLAGSPrint= -DPRINT=2

# change the flag as you wish
# FLAGSGaussSeidel is the default
# FLAGSJacobi is the Jacobi method
# FLAGSPrint is to print the solution, 1 for raw print, 2 for pretty print
# if FLAGSPrint is not defined, the solution will be saved in solution.txt
# DITERATIONS is the number of iterations
# DPRECONDITIONING is to use preconditioning
# DUSER_INPUT is to use user input for iterations, method and preconditioning instead of the default values
FLAGS= $(CFLAGS) $(FLAGSGaussSeidel) $(FLAGSPrint) #-DUSER_INPUT

SDIR=src
IDIR=include
ODIR=obj
BDIR=bin

_DEPS = functions.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main.o functions.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: $(SDIR)/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) -I$(IDIR) $(FLAGS)

$(BDIR)/main: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) -I$(IDIR) $(FLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o $(BDIR)/main
	if [ -d "./r000hs/" ]; then rm -rf ./r000hs/; fi
	if [ -f solution.txt ]; then rm solution.txt; fi

