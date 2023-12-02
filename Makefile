CC=gcc
CFLAGS=-lm -Ofast
DFLAGS=-Wall -Wextra -W -g -O0 -lm

FLAGSGaussSeidel= -DGAUSS_SEIDEL -DITERATIONS=100 -DPRECONDITIONING
FLAGSJacobi= -DJACOBI -DITERATIONS=100 -DPRECONDITIONING

SDIR=src
IDIR=include
ODIR=obj
BDIR=bin

_DEPS = functions.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main.o functions.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: $(SDIR)/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) -I$(IDIR) $(FLAGSGaussSeidel)

$(BDIR)/main: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) -I$(IDIR) $(FLAGSGaussSeidel)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o $(BDIR)/main
	rm solution.txt
	rm -rf ./r000hs/