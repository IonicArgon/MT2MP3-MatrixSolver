CC=gcc
CFLAGS=-lm -Ofast
DFLAGS=-Wall -Wextra -W -g -O0 -lm

FLAGSSOR= -DSOR -DMAX_ITER=10000 -DTHRESHOLD=1e-14 -DOMEGA=1.5  -DDIAGONAL_CHECK
FLAGSJacobi= -DJACOBI -DMAX_ITER=10000 -DTHRESHOLD=1e-14 -DDIAGONAL_CHECK 
FLAGSPrint= -DPRINT=1

FLAGS= $(CFLAGS) $(FLAGSSOR) # $(FLAGSPrint) #-DUSER_INPUT

SDIR=src
IDIR=include
ODIR=obj
BDIR=bin

_DEPS = functions.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main.o functions.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

all: check_dir $(BDIR)/main

$(ODIR)/%.o: $(SDIR)/%.c $(DEPS)
	$(CC) -c -o $@ $< -I$(IDIR) $(FLAGS)

$(BDIR)/main: $(OBJ)
	$(CC) -o $@ $^ -I$(IDIR) $(FLAGS)


.PHONY: clean check_dir


clean:
	rm -f $(ODIR)/*.o $(BDIR)/main
	if [ -d "./r000hs/" ]; then rm -rf ./r000hs/; fi
	if [ -f solution.txt ]; then rm solution.txt; fi
	if [ -f smvp_output.txt ]; then rm smvp_output.txt; fi

check_dir:
	if [ ! -d "./$(ODIR)/" ]; then mkdir $(ODIR); fi
	if [ ! -d "./$(BDIR)/" ]; then mkdir $(BDIR); fi

