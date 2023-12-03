
CC=gcc
CFLAGS=-lm -Ofast
DFLAGS=-Wall -Wextra -W -g -O0 -lm


#FLAGSGaussSeidel= -DGAUSS_SEIDEL -DMAX_ITER=10000 -DTHRESHOLD=1e- #-DPRECONDITIONING
FLAGSJacobi= -DJACOBI -DMAX_ITER=10000 -DTHRESHOLD=1e-7 #-DPRECONDITIONING
FLAGSPrint= -DPRINT=2

FLAGS= $(CFLAGS) $(FLAGSJacobi) # $(FLAGSPrint) #-DUSER_INPUT

SDIR=src
IDIR=include
ODIR=obj
BDIR=bin

_DEPS = functions.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main.o functions.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: $(SDIR)/%.c $(DEPS)
	$(CC) -c -o $@ $< -I$(IDIR) $(FLAGS)

$(BDIR)/main: $(OBJ)
	$(CC) -o $@ $^ -I$(IDIR) $(FLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o $(BDIR)/main
	if [ -d "./r000hs/" ]; then rm -rf ./r000hs/; fi
	if [ -f solution.txt ]; then rm solution.txt; fi
	if [ -f smvp_output.txt ]; then rm smvp_output.txt; fi

