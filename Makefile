CC=gcc
CFLAGS=-lm -Ofast
DFLAGS=-Wall -Wextra -W -g -O0 -lm

SDIR=src
IDIR=include
ODIR=obj
BDIR=bin

_DEPS = functions.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main.o functions.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: $(SDIR)/%.c $(DEPS)
	$(CC) -c -o $@ $< $(DFLAGS) -I$(IDIR)

$(BDIR)/main: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) -I$(IDIR)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o $(BDIR)/main