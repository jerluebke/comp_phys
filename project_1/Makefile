CC = gcc
CFLAGS 	= -Wall -Wextra -pedantic -std=c11 -I./include
COFLAGS = -O3 -DNDEBUG
CGFLAGS = $(shell pkg-config libgvc --cflags)
LDFLAGS = $(shell pkg-config libgvc --libs)
TFLAGS  = -I./tests
DFLAGS  = -Og -g -pg

SRCDIR 		:= src/
TESTDIR 	:= tests/
OBJDIR		:= obj/
BINDIR 		:= bin/

SRCFILES 	:= $(shell find $(SRCDIR) \( -name "*.c" -not -name "cvisualise.c" \))
TESTFILES 	:= $(shell find $(TESTDIR) \( -name "*.c" -not -name "mem_test.c" -not -name "timeit.c" \))
TIMEITFILES := $(TESTDIR)timeit.c
ALLSRC		:= $(SRCFILES) $(TESTFILES)

OBJSRC 		:= $(patsubst $(SRCDIR)%.c, $(OBJDIR)%.o, $(SRCFILES))
OBJTEST 	:= $(patsubst $(TESTDIR)%.c, $(OBJDIR)%.o, $(TESTFILES))
OBJTIME		:= $(OBJDIR)timeit.o
ALLOBJ		:= $(OBJSRC) $(OBJTEST) $(OBJMEM) $(OBJTIME)

TEST		:= $(BINDIR)test.out
TIMEIT 		:= $(BINDIR)timeit.out
ALLBIN 		:= $(TEST) $(MEMTEST) $(TIMEIT)


test: $(OBJSRC) $(OBJTEST)
	$(CC) $^ -lm -o $(TEST)

timeit: $(OBJSRC) $(OBJTIME)
	$(CC) $^ -lm -o $(TIMEIT)


$(OBJTIME): $(TIMEITFILES)
	$(CC) $(CFLAGS) $(COFLAGS) -c $< -o $@

$(OBJDIR)%.o: $(SRCDIR)%.c
	$(CC) $(CFLAGS) $(DFLAGS) -c $< -o $@

$(OBJDIR)%.o: $(TESTDIR)%.c
	$(CC) $(CFLAGS) $(DFLAGS) $(TFLAGS) -c $< -o $@


clean:
	rm -f $(ALLOBJ) $(ALLBIN) # data/*
