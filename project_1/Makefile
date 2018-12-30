CC = gcc
CFLAGS 	= -Wall -Wextra -pedantic -std=c11 -I./include
COFLAGS = -O3
CGFLAGS = $(shell pkg-config libgvc --cflags)
LDFLAGS = $(shell pkg-config libgvc --libs)	-lm
TFLAGS  = -I./tests
DFLAGS  = -g -pg

SRCDIR 		:= src/
CYTDIR 		:= cython_test/
TESTDIR 	:= tests/
OBJDIR		:= obj/
BINDIR 		:= bin/

SRCFILES 	:= $(shell find $(SRCDIR) \( -name "*.c" -not -name "print_tree.c" \))
TESTFILES 	:= $(shell find $(TESTDIR) \( -name "*.c"					$\
										-not -name "mem_test.c" \))
CYTFILES 	:= $(shell find $(CYTDIR) -name "*.c")
MEMFILES 	:= $(TESTDIR)mem_test.c
ALLSRC		:= $(SRCFILES) $(TESTFILES)

OBJSRC 		:= $(patsubst $(SRCDIR)%.c, $(OBJDIR)%.o, $(SRCFILES))
OBJTEST 	:= $(patsubst $(TESTDIR)%.c, $(OBJDIR)%.o, $(TESTFILES))
OBJCYT 		:= $(patsubst $(CYTDIR)%.c, $(OBJDIR)%.o, $(CYTDIR))
OBJMEM   	:= $(OBJDIR)mem_test.o
ALLOBJ		:= $(OBJSRC) $(OBJTEST) $(OBJVISUAL) $(OBJMEM)

TEST		:= $(BINDIR)test.out
CYTEST 		:= $(BINDIR)cytest.out
MEMTEST		:= $(BINDIR)mem_test.out
ALLBIN 		:= $(TEST) $(VISUAL) $(MEMTEST)


test: $(OBJSRC) $(OBJTEST)
	$(CC) $^ -o $(TEST)

cytest: $(OBJSRC) $(OBJCYT)
	$(CC) $^ -o $(VISUAL)

memtest: $(OBJSRC) $(OBJMEM)
	$(CC) $^ -o $(MEMTEST)


$(OBJDIR)%.o: $(CYTDIR)%.c
	$(CC) $(CFLAGS) $(DFLAGS) -c $< -o $@

$(OBJMEM): $(MEMFILES)
	$(CC) $(CFLAGS) $(DFLAGS) -c $< -o $@

$(OBJDIR)%.o: $(SRCDIR)%.c
	$(CC) $(CFLAGS) $(DFLAGS) -c $< -o $@

$(OBJDIR)%.o: $(TESTDIR)%.c
	$(CC) $(CFLAGS) $(DFLAGS) $(TFLAGS) -c $< -o $@


clean:
	rm -f $(ALLOBJ) $(ALLBIN) # data/*
