CC = gcc
CFLAGS 	= -Wall -Wextra -pedantic -std=c11 -I./include
COFLAGS = -O3
CGFLAGS = $(shell pkg-config libgvc --cflags)
LDFLAGS = $(shell pkg-config libgvc --libs)	# -lm
TFLAGS  = -I./tests
DFLAGS  = -g -pg

SRCDIR 		:= src/
VISUALDIR 	:= visual/
TESTDIR 	:= tests/
OBJDIR		:= obj/
BINDIR 		:= bin/

SRCFILES 	:= $(shell find $(SRCDIR) -name "*.c")
TESTFILES 	:= $(shell find $(TESTDIR) \( -name "*.c"					$\
											-not -name "*printing.c"	$\
											-not -name "mem_test.c" \))
VFILES 		:= $(shell find $(VISUALDIR) -name "*.c")
MEMFILES 	:= $(TESTDIR)mem_test.c
ALLSRC		:= $(SRCFILES) $(TESTFILES)

OBJSRC 		:= $(patsubst $(SRCDIR)%.c, $(OBJDIR)%.o, $(SRCFILES))
OBJTEST 	:= $(patsubst $(TESTDIR)%.c, $(OBJDIR)%.o, $(TESTFILES))
OBJVISUAL 	:= $(patsubst $(VISUALDIR)%.c, $(OBJDIR)%.o, $(VFILES))
OBJMEM   	:= $(OBJDIR)mem_test.o
ALLOBJ		:= $(OBJSRC) $(OBJTEST) $(OBJVISUAL) $(OBJMEM)

TEST		:= $(BINDIR)test.out
VISUAL 		:= $(BINDIR)visual.out
MEMTEST		:= $(BINDIR)mem_test.out
ALLBIN 		:= $(TEST) $(VISUAL) $(MEMTEST)


test: $(OBJSRC) $(OBJTEST)
	$(CC) $^ -o $(TEST)

visual: $(OBJSRC) $(OBJVISUAL)
	$(CC) $^ $(LDFLAGS) -o $(VISUAL)

memtest: $(OBJSRC) $(OBJMEM)
	$(CC) $^ -o $(MEMTEST)


$(OBJDIR)%.o: $(VISUALDIR)%.c
	$(CC) $(CFLAGS) $(CGFLAGS) $(DFLAGS) -c $< -o $@

$(OBJMEM): $(MEMFILES)
	$(CC) $(CFLAGS) $(DFLAGS) -c $< -o $@

$(OBJDIR)%.o: $(SRCDIR)%.c
	$(CC) $(CFLAGS) $(DFLAGS) -c $< -o $@

$(OBJDIR)%.o: $(TESTDIR)%.c
	$(CC) $(CFLAGS) $(DFLAGS) $(TFLAGS) -c $< -o $@


clean:
	rm -f $(ALLOBJ) $(ALLBIN) # data/*
