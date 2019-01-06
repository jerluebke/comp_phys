CC = gcc
CFLAGS 	= -Wall -Wextra -pedantic -std=c11 -I./include
COFLAGS = -O3
CGFLAGS = $(shell pkg-config libgvc --cflags)
LDFLAGS = $(shell pkg-config libgvc --libs) -lm
TFLAGS  = -I./tests
DFLAGS  = -g -pg

SRCDIR 		:= src/
TESTDIR 	:= tests/
OBJDIR		:= obj/
BINDIR 		:= bin/

SRCFILES 	:= $(shell find $(SRCDIR) \( -name "*.c" -not -name "print_tree.c" \))
TESTFILES 	:= $(shell find $(TESTDIR) \( -name "*.c" -not -name "mem_test.c" \))
MEMFILES 	:= $(TESTDIR)mem_test.c
ALLSRC		:= $(SRCFILES) $(TESTFILES)

OBJSRC 		:= $(patsubst $(SRCDIR)%.c, $(OBJDIR)%.o, $(SRCFILES))
OBJTEST 	:= $(patsubst $(TESTDIR)%.c, $(OBJDIR)%.o, $(TESTFILES))
OBJMEM   	:= $(OBJDIR)mem_test.o
ALLOBJ		:= $(OBJSRC) $(OBJTEST) $(OBJVISUAL) $(OBJMEM)

TEST		:= $(BINDIR)test.out
MEMTEST		:= $(BINDIR)mem_test.out
ALLBIN 		:= $(TEST) $(MEMTEST)


test: $(OBJSRC) $(OBJTEST)
	$(CC) $^ $(LDFLAGS) -o $(TEST)

memtest: $(OBJSRC) $(OBJMEM)
	$(CC) $^ -lm -o $(MEMTEST)


$(OBJMEM): $(MEMFILES)
	$(CC) $(CFLAGS) $(DFLAGS) -c $< -o $@

$(OBJDIR)%.o: $(SRCDIR)%.c
	$(CC) $(CFLAGS) $(CGFLAGS) $(DFLAGS) -c $< -o $@

$(OBJDIR)%.o: $(TESTDIR)%.c
	$(CC) $(CFLAGS) $(DFLAGS) $(TFLAGS) -c $< -o $@


clean:
	rm -f $(ALLOBJ) $(ALLBIN) # data/*
