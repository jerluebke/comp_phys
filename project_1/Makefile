CC = gcc
CFLAGS = -Wall -Wextra -pedantic -std=c11 -I./include
TFLAGS = -I./tests
DFLAGS = -g -pg

SRCDIR 		:= src/
TESTDIR 	:= tests/
OBJDIR		:= obj/
BINDIR 		:= bin/

SRCFILES 	:= $(shell find $(SRCDIR) -name "*.c")
TESTFILES 	:= $(shell find $(TESTDIR) \( -name "*.c"					$\
											-not -name "*printing.c"	$\
											-not -name "mem_test.c" \))
MEMFILES 	:= $(TESTDIR)mem_test.c
ALLSRC		:= $(SRCFILES) $(TESTFILES)

OBJSRC 		:= $(patsubst $(SRCDIR)%.c, $(OBJDIR)%.o, $(SRCFILES))
OBJTEST 	:= $(patsubst $(TESTDIR)%.c, $(OBJDIR)%.o, $(TESTFILES))
OBJMEM   	:= $(OBJDIR)mem_test.o
ALLOBJ		:= $(OBJSRC) $(OBJTEST)

TEST		:= $(BINDIR)test.out
MEMTEST		:= $(BINDIR)mem_test.out


test: $(ALLOBJ)
	$(CC) $^ -o $(TEST)

memtest: $(OBJSRC) $(OBJMEM)
	$(CC) $^ -o $(MEMTEST)

$(OBJMEM): $(MEMFILES)
	$(CC) $(CFLAGS) $(DFLAGS) -c $< -o $@

$(OBJDIR)%.o: $(SRCDIR)%.c
	$(CC) $(CFLAGS) $(DFLAGS) -c $< -o $@

$(OBJDIR)%.o: $(TESTDIR)%.c
	$(CC) $(CFLAGS) $(DFLAGS) $(TFLAGS) -c $< -o $@

clean:
	rm -f $(ALLOBJ) $(TEST)
