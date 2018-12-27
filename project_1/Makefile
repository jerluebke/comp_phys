CC = gcc
CFLAGS = -Wall -Wextra -pedantic -std=c11 -I./include
LFLAGS = -lm
TFLAGS = -I./tests
DFLAGS = -g -pg

SRCDIR 		:= src/
VISUALDIR 	:= visual/
TESTDIR 	:= tests/
OBJDIR		:= obj/
BINDIR 		:= bin/

SRCFILES 	:= $(shell find $(SRCDIR) -name "*.c")
TESTFILES 	:= $(shell find $(TESTDIR) \( -name "*.c"					$\
											-not -name "*printing.c"	$\
											-not -name "mem_test.c" \))
VFILES		:= $(VISUALDIR)print_tree.c
MEMFILES 	:= $(TESTDIR)mem_test.c
ALLSRC		:= $(SRCFILES) $(TESTFILES)

OBJSRC 		:= $(patsubst $(SRCDIR)%.c, $(OBJDIR)%.o, $(SRCFILES))
OBJTEST 	:= $(patsubst $(TESTDIR)%.c, $(OBJDIR)%.o, $(TESTFILES))
OBJVISUAL	:= $(OBJDIR)print_tree.o
OBJMEM   	:= $(OBJDIR)mem_test.o
ALLOBJ		:= $(OBJSRC) $(OBJTEST)

TEST		:= $(BINDIR)test.out
VISUAL 		:= $(BINDIR)visual.out
MEMTEST		:= $(BINDIR)mem_test.out


test: $(ALLOBJ)
	$(CC) $^ -o $(TEST)

visual: $(OBJSRC) $(OBJVISUAL)
	$(CC) $^ $(LFLAGS) -o $(VISUAL)

memtest: $(OBJSRC) $(OBJMEM)
	$(CC) $^ -o $(MEMTEST)


$(OBJVISUAL): $(VFILES)
	$(CC) $(CFLAGS) $(DFLAGS) -c $< -o $@

$(OBJMEM): $(MEMFILES)
	$(CC) $(CFLAGS) $(DFLAGS) -c $< -o $@

$(OBJDIR)%.o: $(SRCDIR)%.c
	$(CC) $(CFLAGS) $(DFLAGS) -c $< -o $@

$(OBJDIR)%.o: $(TESTDIR)%.c
	$(CC) $(CFLAGS) $(DFLAGS) $(TFLAGS) -c $< -o $@


clean:
	rm -f $(ALLOBJ) $(TEST)
