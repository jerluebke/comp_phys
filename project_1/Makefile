CC = gcc
CFLAGS = -Wall -Wextra -pedantic -std=c11 -I./include
TFLAGS = -I./tests
DFLAGS = -g -pg

SRCDIR 		:= src/
TESTDIR 	:= tests/
OBJDIR		:= obj/
BINDIR 		:= bin/

SRCFILES 	:= $(shell find $(SRCDIR) -name "*.c")
TESTFILES 	:= $(shell find $(TESTDIR) -name "*.c")
ALLSRC		:= $(SRCFILES) $(TESTFILES)

OBJSRC 		:= $(patsubst $(SRCDIR)%.c, $(OBJDIR)%.o, $(SRCFILES))
OBJTEST 	:= $(patsubst $(TESTDIR)%.c, $(OBJDIR)%.o, $(TESTFILES))
ALLOBJ		:= $(OBJSRC) $(OBJTEST)

TEST		:= $(BINDIR)test.out


test: $(ALLOBJ)
	$(CC) $^ -o $(TEST)

$(OBJDIR)%.o: $(SRCDIR)%.c
	$(CC) $(CFLAGS) $(DFLAGS) -c $< -o $@

$(OBJDIR)%.o: $(TESTDIR)%.c
	$(CC) $(CFLAGS) $(DFLAGS) $(TFLAGS) -c $< -o $@

clean:
	rm -f $(ALLOBJ) $(TEST)
