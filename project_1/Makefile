CC = gcc
CFLAGS = -Wall -Wextra -pedantic -std=c99 -I./include
TFLAGS = -I./tests
DFLAGS = -g -pg

SRCDIR 		:= src/
TESTDIR 	:= tests/
OBJDIR 		:= obj/
BINDIR 		:= bin/
SRCFILES 	:= $(shell find $(SRCDIR) -name "*.c")
OBJSRC		:= $(patsubst $(SRCDIR)%.c, $(OBJDIR)%.o, $(SRCFILES))
TESTFILES 	:= $(shell find $(TESTDIR) -name "*c")
OBJTEST		:= $(patsubst $(TESTDIR)%.c, $(OBJDIR)%.o, $(TESTFILES))
ALLFILES	:= $(SRCFILES) $(TESTFILES)
TEST		:= $(BINDIR)test.out

test: $(OBJSRC) $(OBJTEST)
	$(CC) $^ -o $(TEST)

$(OBJDIR)%.o: $(SRCDIR)%.c
	$(CC) $(CFLAGS) $(DFLAGS) -c $< -o $@

$(OBJDIR)%.o: $(TESTDIR)%.c
	$(CC) $(CFLAGS) $(TFLAGS) $(DFLAGS) -c $< -o $@

clean:
	rm -f $(OBJSRC) $(OBJTEST) $(TEST)
