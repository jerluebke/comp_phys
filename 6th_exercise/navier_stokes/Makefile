CFLAGS = -Wall -Wextra -I./include -L./lib/linux_libs/fftw-3.3.8
LDFLAGS = -lfftw3 -lm

ns : ./src/test.c ./src/navier_stokes.c
	gcc $^ -g -o $@ $(CFLAGS) $(LDFLAGS)
