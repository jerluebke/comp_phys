
# run with:
# mpirun -np 4 heat_parallel

MPIF90 = mpifort
FC = gfortran
CXX = g++

LDFLAGS =  -lstdc++

VTKSRC = output.cpp
VTKOBJ = $(VTKSRC:%.cpp=%.o)

all: heat_parallel

heat_parallel: output.o heat_parallel.o $(VTKOBJ)
	$(MPIF90) -o $@ output.o heat_parallel.o $(LDFLAGS)

$(VTKOBJ): %.o: %.cpp
	$(CXX) -c -O3 $<

heat_parallel.o: heat_parallel.F90
	$(MPIF90) -c -O3 $<

clean:
	rm -rf *.o *.mod heat_parallel
