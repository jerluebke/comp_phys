SRC := simu.cc init.cc step.cc output.cc
OBJ := $(SRC:%.cc=%.o)

test:  simu
	./simu

$(OBJ): %.o: %.cc simu.h
	g++ -Wall -std=c++0x -O2 -c -pg $<

simu:  $(OBJ)
	g++ -o simu -pg $^

clean:
	rm -f simu $(OBJ) data/partData.* data/fieldData.*



