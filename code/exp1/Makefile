# This is a makefile.

PROG = openmP
OBJ = main.o OpenMPExamples.o

#flags.
CC = g++ -fopenmp -std=c++11
LFLAG = -g
CFLAG = -c
OFLAG = -o

all: $(PROG)

$(PROG) : $(OBJ)
	$(CC) $(OFLAG) $(PROG) $(OBJ) $(LFLAG)
%.o : %.cpp
	$(CC) -c $(CFLAG) $< -o $@

clean:
	rm -f *.o