# This is a makefile.

PROG = openmP
OBJ = main.o OpenMPExamples.o

#flags.
CC = g++ -fopenmp
LFLAG =
CFLAG = -c
OFLAG = -o

all: $(PROG)

$(PROG) : $(OBJ)
	$(CC) $(OFLAG) $(PROG) $(OBJ) $(LFLAG)
%.o : %.cpp
	$(CC) -c $(CFLAG) $< -o $@

clean:
	rm -f *.o