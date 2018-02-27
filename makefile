NAME = simulation
TARGET = $(NAME)
OBJS = $(NAME).o mc.o
CC = icpc
#CC = g++
DEBUG = -g

# -funroll-loops
#CFLAGS = -Wall -c -std=c++11 -cxxlib=/u/shared/programs/x86_64/gcc/5.4.0/ -marc$
#LFLAGS = -Wall -cxxlib=/u/shared/programs/x86_64/gcc/5.4.0/ $(DEBUG)
CFLAGS = -Wall -c -std=c++11$
LFLAGS = -Wall $(DEBUG)
LIBS =
#LIBS = -I/u/cm/fsurace/armadillo/usr/include/ -L/u/cm/fsurace/armadillo/usr/li$

$(TARGET) : $(OBJS)
	$(CC) $(LFLAGS) $(LIBS) $(OBJS) -o $(TARGET)

$(NAME).o : $(NAME).cpp mc.h
	$(CC) $(CFLAGS) $(LIBS) $(NAME).cpp

mc.o : mc.h mc.cpp
	$(CC) $(CFLAGS) $(LIBS) mc.cpp

clean:
	\rm *.o $(TARGET)


# COMPILARE:
# module load mkl/11.2 intel/18.0
# icpc -Wall -std=c++11 -cxxlib=/u/shared/programs/x86_64/gcc/5.4.0/ -march=nat$
#
# ESEGUIRE:
# export LD_LIBRARY_PATH=/u/shared/programs/x86_64/gcc/5.4.0/lib64/:$LD_LIBRARY$
# ./testing2
# COMPILARE:
# module load mkl/11.2 intel/18.0
# icpc -Wall -std=c++11 -cxxlib=/u/shared/programs/x86_64/gcc/5.4.0/ -march=nat$
#
# ESEGUIRE:
# export LD_LIBRARY_PATH=/u/shared/programs/x86_64/gcc/5.4.0/lib64/:$LD_LIBRARY$
# ./testing2


