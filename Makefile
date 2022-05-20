CC = g++ 
CFLAGS = -Wall -g -std=c++11
LINKFLAGS = -lrt
# if mac use this LINFLAGS
# LINKFLAGS =

all : mpi user tests
clean : 
	rm -f user mpi /dev/shm/simpi_shared_mem

Simpi : Simpi.cpp Simpi.h
	$(CC) $(CFLAGS) -c Simpi.cpp -o Simpi $(LINKFLAGS)
Vector : Vector.cpp Vector.h Simpi Simpi.h
	$(CC) $(CFLAGS) -c Vector.cpp -o Vector $(LINKFLAGS)
Matrix : Matrix.cpp Matrix.h Vector Vector.h Simpi Simpi.h
	$(CC) $(CFLAGS) -c Matrix.cpp -o Matrix $(LINKFLAGS)

mpi : mpi.cpp user  Simpi Simpi.h
	$(CC) $(CFLAGS) mpi.cpp -o mpi Simpi $(LINKFLAGS)
tests : tests.cpp Simpi Simpi.h Vector Vector.h Matrix Matrix.h
	$(CC) $(CFLAGS) tests.cpp -o tests Matrix Vector Simpi $(LINKFLAGS)
user : user.cpp Simpi Simpi.h Vector Vector.h Matrix Matrix.h
	$(CC) $(CFLAGS) user.cpp -o user Matrix Vector Simpi $(LINKFLAGS)