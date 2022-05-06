CC = g++ 
CFLAGS = -Wall -g -std=c++11
LINKFLAGS = -lrt
# if mac use this LINFLAGS
# LINKFLAGS =

all : mpi user
clean : 
	rm -f user mpi /dev/shm/simpi_shared_mem
Simpi : Simpi.cpp Simpi.h
	$(CC) $(CFLAGS) -c Simpi.cpp -o Simpi $(LINKFLAGS)
mpi : mpi.cpp user  Simpi Simpi.h
	$(CC) $(CFLAGS) mpi.cpp -o mpi Simpi $(LINKFLAGS)
user : user.cpp Simpi Simpi.h Matrix Matrix.h
	$(CC) $(CFLAGS) user.cpp -o user Matrix Simpi $(LINKFLAGS)
Matrix : Matrix.cpp Matrix.h Simpi Simpi.h
	$(CC) $(CFLAGS) -c Matrix.cpp -o Matrix $(LINKFLAGS)