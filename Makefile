CC = g++ 
CFLAGS = -Wall -g -std=c++11
LINKFLAGS = -lrt
# if mac use this LINFLAGS
# LINKFLAGS =

all : mpi user test_cases test_times
clean : 
	rm -f user mpi /dev/shm/simpi_shared_mem

Simpi : Simpi.cpp Simpi.h
	$(CC) $(CFLAGS) -c Simpi.cpp -o Simpi $(LINKFLAGS)
Matrix : Matrix.cpp Matrix.h Simpi Simpi.h
	$(CC) $(CFLAGS) -c Matrix.cpp -o Matrix $(LINKFLAGS)

mpi : mpi.cpp user  Simpi Simpi.h
	$(CC) $(CFLAGS) mpi.cpp -o mpi Simpi $(LINKFLAGS)
test_cases : test_cases.cpp Simpi Simpi.h Matrix Matrix.h
	$(CC) $(CFLAGS) test_cases.cpp -o test_cases Matrix Simpi $(LINKFLAGS)
test_times : test_times.cpp Simpi Simpi.h Matrix Matrix.h
	$(CC) $(CFLAGS) test_times.cpp -o test_times Matrix Simpi $(LINKFLAGS)
user : user.cpp Simpi Simpi.h Matrix Matrix.h
	$(CC) $(CFLAGS) user.cpp -o user Matrix Simpi $(LINKFLAGS)