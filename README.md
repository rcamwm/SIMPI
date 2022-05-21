# SIMPI
Simple message passing interface in C++ designed for linux. 

## Overview
SIMPI is a new simple way of executing complex math operations such as matrix inversion or equation solving quickly with the help of parallelization. This is accomplished by having the mpi program calling the user program n-times. In the example of matrix multiplication with two 10x10 matrices, the mpi program could call the user program up to 10 times and have each instance of user calculating one row of the multiplication. This allows matrices to be multiplied roughly 10x faster. This is especially true if SIMPI is implemented on a "supercomputer" and has access to a large number of cores. 

## Main Files 
* user.cpp
* mpi.cpp
* Simpi.cpp
* Matrix.cpp
* tests.cpp

## The User File
The user.cpp file is the file in which a user needs to change to run math calculations. 

## The MPI File
The mpi.cpp file is the file that will run the user file for each process, additionly it is the file that is actully run by a user.

## The SIMPI File
The Simpi.cpp file is the file that contains all of the code that makes up the backbone of the program. 
No interaction with this file is necessary from a user. 

## The Matrix File
The Matrix.cpp file houses all of the mathmatical functions.
No interaction with this file is necessary from a user. 

## The Tests File
The tests.cpp file can optionally be run through MPI, exactly how the user.cpp file would be.
It runs through every major Matrix mathmatical function and confirms it's working correctly.

## Usage 
To run this program first clone this repository to a local linux machine and edit the user.cpp file, then:

Compile with:
../SIMPI$ make

Run with:
../SIMPI$ ./mpi user number_of_parallel_processes_to_use

Run test cases with:
../SIMPI$ ./mpi tests number_of_parallel_processes_to_use