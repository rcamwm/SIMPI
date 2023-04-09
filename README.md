# SIMPI
Simple message passing interface in C++ designed for linux. 

## How to Use 
**SIMPI currently only runs on Linux.**  

To run this program, first clone the repository to a local linux machine and edit the user.cpp file, then:

Compile with:
```
make
```

Run with:
```
./mpi <program> <n>
```

where `<program>` is `user`, `test_cases`, or `test_times`, and `<n>` is the desired number of parallel processes to use.


## Overview
SIMPI is a new simple way of executing complex math operations such as matrix inversion or equation solving quickly with the help of parallelization. This is accomplished by having the mpi program calling the user program n-times. In the example of matrix multiplication with two 10x10 matrices, the mpi program could call the user program up to 10 times and have each instance of user calculating one row of the multiplication. This allows matrices to be multiplied roughly 10x faster. This is especially true if SIMPI is implemented on a "supercomputer" and has access to a large number of cores. 

### Main Files 

#### user.cpp
The file a user needs to change to run math calculations. 

#### mpi.cpp
The file that will run the user file for each process, additionally it is the file that is actually run by a user.

#### Simpi.cpp
The file that contains all the backbone of the program.  
No interaction with this file is necessary from a user. 

#### Matrix.cpp
This file houses all of the mathematical functions.  
No interaction with this file is necessary from a user. 

#### test_cases.cpp and test_times.cpp
These files can optionally be run through `./mpi`, exactly how the user.cpp file would be.
test_cases.cpp runs through every major Matrix mathematical function and confirms it's working correctly.
test_times.cpp measures how long it takes to perform Matrix operations with larger matrices. 
