#include <signal.h>
#include <string.h>

#include "Simpi.h"
#include "Matrix.h"
#include "Vector.h"
#define MATRIX_DIMENSION_X 5
#define MATRIX_DIMENSION_Y 5
using namespace SimpiNS;

int processID;

void testInverse(Simpi *mainSimpi);
void testSolveSystemJacobi(Simpi *mainSimpi); // must be 3x3
void testSolveSystemInverse(Simpi *mainSimpi); // must be 3x3

/*
    userDefinedActivity() is to be filled in by the user of this program.
    All necessary initialization occurs in main() before this function is called.
    Can use Matrix and Vector objects here and perform operations with them.
    Use mainSimpi->sync() to synchronize processes.
*/
void userDefinedActivity(Simpi *mainSimpi) 
{
    testInverse(mainSimpi);
    testSolveSystemJacobi(mainSimpi);
    // testSolveSystemInverse(mainSimpi);
}

void testSolveSystemJacobi(Simpi *mainSimpi)
{
    Matrix A(MATRIX_DIMENSION_X, MATRIX_DIMENSION_Y);
    Vector B(MATRIX_DIMENSION_Y);
    Vector C(MATRIX_DIMENSION_Y);

    mainSimpi->synch();
    
    for (int y = 0; y < MATRIX_DIMENSION_Y; y++)
    {
        for (int x = 0; x < MATRIX_DIMENSION_X; x++)
        {
            A.get(x, y) = rand() % 10 + 1;
        }
        B.set(y, rand() % 10 + 1);
    }
    mainSimpi->synch();

    for (int y = 0; y < MATRIX_DIMENSION_Y; y++)
    {
        for (int x = 0; x < MATRIX_DIMENSION_X; x++)
            if (mainSimpi->getID() == 0 && x == y)
                A.get(x, y) *= 30;
    }
    mainSimpi->synch();

    for (int i = 0; i < MATRIX_DIMENSION_Y; i++)
        C.set(i, 0);
    mainSimpi->synch();

    A.solveSystem(&B, &C);
    if (mainSimpi->getID() == 0) { std::cout << A << "times\n" << C << "equals\n" << B; }
}

void testSolveSystemInverse(Simpi *mainSimpi)
{
    Matrix A(MATRIX_DIMENSION_X, MATRIX_DIMENSION_Y);
    Vector B(3);
    Vector C(3);

    mainSimpi->synch();
    
    A.get(0, 0) = 1;
    A.get(0, 1) = 2;
    A.get(0, 2) = 4;

    A.get(1, 0) = 3;
    A.get(1, 1) = 3;
    A.get(1, 2) = 4;

    A.get(2, 0) = 6;
    A.get(2, 2) = 4;
    A.get(2, 1) = 5;

    mainSimpi->synch();
    
    B.set(0, 5);
    B.set(1, 3);
    B.set(2, 9);

    mainSimpi->synch();

    for (int i = 0; i < 3; i++)
        C.set(i, 0);
    
    mainSimpi->synch();
    A.solveSystem(&B, &C);
    if (mainSimpi->getID() == 0) { std::cout << A << "times\n" << C << "equals\n" << B; }
    
}

void testInverse(Simpi *mainSimpi)
{
    Matrix A(MATRIX_DIMENSION_X, MATRIX_DIMENSION_Y);
    Matrix C(MATRIX_DIMENSION_X, MATRIX_DIMENSION_Y);
    mainSimpi->synch();

    for (int y = 0; y < MATRIX_DIMENSION_Y; y++)
    {
        for (int x = 0; x < MATRIX_DIMENSION_X; x++)
        {
            A.get(x, y) = rand() % 10 + 1;
        }
    }

    mainSimpi->synch();

    A.inverse(&C);
    std::cout << A;
    std::cout << C;
}

void segfault_printer(int dummy)
{
    char buf[20];
    sprintf(buf, "%d: segfaulted\n", processID);
    write(STDOUT_FILENO, buf, strlen(buf));
    exit(1);
}

int main(int argc, char* argv[])
{
    signal(SIGSEGV, segfault_printer);
    processID = atoi(argv[1]); // Global
    Simpi *mainSimpi = new Simpi(processID, atoi(argv[2])); // argv[2]: # of processes being used
    Matrix::setSimpi(mainSimpi);
    Vector::setSimpi(mainSimpi);
    userDefinedActivity(mainSimpi);
}

