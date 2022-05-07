#include <signal.h>
#include <string.h>

#include "Simpi.h"
#include "Matrix.h"
#define MATRIX_DIMENSION_X 10
#define MATRIX_DIMENSION_Y 10

int processID;

Simpi *initializeProcess(char* processIDString, char* processCountString);
void segfault_printer(int dummy);

void testInverse(Simpi *mainSimpi);
void testSolveSystemJacobi(Simpi *mainSimpi); // must be 3x3
void testSolveSystemInverse(Simpi *mainSimpi); // must be 3x3

int main(int argc, char* argv[])
{
    Simpi *mainSimpi = initializeProcess(argv[1], argv[2]);

    testInverse(mainSimpi);
    // testSolveSystemJacobi(mainSimpi);
    // testSolveSystemInverse(mainSimpi);
}

Simpi *initializeProcess(char* processIDString, char* processCountString)
{
    signal(SIGSEGV, segfault_printer);
    processID = atoi(processIDString); // Global

    Simpi *mainSimpi = new Simpi(processID, atoi(processCountString));
    Matrix::setSimpi(mainSimpi);
    return mainSimpi;
}

void testSolveSystemJacobi(Simpi *mainSimpi)
{
    Matrix A(MATRIX_DIMENSION_X, MATRIX_DIMENSION_Y);
    Matrix::Vector B(3);
    Matrix::Vector C(3);

    mainSimpi->synch();
    
    A.get(0, 0) = 10;
    A.get(0, 1) = 2;
    A.get(0, 2) = 4;

    A.get(1, 0) = 3;
    A.get(1, 1) = 12;
    A.get(1, 2) = 4;

    A.get(2, 0) = 6;
    A.get(2, 2) = 15;
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
    std::cout << A << "times\n" << C << "equals\n" << B;
}

void testSolveSystemInverse(Simpi *mainSimpi)
{
    Matrix A(MATRIX_DIMENSION_X, MATRIX_DIMENSION_Y);
    Matrix::Vector B(3);
    Matrix::Vector C(3);

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
    std::cout << A << "times\n" << C << "equals\n" << B;
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