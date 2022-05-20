#include <signal.h>
#include <string.h>

#include "Simpi.h"
#include "Matrix.h"
#include "Vector.h"
#define MATRIX_DIMENSION_X 4
#define MATRIX_DIMENSION_Y 4
using namespace SimpiNS;

int processID;
Simpi *mainSimpi;

void testEquality();
void testMultiplication();
void testDeterminant();
void testAdjoint();
void testInverse();
void testSolveSystemJacobi(); 
void testSolveSystemInverse();

/*
    userDefinedActivity() is to be filled in by the user of this program.
    All necessary initialization occurs in main() before this function is called.
    Can use Matrix and Vector objects here and perform operations with them.
    Use mainSimpi->synch() to synchronize processes.
*/
void userDefinedActivity() 
{
    testEquality();
    // testMultiplication();
    // testDeterminant();
    // testAdjoint();
    // testInverse();
    // testSolveSystemJacobi();
    // testSolveSystemInverse();
}

void testEquality()
{
    Matrix A(4, 4);
    Matrix B(4, 4);
    Matrix C(4, 4);
    for (int row = 0; row < A.getRows(); row++)
    {
        for (int col = 0; col < A.getCols(); col++)
        {
            A.get(row, col) = rand() % 10 + 1;
            B.get(row, col) = A.get(row, col);
            C.get(row, col) = A.get(row, col);
        }
    }
    C.get(3,3) += 0.001;
    bool AB = A == B;
    bool AC = A == C;
    bool AC2 = A != C;
    if (mainSimpi->getID() == 0) { std::cout << "\n1: " << AB << "\n0: " << AC << "\n1: " << AC2 << std::endl; }

    Matrix D(3, 4);
    Matrix E(4, 4);
    bool DE = D == E;
    if (mainSimpi->getID() == 0) { std::cout << "0: " << DE << std::endl; }
}

void testMultiplication()
{
    Matrix A(4, 4);
    for (int row = 0; row < A.getRows(); row++)
    {
        for (int col = 0; col < A.getCols(); col++)
        {
            A.get(row, col) = rand() % 10 + 1;
        }
    }
    Matrix B(4, 4);
    for (int row = 0; row < B.getRows(); row++)
    {
        for (int col = 0; col < B.getCols(); col++)
        {
            B.get(row, col) = rand() % 10 + 1;
        }
    }

    Matrix C = A * B;
    if (mainSimpi->getID() == 0) { std::cout << A << "times" << B << "equals" << C; }

    A *= B;
    if (mainSimpi->getID() == 0) { std::cout << "\nA *= B is now" << A; }
}

void testDeterminant()
{
    Matrix A(MATRIX_DIMENSION_X, MATRIX_DIMENSION_Y);
    for (int y = 0; y < MATRIX_DIMENSION_Y; y++)
    {
        for (int x = 0; x < MATRIX_DIMENSION_X; x++)
        {
            A.get(x, y) = rand() % 10 + 1;
        }
    }
    std::cout << A << "Determinant: " << A.determinant() << std::endl;
}

void testAdjoint()
{
    Matrix A(MATRIX_DIMENSION_X, MATRIX_DIMENSION_Y);
    Matrix adj(MATRIX_DIMENSION_X, MATRIX_DIMENSION_Y);
    mainSimpi->synch();

    for (int y = 0; y < MATRIX_DIMENSION_Y; y++)
    {
        for (int x = 0; x < MATRIX_DIMENSION_X; x++)
        {
            A.get(x, y) = rand() % 10 + 1;
        }
    }

    mainSimpi->synch();
    A.adjoint(&adj);
    if (mainSimpi->getID() == 0) { std::cout << "\nA = " << A << "\nadj(A) = " << adj << std::endl; }
}

void testSolveSystemJacobi()
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

void testSolveSystemInverse()
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
                A.get(x, y) /= 3;
    }
    mainSimpi->synch();

    for (int i = 0; i < MATRIX_DIMENSION_Y; i++)
        C.set(i, 0);
    
    mainSimpi->synch();
    A.solveSystem(&B, &C);
    if (mainSimpi->getID() == 0) { std::cout << A << "times\n" << C << "equals\n" << B; }
    
}

void testInverse()
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
    processID = atoi(argv[1]);
    mainSimpi = new Simpi(processID, atoi(argv[2])); // argv[2]: # of processes being used
    //Matrix::setEqualityPrecision(0.001f);
    Matrix::setSimpi(mainSimpi);
    Vector::setSimpi(mainSimpi);
    userDefinedActivity();
    delete mainSimpi;
}

