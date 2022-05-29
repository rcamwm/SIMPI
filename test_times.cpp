/*
    This file contains tests to check the speed of Matrix functions with different process counts.
    Each function loads Matrix objects with random values using Matrix::fillRandom().
    Except for the column vectors in the failsafe and jacobi functions
    all Matrix objects have the same dimensions: global variables ROWS and COLS.
    Keep in mind that many functions will not run if ROWS and COLS are not equal.
*/

#include <assert.h>
#include <signal.h>
#include <string.h>
#include <time.h>

#include "Simpi.h"
#include "Matrix.h"
using namespace SimpiNS;

int processID;
Simpi *mainSimpi;
clock_t mainClock;

const int ROWS = 10;
const int COLS = 10;

void record(std::string message, bool includeTime)
{
    double total = (double)(std::clock() - mainClock) / CLOCKS_PER_SEC;
    if (mainSimpi->getID() == 0)
    {
        if (includeTime)
            std::cout << std::fixed << std::setprecision(4) << total << " seconds: ";
        std::cout << message << std::endl;;
    }
}

void test_equality()
{
    std::string passMessage = "test_equality()";

    Matrix A(ROWS, COLS);
    A.fillRandom(-50, 50);
    Matrix B = A;

    mainClock = std::clock();
    A == B;
    record(passMessage, true);
}

void test_copy_constructor()
{
    std::string passMessage = "test_copy_constructor()";
    Matrix A(ROWS, COLS);
    A.fillRandom(-50, 50);

    mainClock = std::clock();
    Matrix B = A;
    record(passMessage, true); 
}

void test_matrix_multiplication()
{
    std::string passMessage = "test_matrix_multiplication()";

    Matrix A(ROWS, COLS);
    A.fillRandom(-50, 50);
    Matrix B(ROWS, COLS);
    B.fillRandom(-50, 50);
    
    mainClock = std::clock();
    A * B;
    record(passMessage, true);
}

void test_matrix_scalar_multiplication()
{
    std::string passMessage = "test_matrix_scalar_multiplication()";

    Matrix A(ROWS, COLS);
    A.fillRandom(-50, 50);
    double lambda = A.getVal(0,0) * 1.23;
    
    mainClock = std::clock();
    A * lambda;
    record(passMessage, true);
}

void test_inverse()
{
    std::string passMessage = "test_inverse()";

    Matrix A(ROWS, COLS);
    A.fillRandom(-50, 50);
    Matrix inv(ROWS, COLS);

    mainClock = std::clock();
    A.inverse(inv);
    record(passMessage, true);
}

/***
 * DO NOT RUN AT ALL
*/
void test_determinant()
{
    std::string passMessage = "test_determinant()";

    Matrix A(ROWS, COLS);
    A.fillRandom(-50, 50);
    
    mainClock = std::clock();
    A.determinant();
    record(passMessage, true);
}

void test_adjoint()
{
    std::string passMessage = "test_adjoint()";
   
    Matrix A(ROWS, COLS);
    A.fillRandom(-50, 50);
    Matrix adj(ROWS, COLS);

    mainClock = std::clock();
    A.adjoint(adj);
    record(passMessage, true);
}

/***
 * DO NOT RUN WITH MULTIPLE PROCESSES
*/
void test_solve_system_with_jacobi() 
{
    std::string passMessage = "test_solve_system_with_jacobi()";
    
    Matrix A(ROWS, COLS);
    A.fillRandom(-50, 50);
    for (int i = 0; i < ROWS; i++) // Assuming ROWS and COLS have same value
        A.getRef(i, i) *= 10; // To pass diagonal dominance test

    Matrix B(ROWS, 1);
    B.fillRandom(-50, 50);

    Matrix x(ROWS, 1);

    mainClock = std::clock();
    A.solveSystem(x, B);
    record(passMessage, true);
}

void test_solve_system_with_failsafe()
{
    std::string passMessage = "test_solve_system_with_failsafe()";

    Matrix A(ROWS, COLS);
    A.fillRandom(-50, 50);
    for (int i = 0; i < ROWS; i++) // Assuming ROWS and COLS have same value
        A.getRef(i, i) /= 10; // To fail diagonal dominance test

    Matrix B(ROWS, 1);
    B.fillRandom(-50, 50);

    Matrix x(ROWS, 1);

    mainClock = std::clock();
    A.solveSystem(x, B);
    record(passMessage, true);
}

void segfault_printer(int dummy)
{
    char buf[20];
    sprintf(buf, "%d: segfaulted\n", processID);
    write(STDOUT_FILENO, buf, strlen(buf));
    exit(1);
}

void runTests() 
{
    record("\n---Starting tests for " + std::to_string(mainSimpi->getProcessCount()) + " processes", false);
    record("---Matrix rows = " + std::to_string(ROWS) + ", Matrix cols = " + std::to_string(COLS), false);
    test_equality();
    test_copy_constructor();
    test_matrix_multiplication();
    test_matrix_scalar_multiplication();
    test_inverse();
    // test_determinant(); // DO NOT RUN UNTIL ALGORITHM TAKES ADVANTAGE OF MULTUPLE PROCESSES
    test_adjoint();
    // test_solve_system_with_jacobi(); // DO NOT RUN WITH MULTIPLE PROCESSES
    test_solve_system_with_failsafe();
    record("---Completed tests for " + std::to_string(mainSimpi->getProcessCount()) + " processes", false);
}

int main(int argc, char* argv[])
{
    signal(SIGSEGV, segfault_printer);
    processID = atoi(argv[1]);
    mainSimpi = new Simpi(processID, atoi(argv[2])); // argv[2]: # of processes being used
    Matrix::setSimpi(mainSimpi);
    mainClock = std::clock();
    runTests();
    delete mainSimpi; 
}

