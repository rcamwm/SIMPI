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
const int TESTS = 10;

void record(std::string message, bool includeTime, double averageTime)
{
    if (mainSimpi->getID() == 0)
    {
        if (includeTime)
            std::cout << std::fixed << std::setprecision(4) << averageTime << " seconds: ";
        std::cout << message << std::endl;;
    }
}

void test_equality()
{
    std::string passMessage = "test_equality()";
    double averageTime = 0;
    for (int i = 0; i < TESTS; i++)
    {
        Matrix A(ROWS, COLS);
        A.fillRandom(-50, 50);
        Matrix B = A;

        mainClock = std::clock();
        A == B;
        averageTime += (double)(std::clock() - mainClock);
    }
    averageTime /= (double)(CLOCKS_PER_SEC * TESTS);
    record(passMessage, true, averageTime);
}

void test_copy_constructor()
{
    std::string passMessage = "test_copy_constructor()";
    double averageTime = 0;
    for (int i = 0; i < TESTS; i++)
    {
        Matrix A(ROWS, COLS);
        A.fillRandom(-50, 50);

        mainClock = std::clock();
        Matrix B = A;
        averageTime += (double)(std::clock() - mainClock);
    }
    averageTime /= (double)(CLOCKS_PER_SEC * TESTS);
    record(passMessage, true, averageTime);
}

void test_matrix_multiplication()
{
    std::string passMessage = "test_matrix_multiplication()";
    double averageTime = 0;
    for (int i = 0; i < TESTS; i++)
    {
        Matrix A(ROWS, COLS);
        A.fillRandom(-50, 50);
        Matrix B(ROWS, COLS);
        B.fillRandom(-50, 50);
        
        mainClock = std::clock();
        A * B;
        averageTime += (double)(std::clock() - mainClock);
    }
    averageTime /= (double)(CLOCKS_PER_SEC * TESTS);
    record(passMessage, true, averageTime);
}

void test_matrix_scalar_multiplication()
{
    std::string passMessage = "test_matrix_scalar_multiplication()";
    double averageTime = 0;
    for (int i = 0; i < TESTS; i++)
    {
        Matrix A(ROWS, COLS);
        A.fillRandom(-50, 50);
        double lambda = A.getVal(0,0) * 1.23;
        
        mainClock = std::clock();
        A * lambda;
        averageTime += (double)(std::clock() - mainClock);
    }
    averageTime /= (double)(CLOCKS_PER_SEC * TESTS);
    record(passMessage, true, averageTime);
}

void test_matrix_addition()
{
    std::string passMessage = "test_matrix_addition()";
    double averageTime = 0;
    for (int i = 0; i < TESTS; i++)
    {
        Matrix A(ROWS, COLS);
        A.fillRandom(-50, 50);

        Matrix B(ROWS, COLS);
        A.fillRandom(-50, 50);
        
        
        mainClock = std::clock();
        A + B;
        averageTime += (double)(std::clock() - mainClock);
    }
    averageTime /= (double)(CLOCKS_PER_SEC * TESTS);
    record(passMessage, true, averageTime);
}

void test_matrix_subtraction()
{
    std::string passMessage = "test_matrix_subtraction()";
    double averageTime = 0;
    for (int i = 0; i < TESTS; i++)
    {
        Matrix A(ROWS, COLS);
        A.fillRandom(-50, 50);

        Matrix B(ROWS, COLS);
        A.fillRandom(-50, 50);
        
        
        mainClock = std::clock();
        A - B;
        averageTime += (double)(std::clock() - mainClock);
    }
    averageTime /= (double)(CLOCKS_PER_SEC * TESTS);
    record(passMessage, true, averageTime);
}

void test_transpose()
{
    std::string passMessage = "test_transpose()";
    double averageTime = 0;
    for (int i = 0; i < TESTS; i++)
    {
        Matrix A(ROWS, COLS);
        A.fillRandom(-50, 50);
        
        mainClock = std::clock();
        A.transpose();
        averageTime += (double)(std::clock() - mainClock);
    }
    averageTime /= (double)(CLOCKS_PER_SEC * TESTS);
    record(passMessage, true, averageTime);
}

void test_inverse()
{
    std::string passMessage = "test_inverse()";
    double averageTime = 0;
    for (int i = 0; i < TESTS; i++)
    {
        Matrix A(ROWS, COLS);
        A.fillRandom(-50, 50);

        mainClock = std::clock();
        A.inverse();
        averageTime += (double)(std::clock() - mainClock);
    }
    averageTime /= (double)(CLOCKS_PER_SEC * TESTS);
    record(passMessage, true, averageTime);
}

/***
 * DO NOT RUN AT ALL
*/
void test_determinant()
{
    std::string passMessage = "test_determinant()";
    double averageTime = 0;
    for (int i = 0; i < TESTS; i++)
    {
         Matrix A(ROWS, COLS);
        A.fillRandom(-50, 50);
        
        mainClock = std::clock();
        A.determinant();
        averageTime += (double)(std::clock() - mainClock);
    }
    averageTime /= (double)(CLOCKS_PER_SEC * TESTS);
    record(passMessage, true, averageTime);
}

void test_adjoint()
{
    std::string passMessage = "test_adjoint()";
    double averageTime = 0;
    for (int i = 0; i < TESTS; i++)
    {
        Matrix A(ROWS, COLS);
        A.fillRandom(-50, 50);

        mainClock = std::clock();
        A.adjoint();
        averageTime += (double)(std::clock() - mainClock);
    }
    averageTime /= (double)(CLOCKS_PER_SEC * TESTS);
    record(passMessage, true, averageTime);
}

/***
 * DO NOT RUN WITH MULTIPLE PROCESSES
*/
void test_solve_system_with_jacobi() 
{
    std::string passMessage = "test_solve_system_with_jacobi()";
    double averageTime = 0;
    for (int i = 0; i < TESTS; i++)
    {
        Matrix A(ROWS, COLS);
        A.fillRandom(-50, 50);
        for (int i = 0; i < ROWS; i++) // Assuming ROWS and COLS have same value
            A.getRef(i, i) *= 10; // To pass diagonal dominance test

        Matrix B(ROWS, 1);
        B.fillRandom(-50, 50);

        Matrix x(ROWS, 1);

        mainClock = std::clock();
        A.solveSystem(x, B);
        averageTime += (double)(std::clock() - mainClock);
    }
    averageTime /= (double)(CLOCKS_PER_SEC * TESTS);
    record(passMessage, true, averageTime);
}

void test_solve_system_with_failsafe()
{
    std::string passMessage = "test_solve_system_with_failsafe()";
    double averageTime = 0;
    for (int i = 0; i < TESTS; i++)
    {
        Matrix A(ROWS, COLS);
        A.fillRandom(-50, 50);
        for (int i = 0; i < ROWS; i++) // Assuming ROWS and COLS have same value
            A.getRef(i, i) /= 10; // To fail diagonal dominance test

        Matrix B(ROWS, 1);
        B.fillRandom(-50, 50);

        Matrix x(ROWS, 1);

        mainClock = std::clock();
        A.solveSystem(x, B);
        averageTime += (double)(std::clock() - mainClock);
    }
    averageTime /= (double)(CLOCKS_PER_SEC * TESTS);
    record(passMessage, true, averageTime);
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
    record("\n---Starting tests for " + std::to_string(mainSimpi->getProcessCount()) + " processes", false, 0);
    record("---Matrix rows = " + std::to_string(ROWS) + ", Matrix cols = " + std::to_string(COLS), false, 0);
    record("---Average time for " + std::to_string(TESTS) + " tests:", false, 0);
    test_equality();
    test_copy_constructor();
    test_matrix_multiplication();
    test_matrix_scalar_multiplication();
    test_matrix_addition();
    test_matrix_subtraction();
    test_transpose();
    test_inverse();
    // test_determinant(); // DO NOT RUN UNTIL ALGORITHM TAKES ADVANTAGE OF MULTUPLE PROCESSES
    test_adjoint();
    // test_solve_system_with_jacobi(); // DO NOT RUN WITH MULTIPLE PROCESSES
    test_solve_system_with_failsafe();
    record("---Completed tests for " + std::to_string(mainSimpi->getProcessCount()) + " processes", false, 0);
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

