/*
     
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

const int X = 13;
const int Y = 13;

void record(std::string message, bool includeTime)
{
    clock_t total = (double)(std::clock() - mainClock) / CLOCKS_PER_SEC;
    if (mainSimpi->getID() == 0)
    {
        std::cout << message;
        if (includeTime)
            std::cout << "completed in " << total << " seconds";
        std::cout << std::endl;
    }
}

void test_equality()
{
    std::string passMessage = "test_equality() ";

    Matrix A(X, Y);
    A.fillRandom(-50, 50);
    Matrix B = A;

    mainClock = std::clock();
    A == B;
    record(passMessage, true);
}

void test_matrix_multiplication()
{
    std::string passMessage = "test_matrix_multiplication() ";

    Matrix A(X, Y);
    A.fillRandom(-50, 50);
    Matrix B(X, Y);
    B.fillRandom(-50, 50);
    
    mainClock = std::clock();
    A * B;
    record(passMessage, true);
}

void test_inverse()
{
    std::string passMessage = "test_inverse() ";

    Matrix A(X, Y);
    A.fillRandom(-50, 50);
    Matrix inv(X, Y);

    mainClock = std::clock();
    A.inverse(inv);
    record(passMessage, true);
}

/***
 * DO NOT RUN AT ALL
*/
void test_determinant()
{
    std::string passMessage = "test_determinant() ";

    Matrix A(X, Y);
    A.fillRandom(-50, 50);
    
    mainClock = std::clock();
    A.determinant();
    record(passMessage, true);
}

void test_adjoint()
{
    std::string passMessage = "test_adjoint() ";
   
    Matrix A(X, Y);
    A.fillRandom(-50, 50);
    Matrix adj(X, Y);

    mainClock = std::clock();
    A.adjoint(adj);
    record(passMessage, true);
}

/***
 * DO NOT RUN WITH MULTIPLE PROCESSES
*/
void test_solve_system_with_jacobi() 
{
    std::string passMessage = "test_solve_system_with_jacobi() ";
    
    Matrix A(X, Y);
    A.fillRandom(-50, 50);
    for (int i = 0; i < X; i++) // Assuming X and Y have same value
        A.get(i, i) *= 10; // To pass diagonal dominance test

    Matrix B(X, 1);
    B.fillRandom(-50, 50);

    Matrix x(X, 1);

    mainClock = std::clock();
    A.solveSystem(x, B);
    record(passMessage, true);
}

void test_solve_system_with_failsafe()
{
    std::string passMessage = "test_solve_system_with_failsafe() ";

    Matrix A(X, Y);
    A.fillRandom(-50, 50);
    for (int i = 0; i < X; i++) // Assuming X and Y have same value
        A.get(i, i) /= 10; // To fail diagonal dominance test

    Matrix B(X, 1);
    B.fillRandom(-50, 50);

    Matrix x(X, 1);

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
    test_equality();
    test_matrix_multiplication();
    test_inverse();
    // test_determinant(); // DO NOT RUN UNTIL ALGORITHM TAKES ADVANTAGE OF MULTUPLE PROCESSES
    test_adjoint();
    // test_solve_system_with_jacobi(); // DO NOT RUN WITH MULTIPLE PROCESSES
    test_solve_system_with_failsafe();
    record("\n---Completed tests for " + std::to_string(mainSimpi->getProcessCount()) + " processes", false);
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

