/*
    This file contains test cases to determine if Matrix functions are continuing to give accurate results.
    Test cases use smaller Matrix sizes and have specific hard coded inputs that were decided with a random 
    number generator. These random inputs were then calculated outside of this program and those specific 
    results are compared to the output of each Matrix function.

    These test cases should be run whenever major changes or updates are made to either Matrix.h/cpp or Simpi.h/cpp.
*/

#include <assert.h>
#include <signal.h>
#include <string.h>

#include "Simpi.h"
#include "Matrix.h"
using namespace SimpiNS;

int processID;
Simpi *mainSimpi;

void assertAndPrint(std::string message, bool test)
{
    assert(test);
    if (processID == 0)
        std::cout << message << std::endl;    
}

void test_equality()
{
    std::string passMessage = "test_equality() passes test ";
    int testNo = 1;

    double x[] = {1, 2, 3,
                  4, 5, 6,
                  7, 8, 9};
    Matrix A(3, 3); 
    A.fill(x);
    assertAndPrint(passMessage + std::to_string(testNo++), A == A);

    Matrix B(3, 3); 
    B.fill(x);
    assertAndPrint(passMessage + std::to_string(testNo++), A == B);

    B.getRef(2, 2) += 0.001;
    assertAndPrint(passMessage + std::to_string(testNo++), !(A == B));

    double y[] = {1, 2, 3, 4, 5, 6, 7, 8, 9}; 
    Matrix C(1, 9);
    C.fill(y);
    assertAndPrint(passMessage + std::to_string(testNo++), !(A == C));
}

void test_inequality()
{
    std::string passMessage = "test_inequality_operator() passes test ";
    int testNo = 1;
    
    double x[] = {1, 2, 3,
                  4, 5, 6,
                  7, 8, 9};
    double y[] = {0.1, 2.0, 3.0,
                  4.0, 0.5, 6.0,
                  7.0, 8.0, 0.9};
    Matrix A(3, 3); 
    A.fill(x);
    Matrix B(3, 3); 
    B.fill(y);
    assertAndPrint(passMessage + std::to_string(testNo++), A != B);

    if (processID == 0)
        for (int i = 0; i < B.getCols(); i++) { B.getRef(i, i) *= 10; }
    mainSimpi->synch();
    assertAndPrint(passMessage + std::to_string(testNo++), !(A != B));

    Matrix C(2, 3);
    assertAndPrint(passMessage + std::to_string(testNo++), B != C);
}

void test_copy_constructor()
{
    std::string passMessage = "test_copy_constructor() passes test ";
    int testNo = 1;
    double x[] = {1, 2, 3,
                  4, 5, 6,
                  7, 8, 9};
    Matrix A(3, 3);
    A.fill(x);
    Matrix B = A;
    assertAndPrint(passMessage + std::to_string(testNo++), A == B);

    B.getRef(2,2) += 1;
    assertAndPrint(passMessage + std::to_string(testNo++), A != B);
}

void test_matrix_multiplication_same_size()
{
    std::string passMessage = "test_matrix_multiplication_operator_same_size() passes test ";
    int testNo = 1;
    double x[] = { 4,  7, 8, 6, 
                   4,  6, 7, 3, 
                  10,  2, 3, 8, 
                   1, 10, 4, 7};
    Matrix A(4, 4);
    A.fill(x);
    
    double y[] = {1, 7,  3,  7, 
                  2, 9,  8, 10, 
                  3, 1,  3,  4, 
                  8, 6, 10,  3};
    Matrix B(4, 4);
    B.fill(y);

    Matrix C = A * B;

    double z[] = {90, 135, 152, 148, 
                  61, 107, 111, 125, 
                  87, 139, 135, 126, 
                  89, 143, 165, 144};
    Matrix AB(4, 4);
    AB.fill(z);
    assertAndPrint(passMessage + std::to_string(testNo++), AB == C);

    A *= B;
    assertAndPrint(passMessage + std::to_string(testNo++), A == C);
}

void test_matrix_multiplication_more_rows()
{
    std::string passMessage = "test_matrix_multiplication_more_rows() passes test ";
    int testNo = 1;
    double x[] = { 4,  7,  8, 
                   6,  4,  6, 
                   7,  3,  10, 
                   2,  3,  8, 
                   1, 10,  4,};
    Matrix A(5, 3);
    A.fill(x);
    
    double y[] = {7, 1, 
                  7, 3, 
                  7, 2,};
    Matrix B(3, 2);
    B.fill(y);

    Matrix C = A * B;

    double z[] = {133, 41, 
                  112, 30, 
                  140, 36, 
                   91, 27, 
                  105, 39};
    Matrix AB(5, 2);
    AB.fill(z);
    assertAndPrint(passMessage + std::to_string(testNo++), AB == C);

    A *= B;
    assertAndPrint(passMessage + std::to_string(testNo++), A == C);

}

void test_matrix_multiplication_more_cols()
{
    std::string passMessage = "test_matrix_multiplication_more_cols() passes test ";
    int testNo = 1;
    double x[] = {4, 7, 8, 6, 
                  4, 6, 7, 3};
    Matrix A(2, 4);
    A.fill(x);
    
    double y[] = {10, 2, 3, 8, 1, 
                  10, 4, 7, 1, 7, 
                   3, 7, 2, 9, 8, 
                  10, 3, 1, 3, 4};
    Matrix B(4, 5);
    B.fill(y);

    Matrix C = A * B;

    double z[] = {194, 110, 83, 129, 141, 
                  151,  90, 71, 110, 114,};
    Matrix AB(2, 5);
    AB.fill(z);
    assertAndPrint(passMessage + std::to_string(testNo++), AB == C);

    A *= B;
    assertAndPrint(passMessage + std::to_string(testNo++), A == C);

}

void test_matrix_scalar_multiplication()
{
    std::string passMessage = "test_matrix_scalar_multiplication() passes test ";
    int testNo = 1;
    double x[] = { 1, 2, 3,
                   4, 5, 6,
                   7, 8, 9};
    Matrix A(3, 3);
    A.fill(x);
    
    double xTimes2Point5[] = { 2.5,  5.0,  7.5,
                              10.0, 12.5, 15.0,
                              17.5, 20.0, 22.5};
    Matrix A2Point5(3, 3);
    A2Point5.fill(xTimes2Point5);

    assertAndPrint(passMessage + std::to_string(testNo++), A2Point5 == A * 2.5);

    A.fill(x);
    assertAndPrint(passMessage + std::to_string(testNo++), A2Point5 == 2.5 * A);

    A.fill(x);
    A *= 2.5;
    assertAndPrint(passMessage + std::to_string(testNo++), A2Point5 == A);
}

void test_matrix_addition()
{
    std::string passMessage = "test_matrix_addition() passes test ";
    int testNo = 1;
    double x[] = { 1, 2, 3,
                   4, 5, 6,
                   7, 8, 9};
    Matrix A(3, 3);
    A.fill(x);

    double y[] = { -5.1,  3.2,  6.9,
                    2.4, -5.0, -5.6,
                    7.3, -2.8,  1.9};
    Matrix B(3, 3);
    B.fill(y);

    double z[] = {  -4.1, 5.2,  9.9,
                     6.4, 0.0,  0.4,
                    14.3, 5.2, 10.9};
    Matrix C(3, 3);
    C.fill(z);
    
    Matrix::setEqualityPrecision(0.00000001f);
    assertAndPrint(passMessage + std::to_string(testNo++), A + B == C);

    A.fill(x);
    A += B;
    assertAndPrint(passMessage + std::to_string(testNo++), A == C);

    Matrix::setEqualityPrecision(0.00f);
}

void test_matrix_subtraction()
{
    std::string passMessage = "test_matrix_subtraction() passes test ";
    int testNo = 1;

    double x[] = { -5.1,  3.2,  6.9,
                    2.4, -5.0, -5.6,
                    7.3, -2.8,  1.9};
    Matrix A(3, 3);
    A.fill(x);

    double y[] = {  4.2,  1.2, -9.3,
                    9.2, -2.7,  6.5,
                   10.3,  1.8,  2.1};
    Matrix B(3, 3);
    B.fill(y);

    double z[] = { -9.3,  2.0,  16.2,
                   -6.8, -2.3, -12.1,
                   -3.0, -4.6,  -0.2};
    Matrix C(3, 3);
    C.fill(z);
    Matrix::setEqualityPrecision(0.00000001f);
    assertAndPrint(passMessage + std::to_string(testNo++), A - B == C);

    A.fill(x);
    A -= B;
    assertAndPrint(passMessage + std::to_string(testNo++), A == C);

    Matrix::setEqualityPrecision(0.00f);
}

void test_matrix_transpose_same_size()
{
    std::string passMessage = "test_matrix_transpose_same_size() passes test ";
    int testNo = 1;
    double x[] = { 4,  7, 8, 6, 
                   4,  6, 7, 3, 
                  10,  2, 3, 8, 
                   1, 10, 4, 7};
    Matrix A(4, 4);
    A.fill(x);
    
    double y[] = {4, 4, 10,  1, 
                  7, 6,  2, 10, 
                  8, 7,  3,  4, 
                  6, 3,  8,  7};
    Matrix A_T(4, 4);
    A_T.fill(y);

    assertAndPrint(passMessage + std::to_string(testNo++), A.transpose() == A_T);
}

void test_matrix_transpose_more_rows()
{
    std::string passMessage = "test_matrix_transpose_more_rows() passes test ";
    int testNo = 1;
    double x[] = { 4,  7,  8, 
                   6,  4,  6, 
                   7,  3,  10, 
                   2,  3,  8, 
                   1, 10,  4,};
    Matrix A(5, 3);
    A.fill(x);
    
    double y[] = {4, 6,  7, 2,  1, 
                  7, 4,  3, 3, 10, 
                  8, 6, 10, 8,  4};
    Matrix A_T(3, 5);
    A_T.fill(y);

    assertAndPrint(passMessage + std::to_string(testNo++), A.transpose() == A_T);
}

void test_matrix_transpose_more_cols()
{
    std::string passMessage = "test_matrix_transpose_more_cols() passes test ";
    int testNo = 1;
    double x[] = {4, 7, 8, 6, 
                  4, 6, 7, 3};
    Matrix A(2, 4);
    A.fill(x);
    
    double y[] = {4, 4,
                  7, 6,
                  8, 7,
                  6, 3};
    Matrix A_T(4, 2);
    A_T.fill(y);

    assertAndPrint(passMessage + std::to_string(testNo++), A.transpose() == A_T);
}

void test_inverse()
{
    std::string passMessage = "test_inverse() passes test ";
    int testNo = 1;

    double x[] = {-5, 1,  3,  7, -2, 
                   8, 7, -8, -7,  8, 
                  -8, 7,  1, -6,  1, 
                   9, 1,  7,  5, -1, 
                  -6, 7,  10, 3,  3};
    Matrix A(5, 5);
    A.fill(x);
    Matrix inv = A.inverse();

    double y[] = {-0.0133363,  0.0174213,  0.0078267,  0.0726889, -0.0337270, 
                   0.1305650,  0.0546668,  0.1280080,  0.0901617, -0.0713501, 
                  -0.1119250, -0.0519893, -0.0115341,  0.0244585,  0.0760187, 
                   0.1357830,  0.0388418, -0.0603824, -0.0103498,  0.0036216, 
                  -0.0940235,  0.0417425, -0.1842020, -0.1361780,  0.1753460};
    
    Matrix B(5, 5);
    B.fill(y);
    
    Matrix::setEqualityPrecision(0.000001f);
    assertAndPrint(passMessage + std::to_string(testNo++), inv == B);

    Matrix::setEqualityPrecision(0.0f);
}

void test_determinant()
{
    std::string passMessage = "test_determinant() passes test ";
    int testNo = 1;

    double x[] = {4,  6,  3, 7,  2, 
                  7,  7,  8, 1,  9, 
                  8,  3,  1, 7,  8, 
                  6, 10, 10, 3, 10, 
                  4,  2,  4, 7,  3};
    Matrix A(5, 5);
    A.fill(x);
    assertAndPrint(passMessage + std::to_string(testNo++), (int)A.determinant() == -5399);

    double y[] = {-5.2, 1.0,  3.0,  7.2,
                   8.2, 7.0, -8.0, -7.3,
                  -8.5, 7.0,  1.0, -6.1,
                   9.9, 1.0,  7.0,  5.1};

    Matrix B(4, 4);
    B.fill(y);
    assertAndPrint(passMessage + std::to_string(testNo++), (int)(B.determinant() * 100) == (int)(-11130.32 * 100)); 
}

void test_adjoint()
{
    std::string passMessage = "test_adjoint() passes test ";
    int testNo = 1;
    double x[] = { 8, -7,  4,  6, 7, 
                  -1,  7, -1,  5, 9, 
                  -5, -6, -7, -7, 5, 
                   5,  1,  9, -8, 7, 
                  -7, -7,  1,  4, 8,};
    Matrix A(5, 5);
    A.fill(x);
    Matrix adj = A.adjoint();

    double y[] = { 21253,   2491,   6943,  -1007, -24857, 
                  -11306,  17404,  -6008,   3262,  -8786, 
                   -8627, -10823, -22655,  16369,  19561, 
                    7460,   5546, -13012, -13630,   7292, 
                    6052,  15988,  10156,   6742,   2194};
    
    Matrix B(5, 5);
    B.fill(y);
    
    assertAndPrint(passMessage + std::to_string(testNo++), adj == B);
 
}

/***
 * DO NOT RUN WITH MULTIPLE PROCESSES
*/
void test_solve_system_with_jacobi() 
{
    std::string passMessage = "test_solve_system_with_jacobi() passes test ";
    int testNo = 1;
    double a[] = {40,  7,  1,  3,  3, 
                   7, 30, 10,  7,  1, 
                   8, 10, 40,  2,  3, 
                   6,  2,  7, 90,  4, 
                   4,  3,  1,  8, 80};
    Matrix A(5, 5);
    A.fill(a);

    double b[] = { 6, 
                   8, 
                   7, 
                  10, 
                   6};
    Matrix B(5, 1);
    B.fill(b);

    Matrix x(5, 1);
    A.solveSystem(x, B);

    Matrix::setEqualityPrecision(0.00001f);

    assertAndPrint(passMessage + std::to_string(testNo++), A * x == B);

    Matrix::setEqualityPrecision(0.0f);
}

void test_solve_system_with_failsafe()
{
    std::string passMessage = "test_solve_system_with_failsafe() passes test ";
    int testNo = 1;
    double a[] = {10,  0,   3,  6, -6, 
                  -3, -7,   4, -7,  6, 
                   3,  8,  -2,  4, -8, 
                  -4,  2,  -5,  2,  4, 
                   1,  10, -2,  5, -2};
    Matrix A(5, 5);
    A.fill(a);

    double b[] = { 3, 
                   9, 
                   3, 
                  10, 
                   7};
    Matrix B(5, 1);
    B.fill(b);

    Matrix x(5, 1);
    A.solveSystem(x, B);

    Matrix::setEqualityPrecision(0.00001f);

    assertAndPrint(passMessage + std::to_string(testNo++), A * x == B);

    Matrix::setEqualityPrecision(0.0f);
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
    assertAndPrint("\n---Now running all test cases...", true);
    test_equality();
    test_inequality();
    test_copy_constructor();
    test_matrix_multiplication_same_size();
    test_matrix_multiplication_more_rows();
    test_matrix_multiplication_more_cols();
    test_matrix_scalar_multiplication();
    test_matrix_addition();
    test_matrix_subtraction();
    test_matrix_transpose_same_size();
    test_matrix_transpose_more_rows();
    test_matrix_transpose_more_cols();
    test_inverse();
    test_determinant();
    test_adjoint();
    //test_solve_system_with_jacobi(); // DO NOT RUN WITH MULTIPLE PROCESSES
    test_solve_system_with_failsafe();
    assertAndPrint("---All tests have passed!", true);
}

int main(int argc, char* argv[])
{
    signal(SIGSEGV, segfault_printer);
    processID = atoi(argv[1]);
    mainSimpi = new Simpi(processID, atoi(argv[2])); // argv[2]: # of processes being used
    Matrix::setSimpi(mainSimpi);
    runTests();
    delete mainSimpi; 
}

