#include <assert.h>
#include <signal.h>
#include <string.h>

#include "Simpi.h"
#include "Matrix.h"
#include "Vector.h"
using namespace SimpiNS;

int processID;
Simpi *mainSimpi;

void assertAndPrint(std::string message, int i, bool test)
{
    if (i != 0)
        assert(test);

    if (processID == 0)
    {
        std::cout << message;
        if (i != 0)
            std::cout << i;
        std::cout << std::endl;
    }     
}

void test_matrix_equality()
{
    std::string passMessage = "test_matrix_equality() passes test ";
    int testNo = 1;

    double x[] = {1, 2, 3,
                  4, 5, 6,
                  7, 8, 9};
    Matrix A(3, 3); A.fill(x);
    assertAndPrint(passMessage, testNo++, A == A);

    Matrix B(3, 3); B.fill(x);
    assertAndPrint(passMessage, testNo++, A == B);

    B.get(2, 2) += 0.001;
    assertAndPrint(passMessage, testNo++, !(A == B));

    double y[] = {1, 2, 3, 4, 5, 6, 7, 8, 9}; 
    Matrix C(1, 9);
    C.fill(y);
    assertAndPrint(passMessage, testNo++, !(A == C));
}

void test_matrix_inequality()
{
    std::string passMessage = "test_matrix_inequality_operator() passes test ";
    int testNo = 1;
    
    double x[] = {1, 2, 3,
                  4, 5, 6,
                  7, 8, 9};
    double y[] = {0.1, 2.0, 3.0,
                  4.0, 0.5, 6.0,
                  7.0, 8.0, 0.9};
    Matrix A(3, 3); A.fill(x);
    Matrix B(3, 3); B.fill(y);
    assertAndPrint(passMessage, testNo++, A != B);

    if (processID == 0)
        for (int i = 0; i < B.getCols(); i++) { B.get(i, i) *= 10; }
    mainSimpi->synch();
    assertAndPrint(passMessage, testNo++, !(A != B));

    Matrix C(2, 3);
    assertAndPrint(passMessage, testNo++, B != C);
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
    assertAndPrint(passMessage, testNo++, AB == C);

    A *= B;
    assertAndPrint(passMessage, testNo++, A == C);
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
    assertAndPrint(passMessage, testNo++, AB == C);

    A *= B;
    assertAndPrint(passMessage, testNo++, A == C);

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
    assertAndPrint(passMessage, testNo++, AB == C);

    A *= B;
    assertAndPrint(passMessage, testNo++, A == C);

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

    Matrix inv(5, 5);
    A.inverse(inv);

    double y[] = {-0.0133363,  0.0174213,  0.0078267,  0.0726889, -0.0337270, 
                   0.1305650,  0.0546668,  0.1280080,  0.0901617, -0.0713501, 
                  -0.1119250, -0.0519893, -0.0115341,  0.0244585,  0.0760187, 
                   0.1357830,  0.0388418, -0.0603824, -0.0103498,  0.0036216, 
                  -0.0940235,  0.0417425, -0.1842020, -0.1361780,  0.1753460};
    
    Matrix B(5, 5);
    B.fill(y);
    
    Matrix::setEqualityPrecision(0.000001f);
    assertAndPrint(passMessage, testNo, inv == B);

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
    assertAndPrint(passMessage, testNo++, (int)A.determinant() == -5399);

    double y[] = {-5.2, 1.0,  3.0,  7.2,
                   8.2, 7.0, -8.0, -7.3,
                  -8.5, 7.0,  1.0, -6.1,
                   9.9, 1.0,  7.0,  5.1};

    Matrix B(4, 4);
    B.fill(y);
    assertAndPrint(passMessage, testNo++, (int)(B.determinant() * 100) == (int)(-11130.32 * 100)); 
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

    Matrix adj(5, 5);
    A.adjoint(adj);


    double y[] = { 21253,   2491,   6943,  -1007, -24857, 
                  -11306,  17404,  -6008,   3262,  -8786, 
                   -8627, -10823, -22655,  16369,  19561, 
                    7460,   5546, -13012, -13630,   7292, 
                    6052,  15988,  10156,   6742,   2194};
    
    Matrix B(5, 5);
    B.fill(y);
    
    assertAndPrint(passMessage, testNo, adj == B);
 
}

/***
 * Copy-paste test function for failsafe(),
 * and just make sure that the Matrices are diagonally dominant
 * 
 * No use filling this in until jacobi() has been fixed to work with multiple processes().
 * 
*/
void test_solve_system_with_jacobi() {}

void test_solve_system_with_failsafe()
{
    Matrix A(4, 5);
    Vector B(5);
    Vector C(5);

    mainSimpi->synch();
    
    for (int y = 0; y < 5; y++)
    {
        for (int x = 0; x < 4; x++)
        {
            A.get(x, y) = rand() % 10 + 1;
        }
        B.set(y, rand() % 10 + 1);
    }
    mainSimpi->synch();

    for (int y = 0; y < 5; y++)
    {
        for (int x = 0; x < 4; x++)
            if (mainSimpi->getID() == 0 && x == y)
                A.get(x, y) /= 3;
    }
    mainSimpi->synch();

    for (int i = 0; i < 5; i++)
        C.set(i, 0);
    
    mainSimpi->synch();
    A.solveSystem(&B, &C);
    if (mainSimpi->getID() == 0) { std::cout << A << "times\n" << C << "equals\n" << B; }
    
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
    assertAndPrint("\nNow running all test cases...", 0, 0);
    test_matrix_equality();
    test_matrix_inequality();
    test_matrix_multiplication_same_size();
    test_matrix_multiplication_more_rows();
    test_matrix_multiplication_more_cols();
    test_inverse();
    test_determinant();
    test_adjoint();

    // test_solve_system_with_jacobi(); // Does nothing -> fix jacobi() first
    // test_solve_system_with_failsafe(); // Need to write comparison operators for Vector
    assertAndPrint("All tests have passed!", 0, 0);
}

int main(int argc, char* argv[])
{
    signal(SIGSEGV, segfault_printer);
    processID = atoi(argv[1]);
    mainSimpi = new Simpi(processID, atoi(argv[2])); // argv[2]: # of processes being used
    Matrix::setSimpi(mainSimpi);
    Vector::setSimpi(mainSimpi);
    runTests();
    delete mainSimpi; 
}

