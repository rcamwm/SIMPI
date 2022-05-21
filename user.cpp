#include <signal.h>
#include <string.h>

#include "Simpi.h"
#include "Matrix.h"
using namespace SimpiNS;

int processID;
Simpi *mainSimpi;

/**
 * userDefinedActivity() is to be filled in by the user of this program.
 * All necessary initialization occurs in main() before this function is called.
 * Can use Matrix and Vector objects here and perform operations with them.
 * Use mainSimpi->synch() to synchronize processes.
 */
void userDefinedActivity() 
{
    if (processID == 0) { std::cout << std::endl << "EXAMPLE DEMONSTRATION" << std::endl; }
    double x[] = {6, 1, 3,
                  2, 4, 8,
                  3, 7, 4};
    Matrix A(3,3);
    A.fill(x);
    if (processID == 0) { std::cout << "Now printing Matrix object A:"; }
    std::cout << A; // Doesn't need to be wrapped in if (processID == 0) statement
    if (processID == 0) { std::cout << std::endl; }

    Matrix inverseA(3,3);
    A.inverse(inverseA);
    if (processID == 0) { std::cout << "Now printing the inverse of Matrix object A:" << inverseA << std::endl; }

    Matrix I = A * inverseA;
    if (processID == 0) { 
        std::cout << "A times its inverse is equal to the identity matrix:" << I; 
        std::cout << "But notice that one of the 0's is negative." << std::endl << std::endl; }
    
    Matrix::setPrintPrecision(8);
    if (processID == 0) { std::cout << "Setting the print precision to 8 will show you why:" << I << std::endl; }

    double i[] = {1, 0, 0,
                  0, 1, 0,
                  0, 0, 1};
    Matrix identity(3,3);
    identity.fill(i);

    if (identity != I)
    {
        if (processID == 0) { std::cout << "Unfortunately, the calculated identity is != to the actual identity" << identity << std::endl; }

        Matrix::setEqualityPrecision(0.0000001f);
        if (identity == I && processID == 0)
        {
            std::cout << "But by changing the equality precision to 0.0000001f, " << std::endl;
            std::cout << "the equality operator now returns true for the two Matrix objects" << std::endl;
        }
    }
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
    Matrix::setSimpi(mainSimpi);
    userDefinedActivity();
    delete mainSimpi;
}

