#include <signal.h>
#include <string.h>

#include "Simpi.h"
#include "Matrix.h"
#include "Vector.h"
using namespace SimpiNS;

int processID;
Simpi *mainSimpi;

/*
    userDefinedActivity() is to be filled in by the user of this program.
    All necessary initialization occurs in main() before this function is called.
    Can use Matrix and Vector objects here and perform operations with them.
    Use mainSimpi->synch() to synchronize processes.
*/
void userDefinedActivity() 
{
    
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
    Matrix::setEqualityPrecision(0.00f);
    Matrix::setSimpi(mainSimpi);
    Vector::setSimpi(mainSimpi);
    userDefinedActivity();
    delete mainSimpi;
}

