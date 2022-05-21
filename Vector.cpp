#include "Vector.h"

namespace SimpiNS
{
    Simpi* Vector::mainSimpi;
    double Vector::equalityPrecision;

    Vector::Vector(int a)
    {
        // use simp and init the Matrix for all processes. The id is also in simp
        std::pair<std::string, double*> passBack(mainSimpi->createMatrix(1, a));
        uniqueID = passBack.first;
        arr = passBack.second;
        dim = a;
    }

    Vector::~Vector() 
    {
        // use mainSimpi for getting rid of the mem and unlink stuff
        mainSimpi->freeMatrix(uniqueID);
    }

    void Vector::fill(double *fillArray)
    {
        if (mainSimpi->getID() == 0)
        {
            int i = 0;
            for (int row = 0; row < dim; row++)
            { 
                getRef(row) = fillArray[i];
                i++;
            }
        }
        mainSimpi->synch();
    }

    std::ostream& operator<<(std::ostream& out, const Vector& v)
    {
        
        if (v.getSimpiID() == 0)
        {
            out << std::endl;
            for (int i = 0; i < v.getSize(); i++)
                out << std::fixed << std::setprecision(4) << v.getVal(i) << ", " << std::endl;
        }
        return out;
    }

    bool Vector::equals(Vector &comparand)
    {
        if (dim != comparand.dim)
            return false;
            
        int fd;
        bool *equalityBool = getSharedBool(fd); // Shared between all processes so the same value is always returned
        *equalityBool = true;
        mainSimpi->synch();

        int processCount = mainSimpi->getProcessCount();
        int processID = mainSimpi->getID();

        if (dim <= processCount)
        {
            int start = processID;
            int end = start + 1;
            if (processID < dim)
                determineEquality(comparand, start, end, equalityBool);
        }
        else 
        {
            int work = dim / processCount;
            int start = processID * work;
            int end = start + work;
            determineEquality(comparand, start, end, equalityBool);

            int leftoverWork = dim % processCount;
            if (leftoverWork != 0)
            {
                start = (work * processCount) + processID;
                end = start + 1;
                if (processID < leftoverWork)
                    determineEquality(comparand, start, end, equalityBool);
            }         
        }
        mainSimpi->synch();
        bool eqValue = *equalityBool;
        
        mainSimpi->synch();
        close(fd);
        shm_unlink("TEMP_SHARED_BOOL_FOR_EQUALITY_CHECK");
        munmap(equalityBool, sizeof(bool*));
        
        mainSimpi->synch();
        return eqValue;
    }

    bool operator==(Vector &lhs, Vector &rhs)
    {
        return lhs.equals(rhs);
    }

    bool operator!=(Vector &lhs, Vector &rhs)
    {
        return !(lhs.equals(rhs));
    }

    void Vector::determineEquality(Vector &comparand, int start, int end, bool* eqValue)
    {
        for (int i = start; i < end; i++)
        {
            if (fabs(this->getVal(i) - comparand.getVal(i)) > equalityPrecision)
                *eqValue = false;                    
            if (!*eqValue)
                return;                
        }
    }

    bool* Vector::getSharedBool(int &fd) 
    {
        bool *eq;
        if (mainSimpi->getID() == 0) 
        {
            fd = shm_open("TEMP_SHARED_BOOL_FOR_EQUALITY_CHECK", O_RDWR | O_CREAT, 0777); 
            ftruncate(fd, sizeof(bool*));
            eq = (bool*)mmap(NULL, sizeof(bool*), PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
            mainSimpi->synch();
        }
        else 
        {
            mainSimpi->synch();
            fd = shm_open("TEMP_SHARED_BOOL_FOR_EQUALITY_CHECK", O_RDWR, 0777);
            eq = (bool*)mmap(NULL, sizeof(bool*), PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
        }
        return eq;
    }
}