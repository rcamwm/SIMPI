#include "Vector.h"

namespace SimpiNS
{
    Simpi* Vector::mainSimpi;

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

    std::ostream& operator<<(std::ostream& out, const Vector& v)
    {
        
        if (v.getSimpiID() == 0)
        {
            out << std::endl;
            for (int i = 0; i < v.getSize(); i++)
                out << std::fixed << std::setprecision(2) << v.getVal(i) << ", " << std::endl;
        }
        return out;
    }
}