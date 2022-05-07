#ifndef COL_VECTOR_H
#define COL_VECTOR_H

#include "Simpi.h"

namespace SimpiNS 
{
    class Vector
    {
        private:
            int dim;
            double* arr;
            //Simpi* mysimpi = NULL;  // for later reference
            std::string uniqueID;
            static Simpi *mainSimpi;
            
        public:
            static void setSimpi(Simpi *s) { mainSimpi = s; }
            Vector(int a);
            ~Vector();

            int getSimpiID() const { return mainSimpi->getID(); }
            int getSize() const { return dim; }
            double& getRef(int pos) { return arr[pos]; }
            double getVal(int pos) const { return arr[pos]; }
            void set(int pos, double val) { arr[pos] = val; }
            
            friend std::ostream& operator<<(std::ostream& out, const Vector& v);     
    };
}

#endif