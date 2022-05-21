#ifndef COL_VECTOR_H
#define COL_VECTOR_H

#include <math.h>
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

            static double equalityPrecision;
            void determineEquality(Vector &comparand, int start, int end, bool* eqValue);
            bool *getSharedBool(int &fd);

            
        public:
            static void setSimpi(Simpi *s) { mainSimpi = s; }
            static void setEqualityPrecision(double d) { equalityPrecision = d; }
            Vector(int a);
            ~Vector();
            void fill(double *fillArray);

            int getSimpiID() const { return mainSimpi->getID(); }
            int getSize() const { return dim; }
            double& getRef(int pos) { return arr[pos]; }
            double getVal(int pos) const { return arr[pos]; }
            void set(int pos, double val) { arr[pos] = val; }
            
            friend std::ostream& operator<<(std::ostream& out, const Vector& v); 

            bool equals(Vector &comparand);  
            friend bool operator==(Vector &lhs, Vector &rhs);
            friend bool operator!=(Vector &lhs, Vector &rhs);  
  
    };
}

#endif