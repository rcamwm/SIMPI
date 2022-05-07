#ifndef MATRIX_H
#define MATRIX_H

#include "Simpi.h"

class Matrix
{
    private:
        int xdim, ydim;
        double* arr;
        std::string uniqueID;
        static Simpi* mainSimpi;
        //simpi* mysimpi = NULL;  // for later reference

    public:
        class Vector;
        static void setSimpi(Simpi *s) { mainSimpi = s; }

        Matrix(int x, int y);
        ~Matrix();

        int getSimpiID() const { return mainSimpi->getID(); }
        int getX() { return xdim; }
        int getY() { return ydim; }
        double& get(int x, int y) { return arr[x + y * xdim]; }
        double getAlgbera(int pos) { return arr[pos]; }
        void set(int pos, int val) { arr[pos] = val; }

        int determinant(double* A, int n, int order);
        void adjoint(double* A, double* adj, int order, int processID, int processCount);
        void getCofactor(double* A, double* temp, int p, int q, int n, int order);
        
        void inverse(Matrix* inverse);
        void luDecomposition(Matrix* lower, Matrix* upper);
        void backwardSubstitution(float* b, float* x);
        void forwardSubstitution(float *b, float* x); 

        void solveSystem(Matrix::Vector* constants, Matrix::Vector* solution);
        void jacobi(Matrix::Vector* constants, Matrix::Vector* solution);
        void failSafe(Matrix::Vector* constants, Matrix::Vector* solution);
        bool isDiagonallyDominant();

        friend std::ostream& operator<<(std::ostream& out, const Matrix& m);      
        // overload +, -, *, /, ..... 
};

class Matrix::Vector
{
    private:
        int dim;
        double* arr;
        Simpi* mysimpi = NULL;  // for later reference
        std::string uniqueID;
        
    public:
        Vector(int a);
        ~Vector();

        int getSimpiID() const { return Matrix::mainSimpi->getID(); }
        int getSize() const { return dim; }
        double& getRef(int pos) { return arr[pos]; }
        double getVal(int pos) const { return arr[pos]; }
        void set(int pos, double val) { arr[pos] = val; }
        
        friend std::ostream& operator<<(std::ostream& out, const Matrix::Vector& v);     
};

#endif