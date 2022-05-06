#ifndef MATRIX_H
#define MATRIX_H

#include "Simpi.h"

class Matrix
{
    private:
        int xdim, ydim;
        double* arr;
        std::string unique_id;
        static Simpi* main_simpi;
        //simpi* mysimpi = NULL;  // for later reference

    public:
        class Vector;
        static void setSimpi(Simpi *s) { main_simpi = s; }

        Matrix(int x, int y);
        ~Matrix();

        int getSimpiID() const { return main_simpi->getID(); }
        int get_x() { return xdim; }
        int get_y() { return ydim; }
        double& get(int x, int y) { return arr[x + y * xdim]; }
        double get_algbera(int pos) { return arr[pos]; }
        void set(int pos, int val) { arr[pos] = val; }

        int determinant(double* A, int n, int order);
        void adjoint(double* A, double* adj, int order, int par_id, int par_count);
        void getCofactor(double* A, double* temp, int p, int q, int n, int order);
        
        void inverse(Matrix* inverse);
        void luDecomposition(Matrix* lower, Matrix* upper);
        void backward_substitution(float* b, float* x);
        void forward_substitution(float *b, float* x); 

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
        std::string unique_id;
        
    public:
        Vector(int a);
        ~Vector();

        int getSimpiID() const { return Matrix::main_simpi->getID(); }
        int getSize() const { return dim; }
        double& getRef(int pos) { return arr[pos]; }
        double getVal(int pos) const { return arr[pos]; }
        void set(int pos, double val) { arr[pos] = val; }
        
        friend std::ostream& operator<<(std::ostream& out, const Matrix::Vector& v);     
};

#endif