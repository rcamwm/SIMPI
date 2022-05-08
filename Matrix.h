#ifndef MATRIX_H
#define MATRIX_H

#include "Simpi.h"
#include "Vector.h"

namespace SimpiNS
{
    class Matrix
    {
        private:
            int xdim, ydim;
            double* arr;
            std::string uniqueID;
            static Simpi* mainSimpi;

            // jacobi() helper functions
            void jacobiSaveInputs(int start, int end, Matrix* saveEq, Vector* constants, Vector* saveConst);
            void jacobiSwitchAndDivide(int start, int end, Vector* constants, Vector* solution, Vector *prev, int synchOffset);
            void jacobiFirstIteration(int start, int end, Vector* solution, Vector *prev, int synchOffset);
            void jacobiRemainingIterations(int start, int end, Vector* solution, Vector *prev, int synchOffset);
            void jacobiRestoreInputs(int start, int end, Matrix* saveEq, Vector* constants, Vector* saveConst);

        public:
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

            void solveSystem(Vector* constants, Vector* solution);
            void jacobi(Vector* constants, Vector* solution);
            void failSafe(Vector* constants, Vector* solution);
            bool isDiagonallyDominant();

            friend std::ostream& operator<<(std::ostream& out, const Matrix& m);      
            // overload +, -, *, /, ..... 
    };
}
#endif