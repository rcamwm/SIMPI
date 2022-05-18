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

            void initializeArrayToZero(double *A, int size);

            // determinate() and adjoint() helper functions
            int calculateDeterminant(double* A, int n, int order);
            void allocateAdjointWork(double* A, double* adj, int order);
            void calculateAdjoint(double* A, double* adj, int order, int start, int end);
            void calculateMinor(double* currentArray, double* minorArray, int i, int j, int n, int order);

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
            bool isSquareMatrix() { return getX() == getY(); }

            int determinant();
            void adjoint(Matrix* adj);
            
            void inverse(Matrix* inv);
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