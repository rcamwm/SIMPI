#ifndef MATRIX_H
#define MATRIX_H

#include "Simpi.h"
#include "Vector.h"

namespace SimpiNS
{
    class Matrix
    {
        private:
            int rows, cols;
            double* arr;
            std::string uniqueID;
            static Simpi* mainSimpi;

            void initializeArrayToZero(double *A, int size);

            // multiply() helper function
            void calculateProduct(Matrix &B, Matrix* C, int start, int end, bool rowGreaterThanCol);

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

            Matrix(int rowCount, int colCount);
            ~Matrix();

            int getSimpiID() const { return mainSimpi->getID(); }
            int getRows() { return rows; }
            int getCols() { return cols; }
            double& get(int row, int col) { return arr[row + (col * rows)]; }
            bool isSquareMatrix() { return rows == cols; }

            int determinant();
            void adjoint(Matrix* adj);
            
            void inverse(Matrix* inv);
            void luDecomposition(Matrix* lower, Matrix* upper);
            void backwardSubstitution(float* b, float* x);
            void forwardSubstitution(float *b, float* x); 

            void solveSystem(Vector* constants, Vector* solution);
            void jacobi(Vector* constants, Vector* solution); // TODO: fix, synch() getting stuck in jacobiRemainingIterations()
            void failSafe(Vector* constants, Vector* solution);
            bool isDiagonallyDominant();

            friend std::ostream& operator<<(std::ostream& out, const Matrix& m);

            // TODO
            Matrix &multiply(Matrix &operand);
            friend Matrix &operator*(Matrix &lhs, Matrix &rhs);
            friend void operator*=(Matrix &lhs, Matrix &rhs);
            
            friend Matrix &operator*(int lhs, Matrix &rhs);
            friend Matrix &operator*(Matrix &lhs, int rhs);
            friend void operator*=(Matrix &lhs, int rhs);

            friend Matrix &operator+(Matrix &lhs, Matrix &rhs);
            friend void operator+=(Matrix &lhs, Matrix &rhs);
            
            friend Matrix &operator-(Matrix &lhs, Matrix &rhs);
            friend void operator-=(Matrix &lhs, Matrix &rhs);
            
            
            friend bool operator==(Matrix &lhs, Matrix &rhs);
            friend bool operator!=(Matrix &lhs, Matrix &rhs);   
    };
}
#endif