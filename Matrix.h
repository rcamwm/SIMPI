#ifndef MATRIX_H
#define MATRIX_H

#include <iomanip>
#include <random>
#include <math.h>

#include "Simpi.h"

namespace SimpiNS
{
    class Matrix
    {
        private:
            int rows, cols;
            double* arr;
            std::string uniqueID;
            static Simpi* mainSimpi;
            static double equalityPrecision;
            static int printPrecision;

            void initializeArrayToZero(double *A, int size);
            void copyElements(const Matrix &m, int start, int end, bool moreRows);

            // determinate() and adjoint() helper functions
            double calculateDeterminant(double* A, int n, int order);
            void allocateAdjointWork(double* A, double* adj, int order);
            void calculateAdjoint(double* A, double* adj, int order, int start, int end);
            void calculateMinor(double* currentArray, double* minorArray, int i, int j, int n, int order);

            // jacobi() helper functions
            void jacobiSaveInputs(int start, int end, Matrix &saveEq, Matrix &constants, Matrix &saveConst);
            void jacobiSwitchAndDivide(int start, int end, Matrix &constants, Matrix &solution, Matrix &prev, int synchOffset);
            void jacobiFirstIteration(int start, int end, Matrix &solution, Matrix &prev, int synchOffset);
            void jacobiRemainingIterations(int start, int end, Matrix &solution, Matrix &prev, int synchOffset);
            void jacobiRestoreInputs(int start, int end, Matrix &saveEq, Matrix &constants, Matrix &saveConst);

            // algebra helper functions
            void determineEquality(Matrix &comparand, int start, int end, bool moreRows, bool* eqValue);
            bool *getSharedBool(int &fd);
            void calculateProduct(Matrix &B, Matrix* C, int start, int end, bool moreRows);
            void calculateScalarProduct(double lambda, Matrix* product, int start, int end, bool moreRows);
            void calculateSum(const Matrix &B, Matrix* sum, int start, int end, bool moreRows);
            void calculateDifference(const Matrix &B, Matrix* sum, int start, int end, bool moreRows);

        public:
            static void setSimpi(Simpi *s);
            static void setEqualityPrecision(double d);
            static void setPrintPrecision(int p);

            Matrix(int rowCount, int colCount);
            Matrix(const Matrix &m);
            ~Matrix();

            int getSimpiID() const { return mainSimpi->getID(); }
            int getRows() { return rows; }
            int getCols() { return cols; }
            double& getRef(int row, int col) { return arr[row + (col * rows)]; }
            double getVal(int row, int col) const { return arr[row + (col * rows)]; }
            bool isSquareMatrix() { return rows == cols; }
            void fill(double *fillArray);
            void fillRandom(int min, int max);

            double determinant();
            void adjoint(Matrix &adj);
            
            void inverse(Matrix &inv);
            void luDecomposition(Matrix* lower, Matrix* upper);
            void backwardSubstitution(float* b, float* x);
            void forwardSubstitution(float *b, float* x); 

            void solveSystem(Matrix &solution, Matrix &constants);
            void jacobi(Matrix &solution, Matrix &constants); // TODO: fix, synch() getting stuck in jacobiRemainingIterations()
            void failSafe(Matrix &solution, Matrix &constants);
            bool isDiagonallyDominant();

            friend std::ostream& operator<<(std::ostream& out, const Matrix& m);

            bool equals(Matrix &comparand);
            friend bool operator==(Matrix &lhs, Matrix &rhs);
            friend bool operator!=(Matrix &lhs, Matrix &rhs);  

            Matrix &multiply(Matrix &operand);
            friend Matrix &operator*(Matrix &lhs, Matrix &rhs);
            friend void operator*=(Matrix &lhs, Matrix &rhs);
            
            Matrix &scalarMultiply(double operand);
            friend Matrix &operator*(double lhs, Matrix &rhs);
            friend Matrix &operator*(Matrix &lhs, double rhs);
            friend void operator*=(Matrix &lhs, double rhs);

            Matrix &add(Matrix &operand);
            friend Matrix &operator+(Matrix &lhs, Matrix &rhs);
            friend void operator+=(Matrix &lhs, Matrix &rhs);
            
            Matrix &subtract(Matrix &operand);
            friend Matrix &operator-(Matrix &lhs, Matrix &rhs);
            friend void operator-=(Matrix &lhs, Matrix &rhs); 
    };
}
#endif