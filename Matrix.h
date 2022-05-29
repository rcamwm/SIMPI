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
            void singleCellWorkDivision(int cellCount, int startEnd[]);

            // determinate() and adjoint() helper functions
            double calculateDeterminant(double* A, int n, int order);
            void allocateAdjointWork(double* A, double* adj, int order);
            void calculateAdjoint(double* A, double* adj, int order, int start, int end);
            void calculateMinor(double* currentArray, double* minorArray, int i, int j, int n, int order);

            // jacobi() helper functions
            void jacobiSaveInputs(int start, int end, Matrix &saveEq, Matrix &B, Matrix &saveB);
            void jacobiSwitchAndDivide(int start, int end, Matrix &B, Matrix &x, Matrix &prev, int synchOffset);
            void jacobiFirstIteration(int start, int end, Matrix &x, Matrix &prev, int synchOffset);
            void jacobiRemainingIterations(int start, int end, Matrix &x, Matrix &prev, int synchOffset);
            void jacobiRestoreInputs(int start, int end, Matrix &saveEq, Matrix &B, Matrix &saveB);

            // algebra helper functions
            void determineEquality(Matrix &B, int start, int end, bool* eqValue);
            bool *getSharedBool(int &fd, std::string sharedMemoryName);
            void calculateProduct(Matrix &B, Matrix* C, int start, int end);
            void calculateTranspose(Matrix* A_T, int start, int end);

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
            int getRow(int index) { return index % rows; }
            int getCol(int index) { return index / rows; }
            bool isSquareMatrix() { return rows == cols; }
            void fill(double *fillArray);
            void fillRandom(int min, int max);

            double determinant();
            void adjoint(Matrix &adj);
            
            void inverse(Matrix &inv);
            void luDecomposition(Matrix* L, Matrix* U);
            void backwardSubstitution(float* b, float* x);
            void forwardSubstitution(float *b, float* x); 

            void solveSystem(Matrix &x, Matrix &B);
            void jacobi(Matrix &x, Matrix &B); // TODO: fix, synch() getting stuck in jacobiRemainingIterations()
            void failSafe(Matrix &x, Matrix &B);
            bool isDiagonallyDominant();

            friend std::ostream& operator<<(std::ostream& out, const Matrix& m);

            bool equals(Matrix &B);
                friend bool operator==(Matrix &lhs, Matrix &rhs);
                friend bool operator!=(Matrix &lhs, Matrix &rhs);  
            Matrix &multiply(Matrix &B);
                friend Matrix &operator*(Matrix &lhs, Matrix &rhs);
                friend void operator*=(Matrix &lhs, Matrix &rhs);
            Matrix &scalarMultiply(double lambda);
                friend Matrix &operator*(double lhs, Matrix &rhs);
                friend Matrix &operator*(Matrix &lhs, double rhs);
                friend void operator*=(Matrix &lhs, double rhs);
            Matrix &add(Matrix &B);
                friend Matrix &operator+(Matrix &lhs, Matrix &rhs);
                friend void operator+=(Matrix &lhs, Matrix &rhs);
            Matrix &subtract(Matrix &B);
                friend Matrix &operator-(Matrix &lhs, Matrix &rhs);
                friend void operator-=(Matrix &lhs, Matrix &rhs); 
            Matrix &transpose();
    };
}
#endif