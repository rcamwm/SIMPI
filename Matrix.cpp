#include "Matrix.h"

namespace SimpiNS
{
    // Must be set before any Matrix objects are created
    Simpi* Matrix::mainSimpi;

    Matrix::Matrix(int x, int y)
    {
        // use mainSimpi to init the Matrix for all processes. The id is also in simp
        std::pair<std::string, double*> passBack(mainSimpi->createMatrix(x, y));
        uniqueID = passBack.first;
        arr = passBack.second;
        xdim = x;
        ydim = y;
    }

    Matrix::~Matrix()
    {
        mainSimpi->freeMatrix(uniqueID); // frees and unlinks memory
    }

    std::ostream& operator<<(std::ostream& out, const Matrix& m)
    {
        if (m.getSimpiID() == 0)
        {
            for (int i = 0; i < m.xdim; i++)
            {
                out << std::endl;
                for (int j = 0; j < m.ydim; j++)
                {
                    out << std::fixed << std::setprecision(2) << m.arr[i + j * m.xdim];
                    out << ", ";
                }
            }
            out << std::endl;
        }
        return out;
    }

    int Matrix::determinant()
    {
        if (!isSquareMatrix()) 
        {
            std::cout << "Invalid Matrix: Must be square matrix" << std::endl;
            exit(1);
        }
        return calculateDeterminant(this->arr, getX(), getY());
    }

   /**
    * Helper function to calculate the determinant of a Matrix.
    * Initially called using this->arr.
    * Computes by recursively calling itself to find minors and calculate cofactors
    * using calculateMinor()
    * 
    * @param A Array Matrix values whose determinant is being evaluated
    * @param n Order of A 
    * @param order Order of *this Matrix
    * @return Matrix's determinant
    * */
    int Matrix::calculateDeterminant(double* A, int n, int order)
    {
        if (n == 1) // Base case : if Matrix contains single element
            return A[0];

        int det = 0; // result
        int sign = 1; // sign multiplier
        double m[order * order]; // Array of minors

        // Iterate for each element of top row
        for (int x = 0; x < n; x++) 
        {
            // Getting minor of matrixArray[0][x]
            calculateMinor(A, m, 0, x, n, order);
            det += sign * A[0 + x * order] * calculateDeterminant(m, n - 1, order);

            // terms are to be added with alternate sign
            sign = -sign;
        }
        return det;
    }

    void Matrix::adjoint(Matrix *adj)
    {
        if (!isSquareMatrix() || !adj->isSquareMatrix() || getX() != adj->getX()) 
        {
            std::cout << "Invalid Matrices: Must be equal sized square matrices" << std::endl;
            exit(1);
        }
        allocateAdjointWork(this->arr, adj->arr, getX());
        mainSimpi->synch();
    }

    void Matrix::allocateAdjointWork(double* A, double* adj, int order)
    {
        int processCount = mainSimpi->getProcessCount();
        int processID = mainSimpi->getID();

        if (order <= processCount)
        {
            int start = processID;
            int end = start + 1;
            if (processID < order)
                calculateAdjoint(A, adj, order, start, end);
        }
        else 
        {
            int work = order / processCount;
            int start = processID * work;
            int end = start + work;
            calculateAdjoint(A, adj, order, start, end);

            int leftoverWork = order % processCount;
            if (leftoverWork != 0)
            {
                start = (work * processCount) + processID;
                end = start + 1;
                if (processID < leftoverWork)
                    calculateAdjoint(A, adj, order, start, end);
            }         
        }
    }

    /**
    * Helper function to calculate the adjoint of a Matrix.
    * 
    * @param A Array of matrix to find the adjoint of 
    * @param adj Array of matrix where adjoint will be stored
    * @param order Order of matrices A and adj
    * @param start Starting row index of A to evaluate in this process
    * @param end Ending row index of A to evaluate in this process
    * */
    void Matrix::calculateAdjoint(double* A, double* adj, int order, int start, int end)
    {
        if (order == 1) 
        {
            adj[0] = 1;
            return;
        }

        // m is used to store minors of A[][]
        int sign = 1;
        double m[order * order];

        for (int i = 0; i < order; i++) 
        {
            for (int j = start; j < end; j++) 
            {
                // Get minor of A[i][j]
                calculateMinor(A, m, i, j, order, order);

                // sign of adj[j][i] positive if sum of row
                // and column indexes is even.
                sign = ((i + j) % 2 == 0) ? 1 : -1;

                // Interchanging rows and columns to get the
                // transpose of the cofactor Matrix
                adj[j + i * order] = (sign) * (calculateDeterminant(m, order - 1, order));
            }
        }
    }

    /**
    * Helper function to calculate the minor of a Matrix.
    * 
    * @param A Array of matrix to find the minor of 
    * @param m Array where matrix elements of the minor are placed into
    * @param i The column in currentArray to find the minor of; Range: (0, ..., n)
    * @param j The row in currentArray to find the minor of; Range: (0, ..., n)
    * @param n The order of the elements in A that are being considered (if n < order, elements of higher order are 0)
    * @param order Matrix's total order
    * */
    void Matrix::calculateMinor(double* A, double* m, int i, int j, int n, int order)
    {
        int mRow = 0; // row of m array
        int mCol= 0; // column of m array

        // Looping for each element of the Matrix
        for (int row = 0; row < n; row++) 
        {
            for (int col = 0; col < n; col++) 
            {
                //  Copying into temporary array only those element
                //  which are not in given row and column
                if (row != i && col != j) 
                {
                    m[(mRow) + (mCol++) * order] = A[row + col * order];

                    // Row is filled, so increase row index and
                    // reset col index
                    if (mCol == n - 1) 
                    {
                        mCol = 0;
                        mRow++;
                    }
                }
            }
        }
    }

    /*
    This function calculates the lower and upper triangular matrices of a square nxn Matrix
    It requires an input of 2 empty nxn matrices that are modified 
    A = LU
    */
    void Matrix::luDecomposition(Matrix* lower, Matrix* upper) {

        if (!isSquareMatrix()) 
        {
            std::cout << "Invalid Matrix";
            exit(1);
        }

        for (int i = 0; i < getX(); i++) 
        {
            // Calculate work per parallel process
            // Has to be calculated on every loop iteration as the inner loop is decrementing
            int processCount = mainSimpi->getProcessCount();
            int processID = mainSimpi->getID();
            int total = xdim - i;
            if (processCount > total)
                processCount = total;

            int rpp = total / processCount;
            int start = rpp * mainSimpi->getID() + i;
            int end = start + rpp;

            // Upper Triangular 
            for (int k = start; k < end; k++) 
            {
                if (k >= getX())
                    break;

                // Summation of L(i, j) * U(j, k) 
                float sum = 0;
                for (int j = 0; j < i; j++) 
                    sum += (lower->get(i, j) * upper->get(j, k)); 

                // Evaluating U(i, k) 
                upper->get(i,k) = get(i,k) - sum; 

            }

            // Calculate and execute which processes take the leftover work 
            if (total % processCount != 0) 
            {
                int leftover = total % processCount;
                if (processID < leftover) 
                {
                    processID += (xdim - leftover);
                    int start = processID;
                    int end = start + 1;
                    for (int a = start; a < end; a++) 
                    {
                        // Summation of L(i, j) * U(j, k) 
                        float sum = 0;
                        for (int j = 0; j < i; j++) 
                            sum += (lower->get(i, j) * upper->get(j, a));

                        // Evaluating U(i, k) 
                        upper->get(i, a) = get(i, a) - sum;
                    }
                }
            }

            mainSimpi->synch();

            total = getX() - i;
            processID = mainSimpi->getID();
            processCount = mainSimpi->getProcessCount();

            // Lower Triangular
            for (int k = start; k < end; k++) 
            {
                if (k >= getX())
                    break;
                
                if (i == k)
                    lower->get(i, i) = 1; // Diagonal as 1 
                else 
                {
                    // Summation of L(k, j) * U(j, i)
                    float sum = 0;
                    for (int j = 0; j < i; j++) 
                    sum += (lower->get(k, j) * upper->get(j, i));
                    // Evaluating L(k, i)
                    lower->get(k, i) = ((get(k, i) - sum) / upper->get(i, i));
                }
            }

            // Calculate and execute which processes take the leftover work 
            if (total % processCount != 0) 
            {
                int leftover = total % processCount;
                if (processID < leftover) 
                {
                    processID += (getX() - leftover);
                    int start = processID;
                    int end = start + 1;
                    for (int a = start; a < end; a++) 
                    {
                        if (i == a)
                            lower->get(i, i) = 1; // Diagonal as 1
                        else 
                        {
                            // Summation of L(k, j) * U(j, i) 
                            float sum = 0;
                            for (int j = 0; j < i; j++) 
                            sum += (lower->get(a, j) * upper->get(j, i));

                            // Evaluating L(k, i) 
                            lower->get(a, i) = (get(a, i)-sum) / upper->get(i, i);
                        }
                    }
                }
            }
            mainSimpi->synch();
        }
        return;
    }

    /*
    This method calculates the inverse of a Matrix by using its LU Decomposition
    A = LU
    LZ = B
    LX = Z
    B corresponds to individual columns of an nxn identity Matrix
    X represents each corresponding column of the inverse Matrix
    */
    void Matrix::inverse(Matrix* inv) 
    {
        if (!isSquareMatrix() || !inv->isSquareMatrix() || getX() != inv->getX()) 
        {
            std::cout << "Invalid Matrices: Must be equal sized square matrices" << std::endl;
            exit(1);
        }

        //Solve for lower and upper matrices
        Matrix* upper = new Matrix(getX(), getY());
        Matrix* lower = new Matrix(getX(), getY());
        luDecomposition(lower,upper);
        mainSimpi->synch();

        //Create Identity nxn Matrix
        Matrix* identity = new Matrix(getX(), getY());
        if (mainSimpi->getID() == 0) 
        {
            for (int i = 0; i<getX(); i++) 
            {
                for (int j = 0; j<getX(); j++)            
                    (i == j) ? identity->get(i,j) = 1 : identity->get(i,j) = 0; 
            }
        }

        mainSimpi->synch();

        // Calculate columns per parallel process
        int processCount = mainSimpi->getProcessCount();
        if (processCount > getX()) 
            processCount = getX();
        
        int cpp = getX() / processCount;
        int start = cpp*mainSimpi->getID();
        int end = start + cpp;

        // Initialize necessary arrays for future calculations
        // Each array is local to its own process
        float identityCol[getX()];
        float zCol[getX()];
        float solutionCol[getX()];

        for (int a = start; a<end; a++) 
        {
            //Get individual columns of identity Matrix
            for (int b = 0; b < getX(); b++)         
                identityCol[b] = identity->get(b,a);
            
            //Reset Z column to solve for again
            for (int d = 0; d<getX(); d++)
                zCol[d] = 0;
            
            //Solve LZ = I
            lower->forwardSubstitution(identityCol, zCol);

            //Reset X column to solve for again
            for (int d = 0; d < getX(); d++)
                solutionCol[d] = 0;

            //Solve UX = Z
            upper->backwardSubstitution(zCol, solutionCol);

            //Input X column to corresponding columnn in final inverse Matrix
            for (int c = 0; c < getX(); c++)
                inv->get(c,a) = solutionCol[c];
        }

        // Calculate and execute which processes take the leftover rows 
        // ex. with 3 processes and a 10x10 Matrix:
        // 0-3 is process 0
        // 3-6 is process 1
        // 6-9 is process 2
        // 9 is the leftover column that is taken by process 0
        int processID = mainSimpi->getID();
        if (getX() % processCount != 0) 
        {
            int leftover = getX() % processCount;
            if (processID < leftover) 
            {
                processID += (getX() - leftover);
                int start = processID;
                int end = start + 1;
                for (int a = start; a < end; a++) 
                {
                    //Get individual columns of identity Matrix
                    for (int b = 0; b < getX(); b++)
                        identityCol[b] = identity->get(b, a);
                    
                    //Reset Z column to solve for again
                    for (int d = 0; d < getX(); d++) 
                        zCol[d] = 0;
                    
                    //Solve LZ = I
                    lower->forwardSubstitution(identityCol, zCol);

                    //Reset X column to solve for again
                    for (int d = 0; d < getX(); d++)
                        solutionCol[d] = 0;
                    
                    //Solve UX = Z
                    upper->backwardSubstitution(zCol, solutionCol);

                    //Input X column to corresponding columnn in final inverse Matrix
                    for (int c = 0; c < getX(); c++) 
                        inv->get(c, a) = solutionCol[c];
                    
                }
            }
        }  

        mainSimpi->synch();
        return;
    }

    /*
    This is a helper function to calculate the solutions of a lower triangular Matrix
    */
    void Matrix::forwardSubstitution(float *b, float* x)
    {
        double suma;
        for(int i = 0; i < getX(); i = i + 1)
        {
            suma = 0;
            for(int j = 0; j < i; j = j + 1)
                suma = suma + get(i,j) * x[j];

            x[i] = (b[i] - suma) / get(i, i);
        }
    }

    /*
    This is a helper function to calculate the solutions of an upper triangular Matrix
    */
    void Matrix::backwardSubstitution(float* b, float* x)
    {
        double suma;
        for(int i = getX() - 1; i >= 0; i = i - 1)
        {
            suma=0;
            for(int j = getX() - 1; j > i; j = j - 1)
                suma = suma + get(i, j) * x[j];
    
            x[i] = (b[i] - suma) / get(i, i);
        }
    }

    /*
    Solves a linear system of equations
    If the Matrix is diagonally dominant, jacobi-iterative method is used
    else, the inverse mutliplication method is used
    takes in a Vector of constants that each row is to be solved for and a solution Vector of 0s in which the solution will be written in
    void -> void
    */
    void Matrix::solveSystem(Vector *constants, Vector* solution)
    {
        bool dd = isDiagonallyDominant(); 
        mainSimpi->synch();
        if (dd)
        {
            if (mainSimpi->getID() == 0) { std::cout << "jacobi" << std::endl; }
            mainSimpi->synch();
            jacobi(constants, solution);
        }
        else
        {
            if (mainSimpi->getID() == 0) { std::cout << "failsafe" << std::endl; }        
            mainSimpi->synch();
            failSafe(constants, solution);
        }
        mainSimpi->synch();
    }

    /*
        Solves a linear system of equations in parallel if the Matrix is diagonally dominant
        outputs solution to a Vector of 0s passed into it
        void->void
    */
    void Matrix::jacobi(Vector* constants, Vector* solution) 
    {
        int id = mainSimpi->getID();
        int processCount = mainSimpi->getProcessCount();
        if (processCount > constants->getSize()) 
            processCount = constants->getSize();
        
        Matrix* saveEq = new Matrix(getX(), getY() + 1); // save equations from modification
        Vector* saveConst = new Vector(constants->getSize()); // saves input vector
        Vector *prev = new Vector(constants->getSize()); // shared mem containing a copy of values

        int work = constants->getSize() / processCount; 
        int start = 0;
        int end = 0;
        if (id < processCount) // Extra processes get no work => start and end stay at 0
        {
            start = id * work;
            end = start + work;
        }
        int remainingWork = constants->getSize() % processCount;
        if (id == processCount - 1) // Last active process gets all remaining work
            end += remainingWork;
        // synch() is called (end - start) number of times in several jacobi helper functions
        // synchOffset makes sure that all processes stay in sync even if some processes are doing more work than others
        int synchOffset = (work + remainingWork) - end + start;
        mainSimpi->synch();

        jacobiSaveInputs(start, end, saveEq, constants, saveConst);
        mainSimpi->synch();

        // Switches diagonals with solution elements,
        // then divides each matrix row by their original diagonal element 
        jacobiSwitchAndDivide(start, end, constants, solution, prev, synchOffset);
        mainSimpi->synch();

        // First iteration with substitution 1
        jacobiFirstIteration(start, end, solution, prev, synchOffset);
        mainSimpi->synch();

        // Repeat iterations with calculated results
        jacobiRemainingIterations(start, end, solution, prev, synchOffset);
        mainSimpi->synch();

        jacobiRestoreInputs(start, end, saveEq, constants, saveConst);
        mainSimpi->synch();

        return;
    }

    void Matrix::jacobiSaveInputs(int start, int end, Matrix* saveEq, Vector* constants, Vector* saveConst)
    {
        for (int i = start; i < end; i++) 
        {
            for (int j = 0; j < getY(); j++) 
                saveEq->get(i,j) = get(i, j);

            saveEq->get(i,getY() + 1) = constants->getRef(i);
            saveConst->getRef(i) = constants->getRef(i);
        }
    }

    /*
        In each row, the diagonal element of the matrix is switched with the corresponding element from the solution vector.
        Every element of that matrix row is then divided by the value that was just placed in the solution vector.
        Each element of that matrix row that was not switched is also multiplied by -1.
    */
    void Matrix::jacobiSwitchAndDivide(int start, int end, Vector* constants, Vector* solution, Vector *prev, int synchOffset)
    {
        for (int i = start; i < end; i++) 
        {
            double temp = get(i, i);
            get(i, i) = constants->getRef(i);
            constants->getRef(i) = temp;
            for (int j = 0; j < getY(); j++) 
            {
                if (j != i)
                    get(i, j) *= -1;
        
                get(i, j) /= constants->getRef(i);
                //mainSimpi->synch();
            }
            prev->getRef(i) = 1.0;
            solution->getRef(i) = constants->getRef(i);
            mainSimpi->synch();
        }
        mainSimpi->synchExtraCycles(synchOffset);
    }

    void Matrix::jacobiFirstIteration(int start, int end, Vector* solution, Vector *prev, int synchOffset)
    {
        for (int i = start; i < end; i++) 
        {
            double rowSum = 0;
            for (int j = 0; j < getY(); j++) 
            {
                if (j == i)
                    rowSum += get(i, j);
                else
                    rowSum += (get(i, j) * prev->getRef(j));
                
                //mainSimpi->synch();
            }
            solution->getRef(i) = rowSum;
            mainSimpi->synch();
        }
        mainSimpi->synchExtraCycles(synchOffset);
    }

    void Matrix::jacobiRemainingIterations(int start, int end, Vector* solution, Vector *prev, int synchOffset)
    {
        for (int k = 0; k < 1000; k++)
        {
            for (int i = start; i < end; i++)
            {
                //save prev value for comparision
                prev->getRef(i) = solution->getRef(i);
                mainSimpi->synch();
            }
            mainSimpi->synchExtraCycles(synchOffset);
            for (int i = start; i < end; i++)
            {
                double rowSum = 0;
                for (int j = 0; j < getY(); j++)
                {
                    if (j == i) {
                        rowSum += get(i, j);
                    } else {
                        rowSum += (get(i, j) * prev->getRef(j));
                    }
                }
                solution->getRef(i) = rowSum;
                //mainSimpi->synch();
            }
            //wait at end of each loop for all processes before beginning next iteration
            mainSimpi->synch();
        }
    }

    void Matrix::jacobiRestoreInputs(int start, int end, Matrix* saveEq, Vector* constants, Vector* saveConst)
    {
        for (int i = start; i < end; i++) 
        {
            for (int j = 0; j < getY(); j++) 
                get(i, j) = saveEq->get(i,j);
            
            constants->getRef(i) = saveConst->getRef(i);
        }
    }

    /*
    Method to solve a system of linear equations if the system is not diagonally dominant
    Uses the inverse-mutliplication method. 
    void->void
    */
    void Matrix::failSafe(Vector* constants, Vector* solution)
    {
        Matrix* inv = new Matrix(getX(), getY());
        inverse(inv);
        mainSimpi->synch();

        int processCount = mainSimpi->getProcessCount();
        int id = mainSimpi->getID();
        
        int vectorSize = constants->getSize();
        int work = vectorSize / processCount; 
        int start = 0;
        int end = 0;
        if (id < processCount) // Extra processes get no work => start and end stay at 0
        {
            start = id * work;
            end = start + work;
        }
        int remainingWork = constants->getSize() % processCount;
        if (id == processCount - 1) // Last active process gets all remaining work
            end += remainingWork;
        
        double sol;
        for(int i = start; i < end; i++)
        {
            sol = 0;
            for(int j = 0; j < vectorSize; j++)
            {
                sol += (inv->get(i, j)*constants->getRef(j));
            }
            solution->getRef(i) = sol;
        }
        mainSimpi->synch();
        return;
    }

    /*
    Checks if a square Matrix is diagonally dominant
    (diagonal terms are greater-than or equal to sum of the rest of their row)
    none->bool
    */
    bool Matrix::isDiagonallyDominant()
    {
        if (!isSquareMatrix())
            return false;

        for(int i = 0; i < getX(); i ++)
        {
            double sq;
            double rest = 0;
            for(int j = 0; j < getY(); j++)
                (i == j) ? sq = get(i, j) : rest += get(i,j);
            
            if (sq < rest)
                return false;
        }
        return true;
    }

    void Matrix::initializeArrayToZero(double *A, int size)
    {
        for (int i = 0; i < size; i++)
            A[i] = 0;
    }

}
