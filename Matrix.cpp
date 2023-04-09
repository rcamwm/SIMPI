#include "Matrix.h"

namespace SimpiNS
{
    /**
     * Sets a Simpi object for the entire Matrix class. The Simpi object will keep track of each 
     * Matrix object's shared memory between processes. This function MUST be called before 
     * any Matrix objects are created.
     * 
     * @param s A Simpi object to share between all Matrix objects
     */
    void Matrix::setSimpi(Simpi *s) { mainSimpi = s; }
    Simpi* Matrix::mainSimpi;

    /**
     * Sets precision to specific decimal place to account for potential rounding errors
     * when comparing two Matric objects. Equality between two Matrix objects is determined 
     * by subtracting each same index element. They are considered equal if the absolute 
     * value difference between them is less than the value of d
     * 
     * @param d the new value to compare two Matrix elements against
     */
    void Matrix::setEqualityPrecision(double d) { equalityPrecision = d; }
    double Matrix::equalityPrecision = 0.0f;

    /**
     * Sets the number of digits after the decimal places to show when printing a Matrix object
     * 
     * @param p how many decimal places to show for each element
     */
    void Matrix::setPrintPrecision(int p) { printPrecision = p; }
    int Matrix::printPrecision = 2;

    /**
     * Construct a new Matrix object
     * 
     * @param rowCount the number of rows of the Matrix (stacked vertically)
     * @param colCount the number of columns of the Matrix (placed horizontally)
     */
    Matrix::Matrix(int rowCount, int colCount)
    {
        // use mainSimpi to init the Matrix for all processes. The id is also in simp
        std::pair<std::string, double*> passBack(mainSimpi->createMatrix(rowCount, colCount));
        uniqueID = passBack.first;
        arr = passBack.second;
        rows = rowCount;
        cols = colCount;
    }

    /**
     * Construct a new Matrix object
     * 
     * @param m an existing Matrix whose contents will be copied to the new Matrix object
     */
    Matrix::Matrix(const Matrix &m)
    {
        std::pair<std::string, double*> passBack(mainSimpi->createMatrix(m.rows, m.cols));
        uniqueID = passBack.first;
        arr = passBack.second;
        rows = m.rows;
        cols = m.cols;

        int cellCount = rows * cols;
        int indices[4];
        singleCellWorkDivision(cellCount, indices);

        for (int i = indices[0]; i < indices[1]; i++) // Initial work
            arr[i] = m.arr[i];
        
        for (int i = indices[2]; i < indices[3]; i++) // Leftover work
            arr[i] = m.arr[i];

        mainSimpi->synch();
    }

    /**
     * Frees and unlinks memory that was allocated with Simpi.
     */
    Matrix::~Matrix()
    {
        mainSimpi->freeMatrix(uniqueID);
    }

    /**
     * This function divides work between processes for functions that work on
     * cells in a Matrix individually. This does not work for all functions.
     * Some use algorithms that calculate results for certain elements based on the results
     * of previously calculated elements. This is only intended for use with functions
     * where an element's result is completely independent from all other elements.
     * 
     * @param cellCount The number of elements in the target Matrix object's double array
     * @param indexArray an integer array of size 4 that stores the starting and ending indices 
     *                   of the target Matrix object's double array to work on. 
     *                   [0] and [1] are the starting and ending indices of work to do for this process 
     *                   and will be -1 if there is no work to do. 
     *                   [2] and [3] are the starting and ending indices of leftover work to do for this process 
     *                   and will be -1 if there is no leftover work to do.
     */
    void Matrix::singleCellWorkDivision(int cellCount, int indexArray[])
    {
        int processCount = mainSimpi->getProcessCount();
        int processID = mainSimpi->getID();
        if (cellCount <= processCount) // More processes than cells => each process gets one cell
        {
            if (processID < cellCount) 
            {
                indexArray[0] = processID;
                indexArray[1] = processID + 1;
            }
            else
            {
                indexArray[0] = -1;
                indexArray[1] = -1;
            }
            indexArray[2] = -1;
            indexArray[3] = -1;
        }
        else // Otherwise, evenly distribute initial and leftover work between processes
        {
            int work = cellCount / processCount;
            indexArray[0] = processID * work;
            indexArray[1] = indexArray[0] + work;

            int leftoverWork = cellCount % processCount;
            if (leftoverWork != 0 && processID < leftoverWork)
            {
                indexArray[2] = (work * processCount) + processID;
                indexArray[3] = indexArray[2] + 1;
            }    
            else
            {
                indexArray[2] = -1;
                indexArray[3] = -1;
            }     
        }
    }

    /**
    * Fills the matrix with the values of arr.
    * Each value of arr will be placed into the next open row position.
    * When a row is filled it will move onto the next column.
    * This will continue until all rows and columns are filled.
    * 
    * This function can be optimized in the future, 
    * but there's probably no use for it outside of test cases.
    * 
    * @param arr Values to be entered into this Matrix.
    *            Make sure that arr has at least (this->getRows * this->getCols) elements
    */
    void Matrix::fill(double *fillArray)
    {
        if (mainSimpi->getID() == 0)
        {
            int i = 0;
            for (int row = 0; row < rows; row++)
            {
                for (int col = 0; col < cols; col++)
                {
                    getRef(row, col) = fillArray[i];
                    i++;
                }
            }
        }
        mainSimpi->synch();
    }

    /**
     * Fills the matrix with random doubles between a certain range. 
     * 
     * @param min the minimum possible value for a randomly generated element
     * @param max the maximum possible value for a randomly generated element
     */
    void Matrix::fillRandom(int min, int max)
    {
        if (mainSimpi->getID() == 0)
        {
            std::random_device rd;
            std::default_random_engine eng(rd());
            std::uniform_real_distribution<double> distr(min, max);
            for (int row = 0; row < rows; row++)
            {
                for (int col = 0; col < cols; col++)
                {
                    getRef(row, col) = distr(eng);
                }
            }
        }
        mainSimpi->synch();
    }

    /**
     * The determinant is a scalar value based on the entries of a square matrix.
     * Calling Matrix must have equal row and col count. 
     * 
     * @return the determinant of the Matrix
     */
    double Matrix::determinant()
    {
        if (!isSquareMatrix()) 
        {
            std::cout << "Invalid Matrix: Must be square matrix" << std::endl;
            exit(1);
        }
        mainSimpi->synch();
        return calculateDeterminant(this->arr, getRows(), getCols());
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
    double Matrix::calculateDeterminant(double* A, int n, int order)
    {
        if (n == 1) // Base case : if Matrix contains single element
            return A[0];

        double det = 0; // result
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

    /**
     * The adjoint of a square matrix is the transpose of its cofactor matrix.
     * Calling Matrix must have equal row and col count.
     */
    Matrix &Matrix::adjoint()
    {
        Matrix *A = this;
        if (!A->isSquareMatrix()) 
        {
            std::cout << "Invalid Matrix: Must be square matrix" << std::endl;
            exit(1);
        }
        Matrix *adj = new Matrix(A->rows, A->cols);
        int indices[4];
        singleCellWorkDivision(getRows(), indices);
        calculateAdjoint(this->arr, adj->arr, getRows(), indices[0], indices[1]);
        calculateAdjoint(this->arr, adj->arr, getRows(), indices[2], indices[3]);

        mainSimpi->synch();
        return *adj;
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

    /**
     * This function calculates the lower and upper triangular matrices of a square matrix,
     * such that A = LU.
     * Calling Matrix must have equal row and col count.
     * 
     * @param L pointer to a Matrix of equal size as this Matrix (A) to be overwritten
     * @param U pointer to a Matrix of equal size as this Matrix (A) to be overwritten
     */
    void Matrix::luDecomposition(Matrix* L, Matrix* U) {

        if (!isSquareMatrix()) 
        {
            std::cout << "Invalid Matrix";
            exit(1);
        }

        for (int i = 0; i < getRows(); i++) 
        {
            // Calculate work per parallel process
            // Has to be calculated on every loop iteration as the inner loop is decrementing
            int processCount = mainSimpi->getProcessCount();
            int processID = mainSimpi->getID();
            int total = rows - i;
            if (processCount > total)
                processCount = total;

            int rpp = total / processCount;
            int start = rpp * mainSimpi->getID() + i;
            int end = start + rpp;

            // Upper Triangular 
            for (int k = start; k < end; k++) 
            {
                if (k >= getRows())
                    break;

                // Summation of L(i, j) * U(j, k) 
                float sum = 0;
                for (int j = 0; j < i; j++) 
                    sum += (L->getVal(i, j) * U->getVal(j, k)); 

                // Evaluating U(i, k) 
                U->getRef(i,k) = getVal(i,k) - sum; 

            }

            // Calculate and execute which processes take the leftover work 
            if (total % processCount != 0) 
            {
                int leftover = total % processCount;
                if (processID < leftover) 
                {
                    processID += (rows - leftover);
                    int start = processID;
                    int end = start + 1;
                    for (int a = start; a < end; a++) 
                    {
                        // Summation of L(i, j) * U(j, k) 
                        float sum = 0;
                        for (int j = 0; j < i; j++) 
                            sum += (L->getVal(i, j) * U->getVal(j, a));

                        // Evaluating U(i, k) 
                        U->getRef(i, a) = getVal(i, a) - sum;
                    }
                }
            }

            mainSimpi->synch();

            total = getRows() - i;
            processID = mainSimpi->getID();
            processCount = mainSimpi->getProcessCount();

            // Lower Triangular
            for (int k = start; k < end; k++) 
            {
                if (k >= getRows())
                    break;
                
                if (i == k)
                    L->getRef(i, i) = 1; // Diagonal as 1 
                else 
                {
                    // Summation of L(k, j) * U(j, i)
                    float sum = 0;
                    for (int j = 0; j < i; j++) 
                    sum += (L->getVal(k, j) * U->getVal(j, i));
                    // Evaluating L(k, i)
                    L->getRef(k, i) = ((getVal(k, i) - sum) / U->getVal(i, i));
                }
            }

            // Calculate and execute which processes take the leftover work 
            if (total % processCount != 0) 
            {
                int leftover = total % processCount;
                if (processID < leftover) 
                {
                    processID += (getRows() - leftover);
                    int start = processID;
                    int end = start + 1;
                    for (int a = start; a < end; a++) 
                    {
                        if (i == a)
                            L->getRef(i, i) = 1; // Diagonal as 1
                        else 
                        {
                            // Summation of L(k, j) * U(j, i) 
                            float sum = 0;
                            for (int j = 0; j < i; j++) 
                            sum += (L->getVal(a, j) * U->getVal(j, i));

                            // Evaluating L(k, i) 
                            L->getRef(a, i) = (getVal(a, i)-sum) / U->getVal(i, i);
                        }
                    }
                }
            }
            mainSimpi->synch();
        }
        return;
    }

    /**
     * This method calculates the inverse of a Matrix by using its LU Decomposition,
     * such that A = LU, where LZ = B and LX = Z,
     * and where B corresponds to individual columns of an identity Matrix of the same size
     * and where X represents each corresponding column of the inverse Matrix.
     * Calling Matrix must have equal row and col count
     */
    Matrix &Matrix::inverse() 
    {
        if (!isSquareMatrix()) 
        {
            std::cout << "Invalid Matrix: Must be square matrix" << std::endl;
            exit(1);
        }
        Matrix *inv = new Matrix(rows, cols);

        //Solve for lower and upper matrices
        Matrix* upper = new Matrix(getRows(), getCols());
        Matrix* lower = new Matrix(getRows(), getCols());
        luDecomposition(lower,upper);
        mainSimpi->synch();

        //Create Identity nxn Matrix
        Matrix* identity = new Matrix(getRows(), getCols());
        if (mainSimpi->getID() == 0) 
        {
            for (int i = 0; i<getRows(); i++) 
            {
                for (int j = 0; j<getRows(); j++)            
                    (i == j) ? identity->getRef(i,j) = 1 : identity->getRef(i,j) = 0; 
            }
        }

        mainSimpi->synch();

        // Calculate columns per parallel process
        int processCount = mainSimpi->getProcessCount();
        if (processCount > getRows()) 
            processCount = getRows();
        
        int cpp = getRows() / processCount;
        int start = cpp*mainSimpi->getID();
        int end = start + cpp;

        // Initialize necessary arrays for future calculations
        // Each array is local to its own process
        float identityCol[getRows()];
        float zCol[getRows()];
        float solutionCol[getRows()];

        for (int a = start; a<end; a++) 
        {
            //Get individual columns of identity Matrix
            for (int b = 0; b < getRows(); b++)         
                identityCol[b] = identity->getVal(b,a);
            
            //Reset Z column to solve for again
            for (int d = 0; d<getRows(); d++)
                zCol[d] = 0;
            
            //Solve LZ = I
            lower->forwardSubstitution(identityCol, zCol);

            //Reset X column to solve for again
            for (int d = 0; d < getRows(); d++)
                solutionCol[d] = 0;

            //Solve UX = Z
            upper->backwardSubstitution(zCol, solutionCol);

            //Input X column to corresponding columnn in final inverse Matrix
            for (int c = 0; c < getRows(); c++)
                inv->getRef(c,a) = solutionCol[c];
        }

        // Calculate and execute which processes take the leftover rows 
        // ex. with 3 processes and a 10x10 Matrix:
        // 0-3 is process 0
        // 3-6 is process 1
        // 6-9 is process 2
        // 9 is the leftover column that is taken by process 0
        int processID = mainSimpi->getID();
        if (getRows() % processCount != 0) 
        {
            int leftover = getRows() % processCount;
            if (processID < leftover) 
            {
                processID += (getRows() - leftover);
                int start = processID;
                int end = start + 1;
                for (int a = start; a < end; a++) 
                {
                    //Get individual columns of identity Matrix
                    for (int b = 0; b < getRows(); b++)
                        identityCol[b] = identity->getVal(b, a);
                    
                    //Reset Z column to solve for again
                    for (int d = 0; d < getRows(); d++) 
                        zCol[d] = 0;
                    
                    //Solve LZ = I
                    lower->forwardSubstitution(identityCol, zCol);

                    //Reset X column to solve for again
                    for (int d = 0; d < getRows(); d++)
                        solutionCol[d] = 0;
                    
                    //Solve UX = Z
                    upper->backwardSubstitution(zCol, solutionCol);

                    //Input X column to corresponding columnn in final inverse Matrix
                    for (int c = 0; c < getRows(); c++) 
                        inv->getRef(c, a) = solutionCol[c];
                    
                }
            }
        }  

        mainSimpi->synch();
        return *inv;
    }

    /*
    This is a helper function to calculate the solutions of a lower triangular Matrix
    */
    void Matrix::forwardSubstitution(float *b, float* x)
    {
        double suma;
        for(int i = 0; i < getRows(); i = i + 1)
        {
            suma = 0;
            for(int j = 0; j < i; j = j + 1)
                suma = suma + getVal(i,j) * x[j];

            x[i] = (b[i] - suma) / getVal(i, i);
        }
    }

    /*
    This is a helper function to calculate the solutions of an upper triangular Matrix
    */
    void Matrix::backwardSubstitution(float* b, float* x)
    {
        double suma;
        for(int i = getRows() - 1; i >= 0; i = i - 1)
        {
            suma=0;
            for(int j = getRows() - 1; j > i; j = j - 1)
                suma = suma + getVal(i, j) * x[j];
    
            x[i] = (b[i] - suma) / getVal(i, i);
        }
    }

    /**
     * Solves a linear system of equations Ax=B.
     * If the Matrix is diagonally dominant, the jacobi-iterative method is used,
     * if not, the inverse mutliplication method is used.
     * 
     * Calling Matrix must have equal row and col count.
     * 
     * @param x a single column Matrix with equal rows to calling Matrix that the solution will be written into
     * @param B a single column Matrix of constants with equal rows to calling Matrix that each row is to be solved for
     */
    void Matrix::solveSystem(Matrix &x, Matrix &B)
    {
        if (x.cols != 1 || x.rows != rows || B.cols != 1 || B.rows != rows)
        {
            std::cout << "Cannot solve linear system of equations for x(" << x.rows << ", " << x.cols << ") ";
            std::cout << "and B(" << B.rows << ", " << B.cols << ")" << std::endl;
            std::cout << "Parameters must be size x(1, " << rows << ") and B(1, " << rows << std::endl;
            exit(1);
        }
        bool dd = isDiagonallyDominant(); 
        mainSimpi->synch();

        if (dd)
        {
            mainSimpi->synch();
            jacobi(x, B);
        }
        else
        {
            mainSimpi->synch();
            failSafe(x, B);
        }
        mainSimpi->synch();
    }

    /**
     * Solves a linear system of equations Ax=B in parallel if the Matrix is diagonally dominant.
     * Calling Matrix must have equal row and col count.
     * 
     * @param x a single column Matrix that the solution will be written into
     * @param B a single column Matrix of constants that each row is to be solved for
    */
    void Matrix::jacobi(Matrix &x, Matrix &B) 
    {
        int id = mainSimpi->getID();
        int processCount = mainSimpi->getProcessCount();
        if (processCount > B.getRows()) 
            processCount = B.getRows();
        
        Matrix saveEq(getRows(), getCols() + 1); // save calling Matrix and B
        Matrix saveB(B.getRows(), 1); // saves B
        Matrix prev(B.getRows(), 1); // shared mem containing a copy of values

        int work = B.getRows() / processCount; 
        int start = 0;
        int end = 0;
        if (id < processCount) // Extra processes get no work => start and end stay at 0
        {
            start = id * work;
            end = start + work;
        }
        int remainingWork = B.getRows() % processCount;
        if (id == processCount - 1) // Last active process gets all remaining work
            end += remainingWork;
        // synch() is called (end - start) number of times in several jacobi helper functions
        // synchOffset makes sure that all processes stay in sync even if some processes are doing more work than others
        int synchOffset = (work + remainingWork) - end + start;
        mainSimpi->synch();

        jacobiSaveInputs(start, end, saveEq, B, saveB);
        mainSimpi->synch();
        std::cout << "ID " << id << ": saved" << std::endl;

        // Switches diagonals with x elements,
        // then divides each matrix row by their original diagonal element 
        jacobiSwitchAndDivide(start, end, B, x, prev, synchOffset);
        mainSimpi->synch();
        std::cout << "ID " << id << ": switched and divided" << std::endl;

        // First iteration with substitution 1
        jacobiFirstIteration(start, end, x, prev, synchOffset);
        mainSimpi->synch();
        std::cout << "ID " << id << ": iterated first" << std::endl;

        // Repeat iterations with calculated results
        jacobiRemainingIterations(start, end, x, prev, synchOffset);
        mainSimpi->synch();
        std::cout << "ID " << id << ": iterated remaining" << std::endl;

        jacobiRestoreInputs(start, end, saveEq, B, saveB);
        mainSimpi->synch();
        std::cout << "ID " << id << ": restored" << std::endl;

        return;
    }

    /**
     * Saves inputs to jacobi() to restore later.
     * 
     * @param start the first index for this process to work on
     * @param end the last index for this process to work on
     * @param saveEq A Matrix that the calling Matrix Object (A) as well as B is being saved to
     * @param B a single column Matrix of constants that jacobi() is solving for
     * @param saveB A single column Matrix that B is being saved to
     */
    void Matrix::jacobiSaveInputs(int start, int end, Matrix &saveEq, Matrix &B, Matrix &saveB)
    {
        for (int i = start; i < end; i++) 
        {
            for (int j = 0; j < getCols(); j++) 
                saveEq.getRef(i,j) = getVal(i, j);

            saveEq.getRef(i,getCols() + 1) = B.getVal(i, 0);
            saveB.getRef(i, 0) = B.getVal(i, 0);
        }
    }

    /**
     * In each row, the diagonal element of the matrix is switched with the corresponding element from x. 
     * Every element of that matrix row is then divided by the value that was just placed in x.
     * Each element of that matrix row that was not switched is also multiplied by -1.
     * 
     * @param start the first index for this process to work on
     * @param end the last index for this process to work on
     * @param B a single column Matrix of constants that each row is to be solved for
     * @param x a single column Matrix that the solution will be written into
     * @param prev shared memory to be written to for later
     * @param synchOffset number of times for this process to call synch if it's not doing the maximum amount of work
     */
    void Matrix::jacobiSwitchAndDivide(int start, int end, Matrix &B, Matrix &x, Matrix &prev, int synchOffset)
    {
        for (int i = start; i < end; i++) 
        {
            double temp = getVal(i, i);
            getRef(i, i) = B.getVal(i, 0);
            B.getRef(i, 0) = temp;
            for (int j = 0; j < getCols(); j++) 
            {
                if (j != i)
                    getRef(i, j) *= -1;
        
                getRef(i, j) /= B.getVal(i, 0);
            }
            prev.getRef(i, 0) = 1.0;
            x.getRef(i, 0) = B.getVal(i, 0);
            mainSimpi->synch();
        }
        mainSimpi->synchExtraCycles(synchOffset);
    }

    void Matrix::jacobiFirstIteration(int start, int end, Matrix &x, Matrix &prev, int synchOffset)
    {
        for (int i = start; i < end; i++) 
        {
            double rowSum = 0;
            for (int j = 0; j < getCols(); j++) 
            {
                if (j == i)
                    rowSum += getVal(i, j);
                else
                    rowSum += (getVal(i, j) * prev.getVal(j, 0));
                
            }
            x.getRef(i, 0) = rowSum;
            mainSimpi->synch();
        }
        mainSimpi->synchExtraCycles(synchOffset);
    }

    void Matrix::jacobiRemainingIterations(int start, int end, Matrix &x, Matrix &prev, int synchOffset)
    {
        for (int k = 0; k < 1000; k++)
        {
            for (int i = start; i < end; i++)
            {
                //save prev value for comparision
                prev.getRef(i, 0) = x.getVal(i, 0);
                mainSimpi->synch();
            }
            mainSimpi->synchExtraCycles(synchOffset);
            for (int i = start; i < end; i++)
            {
                double rowSum = 0;
                for (int j = 0; j < getCols(); j++)
                {
                    if (j == i) {
                        rowSum += getVal(i, j);
                    } else {
                        rowSum += (getVal(i, j) * prev.getVal(j, 0));
                    }
                }
                x.getRef(i, 0) = rowSum;
                //mainSimpi->synch();
            }
            //wait at end of each loop for all processes before beginning next iteration
            mainSimpi->synch();
        }
    }

    /**
     * Restores saved initial inputs to.
     * 
     * @param start the first index for this process to work on
     * @param end the last index for this process to work on
     * @param saveEq A Matrix that the calling Matrix Object (A) as well as B was saved to
     * @param B a single column Matrix of constants that jacobi() solved for
     * @param saveB A single column Matrix that B was saved to
     */
    void Matrix::jacobiRestoreInputs(int start, int end, Matrix &saveEq, Matrix &B, Matrix &saveB)
    {
        for (int i = start; i < end; i++) 
        {
            for (int j = 0; j < getCols(); j++) 
                getRef(i, j) = saveEq.getVal(i,j);
            
            B.getRef(i, 0) = saveB.getVal(i, 0);
        }
    }

    /**
     * Solves a linear system of equations Ax=B in parallel if the Matrix is not diagonally dominant.
     * Uses the inverse-mutliplication method. 
     * Calling Matrix must have equal row and col count.
     * 
     * @param x a single column Matrix that the solution will be written into
     * @param B a single column Matrix of constants that each row is to be solved for
    */
    void Matrix::failSafe(Matrix &x, Matrix &B)
    {
        Matrix inv = inverse();
        mainSimpi->synch();

        int processCount = mainSimpi->getProcessCount();
        int id = mainSimpi->getID();
        
        int vectorSize = B.getRows();
        int work = vectorSize / processCount; 
        int start = 0;
        int end = 0;
        if (id < processCount) // Extra processes get no work => start and end stay at 0
        {
            start = id * work;
            end = start + work;
        }
        int remainingWork = B.getRows() % processCount;
        if (id == processCount - 1) // Last active process gets all remaining work
            end += remainingWork;
        
        double sol;
        for(int i = start; i < end; i++)
        {
            sol = 0;
            for(int j = 0; j < vectorSize; j++)
            {
                sol += (inv.getVal(i, j) * B.getVal(j, 0));
            }
            x.getRef(i, 0) = sol;
        }
        mainSimpi->synch();
        return;
    }

    /**
     * Checks if a Matrix is diagonally dominant.
     * A squre matrix is diagonally dominant if the diagonal terms 
     * are greater-than or equal to sum of the rest of their row.
     * Calling Matrix must have equal row and col count.
     */
    bool Matrix::isDiagonallyDominant()
    {
        if (!isSquareMatrix())
            return false;

        for(int i = 0; i < getRows(); i ++)
        {
            double sq;
            double rest = 0;
            for(int j = 0; j < getCols(); j++)
                (i == j) ? sq = getVal(i, j) : rest += getVal(i,j);
            
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

    /**
     * A == B
     * Checks if another Matrix has equal values to this Matrix.
     * Returns false if dimensions are different,
     * of if any of their values at corresponding (row, col) are unequal.
     * 
     * Due to Matrix values being doubles, 
     * it's possible that stray bits will cause undesired results.
     * This can be adjusted by setting Matrix::setEqualityPrecision() to a higher value.
     * 
     * @param B another Matrix to compare
     */
    bool Matrix::equals(Matrix &B)
    {
        Matrix *A = this;
        if (A->rows != B.rows || A->cols != B.cols)
            return false;
            
        int fd;
        std::string sharedMemoryName = "TEMP_SHARED_BOOL_FOR_EQUALITY_CHECK";
        bool *equalityBool = getSharedBool(fd, sharedMemoryName); // Shared between all processes so the same value is always returned
        *equalityBool = true;
        mainSimpi->synch();

        int cellCount = A->rows * A->cols;
        int indices[4];
        singleCellWorkDivision(cellCount, indices);
        determineEquality(B, indices[0], indices[1], equalityBool); // Initial work
        determineEquality(B, indices[2], indices[3], equalityBool); // Leftover work
        mainSimpi->synch();

        bool eqValue = *equalityBool;
        mainSimpi->synch();

        close(fd);
        shm_unlink(sharedMemoryName.c_str());
        munmap(equalityBool, sizeof(bool*));
        
        mainSimpi->synch();
        return eqValue;
    }

    /**
     * Determines if A == B
     * Compares all corresponding Matrix elements that have the same index.
     * If an unequal set of elements is found the function will return false
     * and force all other active processes to return false as well.
     * 
     * @param B the other Matrix being compared.
     * @param start the first double array index of A and B for this process to work on
     * @param end the last double array index of A and B for this process to work on
     * @param eqValue shared boolean value between all active processes
     */
    void Matrix::determineEquality(Matrix &B, int start, int end, bool* eqValue)
    {
        Matrix *A = this;
        for (int i = start; i < end; i++)
        {
            if (fabs(A->arr[i] - B.arr[i]) > equalityPrecision)
                *eqValue = false;                    
            if (!*eqValue)
                return;
        }
    }

    bool operator==(Matrix &lhs, Matrix &rhs)
    {
        return lhs.equals(rhs);
    }

    bool operator!=(Matrix &lhs, Matrix &rhs)
    {
        return !lhs.equals(rhs);
    }

    /**
     * Returns a pointer to a bool value shared between all active processes.
     * Two Matrix objects are not equal if they have different elements.
     * While processes work in parallel, 
     * they need a way to communicate if they've found an unequal set of elements.
     * 
     * @param fd file descriptor for shared boolean pointer
     */
    bool* Matrix::getSharedBool(int &fd, std::string sharedMemoryName)
    {
        bool *eq;
        if (mainSimpi->getID() == 0) 
        {
            fd = shm_open(sharedMemoryName.c_str(), O_RDWR | O_CREAT, 0777); 
            ftruncate(fd, sizeof(bool*));
            eq = (bool*)mmap(NULL, sizeof(bool*), PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
            mainSimpi->synch();
        }
        else 
        {
            mainSimpi->synch();
            fd = shm_open(sharedMemoryName.c_str(), O_RDWR, 0777);
            eq = (bool*)mmap(NULL, sizeof(bool*), PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
        }
        return eq;
    }

    /**
     * C = AB
     * Solves matrix multiplication in parallel and outputs the product solution.
     * 
     * @param B matrix that this Matrix (A) is being multiplied with.
     *          B's row count must match this Matrix (A)'s column count.
     * @return the matrix C, the product of this Matrix (A) and B.
     */
    Matrix &Matrix::multiply(Matrix &B)
    {
        Matrix *A = this;
        if (cols != B.rows)
        {
            printf("Cannot multiply matrices of sizes %dx%d and %dx%d", A->rows, A->cols, B.rows, B.cols); 
            exit(1);
        }
        Matrix *C = new Matrix(A->rows, B.cols);
        int cellCount = C->rows * C->cols;
        int indices[4];
        singleCellWorkDivision(cellCount, indices);
        calculateProduct(B, C, indices[0], indices[1]); // Initial work
        calculateProduct(B, C, indices[2], indices[3]); // Leftover work

        mainSimpi->synch();
        return *C;
    }

    /**
     * Calculates C = A * B 
     * 
     * @param B the Matrix being multiplied with this Matrix (A)
     * @param C the resulting Matrix that the solution is being written into
     * @param start the first double array index of C for this process to work on
     * @param end the last double array index of C for this process to work on
     */
    void Matrix::calculateProduct(Matrix &B, Matrix* C, int start, int end)
    {
        Matrix *A = this; // For clarity
        for (int i = start; i < end; i++)
        {
            C->arr[i] = 0;
            int row = C->getRow(i);
            int col = C->getCol(i);
            for (int offset = 0; offset < B.rows; offset++) // A->col and B.rows are the same
            {
                C->getRef(row, col) += A->getVal(row, offset) * B.getVal(offset, col);
            }
        }
    }

    Matrix &operator*(Matrix &lhs, Matrix &rhs)
    {
        return lhs.multiply(rhs);
    }

    void operator*=(Matrix &lhs, Matrix &rhs)
    {
        lhs = lhs.multiply(rhs);
    }

    /**
     * ATimesLambda = A * lambda
     * Solves scalar multiplication in parallel and outputs the product solution.
     * 
     * @param lambda scalar value that this Matrix (A) is being multiplied with.
     * @return the Matrix ATimesLambda, the product of this Matrix (A) and lambda.
     */
    Matrix &Matrix::scalarMultiply(double lambda)
    {
        Matrix *A = this;
        Matrix *ATimesLambda = new Matrix(rows, cols);

        int cellCount = ATimesLambda->rows * ATimesLambda->cols;
        int indices[4];
        singleCellWorkDivision(cellCount, indices);

        for (int i = indices[0]; i < indices[1]; i++) // Initial work
            ATimesLambda->arr[i] = A->arr[i] * lambda;

        for (int i = indices[2]; i < indices[3]; i++) // Leftover work
            ATimesLambda->arr[i] = A->arr[i] * lambda;

        mainSimpi->synch();
        return *ATimesLambda;
    }

    Matrix &operator*(double lhs, Matrix &rhs)
    {
        return rhs.scalarMultiply(lhs);
    }

    Matrix &operator*(Matrix &lhs, double rhs)
    {
        return lhs.scalarMultiply(rhs);
    }

    void operator*=(Matrix &lhs, double rhs)
    {
        lhs = lhs.scalarMultiply(rhs);
    }

    /**
     * C = A + B
     * Solves matrix addition in parallel and outputs the product solution.
     * 
     * @param B a Matrix of equal size to this Matrix (A) 
     * @return the Matrix C, the sum of this Matrix (A) and B 
     */
    Matrix &Matrix::add(Matrix &B)
    {
        Matrix *A = this;
        Matrix *C = new Matrix(rows, cols);

        int cellCount = C->rows * C->cols;
        int indices[4];
        singleCellWorkDivision(cellCount, indices);

        for (int i = indices[0]; i < indices[1]; i++) // Initial work
            C->arr[i] = A->arr[i] + B.arr[i];

        for (int i = indices[2]; i < indices[3]; i++) // Leftover work
            C->arr[i] = A->arr[i] + B.arr[i];

        mainSimpi->synch();
        return *C;
    }

    Matrix &operator+(Matrix &lhs, Matrix &rhs)
    {
        return rhs.add(lhs);
    }

    void operator+=(Matrix &lhs, Matrix &rhs)
    {
        lhs = lhs.add(rhs);
    }

    /**
     * C = A - B
     * Solves matrix subtraction in parallel and outputs the product solution.
     * 
     * @param B a Matrix of equal size to this Matrix (A) 
     * @return the Matrix C, the difference between this Matrix (A) and B 
     */
    Matrix &Matrix::subtract(Matrix &B)
    {
        Matrix *A = this;
        Matrix *C = new Matrix(rows, cols);

        int cellCount = C->rows * C->cols;
        int indices[4];
        singleCellWorkDivision(cellCount, indices);

        for (int i = indices[0]; i < indices[1]; i++) // Initial work
            C->arr[i] = A->arr[i] - B.arr[i];

        for (int i = indices[2]; i < indices[3]; i++) // Leftover work
            C->arr[i] = A->arr[i] - B.arr[i];
            
        mainSimpi->synch();
        return *C;
    }

    Matrix &operator-(Matrix &lhs, Matrix &rhs)
    {
        return lhs.subtract(rhs);
    }

    void operator-=(Matrix &lhs, Matrix &rhs)
    {
        lhs = lhs.subtract(rhs);
    }

    /**
     * A^T 
     * Returns the transpose of this Matrix (A).
     * The tranpose of a matrix is one in which its row and col indices are switched,
     * such that (row, col) -> (col, row).
     * Its dimensions are also switched from (row size, col size) to (col size, row size). 
     */
    Matrix &Matrix::transpose()
    {
        Matrix *A = this;
        Matrix *A_T = new Matrix(A->cols, A->rows);

        int cellCount = A->rows * A->cols;
        int indices[4];
        singleCellWorkDivision(cellCount, indices);
        calculateTranspose(A_T, indices[0], indices[1]); // Initial work
        calculateTranspose(A_T, indices[2], indices[3]); // Leftover work

        mainSimpi->synch();
        return *A_T;
    }

    /**
     * Calculates A^T
     * 
     * @param A_T the resulting Matrix that the solution is being written into
     * @param start the first double array index of A_T for this process to work on
     * @param end the last double array index of A_T for this process to work on
     */
    void Matrix::calculateTranspose(Matrix* A_T, int start, int end)
    {
        Matrix *A = this;
        for (int i = start; i < end; i++)
        {
            int rowA = A->getRow(i);
            int colA = A->getCol(i);
            A_T->getRef(colA, rowA) = A->getVal(rowA, colA);
        }
    }

    std::ostream& operator<<(std::ostream& out, const Matrix& m)
    {
        if (m.getSimpiID() == 0)
        {
            for (int i = 0; i < m.rows; i++)
            {
                out << std::endl;
                for (int j = 0; j < m.cols; j++)
                {
                    out << std::fixed << std::setprecision(Matrix::printPrecision) << m.arr[i + j * m.rows];
                    out << ", ";
                }
            }
            out << std::endl;
        }
        return out;
    }
}
