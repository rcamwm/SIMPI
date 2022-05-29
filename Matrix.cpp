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

    Matrix::Matrix(int rowCount, int colCount)
    {
        // use mainSimpi to init the Matrix for all processes. The id is also in simp
        std::pair<std::string, double*> passBack(mainSimpi->createMatrix(rowCount, colCount));
        uniqueID = passBack.first;
        arr = passBack.second;
        rows = rowCount;
        cols = colCount;
    }

    Matrix::Matrix(const Matrix &m)
    {
        std::pair<std::string, double*> passBack(mainSimpi->createMatrix(m.rows, m.cols));
        uniqueID = passBack.first;
        arr = passBack.second;
        rows = m.rows;
        cols = m.cols;

        // Processes divide the rows if there are more rows, and divide the columns if not
        bool moreRows = rows > cols;
        int div = (moreRows) ? rows : cols; // Number of lines to divide between processes

        int processCount = mainSimpi->getProcessCount();
        int processID = mainSimpi->getID();

        if (div <= processCount)
        {
            int start = processID;
            int end = start + 1;
            if (processID < div)
                copyElements(m, start, end, moreRows);
        }
        else 
        {
            int work = div / processCount;
            int start = processID * work;
            int end = start + work;
            copyElements(m, start, end, moreRows);

            int leftoverWork = div % processCount;
            if (leftoverWork != 0)
            {
                start = (work * processCount) + processID;
                end = start + 1;
                if (processID < leftoverWork)
                    copyElements(m, start, end, moreRows);
            }         
        }
        mainSimpi->synch();
    }

    Matrix::~Matrix()
    {
        mainSimpi->freeMatrix(uniqueID); // frees and unlinks memory
    }

    void Matrix::copyElements(const Matrix &m, int start, int end, bool moreRows)
    {
        int rowStart, rowEnd, colStart, colEnd;
        if (moreRows) { rowStart = start, rowEnd = end, colStart = 0, colEnd = cols; }
        else { rowStart = 0, rowEnd = rows, colStart = start, colEnd = end; }

        for (int row = rowStart; row < rowEnd; row++)
        {
            for (int col = colStart; col < colEnd; col++)
            {
                getRef(row, col) = m.getVal(row, col);
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

    void Matrix::adjoint(Matrix &adj)
    {
        if (!isSquareMatrix() || !adj.isSquareMatrix() || getRows() != adj.getRows()) 
        {
            std::cout << "Invalid Matrices: Must be equal sized square matrices" << std::endl;
            exit(1);
        }
        allocateAdjointWork(this->arr, adj.arr, getRows());
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
                    sum += (lower->getVal(i, j) * upper->getVal(j, k)); 

                // Evaluating U(i, k) 
                upper->getRef(i,k) = getVal(i,k) - sum; 

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
                            sum += (lower->getVal(i, j) * upper->getVal(j, a));

                        // Evaluating U(i, k) 
                        upper->getRef(i, a) = getVal(i, a) - sum;
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
                    lower->getRef(i, i) = 1; // Diagonal as 1 
                else 
                {
                    // Summation of L(k, j) * U(j, i)
                    float sum = 0;
                    for (int j = 0; j < i; j++) 
                    sum += (lower->getVal(k, j) * upper->getVal(j, i));
                    // Evaluating L(k, i)
                    lower->getRef(k, i) = ((getVal(k, i) - sum) / upper->getVal(i, i));
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
                            lower->getRef(i, i) = 1; // Diagonal as 1
                        else 
                        {
                            // Summation of L(k, j) * U(j, i) 
                            float sum = 0;
                            for (int j = 0; j < i; j++) 
                            sum += (lower->getVal(a, j) * upper->getVal(j, i));

                            // Evaluating L(k, i) 
                            lower->getRef(a, i) = (getVal(a, i)-sum) / upper->getVal(i, i);
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
    void Matrix::inverse(Matrix &inv) 
    {
        if (!isSquareMatrix() || !inv.isSquareMatrix() || getRows() != inv.getRows()) 
        {
            std::cout << "Invalid Matrices: Must be equal sized square matrices" << std::endl;
            exit(1);
        }

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
                inv.getRef(c,a) = solutionCol[c];
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
                        inv.getRef(c, a) = solutionCol[c];
                    
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

    /*
    Ax=B
    Solves a linear system of equations
    If the Matrix is diagonally dominant, jacobi-iterative method is used
    else, the inverse mutliplication method is used
    takes in a Vector of constants that each row is to be solved for and a solution Vector of 0s in which the solution will be written in
    void -> void
    */
    void Matrix::solveSystem(Matrix &solution, Matrix &constants)
    {
        if (solution.cols != 1 || solution.rows != rows || constants.cols != 1 || constants.rows != rows)
        {
            std::cout << "Cannot solve linear system of equations for x(" << solution.rows << ", " << solution.cols << ") ";
            std::cout << "and B(" << constants.rows << ", " << constants.cols << ")" << std::endl;
            std::cout << "Parameters must be size x(1, " << rows << ") and B(1, " << rows << std::endl;
            exit(1);
        }
        bool dd = isDiagonallyDominant(); 
        mainSimpi->synch();

        if (dd)
        {
            mainSimpi->synch();
            jacobi(solution, constants);
        }
        else
        {
            mainSimpi->synch();
            failSafe(solution, constants);
        }
        mainSimpi->synch();
    }

    /*
        Solves a linear system of equations in parallel if the Matrix is diagonally dominant
        outputs solution to a Vector of 0s passed into it
        void->void
    */
    void Matrix::jacobi(Matrix &solution, Matrix &constants) 
    {
        int id = mainSimpi->getID();
        int processCount = mainSimpi->getProcessCount();
        if (processCount > constants.getRows()) 
            processCount = constants.getRows();
        
        Matrix saveEq(getRows(), getCols() + 1); // save equations from modification
        Matrix saveConst(constants.getRows(), 1); // saves input Matrix
        Matrix prev(constants.getRows(), 1); // shared mem containing a copy of values

        int work = constants.getRows() / processCount; 
        int start = 0;
        int end = 0;
        if (id < processCount) // Extra processes get no work => start and end stay at 0
        {
            start = id * work;
            end = start + work;
        }
        int remainingWork = constants.getRows() % processCount;
        if (id == processCount - 1) // Last active process gets all remaining work
            end += remainingWork;
        // synch() is called (end - start) number of times in several jacobi helper functions
        // synchOffset makes sure that all processes stay in sync even if some processes are doing more work than others
        int synchOffset = (work + remainingWork) - end + start;
        mainSimpi->synch();

        jacobiSaveInputs(start, end, saveEq, constants, saveConst);
        mainSimpi->synch();
        std::cout << "ID " << id << ": saved" << std::endl;

        // Switches diagonals with solution elements,
        // then divides each matrix row by their original diagonal element 
        jacobiSwitchAndDivide(start, end, constants, solution, prev, synchOffset);
        mainSimpi->synch();
        std::cout << "ID " << id << ": switched and divided" << std::endl;

        // First iteration with substitution 1
        jacobiFirstIteration(start, end, solution, prev, synchOffset);
        mainSimpi->synch();
        std::cout << "ID " << id << ": iterated first" << std::endl;

        // Repeat iterations with calculated results
        jacobiRemainingIterations(start, end, solution, prev, synchOffset);
        mainSimpi->synch();
        std::cout << "ID " << id << ": iterated remaining" << std::endl;

        jacobiRestoreInputs(start, end, saveEq, constants, saveConst);
        mainSimpi->synch();
        std::cout << "ID " << id << ": restored" << std::endl;


        return;
    }

    void Matrix::jacobiSaveInputs(int start, int end, Matrix &saveEq, Matrix &constants, Matrix &saveConst)
    {
        for (int i = start; i < end; i++) 
        {
            for (int j = 0; j < getCols(); j++) 
                saveEq.getRef(i,j) = getVal(i, j);

            saveEq.getRef(i,getCols() + 1) = constants.getVal(i, 0);
            saveConst.getRef(i, 0) = constants.getVal(i, 0);
        }
    }

    // /*
    //     In each row, the diagonal element of the matrix is switched with the corresponding element from the solution vector.
    //     Every element of that matrix row is then divided by the value that was just placed in the solution vector.
    //     Each element of that matrix row that was not switched is also multiplied by -1.
    // */
    void Matrix::jacobiSwitchAndDivide(int start, int end, Matrix &constants, Matrix &solution, Matrix &prev, int synchOffset)
    {
        for (int i = start; i < end; i++) 
        {
            double temp = getVal(i, i);
            getRef(i, i) = constants.getVal(i, 0);
            constants.getRef(i, 0) = temp;
            for (int j = 0; j < getCols(); j++) 
            {
                if (j != i)
                    getRef(i, j) *= -1;
        
                getRef(i, j) /= constants.getVal(i, 0);
                //mainSimpi->synch();
            }
            prev.getRef(i, 0) = 1.0;
            solution.getRef(i, 0) = constants.getVal(i, 0);
            mainSimpi->synch();
        }
        mainSimpi->synchExtraCycles(synchOffset);
    }

    void Matrix::jacobiFirstIteration(int start, int end, Matrix &solution, Matrix &prev, int synchOffset)
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
                
                //mainSimpi->synch();
            }
            solution.getRef(i, 0) = rowSum;
            mainSimpi->synch();
        }
        mainSimpi->synchExtraCycles(synchOffset);
    }

    void Matrix::jacobiRemainingIterations(int start, int end, Matrix &solution, Matrix &prev, int synchOffset)
    {
        for (int k = 0; k < 1000; k++)
        {
            for (int i = start; i < end; i++)
            {
                //save prev value for comparision
                prev.getRef(i, 0) = solution.getVal(i, 0);
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
                solution.getRef(i, 0) = rowSum;
                //mainSimpi->synch();
            }
            //wait at end of each loop for all processes before beginning next iteration
            mainSimpi->synch();
        }
    }

    void Matrix::jacobiRestoreInputs(int start, int end, Matrix &saveEq, Matrix &constants, Matrix &saveConst)
    {
        for (int i = start; i < end; i++) 
        {
            for (int j = 0; j < getCols(); j++) 
                getRef(i, j) = saveEq.getVal(i,j);
            
            constants.getRef(i, 0) = saveConst.getVal(i, 0);
        }
    }

    /*
    Method to solve a system of linear equations if the system is not diagonally dominant
    Uses the inverse-mutliplication method. 
    void->void
    */
    void Matrix::failSafe(Matrix &solution, Matrix &constants)
    {
        Matrix* inv = new Matrix(getRows(), getCols());
        inverse(*inv);
        mainSimpi->synch();

        int processCount = mainSimpi->getProcessCount();
        int id = mainSimpi->getID();
        
        int vectorSize = constants.getRows();
        int work = vectorSize / processCount; 
        int start = 0;
        int end = 0;
        if (id < processCount) // Extra processes get no work => start and end stay at 0
        {
            start = id * work;
            end = start + work;
        }
        int remainingWork = constants.getRows() % processCount;
        if (id == processCount - 1) // Last active process gets all remaining work
            end += remainingWork;
        
        double sol;
        for(int i = start; i < end; i++)
        {
            sol = 0;
            for(int j = 0; j < vectorSize; j++)
            {
                sol += (inv->getVal(i, j) * constants.getVal(j, 0));
            }
            solution.getRef(i, 0) = sol;
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

    bool Matrix::equals(Matrix &comparand)
    {
        if (rows != comparand.rows || cols != comparand.cols)
            return false;
            
        int fd;
        bool *equalityBool = getSharedBool(fd); // Shared between all processes so the same value is always returned
        *equalityBool = true;
        mainSimpi->synch();

        // Processes divide the rows if matrices have more rows, and divide the columns if not
        bool moreRows = rows > cols;
        int div = (moreRows) ? rows : cols; // Number of lines to divide between processes

        int processCount = mainSimpi->getProcessCount();
        int processID = mainSimpi->getID();

        if (div <= processCount)
        {
            int start = processID;
            int end = start + 1;
            if (processID < div)
                determineEquality(comparand, start, end, moreRows, equalityBool);
        }
        else 
        {
            int work = div / processCount;
            int start = processID * work;
            int end = start + work;
            determineEquality(comparand, start, end, moreRows, equalityBool);

            int leftoverWork = div % processCount;
            if (leftoverWork != 0)
            {
                start = (work * processCount) + processID;
                end = start + 1;
                if (processID < leftoverWork)
                    determineEquality(comparand, start, end, moreRows, equalityBool);
            }         
        }
        mainSimpi->synch();
        bool eqValue = *equalityBool;
        
        mainSimpi->synch();
        close(fd);
        shm_unlink("TEMP_SHARED_BOOL_FOR_EQUALITY_CHECK");
        munmap(equalityBool, sizeof(bool*));
        
        mainSimpi->synch();
        return eqValue;
    }

    bool operator==(Matrix &lhs, Matrix &rhs)
    {
        return lhs.equals(rhs);
    }

    bool operator!=(Matrix &lhs, Matrix &rhs)
    {
        return !lhs.equals(rhs);
    }

    bool* Matrix::getSharedBool(int &fd)
    {
        bool *eq;
        if (mainSimpi->getID() == 0) 
        {
            fd = shm_open("TEMP_SHARED_BOOL_FOR_EQUALITY_CHECK", O_RDWR | O_CREAT, 0777); 
            ftruncate(fd, sizeof(bool*));
            eq = (bool*)mmap(NULL, sizeof(bool*), PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
            mainSimpi->synch();
        }
        else 
        {
            mainSimpi->synch();
            fd = shm_open("TEMP_SHARED_BOOL_FOR_EQUALITY_CHECK", O_RDWR, 0777);
            eq = (bool*)mmap(NULL, sizeof(bool*), PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
        }
        return eq;
    }

    void Matrix::determineEquality(Matrix &comparand, int start, int end, bool moreRows, bool* eqValue)
    {
        int rowStart, rowEnd, colStart, colEnd;
        if (moreRows) { rowStart = start, rowEnd = end, colStart = 0, colEnd = cols; }
        else { rowStart = 0, rowEnd = rows, colStart = start, colEnd = end; }

        for (int row = rowStart; row < rowEnd; row++)
        {
            for (int col = colStart; col < colEnd; col++)
            {
                if (fabs(this->getVal(row, col) - comparand.getVal(row, col)) > equalityPrecision)
                    *eqValue = false;                    
                if (!*eqValue)
                    return;                
            }
        }
    }
    
    /**
     * C = AB
     * Solves matrix multiplication in parallel and outputs the product solution.
     * 
     * @param B matrix that this Matrix (A) is being being multiplied with.
     *          B's row count must match this Matrix (A)'s column count.
     * @return the matrix C, the product of this Matrix (A) and B.
     */
    Matrix &Matrix::multiply(Matrix &B)
    {
        if (cols != B.rows)
        {
            printf("Cannot multiply matrices of sizes %dx%d and %dx%d", rows, cols, B.rows, B.cols); 
            exit(1);
        }
        Matrix *C = new Matrix(rows, B.cols);

        // Processes divide the rows if C has more rows, and divide the columns if not
        bool moreRows = C->rows > C->cols;
        int div = (moreRows) ? C->rows : C->cols; // Number of lines to divide between processes

        int processCount = mainSimpi->getProcessCount();
        int processID = mainSimpi->getID();

        if (div <= processCount)
        {
            int start = processID;
            int end = start + 1;
            if (processID < div)
                calculateProduct(B, C, start, end, moreRows);
        }
        else 
        {
            int work = div / processCount;
            int start = processID * work;
            int end = start + work;
            calculateProduct(B, C, start, end, moreRows);

            int leftoverWork = div % processCount;
            if (leftoverWork != 0)
            {
                start = (work * processCount) + processID;
                end = start + 1;
                if (processID < leftoverWork)
                    calculateProduct(B, C, start, end, moreRows);
            }         
        }
        mainSimpi->synch();
        return *C;
    }

    Matrix &operator*(Matrix &lhs, Matrix &rhs)
    {
        return lhs.multiply(rhs);
    }

    void operator*=(Matrix &lhs, Matrix &rhs)
    {
        lhs = lhs.multiply(rhs);
    }

    void Matrix::calculateProduct(Matrix &B, Matrix* C, int start, int end, bool moreRows)
    {
        Matrix *A = this; // For clarity

        // Processes divide the rows if C has more rows, and divide the columns if not
        int rowStart, rowEnd, colStart, colEnd;
        if (moreRows) { rowStart = start, rowEnd = end, colStart = 0, colEnd = C->cols; }
        else { rowStart = 0, rowEnd = C->rows, colStart = start, colEnd = end; }

        for (int row = rowStart; row < rowEnd; row++)
        {
            for (int col = colStart; col < colEnd; col++)
            {
                C->getRef(row, col) = 0;
                for (int offset = 0; offset < B.rows; offset++) // A->col == B.rows
                {
                    C->getRef(row, col) += A->getVal(row, offset) * B.getVal(offset, col);
                }
            }
        }
    }

    Matrix &Matrix::scalarMultiply(double operand)
    {
        Matrix *product = new Matrix(rows, cols);

        // Processes divide the rows if there are more rows, and divide the columns if not
        bool moreRows = rows > cols;
        int div = (moreRows) ? rows : cols; // Number of lines to divide between processes

        int processCount = mainSimpi->getProcessCount();
        int processID = mainSimpi->getID();

        if (div <= processCount)
        {
            int start = processID;
            int end = start + 1;
            if (processID < div)
                calculateScalarProduct(operand, product, start, end, moreRows);
        }
        else 
        {
            int work = div / processCount;
            int start = processID * work;
            int end = start + work;
            calculateScalarProduct(operand, product, start, end, moreRows);

            int leftoverWork = div % processCount;
            if (leftoverWork != 0)
            {
                start = (work * processCount) + processID;
                end = start + 1;
                if (processID < leftoverWork)
                    calculateScalarProduct(operand, product, start, end, moreRows);
            }         
        }
        mainSimpi->synch();
        return *product;
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

    void Matrix::calculateScalarProduct(double lambda, Matrix* product, int start, int end, bool moreRows)
    {
        // Processes divide the rows if there are more rows, and divide the columns if not
        int rowStart, rowEnd, colStart, colEnd;
        if (moreRows) { rowStart = start, rowEnd = end, colStart = 0, colEnd = cols; }
        else { rowStart = 0, rowEnd = rows, colStart = start, colEnd = end; }

        for (int row = rowStart; row < rowEnd; row++)
        {
            for (int col = colStart; col < colEnd; col++)
            {
                product->getRef(row, col) = getVal(row, col) * lambda;
            }
        }
    }

    Matrix &Matrix::add(Matrix &operand)
    {
        Matrix *sum = new Matrix(rows, cols);

        // Processes divide the rows if there are more rows, and divide the columns if not
        bool moreRows = rows > cols;
        int div = (moreRows) ? rows : cols; // Number of lines to divide between processes

        int processCount = mainSimpi->getProcessCount();
        int processID = mainSimpi->getID();

        if (div <= processCount)
        {
            int start = processID;
            int end = start + 1;
            if (processID < div)
                calculateSum(operand, sum, start, end, moreRows);
        }
        else 
        {
            int work = div / processCount;
            int start = processID * work;
            int end = start + work;
            calculateSum(operand, sum, start, end, moreRows);

            int leftoverWork = div % processCount;
            if (leftoverWork != 0)
            {
                start = (work * processCount) + processID;
                end = start + 1;
                if (processID < leftoverWork)
                    calculateSum(operand, sum, start, end, moreRows);
            }         
        }
        mainSimpi->synch();
        return *sum;
    }

    Matrix &operator+(Matrix &lhs, Matrix &rhs)
    {
        return rhs.add(lhs);
    }

    void operator+=(Matrix &lhs, Matrix &rhs)
    {
        lhs = lhs.add(rhs);
    }

    void Matrix::calculateSum(const Matrix &B, Matrix* sum, int start, int end, bool moreRows)
    {
        // Processes divide the rows if there are more rows, and divide the columns if not
        int rowStart, rowEnd, colStart, colEnd;
        if (moreRows) { rowStart = start, rowEnd = end, colStart = 0, colEnd = cols; }
        else { rowStart = 0, rowEnd = rows, colStart = start, colEnd = end; }

        for (int row = rowStart; row < rowEnd; row++)
        {
            for (int col = colStart; col < colEnd; col++)
            {
                sum->getRef(row, col) = getVal(row, col) + B.getVal(row, col);
            }
        }
    }

    Matrix &Matrix::subtract(Matrix &operand)
    {
        Matrix *difference = new Matrix(rows, cols);

        // Processes divide the rows if there are more rows, and divide the columns if not
        bool moreRows = rows > cols;
        int div = (moreRows) ? rows : cols; // Number of lines to divide between processes

        int processCount = mainSimpi->getProcessCount();
        int processID = mainSimpi->getID();

        if (div <= processCount)
        {
            int start = processID;
            int end = start + 1;
            if (processID < div)
                calculateDifference(operand, difference, start, end, moreRows);
        }
        else 
        {
            int work = div / processCount;
            int start = processID * work;
            int end = start + work;
            calculateDifference(operand, difference, start, end, moreRows);

            int leftoverWork = div % processCount;
            if (leftoverWork != 0)
            {
                start = (work * processCount) + processID;
                end = start + 1;
                if (processID < leftoverWork)
                    calculateDifference(operand, difference, start, end, moreRows);
            }         
        }
        mainSimpi->synch();
        return *difference;
    }

    Matrix &operator-(Matrix &lhs, Matrix &rhs)
    {
        return lhs.subtract(rhs);
    }

    void operator-=(Matrix &lhs, Matrix &rhs)
    {
        lhs = lhs.subtract(rhs);
    }

    void Matrix::calculateDifference(const Matrix &B, Matrix* difference, int start, int end, bool moreRows)
    {
        // Processes divide the rows if there are more rows, and divide the columns if not
        int rowStart, rowEnd, colStart, colEnd;
        if (moreRows) { rowStart = start, rowEnd = end, colStart = 0, colEnd = cols; }
        else { rowStart = 0, rowEnd = rows, colStart = start, colEnd = end; }

        for (int row = rowStart; row < rowEnd; row++)
        {
            for (int col = colStart; col < colEnd; col++)
            {
                difference->getRef(row, col) = getVal(row, col) - B.getVal(row, col);
            }
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
