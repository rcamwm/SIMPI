#include "Matrix.h"

// Must be set before any Matrix objects are created
Simpi* Matrix::main_simpi;

Matrix::Matrix(int x, int y)  // constructor
{
  // use main_simp init the Matrix for all processes. The id is also in simp
  std::pair<std::string, double*> pass_back(main_simpi->create_matrix(x, y));
  unique_id = pass_back.first;
  arr = pass_back.second;
  xdim = x;
  ydim = y;
}
Matrix::~Matrix()  // destructor
{
  // use main_simpi for getting rid of the mem and unlink stuff
  main_simpi->free_matrix(unique_id);
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
        return out;
    }
    else
    {
        return out;
    }
}

/*
This is a helper function to calculate the determinant of a Matrix
*/
int Matrix::determinant(double* A, int n, int order)
{
  int D = 0;  // Initialize result

  //  Base case : if Matrix contains single element
  if (n == 1)
    return A[0];

  double temp[order * order];  // To store cofactors

  int sign = 1;  // To store sign multiplier

  // Iterate for each element of first row
  for (int f = 0; f < n; f++) {
    // Getting Cofactor of A[0][f]
    getCofactor(A, temp, 0, f, n, order);
    D += sign * A[0 + f * order] * determinant(temp, n - 1, order);

    // terms are to be added with alternate sign
    sign = -sign;
  }

  return D;
}

/*
This is a helper function to calculate the adjoint of a Matrix
*/
void Matrix::adjoint(double* A, double* adj, int order, int par_id, int par_count)
{
  if (order == 1) {
    adj[0] = 1;
    return;
  }

  int rpp = order / par_count;
  int start = par_id * rpp;
  int end = start + rpp;

  // temp is used to store cofactors of A[][]
  int sign = 1;
  double temp[order * order];

  for (int i = 0; i < order; i++) {
    for (int j = start; j < end; j++) {
      // Get cofactor of A[i][j]
      getCofactor(A, temp, i, j, order, order);

      // sign of adj[j][i] positive if sum of row
      // and column indexes is even.
      sign = ((i + j) % 2 == 0) ? 1 : -1;

      // Interchanging rows and columns to get the
      // transpose of the cofactor Matrix
      adj[j + i * order] = (sign) * (determinant(temp, order - 1, order));
    }
  }
}

/*
This is a helper function to calculate the cofactor of a Matrix
*/
void Matrix::getCofactor(double* A, double* temp, int p, int q, int n, int order)
{
  int i = 0, j = 0;

  // Looping for each element of the Matrix
  for (int row = 0; row < n; row++) {
    for (int col = 0; col < n; col++) {
      //  Copying into temporary Matrix only those element
      //  which are not in given row and column
      if (row != p && col != q) {
        temp[(i) + (j++) * order] = A[row + col * order];

        // Row is filled, so increase row index and
        // reset col index
        if (j == n - 1) {
          j = 0;
          i++;
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

  // Check if Matrix is square
  if (get_x() != get_y()) {
    std::cout << "Invalid Matrix";
    exit(1);
  }

  for (int i = 0; i < get_x(); i++) {
    // Calculate work per parallel process
    // Has to be calculated on every loop iteration as the inner loop is decrementing
    int num_processes = main_simpi->getProcessCount();
    int parID = main_simpi->getID();
    int total = xdim - i;
    if (num_processes > total) {
      num_processes = total;
    }
    int rpp = total/num_processes;
    int start = rpp*main_simpi->getID() + i;
    int end = start + rpp;

    // Upper Triangular 
    for (int k = start; k< end; k++) {
      if (k>=get_x()) {
        break;
      }
      // Summation of L(i, j) * U(j, k) 
      float sum = 0;
      for (int j = 0; j < i; j++) 
        sum += (lower->get(i,j)* upper->get(j,k)); 

      // Evaluating U(i, k) 
      upper->get(i,k) = get(i,k) - sum; 

    }

    // Calculate and execute which processes take the leftover work 
    if (total%num_processes != 0) {
      int leftover = total%num_processes;
      if (parID < leftover) {
        parID += (xdim-leftover);
        int start = parID;
        int end = start + 1;
        for (int a = start; a<end; a++) {
          // Summation of L(i, j) * U(j, k) 
          float sum = 0;
          for (int j= 0; j<i; j++) 
            sum+= (lower->get(i,j)*upper->get(j,a));
          // Evaluating U(i, k) 
          upper->get(i,a) = get(i,a) - sum;
        }
      }
    }

    main_simpi->synch();

    total = get_x() - i;
    parID = main_simpi->getID();
    num_processes = main_simpi->getProcessCount();

    // Lower Triangular
    for (int k = start; k<end; k++) {
      if (k>=get_x()) {
        break;
      }
      if (i == k) {
        lower->get(i,i) = 1; // Diagonal as 1 
      }
      else {
        // Summation of L(k, j) * U(j, i)
        float sum = 0;
        for (int j = 0; j<i; j++) 
          sum += (lower->get(k,j) * upper->get(j,i));
        // Evaluating L(k, i)
        lower->get(k,i) = ((get(k,i)-sum) / upper->get(i,i));
      }
    }

    // Calculate and execute which processes take the leftover work 
    if (total%num_processes != 0) {
      int leftover = total%num_processes;
      if (parID < leftover) {
        parID += (get_x()-leftover);
        int start = parID;
        int end = start + 1;
        for (int a = start; a<end; a++) {
          if (i == a)
            lower->get(i,i) = 1; // Diagonal as 1
          else {
            // Summation of L(k, j) * U(j, i) 
            float sum = 0;
            for (int j = 0; j<i; j++) 
              sum+= (lower->get(a,j) * upper->get(j,i));

            // Evaluating L(k, i) 
            lower->get(a,i) = (get(a,i)-sum) / upper->get(i,i);
          }
        }
      }
    }
    main_simpi->synch();

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
void Matrix::inverse(Matrix* inv) {

  //Check if Matrix is square
  if (get_x() != get_y()) {
    std::cout << "Invalid Matrix";
    exit(1);
  }


  //Solve for lower and upper matrices
  Matrix* upper = new Matrix(get_x(), get_y());
  Matrix* lower = new Matrix(get_x(), get_y());
  luDecomposition(lower,upper);
  main_simpi->synch();

  //Create Identity nxn Matrix
  Matrix* identity = new Matrix(get_x(), get_y());
  if (main_simpi->getID() == 0) {
    for (int i = 0; i<get_x(); i++) {
      for (int j = 0; j<get_x(); j++) {
        if (i == j) {
          identity->get(i,j) = 1;
        }
        else {
          identity->get(i,j) = 0;                
        }
      }
    }
  }

  main_simpi->synch();

  // Calculate columns per parallel process
  int processCount = main_simpi->getProcessCount();
  if (processCount > get_x()) {
    processCount = get_x();
  }
  int cpp = get_x()/processCount;
  int start = cpp*main_simpi->getID();
  int end = start + cpp;

  // Initialize necessary arrays for future calculations
  // Each array is local to its own process
  float identity_col[get_x()];
  float z_col[get_x()];
  float soln_col[get_x()];

  for (int a = start; a<end; a++) {

    //Get individual columns of identity Matrix
    for (int b = 0; b<get_x(); b++) {
      identity_col[b] = identity->get(b,a);
    }

    //Reset Z column to solve for again
    for (int d = 0; d<get_x(); d++) {
      z_col[d] = 0;
    }
  
    //Solve LZ = I
    (*lower).forward_substitution(identity_col, z_col);

    //Reset X column to solve for again
    for (int d = 0; d<get_x(); d++) {
      soln_col[d] = 0;
    }

    //Solve UX = Z
    (*upper).backward_substitution(z_col, soln_col);

    //Input X column to corresponding columnn in final inverse Matrix
    for (int c = 0; c<get_x(); c++) {
      inv->get(c,a) = soln_col[c];
    }
  }

  // Calculate and execute which processes take the leftover rows 
  // ex. with 3 processes and a 10x10 Matrix:
  // 0-3 is process 0
  // 3-6 is process 1
  // 6-9 is process 2
  // 9 is the leftover column that is taken by process 0
  int parID = main_simpi->getID();
  if (get_x()% processCount != 0) {
    int leftover = get_x() % processCount;
    if (parID < leftover) {
      parID += (get_x()-leftover);
      int start = parID;
      int end = start + 1;
      for (int a = start; a<end; a++) {
        //Get individual columns of identity Matrix
        for (int b = 0; b<get_x(); b++) {
          identity_col[b] = identity->get(b,a);
        }

        //Reset Z column to solve for again
        for (int d = 0; d<get_x(); d++) {
          z_col[d] = 0;
        }
      
        //Solve LZ = I
        (*lower).forward_substitution(identity_col, z_col);

        //Reset X column to solve for again
        for (int d = 0; d<get_x(); d++) {
          soln_col[d] = 0;
        }

        //Solve UX = Z
        (*upper).backward_substitution(z_col, soln_col);

        //Input X column to corresponding columnn in final inverse Matrix
        for (int c = 0; c<get_x(); c++) {
          inv->get(c,a) = soln_col[c];
        }
      }
    }
  }  

  main_simpi->synch();
  return;
}

/*
This is a helper function to calculate the solutions of a lower triangular Matrix
*/
void Matrix::forward_substitution(float *b, float* x)
{
    double suma;
    for(int i=0; i < get_x(); i=i+1)
    {
        suma = 0;
        for(int j=0;j<i;j=j+1)
            suma = suma+ get(i,j) * x[j];

        x[i] = (b[i]-suma)/get(i,i);
    }
}

/*
This is a helper function to calculate the solutions of an upper triangular Matrix
*/
void Matrix::backward_substitution(float* b, float* x)
{
    double suma;
    for(int i= get_x()-1; i>=0; i=i-1)
    {
        suma=0;
        for(int j = get_x() - 1; j > i; j =j-1)
            suma = suma + get(i,j) * x[j];
 
        x[i]=(b[i] - suma)/get(i,i);
    }
}

/*
Solves a linear system of equations in parallel if the Matrix is diagonally dominant
outputs solution to a Vector of 0s passed into it
void->void
*/

void Matrix::jacobi(Matrix::Vector* constants, Matrix::Vector* solution) {

    int processCount = main_simpi->getProcessCount();
    int id = main_simpi->getID();

    Matrix::Vector *prev = new Matrix::Vector(constants->getSize()); // shared mem containing a copy of values

    Matrix* saveEq = new Matrix(get_x(), get_y() + 1); // save equations from modification
    Matrix::Vector* saveConst = new Matrix::Vector(constants->getSize()); // saves input vector


    int work = constants->getSize() / processCount;

    int i, j, k;
    int start = id * work;
    int end = start + work;

    main_simpi->synch();

    //Save Matrix and Vector
    for (i = start; i < end; i++) {
        for (j = 0; j < get_y(); j++) {
            saveEq->get(i,j) = get(i, j);
        }
        saveEq->get(i,get_y() + 1) = constants->getRef(i);
        saveConst->getRef(i) = constants->getRef(i);
    }

    //synch, wait for all process before solving
    main_simpi->synch();

    //setup, switch var coefficient with row solution, and divide by coefficient
    for (i = start; i < end; i++) {
        double temp = get(i, i);
        get(i, i) = constants->getRef(i);
        constants->getRef(i) = temp;
        for (j = 0; j < get_y(); j++) {
            if (j != i) {
                get(i, j) *= -1;
            }
            get(i, j) /= constants->getRef(i);
            //main_simpi->synch();
        }
        prev->getRef(i) = 1.0;
        solution->getRef(i) = constants->getRef(i);
        main_simpi->synch();
    }

    main_simpi->synch();
    // first iteration by trying substituting 1
    for (i = start; i < end; i++) {
        double rowSum = 0;
        for (j = 0; j < get_y(); j++) {
            if (j == i) {
                rowSum += get(i, j);
            } else {
                rowSum += (get(i, j) * prev->getRef(j));
            }
            //main_simpi->synch();
        }
        solution->getRef(i) = rowSum;
        main_simpi->synch();
    }

    //wait for all processes before repeating iterations with calculated results
    //main_simpi->synch();

    for (k = 0; k < 1000; k++)
    {
        for (i = start; i < end; i++)
        {
            //save prev value for comparision
            prev->getRef(i) = solution->getRef(i);
            main_simpi->synch();
        }
        for (i = start; i < end; i++)
        {
            double rowSum = 0;
            for (j = 0; j < get_y(); j++)
            {
                if (j == i) {
                    rowSum += get(i, j);
                } else {
                    rowSum += (get(i, j) * prev->getRef(j));
                }
            }
            solution->getRef(i) = rowSum;
            //main_simpi->synch();
        }
        //wait at end of each loop for all processes before beginning next iteration
        main_simpi->synch();
    }
    main_simpi->synch();
    //restore original Matrix and Vector
    for (i = start; i < end; i++) {
        for (j = 0; j < get_y(); j++) {
            get(i, j) = saveEq->get(i,j);
        }
        constants->getRef(i) = saveConst->getRef(i);
    }
    //wait for all processes before returning solution Vector
    main_simpi->synch();
    //return solution;
    return;
}


/*
Checks if a square Matrix is diagonally dominant
(diagonal terms are greater-than or equal to sum of the rest of their row)
none->bool
*/
bool Matrix::isDiagonallyDominant()
{
    for(int i = 0; i < get_x(); i ++)
    {
        double sq;
        double rest = 0;
        for(int j = 0; j < get_y(); j++)
        {
            if(i==j)
            {
                sq = get(i,j);
            }
            else
            {
                rest+=get(i,j);
            }
        }
        if (sq < rest)
        {
            return false;
        }
    }
    return true;
}


/*
Solves a linear system of equations
If the Matrix is diagonally dominant, jacobi-iterative method is used
else, the inverse mutliplication method is used
takes in a Vector of constants that each row is to be solved for and a solution Vector of 0s in which the solution will be written in
void -> void
*/
void Matrix::solveSystem(Matrix::Vector *constants, Matrix::Vector* solution)
{
    bool dd = isDiagonallyDominant();
    main_simpi->synch();
    if (dd)
    {
        if (main_simpi->getID() == 0) { std::cout << "jacobi" << std::endl; }
        main_simpi->synch();
        jacobi(constants, solution);
    }
    else
    {
        if (main_simpi->getID() == 0) { std::cout << "failsafe" << std::endl; }        
        main_simpi->synch();
        failSafe(constants, solution);
    }
    main_simpi->synch();
}

/*
Method to solve a system of linear equations if the system is not diagonally dominant
Uses the inverse-mutliplication method. 
void->void
*/
void Matrix::failSafe(Matrix::Vector* constants, Matrix::Vector* solution)
{
    Matrix* inv = new Matrix(get_x(), get_y());
    inverse(inv);
    std::cout << "inverse calculated" << std::endl;
    main_simpi->synch();

    int processCount = main_simpi->getProcessCount();
    int id = main_simpi->getID();
    int work = constants->getSize() / processCount;

    int start = id * work;
    int end = start + work;
    int n = constants->getSize();
    
    double sol;
    for(int i = start; i < end; i++)
    {
        sol = 0;
        for(int j = 0; j < n; j++)
        {
            sol += (inv->get(i, j)*constants->getRef(j));
        }
        solution->getRef(i) = sol;
    }
    main_simpi->synch();
    return;
}

Matrix::Vector::Vector(int a)
{
  // use simp and init the Matrix for all processes. The id is also in simp
  std::pair<std::string, double*> pass_back(Matrix::main_simpi->create_matrix(1, a));
  unique_id = pass_back.first;
  arr = pass_back.second;
  dim = a;
}
Matrix::Vector::~Vector() 
{
  // use main_simpi for getting rid of the mem and unlink stuff
  Matrix::main_simpi->free_matrix(unique_id);
}

std::ostream& operator<<(std::ostream& out, const Matrix::Vector& v)
{
    
    if (v.getSimpiID() == 0)
    {
        out << std::endl;
        for (int i = 0; i < v.getSize(); i++)
            out << std::fixed << std::setprecision(2) << v.getVal(i) << ", " << std::endl;
    }
    return out;
}

