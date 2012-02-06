/**
 * \file PardisoLinearSolver.cc
 *
 * \author smith@cs.columbia.edu
 * \date 05/29/2010
 */

#include "PardisoMKLSolver.hh"
#include "BASim/src/Math/EigenSparseMatrix.hh"

#include "mkl_pardiso.h"
#include "mkl_types.h"

using namespace std;

namespace BASim {

PardisoMKLSolver::PardisoMKLSolver( EigenSparseMatrix& A )
: LinearSolverBase(A)
, prdsomat(A)
{
  assert(A.rows() == A.cols());
}

PardisoMKLSolver::~PardisoMKLSolver()
{
}

void PardisoMKLSolver::parsePardisoError( int error ) const 
{
  switch( error )
  {
    case 0:
    {
      //std::cout << "No error." << std::endl;
      break;
    }
    case -1:
    {
      std::cerr << "Input inconsistent. Exiting." << std::endl;
      exit(0);
      break;
    }
    case -2:
    {
      std::cerr << "Not enough memory. Exiting." << std::endl;
      exit(0);
      break;
    }
    case -3:
    {
      std::cerr << "Reordering problem. Exiting." << std::endl;
      exit(0);
      break;
    }
    case -4:
    {
      std::cerr << "Zero pivot, numerical fact. or iterative refinement problem. Exiting." << std::endl;
      exit(0);
      break;
    }
    case -5:
    {
      std::cerr << "Unclassified (internal) error. Exiting." << std::endl;
      exit(0);
      break;
    }
    case -6:
    {
      std::cerr << "Preordering failed (matrix types 11, 13 only). Exiting." << std::endl;
      exit(0);
      break;
    }
    case -7:
    {
      std::cerr << "Diagonal matrix problem. Exiting." << std::endl;
      exit(0);
      break;
    }
    case -8:
    {
      std::cerr << "32-bit integer overflow problem. Exiting." << std::endl;
      exit(0);
      break;
    }
    case -10:
    {
      std::cerr << "No license file pardiso.lic found. Exiting." << std::endl;
      exit(0);
      break;
    }
    case -11:
    {
      std::cerr << "License is expired. Exiting." << std::endl;
      exit(0);
      break;
    }
    case -12:
    {
      std::cerr << "Wrong username or hostname. Exiting." << std::endl;
      exit(0);
      break;
    }
    case -100:
    {
      std::cerr << "Reached maximum number of Krylov-subspace iteration in iterative solver. Exiting." << std::endl;
      exit(0);
      break;
    }
    case -101:
    {
      std::cerr << "No sufficient convergence in Krylov-subspace iteration within 25 iterations. Exiting." << std::endl;
      exit(0);
      break;
    }
    case -102:
    {
      std::cerr << "Error in Krylov-subspace iteration. Exiting." << std::endl;
      exit(0);
      break;
    }
    case -103:
    {
      std::cerr << "Break-Down in Krylov-subspace iteration. Exiting." << std::endl;
      exit(0);
      break;
    }
    default:
    {
      std::cerr << "INVALID ERROR CODE ENCOUNTERED. THIS IS A BUG. EXITING." << std::endl;
      exit(0);
      break;
    }
  }
}

#if defined(MKL_ILP64)
#define MKL_INT long long
#else
#define MKL_INT int
#endif

int PardisoMKLSolver::solve( VecXd& x_in, const VecXd& b_in )
{

  //convert Eigen matrix to compressed sparse row format

  /* Matrix data. */
  Eigen::SparseMatrix<Scalar,Eigen::RowMajor> mat = prdsomat.getEigenMatrix();
  
  /*
  //Test data in PARDISO format
  MKL_INT nT = 8;
  MKL_INT nnzT = 18;
  MKL_INT iaT[ 9 ] = { 1, 5, 8, 10, 12, 15, 17, 18, 19 };
  MKL_INT jaT[18] = { 1, 3, 6, 7,
    2, 3, 5,
    3, 8,
    4, 7,
    5, 6, 7,
    6, 8,
    7,
    8 };
  double aT[18] = { 7.0, 1.0, 2.0, 7.0,
    -4.0, 8.0, 2.0,
    1.0, 5.0,
    7.0, 9.0,
    5.0, 1.0, 5.0,
    -1.0, 5.0,
    11.0,
    5.0 };

  std::cout << "Converting test data to 0-indexed Eigen DynamicSparseMatrix\n";
  Eigen::DynamicSparseMatrix<Scalar, Eigen::RowMajor> mat2(nT,nT);
  int counter = 0;
  for(int row = 0; row < nT; ++row) {
    while(counter < nnzT && counter < iaT[row+1]-1) {
      mat2.coeffRef(row, jaT[counter]-1) = aT[counter];
      ++counter;
    }
  }

  std::cout << "Converting to Eigen SparseMatrix\n";
  //switch to sparsematrx
  Eigen::SparseMatrix<Scalar,Eigen::RowMajor> mat = mat2;
  */

  MKL_INT   n = mat.rows();
  MKL_INT* ia = new MKL_INT[n+1];
  MKL_INT* ja = new MKL_INT[mat.nonZeros()];
  double*   a = new double[mat.nonZeros()];

  //now fill in the data by walking through the eigen sparse matrix
  std::cout << "Filling matrix\n";
  int counter = 0;
  for (int k=0; k<mat.outerSize(); ++k) {
    ia[k] = counter+1;
    for (Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator it(mat,k); it; ++it)
    {
      ja[counter] = it.col()+1;
      a[counter] = it.value();
      ++counter;
    }
  }
  ia[n] = counter+1;
  
  /*
  std::cout << "Verify that the data is consistent.\n";
  std::cout << "ia vector: ";
  for(int i = 0; i < n+1; ++i) {
    std::cout << ia[i] << " ";
  }
  std::cout << "\nja vector: ";
  for(int i = 0; i < mat.nonZeros(); ++i) {
    std::cout << ja[i] << " ";
  }
  std::cout << "\ndata vector: ";
  for(int i = 0; i < mat.nonZeros(); ++i) {
    std::cout << a[i] << " ";
  }  
  */
  std::cout << "Starting up PARDISO\n";
  



  //MKL_INT mtype = -2; /* Real symmetric matrix */
  MKL_INT mtype = 11; /* Real unsymmetric matrix */
  /* RHS and solution vectors. */
  //double b[8], x[8];
  double* b =  new double[b_in.size()];
  for(int i = 0; i < b_in.size(); ++i)
    b[i] = b_in(i);
  double* x = new double[b_in.size()];

  MKL_INT nrhs = 1; /* Number of right hand sides. */
  /* Internal solver memory pointer pt, */
  /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
  /* or void *pt[64] should be OK on both architectures */
  void *pt[64];
  /* Pardiso control parameters. */
  MKL_INT iparm[64];
  MKL_INT maxfct, mnum, phase, error, msglvl;
  /* Auxiliary variables. */
  MKL_INT i;
  double ddum; /* Double dummy */
  MKL_INT idum; /* Integer dummy. */
  /* -------------------------------------------------------------------- */
  /* .. Setup Pardiso control parameters. */
  /* -------------------------------------------------------------------- */
  for (i = 0; i < 64; i++) {
    iparm[i] = 0;
  }
  iparm[0] = 1; /* No solver default */
  iparm[1] = 2; /* Fill-in reordering from METIS */
  /* Numbers of processors, value of OMP_NUM_THREADS */
  iparm[2] = 1;
  iparm[3] = 0; /* No iterative-direct algorithm */
  iparm[4] = 0; /* No user fill-in reducing permutation */
  iparm[5] = 0; /* Write solution into x */
  iparm[6] = 0; /* Not in use */
  iparm[7] = 2; /* Max numbers of iterative refinement steps */
  iparm[8] = 0; /* Not in use */
  iparm[9] = 13; /* Perturb the pivot elements with 1E-13 */
  iparm[10] = 1; /* Use nonsymmetric permutation and scaling MPS */
  iparm[11] = 0; /* Not in use */
  iparm[12] = 1; /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
  iparm[13] = 0; /* Output: Number of perturbed pivots */
  iparm[14] = 0; /* Not in use */
  iparm[15] = 0; /* Not in use */
  iparm[16] = 0; /* Not in use */
  iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
  iparm[18] = -1; /* Output: Mflops for LU factorization */
  iparm[19] = 0; /* Output: Numbers of CG Iterations */
  maxfct = 1; /* Maximum number of numerical factorizations. */
  mnum = 1; /* Which factorization to use. */
  msglvl = 1; /* Print statistical information in file */
  error = 0; /* Initialize error flag */
  /* -------------------------------------------------------------------- */
  /* .. Initialize the internal solver memory pointer. This is only */
  /* necessary for the FIRST call of the PARDISO solver. */
  /* -------------------------------------------------------------------- */
  for (i = 0; i < 64; i++) {
    pt[i] = 0;
  }
  /* -------------------------------------------------------------------- */
  /* .. Reordering and Symbolic Factorization. This step also allocates */
  /* all memory that is necessary for the factorization. */
  /* -------------------------------------------------------------------- */
  std::cout << "Reorder/Factor\n";
  phase = 11;
  PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
    &n, a, ia, ja, &idum, &nrhs,
    iparm, &msglvl, &ddum, &ddum, &error);
  if (error != 0) {
    printf("\nERROR during symbolic factorization: %d", error);
    exit(1);
  }
  printf("\nReordering completed ... ");
  printf("\nNumber of nonzeros in factors = %d", iparm[17]);
  printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
  /* -------------------------------------------------------------------- */
  /* .. Numerical factorization. */
  /* -------------------------------------------------------------------- */
  phase = 22;
  PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
    &n, a, ia, ja, &idum, &nrhs,
    iparm, &msglvl, &ddum, &ddum, &error);
  if (error != 0) {
    printf("\nERROR during numerical factorization: %d", error);
    exit(2);
  }
  printf("\nFactorization completed ... ");
  /* -------------------------------------------------------------------- */
  /* .. Back substitution and iterative refinement. */
  /* -------------------------------------------------------------------- */
  phase = 33;
  iparm[7] = 2; /* Max numbers of iterative refinement steps. */
  PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
    &n, a, ia, ja, &idum, &nrhs,
    iparm, &msglvl, const_cast<double *>(b_in.data()), x_in.data(), &error);
  if (error != 0) {
    printf("\nERROR during solution: %d", error);
    exit(3);
  }
  printf("\nSolve completed ... ");
  
  ////copy data out.
  //for(int i = 0; i < n; ++i)
  //  x_in(i) = x[i];

  /* -------------------------------------------------------------------- */
  /* .. Termination and release of memory. */
  /* -------------------------------------------------------------------- */
  phase = -1; /* Release internal memory. */
  PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
    &n, &ddum, ia, ja, &idum, &nrhs,
    iparm, &msglvl, &ddum, &ddum, &error);
  
  //delete[] b;
  //delete[] x;
  delete[] ia;
  delete[] ja;
  delete[] a;

  return 0;
}

} // namespace BASim
