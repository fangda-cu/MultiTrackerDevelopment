#include "MKLMatrix.hh"
#include "mkl_cblas.h"


inline MKLMatrix::MKLMatrix(int s)
  : MatrixBase(s, s)
{
  m_data.resize( m_rows * m_cols );
}

inline MKLMatrix::MKLMatrix(int r, int c, int nnz)
  : MatrixBase(r, c)
{
  m_data.resize( m_rows * m_cols );
}

inline MKLMatrix::MKLMatrix(const MKLMatrix& M)
  : MatrixBase(0, 0)
{
  m_rows = M.numRows();
  m_cols = M.numColumns();
  m_data = M.copyOfData();
}

inline MKLMatrix::~MKLMatrix()
{
  m_data.clear();
}

inline Scalar MKLMatrix::operator() (int r, int c) const
{
  return (Scalar)m_data[ r * m_cols + c ];
}

inline int MKLMatrix::set(int r, int c, Scalar val)
{
  m_data[ r * m_cols + c ] = val;
  return 0;
}

inline int MKLMatrix::add(int r, int c, Scalar val)
{
  m_data[ r * m_cols + c ] += val;
  return 0;
}

inline int MKLMatrix::add(const IntArray& rowIdx, const IntArray& colIdx,
                            const MatXd& values)
{
  int numRows = rowIdx.size();
  int numCols = colIdx.size();
  
  for ( size_t c=0; c<numCols; c++ )
    for ( size_t r=0; r<numRows; r++ )
      m_data[ rowIdx[r] * m_cols + colIdx[c] ] +=  values(r,c);
  
    
  /*int ierr;
  ierr = MatSetValues(m_M, rowIdx.size(), &(rowIdx[0]),
                      colIdx.size(), &(colIdx[0]),
                      values.data(), ADD_VALUES);
  CHKERRQ(ierr);*/
  return 0;
}

inline int MKLMatrix::add(const IndexArray& rowIdx, const IndexArray& colIdx,
          const MatXd& values)
{
  int numRows = rowIdx.size();
  int numCols = colIdx.size();
  
  for ( size_t c=0; c<numCols; c++ )
    for ( size_t r=0; r<numRows; r++ )
      m_data[ rowIdx.data()[r] * m_cols + colIdx.data()[c] ] +=  values(r,c);
  
    
  /*int ierr;
  ierr = MatSetValues(m_M, rowIdx.size(), rowIdx.data(),
                      colIdx.size(), colIdx.data(),
                      values.data(), ADD_VALUES);
 ;*/
  return 0;
}

inline int MKLMatrix::scale(Scalar val)
{
  size_t numEntries = m_data.size();
  for ( size_t i=0; i< numEntries; i++ )
    m_data[i] *= val;
  return 0;
}

// y = s*M*x
inline int MKLMatrix::multiply(VecXd& y, Scalar s, const VecXd& x)
{
    cerr << "multiple s*M*x\n";
    
  assert(cols() == x.size());
  assert(rows() == y.size());

  /*PetscUtils::copyToPetscVector(m_v, x);

  MatMult(m_M, m_v, m_w);

  PetscUtils::addFromPetscVector(y, s, m_w);*/
    
  vector<Scalar> a;
  vector<Scalar> b;
  a.resize(m_cols);
  b.resize(m_rows);
  
  for (size_t c=0; c<m_cols; c++)
    a[c] = x[c];
  
  // FIXME:
  // Are Eigen vectors guaranteed to eb contiguos in memory, like STL vectors? If so
  // then we can just pass them into dgemv. But to be safe we'll stick them in a vector
  // for just now.
  
   MKL_INT m = m_rows;
   MKL_INT n = m_cols;
   double alpha = s;
   double beta = 1;
   
   // If we are calling this with a specific matrix type then there are more efficient vector/matrix
   // multiplaction functions.

   cblas_dgemv( CblasRowMajor, CblasNoTrans, m, n, alpha, &m_data[0], m_cols, 
                &a[0], 1, beta, &b[0], 1 );
               
   for (size_t r=0; r<m_cols; r++)
     y[r] = b[r];
  
                                      
   /*PetscUtils::copyPetscVec(a, _v);
   PetscUtils::createPetscVec(b, _w);

   // multiply
   MatMult(_M, _v, _w);

   // add the result to w
   Real* wDat;
   int ierr = VecGetArray(_w, &wDat);
   CHKERRQ(ierr);
   for (uint i = 0; i < b.size(); i++) b[i] += r * wDat[i];
   ierr = VecRestoreArray(_w, &wDat);
   CHKERRQ(ierr);*/

   return 0;
}

  
inline int MKLMatrix::setZero()
{
  std::fill(m_data.begin(), m_data.end(), 0); 
  return 0;
}

inline int MKLMatrix::zeroRows(const IntArray& idx, Scalar diag)
{
  for ( size_t r=0; r<idx.size(); r++ )
    for ( size_t c=0; c<m_cols; c++ )
    {
      if ( idx[r] == c )
        m_data[ idx[r] * m_cols + c ] = diag;
      else
        m_data[ idx[r] * m_cols + c ] = 0.0;
    }
  return 0;
}

inline int MKLMatrix::finalize()
{
  /*int ierr = MatAssemblyBegin(m_M, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatAssemblyEnd(m_M, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatSetOption(m_M, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);
  CHKERRQ(ierr);
  ierr = MatSetOption(m_M, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
  CHKERRQ(ierr);*/
  return 0;
}

inline void MKLMatrix::transpose()
{
    vector<double> copy = copyOfData();

    for (int i_row=0; i_row<m_rows; ++i_row)
        for (int i_column=0; i_column<m_cols; ++i_column)
            m_data[i_column * m_rows + i_row] = copy[i_row * m_cols + i_column];

    std::swap(m_rows, m_cols);
}

