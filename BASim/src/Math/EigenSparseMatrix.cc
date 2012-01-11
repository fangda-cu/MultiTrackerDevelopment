#include "EigenSparseMatrix.hh"

namespace BASim {

int EigenSparseMatrix::setZero()
{
  //iterate through the matrix and explicitly zero the entries (this should leave the structure intact)
  for (int k=0; k<m_dynamic.outerSize(); ++k) {
    for (Eigen::DynamicSparseMatrix<Scalar,Eigen::RowMajor>::InnerIterator it(m_dynamic,k); it; ++it)
    {
      m_dynamic.coeffRef(it.row(), it.col()) = 0;
    }
  }

  return 0;
}

int EigenSparseMatrix::resetNonzeros() {
   
  m_dynamic = Eigen::DynamicSparseMatrix<Scalar,Eigen::RowMajor>(m_rows, m_cols);

   return 0;
}

int EigenSparseMatrix::zeroRows(const IntArray& idx, Scalar diag)
{
  //NOTE: This zeros rows AND symmetric columns

  for(unsigned int i = 0; i < idx.size(); ++i) {
    int row = idx[i];
   
    //iterate across the row zeroing entries
    for (Eigen::DynamicSparseMatrix<Scalar,Eigen::RowMajor>::InnerIterator it(m_dynamic,row); it; ++it)
    {
      int col = it.col();
      if(col != row) {
        m_dynamic.coeffRef(row, col) = 0; //zero it out
        m_dynamic.coeffRef(col, row) = 0; //zero the opposite assuming symmetry, since doing zeroCols properly is slower.
      }
      else
        m_dynamic.coeffRef(row, col) = diag;
    }
    
  }
  return 0;
}


int EigenSparseMatrix::finalize()
{
   return 0;
}

int EigenSparseMatrix::finalizeNonzeros() {
  //Convert over to a regular sparse matrix for
  //doing fast multiplication
  m_static = Eigen::SparseMatrix<Scalar>(m_dynamic);
  return 0;
}


int EigenSparseMatrix::scale(Scalar val)
{
  m_dynamic *= val;

  return 0;
}

int EigenSparseMatrix::multiply(VecXd& y, Scalar s, const VecXd& x) const
{
  //compute y += s*A*x
  assert(cols() == x.size());
  assert(rows() == y.size());
  
  //Use the non-dynamic sparse matrix for multiplying (i.e. within conjugate gradient)
  //since it's expected to be faster. Assumes finalizeNonZeros was called already.
  //y += s*(m_dynamic*x); 
  y += s*(m_static*x);

  return 0;
}

int EigenSparseMatrix::zeroCols(const IntArray& idx, Scalar diag) 
{
   //iterate over the list of columns to clear. 
   //This is brute-forcey, but I don't know that Eigen (and the matrix structure itself) offers us any choice

  /*
  //Iterate through the list of colums
   for(unsigned int i = 0; i < idx.size(); ++i) {

      //go through all the rows, killing the appropriate column
      for(int row = 0; row < rows(); ++row) {
         
         if(row == idx[i]) //set diagonal, if this is the diagonal
            m_dynamic.coeffRef(row,idx[i]) = diag; 
         else { //otherwise zero the entry, if it's not zero
           if(m_dynamic.coeff(row,idx[i]) != 0)
             m_dynamic.coeffRef(row,idx[i]) = 0;
         }
      }
   }
   */ 

   return 0;
}


}