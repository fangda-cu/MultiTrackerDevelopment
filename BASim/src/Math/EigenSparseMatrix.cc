#include "EigenSparseMatrix.hh"

namespace BASim {

int EigenSparseMatrix::setZero()
{
  //iterate through the matrix and explicitly zero the entries (this should leave the structure intact?)
  /*for (int k=0; k<m_dynamic.outerSize(); ++k) {
    for (Eigen::SparseMatrix<Scalar,Eigen::RowMajor>::InnerIterator it(m_dynamic,k); it; ++it)
    {
      m_dynamic.coeffRef(it.row(), it.col()) = 0;
    }
  }*/
  //TODO Maybe do something else here?
  m_triplets.clear();

  return 0;
}

int EigenSparseMatrix::resetNonzeros() {
   
  //m_dynamic = Eigen::SparseMatrix<Scalar,Eigen::RowMajor>(m_rows, m_cols);
  m_triplets.clear();
   return 0;
}

int EigenSparseMatrix::zeroRows(const IntArray& idx, Scalar diag)
{
  
  //NOTE: This zeros rows AND symmetric columns, since that seems quickest.
  for(unsigned int i = 0; i < idx.size(); ++i) {
    int row = idx[i];
    //iterate across the row zeroing entries
    for (Eigen::SparseMatrix<Scalar,Eigen::RowMajor>::InnerIterator it(m_dynamic,row); it; ++it)
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
  /*
  //NOTE: This zeros rows AND symmetric columns
  std::vector< Eigen::Triplet<Scalar> > tripletsNew;
  tripletsNew.reserve(m_triplets.size());

  for(unsigned int i = 0; i < m_triplets.size(); ++i) {
    int row = m_triplets[i].row();
    int col = m_triplets[i].col();
    if(std::find(idx.begin(), idx.end(), row) != idx.end() || std::find(idx.begin(), idx.end(), col) != idx.end()) {
      if(col == row)
        tripletsNew.push_back(Eigen::Triplet<Scalar>() );
    }
    else {
      tripletsNew.push_back(m_triplets[i]);
    }
  }
 
 m_triplets = tripletsNew;
 */

  return 0;
}


int EigenSparseMatrix::finalize()
{
   return 0;
}

int EigenSparseMatrix::finalizeNonzeros() {
  
  m_dynamic.setFromTriplets(m_triplets.begin(), m_triplets.end());
  return 0;
}


int EigenSparseMatrix::scale(Scalar val)
{
  //replace all values with scaled versions
  for(unsigned int i = 0; i < m_triplets.size(); ++i) {
    int row = m_triplets[i].row();
    int col = m_triplets[i].col();
    Scalar value = m_triplets[i].value();
    m_triplets[i] = Eigen::Triplet<Scalar>(row,col,value*val);
  }

  return 0;
}

int EigenSparseMatrix::multiply(VecXd& y, Scalar s, const VecXd& x) const
{
  //compute y += s*A*x
  assert(cols() == x.size());
  assert(rows() == y.size());
  
  //Use the non-dynamic sparse matrix for multiplying (i.e. within conjugate gradient)
  //since it's expected to be faster. Assumes finalizeNonZeros was called already.
  y += s*(m_dynamic*x); 
  
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