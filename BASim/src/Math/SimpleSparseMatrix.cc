#include "SimpleSparseMatrix.hh"

namespace BASim {

int SimpleSparseMatrix::setZero()
{
   for(int i = 0; i < rows(); ++i) {
      m_data[i].clear();
   }

  return 0;
}

int SimpleSparseMatrix::resetNonzeros() {
   return 0;
}

int SimpleSparseMatrix::zeroRows(const IntArray& idx, Scalar diag)
{

  for(unsigned int i = 0; i < idx.size(); ++i) {
    int row = idx[i];
    m_data[row].clear();
    m_data[row][row] = diag;
  }
  return 0;
}


int SimpleSparseMatrix::finalize()
{
   return 0;
}

int SimpleSparseMatrix::finalizeNonzeros() {
   return 0;
}


int SimpleSparseMatrix::scale(Scalar val)
{
  for(int r = 0; r < rows(); ++r) {
    std::map<int,Scalar>::iterator it = m_data[r].begin();
    for(; it != m_data[r].end(); ++it) {
      it->second *= val;
    }
  }

  return 0;
}

int SimpleSparseMatrix::multiply(VecXd& y, Scalar s, const VecXd& x) const
{
  assert(cols() == x.size());
  assert(rows() == y.size());

  for(int row = 0; row < x.size(); ++row)  {
    Scalar sum = 0;
    for(std::map<int,Scalar>::const_iterator it = m_data[row].begin(); it != m_data[row].end(); ++it) {
      int col = it->first;
      sum += it->second*x[col];
    }
    y[row] += s*sum;
  }

  return 0;
}

int SimpleSparseMatrix::zeroCols(const IntArray& idx, Scalar diag) 
{
   //iterate over the list of columns to clear
   for(unsigned int i = 0; i < idx.size(); ++i) {

      //go through all the rows, killing the appropriate column
      for(int row = 0; row < rows(); ++row) {
         
         if(row == idx[i]) //set diagonal, if this is the diagonal
            m_data[row][idx[i]] = diag;
         else { //otherwise delete the entry, if it exists
            std::map<int,Scalar>::iterator it=m_data[row].find(idx[i]);
            if(it != m_data[row].end())
               m_data[row].erase(it); 
         }
      }
   }

   return 0;
}


}