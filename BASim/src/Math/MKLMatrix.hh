/**
 * \file MKLMatrix.hh
 *
 * \author acoull@wetafx.co.nz
 * \date 11/10/2009
 */

#ifndef MKL_MATRIX_HH
#define MKL_MATRIX_HH

#include <vector>

using namespace std;

namespace BASim {

class MKLMatrix : public MatrixBase
{
public:
  MKLMatrix(int s);
  MKLMatrix(int r, int c, int nnz = 1);
  MKLMatrix(const MKLMatrix& M);
  ~MKLMatrix();

  Scalar operator() (int i, int j) const;
  int set(int r, int c, Scalar val);
  int add(int r, int c, Scalar val);
  int add(const IntArray& rowIdx, const IntArray& colIdx, const MatXd& values);
  int add(const IndexArray& rowIdx, const IndexArray& colIdx,
          const MatXd& values);
  int scale(Scalar val);
  int setZero();
  int zeroRows(const IntArray& idx, Scalar diag = 1.0);
  int multiply(VecXd& y, Scalar s, const VecXd& x);
    
  int finalize();
  void transpose();
  
  //const Mat& getMKLMatrix() const { return m_M; }
  //Mat& getMKLMatrix() { return m_M; }

  size_t numRows() const { return m_rows; }
  size_t numColumns() const { return m_cols; }
      
  vector<Scalar>& getData() { return m_data; }
  vector<Scalar> copyOfData() const { return m_data; }
  
protected:

  vector<Scalar> m_data;
  //Mat m_M;       // the underlying Petsc matrix
  //Vec m_v, m_w;  // m_v is size rows(), m_w is size cols()
};

#include "MKLMatrix.inl"

} // namespace BASim

#endif

