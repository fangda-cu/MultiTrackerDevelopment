/**
 * \file SimpleSparseMatrix.hh
 *
 * \author batty@cs.columbia.edu
 * \date 12/10/2011
 */

#ifndef SIMPLESPARSEMATRIX_HH
#define SIMPLESPARSEMATRIX_HH

#include "BASim/src/Math/MatrixBase.hh"

namespace BASim {

class SimpleSparseMatrix : public MatrixBase
{
   
public:

  SimpleSparseMatrix(int s);
  SimpleSparseMatrix(int r, int c, int nnz = 1);
  SimpleSparseMatrix(const SimpleSparseMatrix& M);
  ~SimpleSparseMatrix();

  Scalar operator() (int i, int j) const;
  int set(int r, int c, Scalar val);
  int add(int r, int c, Scalar val);
  int add(const IntArray& rowIdx, const IntArray& colIdx, const MatXd& values);
  int add(const IndexArray& rowIdx, const IndexArray& colIdx,
          const MatXd& values);
  int scale(Scalar val);
  int setZero();
  int zeroRows(const IntArray& idx, Scalar diag = 1.0);
  int multiply(VecXd& y, Scalar s, const VecXd& x) const;

  //junk from the rods side of things - stubs for now
  void vertexStencilAdd(int start, const Eigen::Matrix<Scalar, 11, 11>& localJ) {}
  void edgeStencilAdd(int start, const Eigen::Matrix<Scalar, 6, 6>& localJ) {}
  void pointStencilAdd(int start, const Eigen::Matrix<Scalar, 3, 3>& localJ) {}

  int zeroCols(const IntArray& idx, Scalar diag);
  bool isApproxSymmetric( Scalar eps ) const;
  std::string name() const;

  int resetNonzeros();


  int finalize();
  int finalizeNonzeros();

protected:

  //A barebones, naive matrix structure.
  std::vector< std::map<int, Scalar> > m_data;
};

#include "SimpleSparseMatrix.inl"

} // namespace BASim

#endif // SIMPLESPARSEMATRIX_HH
