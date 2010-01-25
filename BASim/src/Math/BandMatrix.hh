#ifndef BANDMATRIX_HH
#define BANDMATRIX_HH

/**
 * \file BandMatrix.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 11/16/2009
 */

#include "MatrixBase.hh"

namespace BASim {

/** Class for storing band matrices in a two-dimensional array. */
class BandMatrix : public MatrixBase
{
public:

  /** Creates an m-by-n band matrix with kl non-zero sub-diagonals and
      ku non-zero super-diagonals. */
  BandMatrix(int m, int n, int kl, int ku)
    : MatrixBase(m, n)
    , m_kl(kl)
    , m_ku(ku)
    , m_size((kl + ku + 1) * n)
  {
    m_data = new Scalar[m_size];
    setZero();

    m_lower = new int[m_rows];
    m_upper = new int[m_rows];
    for (int i = 0; i < m_rows; ++i) {
      m_lower[i] = std::max(i - m_kl, 0);
      m_upper[i] = std::min(i + m_ku + 1, m_cols);
    }
    /*
      m_lower = new int[m_cols];
      m_upper = new int[m_cols];
      for (int j = 0; j < m_cols; ++j) {
      m_lower[j] = std::max(j - m_ku, 0);
      m_upper[j] = std::min(j + m_kl + 1, m_rows);
      }
    */
  }

  virtual ~BandMatrix()
  {
    delete [] m_data;
    delete [] m_lower;
    delete [] m_upper;
  }

  virtual Scalar operator() (int i, int j) const
  {
    if (j - i > m_ku) return 0;
    if (i - j > m_kl) return 0;
    return m_data[(m_ku + i - j) * cols() + j];
  }

  Scalar& operator() (int i, int j)
  {
    assert(j - i <= m_ku);
    assert(i - j <= m_kl);
    return m_data[(m_ku + i - j) * cols() + j];
  }

  virtual int set(int i, int j, Scalar val)
  {
    (*this)(i, j) = val;
    return 0;
  }

  virtual int add(int i, int j, Scalar val)
  {
    (*this)(i, j) += val;
    return 0;
  }

  virtual int add(const IntArray& rowIdx, const IntArray& colIdx,
                  const MatXd& values)
  {
    for (size_t i = 0; i < rowIdx.size(); ++i) {
      for (size_t j = 0; j < colIdx.size(); ++j) {
        (*this)(rowIdx[i], colIdx[j]) += values(i, j);
      }
    }
    return 0;
  }

  virtual int add(const IndexArray& rowIdx, const IndexArray& colIdx,
                  const MatXd& values)
  {
    for (int i = 0; i < rowIdx.size(); ++i) {
      for (int j = 0; j < colIdx.size(); ++j) {
        (*this)(rowIdx[i], colIdx[j]) += values(i, j);
      }
    }
    return 0;
  }

  virtual int scale(Scalar val)
  {
    for (int i = 0; i < m_size; ++i) m_data[i] *= val;
    return 0;
  }

  virtual int setZero()
  {
    for (int i = 0; i < m_size; ++i) m_data[i] = 0;
    return 0;
  }

  virtual int zeroRows(const IntArray& idx, Scalar diag = 1.0)
  {
    for (int i = 0; i < (int) idx.size(); ++i) {
      int r = idx[i];
      assert(r < m_rows);
      int lower = m_lower[r];
      int upper = m_upper[r];
      for (int j = lower; j < upper; ++j) {
        (*this)(r, j) = 0;
      }
      (*this)(r, r) = diag;
    }
    return 0;
  }

  virtual int multiply(VecXd& y, Scalar s, const VecXd& x) const
  {
    for (int i = 0; i < m_rows; ++i) {
      int lower = m_lower[i];
      int upper = m_upper[i];
      const Scalar* val = &m_data[(m_ku + i - lower) * m_cols + lower];
      Scalar sum = 0;
      for (int j = lower; j < upper; ++j, val += (1 - m_cols)) {
        sum += (*val) * x[j];
      }
      y[i] += s * sum;
    }
    return 0;
  }

  virtual void print() const
  {
    std::cout << "[";
    for (int i = 0; i < m_rows; ++i) {
      for (int j = 0; j < m_cols; ++j) {
        std::cout << (*this)(i, j);
        if (j < m_cols - 1) std::cout << ", ";
      }
      if (i < m_rows - 1) std::cout << "; ";
    }
    std::cout << "]" << std::endl;
  }

  /** Number of sub-diagonals */
  int kl() const { return m_kl; }

  /** Number of super-diagonals */
  int ku() const { return m_ku; }

  /** Returns underlying array of values */
  Scalar* data() { return m_data; }

protected:

  int m_kl; ///< Number of sub-diagonals
  int m_ku; ///< Number of super-diagonals
  int m_size; ///< Size of the array holding the non-zero entries
  Scalar* m_data; ///< Array that stores the non-zero entries
  int* m_lower; ///< For each row, stores the smallest valid column index
  int* m_upper; ///< For each row, stores the largest valid column index
};

} // namespace BASim

#endif // BANDMATRIX_HH
