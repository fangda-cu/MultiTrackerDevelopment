#ifndef BANDMATRIX_HH
#define BANDMATRIX_HH

/**
 * \file BandMatrix.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 11/16/2009
 * Modified by smith@cs.columbia.edu
 *   06/26/2010
 */

#include "MatrixBase.hh"
#include "../Core/Util.hh"

namespace BASim
{

/** Class for storing band matrices in a two-dimensional array. */
class BandMatrix: public MatrixBase
{
public:

    /** Creates an m-by-n band matrix with kl non-zero sub-diagonals and
     ku non-zero super-diagonals. */
    BandMatrix(int m, int n, int kl, int ku) :
        MatrixBase(m, n), m_kl(kl), m_ku(ku), m_size((kl + ku + 1) * n)
    {
        m_data = new Scalar[m_size];
        setZero();

        m_lower = new int[MatrixBase::m_rows];
        m_upper = new int[MatrixBase::m_rows];
        for (int i = 0; i < MatrixBase::m_rows; ++i)
        {
            m_lower[i] = std::max(i - m_kl, 0);
            m_upper[i] = std::min(i + m_ku + 1, MatrixBase::m_cols);
        }
    }

    virtual ~BandMatrix()
    {
        if (m_data != NULL)
            delete[] m_data;
        if (m_lower != NULL)
            delete[] m_lower;
        if (m_upper != NULL)
            delete[] m_upper;
    }

    virtual Scalar operator()(int i, int j) const
    {
        if (!indicesValid(i, j))
            return 0;
        return m_data[(m_ku + i - j) * MatrixBase::cols() + j];
    }

    Scalar& operator()(int i, int j)
    {
        assert(indicesValid(i, j));
        return m_data[(m_ku + i - j) * MatrixBase::cols() + j];
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

    virtual int add(const IntArray& rowIdx, const IntArray& colIdx, const MatXd& values)
    {
        int nr = (int) rowIdx.size();
        int nc = (int) colIdx.size();

        for (int i = 0; i < nr; ++i)
        {
            int r = rowIdx[i];
            assert(r >= 0);
            assert(r < rows());
            Scalar* val = &m_data[(m_ku + r) * MatrixBase::m_cols];
            for (int j = 0; j < nc; ++j)
            {
                assert(indicesValid(r, colIdx[j]));
                *(val + colIdx[j] * (1 - MatrixBase::m_cols)) += values(i, j);
            }
        }

        return 0;
    }

    virtual int add(const IndexArray& rowIdx, const IndexArray& colIdx, const MatXd& values)
    {
        int nr = rowIdx.size();
        int nc = colIdx.size();

        for (int i = 0; i < nr; ++i)
        {
            int r = rowIdx[i];
            Scalar* val = &m_data[(m_ku + r) * MatrixBase::m_cols];
            for (int j = 0; j < nc; ++j)
            {
                assert(indicesValid(r, colIdx[j]));
                *(val + colIdx[j] * (1 - MatrixBase::m_cols)) += values(i, j);
            }
        }

        return 0;
    }

    // Add the local Jacobian to the banded matrix following the edge stencil.
    // A zero row and column are inserted in the middle of localJ before adding it at "start" position on the diagonal.
    virtual void edgeStencilAdd(int start, const Eigen::Matrix<Scalar, 6, 6>& localJ)
    {
        static const int size = 6;

        const int n = MatrixBase::m_cols;
        start += n * m_ku;

        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                m_data[start] += localJ(i, j);
                start += 1 - n;
            }
            start += 1 - n;

            for (int j = 3; j < 6; ++j)
            {
                m_data[start] += localJ(i, j);
                start += 1 - n;
            }
            start += (size + 1) * (n - 1) + n;
        }
        start += n;

        for (int i = 3; i < 6; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                m_data[start] += localJ(i, j);
                start += 1 - n;
            }
            start += 1 - n;

            for (int j = 3; j < 6; ++j)
            {
                m_data[start] += localJ(i, j);
                start += 1 - n;
            }
            start += (size + 1) * (n - 1) + n;
        }
    }

    // Add the local Jacobian to the banded matrix following the vertex stencil.
    // LocaJ is directly added at "start" position on the diagonal.
    virtual void vertexStencilAdd(int start, const Eigen::Matrix<Scalar, 11, 11>& localJ)
    {
        static const int size = 11;

        const int n = MatrixBase::m_cols;
        start += n * m_ku;

        for (int i = 0; i < size; ++i)
        {
            for (int j = 0; j < size; ++j)
            {
                m_data[start] += localJ(i, j);
                start += 1 - n;
            }
            start += size * n + n - size;
        }
    }

    virtual void pointStencilAdd(int start, const Eigen::Matrix<Scalar, 3, 3>& localJ)
    {
        static const int size = 3;

        const int n = MatrixBase::m_cols;
        start += n * m_ku;

        for (int i = 0; i < size; ++i)
        {
            for (int j = 0; j < size; ++j)
            {
                m_data[start] += localJ(i, j);
                start += 1 - n;
            }
            start += size * n + n - size;
        }
    }

    virtual int scale(Scalar val)
    {
        for (int i = 0; i < m_size; ++i)
            m_data[i] *= val;
        return 0;
    }

    virtual int setZero()
    {
        for (int i = 0; i < m_size; ++i)
            m_data[i] = 0;
        return 0;
    }

    virtual int zeroRows(const IntArray& idx, Scalar diag = 1.0)
    {
        for (int i = 0; i < (int) idx.size(); ++i)
        {
            int r = idx[i];
            assert(r >= 0 && r < MatrixBase::m_rows);
            int lower = m_lower[r];
            int upper = m_upper[r];
            for (int j = lower; j < upper; ++j)
            {
                (*this)(r, j) = 0;
            }
            (*this)(r, r) = diag;
        }
        return 0;
    }

    int computeStartOfCol(int col)
    {
        assert(col >= 0);
        assert(col < MatrixBase::cols());
        return std::max(col - ku(), 0);
    }

    int computeEndOfCol(int col)
    {
        assert(col >= 0);
        assert(col < MatrixBase::cols());
        return std::min(col + kl(), MatrixBase::rows() - 1);
    }

    virtual int zeroCols(const IntArray& idx, Scalar diag)
    {
#ifdef DEBUG
        for( int i = 0; i < (int) idx.size(); ++i ) assert( idx[i] >= 0 );
        for( int i = 0; i < (int) idx.size(); ++i ) assert( idx[i] < MatrixBase::cols() );
#endif

        // For each column the user provided
        for (int i = 0; i < (int) idx.size(); ++i)
        {
            int col = idx[i];
            for (int row = computeStartOfCol(col); row <= computeEndOfCol(col); ++row)
            {
                assert(indicesValid(row, col));
                if (row != col)
                    (*this)(row, col) = 0.0;
                else
                    (*this)(row, col) = diag;
            }
        }

        return 0;
    }

    virtual int multiply(VecXd& y, Scalar s, const VecXd& x) const
    {
        assert(y.size() == MatrixBase::m_rows);
        assert(x.size() == MatrixBase::m_cols);

        for (int i = 0; i < MatrixBase::m_rows; ++i)
        {
            int lower = m_lower[i];
            int upper = m_upper[i];
            const Scalar* val = &m_data[(m_ku + i - lower) * MatrixBase::m_cols + lower];
            Scalar sum = 0;
            for (int j = lower; j < upper; ++j, val += (1 - MatrixBase::m_cols))
            {
                sum += (*val) * x[j];
            }
            y[i] += s * sum;
        }
        return 0;
    }

    virtual void print() const
    {
        std::cout << "[";
        for (int i = 0; i < MatrixBase::m_rows; ++i)
        {
            for (int j = 0; j < MatrixBase::m_cols; ++j)
            {
                std::cout << (*this)(i, j);
                if (j < MatrixBase::m_cols - 1)
                    std::cout << ", ";
            }
            if (i < MatrixBase::m_rows - 1)
                std::cout << "; ";
        }
        std::cout << "]" << std::endl;
    }

    // NOTE: Assumes not symmetric if lower band is not same size as upper (could be a bunch of zeros in bigger band)
    virtual bool isApproxSymmetric(Scalar eps) const
    {
        if (MatrixBase::rows() != MatrixBase::cols())
            return false;
        if (kl() != ku())
            return false;
        for (int i = 0; i < MatrixBase::rows(); ++i)
        {
            for (int j = i + 1; j <= i + ku(); ++j)
            {
                if (!approxEq((*this)(i, j), (*this)(j, i), eps))
                    return false;
            }
        }
        return true;
    }

    /** Number of sub-diagonals */
    inline int kl() const
    {
        return m_kl;
    }

    /** Number of super-diagonals */
    inline int ku() const
    {
        return m_ku;
    }

    /** Returns underlying array of values */
    Scalar* data()
    {
        return m_data;
    }

    std::string name() const
    {
        return "BandMatrix";
    }

protected:

    bool indicesValid(int r, int c) const
    {
        return ((r >= 0) && (r < rows()) && (c >= 0) && (c < cols()) && (r - c <= m_kl) && (c - r <= m_ku));
    }

    int m_kl; ///< Number of sub-diagonals
    int m_ku; ///< Number of super-diagonals
    int m_size; ///< Size of the array holding the non-zero entries
    Scalar* m_data; ///< Array that stores the non-zero entries
    int* m_lower; ///< For each row, stores the smallest valid column index
    int* m_upper; ///< For each row, stores the largest valid column index
};

} // namespace BASim

#endif // BANDMATRIX_HH
