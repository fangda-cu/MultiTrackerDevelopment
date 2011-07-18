/*
 * BandMatrix.hh
 *
 *  Created on: 8/07/2011
 *      Author: jaubry
 */

#ifndef BANDMATRIX_HH_
#define BANDMATRIX_HH_

#include "Definitions.hh"

namespace strandsim
{

template<typename ScalarT, IndexType kl, IndexType ku>
class BandMatrix
{
public:
    explicit BandMatrix(const IndexType rows = 0, const IndexType cols = 0) :
        m_rows(rows), m_cols(cols), m_size((kl + ku + 1) * cols)
    {
        m_data.resize(m_size);
        m_lower.resize(m_rows);
        m_upper.resize(m_rows);

        for (IndexType i = 0; i < m_rows; ++i)
        {
            m_lower[i] = std::max(i - kl, 0);
            m_upper[i] = std::min(i + ku + 1, static_cast<int> (m_cols));
        }

        setZero();
    }

    ~BandMatrix()
    {
    }

    void resize(const IndexType rows, const IndexType cols)
    {
        m_rows = rows;
        m_cols = cols;
        m_size = (kl + ku + 1) * cols;

        m_data.resize(m_size);
        m_lower.resize(m_rows);
        m_upper.resize(m_rows);

        for (IndexType i = 0; i < m_rows; ++i)
        {
            m_lower[i] = std::max(i - kl, 0);
            m_upper[i] = std::min(i + ku + 1, static_cast<int> (m_cols));
        }
    }

    const ScalarT operator()(const IndexType i, const IndexType j) const
    {
        assert(indicesValid(i, j));

        return m_data[(ku + i - j) * m_cols + j];
    }

    ScalarT& operator()(IndexType i, IndexType j)
    {
        assert(indicesValid(i, j));

        return m_data[(ku + i - j) * m_cols + j];
    }

    // Add the local Jacobian to the banded matrix, top left corner at "start" position on the diagonal.
    template<IndexType localSize>
    void localStencilAdd(int start, const Eigen::Matrix<ScalarT, localSize, localSize>& localJ)
    {
        start += m_cols * ku;

        for (int i = 0; i < localSize; ++i)
        {
            for (int j = 0; j < localSize; ++j)
            {
                m_data[start] += localJ(i, j);
                start += 1 - m_cols;
            }
            start += localSize * (m_cols - 1) + m_cols;
        }
    }

    // Add the local Jacobian, with a middle row and a middle column of zeros inserted,
    // to the banded matrix, top left corner at "start" position on the diagonal.
    // Parameter localSize must be an even number.
    template<IndexType localSize>
    void edgeStencilAdd(IndexType start, const Eigen::Matrix<ScalarT, localSize, localSize>& localJ)
    {
        start += m_cols * ku;

        for (int i = 0; i < localSize / 2; ++i)
        {
            for (int j = 0; j < localSize / 2; ++j)
            {
                m_data[start] += localJ(i, j);
                start += 1 - m_cols;
            }
            start += 1 - m_cols;

            for (int j = localSize / 2; j < localSize; ++j)
            {
                m_data[start] += localJ(i, j);
                start += 1 - m_cols;
            }
            start += (localSize + 1) * (m_cols - 1) + m_cols;
        }
        start += m_cols;

        for (int i = localSize / 2; i < localSize; ++i)
        {
            for (int j = 0; j < localSize / 2; ++j)
            {
                m_data[start] += localJ(i, j);
                start += 1 - m_cols;
            }
            start += 1 - m_cols;

            for (int j = localSize / 2; j < localSize; ++j)
            {
                m_data[start] += localJ(i, j);
                start += 1 - m_cols;
            }
            start += (localSize + 1) * (m_cols - 1) + m_cols;
        }
    }

    BandMatrix<ScalarT, kl, ku>& operator*=(const ScalarT multiplier)
    {
        for (std::vector<ScalarT>::iterator i = m_data.begin(); i != m_data.end(); ++i)
            *i *= multiplier;

        return *this;
    }

    void setZero()
    {
        for (int i = 0; i < m_size; ++i)
            m_data[i] = 0;
    }

    void zeroRows(const IntArray& idx, ScalarT diag = 1.0)
    {
        for (int i = 0; i < idx.size(); ++i)
        {
            int r = idx[i];
            assert(r >= 0 && r < MatrixBase::m_rows);
            int lower = m_lower[r];
            int upper = m_upper[r];

            for (int j = lower; j < upper; ++j)
                (*this)(r, j) = 0;

            (*this)(r, r) = diag;
        }
    }

    int computeStartOfCol(int col)
    {
        assert(col >= 0);
        assert(col < m_cols);

        return std::max(col - ku, 0);
    }

    int computeEndOfCol(int col)
    {
        assert(col >= 0);
        assert(col < m_cols);

        return std::min(col + kl, m_rows - 1);
    }

    void zeroCols(const IntArray& idx, ScalarT diag)
    {
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
    }

    void multiply(VecXd& y, ScalarT s, const VecXd& x) const
    {
        assert(y.size() == m_rows);
        assert(x.size() == m_cols);

        for (int i = 0; i < m_rows; ++i)
        {
            int lower = m_lower[i];
            int upper = m_upper[i];
            const ScalarT* val = &m_data[(ku + i - lower) * m_cols + lower];
            ScalarT sum = 0.0;
            for (int j = lower; j < upper; ++j, val += (1 - m_cols))
            {
                sum += (*val) * x[j];
            }
            y[i] += s * sum;
        }
    }

    void print() const
    {
        std::cout << "[";
        for (int i = 0; i < m_rows; ++i)
        {
            for (int j = 0; j < m_cols; ++j)
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

    IndexType rows() const
    {
        return m_rows;
    }

    IndexType cols() const
    {
        return m_cols;
    }

    // NOTE: Assumes not symmetric if lower band is not same size as upper (could be a bunch of zeros in bigger band)
    bool isApproxSymmetric(ScalarT eps) const
    {
        if (m_rows != m_cols || kl != ku)
            return false;

        for (int i = 0; i < m_rows; ++i)
            for (int j = i + 1; j <= i + ku(); ++j)
                if (fabs((*this)(i, j) - (*this)(j, i)) > eps)
                    return false;

        return true;
    }

private:
    bool indicesValid(IndexType r, IndexType c) const
    {
        return r < m_rows && c < m_cols && r - c <= kl && c - r <= ku;
    }

    IndexType m_rows;
    IndexType m_cols;
    IndexType m_size; // Number of non-zero entries
    std::vector<ScalarT> m_data; // Storage of non-zero entries
    std::vector<IndexType> m_lower; // For each row, stores the smallest valid column index
    std::vector<IndexType> m_upper; // For each row, stores the largest valid column index
};

}

#endif /* BANDMATRIX_HH_ */
