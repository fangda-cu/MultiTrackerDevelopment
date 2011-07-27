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
    explicit BandMatrix( const IndexType rows = 0, const IndexType cols = 0 ) :
        m_rows( rows ), m_cols( cols ), m_size( ( kl + ku + 1 ) * cols )
    {
        m_data.resize( m_size );
    }

    ~BandMatrix()
    {
    }

    void resize( const IndexType rows, const IndexType cols )
    {
        m_rows = rows;
        m_cols = cols;
        m_size = ( kl + ku + 1 ) * cols;

        m_data.resize( m_size );
    }

    const ScalarT operator()( const IndexType i, const IndexType j ) const
    {
        if ( !indicesValid( i, j ) )
            return 0.0;

        return m_data[( ku + i - j ) * m_cols + j];
    }

    ScalarT& operator()( IndexType i, IndexType j )
    {
        assert( indicesValid( i, j ) );

        return m_data[( ku + i - j ) * m_cols + j];
    }

    // Add lambda * Identity
    void addConstantDiagonal( const Scalar lambda )
    {
        assert( m_cols == m_rows );

        for ( int i = 0; i < m_cols; i++ )
            m_data[ku * m_cols + i] += lambda;
    }

    template<int nfixed>
    void fixFirstDOFs()
    {
        assert( nfixed <= m_rows && nfixed <= m_cols );

        for ( int i = 0; i < nfixed; i++ )
        {
            for ( int j = std::max( 0, i - kl ); j <= std::min( ( int ) m_cols - 1, i + ku ); ++j )
                ( *this )( i, j ) = 0.0;
            for ( int k = std::max( 0, i - ku ); k <= std::min( ( int ) m_rows - 1, i + kl ); k++ )
                ( *this )( k, i ) = 0.0;
            ( *this )( i, i ) = 1.0;

        }
    }

    // Add the local Jacobian to the banded matrix, top left corner at "start" position on the diagonal.
    template<IndexType localSize>
    void localStencilAdd( int start, const Eigen::Matrix<ScalarT, localSize, localSize>& localJ )
    {
        start += m_cols * ku;

        for ( int i = 0; i < localSize; ++i )
        {
            for ( int j = 0; j < localSize; ++j )
            {
                m_data[start] += localJ( i, j );
                start += 1 - m_cols;
            }
            start += localSize * ( m_cols - 1 ) + m_cols;
        }
    }

    // Add the local Jacobian, with a middle row and a middle column of zeros inserted,
    // to the banded matrix, top left corner at "start" position on the diagonal.
    // Parameter localSize must be an even number.
    template<IndexType localSize>
    void edgeStencilAdd( IndexType start,
            const Eigen::Matrix<ScalarT, localSize, localSize>& localJ )
    {
        start += m_cols * ku;

        for ( int i = 0; i < localSize / 2; ++i )
        {
            for ( int j = 0; j < localSize / 2; ++j )
            {
                m_data[start] += localJ( i, j );
                start += 1 - m_cols;
            }
            start += 1 - m_cols;

            for ( int j = localSize / 2; j < localSize; ++j )
            {
                m_data[start] += localJ( i, j );
                start += 1 - m_cols;
            }
            start += ( localSize + 1 ) * ( m_cols - 1 ) + m_cols;
        }
        start += m_cols;

        for ( int i = localSize / 2; i < localSize; ++i )
        {
            for ( int j = 0; j < localSize / 2; ++j )
            {
                m_data[start] += localJ( i, j );
                start += 1 - m_cols;
            }
            start += 1 - m_cols;

            for ( int j = localSize / 2; j < localSize; ++j )
            {
                m_data[start] += localJ( i, j );
                start += 1 - m_cols;
            }
            start += ( localSize + 1 ) * ( m_cols - 1 ) + m_cols;
        }
    }

    BandMatrix<ScalarT, kl, ku>& operator*=( const ScalarT multiplier )
    {
        for ( typename std::vector<ScalarT>::iterator i = m_data.begin(); i != m_data.end(); ++i )
            *i *= multiplier;

        return *this;
    }

    void setZero()
    {
        for ( int i = 0; i < m_size; ++i )
            m_data[i] = 0;
    }

    void multiply( VecXd& y, ScalarT s, const VecXd& x ) const
    {
        assert( y.size() == m_rows );
        assert( x.size() == m_cols );

        for ( int i = 0; i < m_rows; ++i )
        {
            int lower = m_lower[i];
            int upper = m_upper[i];
            const ScalarT* val = &m_data[( ku + i - lower ) * m_cols + lower];
            ScalarT sum = 0.0;
            for ( int j = lower; j < upper; ++j, val += ( 1 - m_cols ) )
            {
                sum += ( *val ) * x[j];
            }
            y[i] += s * sum;
        }
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
    bool isApproxSymmetric( ScalarT eps ) const
    {
        if ( m_rows != m_cols || kl != ku )
            return false;

        for ( int i = 0; i < m_rows; ++i )
            for ( int j = i + 1; j <= i + ku(); ++j )
                if ( fabs( ( *this )( i, j ) - ( *this )( j, i ) ) > eps )
                    return false;

        return true;
    }

    const std::vector<ScalarT>& getData() const
    {
        return m_data;
    }

private:
    bool indicesValid( IndexType r, IndexType c ) const
    {
        return r < m_rows && c < m_cols && r - c <= kl && c - r <= ku;
    }

    IndexType m_rows;
    IndexType m_cols;
    IndexType m_size; // Number of non-zero entries
    std::vector<ScalarT> m_data; // Storage of non-zero entries
};

template<typename ScalarT, IndexType kl, IndexType ku>
std::ostream& operator<<( std::ostream& os, const BandMatrix<ScalarT, kl, ku>& M )
{
    os << '{';
    for ( int i = 0; i < M.rows() - 1; i++ )
    {
        os << '{';
        for ( int j = 0; j < M.cols() - 1; j++ )
            os << M( i, j ) << ", ";
        os << M( i, M.cols() - 1 ) << "}, ";
    }
    os << '{';
    for ( int j = 0; j < M.cols() - 1; j++ )
        os << M( M.rows() - 1, j ) << ", ";
    os << M( M.rows() - 1, M.cols() - 1 ) << "}}";

    return os;
}

}

#endif /* BANDMATRIX_HH_ */
