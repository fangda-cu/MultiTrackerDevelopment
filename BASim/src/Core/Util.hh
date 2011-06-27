/**
 * \file Util.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 08/30/2009
 */

#ifndef UTIL_HH
#define UTIL_HH

#include "Definitions.hh"

namespace BASim
{

static const Scalar SMALL_NUMBER = 1e-12; // std::numeric_limits<Scalar>::epsilon();

void _error(const char* file, const char* function, int line, const char* message) __attribute__((noreturn));

#define BA_ERROR(message) _error(__FILE__, __FUNCTION__, __LINE__, (message))

inline std::ostream& operator<<(std::ostream& os, const std::vector<VecXd>& normals)
{
    for (std::vector<VecXd>::const_iterator normal = normals.begin(); normal != normals.end(); ++normal)
        os << *normal << '\n';

    return os;
}

/** Converts from a string to the given type */
template<class T>
inline void fromString(T& t, const std::string& str)
{
    std::stringstream ss(str);
    ss >> t;
}

/** Converts the given input to a string */
template<class T>
inline std::string toString(const T& t)
{
    std::stringstream ss;
    ss << t;
    return ss.str();
}

/** Checks if the two matrices are approximately equal */
template<typename Derived>
inline bool approxEq(const Eigen::MatrixBase<Derived>& a, const Eigen::MatrixBase<Derived>& b, Scalar tol = SMALL_NUMBER)
{
    return (a - b).isZero(tol);
}

/** Checks if the two Scalars are equal within the given a
 tolerance */
inline bool approxEq(const Scalar& a, const Scalar& b, Scalar tol = SMALL_NUMBER)
{
    return (fabs(a - b) < tol);
}

/** Checks if a matrix is approximately symmetric */
template<typename MatrixT>
inline bool approxSymmetric(const MatrixT& A, Scalar tol = SMALL_NUMBER)
{
    for (int i = 0; i < A.rows(); ++i)
        for (int j = i + 1; j < A.cols(); ++j)
            if (!approxEq(A(i, j), A(j, i), tol))
            {
                std::cerr << "approxSymmetric failing by " << fabs(A(i, j) - A(j, i)) << '\n';
                return false;
            }
    return true;
}

template<typename MatrixT>
Scalar maxAsymmetry(const MatrixT& A)
{
    Scalar maxDev = 0.0;
    for (int i = 0; i < A.rows(); ++i)
        for (int j = i + 1; j < A.cols(); ++j)
            maxDev = max(maxDev, fabs(A(i, j) - A(j, i)));

    return maxDev;
}

template<typename MatrixT>
inline bool isSymmetric(const MatrixT& A)
{
    for (int i = 0; i < A.rows(); ++i)
        for (int j = i + 1; j < A.cols(); ++j)
            if (A(i, j) != A(j, i))
            {
                std::cerr << "isSymmetric failing by " << fabs(A(i, j) - A(j, i)) << '\n';
                return false;
            }
    return true;
}

// Computes u B v^T, assuming B is symmetric.
inline Scalar BProduct(const Mat2d& B, const Vec2d& u, const Vec2d& v)
{
    assert(isSymmetric(B));

    return B(0, 0) * u[0] * v[0] + B(0, 1) * (u[0] * v[1] + u[1] * v[0]) + B(1, 1) * u[1] * v[1];
}

// Computes Q B Q^T, assuming B is symmetric. The result is then (exactly) symmetric.
template<int n>
inline void symBProduct(Eigen::Matrix<Scalar, n, n>& result, const Mat2d& B, const MatXd& Q)
{
    assert(isSymmetric(B));
    assert(Q.rows() == n);
    assert(Q.cols() == 2);

    for (int i = 0; i < n; ++i)
    {
        result(i, i) = BProduct(B, Q.row(i), Q.row(i));
        for (int j = 0; j < i; ++j)
            result(i, j) = result(j, i) = BProduct(B, Q.row(i), Q.row(j));
    }
}

/** Uses dynamic_cast if debugging is turned on and static_cast for
 efficiency when debugging is turned off (by defining NDEBUG. */
template<typename target_ptr, typename source>
inline target_ptr smart_cast(source* s)
{
#ifdef NDEBUG

    return static_cast<target_ptr>(s);

#else

    target_ptr t = dynamic_cast<target_ptr> (s);
    assert(t != NULL);
    return t;

#endif

}

/** Uses dynamic_cast if debugging is turned on and static_cast for
 efficiency when debugging is turned off (by defining NDEBUG. */
template<typename target_ref, typename source>
inline target_ref smart_cast(source& s)
{

#ifdef NDEBUG

    return static_cast<target_ref>(s);

#else

    try
    {
        target_ref t = dynamic_cast<target_ref> (s);
        return t;
    } catch (...)
    {
        assert(!"bad cast");
    }

    return dynamic_cast<target_ref> (s);

#endif

}

namespace Util
{

/** Replacement for std::pair that handles alignment automatically
 based on the type of data it is storing. */
template<class T1, class T2> class pair
{
    enum
    {
        NeedsToAlign = ((sizeof(T1) + sizeof(T2)) % 16) == 0
    };

public:

    typedef T1 first_type;
    typedef T2 second_type;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(NeedsToAlign)

    T1 first;
    T2 second;
    pair() :
        first(T1()), second(T2())
    {
    }
    pair(const T1& x, const T2& y) :
        first(x), second(y)
    {
    }
    template<class U, class V>
    pair(const pair<U, V>& p) :
        first(p.first), second(p.second)
    {
    }
};

} // namespace Util

} // namespace BASim

#endif // UTIL_HH
