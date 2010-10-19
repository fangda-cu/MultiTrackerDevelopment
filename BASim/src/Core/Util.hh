/**
 * \file Util.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 08/30/2009
 */

#ifndef UTIL_HH
#define UTIL_HH

#include "Definitions.hh"

namespace BASim {

void _error(const char* file, const char* function,
                 int line, const char* message) __attribute__((noreturn));

#define BA_ERROR(message) _error(__FILE__, __FUNCTION__, __LINE__, (message))

/** Converts from a string to the given type */
template <class T>
inline void fromString(T& t, const std::string& str)
{
  std::stringstream ss(str);
  ss >> t;
}

/** Converts the given input to a string */
template <class T>
inline std::string toString(const T& t)
{
  std::stringstream ss;
  ss << t;
  return ss.str();
}

/** Checks if the two matrices are approximately equal */
template <typename Derived>
inline bool approxEq(const Eigen::MatrixBase<Derived>& a,
                     const Eigen::MatrixBase<Derived>& b,
                     Scalar tol = SMALL_NUMBER)
{
  return (a - b).isZero(tol);
}

/** Checks if the two Scalars are equal within the given a
    tolerance */
inline bool approxEq(const Scalar& a, const Scalar& b,
                     Scalar tol = SMALL_NUMBER)
{
  return (fabs(a - b) < tol);
}

/** Chekcs if a matrix is approximately symmetric */
template <typename Derived>
inline bool approxSymmetric(const Eigen::MatrixBase<Derived>& A,
                     Scalar tol = SMALL_NUMBER)
{
  for( int i = 0; i < A.rows(); ++i ) for( int j = i+1; j < A.cols(); ++j ) 
    if( !approxEq(A(i,j),A(j,i),tol) ) return false;
  
  return true;
}

template <typename Ref> struct referenceToPointer;
template <typename T> struct referenceToPointer<T&>
{
    typedef T* type;
};

/** Uses dynamic_cast if debugging is turned on and static_cast for
    efficiency when debugging is turned off (by defining NDEBUG. */
template <typename target_ref, typename source>
inline target_ref smart_cast(source& s)
{
  typedef typename referenceToPointer<target_ref>::type target_ptr;

#ifdef NDEBUG

  return static_cast<target_ref>(s);

#else

  target_ptr t = dynamic_cast<target_ptr>(&s);
  assert(t);
  return *t;

#endif

}

namespace Util {

/** Replacement for std::pair that handles alignment automatically
    based on the type of data it is storing. */
template <class T1, class T2> class pair
{
  enum { NeedsToAlign = ((sizeof(T1) + sizeof(T2)) % 16) == 0 };

public:

  typedef T1 first_type;
  typedef T2 second_type;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(NeedsToAlign)

  T1 first;
  T2 second;
  pair() : first(T1()), second(T2()) {}
  pair(const T1& x, const T2& y) : first(x), second(y) {}
  template <class U, class V>
  pair(const pair<U,V>& p) : first(p.first), second(p.second) {}
};

} // namespace Util

} // namespace BASim

#endif // UTIL_HH
