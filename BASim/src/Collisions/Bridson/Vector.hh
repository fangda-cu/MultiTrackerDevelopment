#ifndef __VECTOR_HH__
#define __VECTOR_HH__

#include <cassert>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

//#include "common/Defines.hh"

typedef double Real;

// TODO: Test Vec2
template <class T>
class Vec2
{
public:

   enum { X = 0, Y = 1 };

   Vec2()
   {
      v[0] = v[1] = 0;
   }
   Vec2(const Vec2& c)
   {
      for (uint i = 0; i < 2; ++i) v[i] = c.v[i];
   }
   Vec2(const T& x, const T& y)
   {
      v[0] = x;
      v[1] = y;
   }
   ~Vec2() {}

   Vec2& operator= (const Vec2& c)
   {
      for (uint i = 0; i < 2; ++i) v[i] = c(i);
      return (*this);
   }

   T& operator[] (uint i)
   {
      assert(i < 2);
      return v[i];
   }
   const T& operator[] (uint i) const
   {
      assert(i < 2);
      return v[i];
   }

   uint size() const
   {
      return 2;
   }

   T norm() const
   {
      T t = 0;
      for (uint i = 0; i < 2; ++i) t += v[i]*v[i];
      return sqrt(t);
   }

   T& operator() (uint i)
   {
      assert(i < 2);
      return v[i];
   }
   const T& operator() (uint i) const
   {
      assert(i < 2);
      return v[i];
   }

   T& x()
   {
      return v[X];
   }
   const T& x() const
   {
      return v[X];
   }
   T& y()
   {
      return v[Y];
   }
   const T& y() const
   {
      return v[Y];
   }

   friend Vec2<T> operator+ (const Vec2<T>& a, const Vec2<T>& b)
   {
      Vec2<T> c(a);
      for (uint i = 0; i < 2; ++i) c.v[i]+=b.v[i];
      return c;
   }

   friend Vec2<T> operator- (const Vec2<T>& a, const Vec2<T>& b)
   {
      Vec2<T> c(a);
      for (uint i = 0; i < 2; ++i) c.v[i]-=b.v[i];
      return c;
   }

   friend T norm(const Vec2<T>& a)
   {
      return a.norm();
   }

protected:
   T v[2];
};

template <class T>
Vec2<T> operator- (const Vec2<T>& v)
{
   return Vec2<T>(-v(0), -v(1));
}

template <class T>
Vec2<T> operator- (const Vec2<T>& v1, const Vec2<T>& v2)
{
   return Vec2<T>( v1(0)-v2(0), v1(1)-v2(1) );
}

template <class T>
Vec2<T> operator+ (const Vec2<T>& v1, const Vec2<T>& v2)
{
   return Vec2<T>( v1(0)+v2(0), v1(1)+v2(1) );
}

template <class T>
T dot (const Vec2<T>& v1, const Vec2<T>& v2)
{
   return ( v1(0)*v2(0) + v1(1)*v2(1) );
}

template <class T>
Vec2<T> operator* (const T& t, const Vec2<T>& v)
{
   return Vec2<T>(t*v(0), t*v(1));
}




template <class T>
class Mat2
{
public:
   Mat2()
   {
      for (uint i = 0; i < 4; ++i) v[i] = 0;
   }
   Mat2(const Mat2& c)
   {
      for (uint i = 0; i < 4; ++i) v[i] = c.v[i];
   }
   Mat2(const T& v00, const T& v01, const T& v10, const T& v11)
   {
      v[0] = v00;
      v[1] = v01;
      v[2] = v10;
      v[3] = v11;
   }
   ~Mat2() {}

   template <class G>
   Mat2& operator= (const Mat2<G>& m)
   {
      v[0] = m(0,0);
      v[1] = m(0,1);
      v[2] = m(1,0);
      v[3] = m(1,1);
      return (*this);
   }

   uint nr() const
   {
      return 2;
   }
   uint nc() const
   {
      return 2;
   }

   T& operator() (uint r, uint c)
   {
      assert(r < 2);
      assert(c < 2);
      return v[2*r+c];
   }

   const T& operator() (uint r, uint c) const
   {
      assert(r < 2);
      assert(c < 2);
      return v[2*r+c];
   }

protected:
   T v[4];
};


template <class T, class G>
Vec2<T> operator* (const Vec2<T>& v, const Mat2<G>& m)
{
   return Vec2<T>( v(0)*m(0,0) + v(1)*m(1,0) ,
                   v(0)*m(0,1) + v(1)*m(1,1) );
}

template <class T, class G>
Vec2<T> operator* (const Mat2<G>& m, const Vec2<T>& v)
{
   return Vec2<T>( m(0,0)*v(0) + m(0,1)*v(1) ,
                   m(1,0)*v(0) + m(1,1)*v(1) );
}

template <class T>
Mat2<T> operator* (const Mat2<T>& m1, const Mat2<T>& m2)
{
   return Mat2<T>( m1(0,0)*m2(0,0) + m1(0,1)*m2(1,0) ,
                   m1(0,0)*m2(0,1) + m1(0,1)*m2(1,1) ,
                   m1(1,0)*m2(0,0) + m1(1,1)*m2(1,0) ,
                   m1(1,0)*m2(0,1) + m1(1,1)*m2(1,1) );
}

template <class T>
std::ostream& operator<< (std::ostream& os, const Vec2<T>& a)
{
   os << "[" << a(0) << ", " << a(1) << "]";
   return os;
}

template <class T>
std::ostream& operator<< (std::ostream& os, const Mat2<T>& m)
{
   os << "["
   << m(0,0) << ", " << m(0,1) << "; "
   << m(1,0) << ", " << m(1,1) << "]";
   return os;
}





template <class T> class MMat;


template <class T>
class Vec3
{
public:
   enum { X = 0, Y = 1, Z = 2 };

   Vec3()
   {
      v[0] = v[1] = v[2] = 0;
   }

   Vec3(const Vec3<T>& c)
   {
      for (uint i = 0; i < 3; ++i) v[i] = c.v[i];
   }

   Vec3(const T& a, const T& b, const T& c)
   {
      v[X] = a;
      v[Y] = b;
      v[Z] = c;
   }

   ~Vec3() {}

   Vec3<T>& operator= (const Vec3<T>& c)
   {
      for (uint i = 0; i < 3; ++i) v[i] = c.v[i];
      return *this;
   }

   void scale(const T& a)
   {
      for (uint i = 0; i < 3 ; ++i) v[i] *= a;
   }

   void diagonalScale(const Vec3<T>& a)
   {
      for (uint i = 0; i < 3; ++i) v[i] *= a.v[i];
   }

   void add(const uint& i, const T& val)
   {
      assert( X==i || Y==i || Z==i );
      v[i] += val;
   }

   // this += t*a
   void add(const T& t, const Vec3<T>& a)
   {
      for (uint i = 0; i < 3; ++i) v[i] += t*a.v[i];
   }

   T dot(const Vec3<T>&a) const
   {
      T t = 0;
      for (uint i = 0; i < 3; ++i) t += v[i]*a.v[i];
      return t;
   }

   T norm() const
   {
      T t = 0;
      for (uint i = 0; i < 3; ++i) t += v[i]*v[i];
      return sqrt(t);
   }

   T max() const
   {
       T m = v[0];
       for (uint i = 1; i < 3; ++i) m = max(m, v[i]);
       return m;
   }

   T min() const
   {
       T m = v[0];
       for (uint i = 1; i < 3; ++i) m = min(m, v[i]);
       return m;
   }

   void clear()
   {
      for (uint i = 0; i < 3; ++i) v[i] = 0;
   }

   T* data()
   {
      return v;
   }
   const T* data() const
   {
      return v;
   }

   T& x()
   {
      return v[X];
   }
   const T& x() const
   {
      return v[X];
   }
   T& y()
   {
      return v[Y];
   }
   const T& y() const
   {
      return v[Y];
   }
   T& z()
   {
      return v[Z];
   }
   const T& z() const
   {
      return v[Z];
   }

   T& operator() (const uint& i)
   {
      assert(i < 3);
      return v[i];
   }
   const T& operator() (const uint& i) const
   {
      assert(i < 3);
      return v[i];
   }

   T& operator[] (const uint& i)
   {
      assert(i < 3);
      return v[i];
   }
   const T& operator[] (const uint& i) const
   {
      assert(i < 3);
      return v[i];
   }

   Vec3<T>& operator*= (const T& a)
   {
      for (uint i = 0; i < 3; ++i) v[i] *= a;
      return *this;
   }

   Vec3<T>& operator/= (const T& a)
   {
      if( fabs(a) <= 1.0e-10 )
      {
         std::cout << "This: " << *this << std::endl;
         std::cout << "a: " << a << std::endl;
      }

      assert( fabs(a) > 1.0e-10 );
      for (uint i = 0; i < 3; ++i) v[i] /= a;
      return *this;
   }

   Vec3<T>& operator+= (const Vec3<T>& a)
   {
      for (uint i = 0; i < 3; ++i) v[i] += a.v[i];
      return *this;
   }

   Vec3<T>& operator-= (const Vec3<T>& a)
   {
      for (uint i = 0; i < 3; ++i) v[i] -= a.v[i];
      return *this;
   }

   Vec3<T> operator/ (const T& a) const
   {
      assert( fabs(a) > 1.0e-10 );
      Vec3<T> c(*this);
      c /= a;
      return c;
   }
 
   Vec3<T> operator/(const Vec3<T> &w) const
   {
      Vec3<T> a;
      for (unsigned int i=0; i<3; ++i)
        a[i] = v[i] / w.v[i];

      return a;
   }

   Vec3<T> operator+=(const T &w)
   {
      for (unsigned int i=0; i<3; ++i)
        v[i] += w;

      return *this;
   }

   Vec3<T> operator-=(const T &w)
   {
      for (unsigned int i=0; i<3; ++i)
        v[i] -= w;

      return *this;
   }

   void negate()
   {
      for (uint i = 0; i < 3; ++i) v[i] = -v[i];
   }

   friend void normalize(Vec3<T>& a)
   {
      T norm = a.norm();
      for (uint i = 0; i < 3; ++i) a.v[i] /= norm;
   }

   friend Vec3<T> cross(const Vec3<T>& a, const Vec3<T>& b)
   {
      Vec3<T> c;
      c.v[X] = a.v[Y]*b.v[Z] - a.v[Z]*b.v[Y];
      c.v[Y] = a.v[Z]*b.v[X] - a.v[X]*b.v[Z];
      c.v[Z] = a.v[X]*b.v[Y] - a.v[Y]*b.v[X];
      return c;
   }

   friend T dot(const Vec3<T>& a, const Vec3<T>& b)
   {
      T d = 0;
      for (uint i = 0; i < 3; ++i) d += a.v[i]*b.v[i];
      return d;
   }

   friend T triple(const Vec3<T>& a, const Vec3<T>& b, const Vec3<T>& c)
   {
        return (a.v[0]*(b.v[1]*c.v[2]-b.v[2]*c.v[1])
               +a.v[1]*(b.v[2]*c.v[0]-b.v[0]*c.v[2])
               +a.v[2]*(b.v[0]*c.v[1]-b.v[1]*c.v[0]));
   }

   friend Vec3<T> operator+ (const Vec3<T>& a, const Vec3<T>& b)
   {
      Vec3<T> c(a);
      for (uint i = 0; i < 3; ++i) c.v[i]+=b.v[i];
      return c;
   }

   friend Vec3<T> operator- (const Vec3<T>& a, const Vec3<T>& b)
   {
      Vec3<T> c(a);
      for (uint i = 0; i < 3; ++i) c.v[i]-=b.v[i];
      return c;
   }

   friend Vec3<T> operator* (const T& a, const Vec3<T>& b)
   {
      Vec3<T> c(b);
      c *= a;
      return c;
   }

   friend Vec3<T> operator* (const Vec3<T>& b, const T& a)
   {
      Vec3<T> c(b);
      c *= a;
      return c;
   }

   friend T norm(const Vec3<T>& a)
   {
      return a.norm();
   }

   friend T max(const Vec3<T>& a)
   {
       return a.max();
   }

   friend T min(const Vec3<T>& a)
   {
       return a.min();
   }

   uint size() const
   {
      return 3;
   }

   bool isValid() const
   {
      bool valid = true;
      for (uint i = 0; i < 3; ++i)
      {
         valid = valid && finite(v[i]);
      }
      return valid;
   }

   friend std::ostream& operator<< (std::ostream& os, const Vec3<T>& a)
   {
      os << "[";
      for (uint i = 0; i < 3; ++i)
      {
         //if (fabs(a.v[i]) < 1.0e-10) os << 0;
         //else os << a.v[i];
         os << a.v[i];
         if (i < 2) os << ", ";
      }
      os << "]";
      return os;
   }

   static MMat<T> crossMat(const Vec3<T>& a)
   {
      return MMat<T>(   0, -a[2],  a[1],
                        a[2],     0, -a[0],
                        -a[1],  a[0],     0);
   }

   static void crossMat(MMat<T>& M, const Vec3<T>& a)
   {
      M(0,0) = M(1,1) = M(2,2) = 0;
      M(0,1) = -a(2);
      M(1,0) =  a(2);
      M(0,2) =  a(1);
      M(2,0) = -a(1);
      M(1,2) = -a(0);
      M(2,1) =  a(0);
   }

   static MMat<T> outerProd(const Vec3<T>& a, const Vec3<T>& b)
   {
      return MMat<T>( a(0)*b(0), a(0)*b(1), a(0)*b(2) ,
                      a(1)*b(0), a(1)*b(1), a(1)*b(2) ,
                      a(2)*b(0), a(2)*b(1), a(2)*b(2) );
   }

private:
   T v[3];

};


template <class T>
Vec3<T> operator- (const Vec3<T>& v)
{
   return Vec3<T>(-v(0), -v(1), -v(2));
}

template<class T> 
bool operator<(const Vec3<T> &a, const Vec3<T> &b)
{ 
   for(unsigned int i = 0; i < 3; i++)
      if(a[i] >= b[i])
         return false;
   return true;
}

template<class T> 
bool operator>(const Vec3<T> &a, const Vec3<T> &b)
{ 
   for(unsigned int i = 0; i < 3; i++)
      if(a[i] <= b[i])
         return false;
   return true;
}

template<class T> 
bool operator<=(const Vec3<T> &a, const Vec3<T> &b)
{ 
   for(unsigned int i = 0; i < 3; i++)
      if(a[i] > b[i])
         return false;
   return true;
}

template<class T> 
bool operator>=(const Vec3<T> &a, const Vec3<T> &b)
{ 
   for(unsigned int i = 0; i < 3; i++)
      if(a[i] < b[i])
         return false;
   return true;
}

template<class T>
inline Vec3<int> round(const Vec3<T> &a)
{ 
    Vec3<int> rounded;
    for (unsigned int i=0; i<3; ++i)
    {
        if (a[i] > 0)
            rounded[i] = (a[i]-std::floor(a[i]) < 0.5) ? (int)std::floor(a[i]) : (int)std::ceil(a[i]);
        else
            rounded[i] = (a[i]-std::floor(a[i])<=0.5) ? (int)std::floor(a[i]) : (int)std::ceil(a[i]);
    }

    return rounded; 
}

template<class T>
inline Vec3<int> floor(const Vec3<T> &a)
{ 
    Vec3<int> rounded;
    for (unsigned int i=0; i<3; ++i)
        rounded[i] = (int)std::floor(a[i]);

    return rounded; 
}

template<class T>
inline void minmax(const Vec3<T> &x0, const Vec3<T> &x1, Vec3<T> &xmin, Vec3<T> &xmax)
{
    for(unsigned int i=0; i<3; ++i)
    {
        if (x0[i] < x1[i])
        {
            xmin[i] = x0[i];
            xmax[i] = x1[i];
        }
        else
        {
            xmin[i] = x1[i];
            xmax[i] = x0[i];
        }
    }
}

template<class T>
inline void minmax(const Vec3<T> &x0, const Vec3<T> &x1, const Vec3<T> &x2, const Vec3<T> &x3, Vec3<T> &xmin, Vec3<T> &xmax)
{
    for(unsigned int i=0; i<3; ++i)
    {
        xmin[i] = std::min(x0[i],
                  std::min(x1[i],
                  std::min(x2[i], x3[i])));
        
        xmax[i] = std::max(x0[i],
                  std::max(x1[i],
                  std::max(x2[i], x3[i])));
    }
}

template<class T>
inline void minmax(const Vec3<T> &x0, const Vec3<T> &x1, const Vec3<T> &x2, const Vec3<T> &x3, const Vec3<T> &x4,
                   const Vec3<T> &x5, Vec3<T> &xmin, Vec3<T> &xmax)
{
    for(unsigned int i=0; i<3; ++i)
    {
        xmin[i] = std::min(x0[i],
                  std::min(x1[i],
                  std::min(x2[i],
                  std::min(x3[i],
                  std::min(x4[i], x5[i])))));
        
        xmax[i] = std::max(x0[i],
                  std::max(x1[i],
                  std::max(x2[i],
                  std::max(x3[i],
                  std::max(x4[i], x5[i])))));
    }
}

template<class T>
inline void update_minmax(const Vec3<T> &x, Vec3<T> &xmin, Vec3<T> &xmax)
{
    for(unsigned int i=0; i<3; ++i)
    {
        if (x[i] < xmin[i])
            xmin[i] = x[i];
        else if (x[i] > xmax[i])
            xmax[i] = x[i];
    }
}




// TODO: Test Vec4
template <class T>
class Vec4
{
public:
   enum { X = 0, Y = 1, Z = 2, P = 3 };

   Vec4()
   {
      v[0] = v[1] = v[2] = v[3] = 0;
   }

   Vec4(const Vec4<T>& c)
   {
      for (uint i = 0; i < 4; ++i) v[i] = c.v[i];
   }

   Vec4(T a, T b, T c, T d)
   {
      v[X] = a;
      v[Y] = b;
      v[Z] = c;
      v[P] = d;
   }

   ~Vec4() {}

   Vec4<T>& operator= (const Vec4<T>& c)
   {
      for (uint i = 0; i < 4; ++i) v[i] = c.v[i];
      return *this;
   }

   void scale(T a)
   {
      for (uint i = 0; i < 4 ; ++i) v[i] *= a;
   }

   void diagonalScale(Vec4<T>& a)
   {
      for (uint i = 0; i < 4; ++i) v[i] *= a.v[i];
   }

   void add(uint i, T val)
   {
      v[i] += val;
   }

   // this += t*a
   void add(T t, Vec4<T>& a)
   {
      for (uint i = 0; i < 4; ++i) v[i] += t*a.v[i];
   }

   T dot(const Vec4<T>&a) const
   {
      T t = 0;
      for (uint i = 0; i < 4; ++i) t += v[i]*a.v[i];
      return t;
   }

   T norm() const
   {
      T t = 0;
      for (uint i = 0; i < 4; ++i) t += v[i]*v[i];
      return sqrt(t);
   }

   void clear()
   {
      for (uint i = 0; i < 4; ++i) v[i] = 0;
   }

   T* data()
   {
      return v;
   }
   const T* data() const
   {
      return v;
   }

   T& x()
   {
      return v[X];
   }
   const T& x() const
   {
      return v[X];
   }
   T& y()
   {
      return v[Y];
   }
   const T& y() const
   {
      return v[Y];
   }
   T& z()
   {
      return v[Z];
   }
   const T& z() const
   {
      return v[Z];
   }

   T& p()
   {
      return v[P];
   }
   const T& p() const
   {
      return v[P];
   }

   T& operator() (uint i)
   {
      assert(i < 4);
      return v[i];
   }
   const T& operator() (uint i) const
   {
      assert(i < 4);
      return v[i];
   }

   T& operator[] (uint i)
   {
      assert(i < 4);
      return v[i];
   }
   const T& operator[] (uint i) const
   {
      assert(i < 4);
      return v[i];
   }

   Vec4<T>& operator*= (T a)
   {
      for (uint i = 0; i < 4; ++i) v[i] *= a;
      return *this;
   }

   Vec4<T>& operator/= (T a)
   {
      for (uint i = 0; i < 4; ++i) v[i] /= a;
      return *this;
   }

   Vec4<T>& operator+= (const Vec4<T>& a)
   {
      for (uint i = 0; i < 4; ++i) v[i] += a.v[i];
      return *this;
   }

   Vec3<T>& operator-= (const Vec4<T>& a)
   {
      for (uint i = 0; i < 4; ++i) v[i] -= a.v[i];
      return *this;
   }

   Vec3<T> operator/ (T a) const
   {
      Vec4<T> c(*this);
      c /= a;
      return c;
   }

   void negate()
   {
      for (uint i = 0; i < 4; ++i) v[i] = -v[i];
   }

   friend void normalize(Vec4<T>& a)
   {
      T norm = a.norm();
      for (uint i = 0; i < 4; ++i) a.v[i] /= norm;
   }

   friend T dot(const Vec4<T>& a, const Vec4<T>& b)
   {
      T d = 0;
      for (uint i = 0; i < 4; ++i) d += a.v[i]*b.v[i];
      return d;
   }

   friend Vec4<T> operator+ (const Vec4<T>& a, const Vec4<T>& b)
   {
      Vec4<T> c(a);
      c += b;
      return c;
   }

   friend Vec4<T> operator- (const Vec4<T>& a, const Vec4<T>& b)
   {
      Vec4<T> c(a);
      c -= b;
      return c;
   }

   friend Vec4<T> operator* (T a, const Vec4<T>& b)
   {
      Vec4<T> c(b);
      c *= a;
      return c;
   }

   friend Vec4<T> operator* (const Vec4<T>& b, T a)
   {
      Vec4<T> c(b);
      c *= a;
      return c;
   }

   friend std::ostream& operator<< (std::ostream& os, const Vec4<T>& a)
   {
      os << "[";
      for (uint i = 0; i < 4; ++i)
      {
         if (fabs(a.v[i]) < 1.0e-10) os << 0;
         else os << a.v[i];
         if (i < 3) os << ", ";
      }
      os << "]";
      return os;
   }

   static MMat<T> crossMat(const Vec4<T>& a)
   {
      return MMat<T>(   0, -a[2],  a[1],
                        a[2],     0, -a[0],
                        -a[1],  a[0],     0);
   }

   static void crossMat(MMat<T>& M, const Vec4<T>& a)
   {
      M(0,0) = M(1,1) = M(2,2) = 0;
      M(0,1) = -a(2);
      M(1,0) =  a(2);
      M(0,2) =  a(1);
      M(2,0) = -a(1);
      M(1,2) = -a(0);
      M(2,1) =  a(0);
   }

   static MMat<T> outerProd(const Vec4<T>& a, const Vec4<T>& b)
   {
      return MMat<T>( a(0)*b(0), a(0)*b(1), a(0)*b(2) ,
                      a(1)*b(0), a(1)*b(1), a(1)*b(2) ,
                      a(2)*b(0), a(2)*b(1), a(2)*b(2) );
   }

   uint size() const
   {
      return 4;
   }

protected:
   T v[4];
};





template <class T>
class Mat3
{
public:
   Mat3()
   {
      for (uint i = 0; i < 9; ++i) v[i] = 0;
   }
   Mat3(const Mat3& M)
   {
      for (uint i = 0; i < 9; ++i) v[i] = M.v[i];
   }
   ~Mat3() {}

   Mat3& operator= (const Mat3& M)
   {
      for (uint i = 0; i < 9; ++i) v[i] = M.v[i];
      return *this;
   }

   T& operator() (uint i, uint j)
   {
      assert(i<3 && j<3);
      return v[3*i+j];
   }

   const T& operator() (uint i, uint j) const
   {
      assert(i<3 && j<3);
      return v[3*i+j];
   }

   uint nr() const
   {
      return 3;
   }
   uint nc() const
   {
      return 3;
   }

   void clear()
   {
      for (uint i = 0; i < 9; ++i) v[i] = 0;
   }

protected:
   T v[9];
};



template <class T, class G>
Mat3<T>& operator*= (Mat3<T>& m, const G& t)
{
   for (uint i = 0; i < m.nr(); ++i)
      for (uint j = 0; j < m.nc(); ++j)
         m(i,j) *= t;
   return m;
}

template <class T>
Mat3<T> operator* (const T& t, const Mat3<T>& m1)
{
   Mat3<T> m2(m1);
   m2 *= t;
   return m2;
}

template <class T>
Mat3<T> operator+ (const Mat3<T>& m1, const Mat3<T>& m2)
{
   Mat3<T> m3;
   for (uint i = 0; i < m1.nr(); ++i)
      for (uint j = 0; j < m1.nc(); ++j)
         m3(i,j) = m1(i,j) + m2(i,j);
   return m3;
}

template <class T>
Mat3<T>& operator+= (Mat3<T>& m1, const Mat3<T>& m2)
{
   for (uint i = 0; i < 3; ++i)
      for (uint j = 0; j < 3; ++j)
         m1(i,j) += m2(i,j);
   return m1;
}

template <class T>
Mat3<T>& operator-= (Mat3<T>& m1, const Mat3<T>& m2)
{
   for (uint i = 0; i < 3; ++i)
      for (uint j = 0; j < 3; ++j)
         m1(i,j) -= m2(i,j);
   return m1;
}

template <class M, class V>
V dot(const M& m, const V& v)
{
   assert (v.size() == m.nr());

   V r;

   for (uint i = 0; i < m.nr(); ++i)
      for (uint j = 0; j < m.nc(); ++j)
         r[j] += m(i,j) * v(i);

   return r;
}

template <class V, class T>
void crossMat(Mat3<T>& m, const V& v)
{
   m(0,0) = m(1,1) = m(2,2) = 0;
   m(0,1) = -v(2);
   m(1,0) =  v(2);
   m(0,2) =  v(1);
   m(2,0) = -v(1);
   m(1,2) = -v(0);
   m(2,1) =  v(0);
}

template <class T>
Mat3<T> crossMat(const Vec3<T>& v)
{
   Mat3<T> m;
   crossMat(m, v);
   return m;
}

template <class T>
Mat3<T> outerProd(const Vec3<T>& v1, const Vec3<T>& v2)
{
   Mat3<T> m;
   for (uint i = 0; i < 3; ++i)
   {
      for (uint j = 0; j < 3; ++j)
      {
         m(i,j) += v1(i)*v2(j);
      }
   }
   return m;
}

template <class V, class T>
V operator* (const Mat3<T>& m, const V& v1)
{
   V v2;
   for (uint i = 0; i < 3; ++i)
      for (uint j = 0; j < 3; ++j)
         v2(i) += m(i,j)*v1(j);
   return v2;
}


template <class T>
class Vector
{
public:
   enum { X = 0, Y = 1, Z = 2 };

   Vector(uint n = 3) : _n(n)
   {
      v = new T[_n];
      for (uint i = 0; i < _n; ++i) v[i] = 0;
   }

   Vector(const Vector<T>& c) : _n(c._n)
   {
      v = new T[_n];
      for (uint i = 0; i < _n; ++i) v[i] = c.v[i];
   }

   Vector(T a, T b, T c) : _n(3)
   {
      v = new T[_n];
      v[X] = a;
      v[Y] = b;
      v[Z] = c;
   }

   ~Vector()
   {
      delete [] v;
   }

   Vector<T>& operator= (const Vector<T>& c)
   {
      assert(_n == c._n);
      for (uint i = 0; i < _n; ++i) v[i] = c.v[i];
      return *this;
   }

   void scale(T a)
   {
      for (uint i = 0; i < _n ; ++i) v[i] *= a;
   }

   void diagonalScale(Vector<T>& a)
   {
      assert(_n == a._n);

      for (uint i = 0; i < _n; ++i) v[i] *= a.v[i];
   }

   void add(uint i, T val)
   {
      assert(i < _n);
      v[i] += val;
   }

   // this += t*a
   void add(T t, Vector<T>& a)
   {
      assert(_n == a._n);
      for (uint i = 0; i < _n; ++i) v[i] += t*a.v[i];
   }

   T dot(const Vector<T>&a) const
   {
      assert(_n == a._n);
      T t = 0;
      for (uint i = 0; i < _n; ++i) t += v[i]*a.v[i];
      return t;
   }

   T norm() const
   {
      T t = 0;
      for (uint i = 0; i < _n; ++i) t += v[i]*v[i];
      return sqrt(t);
   }

   void clear()
   {
      for (uint i = 0; i < _n; ++i) v[i] = 0;
   }

   T* data()
   {
      return v;
   }
   const T* data() const
   {
      return v;
   }

   T& x()
   {
      return v[X];
   }
   const T& x() const
   {
      return v[X];
   }
   T& y()
   {
      return v[Y];
   }
   const T& y() const
   {
      return v[Y];
   }
   T& z()
   {
      return v[Z];
   }
   const T& z() const
   {
      return v[Z];
   }

   T& operator() (uint i)
   {
      assert(i < _n);
      return v[i];
   }
   const T& operator() (uint i) const
   {
      assert(i < _n);
      return v[i];
   }

   T& operator[] (uint i)
   {
      /*
      if( i >= _n)
      {
         std::cout << i << ", " << _n << std::endl;
      }
      */
      assert(i < _n);
      return v[i];
   }

   const T& operator[] (uint i) const
   {
      assert(i < _n);
      return v[i];
   }

   Vector<T>& operator*= (T a)
   {
      for (uint i = 0; i < _n; ++i) v[i] *= a;
      return *this;
   }

   Vector<T>& operator/= (T a)
   {
      for (uint i = 0; i < _n; ++i) v[i] /= a;
      return *this;
   }

   Vector<T>& operator+= (const Vector<T>& a)
   {
      assert(_n == a._n);
      for (uint i = 0; i < _n; ++i) v[i] += a.v[i];
      return *this;
   }

   Vector<T>& operator-= (const Vector<T>& a)
   {
      assert(_n == a._n);
      for (uint i = 0; i < _n; ++i) v[i] -= a.v[i];
      return *this;
   }

   Vector<T> operator/ (T a) const
   {
      Vector<T> c(*this);
      c /= a;
      return c;
   }

   void negate()
   {
      for (uint i = 0; i < _n; ++i) v[i] = -v[i];
   }

   friend void normalize(Vector<T>& a)
   {
      T norm = a.norm();
      for (uint i = 0; i < a._n; ++i) a.v[i] /= norm;
   }

   friend Vector<T> cross(const Vector<T>& a, const Vector<T>& b)
   {
      assert(a._n == 3);
      assert(b._n == 3);
      Vector<T> c;
      c.v[X] = a.v[Y]*b.v[Z] - a.v[Z]*b.v[Y];
      c.v[Y] = a.v[Z]*b.v[X] - a.v[X]*b.v[Z];
      c.v[Z] = a.v[X]*b.v[Y] - a.v[Y]*b.v[X];
      return c;
   }

   friend T dot(const Vector<T>& a, const Vector<T>& b)
   {
      assert(a._n == b._n);
      T d = 0;
      for (uint i = 0; i < a._n; ++i) d += a.v[i]*b.v[i];
      return d;
   }

   friend Vector<T> operator+ (const Vector<T>& a, const Vector<T>& b)
   {
      assert(a._n == b._n);
      Vector<T> c(a);
      c += b;
      return c;
   }

   friend Vector<T> operator- (const Vector<T>& a, const Vector<T>& b)
   {
      assert(a._n == b._n);
      Vector<T> c(a);
      c -= b;
      return c;
   }

   friend Vector<T> operator* (T a, const Vector<T>& b)
   {
      Vector<T> c(b);
      c *= a;
      return c;
   }

   friend std::ostream& operator<< (std::ostream& os, const Vector<T>& a)
   {
      os << "[";
      for (uint i = 0; i < a._n-1; ++i) os << a.v[i] << ", ";
      os << a.v[a._n-1] << "]";
      return os;
   }

   static MMat<T> crossMat(const Vector<T>& a)
   {
      MMat<T> m(   0, -a[2],  a[1],
                   a[2],     0, -a[0],
                   -a[1],  a[0],     0);
      return m;
   }

   static MMat<T> outerProd(const Vector<T>& a, const Vector<T>& b)
   {
      assert(a._n == b._n);
      MMat<T> M;
      for (uint i = 0; i < a._n; ++i)
      {
         for (uint j = 0; j < b._n; ++j)
         {
            M(i,j) += a[i]*b[j];
         }
      }
      return M;
   }

   uint size() const
   {
      return _n;
   }

protected:
   uint _n;
   T* v;
};


typedef Vec3<Real>          Vec3r;
//typedef Vec3<double>        Vec3d;
typedef Vec3<float>         Vec3f;
typedef Vec3<int>           Vec3i;
typedef Vec3<unsigned int>  Vec3ui;


typedef std::vector<Vec3<Real> > Positions;
typedef std::vector<Vec3<Real> > Velocities;
typedef std::vector<unsigned int>      Indices;
typedef std::vector<Vec3ui> Vec3Indices;

#endif // __VECTOR_HH__
