/*
 * BoundingBox.hh
 *
 *  Created on: 17/03/2011
 *      Author: Jean-Marie Aubry <jaubry@wetafx.co.nz>
 */

#ifndef BOUNDINGBOX_HH_
#define BOUNDINGBOX_HH_

#include <limits>
#include <algorithm>
#include <iostream>
#include "../Core/Definitions.hh"

namespace BASim
{

typedef unsigned int uint;

template<class T>
class Point
{
public:
    enum
    {
        X = 0, Y = 1, Z = 2
    };

    Point()
    {
        v[0] = v[1] = v[2] = 0;
    }

    Point(const Point<T>& c)
    {
        for (uint i = 0; i < 3; ++i)
            v[i] = c.v[i];
    }

    Point(const T& a, const T& b, const T& c)
    {
        v[X] = a;
        v[Y] = b;
        v[Z] = c;
    }

    explicit Point(const Vec3d& eigen_point)
    {
        v[0] = eigen_point[0];
        v[1] = eigen_point[1];
        v[2] = eigen_point[2];
    }

    ~Point()
    {
    }

    Point<T>& operator=(const Point<T>& c)
    {
        for (uint i = 0; i < 3; ++i)
            v[i] = c.v[i];
        return *this;
    }

    void scale(const T& a)
    {
        for (uint i = 0; i < 3; ++i)
            v[i] *= a;
    }

    void diagonalScale(const Point<T>& a)
    {
        for (uint i = 0; i < 3; ++i)
            v[i] *= a.v[i];
    }

    void add(const uint& i, const T& val)
    {
        assert(X == i || Y == i || Z == i);
        v[i] += val;
    }

    // this += t*a
    void add(const T& t, const Point<T>& a)
    {
        for (uint i = 0; i < 3; ++i)
            v[i] += t * a.v[i];
    }

    T dot(const Point<T>&a) const
    {
        T t = 0;
        for (uint i = 0; i < 3; ++i)
            t += v[i] * a.v[i];
        return t;
    }

    T norm() const
    {
        T t = 0;
        for (uint i = 0; i < 3; ++i)
            t += v[i] * v[i];
        return sqrt(t);
    }

    T max() const
    {
        T m = v[0];
        for (uint i = 1; i < 3; ++i)
            m = max(m, v[i]);
        return m;
    }

    T min() const
    {
        T m = v[0];
        for (uint i = 1; i < 3; ++i)
            m = min(m, v[i]);
        return m;
    }

    void clear()
    {
        for (uint i = 0; i < 3; ++i)
            v[i] = 0;
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

    T& operator()(const uint& i)
    {
        assert(i < 3);
        return v[i];
    }
    const T& operator()(const uint& i) const
    {
        assert(i < 3);
        return v[i];
    }

    T& operator[](const uint& i)
    {
        return v[i];
    }
    const T& operator[](const uint& i) const
    {
        return v[i];
    }

    Point<T>& operator*=(const T& a)
    {
        for (uint i = 0; i < 3; ++i)
            v[i] *= a;
        return *this;
    }

    Point<T>& operator/=(const T& a)
    {
        if (fabs(a) <= 1.0e-10)
        {
            std::cout << "This: " << *this << std::endl;
            std::cout << "a: " << a << std::endl;
        }

        assert(fabs(a) > 1.0e-10);
        for (uint i = 0; i < 3; ++i)
            v[i] /= a;
        return *this;
    }

    Point<T>& operator+=(const Point<T>& a)
    {
        for (uint i = 0; i < 3; ++i)
            v[i] += a.v[i];
        return *this;
    }

    Point<T>& operator-=(const Point<T>& a)
    {
        for (uint i = 0; i < 3; ++i)
            v[i] -= a.v[i];
        return *this;
    }

    Point<T> operator/(const T& a) const
    {
        assert(fabs(a) > 1.0e-10);
        Point<T> c(*this);
        c /= a;
        return c;
    }

    Point<T> operator/(const Point<T> &w) const
    {
        Point<T> a;
        for (unsigned int i = 0; i < 3; ++i)
            a[i] = v[i] / w.v[i];

        return a;
    }

    Point<T> operator+=(const T &w)
    {
        for (unsigned int i = 0; i < 3; ++i)
            v[i] += w;

        return *this;
    }

    Point<T> operator-=(const T &w)
    {
        for (unsigned int i = 0; i < 3; ++i)
            v[i] -= w;

        return *this;
    }

    void negate()
    {
        for (uint i = 0; i < 3; ++i)
            v[i] = -v[i];
    }

    friend void normalize(Point<T>& a)
    {
        T norm = a.norm();
        for (uint i = 0; i < 3; ++i)
            a.v[i] /= norm;
    }

    friend Point<T> cross(const Point<T>& a, const Point<T>& b)
    {
        Point<T> c;
        c.v[X] = a.v[Y] * b.v[Z] - a.v[Z] * b.v[Y];
        c.v[Y] = a.v[Z] * b.v[X] - a.v[X] * b.v[Z];
        c.v[Z] = a.v[X] * b.v[Y] - a.v[Y] * b.v[X];
        return c;
    }

    friend T dot(const Point<T>& a, const Point<T>& b)
    {
        T d = 0;
        for (uint i = 0; i < 3; ++i)
            d += a.v[i] * b.v[i];
        return d;
    }

    friend T triple(const Point<T>& a, const Point<T>& b, const Point<T>& c)
    {
        return (a.v[0] * (b.v[1] * c.v[2] - b.v[2] * c.v[1]) + a.v[1] * (b.v[2] * c.v[0] - b.v[0] * c.v[2]) + a.v[2] * (b.v[0]
                * c.v[1] - b.v[1] * c.v[0]));
    }

    friend Point<T> operator+(const Point<T>& a, const Point<T>& b)
    {
        Point<T> c(a);
        for (uint i = 0; i < 3; ++i)
            c.v[i] += b.v[i];
        return c;
    }

    friend Point<T> operator-(const Point<T>& a, const Point<T>& b)
    {
        Point<T> c(a);
        for (uint i = 0; i < 3; ++i)
            c.v[i] -= b.v[i];
        return c;
    }

    friend Point<T> operator*(const T& a, const Point<T>& b)
    {
        Point<T> c(b);
        c *= a;
        return c;
    }

    friend Point<T> operator*(const Point<T>& b, const T& a)
    {
        Point<T> c(b);
        c *= a;
        return c;
    }

    friend T norm(const Point<T>& a)
    {
        return a.norm();
    }

    friend T max(const Point<T>& a)
    {
        return a.max();
    }

    friend T min(const Point<T>& a)
    {
        return a.min();
    }

private:
    T v[3];

};

template<typename T>
Point<T> Min(const Point<T>& a, const Point<T>& b)
{
    return Point<T> (std::min<T>(a.x(), b.x()), std::min<T>(a.y(), b.y()), std::min<T>(a.z(), b.z()));
}

template<typename T>
Point<T> Max(const Point<T>& a, const Point<T>& b)
{
    return Point<T> (std::max<T>(a.x(), b.x()), std::max<T>(a.y(), b.y()), std::max<T>(a.z(), b.z()));
}

template<typename T>
struct BoundingBox
{
    typedef Point<T> PointType;
    PointType min;
    PointType max;

    BoundingBox() :
        min(std::numeric_limits<T>::max(), std::numeric_limits<T>::max(), std::numeric_limits<T>::max()),
                max(-std::numeric_limits<T>::max(), -std::numeric_limits<T>::max(), -std::numeric_limits<T>::max())
    {
    }

    BoundingBox(T minx, T miny, T minz, T maxx, T maxy, T maxz) :
        min(minx, miny, minz), max(maxx, maxy, maxz)
    {
    }

    bool IsValid() const
    {
        return max.x() >= min.x() && max.y() >= min.y() && max.z() >= min.z();
    }

    void Reset()
    {
        *this = BoundingBox();
    }

    // Expand the box to contain the given point
    inline void Insert(const PointType& p)
    {
        min[0] = std::min(min[0], p[0]);
        min[1] = std::min(min[1], p[1]);
        min[2] = std::min(min[2], p[2]);
        max[0] = std::max(max[0], p[0]);
        max[1] = std::max(max[1], p[1]);
        max[2] = std::max(max[2], p[2]);
    }

    // Expand the box to contain the given sphere
    inline void Insert(const PointType& p, const T& radius)
    {
        min[0] = std::min(min[0], p[0] - radius);
        min[1] = std::min(min[1], p[1] - radius);
        min[2] = std::min(min[2], p[2] - radius);
        max[0] = std::max(max[0], p[0] + radius);
        max[1] = std::max(max[1], p[1] + radius);
        max[2] = std::max(max[2], p[2] + radius);
    }

    // Expand the box to contain the given box
    inline void Insert(const BoundingBox& box)
    {
        Insert(box.min);
        Insert(box.max);
    }

    T Volume() const
    {
        return (max.x() - min.x()) * (max.y() - min.y()) * (max.z() - min.z());
    }

};

template<typename T>
bool Intersect(const BoundingBox<T>& bbox_a, const BoundingBox<T>& bbox_b)
{
    // Attempt to find a separating hyperplane
    if ((bbox_a.max[0] < bbox_b.min[0]) || (bbox_b.max[0] < bbox_a.min[0]) || (bbox_a.max[1] < bbox_b.min[1]) || (bbox_b.max[1]
            < bbox_a.min[1]) || (bbox_a.max[2] < bbox_b.min[2]) || (bbox_b.max[2] < bbox_a.min[2]))
        return false;

    // Otherwise the bounding volumes overlap!
    return true;
}

}

#endif /* BOUNDINGBOX_HH_ */
