#include "CollisionUtils.hh"
#include "../Util/TextLog.hh"

namespace BASim
{

// Adapted from Christer Ericson, "Real Time Collision Detection"
Vec3d ClosestPtPointTriangle(const Vec3d& p, const Vec3d& a, const Vec3d& b, const Vec3d& c)
{
    Vec3d result;

    // Check if P in vertex region outside A
    const Vec3d ab = b - a;
    const Vec3d ac = c - a;
    const Vec3d ap = p - a;
    double d1 = ab.dot(ap);
    double d2 = ac.dot(ap);
    if (d1 <= 0.0 && d2 <= 0.0)
    {
        result = a; // barycentric coordinates (1,0,0)
       // TraceStream(g_log, "") << "ClosestPtPointTriangle: CASE 1 p = " << p << " a = " << a << " b = " << b << " c = " << c << " result = " << result << '\n';
        return result;
    }

    // Check if P in vertex region outside B
    const Vec3d bp = p - b;
    double d3 = ab.dot(bp);
    double d4 = ac.dot(bp);
    if (d3 >= 0.0 && d4 <= d3)
    {
        result = b; // barycentric coordinates (0,1,0)
       // TraceStream(g_log, "") << "ClosestPtPointTriangle: CASE 2 p = " << p << " a = " << a << " b = " << b << " c = " << c << " result = " << result << '\n';
        return result;
    }

    // Check if P in edge region of AB, if so return projection of P onto AB
    double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
    {
        double v = d1 / (d1 - d3);
        result = a + v * ab; // barycentric coordinates (1-v,v,0)
       // TraceStream(g_log, "") << "ClosestPtPointTriangle: CASE 3 p = " << p << " a = " << a << " b = " << b << " c = " << c << " result = " << result << '\n';
        return result;
    }

    // Check if P in vertex region outside C
    const Vec3d cp = p - c;
    double d5 = ab.dot(cp);
    double d6 = ac.dot(cp);
    if (d6 >= 0.0 && d5 <= d6)
    {
        result = c; // barycentric coordinates (0,0,1)
       // TraceStream(g_log, "") << "ClosestPtPointTriangle: CASE 4 p = " << p << " a = " << a << " b = " << b << " c = " << c << " result = " << result << '\n';
        return result;
    }

    // Check if P in edge region of AC, if so return projection of P onto AC
    double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
    {
        double w = d2 / (d2 - d6);
        result = a + w * ac; // barycentric coordinates (1-w,0,w)
       // TraceStream(g_log, "") << "ClosestPtPointTriangle: CASE 5 p = " << p << " a = " << a << " b = " << b << " c = " << c << " result = " << result << '\n';
	return result;
    }

    // Check if P in edge region of BC, if so return projection of P onto BC
    double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0)
    {
        double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        result = b + w * (c - b); // barycentric coordinates (0,1-w,w)
       // TraceStream(g_log, "") << "ClosestPtPointTriangle: CASE 6 p = " << p << " a = " << a << " b = " << b << " c = " << c << " result = " << result << '\n';
	return result;
    }

    // P inside face region. Compute Q through its barycentric coordinates (u,v,w)
    double denom = 1.0 / (va + vb + vc);
    double v = vb * denom;
    double w = vc * denom;
    result = a + ab * v + ac * w; // = u*a + v*b + w*c, u = va * denom = 1.0f - v - w

   // TraceStream(g_log, "") << "ClosestPtPointTriangle: CASE 7 p = " << p << " a = " << a << " b = " << b << " c = " << c << " result = " << result << '\n';

    return result;
}

// Adapted from Christer Ericson, "Real Time Collision Detection"
// Computes closest points C1 and C2 of S1(s)=P1+s*(Q1-P1) and
// S2(t)=P2+t*(Q2-P2), returning s and t. Function result is squared
// distance between between S1(s) and S2(t).
// TODO: Explore behavior in degenerate case more closely.
double ClosestPtSegmentSegment(const Vec3d& p1, const Vec3d& q1, const Vec3d& p2, const Vec3d& q2, double& s, double& t,
        Vec3d& c1, Vec3d& c2)
{
    double EPSILON = 1.0e-12;

    Vec3d d1 = q1 - p1; // Direction vector of segment S1
    Vec3d d2 = q2 - p2; // Direction vector of segment S2
    Vec3d r = p1 - p2;
    double a = d1.dot(d1); // Squared length of segment S1, always nonnegative
    double e = d2.dot(d2); // Squared length of segment S2, always nonnegative
    double f = d2.dot(r);

    // Check if either or both segments degenerate into points
    if (a <= EPSILON && e <= EPSILON)
    {
        // Both segments degenerate into points
        s = t = 0.0;
        c1 = p1;
        c2 = p2;
        return (c1 - c2).dot(c1 - c2);
    }
    if (a <= EPSILON)
    {
        // First segment degenerates into a point
        s = 0.0;
        t = f / e; // s = 0 => t = (b*s + f) / e = f / e
        t = clamp(t, 0.0, 1.0);
    }
    else
    {
        double c = d1.dot(r);
        if (e <= EPSILON)
        {
            // Second segment degenerates into a point
            t = 0.0;
            s = clamp(-c / a, 0.0, 1.0); // t = 0 => s = (b*t - c) / a = -c / a
        }
        else
        {
            // The general nondegenerate case starts here
            double b = d1.dot(d2);
            double denom = a * e - b * b; // Always nonnegative

            // If segments not parallel, compute closest point on L1 to L2, and
            // clamp to segment S1. Else pick arbitrary s (here 0)
            if (denom != 0.0)
            {
                s = clamp((b * f - c * e) / denom, 0.0, 1.0);
            }
            else
                s = 0.0;

            // Compute point on L2 closest to S1(s) using
            // t = Dot((P1+D1*s)-P2,D2) / Dot(D2,D2) = (b*s + f) / e
            t = (b * s + f) / e;

            // If t in [0,1] done. Else clamp t, recompute s for the new value
            // of t using s = Dot((P2+D2*t)-P1,D1) / Dot(D1,D1)= (t*b - c) / a
            // and clamp s to [0, 1]
            if (t < 0.0)
            {
                t = 0.0;
                s = clamp(-c / a, 0.0, 1.0);
            }
            else if (t > 1.0)
            {
                t = 1.0;
                s = clamp((b - c) / a, 0.0, 1.0);
            }
        }
    }

    c1 = p1 + d1 * s;
    c2 = p2 + d2 * t;
    return (c1 - c2).dot(c1 - c2);
}

// Adapted from Christer Ericson, "Real Time Collision Detection"
// Compute barycentric coordinates (u, v, w) for 
// point p with respect to triangle (a, b, c)
void Barycentric(const Vec3d& a, const Vec3d& b, const Vec3d& c, const Vec3d& p, double& u, double& v, double& w)
{
    const Vec3d v0 = b - a;
    const Vec3d v1 = c - a;
    const Vec3d v2 = p - a;
    double d00 = v0.dot(v0);
    double d01 = v0.dot(v1);
    double d11 = v1.dot(v1);
    double d20 = v2.dot(v0);
    double d21 = v2.dot(v1);
    double denom = d00 * d11 - d01 * d01;
    v = (d11 * d20 - d01 * d21) / denom;
    w = (d00 * d21 - d01 * d20) / denom;
    u = 1.0 - v - w;
}

// Adapted from some code on Robert Bridson's website, I believe

void addUnique(std::vector<double>& a, double e)
{
    for (unsigned int i = 0; i < a.size(); ++i)
        if (a[i] == e)
            return;
    a.push_back(e);
}

double triple(const Vec3d& a, const Vec3d& b, const Vec3d& c)
{
    return a[0] * (b[1] * c[2] - b[2] * c[1]) + a[1] * (b[2] * c[0] - b[0] * c[2]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
}

double signed_volume(const Vec3d& x0, const Vec3d& x1, const Vec3d& x2, const Vec3d& x3)
{
    // Equivalent to triple(x1-x0, x2-x0, x3-x0), six times the signed volume of the tetrahedron.
    // But, for robustness, we want the result (up to sign) to be independent of the ordering.
    // And want it as accurate as possible...
    // But all that stuff is hard, so let's just use the common assumption that all coordinates are >0,
    // and do something reasonably accurate in fp.

    // This formula does almost four times too much multiplication, but if the coordinates are non-negative
    // it suffers in a minimal way from cancellation error.
    return (x0[0] * (x1[1] * x3[2] + x3[1] * x2[2] + x2[1] * x1[2]) + x1[0] * (x2[1] * x3[2] + x3[1] * x0[2] + x0[1] * x2[2])
            + x2[0] * (x3[1] * x1[2] + x1[1] * x0[2] + x0[1] * x3[2]) + x3[0] * (x1[1] * x2[2] + x2[1] * x0[2] + x0[1] * x1[2]))

            - (x0[0] * (x2[1] * x3[2] + x3[1] * x1[2] + x1[1] * x2[2]) + x1[0]
                    * (x3[1] * x2[2] + x2[1] * x0[2] + x0[1] * x3[2]) + x2[0] * (x1[1] * x3[2] + x3[1] * x0[2] + x0[1] * x1[2])
                    + x3[0] * (x2[1] * x1[2] + x1[1] * x0[2] + x0[1] * x2[2]));
}

// All roots returned in interval [0,1]. Assumed geometry followed a linear
// trajectory between x and xnew. 
void getCoplanarityTimes(const Vec3d& x0, const Vec3d& x1, const Vec3d& x2, const Vec3d& x3, const Vec3d& xnew0,
        const Vec3d& xnew1, const Vec3d& xnew2, const Vec3d& xnew3, std::vector<double>& times, std::vector<double>& errors)
{
    const double tol = 1e-8;
    times.clear();
    errors.clear();

    // cubic coefficients, A*t^3+B*t^2+C*t+D (for t in [0,1])
    Vec3d x03 = x0 - x3;
    Vec3d x13 = x1 - x3;
    Vec3d x23 = x2 - x3;
    Vec3d v03 = (xnew0 - xnew3) - x03;
    Vec3d v13 = (xnew1 - xnew3) - x13;
    Vec3d v23 = (xnew2 - xnew3) - x23;

    double A = triple(v03, v13, v23);
    double B = triple(x03, v13, v23) + triple(v03, x13, v23) + triple(v03, v13, x23);
    double C = triple(x03, x13, v23) + triple(x03, v13, x23) + triple(v03, x13, x23);
    double D = triple(x03, x13, x23);

    const double convergence_tol = tol * (std::fabs(A) + std::fabs(B) + std::fabs(C) + std::fabs(D));

    // find intervals to check, or just solve it if it reduces to a quadratic =============================
    std::vector<double> interval_times;
    double discriminant = B * B - 3 * A * C; // of derivative of cubic, 3*A*t^2+2*B*t+C, divided by 4 for convenience
    if (discriminant <= 0)
    { // monotone cubic: only one root in [0,1] possible
        // so we just
        interval_times.push_back(0);
        interval_times.push_back(1);
    }
    else
    { // positive discriminant, B!=0
        if (A == 0)
        { // the cubic is just a quadratic, B*t^2+C*t+D ========================================
            discriminant = C * C - 4 * B * D; // of the quadratic
            if (discriminant <= 0)
            {
                double t = -C / (2 * B);
                if (t >= -tol && t <= 1 + tol)
                {
                    t = clamp(t, 0., 1.);
                    double val = std::fabs(
                            signed_volume((1 - t) * x0 + t * xnew0, (1 - t) * x1 + t * xnew1, (1 - t) * x2 + t * xnew2,
                                    (1 - t) * x3 + t * xnew3));
                    if (val < convergence_tol)
                    {
                        times.push_back(t);
                    }
                }
            }
            else
            { // two separate real roots
                double t0, t1;
                if (C > 0)
                    t0 = (-C - std::sqrt(discriminant)) / (2 * B);
                else
                    t0 = (-C + std::sqrt(discriminant)) / (2 * B);
                t1 = D / (B * t0);
                if (t1 < t0)
                    std::swap(t0, t1);
                if (t0 >= -tol && t0 <= 1 + tol)
                {
                    times.push_back(clamp(t0, 0., 1.));
                }
                if (t1 >= -tol && t1 <= 1 + tol)
                {
                    addUnique(times, clamp(t1, 0., 1.));
                }
            }

            for (int i = 0; i < (int) times.size(); ++i)
            {
                double ti = times[i];
                double val = std::fabs(
                        signed_volume((1 - ti) * x0 + ti * xnew0, (1 - ti) * x1 + ti * xnew1, (1 - ti) * x2 + ti * xnew2,
                                (1 - ti) * x3 + ti * xnew3));
                errors.push_back(val);
            }

            return;
        }
        else
        { // cubic is not monotone: divide up [0,1] accordingly =====================================
            double t0, t1;
            if (B > 0)
                t0 = (-B - std::sqrt(discriminant)) / (3 * A);
            else
                t0 = (-B + std::sqrt(discriminant)) / (3 * A);
            t1 = C / (3 * A * t0);
            if (t1 < t0)
                std::swap(t0, t1);
            interval_times.push_back(0);
            if (t0 > 0 && t0 < 1)
                interval_times.push_back(t0);
            if (t1 > 0 && t1 < 1)
                interval_times.push_back(t1);

            interval_times.push_back(1);
        }
    }

    // look for roots in indicated intervals ==============================================================
    // evaluate coplanarity more accurately at each endpoint of the intervals
    std::vector<double> interval_values(interval_times.size());
    for (unsigned int i = 0; i < interval_times.size(); ++i)
    {
        double t = interval_times[i];
        interval_values[i] = signed_volume((1 - t) * x0 + t * xnew0, (1 - t) * x1 + t * xnew1, (1 - t) * x2 + t * xnew2,
                (1 - t) * x3 + t * xnew3);
    }
    // first look for interval endpoints that are close enough to zero, without a sign change
    for (unsigned int i = 0; i < interval_times.size(); ++i)
    {
        if (interval_values[i] == 0)
        {
            times.push_back(interval_times[i]);
        }
        else if (std::fabs(interval_values[i]) < convergence_tol)
        {
            if ((i == 0 || (interval_values[i - 1] >= 0 && interval_values[i] >= 0) || (interval_values[i - 1] <= 0
                    && interval_values[i] <= 0)) && (i == interval_times.size() - 1 || (interval_values[i + 1] >= 0
                    && interval_values[i] >= 0) || (interval_values[i + 1] <= 0 && interval_values[i] <= 0)))
            {
                times.push_back(interval_times[i]);
            }
        }
    }
    // and then search in intervals with a sign change
    for (unsigned int i = 1; i < interval_times.size(); ++i)
    {
        double tlo = interval_times[i - 1], thi = interval_times[i], tmid;
        double vlo = interval_values[i - 1], vhi = interval_values[i], vmid;
        if ((vlo < 0 && vhi > 0) || (vlo > 0 && vhi < 0))
        {
            // start off with secant approximation (in case the cubic is actually linear)
            double alpha = vhi / (vhi - vlo);
            tmid = alpha * tlo + (1 - alpha) * thi;
            for (int iteration = 0; iteration < 50; ++iteration)
            {
                vmid = signed_volume((1 - tmid) * x0 + tmid * xnew0, (1 - tmid) * x1 + tmid * xnew1,
                        (1 - tmid) * x2 + tmid * xnew2, (1 - tmid) * x3 + tmid * xnew3);
                if (std::fabs(vmid) < 1e-2 * convergence_tol)
                    break;
                if ((vlo < 0 && vmid > 0) || (vlo > 0 && vmid < 0))
                { // if sign change between lo and mid
                    thi = tmid;
                    vhi = vmid;
                }
                else
                { // otherwise sign change between hi and mid
                    tlo = tmid;
                    vlo = vmid;
                }
                if (iteration % 2)
                    alpha = 0.5; // sometimes go with bisection to guarantee we make progress
                else
                    alpha = vhi / (vhi - vlo); // other times go with secant to hopefully get there fast
                tmid = alpha * tlo + (1 - alpha) * thi;
            }
            times.push_back(tmid);
        }
    }
    std::sort(times.begin(), times.end());

    for (int i = 0; i < (int) times.size(); ++i)
    {
        double ti = times[i];
        double val = std::fabs(
                signed_volume((1 - ti) * x0 + ti * xnew0, (1 - ti) * x1 + ti * xnew1, (1 - ti) * x2 + ti * xnew2,
                        (1 - ti) * x3 + ti * xnew3));
        errors.push_back(val);
    }
}

void getIntersectionPoint(const Vec3d& x_edge_0, const Vec3d& x_edge_1, const Vec3d& x_face_0, const Vec3d& x_face_1,
        const Vec3d& x_face_2, std::vector<double>& times, std::vector<double>& errors)
{
    const double tol = 1e-12;
    times.clear();
    errors.clear();

    Vec3d x03 = x_edge_0 - x_face_2;
    Vec3d x13 = x_face_0 - x_face_2;
    Vec3d x23 = x_face_1 - x_face_2;
    Vec3d v03 = x_edge_1 - x_face_2 - x03;

    double C = triple(v03, x13, x23);
    double D = triple(x03, x13, x23);

    const double convergence_tol = tol * (std::fabs(0) + std::fabs(0) + std::fabs(C) + std::fabs(D));

    // find intervals to check, or just solve it if it reduces to a quadratic =============================
    std::vector<double> interval_times;
    interval_times.push_back(0);
    interval_times.push_back(1);

    // look for roots in indicated intervals ==============================================================
    // evaluate coplanarity more accurately at each endpoint of the intervals
    std::vector<double> interval_values(interval_times.size());
    for (unsigned int i = 0; i < interval_times.size(); ++i)
    {
        double t = interval_times[i];
        interval_values[i] = signed_volume((1 - t) * x_edge_0 + t * x_edge_1, x_face_0, x_face_1, x_face_2);
    }
    // first look for interval endpoints that are close enough to zero, without a sign change
    for (unsigned int i = 0; i < interval_times.size(); ++i)
    {
        if (interval_values[i] == 0)
        {
            times.push_back(interval_times[i]);
        }
        else if (std::fabs(interval_values[i]) < convergence_tol)
        {
            if ((i == 0 || (interval_values[i - 1] >= 0 && interval_values[i] >= 0) || (interval_values[i - 1] <= 0
                    && interval_values[i] <= 0)) && (i == interval_times.size() - 1 || (interval_values[i + 1] >= 0
                    && interval_values[i] >= 0) || (interval_values[i + 1] <= 0 && interval_values[i] <= 0)))
            {
                times.push_back(interval_times[i]);
            }
        }
    }
    // and then search in intervals with a sign change
    for (unsigned int i = 1; i < interval_times.size(); ++i)
    {
        double tlo = interval_times[i - 1], thi = interval_times[i], tmid;
        double vlo = interval_values[i - 1], vhi = interval_values[i], vmid;
        if ((vlo < 0 && vhi > 0) || (vlo > 0 && vhi < 0))
        {
            // start off with secant approximation (in case the cubic is actually linear)
            double alpha = vhi / (vhi - vlo);
            tmid = alpha * tlo + (1 - alpha) * thi;
            for (int iteration = 0; iteration < 50; ++iteration)
            {
                vmid = signed_volume((1 - tmid) * x_edge_0 + tmid * x_edge_1, (1 - tmid) * x_face_0 + tmid * x_face_0,
                        (1 - tmid) * x_face_1 + tmid * x_face_1, (1 - tmid) * x_face_2 + tmid * x_face_2);
                if (std::fabs(vmid) < 1e-2 * convergence_tol)
                    break;
                if ((vlo < 0 && vmid > 0) || (vlo > 0 && vmid < 0))
                { // if sign change between lo and mid
                    thi = tmid;
                    vhi = vmid;
                }
                else
                { // otherwise sign change between hi and mid
                    tlo = tmid;
                    vlo = vmid;
                }
                if (iteration % 2)
                    alpha = 0.5; // sometimes go with bisection to guarantee we make progress
                else
                    alpha = vhi / (vhi - vlo); // other times go with secant to hopefully get there fast
                tmid = alpha * tlo + (1 - alpha) * thi;
            }
            times.push_back(tmid);
        }
    }
    std::sort(times.begin(), times.end());

    for (int i = 0; i < (int) times.size(); ++i)
    {
        double ti = times[i];
        double val = std::fabs(signed_volume((1 - ti) * x_edge_0 + ti * x_edge_1, x_face_0, x_face_1, x_face_2));
        errors.push_back(val);
    }
}

}
