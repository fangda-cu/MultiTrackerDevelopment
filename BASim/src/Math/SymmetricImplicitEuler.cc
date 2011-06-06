/*
 * SymmetricImplicitEuler.cc
 *
 *  Created on: 6/06/2011
 *      Author: jaubry
 */

typedef double Scalar;
#include <string>
#include "SymmetricImplicitEuler.hh"
#include "../Physics/ElasticRods/MultipleRodTimeStepper.hh"
#include "../Physics/ElasticRods/RodTimeStepper.hh"

namespace BASim
{

// Initial guess based on rigid motion of the first two vertices
template<>
void SymmetricImplicitEuler<RodTimeStepper>::generateInitialIterate0(VecXd& dx)
{
    const Vec3d p0 = x0.segment<3> (0);
    const Vec3d p1 = x0.segment<3> (4);
    const Vec3d w0 = v0.segment<3> (0);
    const Vec3d w1 = v0.segment<3> (4);
    const Vec3d q0 = p0 + m_dt * w0;
    const Vec3d q1 = p1 + m_dt * w1;

    const double cosAngle = (p1 - p0).dot(q1 - q0) / ((p1 - p0).norm() * (q1 - q0).norm());
    const double angle = cosAngle >= 1.0 ? 0.0 : cosAngle <= -1.0 ? M_PI : acos(cosAngle);
    assert(!isnan(angle));

    Vec3d normal = (p1 - p0).cross(w1 - w0);
    const double normalNorm = normal.norm();
    if (normalNorm < std::numeric_limits<double>::epsilon())
        normal = Vec3d(1, 0, 0); // The rotation is either identity or central symmetry, any axis will do...
    else
        normal = normal / normalNorm;
    assert(approxEq(normal.norm(), 1.0));
    TraceStream(g_log, "") << "Initial guess by rigid motion: normal = " << normal << " angle = " << angle << '\n';

    const Eigen::AngleAxis<double> rotation(angle, normal);
    for (int i = 0; i < m_ndof; i += 4)
    {
        dx.segment<3> (i) = m_dt * w0 + rotation._transformVector(x0.segment<3> (i) - p0);
        // dx(i + 3) = 0;
    }
}

template<>
void SymmetricImplicitEuler<MultipleRodTimeStepper>::generateInitialIterate0(VecXd& dx)
{
    ErrorStream(g_log, "")
            << "SymmetricImplicitEuler<MultipleRodTimeStepper>::generateInitialIterate0 not implemented yet for this ODE\n";
    assert(0);
}

}
