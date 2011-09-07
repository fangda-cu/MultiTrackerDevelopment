#include <iostream>

#include "BASim/src/Physics/DeformableObjects/Shells/BendingDerivatives.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/DSBendingForce.hh"
#include "BASim/src/Math/Math.hh"

namespace BASim {
/** ComputeAngle
 *
 *         A
 *        /
 *       /
 *      B------C
 *
 *  Returns the angle ABC (between -PI and PI) given the orienting
 *  axis n
 */


//does Eigen not have a 2-parameter version of these? it's way more readable.
static Scalar len(const Vec3d& vec) { return sqrt(vec.dot(vec)); }
static Scalar lenSq(const Vec3d& vec) { return vec.dot(vec); }
static Vec3d cross(const Vec3d& a, const Vec3d& b) { return a.cross(b); }
static Scalar dot(const Vec3d& v1, const Vec3d& v2) { return v1.dot(v2); }
static Vec3d dir(const Vec3d& c) { Scalar l = len(c); assert(l != 0); return c/l; }

static void setTranspose(Mat3d& a, const Mat3d& b) {
  for(int i = 0; i < 3; ++i) for(int j = 0; j < 3; ++j)
    a(i,j) = b(j,i);
}

Scalar DSBendingForce::ComputeAngle(Vec3d BA, Vec3d BC) const
{
    BA = BA / len(BA);
    BC = BC / len(BC);

    Vec3d BA_cross_BC = cross(BA, BC);

    Scalar cosABC = dot(BA, BC);
    Scalar sinABC = len(BA_cross_BC);

    Scalar ABC = atan2(sinABC, cosABC);

    return ABC;
}

// ******************** First Derivative Code **********************



/** ComputeDihedralAngleDerivatives
 *
 *  Derivative of the dihedral angle (theta) wrt to the four points
 *  that form it
 *
 *         q1
 *         /\
 *        /  \
 *      p1----p2
 *        \  /
 *         \/
 *         q2
 */
void DSBendingForce::ComputeDihedralAngleDerivatives(Vec3d& del_p1_theta,
                                                Vec3d& del_p2_theta,
                                                Vec3d& del_q1_theta,
                                                Vec3d& del_q2_theta,
                                                const Vec3d& p1,
                                                const Vec3d& p2,
                                                const Vec3d& q1,
                                                const Vec3d& q2) const
{
    // edge vectors
    Vec3d v11 = q1 - p2;
    Vec3d v12 = q1 - p1;

    Vec3d v21 = q2 - p2;
    Vec3d v22 = q2 - p1;

    Vec3d v = p2 - p1;

    Scalar vNorm2 = dot(v, v);
    Scalar vNorm  = sqrt(vNorm2);

    // normal vectors
    Vec3d n1 = cross(v12, v);
    Vec3d n2 = cross(v, v22);

    Scalar n1Norm2 = dot(n1, n1);
    Scalar n2Norm2 = dot(n2, n2);

    // gradients of theta
    del_q1_theta = vNorm * n1 / n1Norm2;
    del_q2_theta = vNorm * n2 / n2Norm2;

    Vec3d F11 =   dot(v11, v) / (vNorm * n1Norm2) * n1;
    Vec3d F12 =   dot(v21, v) / (vNorm * n2Norm2) * n2;

    Vec3d F21 = - dot(v12, v) / (vNorm * n1Norm2) * n1;
    Vec3d F22 = - dot(v22, v) / (vNorm * n2Norm2) * n2;

    del_p1_theta = F11 + F12;
    del_p2_theta = F21 + F22;

    assert(len(del_p1_theta + del_p2_theta + del_q1_theta + del_q2_theta) < 1e-4);
    assert(len(cross(p1, del_p1_theta) + cross(p2, del_p2_theta) + cross(q1, del_q1_theta) + cross(q2, del_q2_theta)) < 1e-4);
}

// ******************** End First Derivative Code ********************



// ********************* Second Derivative Code **********************

/** Symmetrize
 *
 *  sets m to (m + m^T)
 */
void DSBendingForce::Symmetrize(Mat3d& m) const
{
    m(0, 0) += m(0, 0);
    m(1, 1) += m(1, 1);
    m(2, 2) += m(2, 2);

    const Scalar x01 = m(0, 1) + m(1, 0);

    const Scalar x02 = m(0, 2) + m(2, 0);

    const Scalar x12 = m(1, 2) + m(2, 1);

    m(0, 1) = x01;

    m(1, 0) = x01;

    m(0, 2) = x02;

    m(2, 0) = x02;

    m(1, 2) = x12;

    m(2, 1) = x12;
}

/** ComputeDihedralAngleSecondDerivatives
 *
 *  double derivative of the dihedral angle wrt the four points that
 *  create it
 *
 *         q1
 *         /\
 *        /  \
 *      p1----p2
 *        \  /
 *         \/
 *         q2
 */
void DSBendingForce::ComputeDihedralAngleSecondDerivatives(EnergyHessian& J,
                                                      Scalar Kb,
                                                      Vec3d& p1,
                                                      Vec3d& p2,
                                                      Vec3d& q1,
                                                      Vec3d& q2,
                                                      int p1Index,
                                                      int p2Index,
                                                      int q1Index,
                                                      int q2Index) const
{
    Scalar Dtheta_phi = Kb * 1.;

    // ************** edge vectors **************
    Vec3d v11 = q1 - p2;
    Vec3d v12 = q1 - p1;

    Vec3d v21 = q2 - p2;
    Vec3d v22 = q2 - p1;

    Vec3d v = p2 - p1;

    Scalar vNorm2 = dot(v, v);
    Scalar vNorm  = sqrt(vNorm2);


    // ************* normal vectors **************
    Vec3d n1 = cross(v12, v);
    Vec3d n2 = cross(v, v22);

    Scalar n1Norm2 = dot(n1, n1);
    Scalar n2Norm2 = dot(n2, n2);

    Scalar n1Norm  = sqrt(n1Norm2);
    Scalar n2Norm  = sqrt(n2Norm2);

    Scalar n1Norm3 = n1Norm * n1Norm2;
    Scalar n2Norm3 = n2Norm * n2Norm2;

    Vec3d n1hat = (1. / n1Norm) * n1;
    Vec3d n2hat = (1. / n2Norm) * n2;


    // ******* Derivative computation begins ******
    Vec3d t11 = cross(n1hat, v11);
    Vec3d t12 = cross(v12, n1hat);

    Vec3d t21 = cross(v21, n2hat);
    Vec3d t22 = cross(n2hat, v22);

    Vec3d t1  = cross(n1hat, v);
    Vec3d t2  = cross(v, n2hat);


    Mat3d t11n1Sym = outerProd(t11, n1);
    Mat3d t12n1Sym = outerProd(t12, n1);
    Mat3d t21n2Sym = outerProd(t21, n2);
    Mat3d t22n2Sym = outerProd(t22, n2);
    Symmetrize(t11n1Sym);  // set t11n1Sym to t11^T n1 + n1^T t11
    Symmetrize(t12n1Sym);
    Symmetrize(t21n2Sym);
    Symmetrize(t22n2Sym);

    Mat3d n1v = outerProd(n1, v);
    Mat3d n2v = outerProd(n2, v);

    Mat3d n1t1 = outerProd(n1, t1);
    Mat3d n2t2 = outerProd(n2, t2);

    Mat3d n1t1Sym = n1t1;
    Mat3d n2t2Sym = n2t2;
    Symmetrize(n1t1Sym);
    Symmetrize(n2t2Sym);

    Scalar vNorm_n1Norm2_Inv = (1. / (vNorm * n1Norm2));
    Scalar vNorm_n2Norm2_Inv = (1. / (vNorm * n2Norm2));

    Mat3d Dq1_Dp1_theta =   vNorm_n1Norm2_Inv * n1v - (vNorm / n1Norm3) * t11n1Sym;
    Mat3d Dq1_Dp2_theta = - vNorm_n1Norm2_Inv * n1v - (vNorm / n1Norm3) * t12n1Sym;
    J.Add(q1Index, p1Index, Dtheta_phi * Dq1_Dp1_theta);
    J.Add(q1Index, p2Index, Dtheta_phi * Dq1_Dp2_theta);

    Mat3d Dp1_Dq1_theta, Dp2_Dq1_theta;
    setTranspose(Dp1_Dq1_theta,Dq1_Dp1_theta);
    setTranspose(Dp2_Dq1_theta,Dq1_Dp2_theta);
    J.Add(p1Index, q1Index, Dtheta_phi * Dp1_Dq1_theta);
    J.Add(p2Index, q1Index, Dtheta_phi * Dp2_Dq1_theta);

    Mat3d Dq2_Dp1_theta =   vNorm_n2Norm2_Inv * n2v - (vNorm / n2Norm3) * t21n2Sym;
    Mat3d Dq2_Dp2_theta = - vNorm_n2Norm2_Inv * n2v - (vNorm / n2Norm3) * t22n2Sym;
    J.Add(q2Index, p1Index, Dtheta_phi * Dq2_Dp1_theta);
    J.Add(q2Index, p2Index, Dtheta_phi * Dq2_Dp2_theta);

    Mat3d Dp1_Dq2_theta, Dp2_Dq2_theta;
    setTranspose(Dp1_Dq2_theta,Dq2_Dp1_theta);
    setTranspose(Dp2_Dq2_theta,Dq2_Dp2_theta);
    J.Add(p1Index, q2Index, Dtheta_phi * Dp1_Dq2_theta);
    J.Add(p2Index, q2Index, Dtheta_phi * Dp2_Dq2_theta);

    Mat3d Dq1_Dq1_theta = - (vNorm / n1Norm3) * n1t1Sym;
    Mat3d Dq2_Dq2_theta = - (vNorm / n2Norm3) * n2t2Sym;
    J.Add(q1Index, q1Index, Dtheta_phi * Dq1_Dq1_theta);
    J.Add(q2Index, q2Index, Dtheta_phi * Dq2_Dq2_theta);

    //Mat3d A11 = outerProd(n1, - v  - (n1Norm/vNorm2)*t1);
    //Mat3d A12 = outerProd(n2, - v  - (n2Norm/vNorm2)*t2);
    Mat3d A21 = outerProd(n1,   v  - (n1Norm / vNorm2) * t1);
    Mat3d A22 = outerProd(n2,   v  - (n2Norm / vNorm2) * t2);

    Mat3d Dp1_F11 = vNorm_n1Norm2_Inv * ((-dot(v11, v) / n1Norm) * t11n1Sym - (n1Norm / vNorm2) * n1t1);
    Mat3d Dp1_F12 = vNorm_n2Norm2_Inv * ((-dot(v21, v) / n2Norm) * t21n2Sym - (n2Norm / vNorm2) * n2t2);

    Mat3d Dp1_F21 = vNorm_n1Norm2_Inv * ((dot(v12, v) / n1Norm) * t11n1Sym - A21);
    Mat3d Dp1_F22 = vNorm_n2Norm2_Inv * ((dot(v22, v) / n2Norm) * t21n2Sym - A22);

    //Mat3d Dp2_F11 = vNorm_n1Norm2_Inv * ((-dot(v11,v)/n1Norm)*t12n1Sym - A11);
    //Mat3d Dp2_F12 = vNorm_n2Norm2_Inv * ((-dot(v21,v)/n2Norm)*t22n2Sym - A12);

    Mat3d Dp2_F21 = vNorm_n1Norm2_Inv * ((dot(v12, v) / n1Norm) * t12n1Sym - (n1Norm / vNorm2) * n1t1);
    Mat3d Dp2_F22 = vNorm_n2Norm2_Inv * ((dot(v22, v) / n2Norm) * t22n2Sym - (n2Norm / vNorm2) * n2t2);

    // Dpi_Dpj_theta force acting on pj, derivative w.r.t. pi
    Mat3d Dp1_Dp1_theta = Dp1_F11 + Dp1_F12;
    Mat3d Dp1_Dp2_theta = Dp1_F21 + Dp1_F22;
    Mat3d Dp2_Dp2_theta = Dp2_F21 + Dp2_F22;

    //Mat3d Dp2_Dp1_theta = Dp2_F11 + Dp2_F12;

    Mat3d Dp2_Dp1_theta;
    setTranspose(Dp2_Dp1_theta,Dp1_Dp2_theta);

    J.Add(p1Index, p1Index, Dtheta_phi * Dp1_Dp1_theta);
    J.Add(p2Index, p2Index, Dtheta_phi * Dp2_Dp2_theta);

    J.Add(p1Index, p2Index, Dtheta_phi * Dp2_Dp1_theta);
    J.Add(p2Index, p1Index, Dtheta_phi * Dp1_Dp2_theta);
}

}
// ******************** End Second Derivative Code **********************
