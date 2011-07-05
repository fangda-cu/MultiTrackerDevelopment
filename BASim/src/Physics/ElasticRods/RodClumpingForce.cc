/*
 * RodClumpingForce.cc
 *
 *  Created on: 30/06/2011
 *      Author: Jean-Marie Aubry <jaubry@wetafx.co.nz>
 */

#include "../../Util/TextLog.hh"
#include "../../Collisions/CollisionUtils.hh"
#include "RodClumpingForce.hh"
#include <limits>

namespace BASim
{

RodClumpingForce::RodClumpingForce() :
    q(0.0), r(1.0), rho2(0.01)
{
    m_name = "clumping";
}

RodClumpingForce::~RodClumpingForce()
{
    // TODO Auto-generated destructor stub
}

Vec3d ClosestPointOnRod(const Vec3d& x, const ElasticRod& other)
{
    Scalar mindist = std::numeric_limits<Scalar>::max();
    Vec3d winner;

    for (int othervidx = 0; othervidx < other.nv() - 1; ++othervidx)
    {
        Vec3d y = ClosestPtPointSegment(x, other.getVertex(othervidx), other.getVertex(othervidx + 1));
        Scalar dist = (y - x).norm();
        if (dist < mindist)
        {
            mindist = dist;
            winner = y;
        }
    }

    return winner;
}

Scalar RodClumpingForce::computeEnergy(const ElasticRod& rod) const
{
    const std::vector<ElasticRod*>& neighbours = rod.getNearestRootNeighbours();
    Scalar energy = 0.0;

    for (std::vector<ElasticRod*>::const_iterator neighbour = neighbours.begin(); neighbour != neighbours.end(); ++neighbour)
        energy += computeEnergy(rod, *(*neighbour));

    return energy;
}

Scalar RodClumpingForce::computeEnergy(const ElasticRod& rod, const ElasticRod& other) const
{
    Scalar energy = 0.0;
    for (int vidx = 0; vidx < rod.nv(); ++vidx)
    {
        const Vec3d& x = rod.getVertex(vidx);

        const Vec3d& y = ClosestPointOnRod(x, other);
        const Scalar normep2 = (x[0] - y[0]) * (x[0] - y[0]) + (x[1] - y[1]) * (x[1] - y[1]) + (x[2] - y[2]) * (x[2] - y[2]);
        energy -= q * normep2 / pow(normep2 + rho2, 1 + 0.5 * r);
    }

    return energy;
}

void RodClumpingForce::computeForce(const ElasticRod& rod, VecXd& force) const
{
    const std::vector<ElasticRod*>& neighbours = rod.getNearestRootNeighbours();

    for (std::vector<ElasticRod*>::const_iterator neighbour = neighbours.begin(); neighbour != neighbours.end(); ++neighbour)
        computeForce(rod, *(*neighbour), force);
}

void RodClumpingForce::computeForce(const ElasticRod& rod, const ElasticRod& other, VecXd& force) const
{
    for (int vidx = 0; vidx < rod.nv(); ++vidx)
    {
        const Vec3d& x = rod.getVertex(vidx);
        Vec3d localForceFromOther;
        localForceFromOther.setZero();

        const Vec3d& y = ClosestPointOnRod(x, other);
        const Scalar normep2 = (x[0] - y[0]) * (x[0] - y[0]) + (x[1] - y[1]) * (x[1] - y[1]) + (x[2] - y[2]) * (x[2] - y[2]);

        Vec3d localForceFromOtherVtx = q * (y - x) * pow(normep2 + rho2, -1 - 0.5 * r) * (2 - (2 + r) * normep2 / (normep2
                + rho2));
        localForceFromOther += localForceFromOtherVtx;

        TraceStream(g_log, "") << "Force on vertex " << vidx << " of rod " << &rod << " from rod " << &other << ": "
                << localForceFromOther << '\n';
        force.segment<3> (rod.vertIdx(vidx, 0)) += localForceFromOther;

    }
}

void RodClumpingForce::computeForceEnergy(const ElasticRod& rod, VecXd& force, Scalar& energy) const
{
    computeForce(rod, force);
    energy += computeEnergy(rod);
}

void RodClumpingForce::computeForceDX(int baseindex, const ElasticRod& rod, Scalar scale, MatrixBase& J) const
{
    //TraceStream(g_log, "") << "Computing global clumping Jacobian for rod " << rod.globalRodIndex << ": before " << J << '\n';

    const std::vector<ElasticRod*>& neighbours = rod.getNearestRootNeighbours();

    for (std::vector<ElasticRod*>::const_iterator neighbour = neighbours.begin(); neighbour != neighbours.end(); ++neighbour)
        computeForceDX(baseindex, rod, *(*neighbour), scale, J);

    //    TraceStream(g_log, "") << "Computing global clumping Jacobian for rod " << rod.globalRodIndex << ": after " << J << '\n';
}

void RodClumpingForce::computeForceDX(int baseindex, const ElasticRod& rod, const ElasticRod& other, Scalar scale,
        MatrixBase& J) const
{
    Mat3d localJ;
    for (int vidx = 0; vidx < rod.nv(); ++vidx)
    {
        computeLocalForceDX(rod, vidx, other, localJ);
        //       TraceStream(g_log, "") << "Local clumping force Jacobian for rod " << rod.globalRodIndex << " vertex " << vidx
        //               << " from rod " << other.globalRodIndex << ": " << localJ << '\n';
        J.pointStencilAdd(rod.vertIdx(vidx, 0) + baseindex, localJ);
    }
}

void RodClumpingForce::computeLocalForceDX(const ElasticRod& rod, int vidx, const ElasticRod& other, Mat3d& localJ) const
{
    localJ.setZero();
    const Vec3d& x = rod.getVertex(vidx);
    const Vec3d& y = ClosestPointOnRod(x, other);
    const Scalar normep2 = (x[0] - y[0]) * (x[0] - y[0]) + (x[1] - y[1]) * (x[1] - y[1]) + (x[2] - y[2]) * (x[2] - y[2]);
    //  const Scalar normepr4 = pow(normep2 + 0.1, 0.5 * r + 2);
    const Scalar normep4 = normep2 * normep2;

    for (int i = 0; i < 3; i++)
    {
        const double xx2 = (x[i] - y[i]) * (x[i] - y[i]);
        localJ(i, i) = q * (pow(normep2 + rho2, -3 - r / 2.) * (normep4 * r - 2 * rho2 * (-2 * (2 + r) * xx2 + rho2) - normep2
                * (2 * r * xx2 + r * r * xx2 + 2 * rho2 - r * rho2)));

        for (int j = 0; j < i; j++)
            localJ(i, j) = localJ(j, i) = q * ((-2 - r) * (x[i] - y[i]) * (x[j] - y[j]) * (normep2 * r - 4 * rho2) * pow(
                    normep2 + rho2, -3 - r / 2.));
    }
}

void RodClumpingForce::computeForceDV(int baseindex, const ElasticRod& rod, Scalar scale, MatrixBase& J) const
{
}

}
